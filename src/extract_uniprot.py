import os
import json
import sys
import time
import gzip
import io
import re
import pickle
import pymongo

try:
    import psutil
except:
    pass

import numpy as np

pid = os.getpid()

class data_extractor_base:
    """Base class defining the interface for data extraction from uniprot 
    entries in text format. Children must implement :func:`id` and 
    :func:`extract`.
    """
    def __init__(self):
        pass
    def id(self):
        """Must be implemented by child class and returns a string identifier of
        the extracted data.
        """
        raise NotImplementedError('id not implemented')

    def extract(self, entry):
        """Must be implemented by child class and returns the extracted data,
        None if not available.

        :param entry:  Content of uniprot entry in text format, e.g. the content
                       of `P69892 <https://www.uniprot.org/uniprot/P69892.txt>`_ .
        :type entry:  :class:`list` of :class:`str`
        """
        raise NotImplementedError('extract not implemented')

class uniprot_extractor:
    """This class is capable if iterating over full uniprot releases and 
    extract data using :class:`data_extractor_base` objects which can
    be registered with :class:`register`. Data extraction happens upon
    calling :func:`extract`.
    
    """
    def __init__(self, mongo_host = None, mongo_port = None, name = 'UniProt'):

        if mongo_host is None or mongo_port is None:
            self.data = list()
            self.entries = list()
            self.type = 'List'
        else:
            self.client = pymongo.MongoClient(mongo_host, mongo_port)
            self.db = self.client['Joana'] 
            self.col = self.db[name]
            self.type = 'Mongo'
            self.name = name

        self.data_extractors = list()
        self.ac_extractor = ac_extractor()
        self.coverages_processor = coverages_processor()
        self.n_entries = 0
        self.item_count = 0

    def clear(self):
        """
        Clean the database. It becomes an empty set of lists.
        """
        if self.type == 'Mongo':
            self.col.drop()
            self.col = self.db[self.name]

        elif self.type == 'List':
            self.data = list()
            self.entries = list()

    def register(self, extractor):
        """Register new data extractor. In the data extraction phase, all 
        registered extractors are applied on each processed uniprot entry.

        :param extractor: Extractor to be registered
        :type extractor: data_extractor_base
        """
        self.data_extractors.append(extractor)

    def extract(self, uniprot_data, max_size = 'all', add_if_empty=True, clear = True, targets = 'all', chunk_size = 1000, print_step = 100000, saveto = None, savestep = 100000, interpro_db = None, process_coverages = True):
        """Process all entries in file(s) specified in *uniprot_data* with 
        registered data extractors.
        The data stored for a single entry is a dict with the return value of
        :func:`data_extractor_base.id` as keys and the return value of
        :func:`data_extractor_base.extract` as values. None values are not 
        added. Once data extraction for an entry is complete, storing 
        happens with :func:`store`.

        :param uniprot_data: File(s) containing uniprot entries in text format
                             separated by //. Files ending with .gz are dealt
                             with accordingly. If a list of files is provided,
                             the content of the files is concatenated when
                             processing.
        :param add_if_empty: None values returned by the data extractors are not 
                             added to the per-entry data dicts. Thus, they might
                             be empty. If set to False, these dicts are ignored,
                             i.e. :func:`store` is not called for those.

        :type uniprot_data: :class:`str` / :class:`list` of :class:`str`
        :type add_if_empty: :class:`bool`
        """

        if clear:
            self.clear()

        self.chunk_size = chunk_size
        self.n_chunks = 0
        self.curr_chunk_size = 0
        self.data_dict = list()

        for entry in self._uniprot_entry_iterator(uniprot_data):
            self.item_count += 1
            data = None

            uniprot_ac = self.ac_extractor.extract(entry)

            if uniprot_ac is None:
                raise RuntimeError('Observed None uniprot AC')

            if targets == 'all' or uniprot_ac in targets:
                entry_data = dict()
                for extractor in self.data_extractors:
                    data = extractor.extract(entry)
                    if data is not None:
                        entry_data[extractor.id()] = data
                
                if len(entry_data) > 0 or add_if_empty:    
                    if process_coverages:
                        entry_data = self.coverages_processor.extract(entry_data, uniprot_ac, interpro_db) 
                        
                    self.store(uniprot_ac, entry_data)

                if self.n_entries % print_step == 0:
                    if targets == 'all':
                        print('UNIPROT:', uniprot_ac, self.item_count, self.n_entries, 'RSS memory used (GB):', round(psutil.Process(pid).memory_info().rss/1024/1024/1024, 3))
                    else:
                        print('UNIPROT:', self.n_entries, 'out of {} targets'.format(len(targets)), 'RSS memory used (GB):', round(psutil.Process(pid).memory_info().rss/1024/1024/1024, 3))
                
            if max_size != 'all' and self.item_count >= max_size:
                break
            elif targets != 'all' and self.n_entries == len(targets):
                break

            if saveto is not None and self.item_count % savestep == 0:
                self.save(saveto, self.item_count)

        if saveto is not None:
            self.save(saveto, self.item_count)

        if self.type == 'Mongo':
            self.store(uniprot_ac, data, end = True)

        print('UNIPROT:', self.item_count, self.n_entries, 'RSS memory used (GB):', round(psutil.Process(pid).memory_info().rss/1024/1024/1024, 3))
        
    def store(self, ac, data, end = False):

        """Stores per-entry data in internal data structure. It defines two lists:
        one where the data per entry is stored, and another where the accession codes
        are stored. 
        dict which 
        If you want to use more fancy storage
        systems, e.g. a database, you can subclass :class:`uniparc_extractor`
        and implement your own store.

        :param ac: Uniparc accession code.
        :param data: Data do be stored
        :type ac: :class:`str`
        :type data: :class:`dict`
        """

        if self.type == 'List':
            self.data.append(gzip.compress(json.dumps(data).encode()))
            self.entries.append(ac)

        elif self.type == 'Mongo':
            if not end:
                self.data_dict.append({'_id': ac, 'data': data})
                self.curr_chunk_size += 1

            if self.curr_chunk_size == self.chunk_size or end:
                try:
                    self.col.insert_many(self.data_dict)
                    self.data_dict = list()
                    self.n_chunks += 1
                    self.curr_chunk_size = 0
                except:
                    pass
                    
        self.last_ac = ac
        self.n_entries += 1
    
                
    ### METHODS FOR WHEN THE TYPE OF DB IS MONGO

    def query(self, ac):
        return self.col.find({ '_id': ac })

    def index_db(self):
        self.col.create_index('_id')
        
    ### METHODS FOR WHEN THE TYPE OF DB IS LIST
    
    def save(self, saveto, curr_count, clean = True):

        self.checkpoint = '{}_{}.obj'.format(saveto, curr_count)
        self.saved_index = curr_count

        print('\n ... Saving to: {}'.format(self.checkpoint))
        pickle.dump(self, open('{}'.format(self.checkpoint), 'wb'))

        self.save_index(saveto)

        if clean:
            self.entries = list()
            self.data = list()

    def save_index(self, saveto):

        print(' ... Saving indexes to: {}\n'.format('{}.INDEX'.format(saveto)))
        saveto = '{}.INDEX'.format(saveto)
        if os.path.exists(saveto):
            append_write = 'a' # append if already exists
        else:
            append_write = 'w' # make a new file if not

        with open(saveto, append_write) as outf:
            for uniprot_ac in self.entries:
                outf.write('{}\t{}\n'.format(uniprot_ac, self.saved_index))

        outf.close()
        
    # HELPING METHODS
    
    def _uniprot_entry_iterator(self, uniprot_data, max_entry_size=1048576):
        # check input, uniprot_data must either be string or list of strings 
        # referring to existing files
        files = None
        if isinstance(uniprot_data, str):
            if os.path.isfile(uniprot_data):
                files = [uniprot_data]
            else:
                raise RuntimeError(f"{uniprot_data} does not exist")
        elif isinstance(uniprot_data, list):
            if False in [isinstance(p, str) for p in uniprot_data] or\
                False in [os.path.isfile(p) for p in uniprot_data]:
                raise RuntimeError('Provided list must contain strings refering \
                                    to existing files')
            else:
                files = uniprot_data
        else:
            raise RuntimeError('uniprot_data must either be string or list of \
                                strings referring to existing files')
        # iterate and return entry by entry
        current_entry = list()
        current_entry_size = 0
        for row in self._uniprot_file_iterator(files):
            # if current_entry_size > max_entry_size:
            #     raise RuntimeError('Observed Entry too large to process')
            if row.startswith('//'):
                if len(current_entry) > 0:
                    yield current_entry
                current_entry = list()
                current_entry_size = 0
            else:
                current_entry.append(row)
                current_entry_size += 1
        if len(current_entry) > 0:
            yield current_entry

    def _uniprot_file_iterator(self, files):
        for f in files:
            if f.endswith('.gz'):
                with gzip.open(f, 'rb') as fh:
                    with io.TextIOWrapper(fh, encoding='utf-8') as decoder:
                        for row in decoder:
                            yield row.rstrip() 
            else:
                with open(f) as fh:
                    for row in fh:
                        yield row.rstrip()
    

class cc_extractor(data_extractor_base):
    """Implements interface defined in :class:`data_extractor_base` and extracts
    regions with coiled coil annotations.
    """

    def id(self):
        """Implements functionality defined in :func:`data_extractor_base.id`
        
        :returns: 'CC'
        :rtype: :class:`str`
        """
        return 'CC'
    def extract(self, entry):
        """Implements functionality defined in :func:`data_extractor_base.extract`

        :returns: Annotated coiled coil regions, None if no entries found
        :rtype: :class:`list` of :class:`str`
        """
        cc_intervals = list()
        for line in entry:
                if line.startswith('FT   COILED'):
                    if '?' not in line:
                        line = line.replace(':',' ')
                        interval = line.split()[-1].split('..')
                        interval = [int(i.replace('<','').replace('>', '')) for i in interval]
                        if len(interval) > 1:
                            cc_intervals.append(interval)

        if len(cc_intervals) > 0:
            return tuple(cc_intervals)
        else:
            return None

class tm_extractor(data_extractor_base):
    """Implements interface defined in :class:`data_extractor_base` and extracts
    regions with transmembrane annotations.
    """

    def id(self):
        """Implements functionality defined in :func:`data_extractor_base.id`
        
        :returns: 'TM'
        :rtype: :class:`str`
        """
        return 'TM'
    
    def extract(self, entry):
        """Implements functionality defined in :func:`data_extractor_base.extract`

        :returns: Annotated coiled coil regions, None if no entries found
        :rtype: :class:`list` of :class:`str`
        """
        tm_intervals = list()
        for line in entry:
                if line.startswith('FT   TRANSMEM'):
                    if '?' not in line:
                        line = line.replace(':',' ')
                        interval = line.split()[-1].split('..')
                        interval = [int(i.replace('<','').replace('>', '')) for i in interval]
                        tm_intervals.append(interval)
                    else:
                        tm_intervals.append('UNK')

        if len(tm_intervals) > 0:
            return tuple(tm_intervals)
        else:
            return None

class sp_extractor(data_extractor_base):
    """Implements interface defined in :class:`data_extractor_base` and extracts
    regions with transmembrane annotations.
    """

    def id(self):
        """Implements functionality defined in :func:`data_extractor_base.id`
        
        :returns: 'SP'
        :rtype: :class:`str`
        """
        return 'SP'
    
    def extract(self, entry):
        """Implements functionality defined in :func:`data_extractor_base.extract`

        :returns: Annotated coiled coil regions, None if no entries found
        :rtype: :class:`list` of :class:`str`
        """
        sp_intervals = list()
        for line in entry:
                if line.startswith('FT   SIGNAL'):
                    if '?' not in line:
                        line = line.replace(':',' ')
                        interval = line.split()[-1].split('..')
                        interval = [int(i.replace('<','').replace('>', '')) for i in interval]
                        sp_intervals.append(interval)
                    else:
                        sp_intervals.append('UNK')

        if len(sp_intervals) > 0:
            return tuple(sp_intervals)
        else:
            return None
        
class disorder_extractor(data_extractor_base):
    """Implements interface defined in :class:`data_extractor_base` and extracts
    regions annotated as instrinsically disordered.
    """

    def id(self):
        """Implements functionality defined in :func:`data_extractor_base.id`
        
        :returns: 'IDP'
        :rtype: :class:`str`
        """
        return 'IDP'

    def extract(self, entry):
        """Implements functionality defined in :func:`data_extractor_base.extract`

        :returns: Annotated intrinsically disordered regions, None if no entries found
        :rtype: :class:`list` of :class:`str`
        """
        idp_intervals = list()
        for line in entry:
                if line.startswith('FT   REGION'):
                    interval = line.split()[2].split('..')
                elif line.startswith('FT') and 'Disordered' in line:
                    interval = [int(i.replace('<','').replace('>', '')) for i in interval]
                    if len(interval) > 1:
                        idp_intervals.append(interval)

        if len(idp_intervals) > 0:
            return tuple(idp_intervals)
        else:
            return None


class ac_extractor(data_extractor_base):
    """Implements interface defined in :class:`data_extractor_base` and extracts
    uniprot accession code (AC).
    """
    def id(self):
        """Implements functionality defined in :func:`data_extractor_base.id`
        
        :returns: 'AC'
        :rtype: :class:`str`
        """
        return 'AC'

    def extract(self, entry):
        """Implements functionality defined in :func:`data_extractor_base.extract`

        :returns: The uniprot accession code, None if not found
        :rtype: :class:`str`
        """
        # the following is a regex pattern which identifies uniprot accession 
        # codes stolen from https://www.uniprot.org/help/accession_numbers
        ac_pat = '[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}'
        for line in entry:
            if line.startswith('AC   '):
                ac = line.split()[1].strip(';')
                if(re.match(ac_pat, ac)):
                    return ac
        return None

class seqlen_extractor(data_extractor_base):
    """Implements interface defined in :class:`data_extractor_base` and extracts
    the canonical sequence of an entry.
    """
    def id(self):
        """Implements functionality defined in :func:`data_extractor_base.id`
        
        :returns: 'LEN'
        :rtype: :class:`str`
        """
        return 'LEN'

    def extract(self, entry):
        """Implements functionality defined in :func:`data_extractor_base.extract`

        :returns: The canonical sequence, None if not found
        :rtype: :class:`str`
        """
        n_AA = None
        in_sq = False
        sq_found = False
        canonical_seq = None
        for line in entry:
            if in_sq:
                if line[:2] != '  ':
                    in_seq = False # done reading sequence
                else:
                    canonical_seq += line.replace(' ', '')
            if line.startswith('SQ   '):
                if sq_found:
                    raise RuntimeError('Expect one canonical seq per entry')
                canonical_seq = ''
                in_sq = True
                sq_found = True
                n_AA = int(line.split()[2])
        if canonical_seq and n_AA != len(canonical_seq):
            raise RuntimeError('Invalid seq length observed')
        return n_AA

class seq_extractor(data_extractor_base):
    """Implements interface defined in :class:`data_extractor_base` and extracts
    the canonical sequence of an entry.
    """
    def id(self):
        """Implements functionality defined in :func:`data_extractor_base.id`
        
        :returns: 'LEN'
        :rtype: :class:`str`
        """
        return 'SEQ'

    def extract(self, entry):
        """Implements functionality defined in :func:`data_extractor_base.extract`

        :returns: The canonical sequence, None if not found
        :rtype: :class:`str`
        """
        n_AA = None
        in_sq = False
        sq_found = False
        canonical_seq = None
        for line in entry:
            if in_sq:
                if line[:2] != '  ':
                    in_seq = False # done reading sequence
                else:
                    canonical_seq += line.replace(' ', '')
            if line.startswith('SQ   '):
                if sq_found:
                    raise RuntimeError('Expect one canonical seq per entry')
                canonical_seq = ''
                in_sq = True
                sq_found = True
                n_AA = int(line.split()[2])
        if canonical_seq and n_AA != len(canonical_seq):
            raise RuntimeError('Invalid seq length observed')
        return canonical_seq

class taxid_extractor(data_extractor_base):
    """Implements interface defined in :class:`data_extractor_base` and extracts
    the canonical sequence of an entry.
    """
    def id(self):
        """Implements functionality defined in :func:`data_extractor_base.id`
        
        :returns: 'LEN'
        :rtype: :class:`str`
        """
        return 'TAXID'

    def extract(self, entry):
        """Implements functionality defined in :func:`data_extractor_base.extract`

        :returns: The canonical sequence, None if not found
        :rtype: :class:`str`
        """
        species_name = ''
        taxid = None
        taxa = []
        for line in entry:
            if line.startswith('OX   '):
                taxid = int(line.split()[1].split('_TaxID=')[-1].strip(';'))
            elif line.startswith('OS   '):
                species_name += ' '.join(line.split()[1:]).strip('.')
            elif line.startswith('OC   '):
                taxa += line.strip('OC   ').strip('.').split(';')

        species_name = species_name.replace(',', '.')

        # print([taxid, species_name, taxa])
        return [taxid, species_name, taxa]

class family_extractor(data_extractor_base):
    """Implements interface defined in :class:`data_extractor_base` and extracts
    the canonical sequence of an entry.
    """
    def id(self):
        """Implements functionality defined in :func:`data_extractor_base.id`
        
        :returns: 'LEN'
        :rtype: :class:`str`
        """
        return 'CHAINS'

    def extract(self, entry):
        """Implements functionality defined in :func:`data_extractor_base.extract`

        :returns: The canonical sequence, None if not found
        :rtype: :class:`str`
        """
        families_intervals = []
        found_family = False
        for line in entry:
            if line.startswith('FT   CHAIN') and '..' in line:
                found_family = True
                interval = line.split()[-1].split('..')
                
                if interval[0] == '?':
                    interval[0] = '1'
                if interval[1] == '?':
                    interval[1] = str(seqlen_extractor().extract(entry))
                    
                interval = [int(i.replace('<','').replace('>', '').replace('?', '')) for i in interval]
                title = ''
                
            elif found_family:
                if '/note=' in line:
                    title += line.split('"')[1]
                    if line.strip().endswith('"'):
                        found_family = False
                        families_intervals.append([title, interval])
                else:
                    title += line.split()[-1].replace('"','')
                    found_family = False
                    families_intervals.append([title, interval])
        
        if len(families_intervals) > 0:
            return families_intervals
        else:
            return None
                                               
class name_extractor(data_extractor_base):
    
    """Implements interface defined in :class:`data_extractor_base` and extracts
    the name of an entry.
    """
    def id(self):
        """Implements functionality defined in :func:`data_extractor_base.id`
        
        :returns: 'NAME'
        :rtype: :class:`str`
        """
        return 'NAME'

    def extract(self, entry):
        """Implements functionality defined in :func:`data_extractor_base.extract`

        :returns: The canonical sequence, None if not found
        :rtype: :class:`str`
        """
        name = None
        found_name = False
        is_fragment = False
        source = 'Reviewed'
        for line in entry:
            if line.startswith('DE   RecName') and not found_name:
                found_name = True
                name = line.strip().split('Full=')[1].split(' {')[0]
                if '|' in line:
                    source = line.strip().split('|')[-1].split(':')[0]
            elif line.startswith('DE   SubName') and not found_name:
                found_name = True
                name = line.strip().split('Full=')[1].split(' {')[0]
                if '|' in line:
                    source = line.strip().split('|')[-1].split(':')[0]
            
            if line.startswith('DE   Flags'):
                is_fragment = True
        
        return {'TITLE': name, 'SOURCE': source, 'FRAG_STATUS': is_fragment}

class evidence_extractor(data_extractor_base):
    """Implements interface defined in :class:`data_extractor_base` and extracts
    the canonical sequence of an entry.
    """
    def id(self):
        """Implements functionality defined in :func:`data_extractor_base.id`
        
        :returns: 'EVIDENCE'
        :rtype: :class:`str`
        """
        return 'EVIDENCE'

    def extract(self, entry):
        """Implements functionality defined in :func:`data_extractor_base.extract`

        :returns: The canonical sequence, None if not found
        :rtype: :class:`str`
        """
        
        evidence, evidence_level = 0, 'Undefined'
        
        for line in entry:
            if line.startswith('PE   '):
                evidence = line.strip(';').split(': ')[1]
                evidence_level = int(line.split()[1].replace(':',''))
        
        return {'LEVEL': evidence_level, 'LABEL': evidence}
    

class coverages_processor(data_extractor_base):
    """Implements interface defined in :class:`data_extractor_base` and computes the coverage with annotations
     of an entry.
    """
    def id(self):
        """Implements functionality defined in :func:`data_extractor_base.id`
        
        :returns: 'LEN'
        :rtype: :class:`str`
        """
        return 'ANNOTCOV'

    def extract(self, entry_uniprot, uniprot_ac, interpro_db):
        """Implements functionality defined in :func:`data_extractor_base.extract`

        :returns: The canonical sequence, None if not found
        :rtype: :class:`str`
        """        
        self.cc_cov = np.nan
        self.dis_cov = np.nan
        self.domains_cov = np.nan
        self.families_cov = np.nan
        all_intervals = []
        self.full_coverage = 0

        curr_len = entry_uniprot['LEN']
        if 'CC' in entry_uniprot:
            self.cc_cov, curr_intervals = self._compute_coverage(entry_uniprot['CC'], curr_len)
            all_intervals += curr_intervals
        if 'IDP' in entry_uniprot:
            self.dis_cov, curr_intervals = self._compute_coverage(entry_uniprot['IDP'], curr_len)
            all_intervals += curr_intervals
        if 'CHAINS' in entry_uniprot:
            try:
                self.families_cov, curr_intervals = self._compute_coverage(entry_uniprot['CHAINS'], curr_len, exclude_duf = True)
                all_intervals += curr_intervals
            except IndexError:
                pass
        try:
            entry_interpro = interpro_db.query(uniprot_ac)[0]['data']
            found_in_interpro = True
            self.domains_cov, curr_intervals = self._compute_coverage(entry_interpro, curr_len)
            self.domains_cov_noDuff, curr_intervals = self._compute_coverage(entry_interpro, curr_len, exclude_duf = True)
            all_intervals += curr_intervals
        except IndexError:
            pass

        if len(all_intervals) > 0:
            all_intervals = self._update_intervals(all_intervals, [], mode = 'multiple')
            self.full_coverage, all_intervals =self._compute_coverage(all_intervals, curr_len) # all intervals already excludes dufs
        
        entry_uniprot = self.register(entry_uniprot)
        
        return entry_uniprot
    
    def register(self, entry_uniprot):
        entry_uniprot[self.id()] = {'CC': self.cc_cov, 'IDP': self.dis_cov, 'DOMAINS': self.domains_cov, 'DOMAINS_noDUF': self.domains_cov_noDuff, 'FAMILIES': self.families_cov, 'FULL_noDUF': self.full_coverage}        
        return entry_uniprot
    
    def _compute_coverage(self, intervals, length, exclude_duf = False, ignore = ['Putative','DUF','Uncharacter','pothetical', 'nknown']):
        
        intervals = [anno for anno in intervals if len(anno)>1]
        try:
            if exclude_duf:
                intervals = [anno for anno in intervals if not any(ext in anno[0] for ext in ignore)]
            intervals = [anno[i] for anno in intervals for i in range(len(anno)) if not isinstance(anno[i][0], str)]
            return round(sum([anno[1]-anno[0] for anno in intervals])*100/length, 2), intervals
        except:
            return round(sum([anno[1]-anno[0] for anno in intervals])*100/length, 2), intervals
                
    def _update_intervals(self, interval, data, mode = 'single'):

        interval = [anno for anno in interval if len(anno)>1]
        interval = np.array(interval)
        interval.sort(axis=0)

        new_data = []
        for curr_interval in interval:
            curr_interval = list(curr_interval)
            if len(curr_interval) > 1:
                if len(new_data) == 0:
                    new_data.append(curr_interval)
                else:
                    previous_interval = list(new_data[-1])
                    if curr_interval[0] <= previous_interval[1]:
                        overlap_range = curr_interval+previous_interval
                        new_data[-1] = [min(overlap_range), max(overlap_range)]
                    else:
                        new_data.append(curr_interval)

        return new_data
    
