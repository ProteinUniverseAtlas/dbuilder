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

import urllib.parse
import urllib.request
import subprocess as sp

import numpy as np

pid = os.getpid()

# UNIREF

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

class uniref_extractor:
    """This class is capable of iterating over full uniref releases and 
    extract data using :class:`data_extractor_base` objects which can
    be registered with :class:`register`. Data extraction happens upon
    calling :func:`extract`.
    """

    def __init__(self, mongo_host = None, mongo_port = None, name = 'UniRef50'):

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
        self.darkness_extractor = darkness_extractor()
        self.n_entries = 0
        self.item_count = 0

    def clear(self):
        """
        Clean the database. It becomes an empty set of lists.
        """
        if self.type == 'Mongo':
            print(' ... ... Cleaning DB')
            self.col.drop()
            self.col = self.db[self.name]
            print(' ... ... Done!')

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

    def extract(self, uniref_data, add_if_empty=True, clear = True, max_size = 'all', targets = 'all', chunk_size = 1000, print_step = 100000, saveto = None, savestep = 100000, uniprot_db = None, uniparc_db = None, alphafold_db = None, update_unip = False):
        """Process all entries in file(s) specified in *uniref_data* with 
        registered data extractors.
        The data stored for a single entry is a dict with the return value of
        :func:`data_extractor_base.id` as keys and the return value of
        :func:`data_extractor_base.extract` as values. None values are not 
        added. Once data extraction for an entry is complete, storing 
        happens with :func:`store`.

        :param uniref_data:  File(s) containing uniref entries in xml format. 
                             Files ending with .gz are dealt with accordingly.
                             If a list of files is provided, the content of 
                             the files is concatenated when processing.
        :param add_if_empty: None values returned by the data extractors are not 
                             added to the per-entry data dicts. Thus, they might
                             be empty. If set to False, these dicts are ignored,
                             i.e. :func:`store` is not called for those.
        :param max_size:	 Maximum number of elements to extract. Default: all
        :param targets:	  Target accession codes to extract. Default: all

        :type uniref_data: :class:`str` / :class:`list` of :class:`str`
        :type add_if_empty: :class:`bool`
        :type max_size: :class:`int`
        :type targets: :class:`str` / :class:`list` of :class:`str`
        """

        if clear:
            self.clear()

        self.chunk_size = chunk_size
        self.n_chunks = 0
        self.curr_chunk_size = 0
        self.data_dict = list()

        for entry in self._uniref_entry_iterator(uniref_data):
            self.item_count += 1
            data = None

            uniref_ac = self.ac_extractor.extract(entry)

            # contine only if it is a valid UniRef identifier and it does not exist in the DB already
            if uniref_ac is not None and self.query(uniref_ac).count() == 0: 
                if targets == 'all' or uniref_ac in targets:
                    entry_data = dict()
                    for extractor in self.data_extractors:
                        data = extractor.extract(entry)
                        if data is not None:
                            entry_data[extractor.id()] = data

                    if len(entry_data) > 0 or add_if_empty:	
                        entry_data = self.darkness_extractor.extract(entry_data, uniref_ac, uniprot_db, uniparc_db, alphafold_db)                         
                        if update_unip and uniprot_db is not None and uniparc_db is not None:
                            self.update_uniprot_and_uniparc(uniref_ac, entry_data, uniprot_db, uniparc_db, alphafold_db)
                        
                        self.store(uniref_ac, entry_data)

#                     if self.n_entries % print_step == 0:
#                         if targets == 'all':
#                             print('UNIREF:', self.item_count, self.n_entries, 'RSS memory used (GB):', round(psutil.Process(pid).memory_info().rss/1024/1024/1024, 3))
#                         else:
#                             print('UNIREF:', self.item_count, self.n_entries, 'out of {} targets'.format(len(targets)), 'RSS memory used (GB):', round(psutil.Process(pid).memory_info().rss/1024/1024/1024, 3))
                
            if self.item_count % print_step == 0:
                if targets == 'all':
                    print('UNIREF:', self.item_count, self.n_entries, 'RSS memory used (GB):', round(psutil.Process(pid).memory_info().rss/1024/1024/1024, 3))
                else:
                    print('UNIREF:', self.item_count, self.n_entries, 'out of {} targets'.format(len(targets)), 'RSS memory used (GB):', round(psutil.Process(pid).memory_info().rss/1024/1024/1024, 3))

            if max_size != 'all' and self.item_count >= max_size:
                break
            elif targets != 'all' and self.n_entries == len(targets):
                break

            if saveto is not None and self.item_count % savestep == 0:
                print(uniref_ac)
                self.save(saveto, self.item_count)

        if saveto is not None:
            self.save(saveto, self.item_count)

        if self.type == 'Mongo':
            self.store(uniref_ac, data, end = True)

        print('UNIREF:', self.item_count, self.n_entries, 'RSS memory used (GB):', round(psutil.Process(pid).memory_info().rss/1024/1024/1024, 3))


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
        return self.col.find({ '_id': { "$in": [ac] } })

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
            for uniref_ac in self.entries:
                outf.write('{}\t{}\n'.format(uniref_ac, self.saved_index))

        outf.close()

    def _uniref_entry_iterator(self, uniref_ac, max_entry_size=1048576):
        # check input, uniref_ac must either be string or list of strings 
        # referring to existing files
        files = None
        if isinstance(uniref_ac, str):
            if os.path.isfile(uniref_ac):
                files = [uniref_ac]
            else:
                raise RuntimeError(f"{uniref_ac} does not exist")
        elif isinstance(uniref_ac, list):
            if False in [isinstance(p, str) for p in uniref_ac] or\
               False in [os.path.isfile(p) for p in uniref_ac]:
                raise RuntimeError('Provided list must contain strings refering \
                                    to existing files')
            else:
                files = uniref_ac
        else:
            raise RuntimeError('uniprot_data must either be string or list of \
                                strings referring to existing files')
        # iterate and return entry by entry
        current_entry = list()
        current_entry_size = 0
        for row in self._uniref_file_iterator(files):
            # if current_entry_size > max_entry_size:
            # 	raise RuntimeError('Observed Entry too large to process')
            if row.startswith('</entry>'):
                if len(current_entry) > 0:
                    yield current_entry
                current_entry = list()
                current_entry_size = 0
            else:
                current_entry.append(row)
                current_entry_size += 1
        if len(current_entry) > 0:
            yield current_entry

    def _uniref_file_iterator(self, files):
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

    # functions to update uniprot and uniparc with darkness values
    def update_uniprot_and_uniparc(self, uniref_ac, entry_data, uniprot_db, uniparc_db, alphafold_db):
        
        uniref_data = {'UNIREF_AC': uniref_ac,  
                       'FULL_noDUF': entry_data['DARKNESS']['FULL_noDUF']}
        
        uniprot_task_array = list()
        uniparc_task_array = list()
        alphfol_task_array = list()
        
        for ac in entry_data['ACC']:
            if ac.startswith('UP'):
                uniparc_task_array.append(pymongo.UpdateOne({'_id': ac}, {'$set': {self.name: uniref_data}}))
            else:
                uniprot_task_array.append(pymongo.UpdateOne({'_id': ac}, {'$set': {self.name: uniref_data}}))
                alphfol_task_array.append(pymongo.UpdateOne({'_id': ac}, {'$set': {self.name: uniref_data}}))
                
        if len(uniprot_task_array) > 0 and uniprot_db is not None:
            uniprot_db.col.bulk_write(uniprot_task_array)

        if len(uniparc_task_array) > 0 and uniparc_db is not None:
            uniparc_db.col.bulk_write(uniparc_task_array)
        
        if len(alphfol_task_array) > 0 and alphafold_db is not None:
            alphafold_db.col.bulk_write(alphfol_task_array)
        
        
class ac_extractor(data_extractor_base):
    """Implements interface defined in :class:`data_extractor_base` and extracts
    uniref accession code (AC).
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
        for line in entry:
            if line.startswith('<entry id='):
                ac = line.split('="')[1].split('"')[0]
                return ac
        return None

class entries_extractor(data_extractor_base):
    """Implements interface defined in :class:`data_extractor_base` and extracts
    the uniprotKB and uniparc accession codes of the corresponding uniref cluster.
    """
    def id(self):
        """Implements functionality defined in :func:`data_extractor_base.id`

        :returns: 'ACC'
        :rtype: :class:`str`
        """
        return 'ACC'

    def extract(self, entry):
        """Implements functionality defined in :func:`data_extractor_base.extract`

        :returns: Cluster member entries, None if no entries found
        :rtype: :class:`list` of :class:`str`
        """
        member_entries = list()
        search_uniprotid = False
        for line in entry:
            if line.startswith('<dbReference type='):
                if 'UniProtKB' in line:
                    # member = line.split('id=')[1].split('"')[1].split('_')[0]
                    search_uniprotid = True
                else:
                    member = line.split('id=')[1].split('"')[1].split('"')[0]
                    member_entries.append(member)	

            elif line.startswith('<property type=') and search_uniprotid:
                if 'UniProtKB' in line:
                    member = line.split('value=')[1].split('"')[1]
                    member_entries.append(member)
                    search_uniprotid = False

        if len(member_entries) > 0:
            return list(set(member_entries))
        else:
            return None

class unirefs_extractor(data_extractor_base):
    """Implements interface defined in :class:`data_extractor_base` and extracts
    the other UniRef ids associated with the corresponding uniref cluster.
    """
    def id(self):
        """Implements functionality defined in :func:`data_extractor_base.id`

        :returns: 'UNIREF'
        :rtype: :class:`str`
        """
        return 'UNIREF'

    def extract(self, entry):
        """Implements functionality defined in :func:`data_extractor_base.extract`

        :returns: Cluster UniRef clusters associates, None if no entries found
        :rtype: :class:`list` of :class:`str`
        """
        member_entries = {}
        for line in entry:
            if line.startswith('<property type='):
                if 'UniRef' in line:
                    member = line.split('value=')[1].split('"')[1]
                    if member.split('_')[0] not in member_entries:
                        member_entries[member.split('_')[0]] = []
                        
                    member_entries[member.split('_')[0]].append(member)

        for uniref_level in member_entries:
            member_entries[uniref_level] = list(set(member_entries[uniref_level]))
            
        if len(member_entries) > 0:
            return member_entries
        else:
            return None

class darkness_extractor(data_extractor_base):
    """Implements interface defined in :class:`data_extractor_base` and computes the coverage with annotations
     of an entry.
    """
    def id(self):
        """Implements functionality defined in :func:`data_extractor_base.id`

        :returns: 'LEN'
        :rtype: :class:`str`
        """
        return 'DARKNESS'

    def extract(self, entry_uniref, uniprot_ac, uniprot_db, uniparc_db, alphafold_db):
        """Implements functionality defined in :func:`data_extractor_base.extract`

        :returns: The canonical sequence, None if not found
        :rtype: :class:`str`
        """		
        uniprot_acs = entry_uniref['ACC']

        self.uniprot_data   = uniprot_db.col.find({ '_id': { "$in": uniprot_acs }})
        self.uniparc_data   = uniparc_db.col.find({ '_id': { "$in": uniprot_acs }})  
        self._select_representative()

        if alphafold_db is not None:
            self.alphafold_data = alphafold_db.col.find({ '_id': { "$in": uniprot_acs }})
            self._add_alphafold_confidences()

        entry_uniref = self.register(entry_uniref)

        return entry_uniref

    def _select_representative(self, parameter='FULL_noDUF'):
        self.full_coverage      = 0
        self.representative     = None
        self.length_rep         = 0
        self.is_transmembrane   = False
        self.has_signalpeptide  = False

        for db_data in [self.uniprot_data, self.uniparc_data]:
            for document in db_data:
                ac = document['_id']
                curr_cov = document['data']['ANNOTCOV'][parameter]                
                if curr_cov > self.full_coverage:
                    self.full_coverage = curr_cov
                    self.representative  = ac
                if 'TM' in document['data'] and 'SP' in document['data']:
                    if document['data']['TM'] is not None:
                        self.is_transmembrane = True
                    if document['data']['SP'] is not None:
                        self.has_signalpeptide = True
                        

    def _add_alphafold_confidences(self):
        self.pLDDTs = []
        self.best_af2  = None
        self.worst_af2 = None
            
        for document in self.alphafold_data:
            ac   = document['_id']
            data = document['data']
            
            avgPLDDT = []
            n_res = 0
            for fragment in data:
                avgPLDDT.append(data[fragment]['pLDDT']['avg_pLDDT']*data[fragment]['pLDDT']['Lenght'])
                n_res += data[fragment]['pLDDT']['Lenght']
            
            fullprotein_pLDDT = sum(avgPLDDT)/n_res            
            if len(self.pLDDTs) == 0 or fullprotein_pLDDT > max(self.pLDDTs):
                self.best_af2 = {'ACC': ac, 'LEN': n_res}
            
            if len(self.pLDDTs) == 0 or fullprotein_pLDDT < min(self.pLDDTs):
                self.worst_af2 = {'ACC': ac, 'LEN': n_res}
            
            self.pLDDTs.append(fullprotein_pLDDT)
        
        if self.best_af2 == self.worst_af2:
            self.worst_af2 = None

#     def update_db(self, entry_uniref, db, db_data):

#         for document in db_data:
#             ac   = document['_id']
#             data = document['data']
#             data[self.id()] = {'UNIREF_AC': entry_uniref, 'REP': self.representative, 'FULL_noDUF': self.full_coverage}
#             db.col.update({'_id': ac}, {'$set': {'data': data}})

    def register(self, entry_uniref):
        entry_uniref[self.id()] = {'REP': self.representative, 
                                'FULL_noDUF': self.full_coverage, 
                                'TM': self.is_transmembrane,
                                'SP': self.has_signalpeptide}
        try:
            entry_uniref[self.id()]['pLDDTs']  = self.pLDDTs
            entry_uniref[self.id()]['AF2_REP_best'] = self.best_af2
            entry_uniref[self.id()]['AF2_REP_worst'] = self.worst_af2
        except:
            pass
        
        return entry_uniref
	
	

