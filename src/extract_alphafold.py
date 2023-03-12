import os
import json
import sys
import time
import gzip
import tarfile
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

class alphafold_extractor:
    """This class is capable if iterating over full uniprot releases and 
    extract data using :class:`data_extractor_base` objects which can
    be registered with :class:`register`. Data extraction happens upon
    calling :func:`extract`.
    """
    def __init__(self, mongo_host = None, mongo_port = None):

        if mongo_host is None or mongo_port is None:
            self.data = list()
            self.entries = list()
            self.type = 'List'
        else:
            self.client = pymongo.MongoClient(mongo_host, mongo_port)
            self.db = self.client['Joana'] 
            self.col = self.db['AlphaFold']
            self.type = 'Mongo'

        self.data_extractors = list()
        self.n_entries = 0
        self.item_count = 0

    def clear(self):
        """
        Clean the database. It becomes an empty set of lists.
        """
        if self.type == 'Mongo':
            self.col.drop()
            self.col = self.db['AlphaFold']

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

    def extract(self, alphafold_data, max_size = 'all', add_if_empty=True, clear = True, targets = 'all', chunk_size = 1000, print_step = 100000, saveto = None, savestep = 100000):
        """Process all entries in file(s) specified in *alphafold_data* with 
        registered data extractors.
        The data stored for a single entry is a dict with the return value of
        :func:`data_extractor_base.id` as keys and the return value of
        :func:`data_extractor_base.extract` as values. None values are not 
        added. Once data extraction for an entry is complete, storing 
        happens with :func:`store`.

        :param alphafold_data: Folder containing the alphafold proteomes in tar format
                             separated by //. Files ending with .tar are dealt
                             with accordingly. If a list of files is provided,
                             the content of the files is concatenated when
                             processing.
        :param add_if_empty: None values returned by the data extractors are not 
                             added to the per-entry data dicts. Thus, they might
                             be empty. If set to False, these dicts are ignored,
                             i.e. :func:`store` is not called for those.

        :type alphafold_data: :class:`str` / :class:`list` of :class:`str`
        :type add_if_empty: :class:`bool`
        """

        if clear:
            self.clear()

        self.chunk_size = chunk_size
        self.n_chunks = 0
        self.curr_chunk_size = 0
        self.data_dict = list()

        for uniprot_ac, entry in self._alphafold_entry_iterator(alphafold_data):
            self.item_count += 1
            data = None

            if targets == 'all' or uniprot_ac in targets:
                entry_data = dict()
                for extractor in self.data_extractors:
                    data = extractor.extract(entry)
                    if data is not None:
                        entry_data[extractor.id()] = data
                
                if len(entry_data) > 0 or add_if_empty:
                    self.store(uniprot_ac, entry_data)

                if self.item_count % print_step == 0:
                    if targets == 'all':
                        print('ALPHAFOLD:', self.item_count, self.n_entries, 'RSS memory used (GB):', round(psutil.Process(pid).memory_info().rss/1024/1024/1024, 3))
                    else:
                        print('ALPHAFOLD:', self.n_entries, 'out of {} targets'.format(len(targets)), 'RSS memory used (GB):', round(psutil.Process(pid).memory_info().rss/1024/1024/1024, 3))

            if max_size != 'all' and self.item_count >= max_size:
                break
            elif targets != 'all' and self.n_entries == len(targets):
                break

            if saveto is not None and self.item_count % savestep == 0:
                self.save(saveto, self.item_count)

        if saveto is not None:
            self.save(saveto, self.item_count)

        if self.type == 'Mongo' and data is not None:
            self.store(uniprot_ac, data, end = True)

        print('ALPHAFOLD:', self.item_count, self.n_entries, 'RSS memory used (GB):', round(psutil.Process(pid).memory_info().rss/1024/1024/1024, 3))

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
                self._reorganize_data_dict()
                if len(self.data_dict) > 0:
                    self.col.insert_many(self.data_dict)
                    self.data_dict = list()
                    self.n_chunks += 1
                    self.curr_chunk_size = 0

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

    def _alphafold_entry_iterator(self, db_path):

        db_path = '{}/proteomes'.format(db_path)
        for f in os.scandir(db_path): 
            if f.path.endswith('.tar'):
                tar = tarfile.open(f.path)
                for fh in tar.getmembers():
                    if '.cif' in fh.name :
                        uniprot_ac = '-'.join(fh.name.split('-')[1:3])
                        data = {'coords': fh.name}
                    elif 'confidence' in fh.name and uniprot_ac in fh.name:
                        conf_data = tar.extractfile(fh)
                        conf_data = json.load(gzip.open(conf_data, 'rb'))
                        data['confidence'] = conf_data
                        
                        yield uniprot_ac, data
    
    def _reorganize_data_dict(self):
        
        new_dict = {}
        for entry in self.data_dict:
            _id  = entry['_id']
            data = entry['data']
            
            new_id = _id.split('-')[0]
            frg_id = _id.split('-')[1]
            
            if new_id not in new_dict:
                new_dict[new_id] = {}
            new_dict[new_id][frg_id] = data
        
        self.data_dict = list()
        for ac in new_dict:
            data = new_dict[ac]
            self.data_dict.append({'_id': ac, 'data': data})

class pLDDT_extractor(data_extractor_base):
    """Implements interface defined in :class:`data_extractor_base` and extracts
    the average pLDDT of the model.
    """

    def id(self):
        """Implements functionality defined in :func:`data_extractor_base.id`

        :returns: 'pLDDT'
        :rtype: :class:`str`
        """
        return 'pLDDT'
    
    def extract(self, entry):
        """Implements functionality defined in :func:`data_extractor_base.extract`

        :returns: average pLDDT of the model, None if no entries found
        :rtype: :class:`float`
        """        
        plddt = np.mean(entry['confidence']['confidenceScore'])
        confcat_freq = {}
        for categ in {'M', 'D', 'H', 'L'}:
            confcat_freq[categ] = entry['confidence']['confidenceCategory'].count(categ)*100/len(entry['confidence']['confidenceCategory'])
        
        return {'avg_pLDDT': plddt, 'CategoriesFreq': confcat_freq, 'Lenght': len(entry['confidence']['confidenceScore'])}

    