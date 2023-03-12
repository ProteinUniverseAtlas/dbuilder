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

# Interpro class

class interpro_db_diggested:
	"""Sets up a dictionary database, that can be filled from a text file (the protein2ipr.dat.gz )
	available from Interpro, or a prefilled json file.
	"""
	def __init__(self, mongo_host = None, mongo_port = None):

		if mongo_host is None or mongo_port is None:
			self.data = list()
			self.entries = list()
			self.type = 'List'
		else:
			self.client = pymongo.MongoClient(mongo_host, mongo_port)
			self.db = self.client['Joana'] 
			self.col = self.db['Interpro']
			self.type = 'Mongo'

		self.n_entries = 0
		self.item_count = 0

	def clear(self):
		"""
		Clean the database. It becomes an empty set of lists.
		"""
		if self.type == 'Mongo':
			self.col.drop()
			self.col = self.db['Interpro']

		elif self.type == 'List':
			self.data = list()
			self.entries = list()
	
	def save_to_pickle(self, out_file):
		"""
		Save the database to pickle
		
		:param out_file: The file where to save.
		:type out_file: :class:`str`
		"""
		pickle.dump(self, open(out_file, 'wb'))

		
	def fill_from_file(self, file_path, max_size = 'all', print_step = 100000, targets = 'all', saveto = None, clear = True, savestep = 100000, chunk_size = 1000, simplify = False, exclude_duf = False, ignore = ['Putative','DUF','Uncharacter','Putative', 'nknown']):
		"""Clears content from interpro collection in db and refills it with the
		content from *file_path*.

		:param file_path:   Path to file containing InterPro data in a 
							tab-delimited format. Can be downloaded from 
							<https://www.ebi.ac.uk/interpro/download/>_ and
							represents the mapping from UniProtKB proteins 
							to InterPro entries. At the time of writing this
							was called protein2ipr.dat.gz.
		:param max_size:	Maximum number of elements in the database. Useful
							for testing, when we want to include a maximum of 
							entries in the database. Default: None					   
		:param print_step:  Step count to print progress. Default: 100000 

		:type file_path: :class:`str`
		:type max_size: :class:`str` 
		:type print_step: :class:`str` 
		"""

		# clean the mongo collection
		if clear:
			self.clear()

		added_domains = 0
		previous_entry = None
		self.chunk_size = chunk_size
		self.n_chunks = 0
		self.curr_chunk_size = 0
		self.data_dict = list()

		n_targets = len(targets)
		for item in self._file_iterator(file_path):
			item = item.split('\t')
			uniprot_ac = item[0]
			if uniprot_ac != previous_entry:
				if previous_entry is not None:
					self.item_count += 1
					if previous_entry in targets or targets == 'all':

						id_data = self._format_annotations(id_data, simplify, exclude_duf, ignore = ignore)
						self.store(previous_entry, id_data)
						added_domains += len(id_data)
											
						if self.item_count%print_step == 0:
							if targets == 'all':
								print('INTERPRO:', self.item_count, self.n_entries, 'RSS memory used (GB):', round(psutil.Process(pid).memory_info().rss/1024/1024/1024, 3))
							else:
								print('INTERPRO:', self.item_count, self.n_entries, 'out of {} targets'.format(n_targets), 'RSS memory used (GB):', round(psutil.Process(pid).memory_info().rss/1024/1024/1024, 3))
							
					if max_size != 'all' and self.item_count >= max_size :
						break
					elif targets != 'all' and self.n_entries == n_targets:
						break

					if self.type != 'Mongo' and saveto is not None and self.item_count % savestep == 0:
						self.save(saveto, self.item_count)

				id_data = {'Annotation': [item[2]], 'Interval': [[int(item[4]), int(item[5])]]}
			else:
				curr_interval = [int(item[4]), int(item[5])]
				curr_annotation = item[2]
				
				id_data = self._update_intervals(curr_interval, curr_annotation, id_data, ignore = ignore)
				
			previous_entry = uniprot_ac

		if self.type != 'Mongo' and saveto is not None:
			self.save(saveto, self.item_count)

		if self.type == 'Mongo':
			self.store(uniprot_ac, id_data, end = True)

		print('INTERPRO:', self.item_count, self.n_entries, 'RSS memory used (GB):', round(psutil.Process(pid).memory_info().rss/1024/1024/1024, 3))

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
		
	def _file_iterator(self, f):
		if f.endswith('.gz'):
			with gzip.open(f, 'rb') as fh:
				with io.TextIOWrapper(fh, encoding='utf-8') as decoder:
					for row in decoder:
						yield row.rstrip() 
		else:
			with open(f) as fh:
				for row in fh:
					yield row.rstrip()
   
	def _update_intervals(self, interval, annotation, data, ignore):

		new_data = data.copy()
		overlaps_with = []
		
		for i in range(len(data['Interval'])):
			curr_interval = data['Interval'][i]
			overlap_range = curr_interval+interval
			if (curr_interval[1]-curr_interval[0])+(interval[1]-interval[0]) > max(overlap_range)-min(overlap_range):
				# they overlap
				overlaps_with.append(i)
				
				if curr_interval[1]-curr_interval[0] > interval[1]-interval[0]:
					if not any(ext in data['Annotation'][i] for ext in ignore):
					   annotation = data['Annotation'][i]
					
				interval = [min(overlap_range), max(overlap_range)]
			
		new_data = {key: [data[key][i] for i in range(len(data[key])) if i not in overlaps_with] for key in data}
		new_data['Annotation'].append(annotation)
		new_data['Interval'].append(interval)

		return new_data
	
	def _format_annotations(self, data, simplify, exclude_duf, ignore):

		if len(data['Annotation']) > 0:
			annotations = []
			for annotation in set(data['Annotation']):
				if not exclude_duf or not any(ext in annotation for ext in ignore):
					if simplify:
						curr_lst = tuple(tuple(data['Interval'][i]) for i in range(len(data['Annotation'])) if data['Annotation'][i] == annotation)
					else:
						curr_lst = tuple([annotation] + [tuple(data['Interval'][i]) for i in range(len(data['Annotation'])) if data['Annotation'][i] == annotation])
					annotations.append(curr_lst)
			return tuple(annotations)
		else:
			return None