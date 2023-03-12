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

# UniParc class and extractors

class data_extractor_base:
	"""Base class defining the interface for data extraction from uniparc 
	entries in xml format. Children must implement :func:`id` and 
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

		:param entry:  Content of uniparc entry in xml format
		:type entry:  :class:`list` of :class:`str`
		"""
		raise NotImplementedError('extract not implemented')

class uniparc_extractor:
	"""This class is capable of iterating over full uniparc releases and 
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
			self.col = self.db['UniParc']
			self.type = 'Mongo'

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
			print(' ... ... Cleaning DB')
			self.col.drop()
			self.col = self.db['UniParc']
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

	def extract(self, uniparc_data, add_if_empty=True, clear = True, max_size = 'all', targets = 'all', chunk_size = 1000, print_step = 10000, saveto = None, savestep = 100000, min_len = 0, simplify_anno = False, exclude_duf = False, ignore = ['Putative','DUF','Uncharacter','Putative', 'nknown']):
		"""Process all entries in file(s) specified in *uniparc_data* with 
		registered data extractors. 

		>>> All entries with a UniProtKB ID are ignored!
		
		The data stored for a single entry is a dict with the return value of
		:func:`data_extractor_base.id` as keys and the return value of
		:func:`data_extractor_base.extract` as values. None values are not 
		added. Once data extraction for an entry is complete, storing 
		happens with :func:`store`.

		:param uniparc_data:  File(s) containing uniparc entries in xml format. 
							 Files ending with .gz are dealt with accordingly.
							 If a list of files is provided, the content of 
							 the files is concatenated when processing.
		:param add_if_empty: None values returned by the data extractors are not 
							 added to the per-entry data dicts. Thus, they might
							 be empty. If set to False, these dicts are ignored,
							 i.e. :func:`store` is not called for those.
		:param max_size:	 Maximum number of elements to extract. Default: all
		:param targets:	  Target accession codes to extract. Default: all

		:type uniparc_data: :class:`str` / :class:`list` of :class:`str`
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

		for entry in self._uniparc_entry_iterator(uniparc_data):
			self.item_count += 1

			uniparc_ac = self.ac_extractor.extract(entry)
			if uniparc_ac is not None:
				# raise RuntimeError('Observed None uniparc AC')
			
				if targets == 'all' or uniparc_ac in targets:
					entry_data = dict()
					for extractor in self.data_extractors:
						if extractor.id() == 'ANNO':
							data = extractor.extract(entry, simplify = simplify_anno, exclude_duf = exclude_duf, ignore = ignore)
						else:
							data = extractor.extract(entry)

						if data is not None:
							entry_data[extractor.id()] = data

					if len(entry_data) > 0 or add_if_empty:
						if 'LEN' not in entry_data or entry_data['LEN'] >= min_len:
							entry_data = self.coverages_processor.extract(entry_data)
							self.store(uniparc_ac, entry_data)

					if self.n_entries % print_step == 0 and self.n_entries > 0:
						if targets == 'all':
							print('UNIPARC:', self.item_count, self.n_entries, 'RSS memory used (GB):', round(psutil.Process(pid).memory_info().rss/1024/1024/1024, 3))
						else:
							print('UNIPARC:', self.item_count, self.n_entries, 'out of {} targets'.format(len(targets)), 'RSS memory used (GB):', round(psutil.Process(pid).memory_info().rss/1024/1024/1024, 3))
					
				if max_size != 'all' and self.item_count >= max_size:
					break
				elif targets != 'all' and self.n_entries == len(targets):
					break

				if saveto is not None and self.item_count % savestep == 0:
					self.save(saveto, self.item_count)

		if saveto is not None:
			self.save(saveto, self.item_count)

		if self.type == 'Mongo':
			self.store(uniparc_ac, data, end = True)

		print('UNIPARC:', self.item_count, self.n_entries, 'RSS memory used (GB):', round(psutil.Process(pid).memory_info().rss/1024/1024/1024, 3))

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
			for uniparc_ac in self.entries:
				outf.write('{}\t{}\n'.format(uniparc_ac, self.saved_index))

		outf.close()
	

	def _uniparc_entry_iterator(self, uniparc_ac, max_entry_size=1048576):
		# check input, uniparc_ac must either be string or list of strings 
		# referring to existing files
		files = None
		if isinstance(uniparc_ac, str):
			if os.path.isfile(uniparc_ac):
				files = [uniparc_ac]
			else:
				raise RuntimeError(f"{uniparc_ac} does not exist")
		elif isinstance(uniparc_ac, list):
			if False in [isinstance(p, str) for p in uniparc_ac] or\
			   False in [os.path.isfile(p) for p in uniparc_ac]:
				raise RuntimeError('Provided list must contain strings refering \
									to existing files')
			else:
				files = uniparc_ac
		else:
			raise RuntimeError('uniparc_data must either be string or list of \
								strings referring to existing files')
		# iterate and return entry by entry
		current_entry = list()
		current_entry_size = 0
		for row in self._uniparc_file_iterator(files):
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

	def _uniparc_file_iterator(self, files):
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
		for line in entry:
			if '<accession>U' in line:
				ac = line.split('>')[1].split('<')[0]
			elif 'type="UniProtKB' in line:
				return None
			elif '<sequence' in line:
				return ac
		return None

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
			if 'NCBI_taxonomy_id' in line:
				taxid = int(line.split('=')[-1].split('"')[1])

		# print([taxid, species_name, taxa])
		return [taxid, species_name, taxa]

class seqlen_extractor(data_extractor_base):
	"""Implements interface defined in :class:`data_extractor_base` and extracts
	the canonical sequence len of an entry.
	"""
	def id(self):
		"""Implements functionality defined in :func:`data_extractor_base.id`
		
		:returns: 'LEN'
		:rtype: :class:`str`
		"""
		return 'LEN'

	def extract(self, entry):
		"""Implements functionality defined in :func:`data_extractor_base.extract`

		:returns: The canonical sequence len, None if not found
		:rtype: :class:`str`
		"""
		n_AA = None
		in_sq = False
		sq_found = False
		canonical_seq = None
		for line in entry:
			if '<sequence' in line:
				if sq_found:
					raise RuntimeError('Expect one canonical seq per entry')
				canonical_seq = line.split('>')[1].split('<')[0]
				in_sq = True
				sq_found = True
				n_AA = int(line.split('"')[1].strip('"'))
		if canonical_seq and n_AA != len(canonical_seq):
			raise RuntimeError('Invalid seq length observed')
		return n_AA

class seq_extractor(data_extractor_base):
	"""Implements interface defined in :class:`data_extractor_base` and extracts
	the canonical sequence len of an entry.
	"""
	def id(self):
		"""Implements functionality defined in :func:`data_extractor_base.id`
		
		:returns: 'LEN'
		:rtype: :class:`str`
		"""
		return 'SEQ'

	def extract(self, entry):
		"""Implements functionality defined in :func:`data_extractor_base.extract`

		:returns: The canonical sequence len, None if not found
		:rtype: :class:`str`
		"""
		n_AA = None
		in_sq = False
		sq_found = False
		canonical_seq = None
		for line in entry:
			if '<sequence' in line:
				if sq_found:
					raise RuntimeError('Expect one canonical seq per entry')
				canonical_seq = line.split('>')[1].split('<')[0]
				in_sq = True
				sq_found = True
				n_AA = int(line.split('"')[1].strip('"'))
		if canonical_seq and n_AA != len(canonical_seq):
			raise RuntimeError('Invalid seq length observed')
		return canonical_seq
	
class annotations_extractor(data_extractor_base):
	"""Implements interface defined in :class:`data_extractor_base` and extracts
	the annotations of an entry.
	"""
	def id(self):
		"""Implements functionality defined in :func:`data_extractor_base.id`
		
		:returns: 'ANNO'
		:rtype: :class:`str`
		"""
		return 'ANNO'

	def extract(self, entry, diggest = True, simplify = False, exclude_duf = False, ignore = ['Putative','DUF','Uncharacter','Putative', 'nknown']):
		"""Implements functionality defined in :func:`data_extractor_base.extract`
		
		:param diggest: Boolean statement to merge overlapping annotations. 
						Default: False
		:type diggest: :class:`bool`
		
		:returns: The sequence signatures/matches with corresponding intervals, None if not found. Overlapping annotations are merged if diggest set to True.
		:rtype: :class:`str`
		"""
		annotations = {'Annotation': [], 'Interval': []}
		found_annotation = False
		for line in entry:
			if '<signatureSequenceMatch database="' in line:
				found_annotation = True
			elif '</signatureSequenceMatch>' in line:
				found_annotation = False
			elif found_annotation:
				if '<ipr name=' in line:
					annotation = line.split('=')[1].split('"')[1].strip('"')
				elif '<lcn start=' in line:
					interval = [int(line.split('=')[1].split()[0].strip('"')), int(line.split('=')[2].split('/')[0].strip('"'))]
					
					if diggest:
						annotations = self._update_intervals(interval, annotation, annotations, ignore = ignore)
					else:
						annotations['Annotation'].append(annotation)
						annotations['Interval'].append(interval)
		
		return self._format_annotations(annotations, simplify, exclude_duf, ignore = ignore)
		
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

	def extract(self, entry_uniparc):
		"""Implements functionality defined in :func:`data_extractor_base.extract`

		:returns: The canonical sequence, None if not found
		:rtype: :class:`str`
		"""		
		self.cc_cov = np.nan
		self.dis_cov = np.nan
		self.domains_cov = np.nan
		self.domains_cov_noDuff = np.nan
		self.families_cov = np.nan
		all_intervals = []
		self.full_coverage = 0

		if 'ANNO' in entry_uniparc:
			curr_len = entry_uniparc['LEN']
			self.domains_cov, curr_intervals = self._compute_coverage(entry_uniparc['ANNO'], curr_len)
			self.domains_cov_noDuff, curr_intervals = self._compute_coverage(entry_uniparc['ANNO'], curr_len, exclude_duf = True)
			all_intervals += curr_intervals

		if len(all_intervals) > 0:
			all_intervals = self._update_intervals(all_intervals, [], mode = 'multiple')
			self.full_coverage, all_intervals = self._compute_coverage(all_intervals, curr_len) # all intervals already excludes dufs

		entry_uniparc = self.register(entry_uniparc)
		
		return entry_uniparc
	
	def register(self, entry_uniprot):
		entry_uniprot[self.id()] = {'CC': self.cc_cov, 'IDP': self.dis_cov, 'DOMAINS': self.domains_cov, 'DOMAINS_noDUF': self.domains_cov_noDuff, 'FAMILIES': self.families_cov, 'FULL_noDUF': self.full_coverage}		
		return entry_uniprot
	
	def _compute_coverage(self, intervals, length, exclude_duf = False, ignore = ['Putative','DUF','Uncharacter','pothetical', 'nknown']):

		try:
			if exclude_duf:
				intervals = [anno for anno in intervals if not any(ext in anno[0] for ext in ignore)]
			intervals = [anno[i] for anno in intervals for i in range(len(anno)) if not isinstance(anno[i][0], str)]
			return round(sum([anno[1]-anno[0] for anno in intervals])*100/length, 2), intervals
		except:
			return round(sum([anno[1]-anno[0] for anno in intervals])*100/length, 2), intervals
	
	def _update_intervals(self, interval, data, mode = 'single'):

		interval = np.array(interval)
		interval.sort(axis=0)

		new_data = []
		for curr_interval in interval:
			curr_interval = list(curr_interval)
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