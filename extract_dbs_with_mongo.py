import os
import json
import sys
import time
import gzip
import io
import re
import pickle
import pymongo
import numpy as np

try:
    import psutil
except:
    pass

from datetime import date

pid = os.getpid()
start_memory = psutil.Process(pid).memory_info().rss
print('\nSTART RSS memory use (GB): {}'.format(round(start_memory/1024/1024/1024,3)))

# import my classes from src
from src import extract_interpro   as interpro
from src import extract_uniparc    as uniparc
from src import extract_uniprot    as uniprot
from src import extract_uniref     as uniref
from src import extract_alphafold  as alphafold

### RUN FIRST: DOWNLOAD THE FILES BELOW AND SAVE TO A FOLDER NAMED 'databases'

# INPUTS:

uniprot_files  = ['databases/uniprot_sprot_2022-08-03.dat.gz', 'databases/uniprot_trembl_2022-08-03.dat.gz']
uniparc_file   =  'databases/uniparc_all_2022-08-03.xml.gz'
uniref50_file    =  'databases/uniref50_2022-08-03.xml.gz'
uniref90_file    =  'databases/uniref90_2022-08-03.xml.gz'
interpro_file  =  'databases/protein2ipr_2022-08-03.dat.gz'
alphafold_path =  # path to AFDBv4

uniprot_files_updated_names  = ['databases/uniprot_sprot_2022-12-14.dat.gz', 'databases/uniprot_trembl_2022-12-14.dat.gz']

n_to_extract = 'all'
chunk_size = 100000

MONGO_HOST = None # insert your previously defined mongo host, e.g."10.1.0.202"
MONGO_PORT = None # insert your previously defined mongo port, e.g. 30077

# Define automatic logger

class Logger(object):
	def __init__(self, name):
		self.terminal = sys.stdout
		self.log = open(name, "w")

	def write(self, message):
		self.terminal.write(message)
		self.log.write(message)  

	def flush(self):
		#this flush method is needed for python 3 compatibility.
		#this handles the flush command by doing nothing.
		#you might want to specify some extra behavior here.
		pass	

logfile = 'databases_extraction_{}_{}.log'.format(n_to_extract, date)
sys.stdout = Logger(logfile)


# START EXTRACTION:

# 1. First extract interpro 

print('\n1. TAKING CARE OF INTERPRO: {}\n'.format(interpro_file))

start = time.time()
day0  = date.today()

interpro_db = interpro.interpro_db_diggested(mongo_host = MONGO_HOST, mongo_port = MONGO_PORT)

print(' ... Extracting!\n')
start_memory = psutil.Process(pid).memory_info().rss
interpro_db.fill_from_file(interpro_file, print_step = 10000, chunk_size = chunk_size)

numb_seconds = time.time() - start
numb_days    = date.today() - day0
print('\n ... Time to fill interpro_db: {} days {}'.format(numb_days.days, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))
print(' ... Number of entries collected:', interpro_db.n_entries)
print(' ... Total RSS memory use so far (GB):',  round(psutil.Process(pid).memory_info().rss/1024/1024/1024, 3))
print(' ... Total RSS memory use for this step (GB):',  round((psutil.Process(pid).memory_info().rss - start_memory)/1024/1024/1024,3))

print(' ... Indexing!\n')
start = time.time()
interpro_db.index_db()
numb_seconds = time.time() - start
print(' ... Time to index interpro_db: {}'.format(time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))
print(' ... Total RSS memory use so far (GB):',  round(psutil.Process(pid).memory_info().rss/1024/1024/1024, 3))
print(' ... Total RSS memory use for this step (GB):',  round((psutil.Process(pid).memory_info().rss - start_memory)/1024/1024/1024,3))

# 2. Now extract target uniprots from uniprot (sprot and trembl) and add annotations coverages

print('\n2. TAKING CARE OF UNIPROTS: {}\n'.format(uniprot_files))

start = time.time()
day0  = date.today()

uniprot_db = uniprot.uniprot_extractor(mongo_host = MONGO_HOST, mongo_port = MONGO_PORT)
uniprot_db.register(uniprot.seqlen_extractor())
uniprot_db.register(uniprot.disorder_extractor())
uniprot_db.register(uniprot.cc_extractor())
uniprot_db.register(uniprot.tm_extractor())
uniprot_db.register(uniprot.sp_extractor())
uniprot_db.register(uniprot.seq_extractor())
uniprot_db.register(uniprot.taxid_extractor())
uniprot_db.register(uniprot.family_extractor())
uniprot_db.register(uniprot.name_extractor())
uniprot_db.register(uniprot.evidence_extractor())

start_memory = psutil.Process(pid).memory_info().rss

clear = True
for uniprot_file in uniprot_files:

    print(' ... Extracting {}\n'.format(uniprot_file))
    uniprot_db.extract(uniprot_file, clear = clear, print_step = 10000, chunk_size = chunk_size, interpro_db = interpro_db)

    numb_seconds = time.time() - start
    numb_days    = date.today() - day0
    print(' ... ... Time to fill uniprot_db with {}: {} days {}'.format(uniprot_file, numb_days.days, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))
    print(' ... ... Number of entries collected:',uniprot_db.n_entries)
    print(' ... ... Total RSS memory use so far (GB):',  round(psutil.Process(pid).memory_info().rss/1024/1024/1024, 3))
    print(' ... ... Total RSS memory use for this step (GB):',  round((psutil.Process(pid).memory_info().rss - start_memory)/1024/1024/1024, 3))
    clear = False

print(' ... Total RSS memory use so far (GB):', psutil.Process(pid).memory_info().rss/1024/1024/1024)

print(' ... Indexing!\n')
start = time.time()
uniprot_db.index_db()
numb_seconds = time.time() - start
print(' ... Time to index uniprot_db: {}'.format(time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))
print(' ... Total RSS memory use for this step (GB):',  round((psutil.Process(pid).memory_info().rss - start_memory)/1024/1024/1024,3))

# 3. Now extract uniparcs from uniparc and add annotations coverages

print('\n3. TAKING CARE OF UNIPARCS: {}\n'.format(uniparc_file))

# start = time.time()
# day0  = date.today()

uniparc_db = uniparc.uniparc_extractor(mongo_host = MONGO_HOST, mongo_port = MONGO_PORT)
uniparc_db.register(uniparc.annotations_extractor())
uniparc_db.register(uniparc.seqlen_extractor())
uniparc_db.register(uniparc.seq_extractor())
uniparc_db.register(uniparc.taxid_extractor())

print(' ... Extracting!\n')
start_memory = psutil.Process(pid).memory_info().rss
uniparc_db.extract(uniparc_file, print_step = 10000, chunk_size = chunk_size)

numb_seconds = time.time() - start
numb_days    = date.today() - day0
print('\n ... Time to fill uniparc_db: {} days {}'.format(numb_days.days, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))
print(' ... Number of entries collected:', uniparc_db.n_entries)
print(' ... Total RSS memory use so far (GB):',  round(psutil.Process(pid).memory_info().rss/1024/1024/1024, 3))
print(' ... Total RSS memory use for this step (GB):',  round((psutil.Process(pid).memory_info().rss - start_memory)/1024/1024/1024,3))

print(' ... Indexing!\n')
start = time.time()
uniparc_db.index_db()
numb_seconds = time.time() - start
print(' ... Time to index uniparc_db: {}'.format(time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))
print(' ... Total RSS memory use for this step (GB):',  round((psutil.Process(pid).memory_info().rss - start_memory)/1024/1024/1024,3))

# 4. Now extract the AlphaFold database information

print('\n4. TAKING CARE OF {}\n'.format(alphafold_path))

start = time.time()
day0  = date.today()

alphafold_db = alphafold.alphafold_extractor(mongo_host = MONGO_HOST, mongo_port = MONGO_PORT)
alphafold_db.register(alphafold.pLDDT_extractor())

print('\n ... Extracting!\n')
alphafold_db.extract(alphafold_path, max_size = n_to_extract, print_step = 10000, chunk_size = 10000)

numb_seconds = time.time() - start
numb_days    = date.today() - day0
print('\n ... Time to fill Alphafold db: {} days {}'.format(numb_days.days, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))
print(' ... Total RSS memory use so far (GB):', round(psutil.Process(pid).memory_info().rss/1024/1024/1024, 3))

print(' ... Indexing!\n')
start = time.time()
alphafold_db.index_db()
numb_seconds = time.time() - start
print(' ... Time to index alphafold_db: {}'.format(time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))
print(' ... Total RSS memory use for this step (GB):',  round((psutil.Process(pid).memory_info().rss - start_memory)/1024/1024/1024,3))

# 5. Now extract uniref50, add darkness and select the darkness representatives

print('\n5. TAKING CARE OF {}\n'.format(uniref50_file))

start = time.time()
day0  = date.today()

uniref50_db = uniref.uniref_extractor(mongo_host = MONGO_HOST, mongo_port = MONGO_PORT, name = 'UniRef50')
uniref50_db.register(uniref.entries_extractor())
uniref50_db.register(uniref.unirefs_extractor())

print('\n ... Extracting!\n')

clear = True
uniref50_db.extract(uniref50_file, clear = clear, max_size = n_to_extract, print_step = 1000, savestep = 5000, chunk_size = 10000, uniprot_db = uniprot_db, uniparc_db = uniparc_db, alphafold_db = alphafold_db, update_unip = True)

numb_seconds = time.time() - start
numb_days    = date.today() - day0

print('\n ... Time to fill uniref50 db: {} days {}'.format(numb_days.days, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))
print(' ... Total RSS memory use so far (GB):', round(psutil.Process(pid).memory_info().rss/1024/1024/1024, 3))

print(' ... Indexing!\n')
start = time.time()
uniref50_db.index_db()
numb_seconds = time.time() - start
print(' ... Time to index uniref50_db: {}'.format(time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))
print(' ... Total RSS memory use for this step (GB):',  round((psutil.Process(pid).memory_info().rss - start_memory)/1024/1024/1024,3))

## 6. Now extract the updated names of the uniprots after Google naming

print('\n6. TAKING CARE OF UPDATED UNIPROTS: {}\n'.format(uniprot_files_updated_names))

start = time.time()
day0  = date.today()

uniprot_db = uniprot.uniprot_extractor(mongo_host = MONGO_HOST, mongo_port = MONGO_PORT, name='UniProt_updated_names')
uniprot_db.register(uniprot.name_extractor())

start_memory = psutil.Process(pid).memory_info().rss

clear = False
for uniprot_file in uniprot_files_updated_names:

    print(' ... Extracting {}\n'.format(uniprot_file))
    uniprot_db.extract(uniprot_file, clear = clear, print_step = 10000, chunk_size = chunk_size, process_coverages = False)

    numb_seconds = time.time() - start
    numb_days    = date.today() - day0
    print(' ... ... Time to fill uniprot_db with {}: {} days {}'.format(uniprot_file, numb_days.days, time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))
    print(' ... ... Number of entries collected:',uniprot_db.n_entries)
    print(' ... ... Total RSS memory use so far (GB):',  round(psutil.Process(pid).memory_info().rss/1024/1024/1024, 3))
    print(' ... ... Total RSS memory use for this step (GB):',  round((psutil.Process(pid).memory_info().rss - start_memory)/1024/1024/1024, 3))
    clear = False

print(' ... Total RSS memory use so far (GB):', psutil.Process(pid).memory_info().rss/1024/1024/1024)

print(' ... Indexing!\n')
start = time.time()
uniprot_db.index_db()
numb_seconds = time.time() - start
print(' ... Time to index uniprot_db: {}'.format(time.strftime('%H hours %M min %S sec', time.gmtime(numb_seconds))))
print(' ... Total RSS memory use for this step (GB):',  round((psutil.Process(pid).memory_info().rss - start_memory)/1024/1024/1024,3))
