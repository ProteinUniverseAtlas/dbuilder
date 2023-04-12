# dbuilder

This repository contains the code to extract data from multiple protein databases and organise it into a MongoDB. 

The main script is `extract_dbs_with_mongo.py`, which given a set of DB locations (which must be defined in the script), will generate a big, juicy mongoDB with multiple collections containing data in UniProtKB (Swiss-Prot + TrEMBL), InterPro, UniRef, UniParc and AlphaFoldDB.

The script relies on multiple modules found in `src`, which can be also used independently to access any data in the DB in follow-up scripts and jupyter notebooks (as those in the repository `ProteinUniverseAtlas/AFDB90v4`).

### How to use this repo

Just download the full content of the repo and install the dependencies below with `pip`.

To generate your own local database, make sure you have a Mongo server running and then run the `extract_dbs_with_mongo.py` script. It assumes there is a folder `databases` with the `.gz` files from uniprot and a full copy of the AlphaFold database. 

To access data in your database, check the examples in the `ProteinUniverseAtlas/AFDB90v4` repository.

### Dependencies apart from standard python libraries

The code was written in Python 3.6 and standard python modules were used. Extra modules required are:

- pymongo
- urllib
- numpy
- pickle
- gzip