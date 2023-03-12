# dbuilder

This repository contains the code to extract data from multiple protein databases and organise it into a MongoDB. 

The main script is `extract_dbs_with_mongo.py`, which given a set of DB locations (which must be defined in the script), will generate a big, juicy mongoDB with multiple collections containing data in UniProtKB (Swiss-Prot + TrEMBL), InterPro, UniRef, UniParc and AlphaFoldDB.

The script relies on multiple modules found in `src`, which can be also used independently to access any data in the DB in follow-up scripts and jupyter notebooks (as those in the repository `ProteinUniverseAtlas/AFDB90v4`.


### Dependencies apart from standard python libraries

- pymongo
- urllib
- numpy
- pickle
- gzip
