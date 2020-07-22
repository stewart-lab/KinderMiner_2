# **KinderMiner 2.0**

This repository provides the python script for KinderMiner 2.0, a general text mining system to find association between any two terms in ~30 million PubMed articles. This project is done by Stewart Computational Biology Group (https://morgridge.org/research/regenerative-biology/bioinformatics/) within Thomson Lab (https://morgridge.org/research/regenerative-biology/thomson-lab/) at Morgridge Institute for Research, Madison, WI, USA.

A local PubMed database version is required for executing KinderMiner 2.0. Please see the github project https://github.com/iross/km_indexer/releases/tag/v1.2 for details.  

---- COMPILE AND RUN ON THE COMMAND LINE ----

$ python kinderminer2.py TARGETS_FILE KEYPHRASE_FILE -o OUTPUT_DIRECTORY

---- EVALUATE AND RANK TARGETS ----

To evaluate at default Fisher exact test (FET) p-value, 1.0E-05:
$ python evaluate_fisher_exact_fetpvalue_and_ratio_sorted.py OUTPUT_DIRECTORY/OUTPUT_FILE OUTPUT_DIRECTORY/OUTPUT_FILE_EVALUATION_RESULT

To evaluate at cutomized FET p-value (ex. 0.05):
$ python evaluate_fisher_exact_fetpvalue_and_ratio_sorted.py OUTPUT_DIRECTORY/OUTPUT_FILE OUTPUT_DIRECTORY/OUTPUT_FILE_EVALUATION_RESULT 0.05
