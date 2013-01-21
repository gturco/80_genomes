80 genomes
==============
Python programs written to find gene loss and gain between genomes


Requirements
==============

  - Python 2.7 or higher  `pythonbrew <https://github.com/utahta/pythonbrew/>`_
  - `FuzzyWuzzy <https://github.com/seatgeek/fuzzywuzzy>`_

Run
============
python create_clusters.py

Searches annotation files for a term/regex expression and 
returns protein cluster files that contain the term along·
with any additional terms for that cluster.

Returns summary files containing clade information on presence or abasance of genes 
  - cluster_summary.txt
  - clade_summary.txt


--term: grabs all annotations containing this term
--cluster path to cluster file
--anno directory containing annotation files


Example::

  python create_cluster_files.py --term 'signal transduction' --cluster /EightyHalophiles/FilteredHomology/HomologyClusters_I20.txt --anno /EightyHalophiles/Annotations/





####finds localdups of a gene looks up and downstream for genes in same
cluster on same scaffold and within a certain dist from on another
--cluster", dest="cluster_data", help="path to cluster_data_dir"
--pad", dest="padding", default=2000, help="bps to look up and
downstream gene for local dups
python find_local_dups --cluster "cluster_data/" --pad 2000
