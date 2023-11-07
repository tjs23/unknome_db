unknome_db
----------

These scripts create an up-to-date version of the Unknome database, as a relational database in SQLite3 file format.
A helper script is included to create summary files and charts from a database file, as used on the [Unknome web site](https://unknome.mrc-lmb.cam.ac.uk). These scripts compile the data that is made available on the Unknome web site, but remain separate from the site itself.

The overall principle of the Unknome database is to assign a knownness score to proteins, representing the amount of human knowledge they have.
Proteins from the [UniProt](https://www.uniprot.org/) database are placed in a cluster of orthologues based on the [Panther database](https://pantherdb.org/).
The knowness score is defined as the largest number of [Gene Ontology](http://geneontology.org/) terms that has been assigned to a member of that cluster. Because GO annotations vary in confidence and relevance to function, different types of evidence aree assigned a different weight when calculating the score. 
  
Installation
------------

These scripts work with vanilla Python (v3) and have no requirements outside of the standard library.
They have no special in stallation requirements and can be run directly from a cloned repository.

If you move the `setup_unknome_db.py` Python script the `constants.py` and `sql_tables.py` modules must be
on the import PYTHONPATH.


Creating an Unkome Database
---------------------------

The `setup_unknome_db.py` is simply run from the command line using a Python3 interpreter, including the working directory (for storing GO and PANTHER data files etc) and any other command line options (see below).

Typical use:

  `python3 setup_unknome_db.py /save_dir_path/`


Input data
-----------

The input data, from which the Unknome database is created, will be downloaded from various Internet resources. 

These data, with their default locations, are as follows:

* [Gene Ontology evidence codes specifications](https://raw.githubusercontent.com/evidenceontology/evidenceontology/master/eco.obo)

* [UniProt-GO GOA .gaf mapping file](https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf)

* [PANTHER summery of reference proteomes](https://www.pantherdb.org/panther/summaryStats.jsp)

* [PANTHER classifications release](http://data.pantherdb.org/ftp/hmm_classifications/current_release/)

* [NCBI taxonomy .dmp files](https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.Z)

As time passes the default values for these URLs may no longer work.
However, all the URLs may be specified as command line options in the event that they change. 

Output data
-----------

The compiled database is generated as a single SQLite3 format file, which may be interrogated with the standard [sqlite3](https://docs.python.org/3/library/sqlite3.html) Python module. See the `sql_tables.py` file from an indication fo the SQL schema.

Command line options
--------------------

```
usage: setup_unknome_db [-h] [-o OUT_DB_FILE] [-f] [-eco FILE_URL]
                        [-ps WEB_URL] [-pd DIR_URL] [-gaf FILE_URL]
                        [-tax FILE_URL]
                        WORKING_DIR
```

```
positional arguments:
  WORKING_DIR     The working directory, into which downloaded data files and
                  the databse .sqlite file will be stored.
```

```
optional arguments:
  -h, --help      show this help message and exit
  -o OUT_DB_FILE  Optional output file path to write the SQLite3 database file
                  to. The default is to create a file in the working
                  directorly labelled by creation date, of the form
                  "unknome_db_{date}.sqlite".
  -f, --force     Whether to force aa re-download of all data files, if
                  already present.
  -eco FILE_URL   URL to download evidence code definitions, as used by GO, as
                  an "eco.obo" file. Default=https://raw.githubusercontent.com
                  /evidenceontology/evidenceontology/master/eco.obo
  -ps WEB_URL     URL for PANTHER proteome summary table, to obtain a list of
                  species.
                  Default=https://www.pantherdb.org/panther/summaryStats.jsp
  -pd DIR_URL     URL for the DIRECTORY containing the latest PANTHER HMM
                  classifications file. Default=http://data.pantherdb.org/ftp/
                  hmm_classifications/current_release/
  -gaf FILE_URL   URL for GeneOntology to UniProt linking GOA .gaf file. Defau
                  lt=https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_un
                  iprot_all.gaf
  -tax FILE_URL   URL for TAR GZIPPEd archive containing NCBI taxonomy .dmp
                  files. Default=https://ftp.ncbi.nih.gov/pub/taxonomy/new_tax
                  dump/new_taxdump.tar.Z
```

Reference and License
---------------------

If you use these scripts in your research please cite:

**"Functional unknomics: Systematic screening of conserved genes of unknown function"**
Joao Rocha, Satish Arcot Jayaram, Tim J Stevens, Nadine Muschalik,
Rajen D Shah, Sahar Emran, Cristina Robles, Matthew Freeman, Sean Munro
[PLoS Biol. 2023 Aug; 21(8): e3002222](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3002222) PMID: [37552676](https://pubmed.ncbi.nlm.nih.gov/37552676/)

We have chosen to apply the Creative Commons Attribution 4.0 International (CC BY 4.0) License to any and all copyrightable parts of the Unknome database. We make no warranties regarding the correctness of the data presented, and disclaim liability for any damages that may result from its use. Users of the data are solely responsible for compliance with any copyright restrictions, patents or other rights. All data is provided “as-is” without any warranty, expressed or implied. 

