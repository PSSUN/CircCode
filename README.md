# CircCode

### Introduction

CircCode is a Python3-base pipeline for translated circular RNA identification. It automatically tandem links sequence in series and processes a given ribosome profile data (including quality control, filtering and alignment). Finally, based on random forest and J48 classification, the final translated circular RNA was predicted. The user only needs to ***fill in the given configuration file*** and ***run the python scripts*** to get the predicted translated circular RNA.

### Requirement

#### Data:

- Genome sequence (fasta format)
- Candidate circRNA sequence (fasta/bed format)
- rRNA sequence (fasta format)
- Adapter sequence (fasta format)
- Ribosome profiling data (sra format)
- Coding and non-coding sequence (fasta format)
#### Software:

- bedtools (v.2.26.0+): (https://bedtools.readthedocs.io/en/latest/)
- bowtie (v.1.2.2+): (http://bowtie-bio.sourceforge.net/index.shtml)
- STAR (v.2.7.1+): (https://github.com/alexdobin/STAR)
- Python3 (v.3.6.5+): (https://www.python.org/)
- R language (v.3.4.4+): (https://www.r-project.org/)

#### python3 package:

- Biopython (v.1.72+): (https://pypi.org/project/biopython/)
- Pandas (v.0.23.3+): (https://pypi.org/project/pandas/)

#### R package:

- BASiNET: (https://github.com/cran/BASiNET)
- Seqinr: (https://cran.r-project.org/web/packages/seqinr/index.html)

### Usage

##### - You can run all CircCode pipeline by one script

1. Fill the config file (https://github.com/Sunpeisen/CircCode/blob/master/config.yaml), input full path of each required file.

2. Run bash script on command line  with your config file.

 ```bash
   sh one_step_script.sh config.yaml
 ```

##### - Or you can run CircCode step by step

  1. Fill the config file (https://github.com/Sunpeisen/CircCode/blob/master/config.yaml), input full path of each required file.

  2. Making virtual genomes

  ```python
   python3 make_virtual_genomes.py -y config.yaml
  ```
  3. Filter reads and compare to virtual genomes

  ```python
   python3 map_to_virtual_genomes.py -y config.yaml
  ```
  4. Find RPF-covered region on junction (RCRJ) and classification of RCRJ by sequence features

  ```python
   python3 find_RCRJ_and_classify.py -y config.yaml
  ```
  5. Find longest peptide of translated circRNA

 ```python
  python3 find_longest_pep.py -y config.yaml
 ```

### Run example
You can downlad the required sra file from [NCBI-SRA](https://www.ncbi.nlm.nih.gov/sra/SRR3495992), we also provide the other required files (includes genome.fa, genome.gtf etc.) in [example.tar.xz](https://github.com/PSSUN/CircCode/blob/master/example.tar.xz). Fill in the path of the corresponding file into the project corresponding to config.yaml. Then follow the steps mentioned above to run each script.

**NOTE**：The test file is only used to test whether the software can run smoothly and does not represent the actual research results.



### Contact us

If you encounter any problems while using CircCode, please send an email (sps@snnu.edu.cn / glli@snnu.edu.cn) or submit the issues on GitHub (https://github.com/Sunpeisen/circCode/issues) and we will resolve it as soon as possible.
