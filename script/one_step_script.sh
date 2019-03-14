#!/bin/bash

yaml = $1

python3 make_virtual_genomes.py -y ${yaml}
python3 map_to_virtual_genomes.py -y ${yaml}
python3 find_RCRJ_and_classify.py -y ${yaml}
python3 find_longest_pep.py -y ${yaml}
