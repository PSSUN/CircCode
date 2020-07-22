#!usr/bin/bash

do
python3 make_virtual_genomes.py -y config.yaml
python3 map_to_virtual_genomes.py -y config.yaml
python3 find_RCRJ_and_classify.py -y config.yaml
done
