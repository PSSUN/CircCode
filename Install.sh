# install python package
pip3 install pandas
pip3 install Biopython

# install bedtools
tar -zxvf bedtools-2.26.0.tar.gz
cd bedtools-2.26.0
make
cd bin/
export PATH=$PWD:$PATH
