# install python package
pip3 install -r requirements.txt

# install bedtools
tar -zxvf bedtools-2.26.0.tar.gz
cd bedtools-2.26.0
make
cd bin/
export PATH=$PWD:$PATH

# Give the Permission
chmod a+x ./requiredSoft/STAR
chmod a+x ./requiredSoft/FragGeneScan

