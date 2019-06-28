
# bedtools
tar -zxvf bedtools-2.26.0.tar.gz
cd bedtools2
make
cd bin/
export PATH=$PWD:$PATH

# STAR
tar -xzf 2.7.1a.tar.gz
cd STAR-2.7.1a
# Compile
cd STAR/source
make STAR
