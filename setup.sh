#!/usr/bin/env bash

# setup.sh
# 
# Copyright (C) 2019 University of Glasgow
# 
# Author: Dario Beraldi <dario.beraldi@glasgow.ac.uk>
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

VERSION=0.1.0

set -e
set -o pipefail

# Parse arguments
# ===============

PG=`basename "$0"`
bin_dir=${HOME}/bin

# From https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -b|--bin_dir)
        bin_dir="$2"
        shift # past argument
        shift # past value
    ;;
    -v|--version)
        version=1
        shift # past argument
    ;;
    -h|--help)
        help=1
        shift # past argument
    ;;
    *)
        echo "Unknown option in $@"  
        exit 1
    shift # past argument
    ;;
esac
done

if [[ $help == 1 ]]
then
cat <<EOF
DESCRIPTION
Installer of several required programs. Only programs not found on PATH will be
installed. See code in this script for programs and versions.

-b|--bin_dir  Install missing programs here. This dir should writable and on
              your PATH. Default $bin_dir
-v|--version  Show version
-h|--help     Show help

USAGE EXAMPLE
bash setup.sh -b $bin_dir

Version $VERSION
EOF
exit 0
fi

if [[ $version == 1 ]]
then
    echo "$PG $VERSION"
    exit 0
fi
# End argument parsing
# ====================

cwd=`pwd`
mkdir -p downloads
mkdir -p bin

# pindel (?)
# pysam <- Rewrite requesting package(s) in java?

function install_htslib(){
    # Download and install htslib. Compiled stuff is in `pwd`/htslib 
    pushd .
    rm -f htslib-1.8.tar.bz2
    wget https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2
    tar xf htslib-1.8.tar.bz2
    rm htslib-1.8.tar.bz2
    mv htslib-1.8 htslib
    cd htslib
    ./configure --prefix=`pwd`
    make -j 4
    make install
    popd 
}

function check_bin_dir(){
    # USAGE: check_bin_dir /foo/bar
    # 
    # Check `/foo/bar` is on PATH. If it is on PATH but it does not exist,
    # try to create it
    python3 -c "import os, sys
PATH= os.environ['PATH'].split(os.path.pathsep)
PATH= [os.path.abspath(x) for x in PATH]
bin= os.path.abspath('$1')
if bin not in PATH:
    sys.stderr.write('\n\033[31mError: bin directory requested by -b/--bin_dir \'$1\' is not on PATH\033[0m\n')
    sys.exit(1)
os.makedirs(bin, exist_ok=True) 
"
}
check_bin_dir ${bin_dir}

# Start installing missing bits
# =============================

set -x 

# HTSLIB (tabix & bgzip)
found=`command -v tabix` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    install_htslib
    cp htslib/tabix ${bin_dir}/
    cp htslib/bgzip ${bin_dir}/
fi
command -v tabix
tabix --version
command -v bgzip
bgzip --version

# filterSomaticPindel
chmod a+x bin/filterSomaticPindel.py
rsync --update bin/filterSomaticPindel.py ${bin_dir}
filterSomaticPindel.py --version

# Picard. 
# Note that we search for it on PATH and we put it in bin if not found. 
path=`echo $PATH | tr ':' ' '`
found=`find $path -readable -name 'picard.jar' | head -n 1` || true
if [[ -z $found ]]
then
    rm -f picard.jar
    wget https://github.com/broadinstitute/picard/releases/download/2.18.2/picard.jar
    mv picard.jar ${bin_dir}/
fi
found=`find $path -readable -name 'picard.jar' | head -n 1` || true
java -jar $found MarkDuplicates --version || true

# IGVTools
found=`command -v igvtools` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    rm -f igvtools_2.3.98.zip
    wget http://data.broadinstitute.org/igv/projects/downloads/2.3/igvtools_2.3.98.zip
    unzip -q igvtools_2.3.98.zip
    rm igvtools_2.3.98.zip
    cp IGVTools/igvtools ${bin_dir}/
    cp IGVTools/igvtools.jar ${bin_dir}/
fi
command -v igvtools
igvtools help

# Manta
found=`command -v ${cwd}/bin/manta/bin/configManta.py` || true
if [[ -z $found ]]
then
    cd ${cwd}/bin
    rm -rf manta-1.3.2.centos6_x86_64.tar.bz2 manta-1.3.2.centos6_x86_64
    wget https://github.com/Illumina/manta/releases/download/v1.3.2/manta-1.3.2.centos6_x86_64.tar.bz2
    tar xf manta-1.3.2.centos6_x86_64.tar.bz2
    rm manta-1.3.2.centos6_x86_64.tar.bz2
    mv manta-1.3.2.centos6_x86_64 manta
fi
${cwd}/bin/manta/bin/configManta.py -h

# BBDuk
found=`command -v bbduk.sh` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    rm -f BBMap_37.98.tar.gz
    wget --no-check-certificate https://sourceforge.net/projects/bbmap/files/BBMap_37.98.tar.gz
    tar xf BBMap_37.98.tar.gz
    rm BBMap_37.98.tar.gz
    ln -s `pwd`/bbmap/bbduk.sh ${bin_dir}/
fi
command -v bbduk.sh
bbduk.sh --help

# FastQC
found=`command -v fastqc` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    rm -f fastqc_v0.11.7.zip
    wget --no-check-certificate https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
    unzip -q -o fastqc_v0.11.7.zip
    rm fastqc_v0.11.7.zip
    ln -s `pwd`/FastQC/fastqc ${bin_dir}/
    chmod a+x `pwd`/FastQC/fastqc
fi
command -v fastqc
fastqc --version

# BWA
found=`command -v bwa` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    rm -f v0.7.17.tar.gz
    wget https://github.com/lh3/bwa/archive/v0.7.17.tar.gz
    tar xf v0.7.17.tar.gz
    rm v0.7.17.tar.gz
    cd bwa-0.7.17
    make -j 4
    cp bwa ${bin_dir}/
fi
command -v bwa
bwa || true # bwa doesn't have a --version or --help option

# GATK4
found=`command -v gatk` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    rm -f gatk-4.0.4.0.zip
    wget https://github.com/broadinstitute/gatk/releases/download/4.0.4.0/gatk-4.0.4.0.zip
    unzip -q gatk-4.0.4.0.zip
    rm gatk-4.0.4.0.zip
    ln -s `pwd`/gatk-4.0.4.0/gatk ${bin_dir}/
fi
command -v gatk
gatk Mutect2 --version

# Snakemake
found=`command -v snakemake` || true
if [[ -z $found ]]
then
    pip3 install --user 'snakemake==4.8.0'
fi
command -v snakemake
snakemake --version

# Python/Pandas
found=`python3 -c "import pandas"` || true
if [[ ! -z $found ]]
then
    pip3 install --user pandas
fi
python3 -c "import pandas; print(pandas.__version__)"

# VEP
found=`command -v vep` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    export PERL_MM_USE_DEFAULT=1
    rm -f 92.1.tar.gz
    wget https://github.com/Ensembl/ensembl-vep/archive/release/92.1.tar.gz
    tar xf 92.1.tar.gz
    rm 92.1.tar.gz
    cd ensembl-vep-release-92.1
    perl INSTALL.pl --NO_UPDATE --NO_HTSLIB --AUTO a 
    ln -s `pwd`/vep ${bin_dir}/
fi
command -v vep
vep --help

# Samtools
found=`command -v samtools` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    rm -f samtools-1.8.tar.bz2
    wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2
    tar xf samtools-1.8.tar.bz2
    rm samtools-1.8.tar.bz2
    cd samtools-1.8
    ./configure --prefix=`pwd`
    make -j 4
    cp samtools ${bin_dir}/
fi
command -v samtools
samtools --version


# bcftools
found=`command -v bcftools` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    rm -f bcftools-1.8.tar.bz2
    wget https://github.com/samtools/bcftools/releases/download/1.8/bcftools-1.8.tar.bz2
    tar xf bcftools-1.8.tar.bz2
    rm bcftools-1.8.tar.bz2
    cd bcftools-1.8
    ./configure --prefix=`pwd`
    make -j 4
    cp bcftools ${bin_dir}/
fi
command -v bcftools
bcftools --version

# BEDTOOLS
#found=`command -v bedtools` || true
#if [[ -z $found ]]
#then
#    cd ${cwd}/downloads
#    rm -f bedtools-2.27.1.tar.gz
#    wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz 
#    tar xf bedtools-2.27.1.tar.gz 
#    rm bedtools-2.27.1.tar.gz
#    cd bedtools2
#    make -j 4
#    cp bin/* ${bin_dir}/
#fi
#command -v bedtools
#bedtools --version

# facets
found=`command -v cnv_facets.R` || true
if [[ -z $found ]]
then
    cd ${cwd}/downloads
    rm -rf cnv_facets
    git clone https://github.com/<FIXME>/cnv_facets.git OK
    cd cnv_facets
    bash setup.sh --bin_dir ${bin_dir}
fi
command -v cnv_facets.R
cnv_facets.R --help

set +x
echo -e "\n\033[32mSetup successfully completed\033[0m\n"
