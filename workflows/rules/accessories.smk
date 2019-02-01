# Rules to make files that typically are prepared just once. These rules are
# not expected to contain wildcards relative to libraries and the manifest file
# should not be read at all.
# 
# The output files may go in the same dir as the input. Make sure that such dir
# is writable.
#
# accessories.smk
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

localrules: createSequenceDictionary
rule createSequenceDictionary: 
    input: 
        config['ref'] 
    output: 
        os.path.splitext(config['ref'])[0] + '.dict', 
    shell: 
        """ 
        gatk CreateSequenceDictionary -R={input} -O={output}
        """ 
 
localrules: simpleRepeat 
rule simpleRepeat: 
    output: 
        bed= 'pindel/pon/simpleRepeat.bed.gz', 
        tbi= 'pindel/pon/simpleRepeat.bed.gz.tbi', 
    shell: 
        """ 
        curl -s -S http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz \\ 
        | zcat \\ 
        | awk '$6 <= 6' \\ 
        | cut -f 2- \\ 
        | bgzip > {output.bed} 
        tabix -f {output.bed} 
        """ 

localrules: igvChromSizes  
rule igvChromSizes: 
    input: 
        config['ref'] + '.fai' 
    output: 
        re.sub('\.fa$|\.fasta$', '.chrom.sizes', config['ref'])
    shell: 
        """ 
        cp {input} {output} 
        """ 
 
localrules: faidx 
rule faidx: 
    input: 
        config['ref'] 
    output: 
        config['ref'] + '.fai' 
    shell: 
        """ 
        samtools faidx {input} 
        """ 
 
rule bwaIndex: 
    input: 
        ancient(config['ref']) 
    output: 
        config['ref'] + '.bwt' 
    params: 
        log_dir= '.',
        jobname= 'bwaIndex', 
        cpus_per_task= 2,
        mem= 16000, 
    shell: 
        """ 
        bwa index {input} 
        """ 

