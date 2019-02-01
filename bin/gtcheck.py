#!/usr/bin/env python3

# gtcheck.py
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

import os
import pysam
from pysam import bcftools
import argparse
import tempfile
import re
import subprocess
import sys
import pandas
from collections import OrderedDict
import shutil
import atexit
import datetime 
import statistics
import multiprocessing

parser = argparse.ArgumentParser(description= """
DESCRIPTION
Check genotype consistency between two input bam files
""", formatter_class= argparse.RawTextHelpFormatter, prog= os.path.basename(__file__))

inputArgs= parser.add_argument_group('Input arguments')
snpArgs= parser.add_argument_group('Options to filter reference SNPs')
gtArgs= parser.add_argument_group('Options to filter genoptypes')
miscArgs= parser.add_argument_group('Miscellanea arguments')

inputArgs.add_argument('--vcf', '-V',
                   required= True,
                   help='''\
VCF file of SNPs to use for genotyping. Use '-' to read from stdin
                   ''')

inputArgs.add_argument('--bams', '-b',
                   required= True,
                   nargs= 2,
                   help='''\
The pair of BAM files to check for consistency
                   ''')

inputArgs.add_argument('--ref', '-r',
                   required= True,
                   help='''\
Indexed fasta reference file. It must be consistent with the reference file
used for alignment
                   ''')

gtArgs.add_argument('--mapq', '-q',
                   type= int,
                   default= 30,
                   help='''\
Skip alignments with mapq smaller than MAPQ. Default: %(default)s 
                   ''')

gtArgs.add_argument('--min-gq', '-gq',
                   type= int,
                   default= 30,
                   help='''\
Skip genotype calls with quality less the MIN_GQ. Default %(default)s
                   ''')

gtArgs.add_argument('--min-dp4', '-dp4',
                   type= int,
                   default= 10,
                   help='''\
Skip genotype calls supported by fewer than MIN_DP4 high-quality reads. Default %(default)s
                   ''')

snpArgs.add_argument('--af-tag', '-af',
                   type= str,
                   help='''\
Tag in the VCF file containing the allele frequency. This tag is usually the
'AF' (e.g in gnomAD) or 'CAF' (dbSNP). If not given, SNPs will not be filtered
for allele frequency (not recommended)
                   ''')

snpArgs.add_argument('--af-range', '-ar',
                   type= float,
                   nargs= 2,
                   default= [0, 1],
                   help='''\
Require the allele frequency to be within this range. Default %(default)s
                   ''')

miscArgs.add_argument('--keeptmp', '-k',
                   action= 'store_true',
                   help='''\
Keep the temporary directory and its files created during the execution instead
of deleting it. (Only use for debugging)
                   ''')

miscArgs.add_argument('--version', '-v', action='version', version='%(prog)s 0.2.0')

args= parser.parse_args()

# ------------------------------------

def is_valid_variant(variantRecord, af_tag, min_af, max_af):
    """True if the variant record is suitable for genotyping.
    
    af_tag: The tag to use for the population allele frequency of the REF
        allele
    
    min_af, max_af: Discard variant if allele ferquency is below or above
        these values.
    """
    if len(variantRecord.alleles) != 2:
        return False
    if len(variantRecord.alleles[0]) != 1 or len(variantRecord.alleles[1]) != 1:
        # Only keep single nucleotide variants
        return False
    if af_tag is not None:
        if af_tag not in variantRecord.info:
            return False
        af= float(variantRecord.info[af_tag][0])
        if af < min_af or af > max_af:
            return False
    return True

def genotype_bam(argd):
    """Call genotypes in bam file at the given SNPs. argd is a dictionary of
    arguments
    """
    # Extract reads in panel SNPs
    aln_f= pysam.AlignmentFile(argd['bam'])
    panel_bam= tempfile.NamedTemporaryFile(prefix= re.sub('\.bam$', '', os.path.basename(argd['bam'])) + '.', suffix= '.bam', delete= False, mode= 'w', dir= argd['tmpdir'])
    bname= os.path.splitext(panel_bam.name)[0]
    aln_out_us= pysam.AlignmentFile(panel_bam, 'wb', header= aln_f.header)

    panel= pysam.VariantFile(argd['snp_panel'])
    prev_reads= []
    for snp in panel:
        reads= []
        for read in aln_f.fetch(snp.chrom, snp.start, snp.stop):
            reads.append(read)
        # Remove reads already written
        outreads= [x for x in reads if x not in prev_reads]
        for read in outreads:
            aln_out_us.write(read)
        prev_reads= reads
    aln_f.close()
    aln_out_us.close()
    panel.close()

    # Pileup
    pysam.mpileup('-f', args.ref, '-g', '--min-MQ', str(args.mapq), '-o', bname + '.all.bcf', bname + '.bam', catch_stdout= False)
    bcftools.index(bname + '.all.bcf', catch_stdout= False)
    
    # Call genotypes
    # We use subprocess instead of pysam because of https://github.com/pysam-developers/pysam/issues/693
    cmd= ['bcftools', 'call', '-T', snp_panel, '-m', '--skip-variants', 'indels', '-O', 'z', '-o', bname + '.calls.vcf.gz', bname + '.all.bcf']
    sp= subprocess.Popen(cmd, stderr= subprocess.PIPE)
    stdout, stderr= sp.communicate()
    stderr= stderr.decode().split('\n')
    stderr= '\n'.join([x for x in stderr if 'assuming all sites are diploid' not in x]) # We ignore this warning
    if stderr != '':
        sys.stderr.write(stderr)
    
    if sp.returncode != 0:
        raise Exception('\n%s exit code from \n\n%s\n' % (sp.returncode, ' '.join(cmd)))
    
    # Make tabular format
    with open(bname + '.txt', 'w') as txt:
        txt.write('\t'.join(['chrom', 'pos', 'alt', 'gt', 'qual']) + '\n')
        calls= pysam.VariantFile(bname + '.calls.vcf.gz')
        for x in calls:
            if x.qual < argd['min_gq'] or sum(x.info['DP4']) < argd['min_dp4']:
                continue
            if x.alts is None:
                alt= '.'
            else:
                alt= x.alts[0]
            gt= x.samples[0]['GT']
            if gt[0] is None or gt[1] is None:
                continue
            line= [x.chrom, str(x.pos), alt, str(gt[0]) + "/" + str(gt[1]), str(round(x.qual, 1))]
            txt.write('\t'.join(line) + '\n')
        calls.close()    
    return {'table': bname + '.txt', 'bam': argd['bam']}

def med(x):
    if len(x) == 0:
        return float('nan')
    return statistics.median(x)

# -----------------------------------

dt= datetime.datetime.now().strftime('%y%m%d-%H%M%S-%f')[:-3]

tmpdir= tempfile.mkdtemp(prefix= 'gtc_%s_' % dt, suffix= '', dir= '.')
if not args.keeptmp:
    atexit.register(shutil.rmtree, tmpdir)
else:
    sys.stderr.write('Temp directory: %s\n' % tmpdir)

# Prepare panel VCF file
snp_panel= os.path.join(tmpdir, 'panel.vcf')
panel= pysam.VariantFile(args.vcf)
tmp= pysam.VariantFile(snp_panel, 'w', header= panel.header)
chrom= ''
pos= 0
nsnp= 0
for snp in panel:
    if not is_valid_variant(snp, af_tag= args.af_tag, min_af= args.af_range[0], max_af= args.af_range[1]):
        continue
    if chrom == snp.chrom and pos > snp.pos:
        raise Exception('VCF file %s is not sorted' % args.vcf)
    nsnp += 1
    tmp.write(snp)
    chrom= snp.chrom
    pos= snp.pos
tmp.close()

# Genotype
with multiprocessing.Pool(2) as p:
    genotypes= p.map(genotype_bam, [{'bam': args.bams[0], 'snp_panel': snp_panel, 'min_gq': args.min_gq, 'min_dp4': args.min_dp4, 'tmpdir': tmpdir}, 
                                    {'bam': args.bams[1], 'snp_panel': snp_panel, 'min_gq': args.min_gq, 'min_dp4': args.min_dp4, 'tmpdir': tmpdir}])

# Compare calls
df1= pandas.read_csv(genotypes[0]['table'], sep= '\t')
df2= pandas.read_csv(genotypes[1]['table'], sep= '\t')

# Merge throws a harmless warning if both dataframes are empty (https://github.com/pandas-dev/pandas/issues/17776)
mrg= df1.merge(df2, how= 'inner', on= ['chrom', 'pos'], suffixes= ['_1', '_2'])
mrg.to_csv(os.path.join(tmpdir, 'calls.mrg.tsv'), sep= '\t', index= False)

n_disc= len(mrg[mrg['gt_1'] != mrg['gt_2']])
if len(mrg) == 0:
    pct_disc= 'nan'
else:
    pct_disc= round(n_disc / len(mrg) * 100.0, 2)

med_match_1= med(mrg[mrg['gt_1'] == mrg['gt_2']]['qual_1'])
med_match_2= med(mrg[mrg['gt_1'] == mrg['gt_2']]['qual_2'])
med_disc_1= med(mrg[mrg['gt_1'] != mrg['gt_2']]['qual_1'])
med_disc_2= med(mrg[mrg['gt_1'] != mrg['gt_2']]['qual_2'])

# Report
sys.stdout.write('N_snp\t%s\t# SNPs in %s\n' % (nsnp, args.vcf))
sys.stdout.write('N_snp_1\t%s\t# SNPs called in %s\n' % (len(df1), genotypes[0]['bam']))
sys.stdout.write('N_snp_2\t%s\t# SNPs called in %s\n' % (len(df2), genotypes[1]['bam']))
sys.stdout.write('N_both\t%s\t# SNPs in both 1 and 2\n' % (len(mrg)))
sys.stdout.write('N_discord\t%s\t# Discordant genotypes\n' % n_disc )
sys.stdout.write('Pct_discord\t%s\t# %% discordant genotypes\n' % pct_disc)
sys.stdout.write('med_match_1\t%.2f\t# Median quality of matching genotypes in 1\n' % med_match_1)
sys.stdout.write('med_match_2\t%.2f\t# Median quality of matching genotypes in 2\n' % med_match_2)
sys.stdout.write('med_disc_1\t%.2f\t# Median quality of discordant genotypes in 1\n' % med_disc_1)
sys.stdout.write('med_disc_2\t%.2f\t# Median quality of discordant genotypes in 2\n' % med_disc_2)

