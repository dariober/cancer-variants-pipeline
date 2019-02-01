# Manifest.py 
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

import re
import os
import pandas
import tempfile
import sys
import yaml
from lib import utils
from lib import errors

# Keep these constants in sync with the docs!
COLUMN_NAMES= {'sample', 'library', 'sample_type', 'fastq_1', 'fastq_2', 'bam', 'tn_pair'}
SAMPLE_TYPE= {'tumor', 'tumour', 'normal', 'pon'}

def _read_manifest(mnf):
    """Read and check manifest file and return a pandas dataframe with
    additional columns.
    """
    # First remove comment lines. We don't use pandas' option `comment`
    # since there may be '#' inside file or sample names (I've seen that!)
    fout= tempfile.NamedTemporaryFile(prefix= os.path.basename(mnf), delete= False, mode= 'w')
    with open(mnf) as fin:
        for line in fin:
            if not line.strip().startswith('#'):
                fout.write(line)                
    fout.close()
    manifest= pandas.read_csv(fout.name, sep= '\t')
    os.remove(fout.name)
    
    col= {}
    for x in COLUMN_NAMES:
        col[x]= x

    for col_name in list(manifest.columns):
        if col_name not in col:
            manifest= manifest.drop(col_name, axis= 1) # Drop all columns we don't recognize

    if col['sample'] not in manifest or not all(manifest['sample'].notna()):
        raise errors.InvalidManifestError('"sample" column cannot contain missing values')

    if col['library'] not in manifest:
        manifest['library']= manifest['sample']
    else:
        for index, row in manifest.iterrows():
            if manifest.at[index, 'library'] is None:
                manifest.at[index, 'library']= row['sample']

    # Columns may be missing but if they are present they must be sane
    if col['sample_type'] in manifest:
        manifest['sample_type']= [x.lower().strip() for x in manifest['sample_type']]
        manifest['sample_type']= [re.sub('^tumor$', 'tumour', x)  for x in manifest['sample_type']]
        for x in set(list(manifest['sample_type'])):
            if x not in SAMPLE_TYPE:
                raise errors.InvalidManifestError('Invalid "sample_type" in manifest file: "%s". Allowed are: [%s]' %(x, ', '.join(SAMPLE_TYPE)))
    
    if (col['fastq_1'] in manifest and col['fastq_2'] not in manifest) or (col['fastq_2'] in manifest and col['fastq_1'] not in manifest):
        sys.stderr.write('Columns fastq_1 and fastq_2 must be both present or absent')
        raise errors.InvalidManifestError()

    return manifest

class Manifest:
    """Class to query the sample data and get information about
    librarie and files.
    manifest: 
        File name of the manifest file
    """
    def __init__(self, manifest= None):
        
        self.manifest= manifest

        if manifest is not None:
            self.df= _read_manifest(manifest)
        else:
            self.df= pandas.DataFrame({'sample': []})

        self.sample= set(self.df['sample'])

        if 'bam' in self.df:
            self.df.bam= self._make_abspath(self.df.bam)
        if 'fastq_1' in self.df:
            self.df.fastq_1= self._make_abspath(self.df.fastq_1)
        if 'fastq_2' in self.df:
            self.df.fastq_2= self._make_abspath(self.df.fastq_2)

        self.fastq= pandas.DataFrame({'pair_id': [], 'sample': [], 'library': [], 'fastq_1': [], 'fastq_2': []})
        self.fastqc= pandas.DataFrame({'fastq': [], 'fastq_md5': [], 'sample': []})
        if 'fastq_1' in self.df:
            self.fastq= self.df[['library', 'sample', 'fastq_1', 'fastq_2']].drop_duplicates()
            self.fastq['pair_id']= None
            self.fastq= self.fastq[(self.fastq.fastq_1.notna()) & (self.fastq.fastq_2.notna())].sort_values(by= ['sample'])
            for index, row in self.fastq.iterrows():
                if row['fastq_1'] == row['fastq_2']:
                    sys.stderr.write('\nFirst-in-pair and second-in-pair fastq file are the same: %s' % row.to_string())
                    raise errors.InvalidManifestError()
                self.fastq.at[index, 'pair_id']= re.sub('\.fastq$|\.fq$', '', re.sub('\.gz$|\.bz2$', '', os.path.basename(row['fastq_1']))) + '.' + str(index)
            # pair_id unique as if a primary key
            assert len(self.fastq.pair_id) == len(set(self.fastq.pair_id))
            
            # Assign output directories for QC
            for index, row in self.fastq.iterrows():
                for fq in ['fastq_1', 'fastq_2']:
                    self.fastqc= self.fastqc.append({'fastq': row[fq], 
                                                     'sample': row['sample'], 
                                                     'fastq_md5': os.path.join(self._getFastqcDir(row[fq]), os.path.basename(row[fq])) + '.md5'}, ignore_index= True)
            # As if a primary key
            assert len(self.fastqc.fastq_md5) == len(set(self.fastqc.fastq_md5))

        self.contrasts= pandas.DataFrame({'tumour_smp': [], 'normal_smp': []})
        if 'sample_type' in self.df and 'tn_pair' in self.df:
            # Match up pairs 
            tn_pairs= set(self.df['tn_pair'])
            tn_pairs= [x for x in tn_pairs if x is not None and x != '' and x != '.']
            for p in tn_pairs:
                pdf= self.df[self.df['tn_pair'] == p][['sample', 'sample_type']].drop_duplicates()
                if len(pdf) > 2:
                    sys.stderr.write('\nToo many samples associated to one contrasts:\n%s\n' % pdf.to_string())
                    raise errors.InvalidManifestError()
                if len(pdf) < 2:
                    continue
                if len(set(pdf['sample_type'])) != 2:
                    sys.stderr.write('\nExpected one "tumour" and one "normal" sample type. Got\n%s\n' % pdf.to_string())
                    raise errors.InvalidManifestError()
                normal_smp= pdf[pdf['sample_type'] == 'normal']['sample'].iloc[0]
                tumour_smp= pdf[pdf['sample_type'] == 'tumour']['sample'].iloc[0]
                self.contrasts= self.contrasts.append({'tumour_smp': tumour_smp, 'normal_smp': normal_smp}, ignore_index=True)

    def __str__(self):
        return '\n'.join(['Samples:', str(sorted(self.sample)), 
                          '\nFastq files:', self.fastq.sort_values(['library', 'sample']).to_string(), 
            '\nFastq quality control:', self.fastqc.sort_values(['sample']).to_string(), 
            '\nTumour-Normal contrasts:', self.contrasts.sort_values(['normal_smp', 'tumour_smp']).to_string()])

    def pair_ids_for_sample(self, sample):
        if self.fastq.empty:
            return []
        else:
            return list(self.fastq[self.fastq['sample'] == sample]['pair_id'])

    def get_fq1(self, pair_id):
        if not self.fastq.empty:
            return list(self.fastq[self.fastq['pair_id'] == pair_id]['fastq_1'])[0]
        else:
            return ''

    def get_fq2(self, pair_id):
        if not self.fastq.empty:
            return list(self.fastq[self.fastq['pair_id'] == pair_id]['fastq_2'])[0]
        else:
            return ''

    def symlink_bam(self, aligner= 'bwa'):
        """For each bam file in manifest, create a symlink to the file that would
        be created by the pipeline given the information in the sample sheet.
        """
        if 'bam' in self.df:
            bamdf= self.df[self.df['bam'].notna()].drop_duplicates()
            check_= bamdf[['sample', 'bam']].drop_duplicates()
            if len(check_['sample']) != len(set(check_['sample'])):
                raise errors.InvalidManifestError('A sample has been assigned to multiple bam files.\nIf this is intentional, please report it as a feature request.')

            for index, row in bamdf.iterrows():
                source= os.path.abspath(row['bam'])
                # dest_dir must be consistent with the output of the pipeline!
                dest_dir= os.path.join(os.getcwd(), row['sample'], aligner)
                os.makedirs(dest_dir, exist_ok= True)
                # Name of the bam file must be consistent with the pipeline naming!
                bamfilename= row['sample'] + '.bam'
                dest= os.path.join(dest_dir, bamfilename)
                utils.symLink(source, dest)
                bai= utils.get_bam_index(source)
                utils.symLink(bai, dest + '.bai')

    def _getFastqcDir(self, fq, root_dir= '.'):
        """Get a suitable output directory for FastQC given this fastq file.
        """
        sample= set(self.fastq[self.fastq['fastq_1'] == fq]['sample'])
        if len(sample) == 0:
            sample= set(self.fastq[self.fastq['fastq_2'] == fq]['sample'])

        if len(sample) != 1:
            raise errors.InvalidManifestError('Exactly 1 sample must be associated to fastq file "%s". Got %s' % (fq, str(list(sample))))
        sample= list(sample)[0]

        # Get all the fastq files associated to this sample
        fqs= set(list(self.fastq[self.fastq['sample'] == sample]['fastq_1']) + list(self.fastq[self.fastq['sample'] == sample]['fastq_2']))
        fqs= [x for x in fqs if x is not None]
        fqs= [os.path.abspath(x) for x in fqs]
        # Get the shortest unique filename. Almost always, this will be the 
        # basename of the fastq file. But we could have:
        # /data/a/reads[.fastq.gz] <- fq
        # /data/b/reads[.fastq.gz]
        # In this case we return 'a'
        sfq= [x.split(os.sep)[::-1] for x in fqs]
        depth= 1
        while True:
            base= [os.sep.join(x[0:depth]) for x in sfq]
            if len(base) == len(set(base)):
                # All unique, this depth is suffient
                break
            depth += 1
        pfq= os.sep.join(os.path.abspath(fq).split(os.sep)[::-1][0:depth][::-1])
        pfq= os.path.dirname(pfq) # We remove the name of the fastq file itself since fastqc will put it.
        fqdir= os.path.join(root_dir, sample, 'fastqc', pfq)
        fqdir= os.path.normpath(fqdir)
        return fqdir

    def _make_abspath(self, files):
        """Make the filenames in files absolute 
        """
        root= os.path.abspath(os.path.dirname(self.manifest))
        absfiles= []
        for x in files:
            if x is None or pandas.isnull(x) or x.strip() == '':
                absfiles.append(x)
            elif os.path.isabs(x):
                absfiles.append(x)
            else:
                absfiles.append(os.path.abspath(os.path.join(root, x)))
        return absfiles
    
    def get_tumour_samples(self):
        """Return the set of samples classified as tumour type
        """
        if 'sample_type' in self.df:
            return set(self.df[self.df.sample_type == 'tumour']['sample'])
        else:
            return set()
