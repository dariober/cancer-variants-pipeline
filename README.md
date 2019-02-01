Pipeline for somatic variant detection from Whole Genome Sequencing 
===================================================================

A comprehensive detection of somatic variants in tumour-normal pairs requires
several steps working in concert in order to cover the spectrum of variant
types. *medusa* combines all these steps starting from raw fastq files to
variant calls in a transparent and customizable way.

Quick start
-----------

* Install with conda/bioconda:

```
conda install medusa
```

* Run pipeline

```
./workflows/medusa.py -m manifest.txt \
    -p \
    -gb hg38 \
    --mutect \
    --ref genome.fa
```

Adding program options
----------------------

The pipeline involves several programs and each program takes several command
line options. Most of these options are not exposed to the user in the form of
configurable items (*i.e.* settable via `-C/--config` option) since including
all of them, probably more then a hundred, would make the command line
interface confusing, difficult to maintain and of limited utility since most of
the options are left as default. However, exposing a program option should be
simple in most cases. There are two steps to follow:

* Register the option as a configuration item in `config_opts.yaml`.

* Edit the rule(s) where the option is to be applied.

Adding analysis module
----------------------

E.g. show how to add another variant caller, like Strelka.

* Edit `tasks.yaml` to add the requested output files. These are just the final
output, not the intermediate files that may be produced. For example, an
alignment module may produced: `unsorted.sam`, then
`sorted.bam`, then `sorted_markdups.bam` with associated index
`sorted_markdups.bam.bai`. If we are interested only in the sorted & marked
BAM and its index, we need to add in `tasks.py` only the files
`sorted_markdups.bam` and `sorted_markdups.bam.bai`. In fact, for brevity we
can add only the name of the index file since producing the index necessarily
produces also the associated BAM.

* Define the rules that will produce the output file(s) requested in `tasks.yaml`.
We create a script in `workflows/rules` that contains the necessary rules. Strictly
speaking, these new rules may be written in already existing rule files and need to 
be inside `workflows/rules`. Of course, if the input of this new task comes from rules already 
defined, it is not necessary to repeat them. For example, the input of a new variant caller
may be BAM files sorted and indexed. Since rules producing sorted & indexed BAMs already exists,
we don't need to worry about them.

* If the newly defined rules require some input from the user, edit
`config_opts.yaml` to add the necessary items. For example, if a tool
requires a VCF of reference SNPs we add a configuration item named 'refSnps'
and we execute snakemake with `snakemake -C refSnps=/my/ref.vcf ...`

* Add `include` statement in `main.smk` to import the new module.

Useful HOW-TOs
--------------

* How to find all the output files produced by an analysis module for all
samples. Use the Unix utilities `find` or `ls`.  For example, assuming the
output directory of the pipeline is `workflow`, get all the VCF files produced
by mutect:

```
ls workflow/*/bwa/mutect/*/*.vcf.gz
```

* How to prepare the cache directory for VEP. 

Reference [vep cache](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#pre)

```
cd /my/ref/vep
curl -O ftp://ftp.ensembl.org/pub/release-92/variation/VEP/homo_sapiens_vep_92_GRCh38.tar.gz
tar xzf homo_sapiens_vep_92_GRCh38.tar.gz
```

As documented [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#fasta)

> The first time you run VEP with a specific FASTA file, an index will be
> built. This can take a few minutes, depending on the size of the FASTA file
> and the speed of your system. On subsequent runs the index does not need to
> be rebuilt (if the FASTA file has been modified, VEP will force a rebuild of
> the index).

Then use `vep_dir_cache=/my/ref/vep` as configuration key to enable annotation.

* How to get reference SNPs for FACETS. VCF from [dbSNP](ftp://ftp.ncbi.nih.gov/snp/organisms) have been used previously.
For example, to get the (current at this writing) reference dbSNP for hg38:

```
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/common_all_20180418.vcf.gz
```
