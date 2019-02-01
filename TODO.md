* Add a comment/header section at the beginning of each `rules/*.smk` file
  so a user opening these file to e.g. modify the programs' options knows what he/she is
  looking at.

* Document wildcards: Create a `wildcards.yaml` file which describe each
  wildcard. Some validation could be applied. E.g. check there are no
  undocumented wildcards in Task.py and the wildcards in `wildcards.yaml` are
  actually used in `Task.py`. Or maybe just add a comment section to `Task.py`
  or a link to a wiki page?

* Add diagnostic plots for mutect and possibly other output like samtools stats.
  It could be done on a sample by sample basis (easier) and/or collectively using
  `ggplot::facets_wrap` (more useful). Diagnostics could be: Histogram of VAF,
  number of variants, ...

* Option to analyse selected genomic regions. This option overrides option
  *chroms*. It can be a single bed file or a space separated list of regions.
  Each region is in the usual format `chr` or `chr:start-end`.

* Test scripts: Add option to skip cluster submission so that if the cluster is
  busy we don't have to wait for the submitted jobs to complete.

* See if supporting single-end reads is relatively easy. You would need to:
    * Parse the manifest to allow for `fastq_2` column to have missing values,
      indicating single-end sequencing. [Easy]
    * Make the `bbduk.sh` and `bwa` command strings conditional on pairing. [Easy]
    * Some tools may require paired-end reads (structural variants callers). Remove 
      from the output files where single-end reads are given. [Possibly tricky]
