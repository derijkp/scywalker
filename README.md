Scywalker
========= 
A program for the analysis of single cell Oxford nanopore long read data
Copyright VIB and University of Antwerp

Scywalker
-------
scywalker is a package designed to analyse single cell (10x) Oxford
nanopore long read data. (without the need for matching short read data).
It provides end-to-end analysis in one command: Starting from fastqs, it
will find and assign cellbarcodes, align reads, and reconstruct (based on
[IsoQuant](https://github.com/ablab/IsoQuant)) and quantify isoforms and
genes, producing both bulk and per cell counts.

If cell specific markersets are provided, it will also assign cell-types
and generate pseudobulk (per cell-type) counts for each sample, and make
count files allowing comparison between samples for these pseudobulk
counts.

scywalker is implemented within
[genomecomb](https://github.com/derijkp/genomecomb/). The scywalker
distribution contains all needed dependencies, including the full
genomecomb. This provides, besides the core scywalker algorithms, also
several commands/tools useful for scywalker set-up, analysis and
downstream analysis.

Installation
------------
Binary packages for Linux can be downloaded from github
([https://github.com/derijkp/scywalker](https://github.com/derijkp/scywalker))

Scywalker is distributed as a portable application directory: A
self-contained directory with the scywalker executable (scywalker) and all
needed depencies compiled in a way that they should work on all (except very
ancient) Linux systems.

Installation of the package is as simple as downloading the
[distribution](https://github.com/derijkp/scywalker/releases/download/v0.108.0/scywalker-0.108.0-linux-x86_64.tar.gz)
from github
([https://github.com/derijkp/scywalker](https://github.com/derijkp/scywalker))
and unpacking it, e.g.: 
```
cd ~/bin
wget https://github.com/derijkp/scywalker/releases/download/v0.108.0/scywalker-0.108.0-linux-x86_64.tar.gz
tar xvzf scywalker-0.108.0-linux-x86_64.tar.gz
rm scywalker-0.108.0-linux-x86_64.tar.gz
```

You can call the executables (scywalker, cg) directly from the directory
using the path (e.g. `~/bin/scywalker-0.108.0-linux-x86_64/scywalker ..`) 
or by placing the directory in the PATH environment variable (e.g. using 
`export PATH=~/bin/:$PATH`)
You can also place soft-links to the executables in a directory already in
the PATH. (remark: The executable itself needs to stay in the application
directory to find it's dependencies), e.g.
```
cd ~/bin
ln -s scywalker-0.108.0-linux-x86_64/scywalker .
ln -s scywalker-0.108.0-linux-x86_64/scywalker_makerefdir .
ln -s scywalker-0.108.0-linux-x86_64/cg .
```

Scywalker is largely implemented within [genomecomb](https://github.com/derijkp/genomecomb), 
and its distribution comes with an appropriate full version of genomecomb,
which can be run using the cg executable, providing  which also provides
multiple usefull extra tools for querying tsv files, etc.

Example/test run
----------------
As an example/test, the following code shows you how to download an example data set and run scywalker on it:
```
# download and unpack test data
wget https://github.com/derijkp/scywalker/releases/download/v0.108.0/scywalker_test.tar.gz
tar xvzf scywalker_test.tar.gz

cd scywalker_test
# make refdir; This test data is limited to chromosome 17, so there are no organelles included
# The command will accept compressed source fasta and gtf files (.gz, .zst, ..)
scywalker_makerefdir -organelles '' g17 genome.fa genes.gtf

# run scywalker using 8 cores on local machine, adapt this number to what you have available 
# on your machine (or use e.g. sge to run on a grid engine cluster)
# Of course we cannot properly determine celltypes on this limited data set.
# The "marker" genes in the included markers_chr17.tsv are not good
# markers, they are just made to allow celltyping to at least run on this specific limited dataset
scywalker -v 1 -d 8 \
	-refdir g17 \
	-sc_expectedcells 183 \
	-cellmarkerfile markers_chr17.tsv \
	-threads 6 \
	test10x
```

Reference data
--------------
A scywalker analysis needs a reference genome and a set of known isoforms
in a specific format, with different types of indexes and supporting
files. These must be provided in a reference directory.

You can use genomecomb reference directories for this; These can be downloaded from
the genomecomb website for a number of species (or created new) as
described in the [genomecomb installation
documentation](https://derijkp.github.io/genomecomb/install.html)

You can also use the included command scywalker_makerefdir to make a reference
directory starting from a fasta and a gtf transcript file.
If there are organelles in the genome sequence, it is important to specify them so the organelle
specific algorithm can be used: The isoquant based code often hangs or crashes on the
very different organelle data.
```
scywalker_makerefdir -organelles organelleslist refdir genomesequence.fasta transcripts.gtf
```
where 
* `organelleslist` is a space separated list of chromosomes that represent organelles, e.g. 'chrM chrPt',
* `genomesequence.fasta` is a multifasta file with the genomesequence and transcripts.gtf
* `transcripts.gtf` is a gtf file with transcripts for the given genome sequence. It is also possible
to give a (genomecomb) gene tsv file here.

Sample data
-----------
You provide sample data to scywalker in the form of a **sample directory**:
This is a directory that has (at least) a subdirectory named fastq. This
fastq directory must contain the fastq files (or softlinks to them) The
sample name is determined by the name of the sample directory. Results of the
analysis specific for the sample will be added in this directory

scywalker can analyse multiple samples in one run by providing it a
**project directory**. This is a directory (name of the dir determines the
name of the run/project) that at least contains a subdirectory named
samples. This samples subdir contains a sample directory from each sample
in the run. On analysis of a projectdir, all samples are analysed
individually, and files providing comparisons of multiple samples will be
made in a subdirectory compar.

The starting project directory should look thus like:

* project_directory/
  * sample1/
    * fastq/
      * file1.fq.gz
      * file2.fq.gz
  * sample2/
    * fastq/
      * file1.fq.gz
      * file2.fq.gz

Running scywalker
-----------------
You can run scywalker using the following command
```
scywalker ?options? sampledir/projectdir
```
This will analyse the sample in **sampledir** or all samples in
**projectdir** using the reference data in **refdir**. Results of the
analysis (or intermediate files) are added in place in the sampledir when
finsihed. If analysis is interupted or has an error, you can continue the
analysis by issuing the same command (after fixing what caused the error)

Options typically included (some required) for basic analysis of a 10x v3 (default)
data set would be
`-refdir`
    reference directory with genomesequence, etc. as described previously. This option is required.

`-sc_expectedcells`
    gives the the number of cells expected. This information is required for filtering cells using emptydrops
    However when processing a projectdir, this number is not always the same for all samples. You
    can give different values for this option by writing a tsv (tab separated value) file named
    options.tsv in the projectdir with the following fields: sample option value
    For each sample (that differs from the general option if given) you add a line with the samplename, 
    the option (sc_expectedcells without the -) and the value (the number of expected cells in this case)

`-cellmarkerfile`
    A tsv file providing genes that are indicative specific cell cell types. If not given, scsorter will not
    be used to determine celltypes and making pseudobulk files
    It can contain the following fields:
    **marker**: gene indicative of celltype (obligatory field)
    **celltype**: celltype for which gene is indicative (obligatory field)
    **tissue**: can be used in combination with the -tissue option; only markers of the given tissue will be used
    **markertype**: is the marker expressed "up" or not expressed "down" in celltype (only used by sctype currently)
    **weight**: weight of the marker (only used by sctype currently)

`-tissue`
    The tissue type of the sample. If cellmarkerfile is given, only markers of the given tissue are used.
    If cellmarkerfile is not given, scsorter is not run; however sctype will use its internal database 
    with the given tissue. If neither is given, no celltyping or pseudobulk generation is done.

`-d`
    By default the command is run using a single core (=slow). Use the `-d` option to specify the
    manner of job distribution/parallelisation. Use a number to specify distribution
    over (max) the given number of cores on the local machine, while "sge" or
    "slurm" will distribute jobs over a Grid Engine or SLURM cluster.
    On a cluster the command will finish after submitting all jobs (with dependencies).
    For distributed runs a tab separated log file is created (in the projectdir) named
    process_project_<projectname>.<starttime>.running
    This contains information on all started jobs, and when all jobs are finished, 
    this log file will be renamed to process_project_<projectname>.<starttime>.finished on success
    or to process_project_<projectname>.<starttime>.error when there was an error
    encountered. In this case, specific jobs that had errors can be found in the 
    logfile.
    More information on options for distribution options can be found in the
    [genomecomb joboptions help](https://derijkp.github.io/genomecomb/joboptions.html)

Other options
-------------
Scywalker defaults to analysis of the 10x v3 protocol. 
The following settings influence how **barcodes and UMIs** are found, and some can be used to analyse 10x v2

`-sc_whitelist`
    Used to provide a file with all possible correct barcodes. You should
    specify a different whitelist for 10x v2 using (e.g.) 
    `-sc_whitelist ~/bin/scywalker-0.108.0-linux-x86_64/whitelists/737K-august-2016.txt.gz`

`-sc_umisize`
    The default UMI size is 12 (v3). 10 v2 has a smaller UMI of 10; you can specify this using 
    `-sc_umisize 10`

`-sc_barcodesize`
    The default UMI size is 16  for both v2 and v3. This should normally not be changed

`-sc_adaptorseq`
    The default adapter sequence used to find the barcode and UMI is CTACACGACGCTCTTCCGATCT.
    This should normally not be changed

Some options influence **distribution** (besides `-d` and other [joboptions](https://derijkp.github.io/genomecomb/joboptions.html))

`-threads`
    The number of threads commands in jobs that support threading will use. On a cluster, this many cores
    will be reserved for the job. As threading in this commands often does not scale very well, keep the number
    typically low (4, max 8)

`-distrreg`
    determines how jobs are distributed over regions. Default is g5000000,
    which will distribute over regions of approximately 5M where splits can
    only occur in larger regions without known genes. Other possibilities are
    chr for per chromosome or 0 for no distribution over regions. It is advised
    to keep the default as larger regions also require more (peak) memory
    and 

`-maxfastqdistr`
    alignment is run per fastq; if there are very many small fastqs the overhead to 
    processes (alignment etc.) them separately (default) can become too large.
    The number given here limits the number of separately processed fastqs:
    if there are more separate input fastqs than the number given, they will be merged 
    before processing.

Various **other options** are

`-v`
    default 0, increase (up to 2) to increase the verbosity level, i.e. how much information 
    starting up jobs, dependencies, etc. is displayed

`-stack`
    set to 1 (default 0) to show an extended stack trace on error (mainly for debugging)

`-aligners`
    The default aligner is minimap2_splice (minimap2 with the splice preset). You could (experimentally)
    try to change this to minimap2_splicesmall to run with settings optimized for finding 
    small exons (but probably making more mistakes elsewhere)

Results
-------
After (successful) analysis, a sampledir contains various results files
following the genomecomb naming conventions (what-methods-samplename.extension). 
Result files are often tab separated files that are zstandard
(http://facebook.github.io/zstd/) compressed (extension .zst)

Compressed tsv files can also be easily read in e.g. R using
```
data=read.table(pipe("zstdcat file.tsv.zst"), sep="\t",header=T)
```
or if zstdcat is not installed, the included genomecomb command for this can be used:
```
data=read.table(pipe("cg zcat file.tsv.zst"), sep="\t",header=T)
```

Result files also usualy have an accompanying file with the .analysisinfo
extension that lists the tools (and their versions) used to generate the
file. 

### Sample Results
The most important result files found in the sample directories are:

#### sc_gene_counts_filtered-isoquant_sc-sminimap2_splice-sample1.tsv.zst
zstandard compressed, tab separated file with gene information and UMI counts per
gene/per cell. There are several counts for different ways of correcting
multimapping reads:
- `count`: UMI count per gene, multimapping reads are weighed (if maps to N genes -> each gets 1/N count)
- `nicount`: same, but intronic reads are not counted towards the gene they are in (as in count)
- `maxcount`: UMI count per gene, multimappers count 1 to each gene
- `uniquecount`: multimappers are not counted

The `cell` field indicates which cell the counts apply to. In this file (_filtered) 
only information on the emptydrops approved cells is given, the file
with _raw provides counts for all detected "cells".

The basic data (only the counts field, less info on the genes) is also supplied 
in the 10x (MEX) format in the directory
sc_gene_counts_filtered-isoquant_sc-sminimap2_splice-sample1.10x

#### sc_isoform_counts_filtered-isoquant_sc-sminimap2_splice-sample1.tsv.zst
zstandard compressed, tab separated file with isoform/transcript information and UMI counts per
isoform/per cell.
Transcripts are described in fields using the genePred convention as described 
in the [genomecomb gene/transcipt format](https://derijkp.github.io/genomecomb/format_gene.html).  
The **cell** field again indicates which cell the counts apply to. In this file (_filtered) 
only information on the emptydrops approved cells is given, the file
with _raw provides counts for all detected "cells".

This file also provides different ways counting/correcting for
reads supporting multiple isoforms (many reads are incomplete and could be derived from multiple isoforms):
- `counts_weighed`: reads supporting multiple (N) transcripts are weighed as 1/N
- `counts_unique`: count only reads uniquely supporting this one transcript
- `counts_strict`: only unique reads that cover >= 90% of the transcript
- `counts_aweighed`, `counts_aunique`, `counts_astrict`: same as above, but only reads with polyA (detected) counted

The basic data (using the counts_weighed field) is also supplied in the
10x (MEX) format in the directory as
`sc_isoform_counts_filtered-isoquant_sc-sminimap2_splice-sample1.weighed_count.10x`

#### sc_cellinfo_raw-isoquant_sc-sminimap2_splice-sample1.tsv.zst
A zstd compressed tsv file containing information on the detected cells
(one line for each cell), wth the following main fields:
- `cell` : cell barcode
- `readcount` : nr of reads assigned to to this cell
- `umicount` : nr of UMIs assigned to to this cell
- `is_cell` : 1 if emptydrops categorized the cell as real, 0 if not (the sc_cellinfo_filtered contains only info on cells where is_cell is 1)
- `nCount_RNA` : Seurat UMI count (can be lower than umicount due to filtering))
- `nFeature_RNA` : number of genes/features detected for the cell
Other information from Seurat, and two doublet finders is also added.

#### sc_group-scsorter-isoquant_sc-sminimap2_splice-sample1.tsv
tab separated file assigning cells to specific groups, in this case 
the celltype as determined by scsorter. More than one groupfile may be
present, e.g. based on sctype analysis (sc_group-sctype-isoquant_sc-sminimap2_splice-sample1.tsv)
Main fields are
- `cell` : cell barcode
- `group` : assigned group/cell type
- `group_filtered` : assigned group/cell type filtering put uncertains (foor tools that support this)
- `score` : asigned by celltyper
- `ncells` : number of cells in the group
- `UMAP_1` : UMAP coordinate 1 (for display)
- `UMAP_2` : UMAP coordinate 2 (for display)

#### pb_gene_counts-scsorter-isoquant_sc-sminimap2_splice-sample1.tsv.zst
zstandard compressed, tab separated pseudobulk gene counts file based on the scsorter
celltyper (there can more such files, e.g. one extra for the sctype celltyper).
This file contains the gene info similar to the sc_gene_counts file, and
provides the same types of counts, but the counts here are in wide format
(counts for all cell types on one line); The fields names indicate which
count and which celltype (and sample) is provided using the following
format: 
`<count_type>-<celltype>-scsorter-isoquant_sc-sminimap2_splice-<sample>`
e.g. `count-A549-scsorter-isoquant_sc-sminimap2_splice-mix1`

#### pb_isoform_counts-scsorter-isoquant_sc-sminimap2_splice-sample1.tsv.zst
zstandard compressed, tab separated pseudobulk isoformq counts file based on the scsorter
celltyper (there can more such files, e.g. one extra for the sctype celltyper).
This file contains the gene info similar to the sc_gene_counts file, and
provides the same types of counts, but the counts here are in wide format
(counts for all cell types on one line); The fields names indicate which
count and which celltype (and sample) is provided using the following
format: 
`<count_type>-<celltype>-scsorter-isoquant_sc-sminimap2_splice-<sample>`
e.g. `count-A549-scsorter-isoquant_sc-sminimap2_splice-mix1`

#### report_scywalker-sample1.html
A html reports summarizing the results, and including various informative graphs.

#### read_assignments-isoquant_sc-sminimap2_splice-sample1.tsv.zst
Gives a line for each assignment of a read to an isoform: Each line
contains information about the read, the alignment, the supported isoform,
etc. (in case of assignments to novel isoforms, the events/differences
are with regards to the closest known isoform, not the novel one)

#### map-sminimap2_splice-sample1.bam
bam file created by aligning the reads of sample1. The read names have
embedded cellbarcode and UMI information in both the read name and
comments

### Multi-sample results

The results covering multiple samples are found in the `compar` subdirectory of the projectdir.

#### pb_isoform_counts-scsorter-project.tsv.zst
A compressed tsv file containing the combination of all per sample
pseudobulk isoform count files. It has the same format as the individual
pb_isoform_count files (with fields like
`<count_type>-<celltype>-scsorter-isoquant_sc-sminimap2_splice-<sample>`),
but is wider (having more than one sample).

Novel transcripts often differ slightly on the ends between different
samples (depending on reads present). Such transcripts are matched over
samples (if the have the same junctions) keeping the most outer ends for
the combination.

#### pb_gene_counts-scsorter-project.tsv.zst
A compressed tsv file containing the combination of all per sample pseudobulk gene count files.
It has the same format as the individual pb_gene_count files (with fields 
like `<count_type>-<celltype>-scsorter-isoquant_sc-sminimap2_splice-<sample>`), but is wider
(having more than one sample). 

Novel genes with slightly different ends are similarly matched as isoforms.

Tools
-----
The version of genomecomb included in the scywalker distribution
provides many tools useful for analysis of scywalker results.
You can call these using `cg toolname ...` or `sw toolname ...` if you want
to specifically use the scywalker version. You can get an overview
of all tools in genomecomb using
```
cg help
```
and help on specific tools using
```
cg toolname -h
```
Following tools are typically useful for scywalker analysis:

```
cg viz_transcripts ?options? isoform_counts_file gene output_file
```
[viz_transcripts](https://derijkp.github.io/genomecomb/cg_viz_transcripts.html)
can be used to create a visual presentation of isoform usage of a given gene. get more help

```
cg sc_pseudobulk scgenefile scisoformfile groupfile
```
[sc_pseudobulk](https://derijkp.github.io/genomecomb/cg_sc_pseudobulk.html)
make pseudobulk files of sc_gene and sc_transcript files based on an sc_group file

```
cg multitranscript ?options? multitranscriptfile transcriptfile transcriptfile ?transcriptfile? ...
```
[multitranscript](https://derijkp.github.io/genomecomb/cg_multitranscript.html)
can be used to combine separate per sample transcript files in one multisample transcript file

```
cg viz file.tsv
```
The [viz](https://derijkp.github.io/genomecomb/cg_viz.html) tool
can be used to browse and query (compressed) tsv result files using a graphical interface

```
cg select ?options? ?datafile? ?outfile?
```
Using the [select](https://derijkp.github.io/genomecomb/cg_select.html) tool,
you can query the tsv result files on the command-line

```
cg zcat file ...
```
[zcat](https://derijkp.github.io/genomecomb/cg_zcat.html)
will concatenate (potentially compressed) files to standard output. It differs from
normal zcat in that it supports multiple compression types (based on the
file extension) including zstandard (.gz .zst .lz4 .rz .bz2)

License
-------
The use of this application is governed by the GPL (license.txt).

How to contact me
-----------------

Peter De Rijk
VIB - UAntwerp Center for Molecular Neurology, Neuromics Support Facility - Bioinformatics
University of Antwerp
Universiteitsplein 1
B-2610 Antwerpen, Belgium

tel.: +32-03-265.10.40
E-mail: Peter.DeRijk@uantwerpen.vib.be

