# DpnII-seq
![snakemake-run](https://img.shields.io/badge/snakemake--run-passing-brightgreen)

DpnII-seq is an experimental assay to measure the genome wide cutting frequency of the restriction enzyme DpnII. 
Originally DpnII-seq was developed alongside [Liquid Chromatin Hi-C](https://www.biorxiv.org/content/10.1101/704957v1)
to control for biases in cutting frequency when estimating chromatin contact stability in K562 cells. 
However, the DpnII-seq assay and analysis pipeline could be modified for analyzing cutting frequency for variable restriction enzymes in 
variable cell types. The DpnII-seq workflow performs mapping, filtering of artifacts, multiple resolution binning, copy number correction
and plotting of quality metrics.

This workflow has options for digesting with the following restriction enzymes:
- [DpnII](https://www.neb.com/products/r0543-dpnii#Product%20Information)
- [HindIII](https://www.neb.com/products/r0104-hindiii#Product%20Information)
- [FatI](https://www.neb.com/products/r0650-fati#Product%20Information)

As K562 cells have a primarily triploid karyotype with regions of variable copy number, 
the analysis workflow corrects coverage tracks to a diploid state genome wide. 
If the user is applying DpnII-seq to cells with variable copy number states 
we provide scripts to correct for this bias using a Gaussian mixture model approach.

## Experimental Protocol
<img src="https://github.com/tborrman/DpnII-seq/blob/master/img/dpnII-seq.PNG" alt="dpn" width=300px>

## Computational Workflow
<img src="https://github.com/tborrman/DpnII-seq/blob/master/img/dag.svg" alt="dag" width=500px>

## Requirements
Install Snakemake via Miniconda [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
##### Set up channels
```bash
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels r
```
##### Set up environment
```bash
conda env create --file envs/DpnII-seq_env.yaml
source activate DpnII-seq
```
##### Set up data directory
```bash
.
├── data
│   ├── binning
│   ├── copy_number
│   ├── dpnII_sites
│   ├── fastq
│   ├── fatI_sites
│   ├── hindIII_sites
│   └── indexes
├── DpnII-seq
```
DpnII-seq workflow requires a data directory with the following files:
- binning/
  - bedfiles denoting binned genome (ex. hg19_40kb.bed)
- copy_number/ 
  - bedfiles denoting binned genome with copy number state (4th column)
- \[X\]\_sites/
  - bedfiles denoting restriction sites of restriction enzyme X for each chromosome (ex. dpnII_sites_chr1.bed)
- fastq/
  - paired end fastq files (\*\_R1.fastq, \*\_R2.fastq)
- indexes/
  - bowtie-build index files for reference genome (\*.ebwt)


## Usage
Currently the DpnII-seq workflow has only been tested in an LSF cluster environment
##### Snake Command
```bash
snakemake -j 10 --latency-wait 60 --cluster-config cluster.json --cluster "bsub -q {cluster.queue} -W {cluster.time} -R {cluster.memory} -n {cluster.cores} -o {cluster.output} -e {cluster.error}" -p
```
## Contributing to DpnII-seq
All contributions, bug reports, bug fixes, documentation improvements, enhancements, and ideas are welcome!

## References
Houda Belaghzal*, Tyler Borrman*, Andrew D. Stephens, Denis L. Lafontaine, Sergey V. Venev, Zhiping Weng, John F. Marko, Job Dekker. (2019). [Compartment-dependent chromatin interaction dynamics revealed by liquid chromatin Hi-C](https://www.biorxiv.org/content/10.1101/704957v1). *bioRxiv*

## License 
[GNU General Public License v3.0](LICENSE)
