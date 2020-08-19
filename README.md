# DpnII-seq
![snakemake-run](https://img.shields.io/badge/snakemake--run-passing-brightgreen)

DpnII-seq is an experimental assay to measure the genome wide cutting frequency of the restriction enzyme DpnII. 
Originally DpnII-seq was developed alongside [Liquid Chromatin Hi-C](https://www.biorxiv.org/content/10.1101/704957v1)
to control for biases in cutting frequency when estimating chromatin contact stability in K562 cells. 
However, the DpnII-seq assay and analysis pipeline could be modified for analyzing cutting frequency for variable restriction enzymes in 
variable cell types. 

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
##### Set up Channels
```bash
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels r
```
##### Set up Environment
```bash
conda env create --file envs/DpnII-seq_env.yaml
source activate DpnII-seq
```
##### Snake Command
```bash
snakemake -j 10 --latency-wait 60 --cluster-config cluster.json --cluster "bsub -q {cluster.queue} -W {cluster.time} -R {cluster.memory} -n {cluster.cores} -o {cluster.output} -e {cluster.error}" -p
```

