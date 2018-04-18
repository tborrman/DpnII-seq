# DpnII-seq
## Requirements
Install Snakemake via Miniconda [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
##### Set up Channels
```bash
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```
##### Set up Environment
```bash
conda env create --file envs/DpnII-seq_env.yaml
source activate DpnII-seq
```
