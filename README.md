# DpnII-seq
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
## DAG
<img src="https://github.com/tborrman/DpnII-seq/blob/master/dag.svg" alt="dag" width=1000px>
