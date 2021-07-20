# Activate virtual environment
#. .venv/bin/activate
source activate myenv

snakemake --rerun-incomplete --latency-wait 1200 --cluster-config cluster.yaml --cluster "qsub -cwd -pe smp {cluster.cores} -l h_vmem={cluster.vmem}" --jobs 32


