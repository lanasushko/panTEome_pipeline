
# submit cluster in the background
nohup snakemake --rerun-incomplete --cluster "qsub -cwd -w e -N {cluster.job_name} -l h_rt={cluster.time} -l h_vmem={cluster.mem} -pe {cluster.pe} {cluster.threads} -o {cluster.o} -e {cluster.e} -j {cluster.join} -M {cluster.email}" --cluster-config cluster_config.yml --jobs 100 > snakemake.log 2>&1 &

snakemake --rerun-incomplete --cluster "qsub -cwd -w e -N {cluster.job_name} -l h_rt={cluster.time} -l h_vmem={cluster.mem} -pe {cluster.pe} {cluster.threads} -o {cluster.o} -e {cluster.e} -j {cluster.join} -M {cluster.email}" --cluster-config cluster_config.yml --jobs 100
# current version with conda
snakemake --use-conda --latency-wait 20 --rerun-incomplete --cluster "qsub -cwd -w e -N {cluster.job_name} -l h_rt={cluster.time} -l h_vmem={cluster.mem} -pe {cluster.pe} {cluster.threads} -o {cluster.o} -e {cluster.e} -j {cluster.join} -M {cluster.email}" --cluster-config cluster_config.yml --jobs 100

# dry run
snakemake --dry-run --rerun-incomplete --cluster "qsub -cwd -w e -N {cluster.job_name} -l h_rt={cluster.time} -l h_vmem={cluster.mem} -pe {cluster.pe} {cluster.threads} -o {cluster.o} -e {cluster.e} -j {cluster.join} -M {cluster.email}" --cluster-config cluster_config.yml
snakemake --dry-run --cluster "qsub -cwd -w e -N {cluster.job_name} -l h_rt={cluster.time} -l h_vmem={cluster.mem} -pe {cluster.pe} {cluster.threads} -o {cluster.o} -e {cluster.e} -j {cluster.join} -M {cluster.email}" --cluster-config cluster_config.yml


# kill snakemake
ps aux | grep snakemake  # Find the process ID (PID)
kill -9 <PID>  # Forcefully kill the Snakemake process
qdel -u <username> # delete all user's jobs

# run on chimi
