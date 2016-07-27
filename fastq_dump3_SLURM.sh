#!/bin/bash
#SBATCH --job-name=fastq-dump # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=END,FAIL           # notifications for job done & fail
#SBATCH --mail-user=timothy.h.webster@asu.edu # send-to address
#SBATCH --mem=4000
#SBATCH -t 96:00:00

cd /scratch/thwebste/Anole/Rupp_etal_Anole_DosageCompensation

source activate anole_dosage



for i in SRR1502174 SRR1502175 SRR1502176 SRR1502177 SRR1502178; do fastq-dump --gzip --outdir fastqs/ --readids --split-files $i; done


