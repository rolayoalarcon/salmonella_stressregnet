#!/bin/bash

#------- Descripción del trabajo -------

#SBATCH --job-name='dgr_modeling'
#SBATCH --comment='dgr_modeling'

#------- Parametrización -------

#SBATCH --requeue
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=80
#SBATCH --partition=bqhn01
#SBATCH --time=10-0:0:0

#------- Entrada/Salida -------

#SBATCH --output="/home/roberto_olayo/salmonella_stressregnet/workflow/logdir/out/%x-%j.out"
#SBATCH --error="/home/roberto_olayo/salmonella_stressregnet/workflow/logdir/err/%x-%j.err"

#------- Comando -------

/opt/bayresq.net/R/R-4.2.0/lib/R/bin/Rscript 01-growth_analysis_dgrowthr.R
