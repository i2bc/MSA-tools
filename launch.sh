#! /bin/sh
#PBS -N generate_msa
#PBS -l mem=8gb
#PBS -l ncpus=1
#PBS -l walltime=02:00:00

#source /home/hugo.pointier/miniconda3/bin/activate orenza_env
cd /data/work/I2BC/hugo.pointier/test_mmseqs
structure=structure.pdb
output=output.fasta

python3 get_sequence_from_structure.py -i $structure -o $output

qsub -v "OUTDIR=/data/work/I2BC/hugo.pointier/test_mmseqs/msa,LOGDIR=/data/work/I2BC/hugo.pointier/test_mmseqs/log,INPUT=/data/work/I2BC/hugo.pointier/test_mmseqs/output.fasta,PROGRESSLOG=/data/work/I2BC/hugo.pointier/test_mmseqs/log/progress.log" template_mmseqs.sh
