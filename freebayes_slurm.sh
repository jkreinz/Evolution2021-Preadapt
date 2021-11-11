#!/bin/bash
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=80
#SBATCH --time=12:00:00
#SBATCH --job-name freebayes_parallel

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=1

module load gnu-parallel/20180322
module load python/2.7.14-anaconda5.1.0

#run on multiple nodes
HOSTS=$(scontrol show hostnames $SLURM_NODELIST | tr '\n' ,)

#generate 100KB chunks across the genome
regionsfile=$($SCRATCH/software/freebayes/scripts/fasta_generate_regions.py $SCRATCH/tuberculatus_final.fasta 100000)


# iterate over regions using gnu parallel to dispatch jobs
cat "$regionsfile" | parallel --joblog slurm-$SLURM_JOBID.log -k -j $SLURM_NTASKS_PER_NODE -S $HOSTS --wd $PWD "/scratch/w/wrighste/kreinerj/software/freebayes/bin/freebayes -f $SCRATCH/tuberculatus_final.fasta --use-best-n-alleles 4  --max-complex-gap 1 --haplotype-length 1 --report-monomorphic /scratch/w/wrighste/kreinerj/final_bams/*.bam > /scratch/w/wrighste/kreinerj/regional_vcfs/{}.vcf " --region {}

#join seperately called regions, taking only first header and removing any duplicates
# cat *vcf| python2.7 /scratch/w/wrighste/kreinerj/software/freebayes/vcflib/scripts/vcffirstheader \
#    | /scratch/w/wrighste/kreinerj/software/freebayes/vcflib/bin/vcfstreamsort -w 1000 | /scratch/w/wrighste/kreinerj/software/freebayes/vcflib/bin/vcfuniq >> allsites.vcf #
