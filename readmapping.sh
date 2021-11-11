##ALIGNMENT
# Software requirements:
#       AdapterRemoval: https://github.com/mikkelschubert/adapterremoval
#       bwa: https://github.com/lh3/bwa
#       samtools: https://www.htslib.org/download/
#       picard: https://github.com/broadinstitute/picard/releases/
 #
# USAGE:
# bash script.sh <reference_genome.fasta> <n_threads> <sample name>
#
# You can redirect all stderr to /dev/null if needed. bash script.sh ...options... 2> /dev/null
# parallelize like so: `cat listofsamples | parallel -j 10 "bash script.sh ref threads {}"

# Change paths accordingly
#bwa=/path/to/bwa
#samtools=/path/to/samtools
fastp=/ohta/julia.kreiner/software/fastp
picard=/ohta/julia.kreiner/software/picard.jar
pathtobams=/ohta/julia.kreiner/waterhemp/commongarden/fastqs

reference=$1
threads=$2
prefx=$3

# Assuming names are indivexmple1_R1.fastq.gz, indivexmple1.R2.fastq.gz
# It keeps whatever is before the first dot '.' as basename for the downstream outputs. You can change if needed

#cd $path/raw_data/$prefx

#cat ${prefx}*_1.fq.gz > ${prefx}.R1.fastq.gz
#cat ${prefx}*_2.fq.gz > ${prefx}.R2.fastq.gz

# 1. Remove adapters, polyQ tails
$fastp --in1 ${prefx}.R1.fastq.gz --in2 ${prefx}.R2.fastq.gz --out1 ${prefx}.R1.unmerged.fastq.gz --out2 ${prefx}.R2.unmerged.fastq.gz

# 2. Map reads to Reference Genome
bwa mem -t $threads -R "@RG\tID:$prefx\tSM:$prefx" $reference ${pathtobams}/${prefx}.R1.unmerged.fastq.gz ${pathtobams}/${prefx}.R1.unmerged.fastq.gz | samtools view -@ $threads -Sbh - > /ohta/julia.kreiner/waterhemp/commongarden/bams/secondhalf/${prefx}.uns.bam

cd $path/bams

#total=$(samtools view -@ $threads -c $prefx.full.uns.bam)
#echo -e "TotalReads\n$total" >> $prefx.log

# 3. Sort bam
sambamba sort -m 15GB --tmpdir $path/bams/tmp -t $threads -o $prefx.sorted.bam ${prefx}.uns.bam
rm ${prefx}.uns.bam
#mapped=$(samtools view -@ $threads -c $prefx.bam
#echo -e "MappedReads\n$mapped" >> ${prefx}.log
#echo "EndogenousDNA" >> ${prefx}.log
#python -c "print(float($mapped)/ $total)" >> ${prefx}.log


# 4. Mark duplicates
java -Xmx5G -Djava.io.tmpdir=~/tmp -jar $picard MarkDuplicates I=${prefx}.sorted.bam O=${prefx}.dd.bam M=${prefx}.dedup.log AS=true # Change here java memory usage with the options -Xmx. Here, I am assuming 10G of free memory
samtools index ${prefx}.dd.bam
#samtools flagstat ${prefx}.dd.bam > ${prefx}.stat
