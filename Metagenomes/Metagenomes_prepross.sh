#							Initial QC							  #
#-----------------------------------------------------------------#

mkdir FastQC_initial
cd FastQC_initial

find *.fq.gz > list
sed 's/.fq.gz//g' list > list2
uniq list2 > sample_list
rm -f list*
sample_list=$(cat sample_list)
echo ${sample_list[@]}


#module load java/1.8.0  fastqc/0.11.8
for a in $sample_list
  do
    fastqc -o FastQC_initial/ "$a".fq.gz -t 16
done  


#-----------------------------------------------------------------#
#			          	Pre-processing, Adaptor removal		       		  #
#-----------------------------------------------------------------#
### Trim adapters
mkdir 1-trimmed
cd 1-trimmed

find *.fq.gz > list
sed 's/_..fq.gz//g' list > list2
uniq list2 > sample_list2
rm -f list*
sample_list=$(cat sample_list2)
echo ${sample_list2[@]}

for a in $sample_list2
  do
    AdapterRemoval --file1 "$a"_1.fq.gz --file2 "$a"_2.fq.gz --basename "$a" --output1 1-trimmed/"$a"_filtered_1.fq.gz --output2 1-trimmed/"$a"_filtered_2.fq.gz --adapter1 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT --adapter2 GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG --minlength 50 --threads 10 --gzip
done


#-----------------------------------------------------------------#
#														Remove duplicates										  #
#-----------------------------------------------------------------#

mkdir 0-raw 
mkdir 1-filtered

mv *_filtered.fq.gz 1-filtered
mv *.fq.gz 0-raw

find *.fq.gz > temp
sed 's/_filtered_[1-2].fq.gz//g' temp > temp2
uniq temp2 > sample_list.txt
rm -f temp*
sample_list=$(cat sample_list.txt)


for a in $sample_list
do
cat "$a"_filtered_1.fq.gz| seqkit rmdup -j 20 -s -o "$a"_1_clean.fastq
cat "$a"_filtered_2.fq.gz| seqkit rmdup -j 20 -s -o "$a"_2_clean.fastq
done

cat D00628_filtered_1.fq.gz | seqkit rmdup -j 20 -s -o D00628_1_clean.fastq
cat D00628_filtered_2.fq.gz | seqkit rmdup -j 20 -s -o D00628_2_clean.fastq


#-----------------------------------------------------------------#
#		         		Re-pair fastq after rm duplicates    	         	  #
#-----------------------------------------------------------------#
### jre/1.8.0 bbmap/38.35
#Used to fix files of paired reads that has become disordered after dereplication

find *.fastq > temp
sed 's/_[1-2]_clean.fastq//g' temp > temp2
uniq temp2 > sample_list.txt
rm -f temp*
sample_list=$(cat sample_list.txt)

for a in $sample_list
do
repair.sh in="$a"_1_clean.fastq in2="$a"_2_clean.fastq out=Repaired/"$a"_1_repaired.fq out2=Repaired/"$a"_2_repaired.fq outsingle=Singletons/"$a"_singletons.fq 
rm "$a"_1_clean.fastq
rm "$a"_2_clean.fastq
done

#-----------------------------------------------------------------#
#		         		Filtering of Host DNA, using BWA	   		      	  #
#-----------------------------------------------------------------#

#	Filter against:
#	Human       HG19
#   Bryozoan GCA_914767715.1

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz -P host_ref/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/914/767/715/GCA_914767715.1_tzMemMemb1.1/GCA_914767715.1_tzMemMemb1.1_genomic.fna.gz -P host_ref/


mkdir 1-Host_Removal
cd 1-Host_removal

# 1) Build minimap Host_DB
module load samtools/1.9 bedtools/2.28.0 minimap2/2.6
minimap2 -d ref.mmi $REF_DIR/*

### Human Filtering with minimap
REF_DIR='/home/perbec/PChem0002/Jyllinghavn_metagenomes2021/host_ref'
FASTA='/home/perbec/PChem0002/Jyllinghavn_metagenomes2021/nohost'
OUT='/home/perbec/PChem0002/Jyllinghavn_metagenomes2021/nohost/Repaired'

cd $FASTA

cd $FASTA/
find *.fq > temp
sed 's/_[1-2]_repaired.fq//g' temp > temp2
uniq temp2 > sample_list.txt
rm -f temp*
sample_list=$(cat sample_list.txt)

#human filtering
Ref='human.fna'
for stuff in $sample_list
  do
    minimap2 -ax sr $REF_DIR/$Ref $FASTA/"$stuff"_1_repaired.fq $FASTA/"$stuff"_2_repaired.fq > $OUT/"$stuff"_mapped_and_unmapped.sam
      samtools view --threads 14 -bS $OUT/"$stuff"_mapped_and_unmapped.sam > $OUT/"$stuff"_mapped_and_unmapped.bam
        rm $OUT/"$stuff"_mapped_and_unmapped.sam
          samtools view --threads 14 -b -f 13 -F 1280 $OUT/"$stuff"_mapped_and_unmapped.bam > $OUT/"$stuff"_unmapped.bam
            rm $OUT/"$stuff"_mapped_and_unmapped.bam
              samtools sort --threads 14 -n $OUT/"$stuff"_unmapped.bam -o $OUT/"$stuff"_sorted.bam
                rm $OUT/"$stuff"_unmapped.bam
                  bedtools bamtofastq -i $OUT/"$stuff"_sorted.bam -fq $OUT/"$stuff"_nohuman_1.fastq -fq2 $OUT/"$stuff"_nohuman_2.fastq
                    pigz -p 20 $OUT/"$stuff"_nohuman_*.fastq                  
done


REF_DIR='/home/perbec/PChem0002/Jyllinghavn_metagenomes2021/host_ref'
FASTA='/home/perbec/PChem0002/Jyllinghavn_metagenomes2021/nohost'
OUT='/home/perbec/PChem0002/Jyllinghavn_metagenomes2021/nohost/nobryozoan'

#Bryozoan filtering
Ref='bryozoan.fna'
for stuff in $sample_list
  do
    minimap2 -ax sr $REF_DIR/$Ref $FASTA/"$stuff"_nohuman_1.fastq.gz $FASTA/"$stuff"_nohuman_2.fastq.gz > $OUT/"$stuff"_mapped_and_unmapped.sam
      samtools view --threads 16 -bS $OUT/"$stuff"_mapped_and_unmapped.sam > $OUT/"$stuff"_mapped_and_unmapped.bam
        rm $OUT/"$stuff"_mapped_and_unmapped.sam
          samtools view --threads 16 -b -f 13 -F 1280 $OUT/"$stuff"_mapped_and_unmapped.bam > $OUT/"$stuff"_unmapped.bam
            rm $OUT/"$stuff"_mapped_and_unmapped.bam
              samtools sort --threads 16 -n $OUT/"$stuff"_unmapped.bam -o $OUT/"$stuff"_sorted.bam
                rm $OUT/"$stuff"_unmapped.bam
                  bedtools bamtofastq -i $OUT/"$stuff"_sorted.bam -fq $OUT/"$stuff"_nohuman_nobryozoan_1.fastq -fq2 $OUT/"$stuff"_nohuman_nobryozoan_2.fastq
                    pigz -p 20 $OUT/"$stuff"_nohuman_nobryozoan_*.fastq                  
done

#-------------------------------------------------------------------#
#			               		Co-Assembly, using MegaHit			      			#
#-------------------------------------------------------------------#

# Run on computerome 2.0: 
# I'm using MegaHit with the 'meta-large' preset, because its is very complex environmental samples, if meta-sensitive is used then it takes forever!!! 
#Keep the minimum length at 1000, else we get very bad binning.
module load megahit/1.1.1 anaconda2/4.4.0
FASTA_DIR='<path/to/FASTQ/>'
WORK_DIR='<path/to/WD/>'
FASTA=$(ls $FASTA_DIR/*.fq | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])')
# Check the FASTA list
echo $FASTA
# run megahit
megahit -r $FASTA --min-contig-len 1000 -t 40 --presets meta-large -o $WORK_DIR


#-------------------------------------------------------------------#
#	  Move to Anvio pipeline for binning and curation			      			#
#-------------------------------------------------------------------#
