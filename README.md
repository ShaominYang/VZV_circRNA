# Core for "Identification and Characterization of Varicella Zoster Virus Circular RNA in Lytic Infection"
# Requirements
- 1.Bash (Ubuntu, version 18.04)
- 2.Perl (https://www.perl.org)
- 3.Java (https://javadl.oracle.com)
- 4.BWA (http://bio-bwa.sourceforge.net)
- 5.Bowtie2 (https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- 6.SAMtools (http://www.htslib.org/)
- 7. Trim Galore (http://www.bioinformatics.bbsrc.ac.uk/projects/trim_galore/)
- 8. Unicycler (https://github.com/rrwick/Unicycler)
- 9.CIRI2(version v2.0.6) (https://sourceforge.net/projects/ciri/files/CIRI2/)
- 10.CIRI-full (version 2.0) (https://sourceforge.net/projects/ciri/files/CIRI-full/)
- 11.vircircRNA (https://github.com/jiwoongbio/vircircRNA)
- 12.find_circ (https://github.com/marvin-jens/find_circ)
- 13.CIRI-long (https://github.com/bioinfo-biols/CIRI-long)
# genome 
- 14.VZV pOka_GFP_Luc strain (https://www.ncbi.nlm.nih.gov/nuccore/PP054841)

# Data availability
All sequencing data used in this study is available via NCBI SRA and Gene Expression Omnibus database, accession numbers PRJNA1056528, GSE223870, GSE252124 and GSE223957. 
1. For whole-genome sequencing
https://www.ncbi.nlm.nih.gov/sra/PRJNA1056528
2. For circRNA short reads sequencing
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223870
3. For circRNA long reads sequencing
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE252124
4. For mRNA sequencing
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223957
## To assemble genomic sequences from whole-genome sequencing data
```Shell
conda activate assembly
time unicycler -t 100 \
-1 trim_galore_out_dir/VZV_BAC_WT_short_1_val_1.fq.gz -2 trim_galore_out_dir/VZV_BAC_WT_short_2_val_2.fq.gz \
-o unicycler_VZV_BAC_WT

##Generating genome contained human and VZV

cat chrALL.fa PP054841.fasta > chrALL.fa
```
## 1.1  Build BWA index
```Shell
bwa index chrALL.fa
```
## 1.2  Running CIRI2 pipeline 
#### note: Please replace "/kydata/ysm" with your work path

```Shell
for i in mock_1 mock_2 mock_3 vzv_24h_1 vzv_24h_2 vzv_24h_3 vzv_48h_1 vzv_48h_2 vzv_48h_3
do
mkdir ${i}_ciri2_output
bwa mem -t 100 chrALL.fa ${i}_1.fq.gz ${i}_2.fq.gz > ${i}_ciri2_output/${i}.sam
perl /kydata/ysm/CIRI/CIRI_v2.0.6/CIRI2.pl -I ${i}_ciri2_output/${i}.sam -O ${i}_ciri2_output/${i}.ciri -F chrALL.fa -A mRNA.gtf -T 56
perl /kydata/ysm/CIRI/CIRI-AS/CIRI_AS_v1.2.pl -S ${i}_ciri2_output/${i}.sam -C ${i}_ciri2_output/${i}.ciri -F chrALL.fa -A mRNA.gtf -O ${i}_ciri2_output/${i} -D yes
java -jar /kydata/ysm/CIRI/CIRI-full_v2.0/CIRI-full.jar RO1 -1 ${i}_1.fastq -2 ${i}_2.fastq -o ${i}_ciri2_output/${i} 
bwa mem -t 42 chrALL.fa ${i}_ciri2_output/${i}_ro1.fq > ${i}_ciri2_output/${i}_ro1.sam
java -jar /kydata/ysm/CIRI/CIRI-full_v2.0/CIRI-full.jar RO2 -r chrALL.fa -s ${i}_ciri2_output/${i}_ro1.sam -l 300 -o ${i}_ciri2_output/${i}RO2
java -jar /kydata/ysm/CIRI/CIRI-full_v2.0/CIRI-full.jar Merge -c ${i}_ciri2_output/${i}.ciri -as ${i}_ciri2_output/${i}_jav.list -ro ${i}_ciri2_output/${i}RO2_ro2_info.list -r chrALL.fa -o ${i}_ciri2_output/${i}
unset DISPLAY
java -jar /kydata/ysm/CIRI/CIRI-full_v2.0/CIRI-vis.jar -o ${i}_ciri2_output/${i}_stdir -i ${i}_ciri2_output/${i}_merge_circRNA_detail.anno -l ${i}_ciri2_output/${i}_library_length.list -a mRNA.gtf -r chrALL.fa -min 1
done
```
## 1.3  Build bowtie2 index
```Shell
bowtie2-build chrALL.fa chrALL.fa
```
## 1.4  Running find_circ pipeline
```Shell
for i in mock_1 mock_2 mock_3 vzv_24h_1 vzv_24h_2 vzv_24h_3 vzv_48h_1 vzv_48h_2 vzv_48h_3
do
bowtie2 -p 150 --very-sensitive --score-min=C,-15,0 --mm -x chrALL.fa -q -1 ${i}_1.fq.gz -2 ${i}_2.fq.gz | /kydata/ysm/samtools-1.13/bin/samtools view -hbuS - | /kydata/ysm/samtools-1.13/bin/samtools sort -o ${i}_output.bam
/kydata/ysm/samtools-1.13/bin/samtools view -hf 4 ${i}_output.bam | /kydata/ysm/samtools-1.13/bin/samtools view -Sb - > unmapped.bam
python2 /kydata/ysm/find_circ-1.2/unmapped2anchors.py unmapped.bam | gzip > anchors.fq.gz
bowtie2 -p 150 --reorder --mm  --score-min=C,-15,0 -q -x chrALL.fa -U anchors.fq.gz | python2 /kydata/ysm/find_circ-1.2/find_circ.py --genome=chrALL.fa --prefix=${i}_ --name=${i} --stats=${i}_stats.txt --reads=${i}_splice_reads.fa > ${i}_spliced_sites.bed
grep CIRCULAR ${i}_spliced_sites.bed | grep -v chrM | gawk '$5>=2' | grep UNAMBIGUOUS_BP | grep ANCHOR_UNIQUE | python2 /kydata/ysm/find_circ-1.2/maxlength.py 100000 > ${i}_find_circ.candidates.bed
done
```
# 2.circRNA identification---> Pooling all samples together

```Shell
cat mock_1_1.fq.gz mock_2_1.fq.gz mock_3_1.fq.gz vzv_24h_1_1.fq.gz vzv_24h_2_1.fq.gz vzv_24h_3_1.fq.gz vzv_48h_1_1.fq.gz vzv_48h_2_1.fq.gz vzv_48h_3_1.fq.gz > all_1.fq.gz
cat mock_1_2.fq.gz mock_2_2.fq.gz mock_3_2.fq.gz vzv_24h_1_2.fq.gz vzv_24h_2_2.fq.gz vzv_24h_3_2.fq.gz vzv_48h_1_2.fq.gz vzv_48h_2_2.fq.gz vzv_48h_3_2.fq.gz > all_2.fq.gz
```
## 2.1 Running CIRI2 pipeline
```Shell
for i in all
do
mkdir ${i}_ciri2_output
bwa mem -t 100 chrALL.fa ${i}_1.fq.gz ${i}_2.fq.gz > ${i}_ciri2_output/${i}.sam
perl /kydata/ysm/CIRI/CIRI_v2.0.6/CIRI2.pl -I ${i}_ciri2_output/${i}.sam -O ${i}_ciri2_output/${i}.ciri -F chrALL.fa -A mRNA.gtf -T 56
perl /kydata/ysm/CIRI/CIRI-AS/CIRI_AS_v1.2.pl -S ${i}_ciri2_output/${i}.sam -C ${i}_ciri2_output/${i}.ciri -F chrALL.fa -A mRNA.gtf -O ${i}_ciri2_output/${i} -D yes
java -jar /kydata/ysm/CIRI/CIRI-full_v2.0/CIRI-full.jar RO1 -1 ${i}_1.fastq -2 ${i}_2.fastq -o ${i}_ciri2_output/${i} 
bwa mem -t 42 chrALL.fa ${i}_ciri2_output/${i}_ro1.fq > ${i}_ciri2_output/${i}_ro1.sam
java -jar /kydata/ysm/CIRI/CIRI-full_v2.0/CIRI-full.jar RO2 -r chrALL.fa -s ${i}_ciri2_output/${i}_ro1.sam -l 300 -o ${i}_ciri2_output/${i}RO2
java -jar /kydata/ysm/CIRI/CIRI-full_v2.0/CIRI-full.jar Merge -c ${i}_ciri2_output/${i}.ciri -as ${i}_ciri2_output/${i}_jav.list -ro ${i}_ciri2_output/${i}RO2_ro2_info.list -r chrALL.fa -o ${i}_ciri2_output/${i}
unset DISPLAY
java -jar /kydata/ysm/CIRI/CIRI-full_v2.0/CIRI-vis.jar -o ${i}_ciri2_output/${i}_stdir -i ${i}_ciri2_output/${i}_merge_circRNA_detail.anno -l ${i}_ciri2_output/${i}_library_length.list -a mRNA.gtf -r chrALL.fa -min 1
done
```
## 2.2 Running find_circ pipeline
```
for i in all
do
bowtie2 -p 150 --very-sensitive --score-min=C,-15,0 --mm -x chrALL.fa -q -1 ${i}_1.fq.gz -2 ${i}_2.fq.gz | /kydata/ysm/samtools-1.13/bin/samtools view -hbuS - | /kydata/ysm/samtools-1.13/bin/samtools sort -o ${i}_output.bam
/kydata/ysm/samtools-1.13/bin/samtools view -hf 4 ${i}_output.bam | /kydata/ysm/samtools-1.13/bin/samtools view -Sb - > unmapped.bam
python2 /kydata/ysm/find_circ-1.2/unmapped2anchors.py unmapped.bam | gzip > anchors.fq.gz
bowtie2 -p 150 --reorder --mm  --score-min=C,-15,0 -q -x chrALL.fa -U anchors.fq.gz | python2 /kydata/ysm/find_circ-1.2/find_circ.py --genome=chrALL.fa --prefix=${i}_ --name=${i} --stats=${i}_stats.txt --reads=${i}_splice_reads.fa > ${i}_spliced_sites.bed
grep CIRCULAR ${i}_spliced_sites.bed | grep -v chrM | gawk '$5>=2' | grep UNAMBIGUOUS_BP | grep ANCHOR_UNIQUE | python2 /kydata/ysm/find_circ-1.2/maxlength.py 100000 > ${i}_find_circ.candidates.bed
done
```
## 2.3  Running vircircRNA pipeline
```Shell
##Build BWA index
perl vircircRNA_chromosome.pl PP054841.fasta > concatenated_circular_PP054841.fasta
bwa index concatenated_circular_PP054841.fasta
bwa mem -t 100 -Y concatenated_circular_PP054841.fasta all_1.fq.gz all_2.fq.gz > mapped_read.sam
perl vircircRNA_junction.pl -g PP054841.gff3 -A junction.alignment.html mapped_read.sam concatenated_circular_PP054841.fasta > PP054841_junction.txt
perl vircircRNA_diagram.pl -g PP054841.gff3> PP054841_junction.txt concatenated_circular_PP054841.fasta PP054841.fasta > diagram.png
```
## 3 Running CIRI-long pipeline
```Shell
minimap2 -d chrALL.min chrALL.fa
#######call
CIRI-long call -i SY5Y-pokaA.fq.gz \
               -o SY5Y-pokaA_chrALL \
               -r chrALL.fa \
               -p SY5Y-pokaA_chrALL \
               -t 100
			   
#####collapse		   			   
CIRI-long collapse -i ./SY5Y-pokaA_chrALL.lst \
                   -o ./SY5Y-pokaA_chrALL_collpase \
                   -p SY5Y-pokaA_chrALL \
                   -r chrALL.fa \
                   -t 100
python3 misc/convert_bed.py SY5Y-pokaA_chrALL_collpase/SY5Y-pokaA_chrALL.info SY5Y-pokaA_chrALL_circ.bed	
```
# 4. VZV circVLTs mutagenesis RNA-seq
```Shell
##Build STAR index
STAR --runThreadN 42 --runMode genomeGenerate --genomeDir /kydata/ysm/index/star_index/chrALL_149 \
--genomeFastaFiles chrALL.fa \
--sjdbGTFfile mRNA.gtf \
--sjdbOverhang 149
```
## 4.1 Running RNA-seq gene expression quantification
  
```Shell
for i in pOka_WT_1 pOka_WT_2 pOka_WT_3 pOka_M1_1 pOka_M1_3 pOka_M1_3 pOka_M2_1 pOka_M2_2 pOka_M2_3
do
STAR --genomeDir /kydata/ysm/index/star_index/chrALL_149 \
--readFilesIn ${i}_1.fq.gz ${i}_2.fq.gz \
--readFilesCommand gunzip -c \
--twopassMode Basic \
--outSAMtype BAM Unsorted \
--chimSegmentMin 12 \
--chimJunctionOverhangMin 12 \
--alignSJDBoverhangMin 10 \
--alignMatesGapMax 100000 \
--alignIntronMax 100000 \
--chimSegmentReadGapMax 3 \
--alignSJstitchMismatchNmax 5 -1 5 5 \
--runThreadN 100 \
--outSAMstrandField intronMotif \
--chimOutJunctionFormat 1 \
--outFileNamePrefix $i
/kydata/ysm/samtools-1.13/bin/samtools sort ${i}Aligned.out.bam -o ${i}_sorted.bam -@ 40
/kydata/ysm/samtools-1.13/bin/samtools index ${i}_sorted.bam
echo ${i}_mapped
rm ${i}Aligned.out.bam
rm *.out
featureCounts -T 50 -p -t exon -g gene_id -a mRNA.gtf -o counts/${i}.counts.txt ${i}_sorted.bam
echo ${i}_Counts_finished
cut -f 1,7 counts/${i}.counts.txt |grep -v '^#' >feacturCounts/${i}_feacturCounts.txt
echo ${i}_featureCounts_finished
done 
```
## 5 Running coverage statistics pipeline 

### 5.1 PP054841_WT
```Shell
for i in pOka_WT_1 pOka_WT_2 pOka_WT_3
do
bwa mem -t 100 PP054841.fasta ${i}_1.fq.gz ${i}_2.fq.gz > ${i}.sam
#排序
/kydata/ysm/samtools-1.13/bin/samtools view -b ${i}.sam > ${i}.bam
/kydata/ysm/samtools-1.13/bin/samtools sort ${i}.bam -o ${i}_sorted.bam -@ 100
/kydata/ysm/samtools-1.13/bin/samtools index ${i}_sorted.bam
/kydata/ysm/samtools-1.13/bin/samtools depth -a -m 0 ${i}_sorted.bam > ${i}_coverage.txt
echo ${i}_mapped
done
```
### 5.2 PP054841_M1
```Shell
bwa index PP054841_M1.fasta
for i in pOka_M1_1 pOka_M1_3 pOka_M1_3
do
bwa mem -t 100 PP054841_M1.fasta ${i}_1.fq.gz ${i}_2.fq.gz > ${i}.sam
#排序
/kydata/ysm/samtools-1.13/bin/samtools view -b ${i}.sam > ${i}.bam
/kydata/ysm/samtools-1.13/bin/samtools sort ${i}.bam -o ${i}_sorted.bam -@ 100
/kydata/ysm/samtools-1.13/bin/samtools index ${i}_sorted.bam
/kydata/ysm/samtools-1.13/bin/samtools depth -a -m 0 ${i}_sorted.bam > ${i}_coverage.txt
echo ${i}_mapped
done
```
### 5.3 PP054841_M2
```Shell
bwa index PP054841_M2.fasta
for i in pOka_M2_1 pOka_M2_2 pOka_M2_3
do
bwa mem -t 100 PP054841_M2.fasta ${i}_1.fq.gz ${i}_2.fq.gz > ${i}.sam
#排序
/kydata/ysm/samtools-1.13/bin/samtools view -b ${i}.sam > ${i}.bam
/kydata/ysm/samtools-1.13/bin/samtools sort ${i}.bam -o ${i}_sorted.bam -@ 100
/kydata/ysm/samtools-1.13/bin/samtools index ${i}_sorted.bam
/kydata/ysm/samtools-1.13/bin/samtools depth -a -m 0 ${i}_sorted.bam > ${i}_coverage.txt
echo ${i}_mapped
done
```
