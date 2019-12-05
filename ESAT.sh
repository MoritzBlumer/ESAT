#!/bin/bash

#### Set paths paths:
MOBSTERPATH= 			# set path to Mobster .jar file
MELTPATH= 			# path to MELT .jar file
SIMSEQPATH= 			# path to SimSeq .jar file
SAMTOFASTQPATH=			# path to SamToFastq .jar file
MARKDUPLICATESPATH=		# path to MarkDuplicates .jar file
COLLECTINSERTSIZEMETRICSSPATH=	# path to CollectInsertSizeMetrics .jar file (Picardtools)

#### Usage:
## bash ESAT.sh GenomeSequencePath TESequencePath OutputName Coverage NoOfInsertions RunID

#### process input
## input parameters:
Sequence=$1  #original sequence file
Insertion=$2 #insertion file
Name=$3      #output name
convertNumber_of_insertions=$5
ID=$6

#### Simulate reference sequence:
python3.6 scripts/seq_sim.py $Sequence $Insertion $Name $ID --n $Number_of_insertions 
bwa index $Sequence
bedtools sort -i $ID.insertions.bed > $ID.insertions.sorted.bed
rm $ID.insertions.bed

#### Determine read number:
read_number=$(python3.6 scripts/read_number.py "${Name}.fa" $coverage)

#### Run SimSeq:
java -jar -Xmx4g $SIMSEQPATH -1 90 -2 90 --adapter1 "" --adapter2 "" --error input/error_profile.txt --insert_size 500 --read_number $read_number --read_prefix Read_simulation --reference "${Name}.fa" --out $Name.sam
## convert SAM to FASTQ:
java -jar -Xmx4g $SAMTOFASTQPATH INPUT=$Name.sam FASTQ=library.1.fastq SECOND_END_FASTQ=library.2.fastq INCLUDE_NON_PF_READS=true VALIDATION_STRINGENCY=SILENT 

#### Map reads:
## set variables:
READS1=library.1.fastq
READS2=library.2.fastq
SAMPLE=$Name
RG="@RG\tID:GRE01\tSM:GRE01\tPL:ILLUMINA\tLB:EROB117:" 
## map processed reads:
bwa mem -t 8 -R $RG $Sequence $READS1 $READS2 > ${Name}.sam  &&
## convert SAM to BAM:
samtools view -T $Sequence -Sb ${Name}.sam > ${Name}.bam &&

#### Proces BAM files:
## sort:
samtools sort -@ 8 -m 5G ${Name}.bam -o ${Name}.sorted.bam &&
## mark duplicates:
echo " INPUT=$SAMPLE OUTPUT=${SAMPLE/.bam/.mkdup.bam} METRICS_FILE=${SAMPLE/.bam/.mkdup_metrics} "
java -Xmx2g -jar $MARKDUPLICATESPATH INPUT=${Name}.sorted.bam OUTPUT=${Name}.sorted.mkdup.bam METRICS_FILE=${Name}.sorted.mkdup_metrics  &&
## generate flag statistics:
samtools flagstat  ${Name}.sorted.mkdup.bam &&
## index:
samtools index ${Name}.sorted.mkdup.bam 

#### Run Mobster:
## generate Mobster "properties"-file
printf "#####
#Example of Mobster properties file for the Master jar file (Mobster.jar)
#####
#Whether or not to use Picard to estimate the Insert Size metrics.
#This aids the determination of prediction borders especially for events with only 
#support from one side.
USE_PICARD=true
PICARD_COLLECT_INSERT_SIZE_METRICS_JAR=$COLLECTINSERTSIZESPATH
#Below properties for when you do not use picard.
#It is not necessary to provide when you do not use picard, but it is recommended.
#Now commented out because we use picard.
#MEAN_FRAGMENT_LENGTH=311
#SD_FRAGMENT_LENGTH=13
#LENGTH_99PROCENT_OF_FRAGMENTS=384
#The BAM file in which you want to predict MEIs
#For now only Mosaik and BWA bam files are accepted.
#In the near future all .bam files will be accepted.
IN_FILE=${Name}.sorted.mkdup.bam
#The output prefix
OUT_FILE=$Name.mobster
##Properties about the read length of the original bam file. You are not obligated
#to provide this, but it does improve the prediction windows of single cluster predictions
#For now commented out, as we did not use this setting originally on the simulation bam file.
USE_READ_LENGTH=true
READ_LENGTH=90
#Do you want to do multiple sample calling on your BAM file? Possible when there are valid @RG tags in your bam file each with a sample name field SM:
#Multisample calling?
MULTIPLE_SAMPLE_CALLING=false
#Stringent: at least one sample needs to have >= MINIMUM_SUPPORTING_READS
MULTIPLE_SAMPLE_CALLING_STRINGENT=false
#The mapping tool which was used for creating the bam file supplied in IN_FILE
#This value can either be bwa or mosaik
#MAPPING_TOOL=bwa
# ADDITIONS FROM MOBSTER_LATEST
MAPPING_TOOL=unspecified
#MINIMUM_MAPQ_ANCHOR will only be used when MAPPING_TOOL is set to unspecified.
#It specifies at what mapq a read can be considered as an anchor
MINIMUM_MAPQ_ANCHOR=20
#Whether mobster should try and search for clipped reads
USE_SPLIT=true
#When searching for clipped reads, how long should the unmapped sequence at least be (bp):
MINIMUM_CLIP_LENGTH=50
#If a clipped read is found with MINIMUM_CLIP_LENGTH, how long may the other end of the same read
#be MAXIMUM clipped
MAXIMUM_CLIP_LENGTH=7
#If a clipped read is found with MINIMUM_CLIP_LENGTH, what should the minimum
#average base quality be of that read.
#we use such a low threshold now, because the simulation bam file contains low quality bases
MINIMUM_AVG_QUALITY=10
#For the detection of poly A tails in clipped reads.
#How long should the Poly A tail be and how many times may
#there be a different nucleotide than A for polyA stretch
#or different from T for a polyT stretch
MINIMUM_POLYA_LENGTH=9
MAXIMUM_MISMATCHES_POLYA=1
#Number of records which are held in memory.
#Increase this number if you have more RAM available (this will reduce I/O and thus increase speed)
MAX_RECORDS_IN_RAM=1000000
#The command to do the mobiome mapping
#In this example we use Mosaik to do mobiome mapping
#This can also be done with bwa if you have the mobiome reference
#For bwa mobiome mapping make sure the whole command fits one line (use pipes | or && to achieve this)
#Where bwa requires the input file (fastq) and output file:
#Do not use hardcoded paths but replace with (FASTQ) and (OUT_FILE) just like in the line below.
MOBIOME_MAPPING_CMD=MosaikBuild -q (FASTQ) -st illumina -out (DAT_FILE) -quiet && MosaikAligner -in (DAT_FILE) -out (OUT_FILE) -ia data/mobster_library/SINE-tRNA/CHR2B_Ttr.ME_CAT.SINE2-tRNA.dat -hs 9 -mmp 0.1 -act 20 -j data/mobster_library/SINE-tRNA/CHR2B_Ttr.ME_CAT.SINE2-tRNA -p 2 -annpe input/2.1.26.pe.100.0065.ann -annse input/2.1.26.se.100.005.ann -quiet
#Name of the sample. This will be added as column to the prediction file
SAMPLENAME=Uma_simulation
#The tool which was used for mobiome mapping. Can be either mosaik or bwa
MOBIOME_MAPPING_TOOL=mosaik
#Were there paired-end reads in the BAM file
PAIRED_END=true
#For a discordant cluster supporting a MEI event. How many reads should there at least be in the cluster.
READS_PER_CLUSTER=1
#For discordant clusters coming from both the 5' and 3' side of the same insertion.
#How many bp may the clusters overlap to join the clusters? (to allow for target site duplications)
DISCORDANT_CLUSTERS_MAX_OVERLAP=50
#How many bp may the clusters be away from each other to still join the clusters?
DISCORDANT_CLUSTER_MAX_DISTANCE=600
#Clipped reads coming from one side of the insertion: how much spacing in bp may there be in between the clipped alignment ends
#in order to join these reads into one clipped cluster
MAX_SPACING_OF_CLIPPED_READS=15
#For clusters of clipped reads coming from both 5' and 3' side of insertion:
#How many bp may they overlap or be distanced from each other in order to be joined
MAX_OVERLAP_OF_CLIPPED_CLUSTERS=50
MAX_DISTANCE_OF_CLIPPED_CLUSTERS=20
#How far away may discordant reads be in order to be joined into one cluster?
NEIGHBORHOOD_WINDOW_BP=200
#For a prediction to be kept: how many supporting reads should there be in total for the event
#This is discordant supporting reads + clipped reads
MINIMUM_SUPPORTING_READS=5
#How many reads should there be minimum in a clipped cluster
MINIMUM_INITIAL_SPLIT_CLUSTER_READS=2
#Location of the repeat mask file. To filter prediction which are indiccative of a reference MEI.
REPEATMASK_FILE=input/BGI.scaf.rm.scaffold1.txt
#Temporary directory to write temporary files to. You might need this when the default temporary directory is too small, for now commented out:
#TEMP_DIR=/home/djie/data/public_trio/tmp/" >$Name.mob.properties
## MOBSTER prompt:
java -jar $MOBSTERPATH -properties $Name.mob.properties
## append header to VCF-file
grep "^#" $Name.mobster_predictions.txt >$ID.Mobster.out.vcf

##### Run MELT:
## MELT prompt:
java -jar $MELTPATH Single -bamfile ${Name}.sorted.mkdup.bam -c ${coverage} -d 100000 -e 500 -h $Sequence -n data/melt_library/empty.bed -r 90 -t data/melt_library/SINEtRNA_MELT.zip -w .

#### Process MELT output:
## rename output files:
mv SINEtRNA.bed.list $ID.MELT.bed.list
mv SINEtRNA.final_comp.vcf $ID.MELT.out.vcf
mv SINEtRNA.hum.list $ID.MELT.hum.list
mv SINEtRNA.master.bed $ID.MELT.master.bed
mv SINEtRNA.merged.hum_breaks.sorted.bam $ID.MELT.merged.hum_breaks.sorted.bam
mv SINEtRNA.merged.hum_breaks.sorted.bam.bai $ID.MELT.merged.hum_breaks.sorted.bam.bai
mv SINEtRNA.pre_geno.tsv $ID.MELT.pre_geno.tsv
## create "PASS-only"-file:
grep "^#" $ID.MELT.out.vcf > $ID.MELT.PASS_only.vcf
grep -v "^#" $ID.MELT.out.vcf | awk '$7 == "PASS"' >> $ID.MELT.PASS_only.vcf
## create "PASS-only + homozygote-only"-file:
grep "^#" $ID.MELT.out.vcf > $ID.MELT.homozygote_PASS_only.vcf
grep -v "^#" $ID.MELT.PASS_only.vcf | awk '$10 ~ /1\/1/' >> $ID.MELT.homozygote_PASS_only.vcf
## create "PASS-only + heterozygote-only"-file:
grep "^#" $ID.MELT.out.vcf > $ID.MELT.heterozygote_PASS_only.vcf
grep -v "^#" $ID.MELT.PASS_only.vcf | awk '$10 ~ /1\/0/' >> $ID.MELT.heterozygote_PASS_only.vcf
grep -v "^#" $ID.MELT.PASS_only.vcf | awk '$10 ~ /0\/1/' >> $ID.MELT.heterozygote_PASS_only.vcf

#### Process MOBSTER output:
## convert to VCF:
python3.6 scripts/mobster2vcf_python3.py -i $Name.mobster_predictions.txt -o $ID.Mobster.out.vcf -q 0 -s Read_simulation -r $Sequence

## Compare MOBSTER VCF TEs with original insertions:
python3.6 scripts/compare_insertions.py $ID.insertions.sorted.bed $ID.Mobster.out.vcf $ID #$header
awk '{OFS="\t"; match($11, "MEINFO=[^,]+,"); print $2,$3,$5,"'${coverage}'X",substr($11, RSTART+7, RLENGTH-8),"'${ID}'";}' $ID.comparison.vcf >> $Name.analysis.txt

#### Compile statistics file:
## MOBSTER:
# false positives:
bedtools window -header -v -w 50 -a $ID.Mobster.out.vcf -b $ID.insertions.sorted.bed > $ID.mobster.FP.vcf
# false negatives:
bedtools window -header -v -w 50 -a $ID.insertions.sorted.bed -b $ID.Mobster.out.vcf > $ID.mobster.FN.bed
# true positives:
bedtools window -header -w 50 -b $ID.Mobster.out.vcf -a $ID.insertions.sorted.bed > $ID.mobster.TP.bed
## MELT all detections (raw unfiltered output/ignoring PASS flags):
## false positives:
bedtools window -header -v -w 50 -a $ID.MELT.out.vcf -b $ID.insertions.sorted.bed > $ID.MELT.all.FP.vcf
#false negatives:
bedtools window -header -v -w 50 -a $ID.insertions.sorted.bed -b $ID.MELT.out.vcf > $ID.MELT.all.FN.bed
# true positives:
bedtools window -header -w 50 -a $ID.insertions.sorted.bed -b $ID.MELT.out.vcf > $ID.MELT.all.TP.bed
## MELT "PASS only" (insertions flagged as PASS):
## false positives:
bedtools window -header -v -w 50 -a $ID.MELT.PASS_only.vcf  -b $ID.insertions.sorted.bed > $ID.MELT.PASS_only.FP.vcf
#false negatives:
bedtools window -header -v -w 50 -a $ID.insertions.sorted.bed -b $ID.MELT.PASS_only.vcf > $ID.MELT.PASS_only.FN.bed
# true positives:
bedtools window -header -w 50 -a $ID.insertions.sorted.bed -b $ID.MELT.PASS_only.vcf > $ID.MELT.PASS_only.TP.bed
## MELT "PASS only + homozygote" (insertions flagged as PASS and homozygote):
## false positives:
bedtools window -header -v -w 50 -a $ID.MELT.homozygote_PASS_only.vcf -b $ID.insertions.sorted.bed > $ID.MELT.homozygote_PASS_only.FP.vcf
#false negatives:
bedtools window -header -v -w 50 -a $ID.insertions.sorted.bed -b $ID.MELT.homozygote_PASS_only.vcf > $ID.MELT.homozygote_PASS_only.FN.bed
# true positives:
bedtools window -header -w 50 -a $ID.insertions.sorted.bed -b $ID.MELT.homozygote_PASS_only.vcf > $ID.MELT.homozygote_PASS_only.TP.bed
## MELT "PASS only + heterozygote (insertions flagged as PASS and heterozygote):
## false positives:
bedtools window -header -v -w 50 -a $ID.MELT.heterozygote_PASS_only.vcf -b $ID.insertions.sorted.bed > $ID.MELT.heterozygote_PASS_only.FP.vcf
#false negatives:
bedtools window -header -v -w 50 -a $ID.insertions.sorted.bed -b $ID.MELT.heterozygote_PASS_only.vcf > $ID.MELT.heterozygote_PASS_only.FN.bed
# true positives:
bedtools window -header -w 50 -a $ID.insertions.sorted.bed -b $ID.MELT.heterozygote_PASS_only.vcf > $ID.MELT.heterozygote_PASS_only.TP.bed
## Generate summary TXT files:
echo -en "\n################################################### \n Mobster and MELT statistics for run $ID \n################################################### \n\n\n\n" > stats.$ID.txt
echo "######## Mobster ########" >> stats.$ID.txt
echo -en "\n" >> stats.$ID.txt
echo "# detections:" >> stats.$ID.txt
grep -v "^#" $ID.Mobster.out.vcf | wc -l >> stats.$ID.txt
echo "True Positives:" >> stats.$ID.txt
grep -v "^#" $ID.mobster.TP.bed | wc -l >> stats.$ID.txt
echo "False Negatives:" >> stats.$ID.txt
grep -v "^#" $ID.mobster.FN.bed | wc -l >> stats.$ID.txt
echo "False Positives:" >> stats.$ID.txt
grep -v "^#" $ID.mobster.FP.vcf | wc -l >> stats.$ID.txt
echo -en "\n\n\n" >> stats.$ID.txt
echo "######## MELT ########" >> stats.$ID.txt
echo -en "\n" >> stats.$ID.txt
echo -en "\n" >> stats.$ID.txt
echo "#### (1) All detections  ####" >> stats.$ID.txt
echo -en "\n" >> stats.$ID.txt
echo "# detections:" >> stats.$ID.txt
grep -v "^#" $ID.MELT.out.vcf | wc -l >> stats.$ID.txt
echo "True Positives:" >> stats.$ID.txt
grep -v "^#" $ID.MELT.all.TP.bed | wc -l >> stats.$ID.txt
echo "False Negatives:" >> stats.$ID.txt
grep -v "^#" $ID.MELT.all.FN.bed | wc -l >> stats.$ID.txt
echo "False Positives:" >> stats.$ID.txt
grep -v "^#" $ID.MELT.all.FP.vcf | wc -l >> stats.$ID.txt
echo -en "\n\n" >> stats.$ID.txt
echo "#### (2) PASS only  ####" >> stats.$ID.txt
echo -en "\n" >> stats.$ID.txt
echo "# detections left:" >> stats.$ID.txt
grep -v "^#" $ID.MELT.PASS_only.vcf | wc -l >> stats.$ID.txt
echo "True Positives:" >> stats.$ID.txt
grep -v "^#" $ID.MELT.PASS_only.TP.bed | wc -l >> stats.$ID.txt
echo "False Negatives:" >> stats.$ID.txt
grep -v "^#" $ID.MELT.PASS_only.FN.bed | wc -l >> stats.$ID.txt
echo "False Positives:" >> stats.$ID.txt
grep -v "^#" $ID.MELT.PASS_only.FP.vcf | wc -l >> stats.$ID.txt
echo -en "\n\n" >> stats.$ID.txt
echo "#### (3a) PASS only + homozygote ####" >> stats.$ID.txt
echo -en "\n" >> stats.$ID.txt
echo "# detections left:" >> stats.$ID.txt
grep -v "^#" $ID.MELT.homozygote_PASS_only.vcf | wc -l >> stats.$ID.txt
echo "True Positives:" >> stats.$ID.txt
grep -v "^#" $ID.MELT.homozygote_PASS_only.TP.bed | wc -l >> stats.$ID.txt
echo "False Negatives:" >> stats.$ID.txt
grep -v "^#" $ID.MELT.homozygote_PASS_only.FN.bed | wc -l >> stats.$ID.txt
echo "False Positives:" >> stats.$ID.txt
grep -v "^#" $ID.MELT.homozygote_PASS_only.FP.vcf | wc -l >> stats.$ID.txt
echo -en "\n\n" >> stats.$ID.txt
echo "#### (3b) PASS only + heterozygote ####" >> stats.$ID.txt
echo -en "\n" >> stats.$ID.txt
echo "# detections left:" >> stats.$ID.txt
grep -v "^#" $ID.MELT.heterozygote_PASS_only.vcf | wc -l >> stats.$ID.txt
echo "True Positives:" >> stats.$ID.txt
grep -v "^#" $ID.MELT.heterozygote_PASS_only.TP.bed | wc -l >> stats.$ID.txt
echo "False Negatives:" >> stats.$ID.txt
grep -v "^#" $ID.MELT.heterozygote_PASS_only.FN.bed | wc -l >> stats.$ID.txt
echo "False Positives:" >> stats.$ID.txt
grep -v "^#" $ID.MELT.heterozygote_PASS_only.FP.vcf | wc -l >> stats.$ID.txt
echo -en "\n\n\n" >> stats.$ID.txt
echo "#### ESAT insertions ####" >> stats.$ID.txt
grep -v "^#" $ID.insertions.sorted.bed | wc -l >> stats.$ID.txt
echo -en "\n\n" >> stats.$ID.txt

#### Clean up:
## remove all diles except nohup.ID.out + stats.ID.txt:
rm -rf  library*
rm $ID.comparison.vcf
rm $ID.insertions.sorted.bed
rm $ID.MELT.bed.list
rm $ID.MELT.hum.list
rm $ID.MELT.master.bed
rm $ID.MELT.merged.hum_breaks.sorted.bam
rm $ID.MELT.merged.hum_breaks.sorted.bam.bai
rm $ID.MELT.pre_geno.tsv
rm $Name.analysis.txt
rm $Name.bam
rm $Name.fa
rm $Name.mob.properties
rm $Name.mobster_anchors.bai
rm $Name.mobster*
rm $Name.sam
rm -rf $Name.sorted*
mkdir $ID.results
mv $ID.* nohup.$ID.* $ID.results
echo -en "\n\n\n### fertig XY###" >> nohup.out
mv nohup.out nohup.$ID.out

#### Old analysis part, wasn't used anymore for in this ESAT version:
#Average deviation 
##grep -v '^#' analysis.txt | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {(SUM += abs($3-$1))} END{print "\n Average deviation from original insertion point:" SUM/NR "bp"}' >>analysis.txt
#: <<'END'
####uncomment to count all discovered intersections
#HEADER grep "^#" simulated_sequence.mobster_predictions.txt >analysis.txt 
#echo "#Number of originally inserted TEs: $Number_of_insertions" >>analysis.txt
#discovered_insertions=$(grep -c 'scaffold1' Intersections.vcf)
#	echo "#Discovered Insertions:$discovered_insertions" >>analysis.txt
#discovered_SINEs=$(grep -c 'SINE' Intersections.vcf)
#	echo "#Discovered SINEs:$discovered_SINEs" >>analysis.txt
#discovered_LINEs=$(grep -c 'L1' Intersections.vcf)
#	echo "#Discovered LINEs:$discovered_LINEs" >>analysis.txt
#END
