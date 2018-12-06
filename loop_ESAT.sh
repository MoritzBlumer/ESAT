#!/bin/bash

Sequence=$1  #original sequence file
Insertion=$2 #insertion file
Name=$3	     #output name
#coverage=$4
#Number_of_insertions=$5  #optional argument --n number of insertions, default 10
#runs=$6
#create the end-file
#COUNTER=$((COUNTER+1))


#Name=$Name.$COUNTER


#rm $Name.analysis.txt


#echo -e "Ins""\t""End""\t""Mob""\t""Cov""\t""Type""\t""ID" >$Name.analysis.txt



for run in {1..10}
	#for coverage in {}
		do
			COUNTER=$((COUNTER+1))
			ID=$Name.$COUNTER
			#echo $ID
			bash ESAT.sh $Sequence $Insertion $Name 2 100 $ID

		done


