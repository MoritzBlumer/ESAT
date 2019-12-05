# ESAT – Element Simulating Analysis Tool

ESAT is a tool that simulates the insertion of elements in a genomic sequence and evaluates the detection rate with Mobster.


## Requirements

To run phython scripts individually use python3.4

### Programs 
Additional programs that are required to run the pipeline:

* Bedtools (version 2.17, http://bedtools.readthedocs.org/en/latest/)
* Bwa (version 0.7.5a-r405, https://sourceforge.net/projects/bio-bwa/files/)
* Mobster (http://sourceforge.net/projects/mobster)
* MELT (http://melt.igs.umaryland.edu/index.php)
* Mosaik (version 2.2, https://github.com/wanpinglee/MOSAIK.git/)
* Picard (https://broadinstitute.github.io/picard/)
* Simseq (https://github.com/jstjohn/SimSeq)
* Samtools (version 1.3, https://github.com/samtools/samtools/)

### Additional Files

* Error Profile for SimSeq (https://github.com/jstjohn/SimSeq)
* Repetitive elements for Mobster/MELT

### Before first use:
* set paths to java .jar for Mobster, MELT, SimSeq, SamToFastq, MarkDuplicates and CollectInsertSizeMetrics in ESAT.sh


## Run the Pipeline

```bash
bash ESAT.sh <genomic fragment> <TE sequence> <Name of the output files> <Coverage of the read simulation> <Number of insertions> <ID>
```


## The pipeline

### Steps in the pipeline:

* Sequence Simulation sim_seq.py; Output: simulated sequence (fasta), insert breakpoints (.bed)



* bwa indexing of the new sequence
* Sort insert brekpoints file
* Calculation of read number read_number.py

* read Simulation
	 (change arguments -1 and -2 for different read length (set to 90), use error profile here, 	  	change insert size (set to 500))

* Mapping with bwa
* Generate mobster properties file
* Run Mobster and MELT on the simulated sequence
* Generate vcf file  mobster2vcf_python3.py
* Compare real breakpoints with simulated breakpoints compare_insertions.py
* Generate output summary file

### Output:

* for each ESAT run, a file with statistics of the run will be generated in "./statistics_files"














