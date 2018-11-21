# ESAT – Element Simulating Analysis Tool

Tool that simulates the insertion of elements in a genomic sequence and evaluates the detection rate with Mobster.


To run phython scripts individually use python3.4

Additional Programs that are required to run the pipeline:

Bedtools
Bwa
Mobster
Mosaik
Picard
Simseq
Samtools


Additional Files

Error Profile for SimSeq
Repetitive elements for Mobster


Start the pipeline with

bash ESAT.sh <genomic fragment> <TE sequence><Name for the output files><Coverage for the read simulation><Number of insertions><ID>

or

Start loop with

bash loop.sh <genomic fragment> <TE sequence> <Output name>

- loop generates new ID for every individual run but sam, bam, fastq (etc) files will be overwritten in each run!



The pipeline consists of:

- Sequence Simulation sim_seq.py
	output: simulated sequence (fasta)
	      	 insert breakpoints (.bed)



- bwa indexing of the new sequence
- sort insert brekpoints file
- calculation of read number  read_number.py

- read Simulation
	-change arguments -1 and -2 for different read length (set to 90), use error profile here, 	  	change insert size (set to 500)

- mapping with bwa
- generate mobster properties file
- generate vcf file  mobster2vcf_python3.py
- compare real breakpoints with simulated breakpoints compare_insertions.py
- generate output summary file


Output

- Use Name.sorted.mkdup.bam to view mapped reads
- Use ID.insertions.sorted.bed for original breakpoints (insert locations)
- Use Name.analysis.txt output with original and predicted preakpoints (original beakpoint, breakpoint end, Mobster prediction, Coverage, Type, ID)



















