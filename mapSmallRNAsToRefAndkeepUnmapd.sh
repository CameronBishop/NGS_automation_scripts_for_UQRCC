#!/bin/bash

# usage: mapSmallRNAsToRefAndkeepUnmapd.sh <threads(int)> <reference.fa> <reads_suffix>
#
# The script takes raw smallRNAseq files and a reference genome.
# The script trims adapter sequences, generates fastQC reports, maps reads to reference,
# then collects unmapped reads and outputs them as a fasta and a fastq file.
# For use on the UQ RCC HPC (awoonga, tinaroo, etc).

# BEFORE RUNNING:
#
# All fastq filenames must include exactly 1 'dot' and this dot should deliniate the prefix from the suffix.
# (e.g. 'reads.fastq' or 'reads.fq' or 'some_reads.fq')
# The suffix (with no dot) will need to be included as the third argument so that the script knows what files
# to use as input. As long as the read files are all in the current working directory the script will process them.
#
# Place all read files and reference files together in one directory, and call this script from that same directory.
# The script can be in whatever directory you want.
#
# The adapter sequence used here is the TruSeq3-SE.fa 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA' from the NEBNext multiplex smRNA kit.
#
# Create a bowtie2 index for the reference file (note that bowtie2 requires a '.fa' suffix for the reference, otherwise it throws an error).
#
# This script is designed to work on the UQ RCC cluster (awoonga, tinaroo, etc). The 'module load' commands will probably only work
# on this system.
#
# issue the command 'chmod 700 ./smallRNAseqcountV4.sh' from the dir in which the script is located. (chmod will grant executive privileges).


if [[ -n $1 ]] && [[ -n $2 ]] && [[ -n $3 ]]
then
  echo "number of processing cores requested: "$1
  echo "reference genome file: "$2
  echo "using input read files with "$3" suffix"
  echo ""; echo ""

  #load env variables, etc
  module load Java
  module load trimmomatic
  module load fastqc
  module load bowtie2
  module load samtools
  module load bedtools
  module load bbmap


  #make a bowtie2 prefix
  idxprefix=$(echo "$2" | cut -f 1 -d '.')

  #create an array of fastq filenames
  shopt -s nullglob
  inreads=(*.$3)

  #print number of read files detected, and list them.
  echo "found ${#inreads[@]} read files in working directory"
  for i in "${inreads[@]}"; do echo "$i"; done

  #iterate over inreads array
  for i in "${inreads[@]}"
  do
    #make prefix for readfile
    readspref=$(echo $i | cut -f 1 -d ".")

    #make filenames for trimmed fastqs
    trimmed=$readspref".trimmed.fq"

    #Trim reads for adapter sequences
    java -jar /opt/biotools/trimmomatic/trimmomatic-0.35.jar SE -threads $1 -quiet\
    $i \
    $trimmed\
    ILLUMINACLIP:/opt/biotools/trimmomatic/adapters/TruSeq3-SE.fa:2:30:9
    trimx=$?
    echo ""; echo ""

    #filter reads based on size
    if [[ $trimx -eq 0 ]]
    then
        echo "filtering reads with bbmap..."
        reformat.sh in=$trimmed out=$trimmed".filt.fastq" minlength=18 maxlength=24
        filx=$?
        mv $trimmed".filt.fastq" $trimmed
        echo ""; echo ""
    else
    	echo "trimmomatic failed"
        break 2
    fi

    #generate QC report
    if [[ $filx -eq 0 ]]
    then
        echo "running fastqc..."
        fastqc -q $trimmed
        fqcx=$?
        echo ""; echo ""
    else
        echo "BBmap failed"
        break 2
    fi

    #map reads to reference, convert to bam, keep only unmapped reads
    if [[ $fqcx -eq 0 ]]
    then
        echo "generating bam file..."
        btcores=$(expr $1 - 1)
        bowtie2 -p $btcores -D 10 -R 2 -N 0 -L 18 -i S,1,1.75 -x $idxprefix -U $trimmed | samtools view -h -b -@ 1 -f 4 > $readspref"_unm.bam"
        btstx=$?
        echo ""; echo ""
    else
        echo "fastqc failed"
        break 2
    fi

    #convert bam to fasta
    if [[ $btstx -eq 0 ]]
    then
        echo "converting bam to fastq and fasta"
    	bedtools bamtofastq -i $readspref"_unm.bam" -fq $readspref"_unm.fq"; bedx=$?
        sed -n '1~4s/^@/>/p;2~4p' $readspref"_unm.fq" > $readspref"_unm.fa"; sedx=$?
        echo ""; echo ""
    else
    	echo "bowtie2 or samtools failed"
        break 2
    fi

    if [[ $bedx -ne 0 ]] && [[ $sedx -ne 0 ]]
    then
        echo "bedtools or sed failed"
    fi

  done

    #move temp files
	mkdir outputfiles
    mkdir tempfiles
    mv *unm.fq outputfiles
    mv *unm.fa outputfiles
    rm *unm.bam
    mv *trimmed* tempfiles
	echo "job completed"

else
  echo "An argument is missing"
fi
