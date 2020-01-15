#!/bin/bash

# usage: prok_SE_mapper_no_trim_no_QC.sh <threads> <genome.fa> <fastq>

# a splice-unaware RNAseq mapper for mapping single end reads to prokaryotic reference. Outputs bam, fastq and fasta file for each alignment
#
# All input fastq filenames must include exactly 1 'dot' and this dot should delineate the prefix from the suffix.
# (e.g. 'reads.fastq' or 'reads.fq' or 'some_reads.fq')
# The suffix (with no dot) will need to be included as the third argument so that the script knows what files
# to use as input. As long as the read files are all in the current working directory the script will process them.
#
# This script is designed to work on the UQ RCC cluster (awoonga, tinaroo, etc). The 'module load' commands will probably only work
# on this system.

if [[ -n $1 ]] && [[ -n $2 ]] && [[ -n $3 ]]
then
    echo ""; echo ""
    echo "number of processing cores requested: "$1
    echo "reference fasta filename: "$2
    echo "reading files with suffix: "$3
    echo ""; echo ""

    #load env variables, etc
    module load bowtie2
    module load samtools
    module load bedtools

    #make a bowtie2 prefix
    idxprefix=$(echo "$2" | cut -f 1 -d '.')

    #create an array of fastq filenames
    shopt -s nullglob
    inreads=(*.$3)

    #print number of read files detected, and list them.
    echo "found ${#inreads[@]} read files in working directory"
    echo ""
    for i in "${inreads[@]}"; do echo "$i"; done
    echo ""; echo ""

    #iterate over inreads array
    for i in "${inreads[@]}"
    do
        #make prefix for readfile
        readspref=$(echo $i | cut -f 1 -d ".")

        #make filenames for trimmed fastqs
        echo "Processing readfile: "$i
        echo ""; echo ""
        btcores=$(expr $1 - 1)
        bowtie2 -p $btcores --no-unal -D 10 -R 2 -N 0 -L 18 -i S,1,1.75 -x $idxprefix -U $i | samtools view -h -b -@ 1 > $readspref"_to_"$idxprefix".bam"
        btx=$?
        echo $readspref" mapping is complete"; echo ""; echo""

        #convert bam to fasta
        if [[ $btx -eq 0 ]]
        then
            echo "converting bam to fastq and fasta"
            bedtools bamtofastq -i $readspref"_to_"$idxprefix".bam" -fq $readspref"_to_"$idxprefix".fq"
            sed -n '1~4s/^@/>/p;2~4p' $readspref"_to_"$idxprefix".fq" > $readspref"_to_"$idxprefix".fa"
            bedx=$?
        else
            echo "bowtie2 failed"
            break 3
        fi


        if [[ $bedx -eq 0 ]]
        then
            echo "Finished processing "$readspref
            echo ""; echo ""
        else
            echo "bedtools failed"; echo ""; echo""
        fi

    done < "$inreads"

else
  echo "An argument is missing"
fi
