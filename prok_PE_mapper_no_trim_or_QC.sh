#!/bin/bash

# usage: prok_PE_mapper_no_trim_or_QC.sh <threads> <genome.fa> <readfile_list.txt>

# a splice-unaware RNAseq mapper for mapping paired end reads to prokaryotic reference.
# requires a text file containing one line per sample. Each line must contain a sample identifier,
# forward read filename, and reverse read filename.

# for example:
# sample_1 s1_forward.fq s1_reverse.fq
# sample_2 s2_forward.fq s2_reverse.fq

if [[ -n $1 ]] && [[ -n $2 ]] && [[ -n $3 ]]
then
    echo ""; echo ""
    echo "number of processing cores requested: "$1
    echo "reference fasta filename: "$2
    echo "list of input fastq filenames: "$3
    echo ""; echo ""

    #load env variables, etc
    module load bowtie2
    module load samtools

    #take prefixes from reference filnames
    refprefix=$(echo "$2" | cut -f 1 -d '.')
    input=$3
    while IFS= read -r line
    do
        sample=$(echo $line | cut -f 1 -d ' ')
        read1=$(echo $line | cut -f 2 -d ' ')
        read2=$(echo $line | cut -f 3 -d ' ')
        counter=$(expr $counter + 1)
        echo ""; echo ""
        echo "Processing read-pair No."$counter
        echo $read1; echo $read2
        echo ""; echo ""

        bowtie2 -p $(expr $1 - 1) --fast --no-mixed --no-unal -x $refprefix -1 $read1 -2 $read2 | \
        samtools sort -@ 1 > $sample"_to_"$refprefix".bam"
        btstx=$?
        echo ""; echo ""

        if [[ $btstx -eq 0 ]]
        then
            echo $sample" mapping is complete"; echo ""; echo""
        else
            echo $sample" bowtie2/samtools pipe failed"; echo ""; echo""
            break 2
        fi

        echo "Finished processing read-pair No. "$counter" from "$3

    done < "$input"

else
  echo "An argument is missing"
fi
