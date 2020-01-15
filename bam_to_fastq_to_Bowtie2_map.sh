#!/bin/bash

# usage: bam_to_fastq_to_Bowtie2_map.sh <threads> <reference.fa> <list_of_bam_files.txt>

# ensure that the reference file and all bam files contain only one '.' and that this dot preceeds the suffix (e.g. '.fa' or '.bam')

# for example:
# sample_1.bam
# sample_2.bam

if [[ -n $1 ]] && [[ -n $2 ]] && [[ -n $3 ]]
then
    echo ""; echo ""
    echo "number of processing cores requested: "$1
    echo "reference filename: "$2
    echo "list of input bam filenames: "$3
    echo ""; echo ""

    #load env variables, etc
    module load samtools
    module load bedtools
    module load bowtie2

    #take prefixes from reference filnames
    refprefix=$(echo "$2" | cut -f 1 -d '.')
    counter=0
    failed=0

    #for each line of text file, split into 2 strings with whitespace delimiter, and assign each to read1 and read2.
    #then map readfiles relating to these filenames
    input=$3
    while IFS= read -r line
    do
        sample=$(echo $line | cut -f 1 -d '.')
        counter=$(expr $counter + 1)
        echo ""; echo ""
        echo "Processing sample No."$counter", "$sample
        echo ""; echo ""

        #sort bam by query name
        echo "sorting bam file"
        samtools sort -n -@ $1 $line > $sample"_sorted.bam"; stx=$?

        #convert bam to fastq
        if [[ $stx -eq 0 ]]
        then
        echo "converting bam to fastq and fasta"
            bedtools bamtofastq -i $sample"_sorted.bam" -fq $sample"_first.fq" -fq2 $sample"_second.fq"
            bedx=$?
            echo ""; echo ""
        else
            echo "samtools failed"; echo ""; echo ""
            failed=1
            break 2
        fi

        #map reads to reference
        if [[ $bedx -eq 0 ]]
        then
            echo "generating bam file..."
            bt2thread=$(expr $1 - 1)
            bowtie2 -p $bt2thread --fast -x $refprefix -1 $sample"_first.fq" -2 $sample"_second.fq" | \
            samtools view -@ 1 -b -h -o $sample"_to_"$refprefix".bam"
            btstx=$?
            echo ""; echo ""
        else
            echo "bedtools failed"; echo ""; echo ""
            failed=1
            break 2
        fi


        if [[ $btstx -eq 0 ]]
        then
            echo "Finished processing sample No."$counter", "$sample
        else
            echo "bowtie2/samtools pipe failed"; echo ""; echo ""
            failed=1
            break 2
        fi


    done < "$input"

    if [[ $failed -eq 0 ]]
    then
        rm *sorted* *first* *second*
    fi

else
  echo "An argument is missing"
fi
