#!/bin/bash

# usage: bam2consensusV2.0.sh <threads> <reference.fa> <list_of_bam_files.txt> <header_identifier>

# script takes a list of bam files. for each bam, it calls the consensus sequence, and outputs a fasta and fastq file

# note: In order for python to know which lines are headers, the user must specify a 'header_identifier'
# this should be a non-unique string that is common to every header, and be at the beginning of the line
# Simply specifying '@ is not satisfactory, as '@' also occurs in phred scores.

# fastq header example:
# @locus_tag_*******
# in this example the header prefix would be '@locus_tag_'

# list.txt example:
# sample_1.bam
# sample_2.bam

if [[ -n $1 ]] && [[ -n $2 ]] && [[ -n $3 ]] && [[ -n $4 ]]
then
    echo ""; echo ""
    echo "number of processing cores requested: "$1
    echo "reference filename: "$2
    echo "list of input bam filenames: "$3
    echo "header identifier: " $4
    echo ""; echo ""

    #load env variables, etc
    module load samtools
    module load bcftools
    module load python

    #take prefixes from reference filnames
    refprefix=$(echo "$2" | cut -f 1 -d '.')
    counter=0

    #create python script to convert fastq to fasta
    touch fastq2fasta.py
    echo "import sys" >> fastq2fasta.py
    echo "import re" >> fastq2fasta.py
    echo "fasta = open(sys.argv[2], 'w+')" >> fastq2fasta.py
    echo "getlines = False" >> fastq2fasta.py
    echo "array = []" >> fastq2fasta.py
    echo "plus = ('\+$')" >> fastq2fasta.py
    echo "with open(sys.argv[1]) as fastqin:" >> fastq2fasta.py
    echo "     for line in fastqin:" >> fastq2fasta.py
    echo "        if line.startswith('@line'):" >> fastq2fasta.py
    echo "            getlines = True" >> fastq2fasta.py
    echo "        if re.match(plus, line):" >> fastq2fasta.py
    echo "            getlines = False" >> fastq2fasta.py
    echo "        if getlines:" >> fastq2fasta.py
    echo "            array.append(line.strip('\n'))" >> fastq2fasta.py
    echo "temp = []" >> fastq2fasta.py
    echo "for i in array:" >> fastq2fasta.py
    echo "    if i.startswith('@'):" >> fastq2fasta.py
    echo "        temp = []" >> fastq2fasta.py
    echo "        fasta.write(i)" >> fastq2fasta.py
    echo "    if not i.startswith('@'):" >> fastq2fasta.py
    echo "        temp.append(i.strip())" >> fastq2fasta.py
    echo "        if array.index(i) == (len(array) - 1):" >> fastq2fasta.py
    echo "            fasta.write(''.join(temp))" >> fastq2fasta.py
    echo "        elif array[array.index(i) + 1].startswith('@'):" >> fastq2fasta.py
    echo "            fasta.write(''.join(temp))" >> fastq2fasta.py


    #for each line of text file, remove the suffix make a scalar of it for later use
    input=$3
    while IFS= read -r line
    do
        sample=$(echo $line | cut -f 1 -d '.')
        counter=$(expr $counter + 1)
        echo ""; echo ""
        echo "Processing sample No."$counter
        echo $sample
        echo ""; echo ""

        #sort and index
        echo "sorting and indexing bam file"
        samtools sort -@ $1 $line > $sample"_sorted.bam"; samtools index $sample"_sorted.bam"; stx=$?; echo ""; echo ""

        if [[ $stx -eq 0 ]]
        then
            #pipe the vcf to call a consensus and generate a fastq file
            echo "calling consensus sequence"
            samtools mpileup -d 8000 -Auf $2 $sample"_sorted.bam" | bcftools call --threads $(expr $1 - 2) -c | \
            vcfutils.pl vcf2fq > $sample".fq"; bcfx=$?; echo ""; echo ""
        else
            echo "sort | index pipe failed"; echo ""; echo ""
            break 2
        fi


        if [[ $bcfx -eq 0 ]]
        then
            #convert fastq to fasta and remove textwrapping
            echo "generating fasta file"
            python fastq2fasta.py $line $sample".fa" $4
            pyx=$?
            echo ""; echo ""
        else
            echo "mpileup | call pipe failed"; echo ""; echo ""
            break 2
        fi


        if [[ $pyx -eq 0 ]]
        then
            echo "finished generating fasta file"; echo ""; echo ""
        else
            echo "python failed"; echo ""; echo ""
            break 2
        fi

        echo "Finished processing sample No."$counter

    done < "$input"

    rm fastq2fasta.py

else
  echo "An argument is missing"
fi
