#!/bin/bash

# usage: fasta_to_fastq_remove_wrap.sh <list_of_fasta_files.txt> <header_prefix>

# note: In order for python to know which lines are headers, the user must specify a 'header_iq'
# this should be a non-unique string that is common to every header, and be at the beginning of the line
# Simply specifying '@ is not satisfactory, as '@' also occurs in phred scores, so will fuck-up the output.

# fastq header example:
# @locus_tag_*******
# in this example the header prefix would be '@locus_tag_'

# list.txt example:
# sample_1.fa
# sample_2.fa

if [[ -n $1 ]] && [[ -n $2 ]]
then
    echo ""; echo ""
    echo "list of input fasta filenames: "$1
    echo "header identifier: " $2
    echo ""; echo ""

    counter=0

    #load env variables
    module load python

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
    echo "        if line.startswith(sys.argv[3]):" >> fastq2fasta.py
    echo "            getlines = True" >> fastq2fasta.py
    echo "        if re.match(plus, line):" >> fastq2fasta.py
    echo "            getlines = False" >> fastq2fasta.py
    echo "        if getlines:" >> fastq2fasta.py
    echo "            array.append(line.strip('\n'))" >> fastq2fasta.py
    echo "temp = []" >> fastq2fasta.py
    echo "for i in array:" >> fastq2fasta.py
    echo "    if i.startswith('@'):" >> fastq2fasta.py
    echo "        temp = []" >> fastq2fasta.py
    echo "        fasta.write(i.replace('@', '>') + '\n')" >> fastq2fasta.py
    echo "    if not i.startswith('@'):" >> fastq2fasta.py
    echo "        temp.append(i.strip())" >> fastq2fasta.py
    echo "        if array.index(i) == (len(array) - 1):" >> fastq2fasta.py
    echo "            fasta.write(''.join(temp) + '\n')" >> fastq2fasta.py
    echo "        elif array[array.index(i) + 1].startswith('@'):" >> fastq2fasta.py
    echo "            fasta.write(''.join(temp) + '\n')" >> fastq2fasta.py


    #for each line of text file, remove the suffix make a scalar of it for later use
    input=$1
    while IFS= read -r line
    do
        sample=$(echo $line | cut -f 1 -d '.')
        counter=$(expr $counter + 1)
        echo ""; echo ""
        echo "Processing sample No."$counter
        echo $sample
        echo ""; echo ""

        #convert fastq to fasta and remove textwrapping
        echo "generating fasta file"
        python fastq2fasta.py $line $sample".fa" $2
        pyx=$?
        echo ""; echo ""

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
