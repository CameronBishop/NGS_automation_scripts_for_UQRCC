#!/bin/bash

#usage: fastq2counts.sh <threads> <readfile_list.txt>

if [[ -n $1 ]] && [[ -n $2 ]]
then
    echo ""; echo ""
    echo "number of processing cores requested: "$1
    echo "txt file containing filenames of paired input readfiles: "$2
    echo ""; echo ""


    #load env variables, etc
    module load trimmomatic
    module load fastqc
    module load Java
    mkdir fastqc_results

    #assign call and options to variable for use in trimmomatic
    call_trim="java -jar /opt/biotools/trimmomatic/trimmomatic-0.35.jar PE -threads $1 -phred33"
    trimoptions="ILLUMINACLIP:/opt/biotools/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:20" \

    #for each line of text file, split into 2 strings with whitespace delimiter, and assign each to read1 and read2
    input=$2
    while IFS= read -r line; do
        read1=$(echo $line | cut -f 1 -d ' ')
        read2=$(echo $line | cut -f 2 -d ' ')
        readspref=$(echo $read1 | cut -f 1 -d '.')

        #Trim reads for adapter sequences
        echo "trimming....."; echo ""; echo ""
        $call_trim \
        $read1 $read2 \
        $read1".trimmed.paired" $read1".trimmed.unpaired" \
        $read2".trimmed.paired" $read2".trimmed.unpaired" \
        $trimoptions; trimx=$?
        echo ""; echo ""

        #generate QC report for paired reads
        if [ $trimx -eq 0 ]
        then
            echo "running fastqc on "$read1".trimmed.paired"; fastqc -t $1 -o ./fastqc_results -q $read1".trimmed.paired"
            echo "running fastqc on "$read2".trimmed.paired"; fastqc -t $1 -o ./fastqc_results -q $read2".trimmed.paired"
            echo "running fastqc on "$read1".trimmed.unpaired"; fastqc -t $1 -o ./fastqc_results -q $read1".trimmed.unpaired"
            echo "running fastqc on "$read2".trimmed.unpaired"; fastqc -t $1 -o ./fastqc_results -q $read2".trimmed.unpaired"
            fqcx=$?; echo ""; echo ""
        else
            echo "trimmomatic failed"; echo ""; echo ""
        fi

    done < "$input"

else
  echo "An argument is missing"
fi
