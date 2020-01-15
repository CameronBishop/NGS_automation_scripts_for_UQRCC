#!/bin/bash

# usage: prok_PE_mapper_v2.sh <threads> <genome.fa> <readfile_list.txt>

# a splice-unaware RNAseq mapper for trimming and mapping paired end reads to prokaryotic reference.
# requires a text file containing one line per sample. Each line must have contain a sample identifier,
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
    module load trimmomatic
    module load fastqc
    module load Java
    module load bowtie2
    module load samtools

    mkdir fastQC_output

    #take prefixes from reference filnames
    refprefix=$(echo "$2" | cut -f 1 -d '.')

    #assign call and options to variable for use in trimmomatic
    call_trim="java -jar /opt/biotools/trimmomatic/trimmomatic-0.35.jar PE -threads $1 -phred33"
    trimoptions="ILLUMINACLIP:/opt/biotools/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:20" \

    #for each line of text file, split into 2 strings with whitespace delimiter, and assign each to read1 and read2.
    #then map readfiles relating to these filenames
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

        #Trim reads for adapter sequences
        $call_trim \
        $read1 $read2 \
        $read1".paired" $read1".unpaired" \
        $read2".paired" $read2".unpaired" \
        $trimoptions; trimx=$?
        echo ""; echo ""

        #generate QC report for paired reads
        if [ $trimx -eq 0 ]
        then
            echo "running fastqc on "$read1".paired"; fastqc -t $1 -o ./fastQC_output -q $read1".paired"
            echo "running fastqc on "$read2".paired"; fastqc -t $1 -o ./fastQC_output -q $read2".paired"
            echo "running fastqc on "$read1".unpaired"; fastqc -t $1 -o ./fastQC_output -q $read1".unpaired"
            echo "running fastqc on "$read2".unpaired"; fastqc -t $1 -o ./fastQC_output -q $read2".unpaired"
            fqcx=$?; echo ""; echo ""
        else
            echo "trimmomatic failed"; echo ""; echo ""
            break 2
        fi

        #map reads to reference, and sort bam file by coordinate
        if [[ $fqcx -eq 0 ]]
        then
            echo "generating bam file..."
            btcores=$(expr $1 - 4)
            bowtie2 -p $btcores --fast --no-mixed --no-unal -x $refprefix -1 $read1".paired" -2 $read2".paired" -U $read1".unpaired",$read2".unpaired" | \
            samtools sort -@ 4 > $sample"_to_"$refprefix".bam"
            btstx=$?
            echo ""; echo ""
        else
            echo "fastQC failed"; echo ""; echo ""
            break 2
        fi


        if [[ $btstx -eq 0 ]]
        then
            echo $sample" mapping is complete"; echo ""; echo""
        else
            echo $sample" bowtie2/samtools pipe failed"; echo ""; echo""
        fi

        echo "Finished processing read-pair No. "$counter" from "$3

    done < "$input"

else
  echo "An argument is missing"
fi
