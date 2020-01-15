#!/bin/bash

# usage: prokEndosymRNAseqCount_v3.sh <threads> <host_geneome.fa> <symbiont_genome.fa> <symbiont_features.gff> <featuretype> <readfile_list.txt>

# The script takes as input: untrimmed RNAseq files, a host reference genome, an endosymbiont/pathogen reference genome, and an endosymbiont/pathogen annotation file.
# The script trims the reads on quality and length (minimum 20 bases), generates fastQC reports, maps reads to a eukaryotic host reference,
# collects unmapped reads, maps these to prokaryotic endosymbiont/pathogen reference, counts reads mapping to gene features of feature-type specified as an argument,
# and generates a count matrix for use in DESeq2 or similar.
#
# The adapter sequences searched for by trimmomatic are located in the TruSeq3-PE.fa file in /opt/biotools/trimmomatic/adapters/.
#
# This script is designed to work on the UQ RCC cluster (awoonga, tinaroo, etc). The 'module load' commands change envrironmental
# variables according to the requirements of each software package. If using on another server, this must be done manually.


# BEFORE RUNNING:
#
# For this script to work you must first install htseq-count to your home directory "~/.local/bin" because it does not
# exist on the HPC - at least it did not when I wrote this.
# for instructions on downloading and installing htseq (which is a python package), and dependencies see:
# http://www2.rcc.uq.edu.au/hpc/guides/index.html?secure/Installation_of_Packages_for_Scripting_Languages.html#Python-distutils
# and
# https://htseq.readthedocs.io/en/release_0.11.1/install.html
#
# Create a bowtie2 index for both reference files.
#
# issue the command 'chmod 700 ./smallRNAseqcountV4.sh' from the dir in which the script is located. (chmod will grant executive privileges).

if [[ -n $1 ]] && [[ -n $2 ]] && [[ -n $3 ]] && [[ -n $4 ]] && [[ -n $5 ]] && [[ -n $6 ]]
then
    echo ""; echo ""
    echo "number of processing cores requested: "$1
    echo "host reference genome file: "$2
    echo "symbiont reference genome file: "$3
    echo "gene feature file: "$4
    echo "feature type: "$5
    echo "list of input fastq files: "$6
    echo ""; echo ""

    #add htseq-count to $PATH
    PATH=$PATH":~/.local/bin"

    #load env variables, etc
    module load trimmomatic
    module load fastqc
    module load Java
    module load bowtie2
    module load tophat
    module load samtools
    module load python
    module load bedtools

    mkdir temp_files
    mkdir tophat_temp

    #take prefixes from reference filnames
    hostprefix=$(echo "$2" | cut -f 1 -d '.')
    symprefix=$(echo "$3" | cut -f 1 -d '.')

    #assign call and options to variable for use in trimmomatic
    call_trim="java -jar /opt/biotools/trimmomatic/trimmomatic-0.35.jar PE -threads $1 -phred33"
    trimoptions="ILLUMINACLIP:/opt/biotools/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:20" \

    #for each line of text file, split into 2 strings with whitespace delimiter, and assign each to read1 and read2
    counter=0 #for echoing progress through reads.txt
    failed=0 #if this stops being 0, no count matrix will be compiled.
    input=$6
    while IFS= read -r line
    do
        read1=$(echo $line | cut -f 1 -d ' ')
        read2=$(echo $line | cut -f 2 -d ' ')
        readspref=$(echo $read1 | cut -f 1 -d '.')
        counter=$(expr $counter + 1)
        echo ""; echo ""
        echo "Processing read-pair No."$counter
        echo $read1; echo $read2
        echo ""; echo ""

        #Trim reads for adapter sequences
        $call_trim \
        $read1 $read2 \
        $read1".trimmed.paired" $read1".trimmed.unpaired" \
        $read2".trimmed.paired" $read2".trimmed.unpaired" \
        $trimoptions; trimx=$?
        echo ""; echo ""

        #generate QC report for paired reads
        if [ $trimx -eq 0 ]
        then
            echo "running fastqc on "$read1".trimmed.paired"; fastqc -t $1 -o ./temp_files -q $read1".trimmed.paired"
            echo "running fastqc on "$read2".trimmed.paired"; fastqc -t $1 -o ./temp_files -q $read2".trimmed.paired"
            echo "running fastqc on "$read1".trimmed.unpaired"; fastqc -t $1 -o ./temp_files -q $read1".trimmed.unpaired"
            echo "running fastqc on "$read2".trimmed.unpaired"; fastqc -t $1 -o ./temp_files -q $read2".trimmed.unpaired"
            fqcx=$?; echo ""; echo ""
        else
            echo "trimmomatic failed"; echo ""; echo ""
            failed=1
            break 2
        fi

        #map reads to host, collect unmapped reads, and sort by name for fastq conversion
        if [[ $fqcx -eq 0 ]]
        then
            tophat2 -o tophat_temp -p $1 --b2-very-fast --max-intron-length 5000 --max-segment-intron 5000 $hostprefix \
            $read1".trimmed.paired",$read1".trimmed.unpaired" $read2".trimmed.paired",$read2".trimmed.unpaired"; thx=$?
    		samtools sort -n -@ $1 ./tophat_temp/unmapped.bam > $readspref"_unm_sorted.bam"; stx=$?
            mv ./tophat_temp/logs ./tophat_temp/$readspref"_tophat_logs"
            mv ./tophat_temp/$readspref"_tophat_logs" ./temp_files
            echo ""; echo ""
        else
            echo "fastqc failed"; echo ""; echo ""
            failed=1
            break 2
        fi

        #convert bam to fastq
        if [[ $thx -eq 0 ]] && [[ $stx -eq 0 ]]
        then
        echo "converting bam to fastq and fasta"
        	bedtools bamtofastq -i $readspref"_unm_sorted.bam" -fq $readspref"_first_unm.fq" -fq2 $readspref"_second_unm.fq"
        	bedx=$?
            rm $readspref"_unm_sorted.bam"
        	echo ""; echo ""
        else
        	echo "tophat2 or samtools failed"; echo ""; echo ""
            failed=1
            break 2
        fi

        #map reads to symbiont
        if [[ $bedx -eq 0 ]]
        then
            echo "generating bam file..."
            btcores=$(expr $1 - 4)
            bowtie2 -p $btcores --sensitive -x $symprefix -1 $readspref"_first_unm.fq" -2 $readspref"_second_unm.fq" | samtools sort -n -@ 4 > $readspref"_to_"$symprefix".bam"
            btstx=$?
            echo ""; echo ""
        else
            echo "bedtools failed"; echo ""; echo ""
            failed=1
            break 2
        fi

        #create a count file, and add the count file's name to a txt file for use in R
        if [[ $btstx -eq 0 ]]
        then
            echo "counting reads aligning to features..."
            htseq-count -s no -t $5 -i name -q --nonunique none -f bam $readspref"_to_"$symprefix".bam" $4 > $readspref"_counts.txt"
            echo $readspref"_counts.txt" >> counts_list.txt
            hcx=$?
            echo ""; echo ""
        else
            echo "bowtie2/samtools pipe failed"; echo ""; echo ""
            failed=1
            break 2
        fi

        if [[ $hcx -eq 0 ]]
		then
	        echo $readspref" completed"
	        echo ""; echo ""
	    else
	        echo "htseq-count failed"; echo ""; echo ""
		    rm counts_list.txt
            failed=1
            break 2
        fi

        echo "Finished processing read-pair No. "$counter" from "$6

    done < "$input"


    if [[ $failed -eq 0 ]]
    then

        #####write an R script that will combine counts into a matrix
    	touch make_matrix.R

    	#create a chr vector listing the filenames of all count files in the wd, and make a header"
    	echo 'counts = readLines("counts_list.txt")' >> make_matrix.R; echo 'hdr = gsub("_counts.txt", "", counts)' >> make_matrix.R

    	#make a chr vector to be used as row.names"
    	echo 'features = read.delim(counts[1], header = FALSE)[,1]' >> make_matrix.R

    	#bind the columns together, add header, and add row names"(note: vectors must be made into lists)"
    	echo 'df = do.call(cbind,lapply(counts,function(fn)read.table(fn,header=FALSE, sep="\t")[,2]))' >> make_matrix.R
    	echo 'colnames(df) = as.list(hdr)' >> make_matrix.R; echo 'rownames(df) = as.list(as.character(features))' >> make_matrix.R

    	#write the matrix to file."
    	echo 'write.table(df, file = "count_matrix.txt", sep = "\t", quote = FALSE, col.names = NA, row.names = TRUE)' >> make_matrix.R
    	#####R script finished


    	#run the R script
    	module load R
    	chmod 700 make_matrix.R
    	Rscript ./make_matrix.R
    	Rstatus=$?
    	if [ $Rstatus -eq 1 ]
    	then
    		echo "R failed to combine counts"
    		rm make_matrix.R
    	else
    		#move temp files
    		mv *_counts* temp_files
    		mv *.bam temp_files
    		mv *trimmed* temp_files
            mv *unm* temp_files
    		rm make_matrix.R
            rm -r tophat_temp
    		echo 'job completed'
    		echo "Results are in \'count_matrix.txt\'"
    		echo "alignment files and trimmed readfiles are in /'temp_files/' directory"
    	fi


    else
        echo "Job aborted due to non-zero exit status. Refer to STDERR or nohup.out for details"
    fi


else
  echo "An argument is missing"
fi
