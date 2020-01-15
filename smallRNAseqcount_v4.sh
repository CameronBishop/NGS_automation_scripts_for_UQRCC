#!/bin/bash
# The script takes raw smallRNAseq files, a reference genome, and an annotation file.
# The script trims the reads on quality and length, generates fastQC reports, maps reads to reference,
# counts reads overlapping annotated miRNAs, and generates a count matrix for use in DESeq2 or similar.
# For use on the UQ RCC HPC (awoonga, tinaroo, etc).
#
# usage: smallRNAseqcountV4.sh <threads> <min_read_length> <max_read_length> <features.gff> <reference.fa> <fastq>

# BEFORE RUNNING:
#
# For this script to work you must first install htseq-count to your home directory "~/.local/bin" because it does not
# exist on the HPC - at least it did not when I wrote this.
# for instructions on downloading and installing htseq (which is a python package), and dependencies see:
# http://www2.rcc.uq.edu.au/hpc/guides/index.html?secure/Installation_of_Packages_for_Scripting_Languages.html#Python-distutils
# and
# https://htseq.readthedocs.io/en/release_0.11.1/install.html
#
# All fastq filenames must include exactly 1 'dot' and this dot should deliniate the prefix from the suffix.
# (e.g. 'reads.fastq' or 'reads.fq' or 'some_reads.fq')
# The suffix (with no dot) will need to be included as the sixth argument so that the script knows what files
# to use as input. As long as the read files are all in the current working directory the script will process them.
#
# The adapter sequence used here is the TruSeq3-SE.fa 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA' from the NEBNext multiplex smRNA kit.
#
# Create a bowtie2 index for the reference file (note that bowtie2 requires a '.fa' suffix for the reference, otherwise it throws an error).
#
# This script is designed to work on the UQ RCC cluster (awoonga, tinaroo, etc). The 'module load' commands will probably only work
# on this system.
#
# The script assumes that the gff file lists the feature type as 'region' in column 3, and that miR names are given as 'ID=*****' in column 9.
# eg: chromosome_1	miRbase	region	310530476	310530497	.	-	.	ID=aae_mir_10_5p'
#
# issue the command 'chmod 700 ./smallRNAseqcountV4.sh' from the dir in which the script is located. (chmod will grant executive privileges).




if [[ -n $1 ]] && [[ -n $2 ]] && [[ -n $3 ]] && [[ -n $4 ]] && [[ -n $5 ]] && [[ -n $6 ]]
then
	  echo "number of processing cores requested: "$1
	  echo "minimum read length: "$2
	  echo "maximum read length: "$3
	  echo "gene feature file: "$4
	  echo "reference genome file: "$5
	  echo "using input read files with "$6" suffix"
	  echo ""; echo ""

	  #add htseq-count to $PATH
	  PATH=$PATH":~/.local/bin"

	  #load env variables, etc
	  module load trimmomatic
	  module load fastqc
	  module load bbmap
	  module load Java
	  module load bowtie2
	  module load samtools
	  module load python

	  #make a bowtie2 prefix
	  idxprefix=$(echo "$5" | cut -f 1 -d '.')

	  #create an array of fastq filenames
	  shopt -s nullglob
	  inreads=(*.$6)

	  #print number of read files detected, and list them.
	  echo "found ${#inreads[@]} read files in working directory"
	  for i in "${inreads[@]}"; do echo "$i"; done

	  #iterate over inreads array
	  for i in "${inreads[@]}"
	  do
		  #make prefix for readfile
		  readspref=$(echo $i | cut -f 1 -d ".")

		  #make filename for filtered fastqs
		  filtreads=$readspref"_"$2"-"$3".fastq"

		  #make filenames for trimmed fastqs
		  trimmed=$readspref".trimmed.fq"

		  #make prefix for sam, bam and txt files
		  aligned=$readspref"_to_miRs"
		  #run the pipeline

		  #Trim reads for adapter sequences
		  java -jar /opt/biotools/trimmomatic/trimmomatic-0.35.jar SE -threads $1 -quiet\
		  $i \
		  $trimmed\
		  ILLUMINACLIP:/opt/biotools/trimmomatic/adapters/TruSeq3-SE.fa:2:30:9
		  trimx=$?
		  echo ""; echo ""

		  #generate QC report for paired reads
		  if [ $trimx -eq 0 ]
		  then
		    echo "running fastqc..."
		    fastqc $trimmed
		    fqcx=$?
		    echo ""; echo ""
		  else
		    echo "trimmomatic failed"
		  fi

		  #filter reads based on size
		  if [[ $trimx -eq 0 ]] && [[ $fqcx -eq 0 ]]
		  then
		  echo "filtering reads bbmap..."
		  reformat.sh in=$trimmed out=$filtreads minlength=$2 maxlength=$3
		  filx=$?
		  echo ""; echo ""
		  else
		      echo "fastQC failed"
		  fi

		  #map reads to reference, convert sam to bam and sort by query name
		  if [[ $trimx -eq 0 ]] && [[ $fqc -eq 0 ]] && [[ $filx -eq 0 ]]
		  then
		    echo "generating sorted bam file..."
		    btcores=$(expr $1 - 1)
		    bowtie2 -p $btcores --very-sensitive-local -x $idxprefix -U $filtreads | samtools view -h -@ 1 -uT $5 > $aligned.bam
		    samtools sort -n -@ $1 $aligned.bam > $aligned.sorted.bam; rm $aligned.bam
		    btstx=$?
		    echo ""; echo ""
		  else
		    echo "bbmap failed"
		  fi

		  #count reads to features
		  if [[ $trimx -eq 0 ]] && [[ $fqc -eq 0 ]] && [[ $filx -eq 0 ]] && [[ $btstx -eq 0 ]]
		  then
		    echo "counting reads aligning to features..."
		    htseq-count -s no -t region -i ID -q --nonunique none -f bam $aligned.sorted.bam $4 > $aligned"_counts.txt"
			echo $aligned"_counts.txt" >> counts_list.txt

		    hcx=$?
		    echo ""; echo ""
		  else
		    echo "bowtie2/samtools pipe failed"
		  fi

		  if [[ $trimx -eq 0 ]] && [[ $fqc -eq 0 ]] && [[ $filx -eq 0 ]] && [[ $btstx -eq 0 ]] && [[ $hcx -eq 0 ]]
		  then
		    echo $readspref" completed"
		    echo ""; echo ""
		  else
		    echo "htseq-count failed"
			rm counts_list.txt
		  fi

	  done

		#####write an R script that will combine counts into a matrix

		touch make_matrix.R

		#create a chr vector listing the filenames of all count files in the wd, and make a header"
		echo 'counts = readLines("counts_list.txt")' >> make_matrix.R
		echo 'hdr = gsub("_counts.txt", "", counts)' >> make_matrix.R

		#make a chr vector to be used as row.names"
		echo 'miRs = read.delim(counts[1], header = FALSE)[,1]' >> make_matrix.R

		#bind the columns together, add header, and add row names"
		#(note: vectors must be made into lists)"
		echo 'df = do.call(cbind,lapply(counts,function(fn)read.table(fn,header=FALSE, sep="\t")[,2]))' >> make_matrix.R
		echo 'colnames(df) = as.list(hdr)' >> make_matrix.R
		echo 'rownames(df) = as.list(as.character(miRs))' >> make_matrix.R

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
			mkdir temp_files
			mkdir output_files
			mv *_counts* temp_files
			mv *.bam temp_files
			mv *trimmed* temp_files
			mv *18-24* temp_files
			mv *.fai temp_files
			mv counts_list.txt output_files
			mv make_matrix.R temp_files
			echo "job completed"
			echo 'Results are in 'output_files' directory'
			echo 'alignment files and trimmed readfiles are in 'temp_files' directory'
		fi

else
  echo "An argument is missing"
fi
