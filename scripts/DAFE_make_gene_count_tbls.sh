#! /usr/bin/env bash

# this script make the count table foreach genome based on it's genes

# some default values
COUNT_FILE_NAME="htseq.count.txt"
OUT_FILE_NAME="gene_counts_id60.txt"

# usage and parameters
usage="$(basename "$0") [-h] [-b bam_dir] [-i isolate_names] [-s sample_names] [-c count_file_name] [-o out_file_name]

where:
	-h	show this help text
	-b	bam file directory where the input bam files are stored
	-i	a file with all the isolate names (must correspond to the dir names)
	-s	a file with all the sample names (must correspond to the dir names)
	-c	htseq count file name (not full path just name of file with gene counts. DEFAULT: htseq.count.txt)
	-o	out file name (not full path just name of file to put output table. DEFAULT: gene_counts_id60.txt)"

while getopts ':hb:i:s:c:o:' option; do
	case "$option" in
		h)	echo "$usage"
			exit
			;;
		b)	bam_dir=$OPTARG
			;;
		i)	iso_names=$OPTARG
			;;
		s)	sam_names=$OPTARG
			;;
		c)	count_file_name=$OPTARG
			;;
		o)	out_file_name=$OPTARG
			;;
		:)	printf "missing argument for -%s\n" "$OPTARG" >&2
			echo "$usage" >&2
			exit 1
			;;
		\?)	printf "illegal option: -%s\n" "$OPTARG" >&2
			echo "$usage" >&2
			exit 1
			;;
	esac
done

# set the default values if a variable is undef
# use a : to prevent bash from attempting to execute the value of 
# variable when it is returned
: ${count_file_name=$COUNT_FILE_NAME}
: ${out_file_name=$OUT_FILE_NAME}

# MAIN
for i_name in `cat $iso_names`
do
	# make the output file for each genome
	out_file=$bam_dir/$i_name/$out_file_name
	touch $out_file

	# find a file that has names (ie non-zero size)
	first_s_name='' # reset this to make sure it is empty
	for s_name in `cat $sam_names`
	do
		if [[ -s $bam_dir/$i_name/$s_name/$count_file_name ]]; then
			first_s_name=$s_name
			break
		fi
	done
	
	# what if we don't find a file that has any data (ie no reads map to this genome in any sample
	if [[ -z $first_s_name ]]; then
		echo "WARNING: genome $i_name has no mapped reads in any sample" >&2
		continue
	fi

	# add the names column from the file that I found with the names
	awk '{print $1}' $bam_dir/$i_name/$first_s_name/$count_file_name > $out_file

	# get the column count for a file that is not empty
	col_count=`wc -l $bam_dir/$i_name/$first_s_name/$count_file_name | awk '{print $1}'`;

	# get the columns and print the to the output file
	for s_name in `cat $sam_names`
	do
		# some files may be empty meaning they have no reads that map to them
		# I need to make a column of 0's as long as the other files
		if [[ ! -s $bam_dir/$i_name/$s_name/$count_file_name ]]; then
			for ((i=0; i<$col_count; i++)); do
				echo -e "feature\t0" >> $bam_dir/$i_name/$s_name/$count_file_name
			done
		fi
		
		# update the output files by creating a tmp file and then moving it to the output file
		awk '{print $2}' $bam_dir/$i_name/$s_name/$count_file_name | paste $out_file - > tmp.txt;
		mv tmp.txt $out_file;
	done


	# add the column headers
	echo "name" > .tmp.names.txt
	cat $sam_names >> .tmp.names.txt
	paste -s .tmp.names.txt > tmp.txt
	cat tmp.txt $out_file > tmp2.txt
	mv tmp2.txt $out_file
	rm .tmp.names.txt
	rm tmp.txt

done
