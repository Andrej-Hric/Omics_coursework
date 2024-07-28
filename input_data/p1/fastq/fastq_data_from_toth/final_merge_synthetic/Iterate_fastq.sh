#!/bin/bash

# Assign everyting ended in .fq to FILES variable
FILES=*.fq

# Loop over FILES and assing each to file
for file in $FILES
do 
    filename=$(basename "$file")
    extension="${filename##*.}"
    filename="${filename%.*}"

    echo "File on the loop: 	      	$file"	
    echo "Curated base name of file:  	$filename"
    echo "Extension:                  	$extension"

    # Call FastQC quality analysis	
    /s/software/fastqc/v0.11.8/FastQC/fastqc ${file}

    echo -e "######################\n\n"
done

# Load MultiQC
module use -s /s/mm/modules
module load python/v2

# Run MultiQC 
# -f 	Overwrite any existing reports
#  . 	on this directory
echo "Running MultiQC..."
echo "Capturing fastQC statistics"
multiqc -f .


