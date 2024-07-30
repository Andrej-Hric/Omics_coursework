## Andrej Hric 27/07/2023

# Omics Coursework PART 1 detailed analysis of Negative.fq

## Practical Background information
The practical session involving whole-genome sequencing data, the file `Negative.fq`, extracted using specific barcodes, exhibited poor mapping results compared to other samples like `Positive.fq`. The purpose of this coursework was to analyse and explain this poor performance using best practices in NGS data handling and alignment techniques.

## Environment Setup 

For Part 1 of the coursework I first installed Rosetta(as I am using non intel based MAC), conda and all packages required and suggested by the coursework instructions for Part 1.
```bash
softwareupdate --install-rosetta
```

Then I installed Conda, created an environment for this coursework and installed all packages required based on the coursework instructions.

```bash
conda init
conda create --name omics_cw_env python=3.9
conda activate omics_cw_env
conda install -c bioconda cutadapt bowtie2 samtools fastqc multiqc igv
```

I used conda and created new environment to ensure that there are no software dependencies conflict which I found to be issue originally when not using conda. Also having this specific environment I could call it  whenever neccesary in future.

## Issues with original reads and improving of the mapping 

I ran FastQC report on Negative.fq  file which we want to analyze in detail from practical
```bash
cd input_data/p1/fastq
fastqc trimmed_Negative.fq
```
The report from the practical on Negative.fq only showed problems with 1. Per base sequence quality and 2. Per base N content shown bellow. The fastqc report showed issues with bp possitions [1-6] and [56-59] which are the 5' and 3' ends of the reads. The N base content showed that these bases were classified as other than ACTG bases and labelled as N due to low confidence of the base being one of the four main bases.

![](/input_data/p1/fastq/trimmed_Negative_fastqc/Images/per_base_quality.png)

![](/input_data/p1/fastq/trimmed_Negative_fastqc/Images/per_base_n_content.png)

Next I preformed an alignemnt as in the practical with the same settings this was --end-to-end. However this alignment failed as it was unable to find a match to the reference genome. This is why I assumed that this is due to the bases at positions[1-6] and [56-59] as indicated by the fastQC report so the first five on 5'and last four on 3'ends. I then also performed a local alignment instead which could be a solution however the alignment result ended up having 0.00% alignment rate. 

This is why I then performed another alignment however specifying start/end positions, effectively ignoring the the end bp which caused the issue not being able to align with the reference genome.This showed that the new alignment was sucessfully aligned with the reference genome to a 99.97% alignment rate. 

This was quite a significant improvement compared to original 0.00% proving that the alignment problem was casued by the N base read calls.


For all of this I used bowtie2 and commands shown bellow.
```bash
bowtie2 --end-to-end --all -x input/p1/genome/AFPN02.1_merge -q input/p1/fastq/trimmed_Negative.fq -S output/p1/Negative_local.sam >& output/p1/Negative_local_bowtie_output_statistics.txt
```
```bash
bowtie2 --local --all -x input/p1/genome/AFPN02.1_merge -q input/p1/fastq/trimmed_Negative.fq -S output/p1/Negative_local.sam >& output/p1/Negative_local_bowtie_output_statistics.txt
# note : the original alignemnt txt file was overwritten as the results showed no alignment  (can solve this by renaming the output file differently) 

-Results from Negative_local_bowtie_output_statistics.txt file of the these alignments
""" 
    1076320 reads; of these:
        1076320 (100.00%) were unpaired; of these:
            1076320 (100.00%) aligned 0 times
            0 (0.00%) aligned exactly 1 time
            0 (0.00%) aligned >1 times
        0.00% overall alignment rate
""" 
```
```bash
 bowtie2 --end-to-end --all -x input/p1/genome/AFPN02.1_merge -q input/p1/fastq/trimmed_Negative.fq -5 5 -3 4 -S output/p1/Negative_trim.sam >& output/p1/Negative_local_bowtie_output_statistics.txt

-Improved results of specific alignment
""" 
    1076320 reads; of these:
        1076320 (100.00%) were unpaired; of these:
            308 (0.03%) aligned 0 times
            1020237 (94.79%) aligned exactly 1 time
            55775 (5.18%) aligned >1 times
        99.97% overall alignment rate
""" 
```

I then produced a new fastQC report as shown bellow. This Report now has high quality for all positions, compared to the previous fastQC of the untrimmed Negative.fq. Also the summary.txt file provided when unziping the fastQC file showed that all tests were passed. The fastQC report is shown below along with the summary.txt content

![](/input_data/p1/fastq/NEW_base_trimmed_Negative_fastqc/Images/per_base_quality.png)


```bash 
""" 
PASS	Basic Statistics	NEW_base_trimmed_Negative.fq
PASS	Per base sequence quality	NEW_base_trimmed_Negative.fq
PASS	Per sequence quality scores	NEW_base_trimmed_Negative.fq
PASS	Per base sequence content	NEW_base_trimmed_Negative.fq
PASS	Per sequence GC content	NEW_base_trimmed_Negative.fq
PASS	Per base N content	NEW_base_trimmed_Negative.fq
PASS	Sequence Length Distribution	NEW_base_trimmed_Negative.fq
PASS	Sequence Duplication Levels	NEW_base_trimmed_Negative.fq
PASS	Overrepresented sequences	NEW_base_trimmed_Negative.fq
PASS	Adapter Content	NEW_base_trimmed_Negative.fq
""" 

```
With this information I could proceed to a new alignment as all information seemed correct up until this point

The alignment was now done with the new file New_base_trimmed_Negative.fq created in the next step
## Remaping based on N calls being the problem

To replicated the alignment score and prove my hypothesis I produced a new file trimming  5'and 3'ends and also all N calls, effectively to eliminate any issue causing bases.

```bash 
# remove all N calls and first 5 and last 4  possitions
cutadapt --trim-n --trim5 5 --trim3 4 -o input_data/p1/fastq/NEW_base_trimmed_Negative.fq input_data/p1/fastq/trimmed_Negative.fq


# new alignment
bowtie2 --end-to-end --all -x input_data/p1/genomes/AFPN02.1_merge -q input_data/p1/fastq/NEW_base_trimmed_Negative.fq -S output_data/p1/New_Negative_trim.sam >& output_data/p1/New_Negative_trim_bowtie_output_statistics.txt

-Results
""" 
    1076320 reads; of these:
        1076320 (100.00%) were unpaired; of these:
            308 (0.03%) aligned 0 times
            1020237 (94.79%) aligned exactly 1 time
            55775 (5.18%) aligned >1 times
        99.97% overall alignment rate
        """ 
```

## FINAL mapping statistics 

To get the mapping statistics I looked up samtools commands for this.

```bash

#first I converted SAM files to BAM format for more efficient handling and processing as BAM file is smaller and faster and recomented to use for samtool

samtools view -bS output_data/p1/New_Negative_trim.sam | samtools sort -o output_data/p1/New_Negative_trim_sorted.bam 

#then to index the bam file 
samtools index output_data/p1/New_Negative_trim_sorted.bam

#after indexing the bam files, statistics can be extracted and we can look at alignment quality stats.
samtools stats output_data/p1/New_Negative_trim_sorted.bam > output_data/p1/New_Negative_trim_stats.txt

# flagstat analysis used as recommended by the coursework 
samtools flagstat output_data/p1/New_Negative_trim_sorted.bam > output_data/p1/New_Negative_trim_flagstats.txt

# lastly again as recommended I used MultiQC to summarize the statistics and generate a report, this should contain all data from quality control and alignment.
multiqc output_data/p1/ -o output_data/p1/multiqc_report
```
By doing this I could resolve the question wht the original figure in the coursework, of Negative_bowtie_stats and Positive_bowtie_stats were different, and I ended up with exactly the same result of Negative_bowtie_stats as Positive_bowtie_stats in the original Bowtie 2 SE alignemnts scores. This is shown on the image below.

![](/output_data/p1/multiqc_report/bowtie2_se_plot.png)

## Conclusion of coursework P1

In my genomic sequencing analysis, I initially attempted to align the reads using bowtie2 without preprocessing. This direct approach resulted in an alignment rate of 0%, which was quite disappointing. The FastQC analysis had already shown that the extremities of the reads were plagued with ‘N’ bases—ambiguous indicators that represent uncertainty in nucleotide identity. These ‘N’ bases, found primarily at the 5’ and 3’ ends, effectively hindered proper alignment to the reference genome, as they obscured the true biological sequences.

Considering these challenges, I decided to implement a trimming strategy. I trimmed the first 5 bases from the 5’ end and the last 4 bases from the 3’ end of each read, targeting precisely those regions flagged by FastQC. I think this approach was necessary to get rid off the problematic bases which were causing misalignments. This adjustment improved the alignment rate to 99.97%. I consider the trimmed alignment method better and easies as when I removed these problem N sections, the reads aligned closely to the genome, demonstrating the substantial impact that even simple preprocessing steps can have on the quality and utility of NGS data. This experience has reinforced my belief in the critical role of preprocessing in enhancing the accuracy of genomic analyses.

