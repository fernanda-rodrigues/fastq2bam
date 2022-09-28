#!/bin/bash
# ------------------------------------------------------------------
#title                  :fastq2bam.sh
#description            :This script will trim input fastq files for adapter sequences and low quality reads and align such files to the reference genome of interest.
#                        Duplicates will be also marked in the output BAM file.
#                        Intermediate fastq and bam files are deleted after use. The outputs are sorted BAM files with marked duplicates.
#author                 :Fernanda Martins Rodrigues (fernanda@wustl.edu)
#date                   :09272022
#usage                  :bash fastq2bam.sh -1 [R1 fastq file] -2 [R2 fastq file] -s [sample ID] -o [output directory] -C [config.ini file] 
#notes                  :Install samtools, trim-galore, fastQC, bwa and picard to use this script.
#bash_version           :4.2.46(2)-release (x86_64-redhat-linux-gnu)
# ------------------------------------------------------------------

USAGE="Usage: bash fastq2bam.sh -1 [R1 fastq file] -2 [R2 fastq file] -s [sample ID] -o [output directory] -C [config.ini file] "

# --- Options processing -------------------------------------------
if [ "$#" -ne 10 ] ; then
    echo $USAGE
    echo $DETAILS
    exit 1;
fi


while getopts "1:2:s:o:C:" opt; do
    case $opt in
        1)
            FQ1=$OPTARG # input R1 fastq.gz file; please provide full path
            ;;
        2)
            FQ2=$OPTARG # input R2 fastq.gz file; please provide full path
            ;;
        s)
            SAMPLE_ID=$OPTARG # SAMPLE_ID ID which will serve as basename for output files
            ;;
        o)
            OUT_DIRECTORY=$OPTARG # output directory
            ;;
        C)
            CONFIG=$OPTARG
            ;;
        \?)
            echo "Invalid option -$OPTARG"
            echo $USAGE
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument."
            echo $USAGE
            exit 1
            ;;
    esac
done

if [ ! -d ${OUT_DIRECTORY} ]; then
    mkdir -p "${OUT_DIRECTORY}"
fi

if [ ${OUT_DIRECTORY: -1} != "/" ]; then
    OUT_DIRECTORY=${OUT_DIRECTORY}/
fi


# --- Body --------------------------------------------------------

source ${CONFIG}

# STEP 1: TRIM INPUT FASTQs USING TRIM GALORE

echo "STEP 1: TRIMMING ..."
$TRIMGALORE -j 8 --phred33 --fastqc --length 36 -q 20 --no_report_file -o ${OUT_DIRECTORY} --basename ${SAMPLE_ID} --paired ${FQ1} ${FQ2}

# STEP 2: MAP

echo "STEP 2: ALIGNING ..."
$BWA mem -t 8 -M -R "@RG\tID:$SAMPLE_ID\tPL:illumina\tLB:$SAMPLE_ID\tPU:$SAMPLE_ID\tSM:$SAMPLE_ID" ${GENOME} ${OUT_DIRECTORY}${SAMPLE_ID}_val_1.fq.gz ${OUT_DIRECTORY}${SAMPLE_ID}_val_2.fq.gz | samtools view -Shb -o ${OUT_DIRECTORY}${SAMPLE_ID}.bam -

# STEP 3: SORT BAM

if [ -f '${OUT_DIRECTORY}${SAMPLE_ID}.bam' ]
then
    rm -I ${OUT_DIRECTORY}${SAMPLE_ID}_val_1.fq.gz
    rm -I ${OUT_DIRECTORY}${SAMPLE_ID}_val_2.fq.gz
fi

echo "STEP 3: SORTING BAM FILE ..."
$SAMTOOLS sort ${OUT_DIRECTORY}${SAMPLE_ID}.bam -o ${OUT_DIRECTORY}${SAMPLE_ID}.sorted.bam

# STEP 4: MARK DUPLICATES

if [ -f '${OUT_DIRECTORY}${SAMPLE_ID}.sorted.bam' ]
then
    rm -I ${OUT_DIRECTORY}${SAMPLE_ID}.bam
fi

echo "STEP 4: MARKING DUPLICATES ..."
$PICARD MarkDuplicates I=${OUT_DIRECTORY}${SAMPLE_ID}.sorted.bam O=${OUT_DIRECTORY}${SAMPLE_ID}.sorted.markedDup.bam M=${OUT_DIRECTORY}${SAMPLE_ID}.marked_dup_metrics.txt VALIDATION_STRINGENCY=STRICT CREATE_MD5_FILE=true

# STEP 5: INDEX FINAL BAM

if [ -f '${OUT_DIRECTORY}${SAMPLE_ID}.sorted.markedDup.bam' ]
then
    rm -I ${OUT_DIRECTORY}${SAMPLE_ID}.sorted.bam
fi

echo "STEP 5: INDEXING FINAL BAM ..."
$SAMTOOLS index ${OUT_DIRECTORY}${SAMPLE_ID}.sorted.markedDup.bam

echo "DONE!"

# --- End --------------------------------------------------------
