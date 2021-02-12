#!/bin/bash
#SBATCH -N 1 -n 16 --mem 32gb --out logs/bwa_single.%a.log -p intel --time 64:00:00
module load bwa
module load samtools/1.11
module load picard
module load gatk/4
module load java/13

MEM=32g
TMPOUTDIR=tmp
if [ -f config.txt ]; then
  source config.txt
fi
SAMPFILE=samples_single.csv
if [ -z $REFGENOME ]; then
  echo "NEED A REFGENOME - set in config.txt and make sure 00_index.sh is run"
  exit
fi

if [ ! -f $REFGENOME.dict ]; then
  echo "NEED a $REFGENOME.dict - make sure 00_index.sh is run"
fi
mkdir -p $TMPOUTDIR $ALNFOLDER

CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi
N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
  N=$1
fi
if [ -z $N ]; then
  echo "cannot run without a number provided either cmdline or --array in sbatch"
  exit
fi

MAX=$(wc -l $SAMPFILE | awk '{print $1}')
if [ $N -gt $MAX ]; then
  echo "$N is too big, only $MAX lines in $SAMPFILE"
  exit
fi

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read STRAIN FILEBASE ISOLATE BioProject BIOSAMPLE LAT_LONG
do
  PREFIX=$STRAIN
  FINALFILE=$ALNFOLDER/$STRAIN.$HTCEXT
  if [ ! -s $FINALFILE ]; then
    BAMSTOMERGE=()
    for BASE in $(echo $FILEBASE | perl -p -e 's/;/,/')
    do
      # BEGIN THIS PART IS PROBABLY PROJECT SPECIFIC
      # THIS COULD NEED TO BE CHANGED TO R1 R2 or R1_001 and R2_001 etc
      PAIR1=$FASTQFOLDER/${BASE}.$FASTQEXT
      # END THIS PART IS PROBABLY PROJECT SPECIFIC
      echo "STRAIN is $STRAIN $PAIR1"
      TMPBAMFILE=$TEMP/$BASE.unsrt.bam
      SRTED=$TMPOUTDIR/$BASE.srt.bam
      DDFILE=$TMPOUTDIR/$BASE.DD.bam

      READGROUP="@RG\tID:$BASE\tSM:$STRAIN\tLB:$BASE\tPL:illumina\tCN:$RGCENTER"

      if [ ! -s $DDFILE ]; then
        if [ ! -s $SRTED ]; then
          if [ -e $PAIR1 ]; then
            if [ ! -f $SRTED ]; then
              # potential switch this to bwa-mem2 for extra speed
              bwa mem -t $CPU -R $READGROUP $REFGENOME $PAIR1 | samtools sort --threads $CPU -O bam -o $SRTED -T $TEMP -
            fi
          else
            echo "Cannot find $PAIR1, skipping $STRAIN"
            exit
          fi
        fi # SRTED file exists or was created by this block

        time java -jar $PICARD MarkDuplicates I=$SRTED O=$DDFILE \
          METRICS_FILE=logs/$STRAIN.dedup.metrics CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
        if [ -f $DDFILE ]; then
          rm -f $SRTED
        fi
      fi # DDFILE is created after this or already exists
      BAMSTOMERGE+=( $DDFILE )
    done
    samtools merge --write-index -O $HTCFORMAT --threads $CPU --reference $REFGENOME $FINALFILE "${BAMSTOMERGE[@]}"
    if [ -f $FINALFILE ]; then
      rm -f "${BAMSTOMERGE[@]}"
      rm -f $(echo "${BAMSTOMERGE[@]}" | sed 's/bam$/bai/')
    fi
  fi #FINALFILE created or already exists
done
