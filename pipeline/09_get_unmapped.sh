#!/bin/bash
#SBATCH -N 1 -n 4 --mem 2gb --out logs/unmapped.log --time 2:00:00

module load samtools/1.11
module load picard
module load gatk/4
module load java/13

MEM=2g
TMPOUTDIR=tmp
if [ -f config.txt ]; then
  source config.txt
fi
if [ -z $REFGENOME ]; then
  echo "NEED A REFGENOME - set in config.txt and make sure 00_index.sh is run"
  exit
fi
mkdir -p $UNMAPPED

CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi


IFS=,
tail -n +2 $SAMPFILE | while read STRAIN FILEBASE
do
  # BEGIN THIS PART IS PROBABLY PROJECT SPECIFIC
  # THIS COULD NEED TO BE CHANGED TO R1 R2 or R1_001 and R2_001 etc
  PREFIX=$STRAIN
  # END THIS PART IS PROBABLY PROJECT SPECIFIC

  FINALFILE=$ALNFOLDER/$STRAIN.$HTCEXT

  FQ=$(basename $FASTQEXT .gz)
  UMAP=$UNMAPPED/${STRAIN}.$FQ
  UMAPF=$UNMAPPED/${STRAIN}.F.$FQ
  UMAPR=$UNMAPPED/${STRAIN}.R.$FQ

  UMAPSINGLE=$UNMAPPED/${STRAIN}_single.$FQ
  #echo "$UMAP $UMAPSINGLE $FQ"

  if [ ! -f $UMAP.gz ]; then
    module load BBMap
    samtools fastq -f 4 --threads $CPU -N -1 $UMAPF -2 $UMAPR $FINALFILE > $UMAPSINGLE
    if [ -s $UMAP_SINGLE ]; then
    	pigz $UMAPSINGLE
    fi
    repair.sh in=$UMAPF in2=$UMAPR out=$UMAP.gz
    unlink $UMAPF
    unlink $UMAPR
  fi
done

SAMPFILE=samples_single.csv

tail -n +2 $SAMPFILE | while read STRAIN FILEBASE
do
  # BEGIN THIS PART IS PROBABLY PROJECT SPECIFIC
  # THIS COULD NEED TO BE CHANGED TO R1 R2 or R1_001 and R2_001 etc
  PREFIX=$STRAIN
  # END THIS PART IS PROBABLY PROJECT SPECIFIC
  echo "STRAIN is $STRAIN $PAIR1 $PAIR2"

  FINALFILE=$ALNFOLDER/$STRAIN.$HTCEXT

  FQ=$(basename $FASTQEXT .gz)
  UMAP=$UNMAPPED/${STRAIN}.$FQ
  UMAPF=$UNMAPPED/${STRAIN}.F.$FQ
  UMAPR=$UNMAPPED/${STRAIN}.R.$FQ

  UMAPSINGLE=$UNMAPPED/${STRAIN}_single.$FQ
  #echo "$UMAP $UMAPSINGLE $FQ"

  if [ ! -f $UMAP.gz ]; then
    module load BBMap
    samtools fastq -f 4 --threads $CPU -N -1 $UMAPF -2 $UMAPR $FINALFILE > $UMAPSINGLE
    if [ -s $UMAP_SINGLE ]; then
        pigz $UMAPSINGLE
    fi
    repair.sh in=$UMAPF in2=$UMAPR out=$UMAP.gz
    unlink $UMAPF
    unlink $UMAPR
  fi
done
