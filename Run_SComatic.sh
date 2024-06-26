conda activate SComatic
myPath=/media/scPRAD
OUTdir=$myPath/result/OUT_SComatic
cd $myPath
SCOMATIC=SComatic
REF=reference/genome/hg38/hg38.fa
editing=$SCOMATIC/RNAediting/AllEditingSites.hg38.txt
PON=$SCOMATIC/PoNs/PoN.scRNAseq.hg38.tsv

mkdir -p $OUTdir
cd $OUTdir

if [ -f $myPath/SComatic.list ]; then
  rm $myPath/SComatic.list
  ls -d $myPath/data/*_cellranger710/ >> $myPath/SComatic.list
else
  ls -d $myPath/data/*_cellranger710/ >> $myPath/SComatic.list
fi

cat $myPath/SComatic.list | while read id
do
SampleName=$(echo "$id" | cut -d'/' -f7)
echo $SampleName
cellranger_outDir=$id/outs
BAM=${cellranger_outDir}/possorted_genome_bam.bam
OUT=$OUTdir/$SampleName
mkdir -p $OUT
META=$(ls /media/yuansh/BANQ/data/scRNA/$SampleName/*.tsv)
echo "bash /media/scPRAD/script/SComatic.sh -S" $SCOMATIC "-s" $SampleName "-b" $BAM "-m" $META "-r" $REF "-e" $editing "-p" $PON "-o" $OUT >> mySComatic.sh
done

cat mySComatic.sh
