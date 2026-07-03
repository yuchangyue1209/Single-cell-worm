scrna

#run star without gtf and assembly
mkdir -p /work/cyu/scRNA/star_index
14 large recommend 13


#tried   --genomeSAindexNbases 13
mkdir -p /work/cyu/scRNA/star_index_sa13

STAR \
  --runThreadN 48 \
  --runMode genomeGenerate \
  --genomeDir /work/cyu/scRNA/star_index_sa13 \
  --genomeFastaFiles /work/shared/cyu/worm/worm.no_rDNA.fa \
  --genomeSAindexNbases 13

#cr_whitelists
3M-february-2018.txt.gz

#gtf convert
gffread /work/cyu/scRNA/ssol-lift.gff_polished \
  -T \
  -o /work/cyu/scRNA/ssol-lift.gtf

#index
rm -rf /work/cyu/scRNA/star_index_sa13_gtf
mkdir -p /work/cyu/scRNA/star_index_sa13_gtf

STAR \
  --runThreadN 48 \
  --runMode genomeGenerate \
  --genomeDir /work/cyu/scRNA/star_index_sa13_gtf \
  --genomeFastaFiles /work/shared/cyu/worm/worm.no_rDNA.fa \
  --sjdbGTFfile /work/cyu/scRNA/ssol-lift.gtf \
  --sjdbOverhang 149 \
  --genomeSAindexNbases 13

#run star for big worm

rm -rf /work/cyu/tmp_big
rm -rf /work/cyu/scRNA/starsolo_big
mkdir -p /work/cyu/scRNA/starsolo_big

STAR \
  --genomeDir /work/cyu/scRNA/star_index_sa13_gtf \
  --readFilesIn /mnt/spareHD_2/scRNA/Big_S1_L005_R2_001.fastq.gz /mnt/spareHD_2/scRNA/Big_S1_L005_R1_001.fastq.gz \
  --readFilesCommand zcat \
  --runThreadN 32 \
  --outFileNamePrefix /work/cyu/scRNA/starsolo_big/Big_ \
  --outTmpDir /work/cyu/tmp_big \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 \
  --soloCBlen 16 \
  --soloUMIstart 17 \
  --soloUMIlen 12 \
  --soloCBwhitelist /work/cyu/scRNA/cr_whitelists/3M-february-2018.txt \
  --soloFeatures Gene GeneFull \
  --soloUMIdedup 1MM_All \
  --soloCellFilter EmptyDrops_CR \
  --outFilterScoreMin 30 \
  --outSAMtype None

#whitelist error run cellranger v 10.0.0 instead 
cd /work/cyu/scRNA/
cellranger mkref \
  --genome=ssol \
  --fasta=/work/shared/cyu/worm/worm.no_rDNA.fa \
  --genes=/work/cyu/scRNA/ssol-lift.gtf

#run big
cellranger count \
  --id=Big \
  --transcriptome=/work/cyu/scRNA/ssol \
  --fastqs=/mnt/spareHD_2/scRNA \
  --sample=Big \
  --localcores=32 \
  --localmem=200

  