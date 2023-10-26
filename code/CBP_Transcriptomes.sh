#!/usr/bin/bash

##running Trinity v2.15.1 on several species from the Canada Biogenomes Project, starting with arctic surf clam
cd /mnt/sda/CanadaBiogenomes/

#first trim reads with fastp v0.23.2

fastp --in1 SurfClam/RNA/SRR24725107_1.fastq \
--in2 SurfClam/RNA/SRR24725107_2.fastq \
--out1 SurfClam/RNA/SRR24725107_1_trimmed.fastq.gz \
--out2 SurfClam/RNA/SRR24725107_2_trimmed.fastq.gz \
-R 'SurfClam RNA fastp report' --verbose 

#Now run Trinity on our sea lemon trimmed reads to create a de novo assembly
sudo docker run --rm -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq Trinity \
      --seqType fq \
      --left ../../../../mnt/sda/CanadaBiogenomes/SurfClam/RNA/SRR24725107_1_trimmed.fastq.gz \
      --right ../../../../mnt/sda/CanadaBiogenomes/SurfClam/RNA/SRR24725107_2_trimmed.fastq.gz \
      --max_memory 120G --CPU 40 --full_cleanup --output `pwd`/trinitySurfClam

'
Runtime
=======
Start:       Wed Oct 25 16:29:23 UTC 2023
End:         Thu Oct 26 04:11:38 UTC 2023
Trinity   42135 seconds
  Inchworm (phase 1 - read clustering)  2196 seconds
  Chrysalis (phase 1 - read clustering) 33176 seconds
  Rest (phase 2 - parallel assembly)       6763 seconds      
 '     
#now we have a Trinity.fasta transcriptome assembly in our output directory
#The next few commands will be used to assess transcriptome completeness
cd SurfClam/RNA/trinitySurfClam #just make sure you're going to the right directory when all these files are getting made

#1. Assess read content using bowtie
#build a bowtie2 index
bowtie2-build SurfClam.Trinity.fasta SurfClam.Trinity.fasta

#then align reads
bowtie2 -p 40 -q --no-unal -k 20 -x SurfClam.Trinity.fasta \
-1 ../SRR24725107_1_trimmed.fastq.gz \
-2 ../SRR24725107_2_trimmed.fastq.gz \
  2>align.stats.txt | samtools view -@10 -Sb -o SurfClamTransc.bam
  
  cat 2>$1 align_stats.txt
  
  
  #visualize read support with IGV
  samtools sort SurfClamTransc.bam -o SurfClamTrans.coordSorted.bam
  samtools index SurfClamTrans.coordSorted.bam
  samtools faidx SurfClam.Trinity.fasta
  igv.sh -g SurfClam.Trinity.fasta SurfClamTrans.coordSorted.bam
  
#estimating transcript abundance with perl scripts
  sudo docker run -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq \
  /usr/local/bin/trinityrnaseq/util/align_and_estimate_abundance.pl
      
  #get N50 values etc 
  $TRINITY_HOME/util/TrinityStats.pl  Trinity.fasta
  
 #2. Use BUSCO to assess completeness
 #Created a conda environment called SeaLemonBusco when installing newest version of Busco
 conda activate SeaLemonBusco 
 busco -i SurfClam.Trinity.fasta -l metazoa_odb10 -o SurfClamTranscBUSCOout --mode transcriptome
 
conda deactivate
