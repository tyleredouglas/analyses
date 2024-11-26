# Raw data processing

### fastQC/MultiQC

generate fastQC commands
```
for file in *.gz; do TRIMMEDNAME=$(echo $file | sed 's/_001.fastq.gz//'); printf "#! /bin/bash\n# Job name:\n#SBATCH --job-name=${TRIMMEDNAME}_QC\n#\n# Account:\n#SBATCH --account=fc_poison\n#\n# Partition:\n#SBATCH --partition=savio\n#\n# Wall clock limit:\n#SBATCH --time=04:00:00\n#\n## Command(s) to run:\nmodule load fastqc\nfastqc -o fastqc_output $file" > ${TRIMMEDNAME}.sh; done
```

### trimmomatic

generate trimmomatic commands:
```
for r1 in *_R1_*.fastq.gz; do r2=$(echo $r1 | sed 's/_R1_/_R2_/'); name=$(echo $r1 | sed 's/_R1_001.fastq.gz//' | sed 's/raw_data\///'); printf "#! /bin/bash\n# Job name:\n#SBATCH --job-name=${TRIMMEDNAME}_QC\n#\n# Account:\n#SBATCH --account=fc_poison\n#\n# Partition:\n#SBATCH --partition=savio2\n#\n# Wall clock limit:\n#SBATCH --time=24:00:00\n#\n## Command(s) to run:\nmodule load java\nmodule load trimmomatic\njava -jar /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/trimmomatic/0.36/trimmomatic-0.36.jar PE -threads 20 $r1 $r2 -baseout trimmed_reads/$name.fastq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:4:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20  MINLEN:50 >& trim_logs/$name.log” >${name}_trim.sh; done 
```

actual trimmomatic command:
  
  ```
    #! /bin/bash
    # Job name:
    #SBATCH --job-name=G16R17_S20_trim
    #
    # Account:
    #SBATCH --account=fc_poison
    #
    # Partition:
    #SBATCH --partition=savio2
    #
    # Wall clock limit:
    #SBATCH --time=24:00:00
    #
    ## Command(s) to run:
    module load java
    module load trimmomatic
    java -jar /clusterfs/vector/home/groups/software/sl-7.x86_64/modules/trimmomatic/0.36/trimmomatic-0.36.jar PE -threads 20 G16R17_S20_R1_001.fastq.gz G16R17_S20_R2_001.fastq.gz -baseout trimmed_reads/G16R17_S20.fastq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:4:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50 >& trim_logs/G16R17_S20.log[tylerdouglas@ln003 raw_reads]$ 
  ```
### bowtie2 mapping

generate bowtie2 commands:
```
for p1 in *1P.fastq.gz; do p2=$(echo $p1 | sed 's/1P/2P/'); u1=$(echo $p1 | sed 's/1P/1U/'); u2=$(echo $p1 | sed 's/1P/2U/'); name=$(echo $p1 | sed 's/_1P.fastq.gz//'); printf "#! /bin/bash\n# Job name:\n#SBATCH --job-name=${name}_bt\n#\n# Account:\n#SBATCH --account=fc_poison\n#\n# Partition:\n#SBATCH --partition=savio2\n#\n# Wall clock limit:\n#SBATCH --time=24:00:00\n#\n## Command(s) to run:\nmodule load bowtie2\nbowtie2 --very-sensitive --phred64 -p 20 -x dmel_r6_ref_in -1 trimmed_reads/${p1} -2 trimmed_reads/${p2} -U trimmed_reads/${u1},trimmed_reads/${u2} --al-conc ${name}_P_out.fastq.gz --al ${name}_U_out.fastq.gz -S ${name}.sam 2> mapped_reads/map_logs/${name}.log" > ${name}_bowtie2.sh; done
```
* ran the above, and got error: "Saw ASCII character 58 but expected 64-based Phred qual. Try not specifying --solexa1.3-quals --phred64-quals. terminate called after throwing an instance of 'int'"
* attempted to rerun but having removed -phred64 flag
* command now runs with -phred64 removed, likely quality is actually phred33. odd as last time i did mapping i had a phred64 flag and looking back, that data also seems encoded as phred33 so unclear why that worked previously.

actual bowtie2 command:

```
#! /bin/bash
# Job name:
#SBATCH --job-name=G14R28A_S18_bt
#
# Account:
#SBATCH --account=fc_poison
#
# Partition:
#SBATCH --partition=savio2
#
# Wall clock limit:
#SBATCH --time=24:00:00
#
## Command(s) to run:
module load bowtie2
bowtie2 --very-sensitive -p 20 -x dmel_r6_ref_in -1 trimmed_reads/G14R28A_S18_1P.fastq.gz -2 trimmed_reads/G14R28A_S18_2P.fastq.gz -U trimmed_reads/G14R28A_S18_1U.fastq.gz,trimmed_reads/G14R28A_S18_2U.fastq.gz --al-conc G14R28A_S18_P_out.fastq.gz --al G14R28A_S18_U_out.fastq.gz -S G14R28A_S18.sam 2> mapped_reads/map_logs/G14R28A_S18.log
```


### Re-mapping reads using bwa mem as bowtie2 performance with pool-seq data is not optimal

`for p1 in *1P.fastq.gz; do p2=$(echo $p1 | sed 's/1P/2P/'); name=$(echo $p1 | sed 's/_1P.fastq.gz//'); printf "#! /bin/bash\n# Job name:\n#SBATCH --job-name=${name}_bt\n#\n# Account:\n#SBATCH --account=fc_poison\n#\n# Partition:\n#SBATCH --partition=savio2\n#\n# Wall clock limit:\n#SBATCH --time=24:00:00\n#\n## Command(s) to run:\nmodule load bwa/0.7.17\nmodule load samtools\nbwa mem -t 20 -M dmel_r6_ref.fna.gz trimmed_reads/${p1} trimmed_reads/${p2} | samtools sort -o bwa_mapped/${name}.sorted.bam - 2> bwa_mapped/map_logs/${name}.log; samtools index bwa_mapped/${name}.sorted.bam" > ${name}_bwa.sh; done`

```
#! /bin/bash
# Job name:
#SBATCH --job-name=G13R30_S17_bt
#
# Account:
#SBATCH --account=fc_poison
#
# Partition:
#SBATCH --partition=savio2
#
# Wall clock limit:
#SBATCH --time=24:00:00
#
## Command(s) to run:
module load bwa/0.7.17
module load samtools
bwa mem -t 20 -M dmel_r6_ref.fna.gz trimmed_reads/G13R30_S17_1P.fastq.gz trimmed_reads/G13R30_S17_2P.fastq.gz | samtools sort -o bwa_mapped/G13R30_S17.sorted.bam - 2> bwa_mapped/map_logs/G13R30_S17.log; samtools index bwa_mapped/G13R30_S17.sorted.bam[tylerdouglas@ln001 g1_g20_data]$ 

```



