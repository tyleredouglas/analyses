# Variant calling on populations, generating table of haplotype frequencies for downstream analysis

### First quantify coverage of populations samples, did so by sampling coverage at every 1000th position and writing to csv for each population
```
for file in *.sorted.bam; do name=$(echo $file | sed 's/.sorted.bam//'); samtools depth $file | awk 'NR % 1000 == 0' > ${name}.cov; done
```

### add population names to read depth files and make density plots

` 
awk -v OFS=", " '
    NR == 1 {print $0, "filename"}
    FNR == 1 {next}
    {printf "%s\t%s\n", $0, FILENAME}
' *.cov > all_pops.cov
`

### Remove duplicates and re-quantify coverage (samtools markdup). Requires bam to first be sorted by read names, then fill in mate coordinates and insert sizes, then sort by coordinates. Testing for on one population (G13R26)

```
samtools sort -n -o G13R26_namesort.bam G13R26_S58.sorted.bam
samtools fixmate -m G13R26_namesort.bam G13R26_S58.fix.bam
samtools sort -o G13R26_positsort.bam G13R26_S58.fix.bam
samtools markdup -r -s -f test_dup.txt G13R26_positsort.bam G13R26_nodup.bam  
```
### The above worked, so creating slurm script for each population for job submission:

```
for file in *.sorted.bam; do name=$(echo $file | sed 's/.sorted.bam//'); printf "#! /bin/bash\n# Job name:\n#SBATCH --job-name=${name}_dup\n#\n# Account:\n#SBATCH --account=fc_poison\n#\n# Partition:\n#SBATCH --partition=savio2\n#\n#Wall clock limit:\n#SBATCH --time=24:00:00\n#\n## Command(s) to run:\nmodule load samtools\nsamtools sort -@ 20 -n -o ${name}_namesort.bam $file\nsamtools fixmate -@ 20 -m ${name}_namesort.bam ${name}.fix.bam\nsamtools sort -@ 20 -o ${name}_positsort.bam ${name}.fix.bam\nsamtools markdup -@ 20-r -s -f ${name}_dup.log ${name}_positsort.bam ${name}_nodup.bam" > ${name}_dup.sh; done
```

### Downsample 1000x R14G1 and R14G16 files to 100x to include in analysis with other populations

```
samtools view -@ 20 -s 0.1 -b G16R14_S72_nodup.bam > G16_R14_100X_nodup.bam  
samtools view -@ 20 -s 0.1 -b G1R14_S65_nodup.bam > G1_R14_100X_nodup.bam  
```

### re-index bams after duplicate removal
```
for file in *.bam; do samtools index -@ 20 $file; donea
```

### create list of all indexed population and founder bams for mpileup input:

```
for file in *.bam; do printf "$file\n">> bam_names.txt; done
```

### Run mpileup by chromosome:

NC_004354.4 = X
NT_033779.5 = 2L
NT_033778.4 = 2R
NT_037436.4 = 3L
NT_033777.3 = 3R
NC_004353.4 = 4
NC_024512.1 = Y

```
#! /bin/bash
# Job name:
#SBATCH --job-name=chrX_mpile
#
# Account:
#SBATCH --account=fc_poison
#
# Partition:
#SBATCH --partition=savio2
#
# Wall clock limit:
#SBATCH --time=72:00:00
#
## Command(s) to run:
module load bcftools
bcftools mpileup -I -d 1000 --threads 20 -t NC_004354.4 -a "FORMAT/AD,FORMAT/DP" -f dmel_r6.fna -b bam_names.txt -o chrX_mpile.txt
```

### Call SNPs from mpileup files

```
#! /bin/bash
# Job name:
#SBATCH --job-name=chr2L_snps
#
# Account:
#SBATCH --account=fc_poison
#
# Partition:
#SBATCH --partition=savio2
#
# Wall clock limit:
#SBATCH --time=48:00:00
#
## Command(s) to run:
module load perl
module load bcftools
bcftools call chr2L_mpile.txt -mv -Ob > chr2L_calls.bcf

```

### Filtering spurious SNPs, identified by any SNPs that are heterozygous in founders

```
for file in *.vcf; do name=$(echo $file | sed 's/_calls.vcf//'); cat $file | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24}' | grep -e 0/1 -e 1/0 | sort >> ${name}_het.vcf; awk 'FNR==NR{map[$2]=1;next;}map[$2]==""{print;}' ${name}_het.vcf $file >> ${name}_snps.vcf; done
```

*validated filtering by running grep '0/1' on founder columns of ${name}_snps.vcf files, all came up blank as expected (no heterozygotes)

### Create table of SNP counts and frequencies

```
bcftools query -e 'GT ="./."' -e'QUAL<60' -f'%CHROM %POS %REF %ALT [ %AD{0} %AD{1}] [%GT]\n' chr2L_calls.bcf | sed 's/NT_033779.5/chr2L/' | grep -v '\.' | perl /global/home/users/tylerdouglas/modules/tdlong/accuracy.freqtab.pl > temp.chr2L.txt

bcftools query -e 'GT ="./."' -e'QUAL<60' -f '%CHROM %POS %REF %ALT [ %AD{0} %AD{1}] [%GT]\n' chr2L_calls.bcf | sed 's/NT_033779.5/chr2L/' |grep -v '\.' | perl /global/home/users/tylerdouglas/modules/tdlong/accuracy.counttab.pl > temp.count.chr2L.txt
```

*Ran the above and got this error for chr3L: "[W::bcf_sr_add_reader] No BGZF EOF marker; file 'chr3L_calls.bcf' may be truncated
[E::bgzf_read] Read block operation failed with error -1 after 728 of 1365 bytes"

*And this one for chr3R: "Number of columns at NT_033777.3:26961056 does not match the number of samples (86 vs 123)"

*mpileup ran without throwing an error but it seems to have terminated the jobs for chr3L/chr3R, both ended within a second of each other and the last line on each mpileup for shows a position a few million bp  short of the entire contig. this happened despite not hitting the wallclock limit. validated that chr2L, 2R, and X all ran completely. will try to rerun mpileup on only the outstanding pieces of 3L and 3R

*to that end, deleted last line (truncated line) of 3R and 3L mpileups `sed -i '$d' chr3L_mpile_fixed.txt`. 3L got up to 26595417, 3R got up to 26961055

*rerun of mpileup on chr3R and chr3L ran without error

### Generating haplotype allele frequencies
```
bcftools query -e 'GT ="./."' -e'QUAL<60' -f'%CHROM %POS %REF %ALT [ %AD{0} %AD{1}] [%GT]\n' chr2L_calls.bcf | sed 's/NT_033779.5/chr2L/' | grep -v '\.' | perl /global/home/users/tylerdouglas/modules/tdlong/accuracy.freqtab.pl > temp.chr2L.txt
```
*generated blank files for all chromosomes except for chr2L. 
*BCF files were generated successfully, converted to vcf and checked each one. 
*contig names were just incorrect in the above code, went back and fixed them

### Create table of all SNP frequencies and counts for bam files in mpileup:
```
`for line in $(cat pool_names.txt); do population=${line%_S*}; printf "N_$population\t" >> SNP.count.txt; done`
`sed -i 's/.bam//g' SNP.count.txt `
`sed -i 's/_100X_nodup//g' SNP.count.txt`
`cat temp.count.chrX.txt >> SNP.count.txt`
```
*when creating table in R with `temp = read.table("SNP.freq.txt", header=TRUE, as.is=TRUE, strip.white=TRUE)` got error "Error in scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  : 
  line 1 did not have 248 elements"

*running `awk '{a[NF]++}END{for(k in a)print k,a[k]}' SNP.count.txt | sort > elements.txt` shows one line in file has 248 elements instead of expected 123
*errant line is at chrX 131074, seems to have recopied population names which should only appear on line one
*that is actually just the first line of the text file, needs new line character after population names or the snp counts/freqs will be appended to line 1 of the file (123 names + 123 values + chrom name + position = 248), running successfully with change


