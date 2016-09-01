# Anole X/A diversity analyses from Rupp et al.

This repository contains scripts and information related to the genetic diversity analyses in Rupp et al (_In review_). Evolution of dosage compensation in _Anolis carolinensis_, a reptile with XX/XY chromosomal sex determination.  Scripts related to other parts of the study (e.g. identifying X-linked scaffolds, differential expression analyses, and Ka/Ks calculations) can be found in [another repository](https://github.com/WilsonSayresLab/Anole_expression).


## QUICK START

This section describes the (more or less) push-button replication of the transciptome assembly, variant calling, and diversity analyses from this paper using [snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home).  Each step in the pipeline is described in greater detail below.

To run the pipeline (has only been tested on Mac and Linux machines):

1) Clone and then enter this repository:
```
git clone https://github.com/thw17/Rupp_etal_Anole_DosageCompensation
cd Rupp_etal_Anole_DosageCompensation
```
This repository contains all scripts necessary to reproduce our analyses, as well as the directory structure for the snakemake pipeline.

2) Set up Anaconda environments (one with Python 3 for snakemake and one with Python 2 for the diversity script).  If you don't already have Anaconda installed, it can be obtained free [from here](https://www.continuum.io/downloads) and you can find more information [here](http://conda.pydata.org/docs/index.html).  You can alternatively install [Miniconda](http://conda.pydata.org/docs/install/quick.html), a lightweight version of Anaconda.  The following commands assume that anaconda has been successfully installed and is in your PATH (it will do this automatically if you allow it):
```
conda config --add channels bioconda
conda env create -f env/anole_dosage.yml
conda env create -f env/diversity_script.yml
```

3) Activate the anole_dosage anaconda environment:
```
source activate anole_dosage
```

4) Download the fastq files from SRA into the "fastqs" directory. A straightforward, but somewhat inefficient way to obtain and compress all of the fastq files would be:
```
cd fastqs
for i in SRR1502164 SRR1502165 SRR1502166 SRR1502167 SRR1502168 SRR1502169 SRR1502170 SRR1502171 SRR1502172 SRR1502173 SRR1502174 SRR1502175 SRR1502176 SRR1502177 SRR1502178 SRR1502179 SRR1502180 SRR1502181 SRR1502182 SRR1502183
do
fastq-dump --gzip --outdir fastqs/ --readids --split-files $i
done
```
This can be sped up significantly by running indepentdent, parallel jobs with 1-3 ids each.

5) Download the reference genome into the "reference" directory and create relvant dictionaries and indices.  Note that STAR might require quite a bit of memory to create the reference index (picard and samtools, on the otherhand, will not).

```
cd reference
wget ftp://ftp.ensembl.org/pub/release-85/fasta/anolis_carolinensis/dna/Anolis_carolinensis.AnoCar2.0.dna.toplevel.fa.gz
gunzip AnoCar2.0.fa.gz
STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {path/to/reference} --genomeFastaFiles AnoCar2.0.fa
samtools faidx AnoCar2.0.fa
picard CreateSequenceDictionary R=AnoCar2.0.fa O=AnoCar2.0.dict
```

6) Edit anoles.config.json with the path to your GATK (if you haven't downloaded it, [you can here](https://software.broadinstitute.org/gatk/download/) .  We used version 3.5 - other versions will give very slightly different results), Snpsift (you can download it [here](http://snpeff.sourceforge.net/) ; we used version 4.2), and where you'd like temporary files to go. 

7) Once all of the fastq files have successfull downloaded, you can run the rest of the pipeline by typing:
```
snakemake -s snakefile -c <number of cores>
```
Snakemake will distribute jobs across nodes. So to speed up the process, have a look at the [documentation](https://bitbucket.org/snakemake/snakemake/wiki/Documentation) and the [tutorial](http://snakemake.bitbucket.org/snakemake-tutorial.html).

For example, our HPC at ASU uses sbatch/slurm and an example command for distributing the pipeline across multiple jobs looks something like:
```
snakemake --snakefile snakefile -j <max number of parallel jobs snakemake can submit> --cluster "sbatch -n <number of cores> -t <time limit> --mail-type=END,FAIL --mail-user=<email address> " --cores <number of cores for the main snakemake process>
```

##DETAILED WALKTHROUGH OF ANALYSES

##Data and samples

The samples for this study consist of RNA-seq data from five biological replicates of regenerating tail tissue each from two male (d15 and d26) and two female (d47 and d61) green anoles (_Anolis carolinensis_) from [Hutchins et al (2014)](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0105004).  Fastq files associated with this study have been deposited in the SRA under project accession [PRJNA253971](http://www.ncbi.nlm.nih.gov/bioproject/PRJNA253971). Further details on the collection, sequencing, and processing of data can be found in [Hutchins et al (2014)](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0105004)

##Obtaining fastq files from SRA

All fastq files are available on SRA (associated with the Hutchins et al. 2014 project; see above) and can be easily downloaded using the [fastq-dump](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump) tool, which is part of the [SRA Toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software).

Basic syntax is something like:
```fastq-dump --splitfiles <sra accession number>```
where <sra accession number> is the short read archive number associated with a sample/run.  [Offical documentation](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump) is a bit minimal, but the Edwards lab put together detailed unoffical documentation [here](https://edwards.sdsu.edu/research/fastq-dump/).

The numbers used in this study, with corresponding sample IDs are:

| Accession | ID |
|---|---|
| SRR1502164 | d15s1 |
| SRR1502165 | d15s2 |
| SRR1502166 | d15s3 |
| SRR1502167 | d15s4 |
| SRR1502168 | d15s5 |
| SRR1502169 | d26s1 |
| SRR1502170 | d26s2 |
| SRR1502171 | d26s3 |
| SRR1502172 | d26s4 |
| SRR1502173 | d26s5 |
| SRR1502174 | d47s1 |
| SRR1502175 | d47s2 |
| SRR1502176 | d47s3 |
| SRR1502177 | d47s4 |
| SRR1502178 | d47s5 |
| SRR1502179 | d61s1 |
| SRR1502180 | d61s2 |
| SRR1502181 | d61s3 |
| SRR1502182 | d61s4 |
| SRR1502183 | d61s5 |

In the sample IDs, the d number is the individual ID and the s number is the biological replicate.  For example, d15s1 corresponds to biological replicate 1 of individual d15.

A straightforward, but somewhat inefficient way to obtain and compress all of the fastq files might be:
```
for i in SRR1502164 SRR1502165 SRR1502166 SRR1502167 SRR1502168 SRR1502169 SRR1502170 SRR1502171 SRR1502172 SRR1502173 SRR1502174 SRR1502175 SRR1502176 SRR1502177 SRR1502178 SRR1502179 SRR1502180 SRR1502181 SRR1502182 SRR1502183
do
fastq-dump --gzip --outdir <path to output directory> --readids --split-files $i
done
```
At about 30 minutes per accession (time it took on my machine), the above command would take somewhere around 10 hours to complete.  Alternatively, this can be sped up by splitting the command up into multiple scripts to be run in parallel (e.g., with 2-5 accession numbers at a time).

##Transcriptome Assembly and Variant Calling

Transcriptome assembly and variant calling more or less followed the [GATK Best Practices workflow for RNA-seq data](https://www.broadinstitute.org/gatk/guide/article?id=3891) using the following tools:

[STAR](https://github.com/alexdobin/STAR)

[samtools](http://www.htslib.org/)

[Picard](https://broadinstitute.github.io/picard/)

[Genome Analysis Toolkit (GATK)](https://www.broadinstitute.org/gatk/)

####First pass read mapping with STAR
We used STAR's 2-pass method to map reads.  Assuming the genome ([AnoCar2](ftp://ftp.ensembl.org/pub/release-85/fasta/anolis_carolinensis/dna/Anolis_carolinensis.AnoCar2.0.dna.toplevel.fa.gz)) has been downloaded and properly indexed (see QUICK START above for more information about this), the first pass with STAR is relatively straightforward.  We used the following template command line (run for each sample):

```
STAR --runThreadN {threads} --genomeDir /path/to/reference/directory --readFilesIn sample_1.fastq.gz sample_2.fastq.gz --readFilesCommand zcat  --outFileNamePrefix sample 
```
This step will only take a few minutes per sample and primarily serves to identify potential splice junctions across all samples that will be incorporated in the second pass. 

####Second pass read alignment with STAR
The second pass of STAR is very similar to the first, but will include information about splice junctions identified in the first.  Again, like the first pass, this step will only take a few minutes per sample.

```
STAR --runThreadN {threads} --genomeDir /path/to/reference/directory --readFilesIn sample_1.fastq.gz sample_2.fastq.gz --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix sample --sjdbFileChrStartEnd list_of_all_sjbd_files
```

####Bam Processing
The next series of steps involves processing bam files and includes sorting, adding read groups, marking duplicates, and running GATK's "Split N Cigar" tool.  If we start with the file ```sample.bam```, the commands look something like

```
samtools sort -T path/to/temp/directory -O bam sample.bam > sample.sorted.bam

picard AddOrReplaceReadGroups INPUT=sample.sorted.bam OUTPUT=sample.sorted.rg.bam RGLB=sample RGPL=sample RGPU=illumina RGSM=sample RGID=sample VALIDATION_STRINGENCY=LENIENT

picard MarkDuplicates I=sample.sorted.rg.bam O=sample.sorted.rg.nodups.bam CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT M=sample.metrics.txt

java -Xmx12g -Djava.io.tmpdir=/path/to/temp/directory -jar /path/to/GATK.jar -T SplitNCigarReads -R /path/to/AnoCar2.0.fa -I sample.sorted.rg.nodups.bam -o sample.sorted.rg.nodups.splitncigar.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
```


####Sample variant calling outputting a GVCF
Then, we can call variants using the full processed bamfile (here, it's ```sample.sorted.rg.nodups.splitncigar.bam```).  We're going to run this separately for each sample to generate a gvcf file for each sample.  This will allow for joint genotyping of all samples later.

```
java -Xmx12g -Djava.io.tmpdir=/path/to/temp/directory -jar /path/to/GATK.jar -T HaplotypeCaller -R /path/to/AnoCar2.0.fa -I sample.sorted.rg.nodups.splitncigar.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o sample.g.vcf
```

####Joint genotyping all samples
Once all of the above steps have been run for all samples, we can jointly genotype all samples and then preliminarily filter the vcf using the following commands:

```
java -Xmx12g -Djava.io.tmpdir=/path/to/temp/directory -jar /path/to/GATK.jar -T GenotypeGVCFs -R /path/to/AnoCar2.0.fa --variant sample1.g.vcf --variant sample2.g.vcf --variant sample3.g.vcf -o allsamples.raw.vcf

bcftools view -m2 -M2 -v snps allsamples.raw.vcf | java -jar /path/to/SnpSift.jar filter '(QUAL >= 30) & (MQ >= 30) & (DP >= 40)' > allsamples.biallelicSNPS.QUAL30.MAPQ30.DP40.vcf
```
The first part of the filtering command selects biallelic SNPs only, while the second will only leave sites with QUAL >= 30, MAPQ >= 30, and a total depth >= 40.  Because our diversity ratio calculations will require sites to be callable in all samples within a single replicate, setting the total depth >= 40 will preliminarily remove any site that cannot possibly be included in any replicate.

##Calculating bootstrapped X/A diversity ratios
First, to calculate diversity at all sites, we need information about the number of sites that were callable overall.  This will help give us the number of monomorphic sites in addition to the variant sites in the VCF (monomorphic sites = number of callable sites - number of variant sites).

For each sample, we use GATK's CallableLoci to count the sites callable with a minimum depth of 10 and a minimum mean mapping quality of 30 (these are parameters we'll use on variant sites).  We can then use sed to grab the callable sites only from the output file.

```
java -jar /path/to/GATK.jar -T CallableLoci -R /path/to/AnoCar2.0.fa -I sample.sorted.rg.nodups.splitncigar.bam --minMappingQuality 30 --minDepth 10 --summary sample.summary -o sample.callablesites

sed -e '/CALLABLE/!d' sample.callablesites > sample.ONLYcallablesites.bed
```

Then, FOR EACH REPLICATE, we'll find the sites that are callable in all 4 samples with bedtools:

```bedtools intersect -a sample1.ONLYcallablesites.bed -b sample2.ONLYcallablesites.bed | bedtools intersect -a stdin -b sample3.ONLYcallablesites.bed | bedtools intersect -a stdin -b sample4.ONLYcallablesites.bed > replicate1.sharedcallable.bed```

Finally, now that we have our filtered vcf, and bed files containing callable sites for each replicate, we will use ```Anole_diversity_from_VCF.py``` to calculate diversity ratios between X-linked and autosomal sequences, and X-linked and microchromosome sequences.  The flags are fairly self explanatory and the commands used in this pipeline look like:
```
source activate diversity_script

python ../scripts/Anole_diversity_from_VCF.py --vcf allsamples.biallelicSNPS.QUAL30.MAPQ30.DP40.vcf --outfile x_autosome_ratios.txt --callable_regions replicate1_callable_loci.bed replicate2_callable_loci.bed replicate3_callable_loci.bed replicate4_callable_loci.bed replicate5_callable_loci.bed --bootstrap 1000 --male_list anole_male_list.txt --population_lists replicate1_pop_ids.txt replicate2_pop_ids.txt replicate3_pop_ids.txt replicate4_pop_ids.txt replicate5_pop_ids.txt --autosomal_scaffolds ../data/autosomes.txt --x_linked_scaffolds x_linked_scaffolds.txt --scaffold_sites_filter 250 --min_cov 10 --QD 2 --FS 30 --QUAL 30 --MAPQ 30

python ../scripts/Anole_diversity_from_VCF.py --vcf allsamples.biallelicSNPS.QUAL30.MAPQ30.DP40.vcf --outfile x_microchromosome_ratio.txt --callable_regions replicate1_callable_loci.bed replicate2_callable_loci.bed replicate3_callable_loci.bed replicate4_callable_loci.bed replicate5_callable_loci.bed --bootstrap 1000 --male_list anole_male_list.txt --population_lists replicate1_pop_ids.txt replicate2_pop_ids.txt replicate3_pop_ids.txt replicate4_pop_ids.txt replicate5_pop_ids.txt --autosomal_scaffolds ../data/microchromsomes_not_x_linked.txt --x_linked_scaffolds x_linked_scaffolds.txt --scaffold_sites_filter 250 --min_cov 10 --QD 2 --FS 30 --QUAL 30 --MAPQ 30

```
The ```source activate diversity_script``` will load the Python2 environment for the script, so that cyvcf is able to load (required for the script).

Note that in this pipeline, the callable loci bed files will be in the callable sites directory, and the population ids, scaffolds, and male list will be in the data directory.  ```--scaffold_sites_filter``` indicates the minimum number of callable sites required for a scaffold to be included, ```min_cov``` is the minimum depth per sample for it to be included, ```--QD``` is a measure of quality controling for quality inflation that can occur with very deep coverage (minimum), ```--FS``` is the Fisher strand score (maximum), ```--QUAL``` is the (minimum) site quality score, and ```--MAPQ``` is the (minimum) site mapping quality.

Also, note that this diversity script is designed specifically for this pipeline and will only work with VCFs output by GATK (it doesn't work with Freebayes, for instance, because of the different annotation by the two different programs).  I'll be hosting a more general script in another repository later (I'll provide a link here when it's up).



