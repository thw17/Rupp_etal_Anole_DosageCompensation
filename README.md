# Anole X/A diversity analyses from Rupp et al.
##CURRENTLY UNDER CONSTRUCTION (7/17/2016)

This repository contains scripts and information related to the genetic diversity analyses in Rupp et al (_In review_). Evolution of dosage compensation in _Anolis carolinensis_, a reptile with XX/XY chromosomal sex determination.  Scripts related to other parts of the study (e.g. identifying X-linked scaffolds, differential expression analyses, and Ka/Ks calculations) can be found in [another repository](https://github.com/WilsonSayresLab/Anole_expression).

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
fastq-dump --split-files --gzip $i
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
We used STAR's 2-pass method to map reads.  Assuming the genome ([AnoCar2](http://hgdownload.cse.ucsc.edu/goldenPath/anoCar2/bigZips/)) has been downloaded and properly indexed, the first pass with STAR is relatively straightforward.  We used the following template command line (run for each sample):

```
STAR=/path/to/STAR
refdir=/path/to/reference/genome
fastqdir=/path/to/fastq/directory
sampleSRA=SRAaccessionForSample
numthread=NumberofThreads
outdir=/path/to/outdir

$STAR --genomeDir $refdir --readFilesIn "$fastqdir""$sampleSRA"_1.fastq "fastqdir""$sampleSRA"_2.fastq --runThreadN $numthread --outFileNamePrefix "$outdir""$sampleSRA"_
```
This step will only take a few minutes per sample and primarily serves to identify potential splice junctions across all samples that will be incorporated in the second pass. 

####Second pass read alignment with STAR
The second p

####Bam Processing

####Sample variant calling outputting a GVCF

####Joint genotyping all samples


##Calculating bootstrapped X/A diversity ratios



