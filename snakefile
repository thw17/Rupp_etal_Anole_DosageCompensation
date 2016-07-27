configfile: "anoles.config.json"

fastq1_list = []
fastq2_list = []
for i in config["samples"]:
	fastq1_list.append(i + "_1.fastq.gz")
	fastq1_list.append(i + "_2.fastq.gz")


rule all:
	input:
		"anole_report.html"
		 
# rule fastq_dump:
#  	input:
#  		a=config["samples"]
#  		lambda wildcards: config["samples"][wildcards.sample]
#  	output:
#  		expand("fastqs/{fq_file}",fq_file=fastq1_list),
#  		expand("fastqs/{fq_file}",fq_file=fastq2_list)
#  	threads: 1
#  	shell:
#  		"fastq-dump --gzip --outdir fastqs --readids --split-files {input}"
# 	run:
# 		for i in config["samples"]:
# 			sra = i
# 			shell("fastq-dump --gzip --outdir fastqs --readids --split-files {sra}")


rule download_genome:
	output:
		"reference/AnoCar2.0.fa",
		"reference/AnoCar2.0.85.gtf"
	shell:
		"wget -O reference/AnoCar2.0.fa.gz ftp://ftp.ensembl.org/pub/release-85/fasta/anolis_carolinensis/dna/Anolis_carolinensis.AnoCar2.0.dna.toplevel.fa.gz  "
		"&& gunzip reference/AnoCar2.0.fa.gz "
		"&& wget -O reference/AnoCar2.0.85.gtf.gz ftp://ftp.ensembl.org/pub/release-85/gtf/anolis_carolinensis/Anolis_carolinensis.AnoCar2.0.85.gtf.gz "
		"&& gunzip reference/AnoCar2.0.85.gtf.gz"

rule prep_genome:
	input:
		ref="reference/AnoCar2.0.fa",
		gtf="reference/AnoCar2.0.85.gtf"
	output:
		"reference/AnoCar2.0.fa.fai"
	threads: 4
	shell:
		"star --runMode genomeGenerate --numThreadN {threads} --genomeDir reference --genomeFastaFiles {input.ref} --sjdbGTFfile {input.gtf} && "
		"samtools faidx {input.ref} "
		"picard CreateSequenceDictionary R={input.ref} O=reference/AnoCar2.0.dict"

rule first_pass_map:
	input:
		fq1="fastqs/{sample}_1.fastq.gz",
		fq2="fastqs/{sample}_2.fastq.gz"		
	output:
		"star_first_pass/{sample}_SJ.out.tab"
	threads: 4
	shell:
		"star ---numThreadN {threads} --genomeDir reference --readFilesIn {input.fq1} {input.fq2} --outFileNamePrefix star_first_pass/{sample}_"
	
rule second_pass_map:
	input:
		sjs=expand("star_first_pass/{sample_list}_SJ.out.tab", sample_list=config["samples"]),
		fq1="fastqs/{sample}_1.fastq.gz",
		fq2="fastqs/{sample}_2.fastq.gz"
	output:
		"mapped_reads/{sample}_2ndpass_Aligned.out.bam"
	threads: 4
	shell:
		"star ---numThreadN {threads} --genomeDir reference --readFilesIn {input.fq1} {input.fq2} --outSAMtype BAM Unsorted --outFileNamePrefix mapped_reads/{sample}_2ndpass_ --sjdbFileChrStartEnd {input.sjs}"

rule samtools_sort:
	input:
		"mapped_reads/{sample}_2ndpass_Aligned.out.bam"
	output:
		temp("mapped_reads/{sample}.sorted.bam")
	params:
		temp_dir=config["temp_directory"]
	threads: 4
	shell:
		"samtools sort -T {params.temp_dir} -O bam {input} > {output}"
		
rule add_rg:
	input:
		"mapped_reads/{sample}.sorted.bam"
	output:
		temp("mapped_reads/{sample}.sorted.rg.bam")
	params:
		RGLB="{sample}",
		RGPL="{sample}",
		RGPU="illumina",
		RGSM="{sample}",
		RGID="{sample}",
		temp_dir=config["temp_directory"],
	threads: 2
	shell:
		"picard AddOrReplaceReadGroups INPUT={input} OUTPUT={output} RGLB={params.RGLB} RGPL={params.RGPL} RGPU={params.RGPU} RGSM={params.RGSM} RGID={params.RGID} VALIDATION_STRINGENCY=LENIENT"

rule remove_dups:
	input:
		"mapped_reads/{sample}.sorted.rg.bam"
	output:
		"mapped_reads/{sample}.sorted.rg.nodups.bam"
	params:
		temp_dir=config["temp_directory"],
		metrics="stats/{sample}.picardmetrics"
	threads: 4
	shell:
		"picard MarkDuplicates I={input} O={output} CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT M={params.metrics}"

rule split_n_cigar:
	input:
		ref="reference/AnoCar2.0.fa",
		bam="mapped_reads/{sample}.sorted.rg.nodups.bam"
	output:
		"processed_bams/{sample}.sorted.rg.nodups.splitcigar.bam"
	threads: 4
	params:
		temp_dir=config["temp_directory"],
		gatk_path=config["GATK"]
	shell:
		"java -Xmx12g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T SplitNCigarReads -R {input.ref} -I {input.bam} -o {output} -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS"

rule call_gvcf:
	input:
		ref="reference/AnoCar2.0.fa",
		bam="processed_bams/{sample}.sorted.rg.nodups.splitcigar.bam"
	output:
		"calls/{sample}.g.vcf"
	params:
		temp_dir=config["temp_directory"],
		gatk_path=config["GATK"]
	threads: 4
	shell:
		"java -Xmx12g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T HaplotypeCaller -R {input.ref} -I {input.bam} -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o calls/{sample}.g.vcf"

rule generate_callable_sites:
	input:
		ref="reference/AnoCar2.0.fa",
		bam="processed_bams/{sample}.sorted.rg.nodups.splitcigar.bam"
	output:
		"callable_sites/{sample}.callablesites"
	params:
		temp_dir=config["temp_directory"],
		gatk_path=config["GATK"],
		summary="stats/{sample}.callable.summary"
	threads: 4
	shell:
		"java -Xmx12g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T CallableLoci -R {input.ref} -I {input.bam} --summary {params.summary} -o {output}"

rule extract_callable_sites:
	input:
		"callable_sites/{sample}.callablesites"
	output:
		"callable_sites/{sample}.ONLYcallablesites.bed"
	shell:
		"sed -e '/CALLABLE/!d' {input} > {output}"
		
rule find_shared_callable_sites:
	input:
		beds=expand("callable_sites/{sample}.callablesites", sample=lambda wildcards: config["groups"][wildcards.group])
	output:
		"callable_sites/{group}.sharedcallable.bed"
	run:
		first = input.beds[0]
		second = input.beds[1]
		third = input.beds[2]
		fourth = input.beds[3]
		fifth = input.beds[4]
		shell("bedtools intersect -a {first} -b {second} | bedtools intersect -a stdin -b {third} | bedtools intersect -a stdin -b {fourth} > {output}")
			
rule genotype_gvcfs:
	input:
		ref="reference/AnoCar2.0.fa",
		vcfs=expand("calls/{sample}.g.vcf", sample=lambda wildcards: config["sra"][wildcards.sra])
	output:
		"vcf/all_anoles.raw.vcf"
	params:
		temp_dir=config["temp_directory"],
		gatk_path=config["GATK"]
	threads: 4
	run:
		for i in input.vcfs:
			i = "--variant " + i
		shell("java -Xmx12g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T GenotypeGVCFs -R {input.ref} {input.vcfs} -o allsamples_anole.raw.vcf")

rule filter_vcf:
	input:
		"vcf/all_anoles.raw.vcf"
	output:
		"vcf/all_anoles.biallelicSNPs.QUAL30.sampleDP10.vcf"
	params:
		Snpsift_path=config["SnpSift"]
	shell:
		"bcftools view -m2 -M2 -v snps {input} | java -jar {params.Snpsift_path} filter '(QUAL >= 30) & (GEN[ALL].DP >= 10)' > {output}"

# rule diversity_analysis:
# 	input:
# 
rule report:
    input:
        "vcf/all_anoles.raw.vcf"
    output:
        "anole_report.html"
    run:
        from snakemake.utils import report
        with open(input[0]) as vcf:
            n_calls = sum(1 for l in vcf if not l.startswith("#"))

        report("""
        Anole transcriptome diversity analyses
        ===================================
		Pipeline to replicate anole transcriptomic diversity analyses from Rupp et al.
		
        Reads were mapped to the AnoCar2
        reference genome using STAR (two pass) 
        and variants called using GATK.

        This resulted in {n_calls} variants (see Table T1_).
        """, output[0], T1=input[0])
	