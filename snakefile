configfile: "anoles.config.json"

SAMPLES = ["SRR1502164","SRR1502165","SRR1502166","SRR1502167","SRR1502168","SRR1502169","SRR1502170","SRR1502171","SRR1502172","SRR1502173","SRR1502174","SRR1502175","SRR1502176","SRR1502177","SRR1502178","SRR1502179","SRR1502180$

rule all:
    input:
		x_a = "results/x_autosome_ratios.txt",
		x_micro = "results/x_microchromosomes_ratios.txt"
                
rule first_pass_map:
    input:
        fq1="fastqs/{sample}_1.fastq.gz",
        fq2="fastqs/{sample}_2.fastq.gz"
    output:
        "star_first_pass/{sample}_SJ.out.tab"
    params:
        prefix="star_first_pass/{sample}_"
    threads: 8
    shell:
            "STAR --runThreadN {threads} --genomeDir reference --readFilesIn {input.fq1} {input.fq2} --readFilesCommand zcat --outFileNamePrefix {params.prefix}"


rule second_pass_map:
	input:
		sjs=expand("star_first_pass/{sample}_SJ.out.tab", sample=SAMPLES),
		fq1="fastqs/{sample}_1.fastq.gz",
		fq2="fastqs/{sample}_2.fastq.gz"
	output:
		"mapped_reads/{sample}_2ndpass_Aligned.out.bam"
	params:
		prefix="mapped_reads/{sample}_2ndpass_"
	threads: 8
	shell:
		"STAR --runThreadN {threads} --genomeDir reference --readFilesIn {input.fq1} {input.fq2} --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix {params.prefix} --sjdbFileChrStartEnd {input.sjs}"

rule samtools_sort:
	input:
		"mapped_reads/{sample}_2ndpass_Aligned.out.bam"
	output:
		"mapped_reads/{sample}.sorted.bam"
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

rule mark_dups:
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
		"vcf/all_anoles.biallelicSNPs.QUAL30.MAPQ30.totalDP40.vcf"
	params:
		Snpsift_path=config["Snpsift"]
	shell:
		"bcftools view -m2 -M2 -v snps {input} | java -jar {params.Snpsift_path} filter '(QUAL >= 30) & (MAPQ >= 30) & (DP >= 40)' > {output}"


rule diversity_analyses:
	input:
		call1 =
		call2 =
		call3 = 
		call4 =
		call5 =
		vcf = "vcf/all_anoles.biallelicSNPs.QUAL30.MAPQ30.totalDP40.vcf"
	output:
		x_a = "results/x_autosome_ratios.txt",
		x_micro = "results/x_microchromosomes_ratios.txt"
	params:
		male_list = "data/anole_male_list.txt",
		population_lists = "data/replicate1_pop_ids.txt data/replicate2_pop_ids.txt data/replicate3_pop_ids.txt data/replicate4_pop_ids.txt data/replicate5_pop_ids.txt",
		auto_list = "data/autosomes.txt",
		micro_list = "data/microchromsomes_not_x_linked.txt",
		x_list = "data/x_linked_scaffolds.txt"
	shell:
		"source activate diversity_script "
		"python scripts/Anole_diversity_from_VCF.py --vcf {input.vcf} --outfile {output.x_a} --callable_regions {input.call1} {input.call2} {input.call3} {input.call4} {input.call5} --bootstrap 1000 --male_list {params.male_list} --population_lists {params.population_lists} --autosomal_scaffolds {params.auto_list} --x_linked_scaffolds {params.x_list} --scaffold_sites_filter 250 --min_cov 10 --QD 2 --FS 30 --QUAL 30 --MAPQ 30"
		"python scripts/Anole_diversity_from_VCF.py --vcf {input.vcf} --outfile {output.x_micro} --callable_regions {input.call1} {input.call2} {input.call3} {input.call4} {input.call5} --bootstrap 1000 --male_list {params.male_list} --population_lists {params.population_lists} --autosomal_scaffolds {params.micro_list} --x_linked_scaffolds {params.x_list} --scaffold_sites_filter 250 --min_cov 10 --QD 2 --FS 30 --QUAL 30 --MAPQ 30"











# fastq1_list = []
# fastq2_list = []
# for i in config["samples"]:
# 	fastq1_list.append(i + "_1.fastq.gz")
# 	fastq1_list.append(i + "_2.fastq.gz")
# 
# 
# rule all:
# 	input:
# 		"anole_report.html"
# 		 
# rule fastq_dump:
#   	input:
#   		a=config["samples"]
#  		lambda wildcards: config["samples"][wildcards.sample]
#   	output:
#   		expand("fastqs/{fq_file}",fq_file=fastq1_list),
#   		expand("fastqs/{fq_file}",fq_file=fastq2_list)
#   	threads: 1
#   	shell:
#   		"fastq-dump --gzip --outdir fastqs --readids --split-files {input}"
#  	run:
#  		for i in config["samples"]:
#  			sra = i
#  			shell("fastq-dump --gzip --outdir fastqs --readids --split-files {sra}")
# 
# 
# rule download_genome:
# 	output:
# 		"reference/AnoCar2.0.fa",
# 		"reference/AnoCar2.0.85.gtf"
# 	shell:
# 		"wget -O reference/AnoCar2.0.fa.gz ftp://ftp.ensembl.org/pub/release-85/fasta/anolis_carolinensis/dna/Anolis_carolinensis.AnoCar2.0.dna.toplevel.fa.gz  "
# 		"&& gunzip reference/AnoCar2.0.fa.gz "
# 		"&& wget -O reference/AnoCar2.0.85.gtf.gz ftp://ftp.ensembl.org/pub/release-85/gtf/anolis_carolinensis/Anolis_carolinensis.AnoCar2.0.85.gtf.gz "
# 		"&& gunzip reference/AnoCar2.0.85.gtf.gz"
# 
# rule prep_genome:
# 	input:
# 		ref="reference/AnoCar2.0.fa",
# 		gtf="reference/AnoCar2.0.85.gtf"
# 	output:
# 		"reference/AnoCar2.0.fa.fai"
# 	threads: 8
# 	shell:
# 		"STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir reference --genomeFastaFiles {input.ref} --sjdbGTFfile {input.gtf} && "
# 		"samtools faidx {input.ref} "
# 		"picard CreateSequenceDictionary R={input.ref} O=reference/AnoCar2.0.dict"
# 
# rule first_pass_map:
# 	input:
# 		fq1="fastqs/{sample}_1.fastq.gz",
# 		fq2="fastqs/{sample}_2.fastq.gz"		
# 	output:
# 		"star_first_pass/{sample}_SJ.out.tab"
# 	threads: 8
# 	shell:
# 		"STAR --runThreadN {threads} --genomeDir reference --readFilesIn {input.fq1} {input.fq2} --outFileNamePrefix star_first_pass/{sample}_"
# 	
# rule second_pass_map:
# 	input:
# 		sjs=expand("star_first_pass/{sample_list}_SJ.out.tab", sample_list=config["samples"]),
# 		fq1="fastqs/{sample}_1.fastq.gz",
# 		fq2="fastqs/{sample}_2.fastq.gz"
# 	output:
# 		"mapped_reads/{sample}_2ndpass_Aligned.out.bam"
# 	threads: 8
# 	shell:
# 		"STAR --runThreadN {threads} --genomeDir reference --readFilesIn {input.fq1} {input.fq2} --outSAMtype BAM Unsorted --outFileNamePrefix mapped_reads/{sample}_2ndpass_ --sjdbFileChrStartEnd {input.sjs}"
# 
# rule samtools_sort:
# 	input:
# 		"mapped_reads/{sample}_2ndpass_Aligned.out.bam"
# 	output:
# 		temp("mapped_reads/{sample}.sorted.bam")
# 	params:
# 		temp_dir=config["temp_directory"]
# 	threads: 4
# 	shell:
# 		"samtools sort -T {params.temp_dir} -O bam {input} > {output}"
# 		
# rule add_rg:
# 	input:
# 		"mapped_reads/{sample}.sorted.bam"
# 	output:
# 		temp("mapped_reads/{sample}.sorted.rg.bam")
# 	params:
# 		RGLB="{sample}",
# 		RGPL="{sample}",
# 		RGPU="illumina",
# 		RGSM="{sample}",
# 		RGID="{sample}",
# 		temp_dir=config["temp_directory"],
# 	threads: 2
# 	shell:
# 		"picard AddOrReplaceReadGroups INPUT={input} OUTPUT={output} RGLB={params.RGLB} RGPL={params.RGPL} RGPU={params.RGPU} RGSM={params.RGSM} RGID={params.RGID} VALIDATION_STRINGENCY=LENIENT"
# 
# rule remove_dups:
# 	input:
# 		"mapped_reads/{sample}.sorted.rg.bam"
# 	output:
# 		"mapped_reads/{sample}.sorted.rg.nodups.bam"
# 	params:
# 		temp_dir=config["temp_directory"],
# 		metrics="stats/{sample}.picardmetrics"
# 	threads: 4
# 	shell:
# 		"picard MarkDuplicates I={input} O={output} CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT M={params.metrics}"
# 
# rule split_n_cigar:
# 	input:
# 		ref="reference/AnoCar2.0.fa",
# 		bam="mapped_reads/{sample}.sorted.rg.nodups.bam"
# 	output:
# 		"processed_bams/{sample}.sorted.rg.nodups.splitcigar.bam"
# 	threads: 4
# 	params:
# 		temp_dir=config["temp_directory"],
# 		gatk_path=config["GATK"]
# 	shell:
# 		"java -Xmx12g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T SplitNCigarReads -R {input.ref} -I {input.bam} -o {output} -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS"
# 
# rule call_gvcf:
# 	input:
# 		ref="reference/AnoCar2.0.fa",
# 		bam="processed_bams/{sample}.sorted.rg.nodups.splitcigar.bam"
# 	output:
# 		"calls/{sample}.g.vcf"
# 	params:
# 		temp_dir=config["temp_directory"],
# 		gatk_path=config["GATK"]
# 	threads: 4
# 	shell:
# 		"java -Xmx12g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T HaplotypeCaller -R {input.ref} -I {input.bam} -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o calls/{sample}.g.vcf"
# 
# rule generate_callable_sites:
# 	input:
# 		ref="reference/AnoCar2.0.fa",
# 		bam="processed_bams/{sample}.sorted.rg.nodups.splitcigar.bam"
# 	output:
# 		"callable_sites/{sample}.callablesites"
# 	params:
# 		temp_dir=config["temp_directory"],
# 		gatk_path=config["GATK"],
# 		summary="stats/{sample}.callable.summary"
# 	threads: 4
# 	shell:
# 		"java -Xmx12g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T CallableLoci -R {input.ref} -I {input.bam} --summary {params.summary} -o {output}"
# 
# rule extract_callable_sites:
# 	input:
# 		"callable_sites/{sample}.callablesites"
# 	output:
# 		"callable_sites/{sample}.ONLYcallablesites.bed"
# 	shell:
# 		"sed -e '/CALLABLE/!d' {input} > {output}"
# 		
# rule find_shared_callable_sites:
# 	input:
# 		beds=expand("callable_sites/{sample}.callablesites", sample=lambda wildcards: config["groups"][wildcards.group])
# 	output:
# 		"callable_sites/{group}.sharedcallable.bed"
# 	run:
# 		first = input.beds[0]
# 		second = input.beds[1]
# 		third = input.beds[2]
# 		fourth = input.beds[3]
# 		fifth = input.beds[4]
# 		shell("bedtools intersect -a {first} -b {second} | bedtools intersect -a stdin -b {third} | bedtools intersect -a stdin -b {fourth} > {output}")
# 			
# rule genotype_gvcfs:
# 	input:
# 		ref="reference/AnoCar2.0.fa",
# 		vcfs=expand("calls/{sample}.g.vcf", sample=lambda wildcards: config["sra"][wildcards.sra])
# 	output:
# 		"vcf/all_anoles.raw.vcf"
# 	params:
# 		temp_dir=config["temp_directory"],
# 		gatk_path=config["GATK"]
# 	threads: 4
# 	run:
# 		for i in input.vcfs:
# 			i = "--variant " + i
# 		shell("java -Xmx12g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T GenotypeGVCFs -R {input.ref} {input.vcfs} -o allsamples_anole.raw.vcf")
# 
# rule filter_vcf:
# 	input:
# 		"vcf/all_anoles.raw.vcf"
# 	output:
# 		"vcf/all_anoles.biallelicSNPs.QUAL30.sampleDP10.vcf"
# 	params:
# 		Snpsift_path=config["Snpsift"]
# 	shell:
# 		"bcftools view -m2 -M2 -v snps {input} | java -jar {params.Snpsift_path} filter '(QUAL >= 30) & (MAPQ >= 30) & (GEN[ALL].DP >= 10)' > {output}"
# 
#  rule diversity_analysis:
#  	input:
#  
# rule report:
#     input:
#         "vcf/all_anoles.raw.vcf"
#     output:
#         "anole_report.html"
#     run:
#         from snakemake.utils import report
#         with open(input[0]) as vcf:
#             n_calls = sum(1 for l in vcf if not l.startswith("#"))
# 
#         report("""
#         Anole transcriptome diversity analyses
#         ===================================
# 		Pipeline to replicate anole transcriptomic diversity analyses from Rupp et al.
# 		
#         Reads were mapped to the AnoCar2
#         reference genome using STAR (two pass) 
#         and variants called using GATK.
# 
#         This resulted in {n_calls} variants (see Table T1_).
#         """, output[0], T1=input[0])
	