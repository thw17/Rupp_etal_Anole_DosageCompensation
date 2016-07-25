# configfile: "anoles.config.json"
# 
# rule all:
# 	input:
# 		"anole_report.html"
		
# rule fastq_dump:
# 	input:
# 		lambda wildcards: config["sra"][wildcards.sra]
# 	output:
# 		"fastqs/{sample}_1.fastq.gz",
# 		"fastqs/{sample}_2.fastq.gz"
# 	threads: 1
# 	shell:
# 		"fastq-dump --gzip --outdir fastqs --readids --split-files {input}"

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
#	params:
#		picardPath=config["picard_path"]
	shell:
		"star --runMode genomeGenerate --genomeDir reference --genomeFastaFiles {input.ref} --sjdbGTFfile {input.gtf} && "
		"samtools faidx {input.ref} "
#		"java -Djava.io.tmpdir={params.temp_dir} -jar -Xmx2g {params.picardPath} CreateSequenceDictionary R={input.ref} O=reference/AnoCar2.0.dict"

# rule first_pass_map:
# 	input:
# 	output:
# 	threads:
# 	shell:
# 	
# rule second_pass_map:
# 
# rule samtools_sort:
# 	input:
# 		"mapped_reads/{sample}.bam"
# 	output:
# 		temp("mapped_reads/{sample}.sorted.bam")
# 	threads: 8
# 	params:
# 		temp_dir=config["temp_directory"]
# 	shell:
# 		"module add java/latest && samtools sort -@ 6 -m 4G -T {params.temp_dir} -O bam {input} > {output}"
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
# 		picardPath=config["picard_path"]
# 	threads: 1
# 	shell:
# 		"module add java/latest && java -Djava.io.tmpdir={params.temp_dir} -jar -Xmx2g {params.picardPath} AddOrReplaceReadGroups INPUT={input} OUTPUT={output} RGLB={params.RGLB} RGPL={params.RGPL} RGPU={params.RGPU} RGSM={params.RGSM} VALIDATION_STRINGENCY=LENIENT"
# 
# rule remove_dups:
# 	input:
# 		"mapped_reads/{sample}.sorted.rg.bam"
# 	output:
# 		protected("processed_bams/{sample}.sorted.rg.nodups.bam")
# 	params:
# 		temp_dir=config["temp_directory"],
# 		picardPath=config["picard_path"],
# 		metrics="stats/{sample}.picardmetrics"
# 	threads: 4
# 	shell:
# 		"module add java/latest && java -Djava.io.tmpdir={params.temp_dir} -jar -Xmx12g {params.picardPath} MarkDuplicates I={input} O={output} CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT M={params.metrics}"
# 
# rule split_n_cigar:
# 
# rule call_gvcf:
# 
# rule generate_callable_sites:
# 	input:
# 		"processed_bams/{sample}.sorted.rg.nodups.bam"
# 	output:
# 		"callable_sites/{sample}.callablesites"
# 	params:
# 		temp_dir=config["temp_directory"]
# 		GATKPath=config["gatk_path"]
# 		ref=config["genome_path"]
# 		summary="stats/{sample}.callable.summary"
# 	threads: 4
# 	shell:
# 		"java -Djava.io.tmpdir={params.temp_dir} -jar -Xmx12g {GATKPath} -T CallableLoci -R {params.ref} -I {input} --summary {summary} -o {output}"
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
# 
# rule genotype_gvcfs:
# 	input:
# 	output:
# 		"vcf/all_anoles.raw.vcf"
# 
# rule diversity_analysis:
# 	input:
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
# 	