configfile: "anoles.config.json"

rule all:
	input:
		"anole_report.html"
		
rule fastq_dump:
	input:
		lambda wildcards: config["sra"][wildcards.sra]
	output:
		"fastqs/{sample}_1.fastq.gz",
		"fastqs/{sample}_2.fastq.gz"
	threads: 1
	shell:
		"fastq-dump --gzip --outdir fastqs --readids --split-files {input}"

rule download_genome:

rule prep_genome:

rule first_pass_map:
	input:
	output:
	threads:
	shell:
	
rule second_pass_map:

rule sort_bam:

rule add_rg:

rule mark_dups:

rule split_n_cigar:

rule call_gvcf:

rule generate_callable_sites:
	input:
		"processed_bams/{sample}.sorted.rg.nodups.bam"
	output:
		"callable_sites/{sample}.callablesites"
	params:
		temp_dir=config["temp_directory"]
		GATKPath=config["gatk_path"]
		ref=config["genome_path"]
		summary="stats/{sample}.callable.summary"
	threads: 4
	shell:
		"java -Djava.io.tmpdir={params.temp_dir} -jar -Xmx12g {GATKPath} -T CallableLoci -R {params.ref} -I {input} --summary {summary} -o {output}"

rule extract_callable_sites:
	input:
		"callable_sites/{sample}.callablesites"
	output:
		"callable_sites/{sample}.ONLYcallablesites.bed"
	shell:
		"sed -e '/CALLABLE/!d' {input} > {output}"
		
rule find_shared_callable_sites:

rule genotype_gvcfs:
	input:
	output:
		"vcf/all_anoles.raw.vcf"

rule diversity_analysis:
	input:

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
	