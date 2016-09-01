configfile: "anoles.config.json"

SAMPLES = ["SRR1502164","SRR1502165","SRR1502166","SRR1502167","SRR1502168","SRR1502169","SRR1502170","SRR1502171","SRR1502172","SRR1502173","SRR1502174","SRR1502175","SRR1502176","SRR1502177","SRR1502178","SRR1502179","SRR1502180","SRR1502181","SRR1502182","SRR1502183"]
replicates = ["replicate1","replicate2","replicate3","replicate4","replicate5"]
group1 = ["SRR1502164","SRR1502169","SRR1502174","SRR1502179"]
group2 = ["SRR1502165","SRR1502170","SRR1502175","SRR1502180"]
group3 = ["SRR1502166","SRR1502171","SRR1502176","SRR1502181"]
group4 = ["SRR1502167","SRR1502172","SRR1502177","SRR1502182"]
group5 = ["SRR1502168","SRR1502173","SRR1502178","SRR1502183"]


rule all:
	input:
		"results/x_autosome_ratios.txt",
		"results/x_microchromosomes_ratios.txt"
		
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
		temp_dir=config["temp_directory"]
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
		"java -Xmx12g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T HaplotypeCaller -R {input.ref} -I {input.bam} -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o {output}"

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
		"java -Xmx12g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T CallableLoci -R {input.ref} -I {input.bam} --minMappingQuality 30 --minDepth 10 --summary {params.summary} -o {output}"

rule extract_callable_sites:
	input:
		"callable_sites/{sample}.callablesites"
	output:
		"callable_sites/{sample}.ONLYcallablesites.bed"
	shell:
		"sed -e '/CALLABLE/!d' {input} > {output}"

rule find_shared_callable_sites1:
	input:
		beds=expand("callable_sites/{id}.ONLYcallablesites.bed", id=group1)
	output:
		"callable_sites/replicate1.sharedcallable.bed"
	run:
		first = input.beds[0]
		second = input.beds[1]
		third = input.beds[2]
		fourth = input.beds[3]
		shell("bedtools intersect -a {first} -b {second} | bedtools intersect -a stdin -b {third} | bedtools intersect -a stdin -b {fourth} > {output}")
		
rule find_shared_callable_sites2:
	input:
		beds=expand("callable_sites/{id}.ONLYcallablesites.bed", id=group2)
	output:
		"callable_sites/replicate2.sharedcallable.bed"
	run:
		first = input.beds[0]
		second = input.beds[1]
		third = input.beds[2]
		fourth = input.beds[3]
		shell("bedtools intersect -a {first} -b {second} | bedtools intersect -a stdin -b {third} | bedtools intersect -a stdin -b {fourth} > {output}")

rule find_shared_callable_sites3:
	input:
		beds=expand("callable_sites/{id}.ONLYcallablesites.bed", id=group3)
	output:
		"callable_sites/replicate3.sharedcallable.bed"
	run:
		first = input.beds[0]
		second = input.beds[1]
		third = input.beds[2]
		fourth = input.beds[3]
		shell("bedtools intersect -a {first} -b {second} | bedtools intersect -a stdin -b {third} | bedtools intersect -a stdin -b {fourth} > {output}")

rule find_shared_callable_sites4:
	input:
		beds=expand("callable_sites/{id}.ONLYcallablesites.bed", id=group4)
	output:
		"callable_sites/replicate4.sharedcallable.bed"
	run:
		first = input.beds[0]
		second = input.beds[1]
		third = input.beds[2]
		fourth = input.beds[3]
		shell("bedtools intersect -a {first} -b {second} | bedtools intersect -a stdin -b {third} | bedtools intersect -a stdin -b {fourth} > {output}")
		
rule find_shared_callable_sites5:
	input:
		beds=expand("callable_sites/{id}.ONLYcallablesites.bed", id=group5)
	output:
		"callable_sites/replicate5.sharedcallable.bed"
	run:
		first = input.beds[0]
		second = input.beds[1]
		third = input.beds[2]
		fourth = input.beds[3]
		shell("bedtools intersect -a {first} -b {second} | bedtools intersect -a stdin -b {third} | bedtools intersect -a stdin -b {fourth} > {output}")
			
rule genotype_gvcfs:
	input:
		ref="reference/AnoCar2.0.fa",
		vcfs=expand("calls/{sample}.g.vcf", sample=SAMPLES)
	output:
		"vcf/all_anoles.raw.vcf"
	params:
		temp_dir=config["temp_directory"],
		gatk_path=config["GATK"]
	threads: 4
	run:
		variant_files = []
		for i in input.vcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell("java -Xmx12g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk_path} -T GenotypeGVCFs -R {input.ref} {variant_files} -o {output}")

rule filter_vcf:
	input:
		"vcf/all_anoles.raw.vcf"
	output:
		"vcf/all_anoles.biallelicSNPs.QUAL30.MAPQ30.totalDP40.vcf"
	params:
		Snpsift_path=config["Snpsift"]
	shell:
		"bcftools view -m2 -M2 -v snps {input} | java -jar {params.Snpsift_path} filter '(QUAL >= 30) & (MQ >= 30) & (DP >= 40)' > {output}"
		
rule diversity_analyses:
	input:
		callable = expand("callable_sites/{id}.sharedcallable.bed", id=replicates),
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
		"source activate diversity_script && python Scripts/Anole_diversity_from_VCF.py --vcf {input.vcf} --outfile {output.x_a} --callable_regions {input.callable} --bootstrap 1000 --male_list {params.male_list} --population_lists {params.population_lists} --autosomal_scaffolds {params.auto_list} --x_linked_scaffolds {params.x_list} --scaffold_sites_filter 250 --min_cov 10 --QD 2 --FS 30 --QUAL 30 --MAPQ 30 && python Scripts/Anole_diversity_from_VCF.py --vcf {input.vcf} --outfile {output.x_micro} --callable_regions {input.callable} --bootstrap 1000 --male_list {params.male_list} --population_lists {params.population_lists} --autosomal_scaffolds {params.micro_list} --x_linked_scaffolds {params.x_list} --scaffold_sites_filter 250 --min_cov 10 --QD 2 --FS 30 --QUAL 30 --MAPQ 30"
