import sys

# configfile:'config.yaml'
SAMPLES = ['HBCRACKHiC-K562-DN-TD-R1_GCCAAT_L008_test']
data_dir = '/home/tb37w/project/Research/digest/dpnII'

rule all:
	input:
		#dynamic(expand('split_reads/{sample}/split_{{n}}', sample=SAMPLES))
		#expand('test_{sample}.txt', sample=SAMPLES)
		#dynamic(expand('distances/{sample}/3_distance_{{n}}', sample=SAMPLES))
		#expand('{sample}_3_distance.sam', sample=SAMPLES),
		#expand('{sample}_abs_frag_distances.txt', sample=SAMPLES)
		expand('plots/{sample}_distance_hist.pdf', sample=SAMPLES),
		#expand('{sample}_coverage_500kb.bed', sample=SAMPLES)
		expand('{sample}_copy_correct_coverage_500kb.bed', sample=SAMPLES),
		expand('{sample}_copy_correct_coverage_40kb.bed', sample=SAMPLES)


# Map reads
rule bowtie_map:
	input:
		R1=data_dir + '/data/fastq/21OCT15_PE100_C7MDWAC/HBCRACKHiC-K562-DN-TD-R1/'
			'{sample}_R1.fastq',
		R2=data_dir + '/data/fastq/21OCT15_PE100_C7MDWAC/HBCRACKHiC-K562-DN-TD-R1/'
			'{sample}_R2.fastq'
	output:
		'mapped_reads/{sample}.sam'
	log:
		'logs/bowtie_map/{sample}_bowtie_map.log'
	threads: 10
	shell:
		'bowtie -S -p {threads} -t -m 1 --chunkmbs 1000 -X 10000 '
		+ data_dir + '/data/indexes/hg19_build -1 {input.R1} '
		'-2 {input.R2} {output} 2> {log}'

rule remove_multimappers:
	input:
		'mapped_reads/{sample}.sam'
	output:
		#b='remove_multimappers/{sample}_remove_multi.bam',
		s='remove_multimappers/{sample}_remove_multi.sam'		
	log:
		'logs/remove_multimappers/{sample}_remove_multimappers.log'
	shell:
		'''
		samtools view -bS {input} | samtools view -F 4 -h - > {output.s} 
		'''
rule get_header:
	input:
		s='remove_multimappers/{sample}_remove_multi.sam'
	output:
		h='remove_multimappers/{sample}_remove_multi.header'
	shell:
		'''
		head -27 {input.s} > {output.h}
		'''
		
rule split:
	input:
		s='remove_multimappers/{sample}_remove_multi.sam'
	output:
		dynamic('split_reads/{sample}/split_{n}')
	shell:
		'''
		tail -n+28 {input.s} | split -a 3 -d -l 100 - \
		split_reads/{wildcards.sample}/split_
		'''
rule filter_near_sites:
	input: 
		#a=dynamic('split_reads/{sample}/split_{n}')
		s ='split_reads/{sample}/split_{n}'
	output:
		'distances/{sample}/3_distance_{n}',
		'distances/{sample}/abs_frag_distances_{n}'
	shell:
		'''
		scripts/dpnII_abs_frag_distance.py -s {input.s} \
		-w {wildcards.sample}
		'''

rule merge_sams:
	input:
		d=dynamic('distances/{sample}/3_distance_{n}'),
		h='remove_multimappers/{sample}_remove_multi.header'
	output:
		t='{sample}_3_distance.sam'
	shell:
		'''
		cat {input.h} {input.d} > {output.t}
		'''
rule merge_distances:
	input:
		d=dynamic('distances/{sample}/abs_frag_distances_{n}')
	output:
		t='{sample}_abs_frag_distances.txt'
	shell:
		'''
		cat {input.d} > {output.t}
		'''
rule plot_distance_distribution:
	input:
		d = '{sample}_abs_frag_distances.txt'
	output:
		p = 'plots/{sample}_distance_hist.pdf'
	shell:
		'''
		scripts/dpnII_abs_frag_distance_hist.R \
		-f {input.d} -w {wildcards.sample}
		''' 

rule bin_500kb:
	input:
		r='{sample}_3_distance.sam'
	output:
		c='{sample}_coverage_500kb.bed'
	shell:
		'samtools view -bS {input.r} | bedtools coverage '
		'-abam - -b ' + data_dir + '/data/binning/hg19_500kb.bed | '
		'grep -v chrM |  grep -v chrY | sort -k1,1V -k2,2n '
		'> {output.c}'

rule bin_40kb:
	input:
		r='{sample}_3_distance.sam'
	output:
		c='{sample}_coverage_40kb.bed'
	shell:
		'samtools view -bS {input.r} | bedtools coverage '
		'-abam - -b ' + data_dir + '/data/binning/hg19_40kb.bed | '
		'grep -v chrM |  grep -v chrY | sort -k1,1V -k2,2n '
		'> {output.c}'

rule copy_number_correct_500kb:
	input:
		c='{sample}_coverage_500kb.bed'
	output:
		o='{sample}_copy_correct_coverage_500kb.bed'
	shell:
		'scripts/copy_correct.py -i {input.c} '
		'-c ' + data_dir + '/data/copy_number/K562_copynumber_500kb.bed '
		'> {output.o}'

rule copy_number_correct_40kb:
	input:
		c='{sample}_coverage_40kb.bed'
	output:
		o='{sample}_copy_correct_coverage_40kb.bed'
	shell:
		'scripts/copy_correct.py -i {input.c} '
		'-c ' + data_dir + '/data/copy_number/K562_copynumber_40kb.bed '
		'> {output.o}'








		

		
