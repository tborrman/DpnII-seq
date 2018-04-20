import sys
import math
# configfile:'config.yaml'
SAMPLES = ['HBCRACKHiC-K562-DN-TD-R1_GCCAAT_L008_test']
data_dir = '/home/tb37w/project/Research/digest/dpnII'

rule all:
	input:
		dynamic(expand('split_reads/{sample}/split_{{n}}', sample=SAMPLES))

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
		b='remove_multimappers/{sample}_remove_multi.bam',
		s='remove_multimappers/{sample}_remove_multi.sam'
	
	shell:
		'''
		samtools view -bS {input} | samtools view -F 4 -h - | \
		tee {output.s} | samtools view -bS - > {output.b}
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


# input:
# 	splitty/{sample}/split_{n}
# rule filter_near_sites:
# 	input:
# 		expand('mapped_reads/{{sample}}_{n}', n= ['{:03d}'.format(x) for x in range(1,15)])
# 		#'mapped_reads/{sample}_000'
# 	output:
# 		expand('distances/{{sample}}_dpnII_abs_frag_distances_{n}',n= ['{:03d}'.format(x) for x in range(1,15)]),
# 		expand('distances/{{sample}}_3_distance_{n}',n= ['{:03d}'.format(x) for x in range(1,15)])
# 	threads:14
# 	run:
# 		for i in input:
# 			shell('scripts/dpnII_abs_frag_distance.py -s {i} -w {wildcards.sample}')

		

		
