import sys
import math
# configfile:'config.yaml'

data_dir = '/home/tb37w/project/Research/digest/dpnII'

rule all:
	input:
		'distances/HBCRACKHiC-K562-DN-TD-R1_GCCAAT_L008_test_3_distance_001',
		'distances/HBCRACKHiC-K562-DN-TD-R1_GCCAAT_L008_test_dpnII_abs_frag_distances_001'

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
		b='mapped_reads/{sample}_remove_multi.bam',
		s='mapped_reads/{sample}_remove_multi.sam'
	shell:
		'samtools view -bS {input} | samtools view -F 4 -h - | '
		'tee {output.s} | samtools view -bS - > {output.b}'

# rule count_splitfiles:
# 	input:
# 		'mapped_reads/{sample}_remove_multi.sam'
# 	output:
# 		'mapped_reads/{sample}_numfiles.txt'
# 	message:
# 		'FUCK'
# 	run:
# 		with open(input[0]) as f:
# 			for i, l in enumerate(f):
# 				pass
# 		numlines = i + 1
# 		if numlines % 1000 == 0:
# 			numfiles = numlines / 1000
# 		else:
# 			numfiles = math.floor(numlines / 1000) + 1
# 		with open(output[0], 'w') as o:
# 			o.write(str(numfiles) + '\n')
# 			o.write(str(numlines) + '\n')

rule splitfiles:
	input:
		s='mapped_reads/{sample}_remove_multi.sam'
	output:
		# Don't know n ahead of time how the fuck do i pass this in?
		expand('mapped_reads/{{sample}}_{n}', n=['{:03d}'.format(x) for x in range(15)])
		#'mapped_reads/{sample}_000'

	shell:
		'tail -n+28 {input.s} | split -a 3 -d -l 100 - '
		'mapped_reads/{wildcards.sample}_'

rule filter_near_sites:
	input:
		expand('mapped_reads/{{sample}}_{n}', n= ['{:03d}'.format(x) for x in range(1,15)])
		#'mapped_reads/{sample}_000'
	output:
		expand('distances/{{sample}}_dpnII_abs_frag_distances_{n}',n= ['{:03d}'.format(x) for x in range(1,15)]),
		expand('distances/{{sample}}_3_distance_{n}',n= ['{:03d}'.format(x) for x in range(1,15)])
	threads:14
	run:
		for i in input:
			shell('scripts/dpnII_abs_frag_distance.py -s {i} -w {wildcards.sample}')

		

		
