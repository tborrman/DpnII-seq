import sys

configfile: 'config.yaml'
#SAMPLES = ['SAMPLE_ID']
#SAMPLES = ['HBCRACKHiC-K562-DN-TD-R1_GCCAAT_L008']
#data_dir = '/home/tb37w/project/Research/digest/dpnII'

print(config['restrict_enzy'])

rule all:
	input:
		expand('plots/{sample}_distance_hist.pdf', sample=config['samples']),
		expand('copy_correct_coverage/500kb/{sample}_copy_correct_coverage_500kb.bed', sample=config['samples']),
		expand('copy_correct_coverage/40kb/{sample}_copy_correct_coverage_40kb.bed', sample=config['samples'])


# Map reads
rule bowtie_map:
	input:
		R1=config['data_dir'] + '/fastq/{sample}_R1.fastq',
		R2=config['data_dir'] + '/fastq/{sample}_R2.fastq'

	output:
		'mapped_reads/{sample}.sam'
	log:
		'logs/bowtie_map/{sample}_bowtie_map.log'
	threads: 10
	shell:
		'bowtie -S -p {threads} -t -m 1 --chunkmbs 1000 -X 10000 '
		+ config['data_dir'] + '/indexes/hg19_build -1 {input.R1} '
		'-2 {input.R2} {output} 2> {log}'

rule remove_multimappers:
	input:
		'mapped_reads/{sample}.sam'
	output:
		s='remove_multimappers/{sample}_remove_multi.sam'		
	shell:
		'''
		samtools view -bS {input} | samtools view -F 4 -h - > {output.s} 
		'''
rule get_header:
	input:
		s='remove_multimappers/{sample}_remove_multi.sam'
	output:
		h=temp('remove_multimappers/{sample}_remove_multi.header')
	shell:
		'''
		head -27 {input.s} > {output.h}
		'''
		
rule split:
	input:
		s='remove_multimappers/{sample}_remove_multi.sam'
	output:
		dynamic('split_reads/{sample}/split_{n}')
		#dynamic('split_reads/{sample}/split_{n}')
	shell:
		'''
		tail -n+28 {input.s} | split -a 4 -d -l 1000000 - \
		split_reads/{wildcards.sample}/split_
		'''
	
if config['restrict_enzy'] == 'DpnII':
	rule filter_near_sites:
		input: 
			s ='split_reads/{sample}/split_{n}'
		output:
			'filtered_sams/{sample}/3_distance_{n}',
			'distances/{sample}/abs_frag_distances_{n}'
		shell:
			'scripts/dpnII_abs_frag_distance.py -s {input.s} '
			'-w {wildcards.sample} -r ' + config['read_length'] +
			' -d ' + config['data_dir']

elif config['restrict_enzy'] == 'HindIII':
	rule filter_near_sites:
		input: 
			s ='split_reads/{sample}/split_{n}'
		output:
			'filtered_sams/{sample}/3_distance_{n}',
			'distances/{sample}/abs_frag_distances_{n}'
		shell:
			'scripts/hindIII_abs_frag_distance.py -s {input.s} '
			'-w {wildcards.sample} -r ' + config['read_length']
else:
	print('No restriction enzyme given in config.yaml')
	sys.exit()

rule merge_sams:
	input:
		d=dynamic('filtered_sams/{sample}/3_distance_{n}'),
		h='remove_multimappers/{sample}_remove_multi.header'
	output:
		t='filtered_sams/{sample}/{sample}_3_distance.sam'
	shell:
		'''
		cat {input.h} {input.d} > {output.t}
		'''

rule merge_distances:
	input:
		d=dynamic('distances/{sample}/abs_frag_distances_{n}')
	output:
		t='distances/{sample}/{sample}_abs_frag_distances.txt'
	shell:
		'''
		cat {input.d} > {output.t}
		'''

if config['restrict_enzy'] == 'DpnII':
	rule plot_distance_distribution:
		input:
			d = 'distances/{sample}/{sample}_abs_frag_distances.txt'
		output:
			p = 'plots/{sample}_distance_hist.pdf'
		shell:
			'''
			scripts/dpnII_abs_frag_distance_hist.R \
			-f {input.d} -w {wildcards.sample}
			''' 

elif config['restrict_enzy'] == 'HindIII':
	rule plot_distance_distribution:
		input:
			d = 'distances/{sample}/{sample}_abs_frag_distances.txt'
		output:
			p = 'plots/{sample}_distance_hist.pdf'
		shell:
			'''
			scripts/hindIII_abs_frag_distance_hist.R \
			-f {input.d} -w {wildcards.sample}
			''' 
else:
	print('No restriction enzyme given in config.yaml')
	sys.exit()

rule bin_500kb:
	input:
		r='filtered_sams/{sample}/{sample}_3_distance.sam'
	output:
		c='coverage/500kb/{sample}_coverage_500kb.bed'
	shell:
		'samtools view -bS {input.r} | bedtools coverage '
		'-abam - -b ' + config['data_dir'] + '/binning/hg19_500kb.bed | '
		'grep -v chrM |  grep -v chrY | sort -k1,1V -k2,2n '
		'> {output.c}'

rule bin_40kb:
	input:
		r='filtered_sams/{sample}/{sample}_3_distance.sam'
	output:
		c='coverage/40kb/{sample}_coverage_40kb.bed'
	shell:
		'samtools view -bS {input.r} | bedtools coverage '
		'-abam - -b ' + config['data_dir'] + '/binning/hg19_40kb.bed | '
		'grep -v chrM |  grep -v chrY | sort -k1,1V -k2,2n '
		'> {output.c}'

rule copy_number_correct_500kb:
	input:
		c='coverage/500kb/{sample}_coverage_500kb.bed'
	output:
		o='copy_correct_coverage/500kb/{sample}_copy_correct_coverage_500kb.bed'
	shell:
		'scripts/copy_correct.py -i {input.c} '
		'-c ' + config['data_dir'] + '/copy_number/COSMIC_K562_copy_number_500kb.bedGraph '
		'> {output.o}'

rule copy_number_correct_40kb:
	input:
		c='coverage/40kb/{sample}_coverage_40kb.bed'
	output:
		o='copy_correct_coverage/40kb/{sample}_copy_correct_coverage_40kb.bed'
	shell:
		'scripts/copy_correct.py -i {input.c} '
		'-c ' + config['data_dir'] + '/copy_number/COSMIC_K562_copy_number_40kb.bedGraph '
		'> {output.o}'




		
