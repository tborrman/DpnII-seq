#configfile:'config.yaml'

data_dir = '/home/tb37w/project/Research/digest/dpnII'

# rule sleepy:
# 	input:
# 		R1=data_dir + '/data/fastq/21OCT15_PE100_C7MDWAC/HBCRACKHiC-K562-DN-TD-R1/'+
# 			'HBCRACKHiC-K562-DN-TD-R1_GCCAAT_L008_R1_test.fastq',
# 		R2=data_dir + '/data/fastq/21OCT15_PE100_C7MDWAC/HBCRACKHiC-K562-DN-TD-R1/'+
# 			'HBCRACKHiC-K562-DN-TD-R1_GCCAAT_L008_R2_test.fastq'
# 	output:
# 		'mapped_reads/test.txt'
# 	shell:
# 		'echo {input.R1} > {output}'


rule bowtie_map:
	input:
		R1=data_dir + '/data/fastq/21OCT15_PE100_C7MDWAC/HBCRACKHiC-K562-DN-TD-R1/'+
			'HBCRACKHiC-K562-DN-TD-R1_GCCAAT_L008_R1_test.fastq',
		R2=data_dir + '/data/fastq/21OCT15_PE100_C7MDWAC/HBCRACKHiC-K562-DN-TD-R1/'+
			'HBCRACKHiC-K562-DN-TD-R1_GCCAAT_L008_R2_test.fastq'
	output:
		'mapped_reads/HBCRACKHiC-K562-DN-TD-R1_GCCAAT_L008_test.sam'
	log:
		'logs/bowtie_map/HBCRACKHiC-K562-DN-TD-R1_GCCAAT_L008_test.log'
	shell:
		'bowtie -S -p 10 -t -m 1 --chunkmbs 1000 -X 10000 '
		+ data_dir + '/data/indexes/hg19_build -1 {input.R1} '
		'-2 {input.R2} {output} 2> {log}'




