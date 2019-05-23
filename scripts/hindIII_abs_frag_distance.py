#!/usr/bin/env python
import argparse
import sys
import os
# Flush STOUT continuously
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

parser = argparse.ArgumentParser(description='Write list of distances to closest HindIII site for fragment '+
	'(distance reported is the minimum distance from either fragment end)')
parser.add_argument('-s', help='SAM file with alignment reads (bowtie_output_mapped.sam)', type=str, required=True)
parser.add_argument('-o', help='output zero-distance and nonzero-distance fragment files', type= bool, default= False)
parser.add_argument('-c', help='cutoff for nucleotide distance to restriction site (ex. 3)', type=int, default=3)
parser.add_argument('-w', help='wildcard for sample (ex. HBCRACKHiC-K562-HN-TD-R1_ACAGTG_L008)', type=str, required=True)
args = parser.parse_args()

def write_distance_files(l1, l2, Zf, nonZf, min_d):
	'''
	Write fragments into files with zero distance to dpnII site
	fragments and nonzero distance to dpnII site fragments
	Args:
		l1: string read 1 (SAM line)
		l2: string read 2 (SAM line)
		Zf: filehandle zero distance fragments
		nonZf: filehandle non zero distance fragments
		min_d: minimum distance to dpnII site for fragment
	'''
	if min_d == 0:
		Zf.write(l1 + l2)
	else:
		nonZf.write(l1 + l2)

def write_cutoff_file(l1, l2, cf, min_d, c):
	'''
	Write fragments into file for fragments which are less than
	or equal to c distance away from a dpnII site
	Args: 
		l1: string read 1 (SAM line)
		l2: string read 2 (SAM line)
		cf: filehandle cutoff distance fragments
		min_d: minimum distance to dpnII site for fragment
		c: cutoff for distance to dpnII site
	'''
	if min_d <= c:
		cf.write(l1 + l2)


def main():

	# Store hindIII coordinates in dictionary
	hindIII_sites = {}
	chroms = map(str, range(1,23)) + ['X','Y']
	for chrom in chroms:
		DPN = open('/home/tb37w/project/Research/digest/dpnII/test/hg19/HN/hindIII/hindIII_sites_chr'+chrom+'.bed', 'r')
		coords = []
		for line in DPN:
			# Change coordinates from BED file format
			coord1 = int(line.split()[1]) + 1
			coord2 = int(line.split()[2])
			# Get coordinates for nucleotides after cut site on 5' end and before cut site on 3' end
			coord1 = coord1 + 1
			coord2 = coord2 - 1 
			coords.append(coord1)
			coords.append(coord2)
		hindIII_sites['chr'+chrom] = coords
		DPN.close()
	print 'Processed coordinate files'

	# Parse SAM file
	SAM = open(args.s, 'r')
	index = args.s[-4:]
	OUT = open('distances/'+args.w + '/abs_frag_distances_'+index,'w')
	if args.o:
		ZERO_DIST_OUT = open('zero_distance_'+index, 'w')
		NONZERO_DIST_OUT = open('nonzero_distance_'+index, 'w')
	CUTOFF_DIST_OUT = open('filtered_sams/'+args.w + '/'+ str(args.c) + '_distance_'+index, 'w')
	counter = 0
	while True:
		line1=SAM.readline()
		line2=SAM.readline()
		if not line2:
			break
		else:
			counter += 1
			if counter%100 == 0:
				print 'On fragment: '+ str(counter)
			# Get chromosome
			chrom1 = line1.split()[2]
			chrom2 = line2.split()[2]
			if chrom1 != chrom2:
				print 'ERROR: wrong chromosomes'
				sys.exit()
			# Ignore mitochondrial DNA
			if chrom1 == 'chrM':
				continue
			# Check alignment has pair
			if (int(line1.split()[8]) < 0) or (int(line2.split()[8]) > 0):
				print 'ERROR: no pair'
				sys.exit()
			# Get start and end of fragment
			start = int(line1.split()[3])
			end = int(line2.split()[3]) + 99
			# Find closest hindIII sites
			closest_start = min(hindIII_sites[chrom1], key=lambda x: abs(start-x))
			closest_end = min(hindIII_sites[chrom1], key=lambda x: abs(end-x))
			# Minimum distance to hindIII sites from either end
			min_distance = min(abs(start-closest_start), abs(end-closest_end))
			if args.o:
				write_distance_files(line1, line2, ZERO_DIST_OUT, NONZERO_DIST_OUT, min_distance)
			OUT.write(str(min_distance) + '\n')
			write_cutoff_file(line1, line2, CUTOFF_DIST_OUT, min_distance, args.c)
	SAM.close()
	OUT.close()
	CUTOFF_DIST_OUT.close()
	if args.o:
		ZERO_DIST_OUT.close()
		NONZERO_DIST_OUT.close()

if __name__ == '__main__':
	main()

