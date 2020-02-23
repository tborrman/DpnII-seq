#!/usr/bin/env python
import argparse
import sys
import os
import timeit
import bisect


# Flush STOUT continuously
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

parser = argparse.ArgumentParser(description='Write list of distances to closest FatI site for fragment '+
	'(distance reported is the minimum distance from either fragment end)')
parser.add_argument('-s', help='SAM file with alignment reads (bowtie_output_mapped.sam)', type=str, required=True)
parser.add_argument('-o', help='output zero-distance and nonzero-distance fragment files', type= bool, default= False)
parser.add_argument('-c', help='cutoff for nucleotide distance to restriction site (ex. 3)', type=int, default=3)
parser.add_argument('-w', help='wildcard for sample (ex. HBCRACKHiC-K562-DN-TD-R1_GCCAAT_L008_test)', type=str, required=True)
parser.add_argument('-r', help='read length (ex. 100)', type=int, required=True)
parser.add_argument('-d', help='path to data directory (ex. /home/tb37w/project/Research/digest/dpnII/data)',
 type=str, required=True)
args = parser.parse_args()

def make_fatI_site_dict(site_path):
	'''
	Make dictionary of fatI sites

	Args: 
		site_path: string path to fatI site data

	Returns:
		fatI_sites: dictionary with chromosomes 
		for keys and sorted list of coordinates of
		start and ends of fatI sites for values

	'''
	fatI_sites = {}
	chroms = map(str, range(1,23)) + ['X','Y']
	for chrom in chroms:
		DPN = open(site_path + '/fatI_sites/FatI_sites_chr' + chrom + '.bed', 'r')
		coords = []
		for line in DPN:
			# Change coordinates from BED file format
			coord1 = int(line.split()[1]) + 1
			coord2 = int(line.split()[2])
			coords.append(coord1)
			coords.append(coord2)
		fatI_sites['chr'+chrom] = coords
		DPN.close()
	#print 'Processed FatI site coordinate files'
	return fatI_sites

def write_distance_files(l1, l2, Zf, nonZf, min_d):
	'''
	Write fragments into files with zero distance to fatI site
	fragments and nonzero distance to fatI site fragments

	Args:
		l1: string read 1 (SAM line)
		l2: string read 2 (SAM line)
		Zf: filehandle zero distance fragments
		nonZf: filehandle non zero distance fragments
		min_d: minimum distance to fatI site for fragment
	'''
	if min_d == 0:
		Zf.write(l1 + l2)
	else:
		nonZf.write(l1 + l2)

def write_cutoff_file(l1, l2, cf, min_d, c):
	'''
	Write fragments into file for fragments which are less than
	or equal to c distance away from a fatI site

	Args: 
		l1: string read 1 (SAM line)
		l2: string read 2 (SAM line)
		cf: filehandle cutoff distance fragments
		min_d: minimum distance to fatI site for fragment
		c: cutoff for distance to fatI site
	'''
	if min_d <= c:
		cf.write(l1 + l2)

def get_nearest_fatI_site(cut_sites, p):
	'''
	Use bisect algorithm to find nearest 
	fatI cut site to p 

	Args:
		cut_sites: sorted list of fatI cut sites 
		p: int base pair position

	Returns:
		n: int nearest fatI cut site to p

	'''
	i = bisect.bisect_left(cut_sites, p)
	if i == len(cut_sites):
		n = cut_sites[-1]
	else:
		n = min(cut_sites[i-1], cut_sites[i], key=lambda x: abs(x-p))
	return n

def main():
	# Make FatI site dictionary 
	fatI_sites = make_fatI_site_dict(args.d)
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
			# if counter%100 == 0:
			# 	print 'On fragment: '+ str(counter)
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
			end = int(line2.split()[3]) + (args.r - 1)
			# Find nearest fatI cut site
			nearest_start = get_nearest_fatI_site(fatI_sites[chrom1], start)
			nearest_end = get_nearest_fatI_site(fatI_sites[chrom1], end)
			# Minimum distance to fatI sites from either end
			min_distance = min(abs(start-nearest_start), abs(end-nearest_end))
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
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print 'Running time: ' +  str(stop - start)

