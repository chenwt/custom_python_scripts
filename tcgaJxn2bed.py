'''This file will convert TCGA RNASeq v2 Junction Data to the junction format for jSplice

    jSplice is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    jSplice is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with jSplice.  If not, see <http://www.gnu.org/licenses/>
    
    (C) copyright ETH Zurich; 2014 Institute for Molecular Health Sciences, Krek group
    Developer: Yann Christinat
    
'''


# Print the mapped reads from TCGA RNASeq v2 to a BED6 format
import sys,argparse,os.path

argparser = argparse.ArgumentParser(description='Converts a TCGA RNASeq v2 junction file into a 6-column BED file.')
argparser.add_argument('-f','--tcga_junction_file',help='STAR junction file',required=True)
argparser.add_argument('-o','--outfile',help='Output file name. (If omitted, the same file name but with a .bed extension will be used.)')

args = argparser.parse_args()

if not os.path.exists(args.tcga_junction_file):
	print 'ERROR file '+args.tcga_junction_file+' does not exists!'
	sys.exit()
	
if not args.outfile:
	l=len(args.tcga_junction_file)
	ext=args.tcga_junction_file[l-4:]
	if ext=='.csv' or ext=='.tab' or ext=='.txt' or ext=='.tsv':
		args.outfile = args.tcga_junction_file[:l-4]+'.bed'
	else:
		args.outfile = args.tcga_junction_file+'.bed'

out=open(args.outfile,'w')
for l in open(args.tcga_junction_file,'r'):
	if 'junction' in l:
		continue
	else:
		frags=l.strip().split('\t')
		chr_parse = frags[0].split(',')
		up_parse = chr_parse[0].split(':')
		down_parse = chr_parse[1].split(':')
		chrom = up_parse[0]
		up_coord = up_parse[1]
		down_coord = down_parse[1]
		cnt = frags[1]
		strand = up_parse[2]
			
		out.write("\t".join([chrom, up_coord, down_coord, 'jxn', cnt, strand]) + '\n')
		
out.close()
print 'Output written in '+args.outfile
