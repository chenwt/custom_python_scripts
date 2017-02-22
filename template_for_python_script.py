#!/usr/bin/env python
desc="""Description"""

#imports
import os, pysam, resource, sys
from datetime import datetime
#from pybedtools.contrib.bigwig import bam_to_bigwig
import pybedtools

#required functions for the script
def sample_function

#main function that will be called when you ask for it
def main():
    #get your parser ready
    import argparse
    parser  = argparse.ArgumentParser(description=desc,formatter_class=argparse.RawTextHelpFormatter)
       
    #now add your arguments   
    parser.add_argument("-i", "--bam",     required=True,     
                        help="BAM file")
    parser.add_argument("-g", "--genome",  required=True,
                        help="genome FASTA file")
    parser.add_argument("-o", "--output",  required=True,
                        help="output stream   [stdout]")
    parser.add_argument("-s", "--strand", default="both", choices=("both", "+","-", "pos", "neg"), 
                        help="report coverage from + or - strand [%(default)s]")
    parser.add_argument("--scaling",       default=True,  action="store_false",
                        help="disable RPM scaling")
    
    #create your argument object
    o = parser.parse_args()
        
    #now run your function with the proper arguments
    #remember that the name will come from the --X in the add_argument function    
    bam2bigwig(o.bam, o.genome, o.output, o.strand, o.scaling, o.verbose)

if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    except IOError as e:
        sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
