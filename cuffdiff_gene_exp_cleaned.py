#!/usr/bin/env python
desc="""Takes the gene_exp.diff file and cleans it up to show only a few columns from cuffdiff """

#imports
from datetime import datetime
import os, pysam, resource, sys
import pandas as pd
import numpy as np
from glob import glob

#required functions for the script
def merge(gene_diff_file,match_df,sample_id):
    dataframe = pd.read_csv(gene_diff_file,sep='\t')
    columns = ['gene_id','gene','locus','value_1','value_2','log2(fold_change)','q_value','biotype']
    merge_df = pd.merge(dataframe,match_df,how='left',on='gene_id')
    merge_df.to_csv(sample_id+'_cuffdiff_with_biotype.txt',sep='\t',na_rep='unknown',columns=columns,index=False)

#main function that will be called when you ask for it
def main():
    #set up any constants in the function
    #gene match list
    header = ['gene_id','external_gene_name','biotype']
    match_df = pd.read_csv('/media/sam/Data1/hnRNPM_iCLIP/iCLIP_analysis/2016_programmatic_iCLIP_analysis/getting_clipper_to_work/clipper_runs_working_with_individual_rt_stops/GO_analysis/ensembl_gene_id_external_gene_names_biotype',sep='\t',names=header)
   
    #get your parser ready
    import argparse
    parser  = argparse.ArgumentParser(description=desc,formatter_class=argparse.RawTextHelpFormatter)
       
    #now add your arguments   
    parser.add_argument("-d", "--diff",required=True,help="cuffdiff gene_exp.diff file")
    parser.add_argument("-s", "--sample",required=True,help="sample_name",type=str)
    
    #create your argument object
    o = parser.parse_args()
        
    #now run your function with the proper arguments
    #remember that the name will come from the --X in the add_argument function    
    merge(o.diff,match_df,o.sample)

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
