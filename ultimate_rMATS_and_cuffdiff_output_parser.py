#!/usr/bin/env python
desc="""The purpose of this script is to take rMATS output and cuffdiff output\
from a two RNA seq datasets and generate a complete output file that contains\
information on read count and gene expression for the rMATS data\
Inputs are the rMATS output folder, whether or not you want reads_on_target (rot)\
or junctions, FDR, dPSI, and average read counts if you want (defaults will be 0)\
You also need to put the sample name since that is not parsed from the rMATS output folder
"""

#imports
import os, sys
from datetime import datetime
import pandas as pd
import numpy as np
from glob import glob
from pandas import ExcelWriter

#for debug
#supp_folder = '/media/sam/Data2/annotations/supp_folder/'
#cuffdiff_output_folder = 'HF_cuffdiff_2_2_1_gencode_v24_lift_37/'
#rMATS_folder = 'HF_rMATS_3_2_5_gencode_v24_list_37_c_05/'

#required functions for the script
def parse_cuffdiff_output(cuffdiff_output_folder,supp_folder,sample_name):
    #put in the files that have the biotype, RNA binding protein, and TF ensemble gene IDs
    biotype_df = pd.read_table(supp_folder+'/biotype.txt',names=['GeneID','biotype'])
    rbp_df = pd.read_table(supp_folder+'/rbp.txt',names=['GeneID','rbp'])
    dbp_df = pd.read_table(supp_folder+'/dbp.txt',names=['GeneID','dbp'])
    
    #now read the gene_exp.diff file from cuffdiff
    dataframe = pd.read_csv(cuffdiff_output_folder+'gene_exp.diff',sep='\t')
    #change gene_id to GeneID and rename other columns as desired
    dataframe.rename(columns={'gene_id':'GeneID','value_1':'fpkm_sample_1','value_2':'fpkm_sample_2'},inplace=True)
    
    #drop any ENSG on GeneID
    dataframe['GeneID'] = dataframe.GeneID.str.split('.').str[0]
    
    merge_df = dataframe.merge(biotype_df,how='left',on='GeneID').merge(rbp_df,how='left',on='GeneID').merge(dbp_df,how='left',on='GeneID')
    
    columns = ['GeneID','gene','locus','fpkm_sample_1','fpkm_sample_2','log2(fold_change)','q_value','biotype','rbp','dbp']
    merge_df.to_csv('cuffdiff_parsing/'+sample_name+'_cuffdiff_with_biotype_rbp_dbp.txt',sep='\t',na_rep='False',columns=columns,index=False)
    merge_df.to_excel('cuffdiff_parsing/'+sample_name+'_cuffdiff_with_biotype_rbp_dbp.xlsx',na_rep='False',columns=columns,index=False)

def parse_rMATS_output(rMATS_folder,rot_or_junction,FDR,dPSI,average_read_count_min,sample_name):
    """Purpose of this is to take the rMATS data and put it indo pandas for filtering as novel and not novel based on the ID from the AS events as well as user called stringency filteres"""
    #first get the rot or junction file list
    data_list = sorted(glob(rMATS_folder+'MATS_output/*Count*'))
    for item in data_list:
        if rot_or_junction == 'rot':
            if 'JunctionCountOnly' in item:
                data_list.remove(item)
        if rot_or_junction == 'junction':
            if 'ReadsOnTarget' in item:
                data_list.remove(item)
    
    #now get the dataframes for the chosen samples
    df_list = []
    counter = 0
    for item in data_list:
        df = pd.read_table(item)
        df_list.append(df)
        counter +=1
        
    #for each of the 5 event types you must create an ID
    #just keep the splice ID at the end... it causes problems
    
    #A3SS
    df_list[0]['splice_id']= 'A3SS:'+df_list[0]['chr'].astype('str')+':'+df_list[0]['longExonStart_0base'].astype('str')+'-'+df_list[0]['longExonEnd'].astype('str')+':'+df_list[0]['shortES'].astype('str')+'-'+df_list[0]['shortEE'].astype('str')+':'+df_list[0]['flankingES'].astype('str')+'-'+df_list[0]['flankingEE'].astype('str')+':'+df_list[0]['strand'].astype('str')+':'+df_list[0]['geneSymbol'].astype('str')
    
    #A5SS
    df_list[1]['splice_id']='A5SS:'+df_list[1]['chr'].astype('str')+':'+df_list[1]['longExonStart_0base'].astype('str')+'-'+df_list[1]['longExonEnd'].astype('str')+':'+df_list[1]['shortES'].astype('str')+'-'+df_list[1]['shortEE'].astype('str')+':'+df_list[1]['flankingES'].astype('str')+'-'+df_list[1]['flankingEE'].astype('str')+':'+df_list[1]['strand'].astype('str')+':'+df_list[1]['geneSymbol'].astype('str')
    
    #MXE
    df_list[2]['splice_id']='MXE:'+df_list[2]['chr'].astype('str')+':'+df_list[2]['1stExonStart_0base'].astype('str')+'-'+df_list[2]['1stExonEnd'].astype('str')+':'+df_list[2]['2ndExonStart_0base'].astype('str')+'-'+df_list[2]['2ndExonEnd'].astype('str')+':'+df_list[2]['upstreamES'].astype('str')+'-'+df_list[2]['upstreamEE'].astype('str')+':'+':'+df_list[2]['downstreamES'].astype('str')+'-'+df_list[2]['downstreamEE'].astype('str')+':'+df_list[2]['strand'].astype('str')+':'+df_list[2]['geneSymbol'].astype('str')
    
    #RI
    df_list[3]['splice_id']='RI:'+df_list[3]['chr'].astype('str')+':'+df_list[3]['riExonStart_0base'].astype('str')+'-'+df_list[3]['riExonEnd'].astype('str')+':'+df_list[3]['upstreamES'].astype('str')+'-'+df_list[3]['upstreamEE'].astype('str')+':'+df_list[3]['downstreamES'].astype('str')+'-'+df_list[3]['downstreamEE'].astype('str')+':'+df_list[3]['strand'].astype('str')+':'+df_list[3]['geneSymbol'].astype('str')
    
    #SE
    df_list[4]['splice_id']='SE:'+df_list[4]['chr'].astype('str')+':'+df_list[4]['exonStart_0base'].astype('str')+'-'+df_list[4]['exonEnd'].astype('str')+':'+df_list[4]['upstreamES'].astype('str')+'-'+df_list[4]['upstreamEE'].astype('str')+':'+df_list[4]['downstreamES'].astype('str')+'-'+df_list[4]['downstreamEE'].astype('str')+':'+df_list[4]['strand'].astype('str')+':'+df_list[4]['geneSymbol'].astype('str')
        
    #now add splicing direction
    def splicing_direction(row):
        if row['IncLevelDifference'] > 0:
            direction = 'skip'
        else:
            direction = 'include'
        return direction
    for item in df_list:
        item['splicing_direction'] = item.apply(splicing_direction,axis=1)
    
    #now add sample 1 and sample2 average read counts
    
    if rot_or_junction == 'rot':
        for item in df_list:
            ic_sample_1 = []
            sc_sample_1 = []
            ic_sample_2 = []
            sc_sample_2 = []
            for row in item['IC_SAMPLE_1']:
                a_list = str(row).split(',')
                a_num = [int(x) for x in a_list]
                mean = np.mean(a_num)
                ic_sample_1.append(mean)
            for row in item['SC_SAMPLE_1']:
                a_list = str(row).split(',')
                a_num = [int(x) for x in a_list]
                mean = np.mean(a_num)
                sc_sample_1.append(mean)
            for row in item['IC_SAMPLE_2']:
                a_list = str(row).split(',')
                a_num = [int(x) for x in a_list]
                mean = np.mean(a_num)
                ic_sample_2.append(mean)
            for row in item['SC_SAMPLE_2']:
                a_list = str(row).split(',')
                a_num = [int(x) for x in a_list]
                mean = np.mean(a_num)
                sc_sample_2.append(mean)

            item['IC_SAMPLE_1_mean'] = ic_sample_1
            item['SC_SAMPLE_1_mean'] = sc_sample_1
            item['SAMPLE_1_mean_total'] = item.IC_SAMPLE_1_mean + item.SC_SAMPLE_1_mean
            item['IC_SAMPLE_2_mean'] = ic_sample_2
            item['SC_SAMPLE_2_mean'] = sc_sample_2
            item['SAMPLE_2_mean_total'] = item.IC_SAMPLE_2_mean + item.SC_SAMPLE_2_mean
    if rot_or_junction == 'junction':
        for item in df_list:
            ic_sample_1 = []
            sc_sample_1 = []
            ic_sample_2 = []
            sc_sample_2 = []
            for row in item['IJC_SAMPLE_1']:
                a_list = str(row).split(',')
                a_num = [int(x) for x in a_list]
                mean = np.mean(a_num)
                ic_sample_1.append(mean)
            for row in item['SJC_SAMPLE_1']:
                a_list = str(row).split(',')
                a_num = [int(x) for x in a_list]
                mean = np.mean(a_num)
                sc_sample_1.append(mean)
            for row in item['IJC_SAMPLE_2']:
                a_list = str(row).split(',')
                a_num = [int(x) for x in a_list]
                mean = np.mean(a_num)
                ic_sample_2.append(mean)
            for row in item['SJC_SAMPLE_2']:
                a_list = str(row).split(',')
                a_num = [int(x) for x in a_list]
                mean = np.mean(a_num)
                sc_sample_2.append(mean)

            item['IC_SAMPLE_1_mean'] = ic_sample_1
            item['SC_SAMPLE_1_mean'] = sc_sample_1
            item['SAMPLE_1_mean_total'] = item.IC_SAMPLE_1_mean + item.SC_SAMPLE_1_mean
            item['IC_SAMPLE_2_mean'] = ic_sample_2
            item['SC_SAMPLE_2_mean'] = sc_sample_2
            item['SAMPLE_2_mean_total'] = item.IC_SAMPLE_2_mean + item.SC_SAMPLE_2_mean
    
    #get rid of everything past the . in the GeneID so you are just searching for ENSG# and not versions
    for item in df_list:
        item['GeneID'] = item.GeneID.str.split('.').str[0]

    #Now you need to calculate the average PSI because it gives it in comma separated values
    #Note - you need to have a way to deal with NA values in some samples. Just throw it out.

    def get_mean_psi_from_dataframe_inc_level_1(row):
        #need to account for non replicated samples - they should be read as float whereas those with multiple samples are strings?
        if type(row['IncLevel1']) != float:
            psi_list = row['IncLevel1'].split(',')
            #need to account for multiple NAs that could show up
            for item in psi_list:
                if 'NA' in psi_list:
                    psi_list.remove('NA')
            psi_list_float = [float(i) for i in psi_list]
            psi_mean = np.mean(psi_list_float)
        else:
            psi_mean = row['IncLevel1']
        return(psi_mean)

    def get_mean_psi_from_dataframe_inc_level_2(row):
        #need to account for non replicated samples - they should be read as float whereas other are strings
        if type(row['IncLevel2']) != float:
            psi_list = row['IncLevel2'].split(',')
            #need to account for multiple NAs that could show up
            for item in psi_list:
                if 'NA' in psi_list:
                    psi_list.remove('NA')
            psi_list_float = [float(i) for i in psi_list]
            psi_mean = np.mean(psi_list_float)
        else:
            psi_mean = row['IncLevel2']
        return(psi_mean)

    for item in df_list:
        item['IncLevel1_mean'] = item.apply(get_mean_psi_from_dataframe_inc_level_1,axis=1)
        item['IncLevel2_mean'] = item.apply(get_mean_psi_from_dataframe_inc_level_2,axis=1)
    
    #Now add in the cuffdiff cleaned file on gene expression of these samples. Label the gene expression,biotype, rbp, or dbp identity
    
    cuffdiff_df = pd.read_table('cuffdiff_parsing/'+sample_name+'_cuffdiff_with_biotype_rbp_dbp.txt')
    
    merged_list = []
    for item in df_list:
        merge_df = item.merge(cuffdiff_df,how='left',on='GeneID')
        merged_list.append(merge_df)
        
    #now filter by FDR and dPSI and readcounts
    #require readcount filter that at least sample 1 AND sample 2 must meet readcount filter because you must be able to calculate a reliable delta psi
    filtered_list = []
    for item in merged_list:
        item = item[item.FDR <= FDR]
        item = item[item.IncLevelDifference.abs() >= dPSI]
        item = item[item.SAMPLE_1_mean_total > average_read_count_min]
        item = item[item.SAMPLE_2_mean_total > average_read_count_min]
        filtered_list.append(item)

    #add a splicing type column to indicate the splicing type
    name_list = ['A3SS','A5SS','MXE','RI','SE']
    counter = 0
    for item in merged_list:
        item['splicing_type'] = name_list[counter]
        counter += 1            
        
    #now get the novel event ids in a list by index
    as_events = sorted(glob(rMATS_folder+'/ASEvents/*'))
    id_list = []
    for item in as_events:
        if 'novel' in item:
            id_list.append(pd.read_table(item,usecols=[0]))
    
    #now filter each dataframe by the novel list
    
    novel_dataframes = []
    annotated_dataframes = []
    counter = 0
    for item in filtered_list:
        novel = item.loc[item.ID.isin(id_list[counter].ID)]
        novel_dataframes.append(novel)       
        annotated = item.loc[-item.ID.isin(id_list[counter].ID)]
        annotated_dataframes.append(annotated)
        counter += 1

    #now combine everything for all, novel, and annotated into one master file and then save it.
    #problem, the columns are all sorted out of wack - I might have to fix this later.
    #all_df = pd.concat(filtered_list,ignore_index=True)
    #all_df.to_csv('rMATS_parsing/'+sample_name+'_all_master'+'_FDR_'+str(FDR)+'_dPSI_'+str(dPSI)+'_read_cutoff_'+str(average_read_count_min)+'.txt',index=False,sep='\t')
    
    #novel_df = pd.concat(novel_dataframes,ignore_index=True)
    #novel_df.to_csv('rMATS_parsing/'+sample_name+'_novel_master'+'_FDR_'+str(FDR)+'_dPSI_'+str(dPSI)+'_read_cutoff_'+str(average_read_count_min)+'.txt',index=False,sep='\t')
    
    #annotated_df = pd.concat(annotated_dataframes,ignore_index=True)
    #annotated_df.to_csv('rMATS_parsing/'+sample_name+'_annotated_master'+'_FDR_'+str(FDR)+'_dPSI_'+str(dPSI)+'_read_cutoff_'+str(average_read_count_min)+'.txt',index=False,sep='\t')
        
    #now save everything as novel and non_novel
    #save as excel as well with different splice type in each sheet
    name_list = ['A3SS','A5SS','MXE','RI','SE']
    counter = 0
    writer = ExcelWriter('rMATS_parsing/'+sample_name+'_all_'+'FDR_'+str(FDR)+'_dPSI_'+str(dPSI)+'_read_cutoff_'+str(average_read_count_min)+'.xlsx')
    for item in filtered_list:
        item.to_csv('rMATS_parsing/'+sample_name+'_all_'+name_list[counter]+'_FDR_'+str(FDR)+'_dPSI_'+str(dPSI)+'_read_cutoff_'+str(average_read_count_min)+'.txt',index=False,sep='\t')
        item.to_excel(writer,index=False,sheet_name=name_list[counter],na_rep='False')
        counter += 1
    counter = 0
    writer = ExcelWriter('rMATS_parsing/'+sample_name+'_novel_'+'FDR_'+str(FDR)+'_dPSI_'+str(dPSI)+'_read_cutoff_'+str(average_read_count_min)+'.xlsx')
    for item in novel_dataframes:
#         item.add_prefix(sample_name+'_')
        item.to_csv('rMATS_parsing/'+sample_name+'_'+'novel_'+name_list[counter]+'_FDR_'+str(FDR)+'_dPSI_'+str(dPSI)+'_read_cutoff_'+str(average_read_count_min)+'.txt',index=False,sep='\t')
        item.to_excel(writer,index=False,sheet_name=name_list[counter],na_rep='False')
        counter += 1
    counter = 0
    writer = ExcelWriter('rMATS_parsing/'+sample_name+'_annotated_'+'FDR_'+str(FDR)+'_dPSI_'+str(dPSI)+'_read_cutoff_'+str(average_read_count_min)+'.xlsx')
    for item in annotated_dataframes:
#         item.add_prefix(sample_name+'_')
        item.to_csv('rMATS_parsing/'+sample_name+'_'+'annotated_'+name_list[counter]+'_FDR_'+str(FDR)+'_dPSI_'+str(dPSI)+'_read_cutoff_'+str(average_read_count_min)+'.txt',index=False,sep='\t')
        item.to_excel(writer,index=False,sheet_name=name_list[counter],na_rep='False')
        counter += 1

    #now move everything into an all, novel, and annotated
    os.system('mkdir rMATS_parsing/'+sample_name+'_all')
    os.system('mv rMATS_parsing/*all*.txt rMATS_parsing/'+sample_name+'_all')
    os.system('mv rMATS_parsing/*all*.xlsx rMATS_parsing/'+sample_name+'_all')
    os.system('mkdir rMATS_parsing/'+sample_name+'_novel')
    os.system('mv rMATS_parsing/*novel*.txt rMATS_parsing/'+sample_name+'_novel')
    os.system('mv rMATS_parsing/*novel*.xlsx rMATS_parsing/'+sample_name+'_novel')
    os.system('mkdir rMATS_parsing/'+sample_name+'_annotated')
    os.system('mv rMATS_parsing/*annotated*.txt rMATS_parsing/'+sample_name+'_annotated')
    #Cant get the last annotated xlsx file into the directory. It doesn't matter. I don't use this file normally 
    #os.system('chmod 775 rMATS_parsing/*xlsx') 
    #os.system('mv rMATS_parsing/*annotated*.xlsx rMATS_parsing/'+sample_name+'_annotated')   


#main function that will be called when you ask for it
def main():
    #get your parser ready
    import argparse
    parser  = argparse.ArgumentParser(description=desc,formatter_class=argparse.RawTextHelpFormatter)
       
    #now add your arguments   
    parser.add_argument("-rm", "--rmats_folder",     required=True,     
                        help="path_to_rMATS_folder")
    parser.add_argument("-cf", "--cuffdiff_folder",  required=True,
                        help="path_to_cuffdiff_folder")
    parser.add_argument("-supp", "--supp_folder",  required=True,
                        help="path to folder where biotype rbp dbp text files are")
    parser.add_argument("-roj", "--rot_or_junction",  required=False,
                        help="type rot or junction",default='junction',type=str)
    parser.add_argument("-fdr", "--fdr",  required=False,
                        help="FDR",default=1,type=float)
    parser.add_argument("-reads", "--average_read_count_min",  required=False,
                        help="Minimum number of inclusion + skipping reads per sample",default=0,type=float)
    parser.add_argument("-dpsi", "--dpsi",  required=False,
                        help="Desired Delta PSI",default=0,type=float)
    parser.add_argument("-s", "--sample_name",  required=False,
                        help="Type the sample name that will append to output",default='sample',type=str)

    #create your argument object
    o = parser.parse_args()
        
    #now run your function with the proper arguments
    #remember that the name will come from the --X in the add_argument function    
    cmd = 'mkdir rMATS_parsing'
    os.system(cmd)
    cmd = 'mkdir cuffdiff_parsing'
    os.system(cmd)

    parse_cuffdiff_output(o.cuffdiff_folder,o.supp_folder,o.sample_name)

    parse_rMATS_output(o.rmats_folder,o.rot_or_junction,o.fdr,o.dpsi,o.average_read_count_min,o.sample_name)

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
