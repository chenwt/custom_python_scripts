#!/bin/bash

#Here I will test different rmats_to_sashimi_plot files

#for the EMT ones, let's test the ones from li's set.

#with bam - try doing a subset of events such as replicates
#try snail
#rmats2sashimiplot --b1 bam_files/SNAIL1Aligned.sortedByCoord.out.bam,bam_files/SNAIL3Aligned.sortedByCoord.out.bam,bam_files/SNAIL5Aligned.sortedByCoord.out.bam --b2 bam_files/SNAIL2Aligned.sortedByCoord.out.bam,bam_files/SNAIL4Aligned.sortedByCoord.out.bam,bam_files/SNAIL6Aligned.sortedByCoord.out.bam -t SE -e emt_select_events_rMATS_output.txt --l1 SNAIL_EPI --l2 SNAIL_MES --exon_s 1 --intron_s 5 -o rmats2sashimiplot_output/snail_emt_events

#The above works ok but requires rMATS output. Maybe better to just use the GFF3 file from Gencode

#try GFF with a specific region of interest - don't worry yet about calculation of PSI
#this works pretty well - it even indicates the junction number
#where can i access the junction number?
rmats2sashimiplot --b1 bam_files/SNAIL1Aligned.sortedByCoord.out.bam,bam_files/SNAIL3Aligned.sortedByCoord.out.bam,bam_files/SNAIL5Aligned.sortedByCoord.out.bam --b2 bam_files/SNAIL2Aligned.sortedByCoord.out.bam,bam_files/SNAIL4Aligned.sortedByCoord.out.bam,bam_files/SNAIL6Aligned.sortedByCoord.out.bam -c chr11:+:129991987:129996905:./annotations/gencode.v24lift37.annotation.transcript_to_mRNA.gff3 --l1 SNAIL_EPI --l2 SNAIL_MES --exon_s 1 --intron_s 5 -o rmats2sashimiplot_output/snail_aplp2_exon_7
