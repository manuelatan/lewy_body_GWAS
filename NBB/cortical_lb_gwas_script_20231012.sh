##### CORTICAL LB GWAS REPLICATION FOR HUW AND LESLEY #####
#Created by: Manuela Tan
#Created: 14/09/2023
#Last updated: 03/10/2023
#Last update: post imputation filtering and run case/control GWAS
#WD: /Users/manuela/Documents/Work/NeuroChip_NBB/LB_GWAS_Lesley
#Aim: run replication for LB GWAS for Lesley and Huw
#Lesley's files were QCed and imputed by GP2 so following the pipeline at: https://github.com/GP2code/GenoTools
#Original files are NeuroChip_NBB_raw.bed bim and fam, copied from original files in raw/
#Note that data is in hg19, whereas GP2 have data directly in hg38 from raw genotype data

### Original files ###
	#490132 variants
	#491 individuals

### QC on non-imputed data ###

	#Following pipeline and filter thresholds at https://github.com/GP2code/GenoTools
	#The functions are in https://github.com/GP2code/GenoTools/blob/main/QC/qc.py
	#The order of QC functions are in https://github.com/GP2code/GenoTools/blob/main/run_qc_pipeline.py

	#CHeck sample genotyping rate and heterozygosity
	/Users/manuela/Documents/Work/software/plink/plink --bfile NeuroChip_NBB_raw \
	--missing \
	--het \
	--out NeuroChip_NBB_raw.sampleqc

	# Sample-level pruning - remove samples with genotyping rate < 98%
	/Users/manuela/Documents/Work/software/plink/plink --bfile NeuroChip_NBB_raw \
	--mind 0.02 \
	--make-bed \
	--out NeuroChip_NBB_raw.sample_0.98
	#110 individuals remaining - 381 removed

	# Sample-level pruning - remove samples with genotyping rate < 97% as otherwise too many samples are removed
	/Users/manuela/Documents/Work/software/plink/plink --bfile NeuroChip_NBB_raw \
	--mind 0.03 \
	--make-bed \
	--out NeuroChip_NBB_raw.sample_0.97
	#No individuals removed

### Sex checking ###

	# Update sex
	/Users/manuela/Documents/Work/software/plink/plink --bfile NeuroChip_NBB_raw.sample_0.97 \
	--keep ./clinical/outputs/NBB_sex.tab \
	--update-sex ./clinical/outputs/NBB_sex.tab \
	--make-bed \
	--out NeuroChip_NBB_raw.sample_0.97.updated_sex
	#377 people remaining
	#Warning: 3554 het. haploid genotypes present (see NeuroChip_NBB_raw.sample_0.97.updated_sex.hh ); many commands treat these as missing.
	#Total genotyping rate in remaining samples is 0.976395.
	#490132 variants and 377 people pass filters and QC.

	# Sex checking stage 1 - on whole genome
	/Users/manuela/Documents/Work/software/plink/plink --bfile NeuroChip_NBB_raw.sample_0.97.updated_sex \
	--check-sex 0.25 0.75 \
	--maf 0.05 \
	--out sex_tmp1
	#--check-sex: 6493 Xchr and 0 Ychr variant(s) scanned, 1 problem detected.

	# Sex checking stage 2 - on chr 23
	/Users/manuela/Documents/Work/software/plink/plink --bfile NeuroChip_NBB_raw.sample_0.97.updated_sex \
	--chr 23 --from-bp 2699520 --to-bp 154931043 \
	--maf 0.05 \
	--geno 0.05 \
	--hwe 1E-5 \
	--check-sex  0.25 0.75 \
	--out sex_tmp2
	#--check-sex: 6319 Xchr and 0 Ychr variant(s) scanned, 1 problem detected.

	# Grab fails from .sexcheck files
	python
	import subprocess
	import argparse
	import pandas as pd
	import numpy as np
	import os
	import glob
	import shutil
	import sys

	sex1 = pd.read_csv(f'sex_tmp1.sexcheck', sep='\s+')
	sex_fail1 = sex1[sex1.STATUS=='PROBLEM']

	sex2 = pd.read_csv(f'sex_tmp2.sexcheck', sep='\s+')
	sex_fail2 = sex2[sex2.STATUS=='PROBLEM']

	# combine and output
	# sex_fail_df = sex_fail1.append(sex_fail2)
	sex_fail_df = pd.concat([sex_fail1, sex_fail2], ignore_index=True)
	sex_fail_ids = sex_fail_df.loc[:,['FID','IID']].drop_duplicates(subset=['FID','IID'])
	sex_fail_count = sex_fail_ids.shape[0]
	sex_fail_ids.to_csv("sex_fails.tab", sep='\t', header=True, index=False)

	quit()

	
	#Remove sex fail samples from geno
	/Users/manuela/Documents/Work/software/plink/plink --bfile NeuroChip_NBB_raw.sample_0.97.updated_sex \
	--remove sex_fails.tab \
	--make-bed \
	--out NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass


### Run ancestry methods ###

	#This differs from the GP2 pipeline - they use prediction to predict ancestry and training/testing models
	#Here I just do the simple way of merging with reference data and removing PC outliers
	#Following GP2 Bioinformatics course pipeline https://github.com/GP2-TNC-WG/GP2-Beginner-Bioinformatics-for-PD-Genetics/blob/master/Module_I.md#4

	#First extract common SNPs between NBB and HapMap - otherwise we get a very low genotyping rate
	/Users/manuela/Documents/Work/software/plink/plink --bfile NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass \
	--extract /Users/manuela/Documents/Work/reference/hapmap/hapmap3r2_CEU.CHB.JPT.YRI.no-at-cg-snps.txt \
	--make-bed \
	--out NBB.hapmap_SNPs
	#197289 SNPs remaining 
	#Using rsIDs

	#Merge NBB data with HapMap by SNP name (rsID)
	/Users/manuela/Documents/Work/software/plink/plink --bfile NBB.hapmap_SNPs \
	--bmerge /Users/manuela/Documents/Work/reference/hapmap/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps \
	--make-bed \
	--out hapmap3_bin_snplis 

	#Flip missnps
	/Users/manuela/Documents/Work/software/plink/plink --bfile NBB.hapmap_SNPs \
	--flip hapmap3_bin_snplis-merge.missnp \
	--make-bed \
	--out NBB.hapmap_SNPs.flipped
	#--flip: 98547 SNPs flipped.

	#Merge attempt2
	/Users/manuela/Documents/Work/software/plink/plink --bfile NBB.hapmap_SNPs.flipped \
	--bmerge /Users/manuela/Documents/Work/reference/hapmap/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps \
	--extract NBB.hapmap_SNPs.bim \
	--make-bed \
	--out hapmap3_bin_snplis 
	#Total genotyping rate is 0.99365.
	#197289 variants and 771 people pass filters and QC.
	#Just merging SNPs that are common across both datasets and merging by SNP name (rsIDs)

	#Prune SNPs
	/Users/manuela/Documents/Work/software/plink/plink --bfile hapmap3_bin_snplis \
	--geno 0.01 \
	--maf 0.05 \
	--indep-pairwise 50 5 0.5 \
	--out pruning

	#Extract pruned SNPs
	/Users/manuela/Documents/Work/software/plink/plink --bfile hapmap3_bin_snplis \
	--extract pruning.prune.in \
	--make-bed \
	--out pruned_data
	#136887 variants and 771 people pass filters and QC.

	#Calculate heterozygosity of pruned SNPs
	/Users/manuela/Documents/Work/software/plink/plink --bfile pruned_data \
	--het \
	--out prunedHet

	#Calculate PCAs
	/Users/manuela/Documents/Work/software/plink/plink --bfile pruned_data \
	--geno 0.01 \
	--pca 10 \
	--out pca
	#Total genotyping rate is 0.998808.
	#0 variants removed due to missing genotype data (--geno).
	#136887 variants and 771 people pass filters and QC.

	#Run PCA_script_CEU.R script
	#This writes a list of individuals who are >6 SDs away from the mean of the HapMap CEU samples for any of the first 10 PCs

	/Users/manuela/Documents/Work/software/plink/plink --bfile NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass \
	--remove PCA_outliers.txt \
	--make-bed \
	--out NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep
	#20 individuals removed
	#356 people remaining

### Relatedness ###

	# crete pfiles
	/Users/manuela/Documents/Work/software/plink2/plink2 --bfile NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep \
	--hwe 0.0001 \
	--mac 2 \
	--make-pgen \
	--out grm1

	# create table of related pairs
    /Users/manuela/Documents/Work/software/plink2/plink2 --pfile grm1 \
    --make-king-table \
    --king-table-filter 0.125 \
    --out NBB_pairs
    #No related individuals, will continue to next step

### Heterozygosity ###
	
	#Prune SNPs
	/Users/manuela/Documents/Work/software/plink2/plink2 --bfile NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep \
	--geno 0.01 \
	--maf 0.05 \
	--indep-pairwise 50 5 0.5 \
	--out het_tmp

	#Extract pruned SNPs
	/Users/manuela/Documents/Work/software/plink2/plink2 --bfile NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep \
	--extract het_tmp.prune.in \
	--make-bed \
	--out het_tmp2
	#138571 variants

	#Calculate heterozygosity on pruned SNPs
	/Users/manuela/Documents/Work/software/plink2/plink2 --bfile het_tmp2 \
	--het \
	--out het_tmp3

	#Write list of heterozygosity outliers
	python

	import pandas as pd
	import numpy as np

	het = pd.read_csv("het_tmp3.het", sep='\s+')
	het_outliers = het[((het.F <= -0.25) | (het.F >= 0.25))]
	outlier_count = het_outliers.shape[0]
	het_outliers.to_csv(f'outliers_out.tab', sep='\t', header=True, index=False)

	quit()

	#There are no heterozygosity outliers but will run command anyway
	/Users/manuela/Documents/Work/software/plink2/plink2 --bfile NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep \
	--remove outliers_out.tab \
	--make-bed \
	--out NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep.het
	#490132 variants and 356 individuals

### Variant QC ###

	/Users/manuela/Documents/Work/software/plink2/plink2 --bfile NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep.het \
	--geno 0.05 \
	--make-bed \
	--out NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep.het.geno
	#15700 variants removed
	#474432 variants remaining

	#Skipped missingness by case/control status

	#Missingness by haplotype using P > 1E-4
	/Users/manuela/Documents/Work/software/plink/plink --bfile NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep.het.geno \
	--test-mishap \
	--out hap_tmp1
	#--test-mishap: 2436 loci checked (471996 skipped).

	# read .missing.hap file and grab flanking snps for P <= 0.0001. write flanking snps to file to exclude w bash
	python
	import pandas as pd
	import numpy as np
	mis_hap = pd.read_csv(f'hap_tmp1.missing.hap', sep='\s+')
	mis_hap_snps = list(mis_hap[mis_hap.P <= 0.0001].loc[:,'FLANKING'].str.split('|'))
	snp_ls_df = pd.DataFrame({'snp':[rsid for ls in mis_hap_snps for rsid in ls]})
	snp_ls_df['snp'].to_csv(f'hap_tmp1.exclude',sep='\t', header=False, index=False)

	#Exclude SNPs by missingness by haplotype
	/Users/manuela/Documents/Work/software/plink2/plink2 --bfile NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep.het.geno \
	--exclude hap_tmp1.exclude \
	--make-bed \
	--out NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep.het.geno.hap
	#474432 variants originally
	#--exclude: 474383 variants remaining.

	#Did not do HWE filtering as this is only done on controls in GenoTools
	#Will apply a less stringent filter later after imputation


### Prepare for imputation hg19 ###

	#356 individuals
	#474383 variants
	#Total genotyping rate 0.999687

	#Run Will Rayne's checking tool
	wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.13.zip

	#Extract files from zip
	unzip HRC-1000G-check-bim-v4.2.13.zip

	#Get allele frequencies
	/Users/manuela/Documents/Work/software/plink/plink --bfile NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep.het.geno.hap \
	--freq \
	--out NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep.het.geno.hap.freq

	#Run checking tool
	perl HRC-1000G-check-bim.pl \
	-b NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep.het.geno.hap.bim \
	-f NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep.het.geno.hap.freq.frq \
	-r /Users/manuela/Documents/Work/reference/HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h

		##Matching to HRC#

	#	Position Matches
	#	 ID matches HRC 328786
	#	 ID Doesn't match HRC 10179
	#	 Total Position Matches 338965
	#	ID Match
	#	 Position different from HRC 7318
	#	No Match to HRC 111953
	#	Skipped (MT) 0
	#	Total in bim file 474383
	#	Total processed 458236#

	#	Indels 15963#

	#	SNPs not changed 45603
	#	SNPs to change ref alt 228695
	#	Strand ok 161200
	#	Total Strand ok 274298#

	#	Strand to change 155624
	#	Total checked 346283
	#	Total checked Strand 316824
	#	Total removed for allele Frequency diff > 0.2 3087
	#	Palindromic SNPs with Freq > 0.4 187#
	#

	#	Non Matching alleles 29272
	#	ID and allele mismatching 5119; where HRC is . 803
	#	Duplicates removed 184

### Liftover to hg38 ###

	#Make liftover file for NeuroChip data
	R
	library(data.table)
	library(dplyr)

	#Read in bim file
	bim <- fread("NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep.het.geno.hap.bim")
	
	#SNPs should be in chrN:start-end formats
	bim <- bim %>%
		mutate(chr = paste("chr", V1, sep = ""),
				start = V4,
				end = V4+1, 
				name = V2)

	liftover <- bim %>%
		select(chr, start, end, name)
	#Disable scientific notation
	options(scipen=999)

	#Save file
	write.table(liftover, "NBB.liftover_hg19_to_hg38.txt", quote = F, row.names = F, col.names = F, sep = "\t")
	q()
	n

	#Run liftOver
	/Users/manuela/Documents/Work/software/liftOver/liftOver NBB.liftover_hg19_to_hg38.txt /Users/manuela/Documents/Work/software/liftOver/hg19ToHg38.over.chain.gz output.bed unlifted.bed


	R
	library(data.table)
	library(dplyr)
	output <- fread("output.bed")
	old_map <- fread("NBB.liftover_hg19_to_hg38.txt")

	#Check whether chromosomes match
	#Merge by SNP name
	merged <- output %>%
		inner_join(old_map, by = "V4")

	merged %>%
		filter(V1.x != V1.y) %>%
		summarise(count = n())
	#57 SNPs with strange chromosome names - just remove these
	#chr0 have been removed

	#Write SNP name and new bp (starting bp)
	export_hg38 <- merged %>%
		filter(V1.x == V1.y) %>%
		select(V4, V2.x)

	write.table(export_hg38, "export_hg38.txt", quote = F, col.names = F, row.names = F)

	#473407 SNPs
	#Original was 474383 SNPs

	q()
	n


	#Update positions with hg38
	/Users/manuela/Documents/Work/software/plink/plink --bfile NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep.het.geno.hap \
	--extract export_hg38.txt \
	--update-map export_hg38.txt \
	--make-bed \
	--out NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep.het.geno.hap.hg38
	#473407 variants 
	#Total genotyping rate is 0.999688

### Prepare for imputation hg38 - to match GP2 and TopMed ###

	#Download updated version of Will Rayne's checking tool
	wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.3.0.zip

	#Extract files from zip
	unzip HRC-1000G-check-bim-v4.3.0.zip

	#Follow instructions at https://www.well.ox.ac.uk/~wrayner/tools/ to download TOPMed reference panel sites and convert to HRC formatted
	#See /Users/manuela/Documents/Work/reference/TOPMed

	#Make frequency file
	#Get allele frequencies
	/Users/manuela/Documents/Work/software/plink/plink --bfile NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep.het.geno.hap.hg38 \
	--freq \
	--out NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep.het.geno.hap.hg38.freq

	#Run checking tool
	#Had to edit HRC-1000G-check-bim.pl script to change gunzip to /usr/bin/gunzip
	perl HRC-1000G-check-bim.pl -b NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep.het.geno.hap.hg38.bim \
	-f NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep.het.geno.hap.hg38.freq.frq \
	-r /Users/manuela/Documents/Work/reference/TOPMed/PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz -h

		#	Matching to HRC#

	#	Position Matches
	#	 ID matches HRC 0
	#	 ID Doesn't match HRC 400962
	#	 Total Position Matches 400962
	#	ID Match
	#	 Position different from HRC 0
	#	No Match to HRC 55997
	#	Skipped (MT) 319
	#	Total in bim file 473407
	#	Total processed 457278#

	#	Indels 15591#

	#	SNPs not changed 44877
	#	SNPs to change ref alt 222626
	#	Strand ok 156211
	#	Total Strand ok 267503#

	#	Strand to change 153285
	#	Total checked 400962
	#	Total checked Strand 309496
	#	Total removed for allele Frequency diff > 0.2 5188
	#	Palindromic SNPs with Freq > 0.4 180#
	#

	#	Non Matching alleles 91286
	#	ID and allele mismatching 91286; where HRC is . 0
	#	Duplicates removed 538

### Make final files for imputation ###

	#Edited Run-plink.sh to refer to plink directory
	sh Run-plink.sh

	#VCFs have already been made

	#Then sort and add chr prefix for hg38 and zip
	sh
	for i in $(seq 1 23)
	do
	vcf-sort NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep.het.geno.hap.hg38-updated-chr$i.vcf | awk -F"\t" '{{if ($0 !~ /^#/) {{print "chr"$0}} else{{print $0}}}}' | bgzip -c > NBB_pre_impute_chr$i.vcf.gz
	done

### Extra QC before GWAS ###
	
	#R2 filtering done in ./TOPMed_imputed/postimpute_bcftools.sh and ./TOPMed_imputed/postimpute_plink.sh
	#Done on imputed data

	#Update sex and pheno
	/Users/manuela/Documents/Work/software/plink/plink --bfile ./TOPMed_imputed/NBB_allchromosomes.converted.R2_0.3.MAF_0.01 \
	--update-sex ./clinical/outputs/NBB_imputed_sex.tab \
	--pheno ./clinical/outputs/NBB_imputed_pheno.tab \
	--make-bed \
	--out LBD.sex.pheno
	#8652116 variants and 356 people pass filters and QC.
	#Among remaining phenotypes, 285 are cases and 71 are controls.

	#Variant level QC
	/Users/manuela/Documents/Work/software/plink/plink --bfile LBD.sex.pheno \
	--maf 0.01 \
	--autosome \
	--exclude /Users/manuela/Documents/Work/exclusion_regions_hg38.txt \
	--make-bed \
	--out LBD.sex.pheno.maf.auto.exclusion
	#8652116 variants and 356 people pass filters and QC.
	#0 variants removed due to minor allele threshold(s)
	#--exclude: 8652116 variants remaining.


	#Filter by genotyping rate
	/Users/manuela/Documents/Work/software/plink/plink --bfile LBD.sex.pheno.maf.auto.exclusion \
	--geno 0.01 \
	--make-bed \
	--out LBD.sex.pheno.maf.auto.exclusion.geno
	#0 variants removed due to missing genotype data (--geno).

	#HWE filtering
	/Users/manuela/Documents/Work/software/plink/plink --bfile LBD.sex.pheno.maf.auto.exclusion.geno \
	--hwe 1e-8 \
	--make-bed \
	--out LBD.sex.pheno.maf.auto.exclusion.geno.hwe
	#--hwe: 0 variants removed due to Hardy-Weinberg exact test.

	# IBD analysis -no duplicates (already removed previously)


### Make PCs from non-imputed data ###

	#This will overwrite the previous PCA files which were with the HapMap samples
 	/Users/manuela/Documents/Work/software/plink/plink --bfile NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep.het.geno.hap.hg38 \
	--indep-pairwise 1000 10 0.02 \
	--autosome \
	--out pruned_data

	/Users/manuela/Documents/Work/software/plink/plink --bfile NeuroChip_NBB_raw.sample_0.97.updated_sex.sexpass.PCA_keep.het.geno.hap.hg38 \
	--extract pruned_data.prune.in \
	--make-bed \
	--out PRUNED
	#10688 variants

	#Calculate PCAs on pruned data
	/Users/manuela/Documents/Work/software/plink/plink --bfile PRUNED \
	--pca \
	--out PCA
	#Total genotyping rate is 0.99965.
	#10688 variants and 356 people pass filters and QC.


### Make PCs from imputed data ###
	#This is what Lesley did - not on directly genotyped data

 	/Users/manuela/Documents/Work/software/plink/plink --bfile LBD.sex.pheno.maf.auto.exclusion.geno.hwe \
	--indep-pairwise 1000 10 0.02 \
	--autosome \
	--out imp.pruned_data

	/Users/manuela/Documents/Work/software/plink/plink --bfile LBD.sex.pheno.maf.auto.exclusion.geno.hwe \
	--extract imp.pruned_data.prune.in \
	--make-bed \
	--out imp.PRUNED
	#17431 variants

	#Calculate PCAs on pruned data
	/Users/manuela/Documents/Work/software/plink/plink --bfile imp.PRUNED \
	--pca \
	--out imp.PCA
	#17431 variants and 356 people pass filters and QC.
	#Scree plot in R script
	#The scree plots and variance explained in the nonimputed vs. imputed PCA eigenvals are very similar


### GWAS ###

	#315 individuals have complete genetic and clinical data
	#Some individuals are missing age at onset which is a covariate
	/Users/manuela/Documents/Work/software/plink2/plink2 --bfile LBD.sex.pheno.maf.auto.exclusion.geno.hwe \
	--ci 0.95 \
	--covar ./clinical/outputs/COVAR.txt \
	--covar-name sex,age_at_death,PC1,PC2,PC3 \
	--covar-variance-standardize \
	--glm hide-covar \
	--out LBD_PATH

	#Make frequency file
	/Users/manuela/Documents/Work/software/plink/plink --bfile LBD.sex.pheno.maf.auto.exclusion.geno.hwe \
	--freq \
	--out LBD.sex.pheno.maf.auto.exclusion.geno.hwe.freq

	#Lift back over to hg19 - because FUMA is not working with hg38 at the moment

	#Make file for liftover
	R
	library(data.table)
	library(dplyr)

	#Read in bim file
	data <- fread("LBD_PATH.PHENO1.glm.logistic.hybrid")
	
	#SNPs should be in chrN:start-end formats
	data <- data %>%
		mutate(chr = paste("chr", `#CHROM`, sep = ""),
				start = POS,
				end = POS+1, 
				name = ID)

	liftover <- data %>%
		select(chr, start, end, name)
	#Disable scientific notation
	options(scipen=999)

	#Save file
	write.table(liftover, "NBB.results.liftover_hg38_to_hg19.txt", quote = F, row.names = F, col.names = F, sep = "\t")
	q()
	n

	#Run liftOver
	/Users/manuela/Documents/Work/software/liftOver/liftOver NBB.results.liftover_hg38_to_hg19.txt /Users/manuela/Documents/Work/software/liftOver/hg38ToHg19.over.chain.gz output.bed unlifted.bed

	R
	library(data.table)
	library(dplyr)
	output <- fread("output.bed")
	results <- fread("LBD_PATH.PHENO1.glm.logistic.hybrid")

	#Check whether chromosomes match
	#Merge by SNP name
	merged <- output %>%
		inner_join(results, by = c("V4" = "ID"))


	#Make new results file for FUMA with hg19 coordinates
	results_hg19 <- merged %>%
		select(CHR = V1,
			BP = V2,
			REF, ALT, A1, OR, `LOG(OR)_SE`, P) %>%
		mutate(CHR = gsub("chr", "", CHR))

	#Remove odd chromosome names
	results_hg19_final <- results_hg19 %>% 
	filter(!is.na(as.numeric(CHR)))


	results_hg19_final %>%
	group_by(CHR) %>%
	summarise(count = n()) %>%
	print(n = 50)

	write.table(results_hg19_final, "results_hg19.txt", quote = F, col.names = T, row.names = F, sep = "\t")

	q()
	n


	#gzip file for fuma
	gzip -c results_hg19.txt > results_hg19.txt.gz

