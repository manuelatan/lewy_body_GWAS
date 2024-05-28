##### CORTICAL LB GWAS REPLICATION FOR HUW AND LESLEY #####
#Created by: Manuela Tan
#Created: 19/01/2024
#Last updated: 19/01/2024
#Last update: started script
#WD: /Users/manuela/Documents/Work/Mayo/LB_GWAS_Lesley
#Aim: run replication for LB GWAS for Lesley and Huw in Mayo data
#Lesley's files were QCed and imputed by GP2 so following the pipeline at: https://github.com/GP2code/GenoTools
#Original files are chr$CHR_IMPUTED_FILTERED_UPDATED_TOASSESS.bed bim and fam
#The data is in hg38, files have already been imputed


### Combine chromosome files ###

	cd ../MayoLBD

	#Unzip all gz files - don't think plink can merge gz files
	gunzip *.gz

	#Make list of files
	#In the list of files, remove .bim suffix - so plink knows to refer to all binary files bed bim and fam
	ls *.bim | sed 's/.bim//g' > list_files.txt

	#Remove chr1 filename from list_files.txt
	#This is because this will be the first file called in plink
	grep -v chr1_IMPUTED list_files.txt > list_files_final.txt

	#Merge in plink
	/Users/manuela/Documents/Work/software/plink/plink --bfile chr1_IMPUTED_FILTERED_UPDATED_TOASSESS \
	--merge-list list_files_final.txt \
	--make-bed \
	--out Mayo.ALL_CHROMOSOMES

### Combined files ###

	#8752750 variants
	#987 individuals
	#Not sure what imputation quality (R2/INFO) filters have been applied after imputation but the number of variants is similar to the NBB data (8652116 variants with R2 > 0.3 and MAF > 0.01)

### Data appears to be in hg38 build ###
	
	#Check gnomAD for the first few variants (by chr:pos) and there are variants with matching alleles in the hg38 (v3.1) of gnomAD but no variants in hg19 (v2.1) 
	#Also searched bim for common APOE SNPs rs7412 and rs429358 which are usually present in imputed datasets
	#These were present when using hg38 coordinates but not hg19 coordinates


### Extra QC before GWAS ###
	
	#Update sex and pheno
	#Sex is already in fam file and concordant with clinical data files
	#Keep only individuals with clinical data 980/987
	/Users/manuela/Documents/Work/software/plink/plink --bfile ../MayoLBD/Mayo.ALL_CHROMOSOMES \
	--pheno ../clinical/outputs/Mayo_imputed_pheno.tab \
	--keep ../clinical/outputs/Mayo_imputed_pheno.tab \
	--make-bed \
	--out Mayo.sex.pheno
	#8752750 variants and 961 people pass filters and QC.
	#Among remaining phenotypes, 356 are cases and 604 are controls. (1 phenotype is missing.)


	#Variant level QC
	/Users/manuela/Documents/Work/software/plink/plink --bfile Mayo.sex.pheno \
	--maf 0.01 \
	--autosome \
	--exclude /Users/manuela/Documents/Work/exclusion_regions_hg38.txt \
	--make-bed \
	--out Mayo.sex.pheno.maf.auto.exclusion
	#--exclude: 8752750 variants remaining.
	#85662 variants removed due to minor allele threshold(s)
	#8667088 variants and 961 people pass filters and QC.

	#Filter by genotyping rate
	/Users/manuela/Documents/Work/software/plink/plink --bfile Mayo.sex.pheno.maf.auto.exclusion \
	--geno 0.01 \
	--make-bed \
	--out Mayo.sex.pheno.maf.auto.exclusion.geno
	#0 variants removed due to missing genotype data (--geno).

	#HWE filtering
	/Users/manuela/Documents/Work/software/plink/plink --bfile Mayo.sex.pheno.maf.auto.exclusion.geno \
	--hwe 1e-8 \
	--make-bed \
	--out Mayo.sex.pheno.maf.auto.exclusion.geno.hwe
	#--hwe: 0 variants removed due to Hardy-Weinberg exact test.

	# IBD analysis - following GP2 Genotools pipeline

	/Users/manuela/Documents/Work/software/plink2/plink2 --bfile Mayo.sex.pheno.maf.auto.exclusion.geno.hwe \
	--hwe 0.0001 \
	--mac 2 \
	--make-pgen psam-cols=fid,parents,sex,pheno1,phenos \
	--out grm1

    # create table of related pairs
    /Users/manuela/Documents/Work/software/plink2/plink2 --pfile grm1 \
    --make-king-table \
    --make-king triangle bin \
    --king-table-filter 0.0884 \
    --out related_pairs

    # see if any samples are related (includes duplicates)
    /Users/manuela/Documents/Work/software/plink2/plink2 --pfile grm1 \
    --king-cutoff related_pairs 0.0884 \
    --out grm2
   	#no related pairs 


### Make PCs from imputed data ###
	#This is what Lesley did - not on directly genotyped data

 	/Users/manuela/Documents/Work/software/plink/plink --bfile Mayo.sex.pheno.maf.auto.exclusion.geno.hwe \
	--indep-pairwise 1000 10 0.02 \
	--autosome \
	--out imp.pruned_data

	/Users/manuela/Documents/Work/software/plink/plink --bfile Mayo.sex.pheno.maf.auto.exclusion.geno.hwe \
	--extract imp.pruned_data.prune.in \
	--make-bed \
	--out imp.PRUNED
	#23910 variants

	#Calculate PCAs on pruned data
	/Users/manuela/Documents/Work/software/plink/plink --bfile imp.PRUNED \
	--pca \
	--out imp.PCA
	#23910 variants and 961 people pass filters and QC.
	#Scree plot in R script


### GWAS ###

	#393 individuals have complete genetic and clinical data
	#Quite a lot of individuals are missing age at onset which is a covariate = only 393 out of 961
	/Users/manuela/Documents/Work/software/plink2/plink2 --bfile Mayo.sex.pheno.maf.auto.exclusion.geno.hwe \
	--ci 0.95 \
	--covar ../clinical/outputs/COVAR.txt \
	--covar-name sex,age_at_death,PC1,PC2,PC3 \
	--covar-variance-standardize \
	--glm hide-covar \
	--out MAYO_PATH

	#Make frequency file
	/Users/manuela/Documents/Work/software/plink/plink --bfile Mayo.sex.pheno.maf.auto.exclusion.geno.hwe \
	--freq \
	--out Mayo.sex.pheno.maf.auto.exclusion.geno.hwe.freq
	#961 individuals

	#Make frequency file but only with individuals included in the GWAS - because a lot were missing AAO
	/Users/manuela/Documents/Work/software/plink/plink --bfile Mayo.sex.pheno.maf.auto.exclusion.geno.hwe \
	--keep ../clinical/outputs/COVAR.txt \
	--freq \
	--out Mayo.sex.pheno.maf.auto.exclusion.geno.hwe.final_individuals.freq
	#393 individuals


	#Run GWAS without filtering for AAO
	#AAO was not included as a covariate in the original GWAS but the covariate file was filtered to drop individuals missing AAO
	#Made new covar file which does not apply this filter
	/Users/manuela/Documents/Work/software/plink2/plink2 --bfile Mayo.sex.pheno.maf.auto.exclusion.geno.hwe \
	--ci 0.95 \
	--covar ../clinical/outputs/COVAR_NOAAO.txt \
	--covar-name sex,age_at_death,PC1,PC2,PC3 \
	--covar-variance-standardize \
	--glm hide-covar \
	--out MAYO_PATH_NOAAO
	#958 individuals
	#2 individuals missing age at death, 1 missing Braak Lewy body stage




