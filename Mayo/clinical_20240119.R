##### LBD GWAS - REPLICATION FOR LESLEY #####
#Created: 19/01/2024
#Last updated: 19/01/2024
#Created by: Manuela Tan
#WD: /Users/manuela/Documents/Work/Mayo/clinical
#Last update: created script



#---Load packages---####

library(tidyverse)
library(data.table)
library(readxl)

#---Read in clinical data---####

#Note that the file called "cohort_data with IDs.xlsx" is the one you want to use - including all 980 samples
clinical <- read_xlsx("cohort_data with IDs.xlsx")

#There is also an additional file with AAO
# I also included a second file including the cases where age at onset is available. 
#In these files LBD-type is the measure of Lewy patholgy, while braak and braak_stage_num is Braak NFT stage. I have used the variable braak_stage_num for my analyses.
extra <- read_xlsx("Norway ID#S File Combinedoar_20221221_onset-FID.xlsx")

#---Read in genetic fam file---####
#This is after imputation performed at Mayo

fam <- read.table("../MayoLBD/Mayo.ALL_CHROMOSOMES.fam")
colnames(fam) <- c("FID", "IID", "m", "f", "sex", "pheno")

#Genders are already in the fam file
#Check concordance between clinical data and fam genders
fam_merged <- fam %>%
  inner_join(clinical, by = c("FID" = "fid"))

#Recode plink fam genders as text
fam_merged <- fam_merged %>%
  mutate(sex_fam = ifelse(sex.x == 1, "Male",
                          ifelse(sex.x == 2, "Female", NA)))

#Check match with clinical data file
fam_merged %>%
  filter(sex_fam !=sex.y)
#No mismatching

#Check any samples missing gender
fam_merged %>% 
  filter(is.na(sex_fam))

fam_merged %>% 
  filter(is.na(sex.y))
#No samples missing sex in either clinical or genetic files



##### Cortical LB GWAS #####


#Cortical LB GWAS
#Inclusion criteria: we are including all autopsy cases with a diagnosis of Lewy bodies with clinical data on sex and age at death
#Analysis: case-control type GWAS with sex, age at death and PC1 - PC3 as covariates 
#In the future, we can re-run this analysis including disease duration as a covariate - this will however decrease sample size due to missing data 
#If it's not too difficult, would it be possible to begin by including iLBD and cases with severe AD co-pathology.
#LBD is defined based on primary or secondary pathology diagnosis (so regardless of clinical diagnosis

#First exclude cases with Braak stage of 0
fam_merged_LB <- fam_merged %>% 
  filter(braak_stage_num > 0)
#19 samples removed

#Case = Braak stage 5-6
#Control = Braak stage 1-4
#Make file for plink Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
#Note there is one sample with braak_stage_num == "NA" as text, recode this as NA
fam_merged_LB <- fam_merged_LB %>% 
  mutate(case_status = ifelse(is.na(braak_stage_num), NA, 
                              ifelse(braak_stage_num == "NA", NA,
                                     ifelse(braak_stage_num >=5, 2,
                                            ifelse(braak_stage_num <5, 1, "check")))))

fam_merged_LB %>% 
  group_by(case_status, braak_stage_num) %>% 
  summarise(count = n())


#Export phenotype for plink 
pheno_export <- fam_merged_LB %>% 
  select(FID, IID, case_status)

write.table(pheno_export, "./outputs/Mayo_imputed_pheno.tab",
            quote = F,
            col.names = F, 
            row.names = F,
            sep = "\t")



#---Make scree plot from imputed data PCA---####

## Make scree plot 
#Read in PCA eigenvalues
eigenval <- read.table("../LB_GWAS_Lesley/imp.pca.eigenval", sep ="\t", header = F, stringsAsFactors = F)

# Update column names
colnames(eigenval)[1] <- "Eigenvalues"
eigenval$PC <- as.numeric(rownames(eigenval))
eigenval$VarianceExplained <- eigenval$Eigenvalues/sum(eigenval$Eigenvalues)*100

# Keeping only the first 10 PCs
eigenval2 <- head(eigenval,10)

# Generating the plot
scree <- ggplot(data = eigenval2, aes(x = PC, y = VarianceExplained)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "Principal Components", breaks = seq(0,10,1), limits = c(NA,10)) +
  scale_y_continuous(name = "Percent (%) Variance Explained", breaks = seq(0,50,5), limits = c(0,50)) +
  ggtitle("Scree Plot: \n LBD") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

scree
ggsave("./plots/scree_PCA_imputed.png")


# Generating the plot - WITHOUT y-axis scaling as it is difficult to see differences
scree <- ggplot(data = eigenval2, aes(x = PC, y = VarianceExplained)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "Principal Components", breaks = seq(0,10,1), limits = c(NA,10)) +
  #scale_y_continuous(name = "Percent (%) Variance Explained", breaks = seq(0,50,5), limits = c(0,50)) +
  ggtitle("Scree Plot: \n LBD") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

scree
ggsave("./plots/scree_PCA_imputed_noscale.png")


#---Make covariates file INCLUDING AAO---####

# Creating covariates file

Eigenvec <- as_tibble(read.table("../LB_GWAS_Lesley/imp.PCA.eigenvec", sep = "")) %>%
  dplyr::rename(FID = V1,
                IID = V2,
                PC1 = V3,
                PC2 = V4,
                PC3 = V5,
                PC4 = V6,
                PC5 = V7,
                PC6 = V8,
                PC7 = V9,
                PC8 = V10,
                PC9 = V11,
                PC10 = V12) %>%
  select(FID:PC3)

#Get age at onset from extra clinical data file
fam_merged_aao <- fam_merged_LB %>% 
  inner_join(extra, by = c("FID" = "fid"))
#Only 393 out of 961 cases have age at onset available

#Make covariate file
#Using PCs from imputed data PCA - following Lesley's pipeline
COVAR <- left_join(fam_merged_aao, Eigenvec, by = c("FID", "IID")) %>%
  select(FID, IID, 
         sex = sex.x,
         age_at_onset = "overall onset",
         age_at_death = age_at_death.x, PC1, PC2,PC3) %>%
  mutate(age_at_death = as.numeric(age_at_death)) %>% 
  na.omit() %>%
  distinct(IID, .keep_all = TRUE)

write.table(COVAR, "./outputs/COVAR.txt", col.names = T, row.names = F, quote = F)

#Check how many case/controls in final file after removing individuals with missing data
case_counts <- fam_merged_aao %>% 
  left_join(Eigenvec, by = c("FID", "IID")) %>%
  select(FID, IID, 
         case_status, braak_stage_num,
         sex = sex.x,
         age_at_onset = "overall onset",
         age_at_death = age_at_death.x) %>% 
  na.omit() %>% 
  group_by(case_status) %>% 
  summarise(count = n())


#Get demographics
demographics <- fam_merged_aao %>% 
  left_join(Eigenvec, by = c("FID", "IID")) %>%
  select(FID, IID, 
         case_status, braak_stage_num,
         sex = sex.x,
         age_at_onset = "overall onset",
         age_at_death = age_at_death.x) %>% 
  na.omit() %>% 
  group_by(case_status) %>% 
  mutate(age_at_death = as.numeric(age_at_death)) %>% 
  summarise(count = n(),
            mean_aao = mean(age_at_onset),
            sd_aao = sd(age_at_onset),
            mean_age_death = mean(age_at_death),
            sd_age_death = sd(age_at_death))


write.table(demographics, "./outputs/demographics_casecontrol.tab", 
            sep = "\t", quote = F, row.names = F, col.names = T)

#Sex counts
sex_counts <- fam_merged_aao %>% 
  left_join(Eigenvec, by = c("FID", "IID")) %>%
  select(FID, IID, 
         case_status, braak_stage_num,
         sex = sex.x,
         age_at_onset = "overall onset",
         age_at_death = age_at_death.x) %>% 
  na.omit() %>% 
  group_by(case_status, sex) %>% 
  summarise(count = n()) %>% 
  group_by(case_status) %>%
  mutate(percentage = count / sum(count) * 100)

write.table(sex_counts, "./outputs/sex_casecontrol.tab",
            quote = F, col.names = T, row.names = F, sep = "\t")

#---Make covariates file EXCLUDING AAO---####

#Make covariate file but without including age at onset
#Although age at onset is not included as a covariate, we filtered the file to drop individuals missing age at onset (following Lesley's script)
#In the Mayo data there are a lot of individuals without age at onset
COVAR_NOAAO <- left_join(fam_merged_LB, Eigenvec, by = c("FID", "IID")) %>%
  select(FID, IID, 
         sex = sex.x,
         age_at_death = age_at_death,
         PC1, PC2,PC3) %>%
  mutate(age_at_death = as.numeric(age_at_death)) %>% 
  na.omit() %>%
  distinct(IID, .keep_all = TRUE)

write.table(COVAR_NOAAO, "./outputs/COVAR_NOAAO.txt", col.names = T, row.names = F, quote = F)



#Check how many case/controls in final file after removing individuals with missing data
case_counts_noaao <- fam_merged_LB %>% 
  left_join(Eigenvec, by = c("FID", "IID")) %>%
  select(FID, IID, 
         case_status, braak_stage_num,
         sex = sex.x,
         age_at_death = age_at_death) %>% 
  na.omit() %>% 
  group_by(case_status) %>% 
  summarise(count = n())


#Get demographics
#Merged with extra clinical file to get AAO but have not removed individuals missing AAO or age at death
demographics_noaao <- fam_merged_LB %>% 
  left_join(Eigenvec, by = c("FID", "IID")) %>%
  left_join(extra, by = c("FID" = "fid")) %>% 
  select(FID, IID, 
         case_status, braak_stage_num,
         sex = sex.x,
         age_at_onset = "overall onset",
         age_at_death = age_at_death.x) %>% 
  mutate(age_at_death = as.numeric(age_at_death)) %>%
  filter(!is.na(case_status)) %>% 
  group_by(case_status) %>% 
  summarise(count = n(),
            mean_aao = mean(age_at_onset, na.rm = TRUE),
            sd_aao = sd(age_at_onset, na.rm = TRUE),
            mean_age_death = mean(age_at_death, na.rm = TRUE),
            sd_age_death = sd(age_at_death, na.rm = TRUE))


write.table(demographics_noaao, "./outputs/demographics_casecontrol_NOAAO.tab", 
            sep = "\t", quote = F, row.names = F, col.names = T)

#Sex counts
sex_counts_noaao <- fam_merged_LB %>% 
  left_join(Eigenvec, by = c("FID", "IID")) %>%
  select(FID, IID, 
         case_status, braak_stage_num,
         sex = sex.x,
         age_at_death = age_at_death) %>% 
  na.omit() %>% 
  group_by(case_status, sex) %>% 
  summarise(count = n()) %>% 
  group_by(case_status) %>%
  mutate(percentage = count / sum(count) * 100)

write.table(sex_counts_noaao, "./outputs/sex_casecontrol_NOAAO.tab",
            quote = F, col.names = T, row.names = F, sep = "\t")


