##### LBD GWAS - REPLICATION FOR LESLEY #####
#Created: 25/09/2023
#Last updated: 05/10/2023
#Created by: Manuela Tan
#WD: /Users/manuela/Documents/Work/NeuroChip_NBB/LB_GWAS_Lesley/clinical
#Last update: run demographics summary for case/control GWAS

#---Load packages---####

library(tidyverse)
library(data.table)

#---Read in clinical data---####

clinical <- read.csv("nbb_demographics_2023_09_19.txt")


#---Read in genetic fam file---####

fam <- read.table("../NeuroChip_NBB_raw.fam")
colnames(fam) <- c("FID", "IID", "m", "f", "sex", "pheno")

#Check match between clinical and genetic IDs
merged <- clinical %>%
  inner_join(fam, by = c("Case_ID" = "FID"))

#Check IDs in clinical data that are not in fam file
clinical %>%
  anti_join(fam, by = c("Case_ID" = "FID"))
#There are only 8 samples with IDs in the clinical file but not in the genetic fam file
#I checked the fam file and there do not seem to be any possible matches for these

#---Export genders for genetic data QC---####

sex_export <- merged %>%
  select(FID = Case_ID, IID, Sex)

write.table(sex_export, "./outputs/NBB_sex.tab",
            quote = F, col.names = F, row.names = F, sep = "\t")

#---Read in updated version of clinical data---####
#Sent by Jon-Anders 25/09/2023
#I have now made a new version of the NBB data where the only filtering is exclusion of cases with too limited clinical data (previously disucssed with Hanneke) and cases without Lewy pathology. 
#We now have 389 cases for the neuropathology GWAS. 


clinical_new <- read.csv("nbb_lp_demograpics_and_neuropathology_2023_09_25.csv")

#---Check diagnosis categories---####

clinical_new %>% 
  group_by(Neuropath_diagnosis) %>% 
  summarise(count = n())

clinical_new %>% 
  group_by(Braak_aSyn_stage) %>% 
  summarise(count = n())

#---Read in imputed fam file and merge---####


#Merge clinical data with fam file
imputed_fam <- read.table("../TOPMed_imputed/NBB_allchromosomes.converted.R2_0.3.MAF_0.01.fam")
colnames(imputed_fam) <- c("FID", "IID", "f", "m", "sex", "pheno")

clinical_fam <- clinical_new %>% 
  mutate(IID = paste(Case_ID, "_", Case_ID, sep = "")) %>% 
  inner_join(imputed_fam, by = c("IID"))

##### Cortical LB GWAS #####

#Cortical LB GWAS
#Inclusion criteria: we are including all autopsy cases with a diagnosis of Lewy bodies with clinical data on sex and age at death
#Analysis: case-control type GWAS with sex, age at death and PC1 - PC3 as covariates 
#In the future, we can re-run this analysis including disease duration as a covariate - this will however decrease sample size due to missing data 
#If it's not too difficult, would it be possible to begin by including iLBD and cases with severe AD co-pathology.
#LBD is defined based on primary or secondary pathology diagnosis (so regardless of clinical diagnosis


#Case = Braak stage 5-6
#Control = Braak stage 1-4
#Make file for plink Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
clinical_fam <- clinical_fam %>% 
  mutate(case_status = ifelse(Braak_aSyn_stage >=5, 2,
                              ifelse(Braak_aSyn_stage <5, 1, NA)))


clinical_fam %>% 
  group_by(case_status, Braak_aSyn_stage) %>% 
  summarise(count = n())

clinical_fam %>% 
  filter(is.na(Braak_aSyn_stage))

#Export sex for plink update sex
sex_export <- clinical_fam %>% 
  select(IID, FID, Sex)

write.table(sex_export, "./outputs/NBB_imputed_sex.tab",
            quote = F,
            col.names = F,
            row.names = F,
            sep = "\t")


#Export phenotype 
pheno_export <- clinical_fam %>% 
  select(IID, FID, case_status)

write.table(pheno_export, "./outputs/NBB_imputed_pheno.tab",
            quote = F,
            col.names = F, 
            row.names = F,
            sep = "\t")



#---Make scree plot from directly genotyped PCA---####

## Make scree plot 
#Read in PCA eigenvalues
eigenval <- read.table("../pca.eigenval", sep ="\t", header = F, stringsAsFactors = F)

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
ggsave("./plots/scree_PCA_nonimputed.png")

#---Make scree plot from imputed data PCA---####

## Make scree plot 
#Read in PCA eigenvalues
eigenval <- read.table("../imp.pca.eigenval", sep ="\t", header = F, stringsAsFactors = F)

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


#---Make covariates file---####

# Creating covariates file

Eigenvec <- as_tibble(read.table("../imp.PCA.eigenvec", sep = "")) %>%
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


#Make age at onset which is age at death minus disease duration
clinical_fam <- clinical_fam %>% 
  mutate(age_at_onset = Age_death - Disease_duration)

#Check individuals missing age at onset
clinical_fam %>% 
  filter(is.na(age_at_onset)) %>% 
  summarise(count = n())
#41 individuals missing age at onset - no other available data that can help to impute missing data

#Make covariate file
#Using PCs from imputed data PCA - following Lesley's pipeline
COVAR <- left_join(clinical_fam, Eigenvec, by = c("FID", "IID")) %>%
  select(FID, IID, 
         sex = Sex,
         age_at_onset,
         age_at_death = Age_death, PC1, PC2,PC3) %>%
  na.omit() %>%
  distinct(IID, .keep_all = TRUE)

write.table(COVAR, "./outputs/COVAR.txt", col.names = T, row.names = F, quote = F)

#Check how many case/controls in final file after removing individuals with missing data
case_counts <- clinical_fam %>% 
  left_join(Eigenvec, by = c("FID", "IID")) %>%
  select(FID, IID, 
         case_status, Braak_aSyn_stage,
         sex = Sex,
         age_at_onset,
         age_at_death = Age_death) %>% 
  na.omit() %>% 
  group_by(case_status) %>% 
  summarise(count = n())

#Get demographics
demographics <- clinical_fam %>% 
  left_join(Eigenvec, by = c("FID", "IID")) %>%
  select(FID, IID, 
         case_status, Braak_aSyn_stage,
         sex = Sex,
         age_at_onset,
         age_at_death = Age_death) %>% 
  na.omit() %>% 
  group_by(case_status) %>% 
  summarise(count = n(),
            mean_aao = mean(age_at_onset),
            sd_aao = sd(age_at_onset),
            mean_age_death = mean(age_at_death),
            sd_age_death = sd(age_at_death))

write.table(demographics, "./outputs/demographics_casecontrol.tab", 
            sep = "\t", quote = F, row.names = F, col.names = T)

#Sex counts
sex_counts <- clinical_fam %>% 
  left_join(Eigenvec, by = c("FID", "IID")) %>%
  select(FID, IID, 
         case_status, Braak_aSyn_stage,
         sex = Sex,
         age_at_onset,
         age_at_death = Age_death) %>% 
  na.omit() %>% 
  group_by(case_status, sex) %>% 
  summarise(count = n()) %>% 
  group_by(case_status) %>%
  mutate(percentage = count / sum(count) * 100)

write.table(sex_counts, "./outputs/sex_casecontrol.tab",
            quote = F, col.names = T, row.names = F, sep = "\t")

#---Plots by clinical diagnosis---####

final <- clinical_fam %>% 
  left_join(Eigenvec, by = c("FID", "IID")) %>%
  select(FID, IID, 
         case_status, Braak_aSyn_stage, Neuropath_diagnosis, Clinical_diagnosis,
         sex = Sex,
         age_at_onset,
         age_at_death = Age_death) %>% 
  na.omit() 



neuropath_dx_summary <- final %>% 
  group_by(Neuropath_diagnosis) %>% 
  summarise(count = n()) %>% 
  arrange(count) %>% 
  mutate(Neuropath_diagnosis=factor(Neuropath_diagnosis, Neuropath_diagnosis)) 

#Line plot of neuropath diagnoses counts
ggplot(data = neuropath_dx_summary, aes(x=Neuropath_diagnosis, y=count)) +
  geom_segment( aes(x=Neuropath_diagnosis, xend=Neuropath_diagnosis, y=0, yend=count), color="grey", size = 3) +
  geom_point(size=5, color="#69b3a2") +
  coord_flip() +
  theme_bw(base_size = 20) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none",
  ) +
  xlab("")
ggsave("./plots/neuropath_dx_counts.png")

clinical_dx_summary <- final %>% 
  group_by(Clinical_diagnosis) %>% 
  summarise(count = n()) %>% 
  arrange(count) %>% 
  mutate(Clinical_diagnosis=factor(Clinical_diagnosis, Clinical_diagnosis)) 

#Line plot of clinical diagnoses counts
ggplot(data = clinical_dx_summary, aes(x=Clinical_diagnosis, y=count)) +
  geom_segment( aes(x=Clinical_diagnosis, xend=Clinical_diagnosis, y=0, yend=count), color="grey", size = 3) +
  geom_point(size=5, color="blue") +
  coord_flip() +
  theme_bw(base_size = 20) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none",
  ) +
  xlab("")
ggsave("./plots/clinical_dx_counts.png")


final %>% 
  group_by(case_status, Braak_aSyn_stage) %>% 
  summarise(count = n()) %>% 
  mutate(Braak_aSyn_stage=factor(Braak_aSyn_stage, Braak_aSyn_stage)) %>% 
  mutate(case_status=factor(case_status)) %>% 
  ggplot(aes(x = case_status, y = count, fill = Braak_aSyn_stage)) +
  geom_col() + 
  theme_bw(base_size = 20) +
  scale_fill_brewer(palette = "Blues")
ggsave("./plots/braak_stage_casestatus.png")

#### Braak stage GWAS #####

#Braak stages GWAS
#Inclusion criteria: LBD cases with Braak stages between 1-6
#Analysis: ordinal GWAS with sex, age at death and PCs as covariates (we can determine number of PCs based on scree plot later)
#Future plans: we can think about including controls with LB Braak stage 0 to answer different questions


