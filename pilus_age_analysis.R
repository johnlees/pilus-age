library(dplyr)
library(tidyr)
library(lme4)

# rrgB COGs are CLS02709, CLS03842 and CLS01960

#
# Maela
#

Drug_resistant <- read.delim("Drug_resistant.txt", stringsAsFactors=FALSE)
epidata <- read.delim("epidata.txt", header=FALSE, stringsAsFactors=FALSE)
colnames(epidata) = c("sample", "host", "date", "longitude", "latitude")


maela <- read.delim("maela.input", stringsAsFactors=FALSE)
maela <- dplyr::as_data_frame(maela[,c("Taxon", "Time", "Serotype", "SC", "CLS02709", "CLS03842", "CLS01960")])
any_pilB = maela$CLS02709 | maela$CLS03842 | maela$CLS01960
lanes = Drug_resistant[match(maela$Taxon, Drug_resistant$SMRU_id), "lane_id"]

# For R
maela_out_R <- maela %>%
                 mutate(lane = lanes, pilus = if_else(any_pilB, 1, 0)) %>%
                 select(lane, pilus, CLS01960, CLS02709, CLS03842)

# For pyseer
maela_out_pyseer <- t(maela %>%
  mutate(Gene = lanes, pilus = if_else(any_pilB, 1, 0)) %>%
  select(Gene, pilus, CLS01960, CLS02709, CLS03842))
write.table(maela_out_pyseer, file="pilus.Rtab", quote = F, sep = "\t", col.names = F)

mother_pheno = data.frame(samples=epidata$sample, mother=ifelse(unlist(strsplit(epidata$host, "_"))[rep(c(FALSE,TRUE))] == "MOTHER", 1, 0))
write.table(mother_pheno, "mother.pheno", quote = F, row.names = F, sep = "\t")

covariates <- read.delim("covariates.txt", header=FALSE, stringsAsFactors=FALSE)
mat_ab_exclude = covariates[exp(covariates[,3]) < 183, 1]
pheno_no_infant = mother_pheno[!(mother_pheno$samples %in% mat_ab_exclude),]
write.table(pheno_no_infant, "mother_no_infant.pheno", quote = F, row.names = F, sep = "\t")

# removing <6 months
maela_df = merge(pheno_no_infant, maela_out_R, by.x = "samples", by.y = "lane")

# Simple association
maela_lm = glm(mother ~ pilus, maela_df, family = binomial())
summary(maela_lm)$coefficients

#              Estimate Std. Error    z value     Pr(>|z|)
# (Intercept) -0.8819167 0.05681796 -15.521794 2.470478e-54
# pilus       -1.1422262 0.13250811  -8.620047 6.692653e-18

# GWAS with pyseer (run in terminal)

# pyseer --lmm --similarity nj_tree_dists.txt --phenotypes mother_no_infant.pheno --pres pilus.Rtab
# Read 2600 phenotypes
# Detected binary phenotype
# Setting up LMM
# Similarity matrix has dimension (3039, 3039)
# Analysing 2497 samples found in both phenotype and similarity matrix
# h^2 = 0.65
# variant	af	filter-pvalue	lrt-pvalue	beta	beta-std-err	variant_h2	notes
# pilus	2.71E-01	8.66E-16	4.98E-01	-2.76E-02	4.07E-02	1.36E-02
# CLS01960	5.89E-02	2.99E-07	4.71E-01	-4.14E-02	5.73E-02	1.44E-02
# CLS02709	1.94E-01	4.09E-08	8.37E-01	1.31E-02	6.40E-02	4.11E-03
# CLS03842	1.80E-02	1.34E-01	3.63E-01	-1.14E-01	1.25E-01	1.82E-02

# LMM with host
carriage_eps <- read.delim("carriage_eps_longest_unique.txt", stringsAsFactors=TRUE, header = F)
colnames(carriage_eps) = c("carriage_episode", "samples")
maela_df_lmm = merge(maela_df, carriage_eps, all.x = T)

# Add single ep ids to missing (mostly mothers, who have short carriage)
single_carriage_ids = seq(from = max(maela_df_lmm$carriage_episode, na.rm = T) + 1,
    to = max(maela_df_lmm$carriage_episode, na.rm = T) + sum(is.na(maela_df_lmm$carriage_episode)))
maela_df_lmm[is.na(maela_df_lmm$carriage_episode), "carriage_episode"] = single_carriage_ids
maela_df_lmm$carriage_episode = as.factor(maela_df_lmm$carriage_episode)

maela_glmm = glmer(formula = mother ~ pilus + (1 | carriage_episode), 
                 data = maela_df_lmm, verbose = 2,
                 family = binomial(link = "logit"), nAGQ = 25)
summary(maela_glmm)

maela_glmm_nopil = glmer(formula = mother ~ (1 | carriage_episode), 
                         data = maela_df_lmm, verbose = 2,
                         family = binomial(link = "logit"), nAGQ = 25)
anova(maela_glmm, maela_glmm_nopil, test="LRT")

# Simple frequency table
table(maela_df_lmm$mother, maela_df_lmm$pilus)

