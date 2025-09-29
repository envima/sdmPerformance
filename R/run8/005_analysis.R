#'@name analysis.R
#'@date 29.09.2025
#'@author Lisa Bald [bald@staff.uni-marburg.de]
#'@description plots the results

# 1 - install and load packages  ####
#-----------------------------------#

# install.packages("ContDataQC")
#remotes::install_github("leppott/ContDataQC", force = TRUE, build_vignettes = F)

library(ContDataQC)
library(ggplot2)
library(tidyverse)
library(hrbrthemes)
library(hexbin)
library(egg)
library(grid)


if (Sys.info()[[4]]=="PC19674") {
  setwd("M:/user/bald/SDM/sdmPerformance")
  nCores=1
} else if (Sys.info()[[4]]=="pc19543") {
  nCores=60
}

# create a unique name for different runs
nameRun <- paste0("run", 8)

source(paste0("R/",nameRun,"/functions/scale_metric.R"))

# 2 - load data ####
#------------------#

data=readRDS(paste0("data/",nameRun,"/results.RDS"))
data$points <-factor(data$points, levels = c(40, 80, 120, 160, 200, 400))

# List of metrics to scale
metric_names <- c("AUC", "COR", "Spec", "Sens", "Kappa", "SBI_m",
                  "PCC", "TSS", "PRG", "SEDI", "trueCor")

# Apply scaling function to each metric and create new "_scaled" columns
data <- data %>%
  mutate(across(all_of(metric_names),
                .fns = ~ scale_metric(., cur_column()),
                .names = "{.col}_scaled"))

if(!dir.exists(paste0("images/",nameRun,"/resultPlots"))) dir.create(paste0("images/",nameRun,"/resultPlots"), recursive=T)


set.seed(123)

dataClean <- data %>%
  group_by(species, model, size, testData, replicate, points) %>%
  filter(n() == 3) %>%
  ungroup()


# 3 - t test ####
#---------------#


metric_names <- c( "AUC_scaled", "COR_scaled", "PRG_scaled", "Spec_scaled", "Sens_scaled", "TSS_scaled", "Kappa_scaled", 
                   "PCC_scaled")#"SEDI_scaled", "SBI_m_scaled"

all_results=list()

for (m in metric_names) {
  print(m)
  
  PBG <- dataClean%>%dplyr::filter(method == "PBG")%>%dplyr::pull(m)
  PAA   <- dataClean%>%dplyr::filter(method == "PAA")%>%dplyr::pull(m)
  PA   <- dataClean%>%dplyr::filter(method == "PA")%>%dplyr::pull(m)
  
  
  # Run one-sided t-test (H0: corrected <= uncorrected; HA: corrected > uncorrected)
  #test_result <- rquery.t.test(
  #  x = PA,
  #  y = PAA,
  #  # alternative = "greater",   # one-sided test+
  #  alternative = "two.sided",  # two-sided test
  #  paired = TRUE              # paired observations
  #)
  
  # Inspect results
  #print(test_result)
  
  diff_PAA <- abs(PA - PAA)
  diff_PBG <- abs(PA - PBG)
  
  # H₀: The median of (diff_PAA−diff_PBG)=0. (p.value less than 0.05)
  # → No difference in closeness of PAA and PBG to PA.
  # H₁: The median of (diff_PAA−diff_PBG)≠0.
  # → One method is closer to PA than the other.
  test_result_twoSided <- wilcox.test(diff_PAA, diff_PBG, paired = TRUE, alternative = "two.sided")
  
  # H₀: The median of (diff_PAA−diff_PBG)≤0. (p.value less than 0.05)
  # → PAA is equally close to or closer than PBG to PA.
  # H₁: The median of (diff_PAA−diff_PBG)>0.
  # → PAA is farther from PA than PBG.
  test_result_greater <- wilcox.test(diff_PAA, diff_PBG, paired = TRUE, alternative = "g")
  
  #H₀: The median of (diff_PAA−diff_PBG)≥0. (p.value less than 0.05)
  #→ PAA is equally far or farther from PA than PBG.
  #H₁: The median of (diff_PAA−diff_PBG)<0.
  #→ PAA is closer to PA than PBG.
  test_result_lesser <- wilcox.test(diff_PAA, diff_PBG, paired = TRUE, alternative = "l")
  
  df=data.frame(metric=m,
                twoSided=round(test_result_twoSided$p.value,2),
                greater=round(test_result_greater$p.value,2),
                lesser = round(test_result_lesser$p.value,2))
  all_results[[m]] <- df
  rm(PAA,PA, PBG,df,test_result_twoSided ,  test_result_greater,  test_result_lesser, diff_PAA, diff_PBG)
};rm(data,scale_metric)


results=do.call(rbind, all_results)
rownames(results) <- 1:nrow(results)
saveRDS(results, paste0("data/",nameRun,"/results_wilcoxTest.RDS"))
write.csv(results, paste0("data/",nameRun,"/results_wilcoxTest.csv"))
rm(all_results,results)


# 5 prevalence ####
#-----------------#


prev=data.frame(species= c("VS01", "VS02", "VS03", "VS04", "VS05", "VS06", "VS07", "VS08", "VS09", "VS10"),
                prevalence= c(0.35,0.34,0.33,0.29,0.26,0.21,0.15,0.12,0.11,0.05),
                class=c("Broad","Broad","Broad","Middle","Middle","Middle","Narrow","Narrow","Narrow","Narrow"))


#VS1–3 were calibrated to be broadly distributed, VS7–10 narrowly distributed and VS4–6 in between
dataClean=merge(dataClean,prev, by.x="species",by.y="species")

all_results=list()

for (m in metric_names) {
  print(m)
  for (prevalence_method in c("Broad", "Narrow", "Middle")){
    
    PBG <- dataClean%>%dplyr::filter(method == "PBG")%>%dplyr::filter(class == prevalence_method)%>%dplyr::pull(m)
    PAA   <- dataClean%>%dplyr::filter(method == "PAA")%>%dplyr::filter(class == prevalence_method)%>%dplyr::pull(m)
    PA   <- dataClean%>%dplyr::filter(method == "PA")%>%dplyr::filter(class == prevalence_method)%>%dplyr::pull(m)
    
    
    diff_PAA <- abs(PA - PAA)
    diff_PBG <- abs(PA - PBG)
    
    # H₀: The median of (diff_PAA−diff_PBG)=0. (p.value less than 0.05)
    # → No difference in closeness of PAA and PBG to PA.
    # H₁: The median of (diff_PAA−diff_PBG)≠0.
    # → One method is closer to PA than the other.
    test_result_twoSided <- wilcox.test(diff_PAA, diff_PBG, paired = TRUE, alternative = "two.sided")
    
    # H₀: The median of (diff_PAA−diff_PBG)≤0. (p.value less than 0.05)
    # → PAA is equally close to or closer than PBG to PA.
    # H₁: The median of (diff_PAA−diff_PBG)>0.
    # → PAA is farther from PA than PBG.
    test_result_greater <- wilcox.test(diff_PAA, diff_PBG, paired = TRUE, alternative = "g")
    
    #H₀: The median of (diff_PAA−diff_PBG)≥0. (p.value less than 0.05)
    #→ PAA is equally far or farther from PA than PBG.
    #H₁: The median of (diff_PAA−diff_PBG)<0.
    #→ PAA is closer to PA than PBG.
    test_result_lesser <- wilcox.test(diff_PAA, diff_PBG, paired = TRUE, alternative = "l")
    
    df=data.frame(metric=m,
                  prevalence=prevalence_method,
                  twoSided=round(test_result_twoSided$p.value,2),
                  greater=round(test_result_greater$p.value,2),
                  lesser = round(test_result_lesser$p.value,2))
    all_results[[paste0(m,"_",prevalence_method)]] <- df
    rm(PAA,PA, PBG,df,test_result_twoSided ,  test_result_greater,  test_result_lesser, diff_PAA, diff_PBG)
  }
};rm(prev)



resultsPrevalence=do.call(rbind, all_results)
rownames(resultsPrevalence) <- 1:nrow(resultsPrevalence)
saveRDS(resultsPrevalence, paste0("data/",nameRun,"/results_wilcoxTest_prevalence.RDS"))
write.csv(resultsPrevalence, paste0("data/",nameRun,"/results_wilcoxTest_prevalence.csv"))
rm(all_results,dataClean,resultsPrevalence)
