#'@author https://doi.org/10.1111/ecog.07218
#'@author Canran Liu, Graeme Newell, Matt White, Josephine Machunter
#'@source https://nsojournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fecog.07218&file=ecog13216-sup-0001-AppendixS1.docx
#'@publication Improving the estimation of the Boyce index using statistical smoothing methods for evaluating species distribution models with presence-only data


#Note: Please load the R package mgcv before using this function
#	Evaluation data: include two vectors (prd1 and prd0), which are suitability values
#		predicted from the species distribution model you want to evaluate using some
#		presence points and some random points (i.e. random background sites).
#	prd1: predictions for presence sites.
#	prd0: predictions for random points.
#	ktry: basis dimension for smoothers, which is a positive integer. 
#	         Generally, the default value 10 is ok.
library(mgcv)

sfbi <- function(prd1, prd0, ktry=10) {
  
  
  
  
  
  p <- c(prd1, prd0)
  n1 <- length(prd1)
  n0 <- length(prd0)
  prd <- seq(min(p), max(p), length=n0)
  oc <- c(rep(1, n1), rep(0, n0))
  
  
  #prd_cr <- tryCatch({
  #  
  #  md_cr = mgcv::gam(oc ~ s(p,bs="cr",k=min(ktry,length(unique(p)))), family=binomial)
  #  prd_cr = predict(md_cr,newdata=data.frame(p=prd),type='response')
  #  
  #}, error = function(e) {
  #  return(NA)
  #})
  prd_cr<-NA
  
  md_tp = mgcv::gam(oc ~ s(p,bs="tp",k=min(ktry,length(unique(p)))), family=binomial)
  prd_tp = predict(md_tp,newdata=data.frame(p=prd),type='response')
  md_bs = mgcv::gam(oc ~ s(p,bs="bs",k=min(ktry,length(unique(p)))), family=binomial)
  prd_bs = predict(md_bs,newdata=data.frame(p=prd),type='response')
  md_ps = mgcv::gam(oc ~ s(p,bs="ps",k=min(ktry,length(unique(p)))), family=binomial)
  prd_ps = predict(md_ps,newdata=data.frame(p=prd),type='response')
  md_ad = mgcv::gam(oc ~ s(p, bs = "ad",k=min(ktry,length(unique(p)))), family=binomial)
  prd_ad = predict(md_ad,newdata=data.frame(p=prd),type='response')
  prd_m = (prd_tp + prd_cr + prd_bs + prd_ps + prd_ad)/5
  SBI_tp <- cor(prd,prd_tp,method="spearman")
  SBI_bs <- cor(prd,prd_bs,method="spearman")
  SBI_ps <- cor(prd,prd_ps,method="spearman")
  SBI_ad <- cor(prd,prd_ad,method="spearman")
  
  SBI_cr <- tryCatch({
    
    SBI_cr <- cor(prd,prd_cr,method="spearman")
    #return(SBI_cr)
  }
  , error = function(e) {return(NA)} 
  )
  
  
  
  SBI_m <- cor(prd,prd_m,method="spearman")
  
  return(c(SBI_tp, SBI_cr, SBI_bs, SBI_ps, SBI_ad, SBI_m))
  
  
}
