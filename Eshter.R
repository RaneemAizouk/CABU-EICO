#################################################################################
# SAMPLE SIZE CALCULATION FOR WATCH ABX JPI AMR
#################################################################################

# Author E van Kleef
# Using methodology:
# 1) Using clusterPower package
# 2) Using EpiDisplay package
# 3) As described in Arnold et al 2011 BMC Medical Research Methodology https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-11-94


rm(list=ls()) 
install.packages("Hmisc")
install.packages("ggplot2")
install.packages("clusterPower")

set.seed(123)
library(Hmisc); library(rms); library(epiDisplay); library("tidyverse"); library("geepack"); library(readxl); library(writexl);
library("ggplot2"); library("clusterPower")

# Read in parameterset for powercalculations
dpower_cp = read.csv("./Data/parameters_clusterPower.csv", sep=";", colClasses = "numeric")
names(dpower_cp)[1] = "nclusters"

#######################################################################################
# 1) clusterPower Package
#######################################################################################

# Difference in Difference approach (https://cran.rstudio.com/web/packages/clusterPower/vignettes/clusterpower.html
# i.e. net difference, see also https://researchmethodsresources.nih.gov/grt-calculator#/Step1)
#######################################################################################

# Simulation approach (takes long time to run, but provides similar results to the above). sigma_b_sq0 = 0.3 ~ icc = 0.027
# The below parameters provide a power of 0.85
did.binary.sim = cps.did.binary(nsim = 100, 
                                nsubjects = 100, 
                                nclusters = 11,
                                p1t0 = 0.28, 
                                p2t0 = 0.28,
                                p1t1 = 0.28, 
                                p2t1 = 0.28*-0.6+0.28, 
                                sigma_b_sq0 = 0.3,
                                sigma_b_sq1 = 0.3, 
                                alpha = 0.05,
                                method ='gee', allSimData = FALSE, seed = 123)

did.binary.sim

# Compare with analytical approach
# The below parameters provides a power of 
base = 0.28
effect = 0.6
p = (base+base*(1-effect))/2
d = base - base*(1-effect) # similar to (C_post - C_base) - (I_post - I_base), see also https://researchmethodsresources.nih.gov/grt-calculator#/Step1

cpa.did.binary(alpha = 0.05, 
               power = NA, 
               nclusters = 10, 
               nsubjects = 100, 
               p = p, 
               d = d, 
               ICC = 0.027, 
               rho_c = 0, 
               rho_s=0)


# Simulation approach - all paramater values
dpower_cp_diff_sim = dpower_cp %>%
  mutate(power = NA,
         sd = case_when(
           icc == 0.01 ~ 0.1,
           icc == 0.03 ~ 0.3,
           icc == 0.05 ~ 0.6),
         icc = NA
  ) %>%
  filter(nsubjects==100) # To limit computational time only take nsubjects = 100

for(i in 1:nrow(dpower_cp_diff_sim)){
  did.binary.sim = cps.did.binary(nsim = 100, 
                                  nsubjects = dpower_cp_diff_sim$nsubjects[i], 
                                  nclusters = dpower_cp_diff_sim$nclusters[i],
                                  p1t0 = dpower_cp_diff_sim$base_prev[i], 
                                  p2t0 = dpower_cp_diff_sim$base_prev[i],
                                  p1t1 = dpower_cp_diff_sim$base_prev[i], 
                                  p2t1 = dpower_cp_diff_sim$base_prev[i]*-dpower_cp_diff_sim$effect[i]+dpower_cp_diff_sim$base_prev[i], 
                                  sigma_b_sq0 = dpower_cp_diff_sim$sd[i],
                                  sigma_b_sq1 = dpower_cp_diff_sim$sd[i], 
                                  alpha = 0.05,
                                  method ='gee', allSimData = FALSE, seed = 123,
                                  lowPowerOverride = TRUE)
  dpower_cp_diff_sim$icc[i] = did.binary.sim$ICC[4]
  dpower_cp_diff_sim$power[i] = did.binary.sim$power[1]
}
dpower_cp_diff_sim$prev_lab = paste0("Base p=",dpower_cp_diff_sim$base_prev)
dpower_cp_diff_sim$effect_lab = paste0("Reduc=",dpower_cp_diff_sim$effect*100,"%")
dpower_cp_diff_sim$power = unlist(dpower_cp_diff_sim$power)


# Analytical approach - all parameter values
dpower_cp_diff = dpower_cp
dpower_cp_diff$power = NULL

for(i in 1:nrow(dpower_cp_diff)){
  base = dpower_cp_diff$base_prev[i]
  effect = dpower_cp_diff$effect[i]
  p = (base+base*(1-effect))/2
  d = base - base*(1-effect)
  power = cpa.did.binary(alpha = 0.05, 
                         power = NA, 
                         nclusters = dpower_cp_diff$nclusters[i], 
                         nsubjects = dpower_cp_diff$nsubjects[i], 
                         p = p, 
                         d = d, 
                         ICC = dpower_cp_diff$icc[i], 
                         rho_c = 0, 
                         rho_s=0)
  dpower_cp_diff$power[i] = power
}
dpower_cp_diff$prev_lab = paste0("Base p=",dpower_cp_diff$base_prev)
dpower_cp_diff$effect_lab = paste0("Reduc=",dpower_cp_diff$effect*100,"%")

# Simple two-arm comparison designs (https://cran.rstudio.com/web/packages/clusterPower/vignettes/clusterpower.html)
# i.e. # Simple difference approach (see also https://researchmethodsresources.nih.gov/grt-calculator#/Step1)
# This approach requires less power
#######################################################################################

# Simulation approach
cps.binary(
  nsim = 100,
  nsubjects = 100,
  nclusters = 11,
  p1 = 0.28,
  p2 = 0.28*-0.6+0.28,
  sigma_b_sq = 0.3,
  sigma_b_sq2 = 0.3,
  alpha = 0.05,
  method = "gee",allSimData = FALSE, seed = 123,
  lowPowerOverride = TRUE
  
)

# Analytical approach - Illustration to compare with difference in difference approach
cpa.binary(alpha = 0.05,
           power = NA,
           nclusters = 11,
           nsubjects = 100,
           CV = 0,
           p1 = 0.28,
           p2 = 0.28*-0.4+0.28,
           ICC = 0.03,
           pooled = FALSE,
           p1inc = TRUE,
           tdist = TRUE,tol = .Machine$double.eps^0.25)

# Simulation approach for all parametersets
dpower_cp_simp_sim = dpower_cp %>%
  mutate(power = NA,
         sd = case_when(
           icc == 0.01 ~ 0.1,
           icc == 0.03 ~ 0.3,
           icc == 0.05 ~ 0.6),
         icc = NA
  ) %>%
  filter(nsubjects==100) # To limit computational time only take nsubjects = 100

for(i in 1:nrow(dpower_cp_simp_sim)){
  cps.binary.sim = cps.binary(
    nsim = 100,
    nsubjects = dpower_cp_simp_sim$nsubjects[i],
    nclusters = dpower_cp_simp_sim$nclusters[i],
    p1 = dpower_cp_simp_sim$base_prev[i],
    p2 = dpower_cp_simp_sim$base_prev[i]*-dpower_cp_simp_sim$effect[i]+dpower_cp_simp_sim$base_prev[i],
    sigma_b_sq = dpower_cp_simp_sim$sd[i],
    sigma_b_sq2 = dpower_cp_simp_sim$sd[i],
    alpha = 0.05,
    method = "gee",allSimData = FALSE, seed = 123,
    lowPowerOverride = TRUE
  )
  dpower_cp_simp_sim$icc[i] = cps.binary.sim$ICC[3]
  dpower_cp_simp_sim$power[i] = cps.binary.sim$power[1]
}
dpower_cp_simp_sim$prev_lab = paste0("Base p=",dpower_cp_simp_sim$base_prev)
dpower_cp_simp_sim$effect_lab = paste0("Reduc=",dpower_cp_simp_sim$effect*100,"%")
dpower_cp_simp_sim$power = unlist(dpower_cp_simp_sim$power)


# Analytical approach for all parametersets
dpower_cp_simp = dpower_cp
dpower_cp_simp$power = NULL

for(i in 1:nrow(dpower_cp_simp)){
  power =  cpa.binary(alpha = 0.05,
                      power = NA,
                      nclusters = dpower_cp_simp$nclusters[i],
                      nsubjects = dpower_cp_simp$nsubjects[i],
                      CV = 0,
                      p1 = dpower_cp_simp$base_prev[i],
                      p2 = dpower_cp_simp$base_prev[i]*-dpower_cp_simp$effect[i]+dpower_cp_simp$base_prev[i],
                      ICC = dpower_cp_simp$icc[i],
                      pooled = FALSE,
                      p1inc = TRUE,
                      tdist = TRUE,tol = .Machine$double.eps^0.25)
  dpower_cp_simp$power[i] = power 
}
dpower_cp_simp$prev_lab = paste0("Base p=",dpower_cp_simp$base_prev)
dpower_cp_simp$effect_lab = paste0("Reduc=",dpower_cp_simp$effect*100,"%")

# Plot results
##############################################################################################
unique(round(dpower_cp_diff_sim$icc,2))

# Difference in difference - simulation approach
p0 = ggplot(dpower_cp_diff_sim, aes(x=effect*100,y=power,col = factor(sd)), group=factor(sd)) + 
  geom_point(size=4) + 
  theme_bw() +
  geom_line(size=1) + 
  xlab("Trial effect (%)") + 
  theme(axis.text=element_text(size=12) ) +
  facet_wrap(~prev_lab, ncol=5)+
  scale_color_manual(name = "SD",
                     values = c("0.1" = "brown2", "0.3" = "darkgoldenrod2", "0.6" = "cyan4"),
                     labels = c("0.1" = paste("0.1 (ICC ~ 0.01)"), "0.3" = paste("0.3 (ICC ~ 0.03)"),"0.6" = paste("0.6 (ICC ~ 0.05)"))) +
  geom_hline(yintercept = 0.8, linetype=2) +
  labs(title = "Power estimation binary outcome - net difference (simulation)",
       subtitle = "Cluster size per arm = 11; nsubjects per cluster = 100")

p0

# Difference in difference - Analytical approach
p1 = ggplot(dpower_cp_diff, aes(x=nsubjects,y=power,col = factor(icc)), group=factor(icc)) + 
  geom_point(size=4) + 
  theme_bw() +
  geom_line(size=1) + 
  xlab("Subjects per cluster") + 
  theme(axis.text=element_text(size=12) ) +
  facet_wrap(~prev_lab+effect_lab, ncol=5)+
  scale_color_manual(name = "ICC",
                     values = c("0.01" = "brown2", "0.03" = "darkgoldenrod2", "0.05" = "cyan4")) +
  geom_hline(yintercept = 0.8, linetype=2) +
  labs(title = "Power estimation binary outcome - net difference (analytical)",
       subtitle = "Cluster size per arm = 11")
p1

# Difference in difference - Analytical approach just 100 participants
p1.1 = dpower_cp_diff %>% filter(nsubjects==100) %>%
  ggplot(., aes(x=effect*100,y=power,col = factor(icc)), group=factor(icc)) + 
  geom_point(size=4) + 
  theme_bw() +
  geom_line(size=1) + 
  xlab("Trial effect (%)") + 
  theme(axis.text=element_text(size=12) ) +
  facet_wrap(~prev_lab, ncol=5)+
  scale_color_manual(name = "SD",
                     values = c("0.01" = "brown2", "0.03" = "darkgoldenrod2", "0.05" = "cyan4")) +
  geom_hline(yintercept = 0.8, linetype=2) +
  labs(title = "Power estimation binary outcome - net difference (analytical)",
       subtitle = "Cluster size per arm = 11; nsubjects per cluster = 100")


# Simple difference - simulation approach
p2 = ggplot(dpower_cp_simp_sim, aes(x=effect*100,y=power,col = factor(sd)), group=factor(sd)) + 
  geom_point(size=4) + 
  theme_bw() +
  geom_line(size=1) + 
  xlab("Trial effect (%)") + 
  theme(axis.text=element_text(size=12) ) +
  facet_wrap(~prev_lab, ncol=5)+
  scale_color_manual(name = "SD",
                     values = c("0.1" = "brown2", "0.3" = "darkgoldenrod2", "0.6" = "cyan4"),
                     labels = c("0.1" = paste("0.1 (ICC ~ 0.01)"), "0.3" = paste("0.3 (ICC ~ 0.03)"),"0.6" = paste("0.6 (ICC ~ 0.05)"))) +
  geom_hline(yintercept = 0.8, linetype=2) +
  labs(title = "Power estimation binary outcome - simple difference (simulation)",
       subtitle = "Cluster size per arm = 11; nsubjects per cluster = 100")

p2


# Simple difference - analytical approach
p3 = ggplot(dpower_cp_simp, aes(x=nsubjects,y=power,col = factor(icc)), group=factor(icc)) + 
  geom_point(size=4) + 
  theme_bw() +
  geom_line(size=1) + 
  xlab("Subjects per cluster") + 
  theme(axis.text=element_text(size=12) ) +
  facet_wrap(~prev_lab+effect_lab, ncol=5)+
  scale_color_manual(name = "ICC",
                     values = c("0.01" = "brown2", "0.03" = "darkgoldenrod2", "0.05" = "cyan4")) +
  geom_hline(yintercept = 0.8, linetype=2) +
  labs(title = "Power estimation binary outcome - simple difference (analytical)",
       subtitle = "Cluster size per arm = 11")
p3

# Simple difference - Analytical approach just 100 participants
p3.1 = dpower_cp_simp %>% filter(nsubjects==100) %>%
  ggplot(., aes(x=effect*100,y=power,col = factor(icc)), group=factor(icc)) + 
  geom_point(size=4) + 
  theme_bw() +
  geom_line(size=1) + 
  xlab("Trial effect (%)") + 
  theme(axis.text=element_text(size=12) ) +
  facet_wrap(~prev_lab, ncol=5)+
  scale_color_manual(name = "SD",
                     values = c("0.01" = "brown2", "0.03" = "darkgoldenrod2", "0.05" = "cyan4")) +
  geom_hline(yintercept = 0.8, linetype=2) +
  labs(title = "Power estimation binary outcome - simple difference (analytical)",
       subtitle = "Cluster size per arm = 11; nsubjects per cluster = 100")


# Save plots
pdf("./Output/Figures/Net_diff_vs_simp_clusterPower_sim.pdf", width=10, height=5)
print(p0)
print(p2)
dev.off()

pdf("./Output/Figures/Net_diff_vs_simp_clusterPower_anal.pdf", width=10, height=5)
print(p1.1)
print(p3.1)
dev.off()

pdf("./Output/Figures/Net_diff_clusterPower_nsubjects_anal.pdf", width=12, height=7)
print(p1)
dev.off()


pdf("./Output/Figures/Simple_diff_clusterPower_nsubjects_anal.pdf", width=12, height=7)
print(p3)
dev.off()

# Write results to Table
write_xlsx(dpower_cp_diff,"./Output/Tables/Net_diff_clusterPower_anal.xlsx")
write_xlsx(dpower_cp_simp,"./Output/Tables/Simp_diff_clusterPower_anal.xlsx")
write_xlsx(dpower_cp_diff_sim,"./Output/Tables/Net_diff_clusterPower_sim.xlsx")
write_xlsx(dpower_cp_simp_sim,"./Output/Tables/Simp_diff_clusterPower_sim.xlsx")

# SAMPSIZE PACKAGE
#######################################################################################

n.for.cluster.2p(p1=.28, p2=0.28*-0.4+0.28, alpha = 0.05, power = 0.80, ratio=1, 
                 mean.cluster.size = 100, max.cluster.size = 100, 
                 min.cluster.size = 100, icc = 0.05)


# Arnold et al 2011 
#################################################################################
#Power simulation for Yij ~ [1 + exp(-[mu + beta*Ai + bi])]^-1
#  (binary outcome  with cluster (i) variability)

#Function parameters:
# mu: mean prevalence of outcome in the control group
# beta: odds ratio of outcome in treatment vs. control
# sdClust: standard deviation of random effect at cluster level
# nPartPerClust: number of interviews per cluster
# nTreatClust: number of treatment clusters
# nControlClust: number of control (untreated) clusters
# nIterations: number of iterations for simulation

#Function for calculating study power for given design and treatment effects

# fnPower <- function(mu,beta,sdClust,nPartPerClust,nTreatClust,nControlClust,nIterations,dots=TRUE){ 
#   start.time <- Sys.time() 
#   
#   if(dots) cat("Simulations (",nIterations,") \n----|--- 1 ---|--- 2 ---|--- 3 ---|--- 4 ---| --- 5 \n",sep="") 
#   
#   #objects to store pvalue, beta, and standard error from each iteration of simulation  
#   pVec <- betaVec <- seVec <- rep(NA,nIterations) 
#   
#   #build design matrix 
#   nClust <- nTreatClust + nControlClust 
#   nObs <- nPartPerClust * nClust 
#   m <- matrix(NA,nrow=nObs,ncol=3) 
#   colnames(m) <- c("cluster","participant","tx") 
#   m[,1] <- rep(1:nClust,each=nPartPerClust) 
#   m[,2] <- 1:nrow(m) 
#   m[,3] <- c(rep(1,nTreatClust*nPartPerClust),rep(0,nTreatClust*nPartPerClust)) 
#   
#   tx <- m[,3]  #treatment dummy  
#   
#   for(i in 1:nIterations){
#     
#     #draw random effects for clusters
#     clustRandEffect <- rep(rnorm(nClust,0,sdClust),each=nPartPerClust)
#     
#     #create outcome
#     prob <- 1/(1 + exp(-(log(mu/(1-mu)) + log(beta)*tx + clustRandEffect)))
#     y <- rbinom(nObs,1,prob)
#     
#     #fit model, store p-value, beta and standard error
#     o <- try(lrm(y~tx,x=TRUE,y=TRUE))
#     #lrm() function occasionally gets stuck in null scenario
#     # (if it does, change initial parameter value as work around)
#     #Note: We have contacted author of rms package, who is correcting this error 
#     # (i.e., we will be able to remove this line before final publication)
#     if(inherits(o,"try-error")==TRUE){ 
#       o <- lrm(y~tx,x=TRUE,y=TRUE,initial=0.5)   
#       }
#     
#     v <- robcov(fit=o,cluster=m[,1]) 
#     pVec[i] <- 2*pnorm(-abs(v$coeff[2]/sqrt(diag(v$var))[2]))
#     betaVec[i] <- v$coef[2]
#     seVec[i] <- sqrt(diag(v$var))[2]
#     if(dots) cat(".",sep="")
#     if(dots && i %% 50 == 0) cat(i,"\n")
# }   
#   if(dots) cat("\nSimulation Run Time:",round(difftime(Sys.time(),start.time,units="hours"),3)," Hours \n")   
#   
#   #calculate power   
#   power <- length(pVec[pVec<0.05])/length(pVec)   
#   
#   return(list(power=power,p=pVec,beta=betaVec,se=seVec)) 
# }
# 
# # ICC = sigma_g^2/(sigma_g^2+sigma_e^2); for log linear model, we only use sdClust
# # sigma_g = between cluster variability = sd(bi), where bi ~ N(0, sigma_g) = sdClust
# # sigma_e = within cluster variability = sd(eij), where eij ~ N(0, sigma_e) = sdResid
# # 
# # For linear model, the below would provide estimate for ICC
# sdClust=0.3; sdResid=1.297  # sdClust=0.3; sdResid=1.297 = ICC of 0.05
# sdClust^2/(sdClust^2+sdResid^2) 
# 
# # # Log odds to probability function
# # logit2prob <- function(logit){
# #   odds <- exp(logit)
# #   prob <- odds / (1 + odds)
# #   return(prob)
# # }
# # logit2prob(-1.1) # = 0.25 probability of WATCH prescribing; want to reduce with 40% = 0.15 probability
# # logit2prob(-1.73) # = 0.15 probability
# # 
# mu = 0.25
# beta = 0.6 # = a Odds ratio of 0.6 (not entirely 40% reduction)
# nPartPerClust=100
# nTreatClust=10
# nControlClust=10
# nIterations = 1000
# 
# # Can take 200 individuals per cluster, as we have two sampling moments.
# #Calibrate by setting beta = 0
# outCalibrate <- fnPower(mu=mu,beta=1,sdClust=sdClust,nPartPerClust=nPartPerClust,nTreatClust=nTreatClust,nControlClust=nControlClust,
#                         nIterations=nIterations)
# outCalibrate$power
# hist(outCalibrate$beta)
# round(quantile(outCalibrate$beta, probs = c(0.25, 0.5, 0.75)), 4)
# 
# hist(outCalibrate$p)
# 
# #Generate power power curve for different number of clusters per arm 
# # (here, assume same number of clusters in treatment and control arms)
# sdClustvec = seq(0.3,0.9,by=0.3)
# betavec = seq(0.6,0.8, by=0.1)
# nPartPerClustvec = seq(50,100,by=25)
# clustersPerArm <- seq(9,13,by=2)
# 
# drows = length(betavec)*length(clustersPerArm)*length(nPartPerClustvec)
# dpower <- as.data.frame(matrix(NA, nrow=drows,ncol=5))
# names(dpower) = c("clusters", "nPartPerClust","sdClust","beta","power")
# 
# dvecpartcl = rep(c(rep(nPartPerClustvec[1],length(clustersPerArm)),rep(nPartPerClustvec[2],length(clustersPerArm)),rep(nPartPerClustvec[3],length(clustersPerArm))),length(nPartPerClustvec))
# dvecsdcl = rep(c(rep(sdClustvec[1],length(clustersPerArm)),rep(sdClustvec[2],length(clustersPerArm)),rep(sdClustvec[3],length(clustersPerArm))),length(nPartPerClustvec))
# 
# dpower = dpower %>%
#   mutate(clusters= rep(clustersPerArm,nrow(dpower)/length(clustersPerArm)),
#          nPartPerClust  = 100,
#        #  sdClust = rep(sdClustvec,nrow(dpower)/length(clustersPerArm)),
#       #   sdClust = sort(sdClust),
#         sdClust = dvecsdcl, 
#         clusters = sort(clusters),
#          beta = rep(betavec,nrow(dpower)/length(clustersPerArm)),
#       beta_labels = paste0("OR(beta)=",beta)
#   )
# 
# # If want to change number of participants
# dpowerp = rbind(dpower,dpower,dpower)
# dpowerp$nPartPerClust = c(rep(nPartPerClustvec[1],nrow(dpower)),rep(nPartPerClustvec[2],nrow(dpower)),rep(nPartPerClustvec[3],nrow(dpower)))
# 
# #st <- c()
# #for(i in 1:length(clustersPerArm)){ 
#   #st[[paste("clust",clustersPerArm[i])]] <- 
# #    fnPower(mu=mu,beta=beta,sdClust=sdClust,nPartPerClust=nPartPerClust,nTreatClust=clustersPerArm[i],
# #            nControlClust=clustersPerArm[i],nIterations=nIterations)
# #}
# #powerVals <- lapply(st,function(x) x[[1]])
# #plot(clustersPerArm,unlist(powerVals),type='l',ylab="Power",xlab="Clusters per arm",ylim=c(0,1), lwd=2, main=paste0("Study power\n (n cluster = i; n per cluster = ", nPartPerClust, "\n mu=",mu, ", SD=",sdClust))
# #points(clustersPerArm,unlist(powerVals),pch=19,cex=2)
# #abline(h = 0.8, lty = 2, col="red", lwd=2)
# 
# #hist(exp(st$`clust 10`$beta))
# #exp(round(quantile(st$`clust 10`$beta, probs = c(0.05, 0.5, 0.95)), 4))
# 
# #png("./ITG/Other_projects/AMR/JPI-AMR/Sample size calculation/power_vs_clustersize.png")
# #plot(clustersPerArm,unlist(powerVals),type='l',ylab="Power",xlab="Clusters per arm",ylim=c(0,1), lwd=2, main=paste0("Study power\n (n cluster = i; n per cluster = ", nPartPerClust, "\n mu=",mu, ", SD=",sdClust))
# #points(clustersPerArm,unlist(powerVals),pch=19,cex=2)
# #abline(h = 0.8, lty = 2, col="red", lwd=2)
# #dev.off()
# 
# for(i in 1:nrow(dpower)){ 
#   print(paste0("Parameterset_", i))
#   p = fnPower(mu=mu,beta=dpower$beta[i],sdClust=dpower$sdClust[i],nPartPerClust=dpower$nPartPerClust[i],nTreatClust=dpower$clusters[i],
#               nControlClust=dpower$clusters[i],nIterations=nIterations)
#   dpower[i,"power"] = p$power
# }
# 
# 
# plot = ggplot(dpower, aes(x=clusters,y=power,col = factor(sdClust)), group=factor(SdClust)) + 
#   geom_point(size=4) + 
#   theme_bw() +
#   geom_line(size=1) + 
#   xlab("Clusters per arm") + 
#   theme(axis.text=element_text(size=12)
#   ) +
#   facet_wrap(~beta_labels)+
#   scale_color_manual(name = "Cluster SD",
#                      values = c("0.3" = "brown2", "0.6" = "darkgoldenrod2", "0.9" = "cyan4")) +
#   geom_hline(yintercept = 0.8, linetype=2) +
#   ggtitle("Sample size estimation for binary outcome  with cluster (i) variability \nYij ~ (1 + exp[-(mu + beta*Intervention_i + bi)])^-1")
# 
# print(plot)
# 
# pdf("./ITG/Other_projects/AMR/JPI-AMR/Sample size calculation/power_vs_clustersize.pdf", width=10,height=4)
# print(plot)
# dev.off()