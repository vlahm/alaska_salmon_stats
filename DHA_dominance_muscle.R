#Jordan Lee - Alaska salmon FA stats
#contact: Mike Vlah (vlahm@uw.edu)
#last edit: 18 Oct 2016

rm(list=ls())

#setup ####

#install any necessary packages you don't already have
package_list <- c('car','lsmeans','devtools','stringr')
new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

#load packages
# install_github('vlahm/manipulateR')
library(manipulateR)
library(lsmeans)
library(devtools)
library(stringr)
library(car)
# library(DirichletReg)

setwd('/home/mike/git/alaska_salmon_stats')
setwd('C:/Users/Mike/git/j_lee_stats_2015')
allfa <- read.csv('DHA_dominance_muscle.csv')


#preprocessing ####

#remove FAs with < 40% representation
    #[NOTICE: take a look at the resulting dataframe (allfa)
    #and make sure none of the FAs that you need to include are missing. If you need
    #them all, change "thresh" to 1 in the line below.]
allfa <- cbind(allfa[,1:3], matfilter(allfa[,-(1:3)], cond='==0', thresh=.6))

#replace 0s with very small values
# allfa[,-(1:3)][allfa[,-(1:3)]==0] <- 0.00000001

#re-calculate compositions
allfa[,-(1:3)] <- t(apply(allfa[,-(1:3)], 1, FUN=function(x) x/sum(x)))

#arc-sine square root transform compositional data for MANOVA
allfa_pre_transform <- allfa
allfa[,-(1:3)] <- asin(sqrt(allfa[,-(1:3)]))

#simplify sample IDs
allfa$Sample.ID <- str_match(allfa$Sample.ID, '(\\d+)-.*')[,2]

#separate creeks
hanson <- allfa[which(allfa$Creek == 'H'),]
rownames(hanson) <- 1:nrow(hanson)
hfa <- as.matrix(hanson[,-(1:3)])
pick <- allfa[which(allfa$Creek == 'P'),]
rownames(pick) <- 1:nrow(pick)
pfa <- as.matrix(pick[,-(1:3)])


#dirichlet regression (may include later if necessary; skip for now) ####

# hfa <- DR_data(hfa, trafo=TRUE) #prepare response dataset for DirichReg
# mod <- DirichReg(hfa ~ hanson$Location)
# summary(mod)

# mod1 <- glm(hfa[,4:6] ~ hanson$Location, family=quasibinomial(link='logit'))


#hanson manova ####

mod1 <- lm(hfa ~ Location, data=hanson)
man1 <- Anova(mod1)
summary(man1) #says there are significant differences somewhere in the dataset

#pairwise comparisons
rg1 <- ref.grid(mod1)
comp <- summary(contrast(rg1, method='pairwise'))

#extract all DHA comparisons (this shows comparisons across time points
    #e.g. DHA-holding vs. EPA-postspawn. if you don't care about these,
    #see next section
dha_ind <- grep('DHA', as.vector(comp$contrast))
(dha_all <- as.vector(comp)[dha_ind,])

#extract all DHA comparisons within time points (corresponds to the results section)
dha_ind2 <- grep('(postspawn|entry|holding).*\\1', as.vector(dha_all[,1]))
(dha_timept <- as.vector(dha_all)[dha_ind2,])


#hanson within-timepoint results ####

#DHA is higher than all other FAs at all timepoints except for the following:
#entry 16:0
#entry & holding 18:1n9 (oleic acid)
#entry & holding 20:1n11
#entry & holding 22:1n11

#means by timepoint for comparison
aggregate(allfa_pre_transform[allfa$Creek=='H',4:ncol(allfa)],
          by=list(hanson$Location), FUN=mean)


#pick manova ####

mod2 <- lm(pfa ~ Location, data=pick)
man2 <- Anova(mod2)
summary(man2) #there are differences in this set too

#pairwise comparisons
rg2 <- ref.grid(mod2)
comp2 <- summary(contrast(rg2, method='pairwise'))

#extract all DHA comparisons (this shows comparisons across time points
    #e.g. DHA-holding vs. EPA-postspawn. if you don't care about these,
    #see next section
dha_ind3 <- grep('DHA', as.vector(comp2$contrast))
(dha_all2 <- as.vector(comp2)[dha_ind3,])

#extract all DHA comparisons within time points
dha_ind4 <- grep('(postspawn|entry|holding).*\\1', as.vector(dha_all2[,1]))
(dha_timept2 <- as.vector(dha_all2)[dha_ind4,])


#pick results ####

#DHA is higher than all other FAs at all timepoints except for the following:
#entry 16:0
#entry & holding 18:1n9 (oleic acid)
#holding 20:1n11
#holding 22:1n11 (entry is very close: p = 0.067)

#means by timepoint for comparison
aggregate(allfa_pre_transform[allfa$Creek=='P',4:ncol(allfa)],
          by=list(pick$Location), FUN=mean)
