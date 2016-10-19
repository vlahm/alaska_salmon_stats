
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


#pick manova ####

mod1 <- lm(hfa ~ Location, data=hanson)
man1 <- Anova(mod1)
summary(man1)

#pairwise comparisons
rg1 <- ref.grid(mod1)
comp <- summary(contrast(rg1, method='pairwise'))

#extract all DHA comparisons (this shows comparisons across time points
    #e.g. DHA-holding vs. EPA-postspawn. if you don't care about this,
    #see next section
dha_ind <- grep('DHA', as.vector(comp$contrast))
(dha_all <- as.vector(comp)[dha_ind,])

#extract all DHA comparisons within time points
grep('[(postspawn)(entry)]', as.vector(dha_all[,1]))
