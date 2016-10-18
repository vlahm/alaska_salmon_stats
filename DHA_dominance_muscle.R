
library(car)
library(lsmeans)
# library(dplyr)
# library(devtools)
# install_github('vlahm/manipulateR')
library(manipulateR)
library(stringr)

setwd('/home/mike/git/alaska_salmon_stats')
allfa <- read.csv('DHA_dominance_muscle.csv')

#remove FAs with < 40% representation
allfa <- cbind(allfa[,1:3], matfilter(allfa[,-(1:3)], cond='==0', thresh=.6))
#re-calculate compositions
allfa[,-(1:3)] <- t(apply(allfa[,-(1:3)], 1, FUN=function(x) x/sum(x)))
#arc-sine square root transform compositional data for MANOVA
allfa[,-(1:3)] <- asin(sqrt(allfa[,-(1:3)]))
#simplify sample IDs
allfa$Sample.ID <- str_match(allfa$Sample.ID, '(\\d+)-.*')[,2]
#separate creeks
hanson <- allfa[which(allfa$Creek == 'H'),]
hfa <- as.matrix(hanson[,-(1:3)])
pick <- allfa[which(allfa$Creek == 'P'),]
pfa <- as.matrix(pick[,-(1:3)])



mod1 <- lm(hfa ~ hanson$Location)
man1 <- Anova(mod1)
summary(man1)
rg1 <- ref.grid(mod1)
contrast(rg1, method='pairwise')

mod_skin_grp_1 <- lm(FA_skin_grouped ~ location, data=skin)
summary(mod_skin_grp_1)
manova_skin_grp_1 <- Anova(mod_skin_grp_1)
summary(manova_skin_grp_1)
rg_skin_grp_1 <-ref.grid(mod_skin_grp_1)
contrast(rg_skin_grp_1, method='pairwise')
