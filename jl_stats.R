#Jordan Lee - Alaska FA stats
#contact: Mike Vlah (vlahm@uw.edu)
#last edit: 16 Oct 2016

#Hi Jordan.  With the updated data, a couple of things have changed.  For muscle-grouped and skin-individual,
#the same models should still be used.  Of course, the p-values will have changed.  For the other
#two categories, the creek and sex terms are no longer significant by themselves.  Although some of the more
#complex models do show significance in these categories, the most parsimonious explanation in each case is that
#location alone is driving variation.  This probably simplifies things quite a bit for your interpretation.

setwd('/home/mike/git/alaska_salmon_stats')
library(car)
library(lsmeans)

#setup_FAgroups_muscle####
muscle <- read.csv('jl_data_muscle.csv')

#arcsin-sqrt transform proportion data for MANOVA
safa_muscle_trans <- asin(sqrt(muscle$safa_rel))
mufa_muscle_trans <- asin(sqrt(muscle$mufa_rel))
pufa_muscle_trans <- asin(sqrt(muscle$pufa_rel))

FA_muscle_grouped <- cbind(safa_muscle_trans, mufa_muscle_trans, pufa_muscle_trans)
colnames(FA_muscle_grouped) <- c('safa', 'mufa', 'pufa')

##MANOVA_FAgroups_muscle####

#I've gone through and sytematically tested the effects of location, creek, and
#sex on SAFA, MUFA, and PUFA, in all their combinations.  Where relevant, I've
#included code for pairwise comparisons.  The first series of operations is annotated.

# mod_musc_grp_1 <- lm(FA_muscle_grouped ~ location, data=muscle) #first let's see if location has an effect on SAFA_muscle, MUFA_muscle, OR PUFA_muscle relative composition.
# summary(mod_musc_grp_1) #use this to see results of univariate tests (only included this time)
# manova_musc_grp_1 <- Anova(mod_musc_grp_1) #this is the Anova function from the car package, creates manova object of class 'Anova.mlm'
# summary(manova_musc_grp_1) #Pillai, Wilks, etc. are different methods for testing significance.  This shit is mad significant across the board.
# rg_musc_grp_1 <-ref.grid(mod_musc_grp_1) #the 'reference grid' object is used by the 'lsmeans' package to perform pairwise comparisons
# contrast(rg_musc_grp_1, method='pairwise')
#
# mod_musc_grp_2 <- lm(FA_muscle_grouped ~ creek, data=muscle)
# manova_musc_grp_2 <- Anova(mod_musc_grp_2)
# summary(manova_musc_grp_2)
# rg_musc_grp_2 <-ref.grid(mod_musc_grp_2)
# contrast(rg_musc_grp_2, method='pairwise')
#
# mod_musc_grp_3 <- lm(FA_muscle_grouped ~ sex, data=muscle)
# manova_musc_grp_3 <- Anova(mod_musc_grp_3)
# summary(manova_musc_grp_3) #sex does not have an effect on FA_muscle composition
#
# mod_musc_grp_4 <- lm(FA_muscle_grouped ~ location + creek, data=muscle)
# manova_musc_grp_4 <- Anova(mod_musc_grp_4)
# summary(manova_musc_grp_4) #note that there are two separate output tables now, one for each predictor
# rg_musc_grp_4 <-ref.grid(mod_musc_grp_4)
# contrast(rg_musc_grp_4, method='pairwise')

mod_musc_grp_5 <- lm(FA_muscle_grouped ~ location + creek + location:creek, data=muscle) #interaction term
manova_musc_grp_5 <- Anova(mod_musc_grp_5)
summary(manova_musc_grp_5) #the interaction term is significant
rg_musc_grp_5 <-ref.grid(mod_musc_grp_5)
contrast(rg_musc_grp_5, method='pairwise') #whew. quite a table. enjoy.

#NOTE: I tried the 'sex' parameter with the others included (in all combinations).
#still not significant.

#setup_indivPUFAs_muscle####

#arcsin-sqrt transform proportion data for MANOVA
X18.2n6c_muscle_trans <- asin(sqrt(muscle$X18.2n6c))
X18.3n3_muscle_trans <- asin(sqrt(muscle$X18.3n3))
X20.4n3_muscle_trans <- asin(sqrt(muscle$X20.4n3))
X20.5n3_muscle_trans <- asin(sqrt(muscle$X20.5n3))
X22.5n3_muscle_trans <- asin(sqrt(muscle$X22.5n3))
X22.6n3_muscle_trans <- asin(sqrt(muscle$X22.6n3))

FA_muscle_indiv <- cbind(X18.2n6c_muscle_trans, X18.3n3_muscle_trans, X20.4n3_muscle_trans,
                         X20.5n3_muscle_trans, X22.5n3_muscle_trans, X22.6n3_muscle_trans)
colnames(FA_muscle_indiv) <- c('X18.2n6c', 'X18.3n3', 'X20.4n3', 'X20.5n3', 'X22.5n3', 'X22.6n3')

##MANOVA_indivPUFAs_muscle####

mod_musc_indiv_1 <- lm(FA_muscle_indiv ~ location, data=muscle)
summary(mod_musc_indiv_1)
manova_musc_indiv_1 <- Anova(mod_musc_indiv_1)
summary(manova_musc_indiv_1)
rg_musc_indiv_1 <-ref.grid(mod_musc_indiv_1)
contrast(rg_musc_indiv_1, method='pairwise')
#
# mod_musc_indiv_2 <- lm(FA_muscle_indiv ~ creek, data=muscle)
# manova_musc_indiv_2 <- Anova(mod_musc_indiv_2)
# summary(manova_musc_indiv_2)
# rg_musc_indiv_2 <-ref.grid(mod_musc_indiv_2)
# contrast(rg_musc_indiv_2, method='pairwise')
#
# mod_musc_indiv_3 <- lm(FA_muscle_indiv ~ sex, data=muscle)
# manova_musc_indiv_3 <- Anova(mod_musc_indiv_3)
# summary(manova_musc_indiv_3)
#
# mod_musc_indiv_4 <- lm(FA_muscle_indiv ~ location + sex, data=muscle)
# manova_musc_indiv_4 <- Anova(mod_musc_indiv_4)
# summary(manova_musc_indiv_4)
# rg_musc_indiv_4 <-ref.grid(mod_musc_indiv_4)
# contrast(rg_musc_indiv_4, method='pairwise')

# mod_musc_indiv_5 <- lm(FA_muscle_indiv ~ location + sex + location:sex, data=muscle)
# manova_musc_indiv_5 <- Anova(mod_musc_indiv_5)
# summary(manova_musc_indiv_5)
# rg_musc_indiv_5 <-ref.grid(mod_musc_indiv_5)
# contrast(rg_musc_indiv_5, method='pairwise')

########################################   NEW   ###########################
#Previously, the best model was the second-to-last (included location and sex).
#Now, the sex term is no longer significant by itself, and neither is the creek term.
#The first model is now the best one.  This makes your interpretation a bit easier.

#setup_FAgroups_skin####
skin <- read.csv('jl_data_skin.csv')

#arcsin-sqrt transform proportion data for MANOVA
safa_skin_trans <- asin(sqrt(skin$safa_rel))
mufa_skin_trans <- asin(sqrt(skin$mufa_rel))
pufa_skin_trans <- asin(sqrt(skin$pufa_rel))

FA_skin_grouped <- cbind(safa_skin_trans, mufa_skin_trans, pufa_skin_trans)
colnames(FA_skin_grouped) <- c('safa', 'mufa', 'pufa')

##MANOVA_FAgroups_skin####

mod_skin_grp_1 <- lm(FA_skin_grouped ~ location, data=skin)
summary(mod_skin_grp_1)
manova_skin_grp_1 <- Anova(mod_skin_grp_1)
summary(manova_skin_grp_1)
rg_skin_grp_1 <-ref.grid(mod_skin_grp_1)
contrast(rg_skin_grp_1, method='pairwise')

# mod_skin_grp_2 <- lm(FA_skin_grouped ~ creek, data=skin)
# manova_skin_grp_2 <- Anova(mod_skin_grp_2)
# summary(manova_skin_grp_2)
# rg_skin_grp_2 <-ref.grid(mod_skin_grp_2)
# contrast(rg_skin_grp_2, method='pairwise')
#
# mod_skin_grp_3 <- lm(FA_skin_grouped ~ sex, data=skin)
# manova_skin_grp_3 <- Anova(mod_skin_grp_3)
# summary(manova_skin_grp_3)
#
# mod_skin_grp_4 <- lm(FA_skin_grouped ~ location + creek, data=skin)
# manova_skin_grp_4 <- Anova(mod_skin_grp_4)
# summary(manova_skin_grp_4)
# rg_skin_grp_4 <-ref.grid(mod_skin_grp_4)
# contrast(rg_skin_grp_4, method='pairwise')
#
# mod_skin_grp_5 <- lm(FA_skin_grouped ~ location + creek + location:creek, data=skin) #interaction term
# manova_skin_grp_5 <- Anova(mod_skin_grp_5)
# summary(manova_skin_grp_5)
# rg_skin_grp_5 <-ref.grid(mod_skin_grp_5)
# contrast(rg_skin_grp_5, method='pairwise')

########################################   NEW   ###########################
#Previously, the best model was the last one (included location, creek, and their interaction).
#Now, the creek term is no longer significant by itself, and neither is the sex term.
#here, too, the first model is now the best one.

#setup_indivPUFAs_skin####

#arcsin-sqrt transform proportion data for MANOVA
X18.2n6c_skin_trans <- asin(sqrt(skin$X18.2n6c))
X20.4n6_skin_trans <- asin(sqrt(skin$X20.4n6))
X20.5n3_skin_trans <- asin(sqrt(skin$X20.5n3))
X22.5n3_skin_trans <- asin(sqrt(skin$X22.5n3))
X22.6n3_skin_trans <- asin(sqrt(skin$X22.6n3))

FA_skin_indiv <- cbind(X18.2n6c_skin_trans, X20.4n6_skin_trans, X20.5n3_skin_trans,
                         X22.5n3_skin_trans, X22.6n3_skin_trans)
colnames(FA_skin_indiv) <- c('X18.2n6c', 'X20.4n6', 'X20.5n3', 'X22.5n3', 'X22.6n3')

##MANOVA_indivPUFAs_skin####

mod_skin_indiv_1 <- lm(FA_skin_indiv ~ location, data=skin)
summary(mod_skin_indiv_1)
manova_skin_indiv_1 <- Anova(mod_skin_indiv_1)
summary(manova_skin_indiv_1)
rg_skin_indiv_1 <-ref.grid(mod_skin_indiv_1)
contrast(rg_skin_indiv_1, method='pairwise')

# mod_skin_indiv_2 <- lm(FA_skin_indiv ~ creek, data=skin)
# manova_skin_indiv_2 <- Anova(mod_skin_indiv_2)
# summary(manova_skin_indiv_2)
# rg_skin_indiv_2 <-ref.grid(mod_skin_indiv_2)
# contrast(rg_skin_indiv_2, method='pairwise')
#
# mod_skin_indiv_3 <- lm(FA_skin_indiv ~ sex, data=skin)
# manova_skin_indiv_3 <- Anova(mod_skin_indiv_3)
# summary(manova_skin_indiv_3)
#
# mod_skin_indiv_4 <- lm(FA_skin_indiv ~ location + creek, data=skin)
# manova_skin_indiv_4 <- Anova(mod_skin_indiv_4)
# summary(manova_skin_indiv_4)
# rg_skin_indiv_4 <-ref.grid(mod_skin_indiv_4)
# contrast(rg_skin_indiv_4, method='pairwise')
#
# mod_skin_indiv_5 <- lm(FA_skin_indiv ~ location + creek + location:creek, data=skin)
# manova_skin_indiv_5 <- Anova(mod_skin_indiv_5)
# summary(manova_skin_indiv_5)
# rg_skin_indiv_5 <-ref.grid(mod_skin_indiv_5)
# contrast(rg_skin_indiv_5, method='pairwise')
