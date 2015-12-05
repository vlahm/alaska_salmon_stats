#Jordan Lee - Alaska FA stats
#contact: Mike Vlah (vlahm@uw.edu)
#last edit: 4 Dec 2015

#NOTE: to test specific hypotheses about parameters, use linearHypothesis()
#details here: https://socserv.socsci.mcmaster.ca/jfox/Books/Companion/appendix/Appendix-Multivariate-Linear-Models.pdf

#NOTES FOR JORDAN:
#use ALT+O to collapse sections.  Click the blue rectangle beside a section to open it.
#I commented out all but the best model for each comparison.  For the first part (MANOVA_FAgroups_muscle), location,
#creek, and the interaction of the two are all significant.  You'll want to acknowledge that in your paper.  if you're 
#not sure about interpretation, talk to me or Gordon.

#For part 2 (MANOVA_indivPUFAs_muscle), there's no significant difference between creeks, but ther IS a difference
#between sexes.  You didn't give me sex data for the skin samples, but we may want to include that too.  You may be
#totally uninterested in the differences between sexes, but otherwise this could be soemthing worth making a graph for too.

#For part 3 (MANOVA_FAgroups_skin), there's an interaction between location and creek again.

#For part 4 (MANOVA_FAindiv_skin), you have a simple situation.  Location is the only source of significant variation.




#setup_FAgroups_muscle####
muscle <- read.csv('jl_data_muscle.csv')
library(car)
library(lsmeans)

#arcsin-sqrt transform proportion data for MANOVA
safa_muscle_trans <- asin(sqrt(muscle$safa_rel))
mufa_muscle_trans <- asin(sqrt(muscle$mufa_rel))
pufa_muscle_trans <- asin(sqrt(muscle$pufa_rel))

FA_muscle_grouped <- cbind(safa_muscle_trans, mufa_muscle_trans, pufa_muscle_trans)
colnames(FA_muscle_grouped) <- c('safa', 'mufa', 'pufa')

##MANOVA_FAgroups_muscle####

###gangster method (ignore unless desperate):

#model1 <- manova(FA ~ location)
#summary.aov(model1) #univariate stats
#summary(model1, test='Pillai')
#pairwise.manova(FA, location, p.method = "fdr") #no longer works?

###the real thing:

#I've gone through and sytematically tested the effects of location, creek, and 
#sex on SAFA, MUFA, and PUFA, in all their combinations.  Where relevant, I've
#included code for pairwise comparisons.  The first series of operations is annotated.

# mod_musc_grp_1 <- lm(FA_muscle_grouped ~ location, data=muscle) #first let's see if location has an effect on SAFA_muscle, MUFA_muscle, OR PUFA_muscle relative composition.
# summary(mod_musc_grp_1) #use this to see results of univariate tests (only included this time)
# manova_musc_grp_1 <- Anova(mod_musc_grp_1) #this is the Anova function from the car package, creates manova object of class 'Anova.mlm'
# summary(manova_musc_grp_1) #Pillai, Wilks, etc. are different methods for testing significance.  This shit is mad significant across the board.
# rg_musc_grp_1 <-ref.grid(mod_musc_grp_1) #the 'reference grid' object is used by the 'lsmeans' package to perform pairwise comparisons
# contrast(rg_musc_grp_1, method='pairwise')

# mod_musc_grp_2 <- lm(FA_muscle_grouped ~ creek, data=muscle)
# manova_musc_grp_2 <- Anova(mod_musc_grp_2)
# summary(manova_musc_grp_2) 
# rg_musc_grp_2 <-ref.grid(mod_musc_grp_2)
# contrast(rg_musc_grp_2, method='pairwise')

# mod_musc_grp_3 <- lm(FA_muscle_grouped ~ sex, data=muscle)
# manova_musc_grp_3 <- Anova(mod_musc_grp_3)
# summary(manova_musc_grp_3) #sex does not have an effect on FA_muscle composition

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

# mod_musc_indiv_1 <- lm(FA_muscle_indiv ~ location, data=muscle)
# summary(mod_musc_indiv_1)
# manova_musc_indiv_1 <- Anova(mod_musc_indiv_1)
# summary(manova_musc_indiv_1)
# rg_musc_indiv_1 <-ref.grid(mod_musc_indiv_1)
# contrast(rg_musc_indiv_1, method='pairwise')

# mod_musc_indiv_2 <- lm(FA_muscle_indiv ~ creek, data=muscle)
# manova_musc_indiv_2 <- Anova(mod_musc_indiv_2)
# summary(manova_musc_indiv_2) 
# rg_musc_indiv_2 <-ref.grid(mod_musc_indiv_2)
# contrast(rg_musc_indiv_2, method='pairwise')

# mod_musc_indiv_3 <- lm(FA_muscle_indiv ~ sex, data=muscle)
# manova_musc_indiv_3 <- Anova(mod_musc_indiv_3)
# summary(manova_musc_indiv_3)

mod_musc_indiv_4 <- lm(FA_muscle_indiv ~ location + sex, data=muscle)
manova_musc_indiv_4 <- Anova(mod_musc_indiv_4)
summary(manova_musc_indiv_4)
rg_musc_indiv_4 <-ref.grid(mod_musc_indiv_4)
contrast(rg_musc_indiv_4, method='pairwise')

# mod_musc_indiv_5 <- lm(FA_muscle_indiv ~ location + sex + location:sex, data=muscle)
# manova_musc_indiv_5 <- Anova(mod_musc_indiv_5)
# summary(manova_musc_indiv_5)
# rg_musc_indiv_5 <-ref.grid(mod_musc_indiv_5)
# contrast(rg_musc_indiv_5, method='pairwise')

#setup_FAgroups_skin####
skin <- read.csv('jl_data_skin.csv')

#arcsin-sqrt transform proportion data for MANOVA
safa_skin_trans <- asin(sqrt(skin$safa_rel))
mufa_skin_trans <- asin(sqrt(skin$mufa_rel))
pufa_skin_trans <- asin(sqrt(skin$pufa_rel))

FA_skin_grouped <- cbind(safa_skin_trans, mufa_skin_trans, pufa_skin_trans)
colnames(FA_skin_grouped) <- c('safa', 'mufa', 'pufa')

##MANOVA_FAgroups_skin####

# mod_skin_grp_1 <- lm(FA_skin_grouped ~ location, data=skin)
# summary(mod_skin_grp_1)
# manova_skin_grp_1 <- Anova(mod_skin_grp_1)
# summary(manova_skin_grp_1)
# rg_skin_grp_1 <-ref.grid(mod_skin_grp_1)
# contrast(rg_skin_grp_1, method='pairwise')
# 
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

mod_skin_grp_5 <- lm(FA_skin_grouped ~ location + creek + location:creek, data=skin) #interaction term
manova_skin_grp_5 <- Anova(mod_skin_grp_5)
summary(manova_skin_grp_5)
rg_skin_grp_5 <-ref.grid(mod_skin_grp_5)
contrast(rg_skin_grp_5, method='pairwise')

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
# 
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
