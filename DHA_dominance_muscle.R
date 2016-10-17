
library(car)
library(lsmeans)
# library(dplyr)
library(devtools)
install_github('vlahm/manipulateR')
library(manipulateR)
library(stringr)

setwd('/home/mike/git/alaska_salmon_stats')
dha <- read.csv('DHA_dominance_muscle.csv')

#remove FAs with < 40% representation
dha[,-(1:3)] <- matfilter(dha[,-(1:3)], cond='==0', thresh=.6)

#simplify "sample.ID" column
dha$Sample.ID <- str_match(dha$Sample.ID, '(\\d+)-.*')[,2]

