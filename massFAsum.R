#change this path so that it points to the folder "j_lee_stats_2015" on your computer
setwd("C:/Users/Mike/git/j_lee_stats_2015")
setwd("/home/mike/git/alaska_salmon_stats")

#run all of this
mass <- read.csv('massFAsum.csv')
mass$Stream <- as.vector(mass$Stream)
mass[mass=='na'] <- NA
mass$Stream <- factor(mass$Stream)

error.bars <- function(x, y, upper, lower=upper, cap.length=0.1, horiz=F,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
        stop("One or more vectors is not the same length")

    if(horiz==F) {
        arrows(x,y+upper, x, y-lower, angle=90, code=3, length=cap.length, ...)
    } else if (horiz==T) {
        arrows(x+upper,y, x-lower, y, angle=90, code=3, length=cap.length, ...)
    }
}

#run this section to see a plot of mass by timepoint
#(make sure the plot window is big enough to see labels)
means <- tapply(mass$Mass.FA.sum, mass$TimePt, mean)
sds <- tapply(mass$Mass.FA.sum, mass$TimePt, sd)
mids <- barplot(means, ylim=c(0, max(means+sds)), ylab='mass FA sum')
error.bars(mids, means, sds)

#run this section to see a plot of mass by stream
#(make sure the plot window is big enough to see labels)
means1 <- tapply(mass$Mass.FA.sum, mass$Stream, mean, na.rm=TRUE)
sds1 <- tapply(mass$Mass.FA.sum, mass$Stream, sd, na.rm=TRUE)
mids1 <- barplot(means1, ylim=c(0, max(means1+sds1)), ylab='mass FA sum')
error.bars(mids1, means1, sds1)

#run this section to see a plot of mass by timepoint*stream
#(make sure the plot window is big enough to see labels)
sites <- as.vector(mass$Stream)
sites[is.na(sites)] <- 'Entry'
sites <- factor(sites)

pdf('massFAsum.pdf', width=6, height=5)
defpar <- par(mar=c(4,4,1,0))
means2 <-aggregate(mass$Mass.FA.sum, list(mass$TimePt, sites), mean)
SEs2 <- aggregate(mass$Mass.FA.sum, list(mass$TimePt, sites),
                  function(i) sd(i)/sqrt(length(i)))
mids2 <- barplot(means2[,3][c(1,2,4,3,5)], ylim=c(0, max(means2[,3]+sds2[,3])),
     # names.arg=paste(means2[,1], means2[,2]),
     ylab='', xlim=c(0,10),
     names.arg='', space=c(.2,.8,.2,.8,.2),
     density=c(-1,-1,30,-1,30), legend.text=c('','Hansen','Pick'),
     width=1.1, col=c('white',rep('gray40',4)), lwd=1.5,
     args.legend=list(x=10.1, y=240, bty='n',
                      border=c('white','black','black')))
error.bars(mids2, means2[,3][c(1,2,4,3,5)], SEs2[,3][c(1,2,4,3,5)],
           cap.length=.05)
axis(1, at=c(.8,3.4,6.7), labels=c('Lake entry', 'Holding', 'Post-spawn'),
     tick=FALSE, line=NA, padj=-1)
mtext('Time', 1, line=2, font=2, cex=1.3, at=4)
mtext('FA %', 2, line=2.4, font=2, cex=1.3)
par(defpar)
dev.off()

#does time alone affect mass?
mod1 <- lm(Mass.FA.sum ~ TimePt, data=mass)
anova(mod1) #no, p=0.1516; F=1.9679; df=2,45

#does stream alone affect mass?
mod2 <- lm(Mass.FA.sum ~ Stream, data=mass)
anova(mod2) #no, p=0.1018; F=2.8193; df=1,36

#does sex alone affect mass?
mod3 <- lm(Mass.FA.sum ~ Sex, data=mass)
anova(mod3) #no, p=0.8678; F=0.028; df=1,46

#what about stream and time?
mod4 <- lm(Mass.FA.sum ~ Stream + TimePt, data=mass)
anova(mod4) #no, let me know if want p/F/df values for this

#stream and sex?
mod5 <- lm(Mass.FA.sum ~ Stream + Sex, data=mass)
anova(mod5) #no, let me know if want p/F/df values for this

#sex and time?
mod6 <- lm(Mass.FA.sum ~ TimePt + Sex, data=mass)
anova(mod6) #no, let me know if want p/F/df values for this

#all three?
mod7 <- lm(Mass.FA.sum ~ ., data=mass)
anova(mod7) #no, let me know if want p/F/df values for this

#testing for interactions:
modx <- lm(Mass.FA.sum ~ Stream * Sex, data=mass) #stream + sex + interaction
anova(modx)
modx <- lm(Mass.FA.sum ~ Stream * TimePt, data=mass) #stream + time + interaction
anova(modx)
modx <- lm(Mass.FA.sum ~ Sex * TimePt, data=mass) #sex + time + interaction
anova(modx)
modx <- lm(Mass.FA.sum ~ . * ., data=mass) #everything plus interactions
anova(modx)

#non significant across the board. let me know if you want p/F/df for any of these
