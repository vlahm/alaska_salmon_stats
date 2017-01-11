#Jordan Lee - Alaska FA plots
#contact: Mike Vlah (vlahm@uw.edu)
#last edit: 2/22/16

setwd('C:/Users/Mike/git/j_lee_stats_2015')
dev.off()
rm(list=ls())
windows(record=T)

#setup####
data <- read.csv('jl_plotdata.csv')
bylocation <- data[1:3,1:7]
byloc_and_creek <- data[1:5,8:14]
colnames(byloc_and_creek) <- colnames(bylocation)
bylocation[,1] <- 1:3
byloc_and_creek[,1] <- c(1,2,2,3,3)

percentifier <- function(x){
  x <- x * 100
}

error.bars <- function(x, y, upper, lower=upper, cap.length=0.05, horiz=F,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("One or more vectors is not the same length")

  if(horiz==F) {
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=cap.length, ...)
  } else if (horiz==T) {
    arrows(x+upper,y, x-lower, y, angle=90, code=3, length=cap.length, ...)
  }
}

bylocation[,2:7] <- sapply(bylocation[,2:7], FUN = percentifier)
byloc_and_creek[,2:7] <- sapply(byloc_and_creek[,2:7], FUN = percentifier)

#Plots####
#Fig 2a
par(mar=c(4,4,4,6))#, xpd=TRUE)
plot(bylocation$location, bylocation$safa_mean, type='b', pch=16, ylim=c(0,80),
     xlab='', ylab='Relative Abundance (%)', xaxt='n', cex=1.5,
     main='FA concentration across tissue and creek')
axis(side=1, at=1:3, labels=c('Lake Entry', 'Holding', 'Post-Spawn'))
# axis(side=4, at=c(0, 0.25, 0.5, 0.75, 1), labels=c('0', '0.25', '0.5', '0.75', '1'))
# legend(x='right', inset=c(-0.25, 0), legend=c('SAFA', 'MUFA', 'PUFA'), pch=c(16, 17, 18))
lines(bylocation$location, bylocation$mufa_mean, type='b', pch=17, cex=1.5)
lines(bylocation$location, bylocation$pufa_mean, type='b', pch=18, cex=1.5)
error.bars(bylocation$location, bylocation$safa_mean, bylocation$safa_sd)
error.bars(bylocation$location, bylocation$mufa_mean, bylocation$mufa_sd)
error.bars(bylocation$location, bylocation$pufa_mean, bylocation$pufa_sd)
legend(x='topleft', legend=c('SAFA', 'MUFA', 'PUFA'), pch=c(16, 17, 18))
par(mar=c(4,4,4,4))#, xpd=FALSE)

#w3:w6 plot
plot(bylocation$location, data$w3w6_mean[1:3], type="l",col="black", xaxt='n',
     xlab='',ylab='w3:w6', lty=1, ylim=c(min(data$w3w6_mean[1:3])-min(data$w3w6_sd[1:3]),
                                    max(data$w3w6_mean[1:3])+max(data$w3w6_sd[1:3])),
     main='w3:w6 across tissue and creek')
axis(side=1, at=1:3, labels=c('Lake Entry', 'Holding', 'Post-Spawn'))
error.bars(bylocation$location, data$w3w6_mean[1:3], data$w3w6_sd[1:3])


#Fig 2a (including creek) [still in progress]
par(mar=c(4.5,4.5,4,4), xpd=TRUE)
plot(byloc_and_creek$location[c(1,2,4)], byloc_and_creek$safa_mean[c(1,2,4)], type='l', col='steelblue3', ylim=c(0,80),
     xlab='', ylab='', xaxt='n', las=2, lwd=2, bty='l', yaxs='i')
axis(side=1, at=1:3, labels=c('Lake Entry', 'Holding', 'Post-Spawn'))
mtext('FA %', 2, line=2.5, font=2, cex=1.3)
# legend(x='right', inset=c(-.16, 0), legend=c('SAFA', 'MUFA', 'PUFA', 'Hanson', 'Pick'), pch=c(16, 17, 18, NA, NA), lty=c(NA, NA, NA, 1, 2))
lines(byloc_and_creek$location[c(1,2,4)], byloc_and_creek$mufa_mean[c(1,2,4)], col='darkorchid4', lwd=2)
lines(byloc_and_creek$location[c(1,2,4)], byloc_and_creek$pufa_mean[c(1,2,4)], col='green4', lwd=2)
lines(byloc_and_creek$location[c(1,3,5)], byloc_and_creek$safa_mean[c(1,3,5)], lty=2, col='steelblue3', lwd=2)
lines(byloc_and_creek$location[c(1,3,5)], byloc_and_creek$mufa_mean[c(1,3,5)], lty=2, col='darkorchid4', lwd=2)
lines(byloc_and_creek$location[c(1,3,5)], byloc_and_creek$pufa_mean[c(1,3,5)], lty=2, col='green4', lwd=2)
error.bars(byloc_and_creek$location[c(1,2,4)], byloc_and_creek$safa_mean[c(1,2,4)],
           byloc_and_creek$safa_sd[c(1,2,4)]*1.96/sqrt(9))
error.bars(byloc_and_creek$location[c(1,2,4)], byloc_and_creek$mufa_mean[c(1,2,4)],
           byloc_and_creek$mufa_sd[c(1,2,4)]*1.96/sqrt(9))
error.bars(byloc_and_creek$location[c(1,2,4)], byloc_and_creek$pufa_mean[c(1,2,4)],
           byloc_and_creek$pufa_sd[c(1,2,4)]*1.96/sqrt(9))
error.bars(byloc_and_creek$location[c(1,3,5)], byloc_and_creek$safa_mean[c(1,3,5)],
           byloc_and_creek$safa_sd[c(1,3,5)]*1.96/sqrt(10))
error.bars(byloc_and_creek$location[c(1,3,5)], byloc_and_creek$mufa_mean[c(1,3,5)],
           byloc_and_creek$mufa_sd[c(1,3,5)]*1.96/sqrt(10))
error.bars(byloc_and_creek$location[c(1,3,5)], byloc_and_creek$pufa_mean[c(1,3,5)],
           byloc_and_creek$pufa_sd[c(1,3,5)]*1.96/sqrt(10))
text(1.5, 14, '***', cex=1.5, pos=1, srt=-10)
text(2.5, 16, '***', cex=1.5, pos=1, srt=15)
text(2.65, 33.5, '***', cex=1.5, pos=1, srt=-29)
text(2.75, 39, '***', cex=1.5, pos=1, srt=-30)
text(2.75, 51.5, '**', cex=1.5, pos=1, srt=22)
text(2.8, 49, '**', cex=1.5, pos=1, srt=22)
par(mar=c(4,4,4,4), xpd=FALSE)
legend(x=.935, y=82.5, legend=c('  SFA', '  MUFA', '  PUFA'),
       fill=c('steelblue3', 'darkorchid4', 'green4'), bty='n')
legend(x=.9, y=69.5, legend=c('Hanson', 'Pick'), lty=1:2, bty='n', lwd=2)
mtext('Time', 1, line=2.5, font=2, cex=1.3)
