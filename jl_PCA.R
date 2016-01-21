#Jordan Lee
#Alaska Salmon FA PCAs 2016
#Last edit: 1/20/16
#Contact: Mike Vlah (vlahm13@gmail.com)

dev.off()
rm(list=ls())
windows(record=T)
setwd("C:/Users/Mike/Desktop/Grad/Projects/for_others/j_lee_stats_2015")
data_full<-read.csv("PCA_data.csv", header=FALSE)
palette(c('deepskyblue', 'green', 'magenta'))

#functions and packages####

library(vioplot)
library(vegan)
library(plyr)
library(plotrix)
library(car)
library(ade4)

replacer <- function(x, cond, replacement, rows=1:length(x[,1]), cols=1:length(x[1,]), asinsqrt=FALSE){
    for (i in rows)
    {
        for (j in cols)
        {
            if (cond=='==NA')
            {
                if (is.na(x[i,j]==TRUE))
                {
                    x[i,j] <- replacement
                }
            } else
            {
                if (eval(parse(text=paste('x[i,j]', cond))))
                {
                    x[i,j] <- replacement
                }
            }
            if (asinsqrt == TRUE)
            {
                x[i,j]<-asin(sqrt(x[i,j]/100))*(2/pi)
            }
        }
    }
    return(x)
}

pca.eigenval <-
    function(x.pca,dim=length(x.pca$sdev),digits=7){

        #check for dim limit
        if(dim>length(x.pca$sdev)){
            cat("Only",length(x.pca$sdev),"axes available\n")
            dim<-length(x.pca$sdev)
        }

        #calculate some variables
        names<-colnames(x.pca$rotation[,1:dim])
        var<-x.pca$sdev^2
        trace<-sum(var)
        prop.var<-var/trace

        #broken-stick distribution
        p<-length(x.pca$sdev)
        y<-rep(0,p)
        for(i in 1:p) y[i]<-1/i
        for(i in 1:p) y[i]<-sum(y[i:p])
        y<-y[1:dim]

        #print results
        cat('Importance of components:\n')
        z<-rbind('Variance(eigenvalue)'=var[1:dim],
                 'Proportion of Variance'=prop.var[1:dim],
                 'Cumulative Proportion'=cumsum(prop.var[1:dim]),
                 'Broken-stick value'=y)
        colnames(z)<-names
        z<-round(z,digits=digits)
        return(z)
    }

pca.eigenvec <-
    function(x.pca,dim=length(x.pca$sdev),
             digits=7,cutoff=0){

        #check for dim limit
        if(dim>ncol(x.pca$rotation)){
            cat("Only",ncol(x.pca$rotation),"axes available\n")
            dim<-ncol(x.pca$rotation)
        }

        #print results
        cat("\nEigenvectors:\n")
        z<-format(round(x.pca$rotation[,1:dim],digits=digits))
        z[abs(x.pca$rotation[,1:dim])<cutoff]<-substring('',1,nchar(z[1,1]))
        z<-as.data.frame(z)
        return(z)
    }

pca.structure <-
    function(x.pca,x,dim=length(x.pca$sdev),
             digits=3,cutoff=0){

        #check for dim limit
        if(dim>length(x.pca$sdev)){
            cat("Only",length(x.pca$sdev),"axes available\n")
            dim<-length(x.pca$sdev)
        }

        #calculate structure correlations
        z<-cor(x,x.pca$x[,1:dim])

        #print results
        cat("\nStructure Correlations:\n")
        z<-round(z,digits=digits)
        z[abs(z)<cutoff]<-substring('',1,nchar(z[1,1]))
        z<-as.data.frame(z)
        return(z)
    }

ordi.monte <-
    function(x,ord,dim=length(x),perm=1000,center=TRUE,
             scale=TRUE,digits=3,plot=TRUE,col.hist='blue',col.line='red',
             lty=2,las=1,lab=c(5,5,4),...){

        p<-length(x)
        if(dim>p){
            cat("Only",p,"axes available\n")
            dim<-p
        }

        if(ord=='pca'){
            z<-prcomp(x,center=center,scale=scale) #prcomp analysis
            z.eig<-z$sdev[1:dim]^2 #compute eigenvalues
            z.teig<-t(z.eig) #transpose eigenvalues
            z.teig<-t(matrix(z.teig,length(z.teig),perm)) #fill matrix with eigenvalues
            write('',file='y.csv') #empty outfile if it exists
            for(i in 1:perm){
                y<-apply(x,2,sample) #permute data matrix
                y<-prcomp(y,center=center,scale=scale) #prcomp on permuted matrix
                y<-y$sdev[1:dim]^2 #compute eigenvalues
                y<-as.data.frame(t(y)) #coerce to data frame and transpose
                write.table(y,file='y.csv',sep=',',append=TRUE,row.names=FALSE,col.names=FALSE)
            }
            y<-read.table('y.csv',sep=',',header=FALSE) #read in permutation results
            p.value<-apply(y>z.teig,2,sum) #compute proportion of random distribution smaller than observed
            p.value<-p.value/perm #compute p-value
            names<-colnames(z$rotation[,1:dim]) #add 'PC#' names
        }

        else if(ord=='ca'){
            library(vegan)
            z<-cca(x) #correspondence analysis
            z.eig<-z$CA$eig[1:dim] #get eigenvalues
            z.teig<-t(z.eig) #transpose eigenvalues
            z.teig<-t(matrix(z.teig,length(z.teig),perm)) #fill matrix with eigenvalues
            write('',file='y.csv') #empty outfile if it exists
            for(i in 1:perm){
                y<-apply(x,2,sample) #permute data matrix
                y<-cca(y) #CA on permuted matrix
                y<-y$CA$eig[1:dim] #get eigenvalues
                y<-as.data.frame(t(y)) #coerce to data frame and transpose
                write.table(y,file='y.csv',sep=',',append=TRUE,row.names=FALSE,col.names=FALSE)
            }
            y<-read.table('y.csv',sep=',',header=FALSE) #read in permutation results
            p.value<-apply(y>z.teig,2,sum) #compute proportion of random distribution smaller than observed
            p.value<-p.value/perm #compute p-value
            names<-names(z$CA$eig[1:dim]) #add 'CA#' names
        }

        else if(ord=='dca'){
            library(vegan)
            if(dim>4){
                cat("Only",4,"axes available\n")
                dim<-4
            }
            z<-decorana(x,...) #detrended correspondence analysis
            z.eig<-z$evals[1:dim] #get eigenvalues
            z.teig<-t(z.eig) #transpose eigenvalues
            z.teig<-t(matrix(z.teig,length(z.teig),perm)) #fill matrix with eigenvalues
            write('',file='y.csv') #empty outfile if it exists
            for(i in 1:perm){
                y<-apply(x,2,sample) #permute data matrix
                y<-decorana(y,...) #DCA on permuted matrix
                y<-y$evals[1:dim] #get eigenvalues
                y<-as.data.frame(t(y)) #coerce to data frame and transpose
                write.table(y,file='y.csv',sep=',',append=TRUE,row.names=FALSE,col.names=FALSE)
            }
            y<-read.table('y.csv',sep=',',header=FALSE) #read in permutation results
            p.value<-apply(y>z.teig,2,sum) #compute proportion of random distribution smaller than observed
            p.value<-p.value/perm #compute p-value
            names<-names(z$eval[1:dim]) #add 'CA#' names
        }

        if(plot==TRUE){
            for(i in 1:dim){
                xmin<-min(min(y[[i]],z.eig[i]))
                xmax<-max(max(y[[i]],z.eig[i]))
                hist(y[[i]],col=col.hist,las=las,lab=lab,
                     xaxs='i',yaxs='i',xlim=c(xmin,xmax),
                     xlab='Eigenvalue',
                     main=paste('Random Permutation Distribution of Eigenvalues for',names[i],sep=' '),...)
                abline(v=z.eig[i],col=col.line,lty=lty,lwd=2,...)
                readline("Press return for next plot ")
            }
        }

        cat('Randomization Test of Eigenvalues:\n')
        z<-rbind('Eigenvalue'=z.eig,'P-value'=p.value)
        colnames(z)<-names
        z<-round(z,digits=digits)
        return(z)
    }

#muscle - individual FAs####

#prepare musc_indiv data
musc_indiv <- data_full[-c(1:6,8,11,18,24,31:35), 1:49]
names <- musc_indiv[,1]
musc_indiv <- apply(musc_indiv[,-1], MARGIN=2, FUN=as.numeric)
rownames(musc_indiv) <- names

#replace 0s with a small value, arcsine square-root transform all data
musc_indiv <- replacer(musc_indiv, cond='==0', replacement=0.0000000001, asinsqrt=TRUE)

#format data for PCA
musc_indiv<-t(as.matrix(musc_indiv))

#perform PCA and associated tests
musc_indiv_pca<-prcomp(musc_indiv, scale=TRUE, scores=TRUE)
#determine eigenvalues
musc_indiv_eigen<-pca.eigenval(musc_indiv_pca)
#see which eigenvalues are significant
screeplot(musc_indiv_pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
musc_indiv_struc<-pca.structure(musc_indiv_pca, musc_indiv, dim=7, cutoff=0.5)
#sample scores
musc_indiv_scores<-musc_indiv_pca$x[,1:7]
#test loadings for significance (two methods)
#testdim(musc_indiv_pca) #gotta figure out how to convert from 'prcomp' class to 'pca' class


#plot PCA
pchs <- factor(as.vector(as.matrix((data_full[4,2:49]))))
levels(pchs) <- c(21,22,24)
pchs <- as.integer(as.vector(pchs))
bgs <- factor(as.vector(as.matrix((data_full[5,2:49]))))
levels(bgs) <- c('black', 'white')
bgs <- as.vector(bgs)
cols <- factor(as.vector(as.matrix(data_full[3,2:length(data_full[1,])])))

#first and second principal components
musc_indiv_12<-ordiplot(musc_indiv_pca, choices=c(1,2), type="none",
                        main="PC 1 and 2", ylim=c(-6.5, 5))
points(musc_indiv_12, "sites", pch=pchs,
       col=cols, bg=bgs)
arrows(0,0,musc_indiv_pca$rotation[-c(1,6,16,17,19),1]*10,
       musc_indiv_pca$rotation[-c(1,6,16,17,19),2]*10, col="black")
text(musc_indiv_pca$rotation[-c(1,6,16,17,19),1]*12,
     musc_indiv_pca$rotation[-c(1,6,16,17,19),2]*12,
     row.names(musc_indiv_pca$rotation), col="black")
legend("bottomleft", legend=c("Entry", "Holding", "Post-spawn", "Pre-stream", "Hansen",
                              "Pick", "Female", "Male"),
       fill=c(1, 2, 3, "white", "white", "white", "white", "white"), border="white",
       pch=c(NA, NA, NA, 21, 22, 24, 21, 21),
       pt.bg=c(NA, NA, NA, "white", "white", "white", "black", "white"))

#first and third PCs
musc_indiv_13<-ordiplot(musc_indiv_pca, choices=c(1,3), type="none", main="PC 1 and 3")
points(musc_indiv_13, "sites", pch=pchs,
       col=cols, bg=bgs)
arrows(0,0,musc_indiv_pca$rotation[-c(1,4,5,9,15,16,19),1]*10,
       musc_indiv_pca$rotation[-c(1,4,5,9,15,16,19),3]*10, col="black")
text(musc_indiv_pca$rotation[-c(1,4,5,9,15,16,19),1]*12,
     musc_indiv_pca$rotation[-c(1,4,5,9,15,16,19),3]*12,
     row.names(musc_indiv_pca$rotation), col="black")
legend("bottomleft", legend=c("Entry", "Holding", "Post-spawn", "Pre-stream", "Hansen",
                              "Pick", "Female", "Male"),
       fill=c(1, 2, 3, "white", "white", "white", "white", "white"), border="white",
       pch=c(NA, NA, NA, 21, 22, 24, 21, 21),
       pt.bg=c(NA, NA, NA, "white", "white", "white", "black", "white"))

#muscle - grouped FAs####

#prepare musc_grouped data
musc_grouped <- data_full[c(31,33:35), 1:49]
names <- musc_grouped[,1]
musc_grouped <- apply(musc_grouped[,-1], MARGIN=2, FUN=as.numeric)
rownames(musc_grouped) <- names

#replace 0s with a small value, arcsine square-root transform all data
musc_grouped <- replacer(musc_grouped, cond='==0', replacement=0.0000000001, asinsqrt=TRUE)

#format data for PCA
musc_grouped<-t(as.matrix(musc_grouped))

#perform PCA and associated tests
musc_grouped_pca<-prcomp(musc_grouped, scale=TRUE, scores=TRUE)
#determine eigenvalues
musc_grouped_eigen<-pca.eigenval(musc_grouped_pca)
#see which eigenvalues are significant
screeplot(musc_grouped_pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
musc_grouped_struc<-pca.structure(musc_grouped_pca, musc_grouped, dim=7, cutoff=0.5)
#sample scores
musc_grouped_scores<-musc_grouped_pca$x[,1:7]
#test loadings for significance (two methods)
#testdim(musc_grouped_pca) #gotta figure out how to convert from 'prcomp' class to 'pca' class

#plot PCA
pchs <- factor(as.vector(as.matrix((data_full[4,2:49]))))
levels(pchs) <- c(21,22,24)
pchs <- as.integer(as.vector(pchs))
bgs <- factor(as.vector(as.matrix((data_full[5,2:49]))))
levels(bgs) <- c('black', 'white')
bgs <- as.vector(bgs)
cols <- factor(as.vector(as.matrix(data_full[3,2:length(data_full[1,])])))

#first and second principal components
musc_grouped_12<-ordiplot(musc_grouped_pca, choices=c(1,2), type="none",
                        main="PC 1 and 2", ylim=c(-2.5, 4))
points(musc_grouped_12, "sites", pch=pchs,
       col=cols, bg=bgs)
arrows(0,0,musc_grouped_pca$rotation[,1]*3, musc_grouped_pca$rotation[,2]*3, col="black")
text(musc_grouped_pca$rotation[,1]*3.5, musc_grouped_pca$rotation[,2]*3.5,
     row.names(musc_grouped_pca$rotation), col="black")
legend("topleft", legend=c("Entry", "Holding", "Post-spawn", "Pre-stream", "Hansen",
                              "Pick", "Female", "Male"),
       fill=c(1, 2, 3, "white", "white", "white", "white", "white"), border="white",
       pch=c(NA, NA, NA, 21, 22, 24, 21, 21),
       pt.bg=c(NA, NA, NA, "white", "white", "white", "black", "white"))

#skin - individual FAs####

#prepare skin_indiv data
skin_indiv <- data_full[-c(1:6,14,18,21,26,31:35), c(1,50:79)]
names <- skin_indiv[,1]
skin_indiv <- apply(skin_indiv[,-1], MARGIN=2, FUN=as.numeric)
rownames(skin_indiv) <- names

#replace 0s with a small value, arcsine square-root transform all data
skin_indiv <- replacer(skin_indiv, cond='==0', replacement=0.0000000001, asinsqrt=TRUE)

#format data for PCA
skin_indiv<-t(as.matrix(skin_indiv))

#perform PCA and associated tests
skin_indiv_pca<-prcomp(skin_indiv, scale=TRUE, scores=TRUE)
#determine eigenvalues
skin_indiv_eigen<-pca.eigenval(skin_indiv_pca)
#see which eigenvalues are significant
screeplot(skin_indiv_pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
skin_indiv_struc<-pca.structure(skin_indiv_pca, skin_indiv, dim=7, cutoff=0.5)
#sample scores
skin_indiv_scores<-skin_indiv_pca$x[,1:7]
#test loadings for significance
#testdim(skin_indiv_pca) #gotta figure out how to convert from 'prcomp' class to 'pca' class


#plot PCA
pchs <- factor(as.vector(as.matrix((data_full[4,2:49]))))
levels(pchs) <- c(21,22,24)
pchs <- as.integer(as.vector(pchs))
bgs <- factor(as.vector(as.matrix((data_full[5,2:49]))))
levels(bgs) <- c('black', 'white')
bgs <- as.vector(bgs)
cols <- factor(as.vector(as.matrix(data_full[3,2:length(data_full[1,])])))

#first and second principal components
skin_indiv_12<-ordiplot(skin_indiv_pca, choices=c(1,2), type="none",
                        main="PC 1 and 2", ylim=c(-3.5, 5))
points(skin_indiv_12, "sites", pch=pchs,
       col=cols, bg=bgs)
arrows(0,0,skin_indiv_pca$rotation[-c(5,10,15,17),1]*8,
       skin_indiv_pca$rotation[-c(5,10,15,17),2]*8, col="black")
text(skin_indiv_pca$rotation[-c(5,10,15,17),1]*8.5,
     skin_indiv_pca$rotation[-c(5,10,15,17),2]*8.5,
     row.names(skin_indiv_pca$rotation), col="black")
legend("topright", legend=c("Entry", "Holding", "Post-spawn", "Pre-stream", "Hansen",
                              "Pick", "Female", "Male"),
       fill=c(1, 2, 3, "white", "white", "white", "white", "white"), border="white",
       pch=c(NA, NA, NA, 21, 22, 24, 21, 21),
       pt.bg=c(NA, NA, NA, "white", "white", "white", "black", "white"))

#first and third PCs
skin_indiv_13<-ordiplot(skin_indiv_pca, choices=c(1,3), type="none",
                        main="PC 1 and 3", ylim=c(-4, 5))
points(skin_indiv_13, "sites", pch=pchs,
       col=cols, bg=bgs)
arrows(0,0,skin_indiv_pca$rotation[-c(1,5,9,15,19),1]*8,
       skin_indiv_pca$rotation[-c(1,5,9,15,19),3]*8, col="black")
text(skin_indiv_pca$rotation[-c(1,5,9,15,19),1]*8.5,
     skin_indiv_pca$rotation[-c(1,5,9,15,19),3]*8.5,
     row.names(skin_indiv_pca$rotation), col="black")
legend("topleft", legend=c("Entry", "Holding", "Post-spawn", "Pre-stream", "Hansen",
                              "Pick", "Female", "Male"),
       fill=c(1, 2, 3, "white", "white", "white", "white", "white"), border="white",
       pch=c(NA, NA, NA, 21, 22, 24, 21, 21),
       pt.bg=c(NA, NA, NA, "white", "white", "white", "black", "white"))

#skin - grouped FAs####

#prepare skin_grouped data
skin_grouped <- data_full[c(31:35), c(1,50:79)]
names <- skin_grouped[,1]
skin_grouped <- apply(skin_grouped[,-1], MARGIN=2, FUN=as.numeric)
rownames(skin_grouped) <- names

#replace 0s with a small value, arcsine square-root transform all data
skin_grouped <- replacer(skin_grouped, cond='==NA', replacement=0.0000000001, asinsqrt=TRUE)

#format data for PCA
skin_grouped<-t(as.matrix(skin_grouped))

#perform PCA and associated tests
skin_grouped_pca<-prcomp(skin_grouped, scale=TRUE, scores=TRUE)
#determine eigenvalues
skin_grouped_eigen<-pca.eigenval(skin_grouped_pca)
#see which eigenvalues are significant
screeplot(skin_grouped_pca, bstick=T)
#see loadings.  square these to get percentage of variance in each original variable
#accounted for by each principal component
skin_grouped_struc<-pca.structure(skin_grouped_pca, skin_grouped, dim=7, cutoff=0.5)
#sample scores
skin_grouped_scores<-skin_grouped_pca$x[,1:7]
#test loadings for significance (two methods)
#testdim(skin_grouped_pca) #gotta figure out how to convert from 'prcomp' class to 'pca' class

#plot PCA
pchs <- factor(as.vector(as.matrix((data_full[4,2:49]))))
levels(pchs) <- c(21,22,24)
pchs <- as.integer(as.vector(pchs))
bgs <- factor(as.vector(as.matrix((data_full[5,2:49]))))
levels(bgs) <- c('black', 'white')
bgs <- as.vector(bgs)
cols <- factor(as.vector(as.matrix(data_full[3,2:length(data_full[1,])])))

#first and second principal components
skin_grouped_12<-ordiplot(skin_grouped_pca, choices=c(1,2), type="none",
                          main="PC 1 and 2", ylim=c(-3.5, 2.5))
points(skin_grouped_12, "sites", pch=pchs,
       col=cols, bg=bgs)
arrows(0,0,skin_grouped_pca$rotation[,1]*3, skin_grouped_pca$rotation[,2]*3, col="black")
text(skin_grouped_pca$rotation[,1]*3.5, skin_grouped_pca$rotation[,2]*3.5,
     row.names(skin_grouped_pca$rotation), col="black")
legend("bottomright", legend=c("Entry", "Holding", "Post-spawn", "Pre-stream", "Hansen",
                           "Pick", "Female", "Male"),
       fill=c(1, 2, 3, "white", "white", "white", "white", "white"), border="white",
       pch=c(NA, NA, NA, 21, 22, 24, 21, 21),
       pt.bg=c(NA, NA, NA, "white", "white", "white", "black", "white"))

#first and third pcs
skin_grouped_12<-ordiplot(skin_grouped_pca, choices=c(1,3), type="none",
                          main="PC 1 and 3")
points(skin_grouped_12, "sites", pch=pchs,
       col=cols, bg=bgs)
arrows(0,0,skin_grouped_pca$rotation[,1]*3, skin_grouped_pca$rotation[,3]*3, col="black")
text(skin_grouped_pca$rotation[,1]*3.5, skin_grouped_pca$rotation[,3]*3.5,
     row.names(skin_grouped_pca$rotation), col="black")
legend("bottomright", legend=c("Entry", "Holding", "Post-spawn", "Pre-stream", "Hansen",
                           "Pick", "Female", "Male"),
       fill=c(1, 2, 3, "white", "white", "white", "white", "white"), border="white",
       pch=c(NA, NA, NA, 21, 22, 24, 21, 21),
       pt.bg=c(NA, NA, NA, "white", "white", "white", "black", "white"))
