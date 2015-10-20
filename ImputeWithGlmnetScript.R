
##
## Imputation script by Martha Hamblin
## Takes a tab delimited file as an argument
## file: col headers = accessions, row headers = markers
## Script transposes the table, ensures numericness and runs the imputation function.
## Feb 2015
##
## Modified 20150320 to remove space from .imputed.txt filename, tranpose after imputation to restore input format
## and skip check.names function while reading input file to avoid changes to variable names.

args <- commandArgs(trailingOnly = TRUE);

writeLines("Loading glmnet library...");
library("glmnet");
library("methods");

chr_file = args[1]


## This is the basic imputation function. Depending on how thw dosage files were made, you may want to add some filters, which are available as ## functions. Make sure you know whether the snps should be by row or by column.

## SNPs are in columns. They come out of perl script in rows.
impute.glmnet <- function(snps){
varRange <- range(0:2)
cvLambda <- exp(-(2:11))
nPred <- min(60, round(ncol(snps) * 0.5))
snpsNoNA <- apply(snps, 2, function(vec){vec[is.na(vec)] <- mean(vec, na.rm=TRUE); return(vec)})
snpsImp <- snps
for(k in 1:ncol(snps)){
	mrkScores <- snps[,k]
	isNA <- which(is.na(mrkScores)==T)
	if (length(isNA) ==0){next}
	if (length(isNA) > (nrow(snps)*0.95)){next}
	if (length(isNA) > 1){
		if (sd(snps[,k], na.rm=TRUE) == 0) {snpsImp[isNA,k] <- snps[-isNA,k][1]} else{
		corMrk <- abs(cov(snps[,k], snps, use="pairwise.complete.obs"))
    	# Retain markers that correlate highly with marker to be imputed
    	predMrk <- setdiff(order(corMrk,decreasing=TRUE)[1:nPred], k)	
		G <- snpsNoNA[,predMrk]
		ans <- cv.glmnet(x = G[-isNA,],y = mrkScores[-isNA]);    
		pred <- predict(ans,s="lambda.min",newx=G[isNA,])
		pred[pred < 0] <- 0
    	pred[pred > 2] <- 2
		snpsImp[isNA,k] <- round(pred,digits = 2)}
		}else{snpsImp[isNA,k] <- mean(snps[,k],na.rm=T)}
		}
	return(snpsImp)	
		}

## nClones is the number of clones for which genotypes were called; it is the max value of N

filterAndImpute <- function(nClones){
	library(glmnet)
	for(chr in 1:19){
	snpstats <- read.table(file=paste("cassava_chr",chr,".snpstats.txt",sep=""),row.names=1)
	snps <- read.table(file=paste("cassava_chr",chr,".snps.txt",sep=""),header=T,row.names=1)
	
	retained <- grep("RETAINED",snpstats[,4])
	snpstats <- snpstats[retained,]
	hiMAF <- which(as.numeric(snpstats[,2]) > 0.5);  ##make all MAFs <= 0.5
	snpstats[hiMAF,2] <- 1 - as.numeric(snpstats[hiMAF,2]);
	
	##filter snps for N as a function of MAF
	p1 <- which(as.numeric(snpstats[,2]) <0.01)
	p2 <- which(as.numeric(snpstats[,2]) >= 0.01 & as.numeric(snpstats[,2]) <0.1)
	lowN1 <- which(as.numeric(snpstats[,1]) < 0.4)
	lowN2 <- which(as.numeric(snpstats[,1]) < 0.3)
	
	snpsToRemove <- c(intersect(p1,lowN1),intersect(p2,lowN2))
	
	snps <- snps[,-snpsToRemove];
	       

	save(snpsImp,file=paste("cassava_chr",chr,".snpsImp.Rdata",sep=""))
	}
}

filterAndImpute1 <- function(nClones){
	library(glmnet)
	for(chr in 1:19){
		snps <- filterByMAF(chr);	
		snps <- as.matrix(snps)
		mode(snps)="numeric"
		snps <- t(snps);
		snpsImp <- impute.glmnet(snps)

		save(snpsImp,file=paste("cassava_chr",chr,".snpsImp1.Rdata",sep=""))
		}
	}
	
filterAndImpute2 <- function(nClones){
	library(glmnet)
	for(chr in 1:19){
		snps <- filterByMAF(chr);	
		mode(snps)="numeric"
		snps <- t(snps);
		snps <- setPoorDosageToNA(snps)
		
		snpsImp <- impute.glmnet(snps)

		save(snpsImp,file=paste("cassava_chr",chr,".snpsImp2.Rdata",sep=""))
		}
	}



filterByMAF <- function(chr){
	
	snpstats <- read.table(file=paste("cassava_chr",chr,".snps.txt.stats",sep=""),row.names=1)
	snps <- read.table(file=paste("cassava_chr",chr,".snps.txt",sep=""),header=T,row.names=1)
	
	retained <- grep("RETAINED",snpstats[,4])
	snpstats <- snpstats[retained,]
	hiMAF <- which(as.numeric(snpstats[,2]) > 0.5);  ##make all MAFs <= 0.5
	snpstats[hiMAF,2] <- 1 - as.numeric(snpstats[hiMAF,2]);
	
	##filter snps for N as a function of MAF
	p1 <- which(as.numeric(snpstats[,2]) <0.01)
	p2 <- which(as.numeric(snpstats[,2]) >= 0.01 & as.numeric(snpstats[,2]) <0.1)
	lowN1 <- which(as.numeric(snpstats[,1]) < 0.4)
	lowN2 <- which(as.numeric(snpstats[,1]) < 0.3)
	
	snpsToRemove <- c(intersect(p1,lowN1),intersect(p2,lowN2))
	
	
	if(length(snpsToRemove) > 0){
		snps <- snps[-snpsToRemove,];
	}
	snpList <- row.names(snps)
	save(snpList,file=paste("cassava_chr",chr,".snplist.Rdata",sep=""))
	return(snps)
}

setPoorDosageToNA <- function(snps){
snpIndex <- c(1:nrow(snps))
for(i in 1:ncol(snps)){
	ref <- which(snps[,i] <= 0.1)
	het <- which(snps[,i] >= 0.9 & snps[,i] <= 1.1)
	alt <- which(snps[,i] >= 1.9)
	setNA <- snpIndex[-c(ref,het,alt)]
	snps[setNA,i] <- NA
	}
	return(snps)
}


mergeDosage <- function(snps1,snps2){
	matchSNPs <- match(colnames(snps1),colnames(snps2))
	snps2 <- snps2[,matchSNPs]
	snps1 <- rbind(snps1,snps2)	
	return(snps1)
}



setPoorDosageToNA2 <- function(snps){
for(j in 1:ncol(snps)){
	a <- which(snps[,j] > 0.1 & snps[,j] < 0.9)
	b <- which(snps[,j] > 1.1 & snps[,j] < 1.9)
	c <- c(a,b)
	if(length(c) > 0){
 		snps[j,c] <- NA
 	}
  }
  return(snps)
 }
 
 
 
 compareImputatedToPoorDosage <- function(snps,snpsImp){
   temp <- array(dim=ncol(snpsImp))
   for(j in 1:ncol(snps)){
	a <- which(snps[,j] > 0.1 & snps[,j] < 0.9)
	b <- which(snps[,j] > 1.1 & snps[,j] < 1.9)
	c <- c(a,b)
	a <- which(snpsImp[c,j] > 0.1 & snpsImp[c,j] < 0.9)
	b <- which(snpsImp[c,j] > 1.1 & snpsImp[c,j] < 1.9)
	temp[j] <- length(c) - (length(a) + length(b))
   }
   return(temp)
 }

# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

chr_file = trim(chr_file)
 
writeLines(paste("Reading file ",chr_file,"..."));

chr <- read.table(chr_file, sep="\t", row.names=1, header=TRUE, check.names=FALSE)

writeLines("Transposing...")

chrtm <- t(chr)

#remove first line - snp_names
#chrtm_wo <- chrtm[2:nrow(chrtm),]

#writeLines(chrtm)
writeLines("Making numeric...")

mode(chrtm) <- "numeric"
 
writeLines("Imputing...");
chrtmi = impute.glmnet(chrtm)

writeLines("Transposing back to original format...")

chri <- t(chrtmi)

writeLines("Writing output...");
write.table(chri, file = paste(chr_file, ".imputed.txt", sep=""), sep="\t", quote=FALSE);# col.names=snp_names);

writeLines("Done.");
