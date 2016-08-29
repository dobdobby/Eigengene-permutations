##########################################################################
##	permute Eigengenes to find sig between-module changes with lifespan	##
##	Adam Dobson, 2016. A.Dobson@ucl.ac.uk								##
##########################################################################
	
	#set path to the eigengenes data file
pathToData <- "./github/eigengenes.csv"

	#open the data
d <- read.csv(pathToData)
head(d)
rownames(d) <- d$X
d <- d[,2:ncol(d)]


	#observed changes in eigengene-eigengene correlations
EAA_SY <- cor(d[d$diet=="EAA",3:ncol(d)], method="spearman") - cor(d[d$diet=="SY",3:ncol(d)], method="spearman")
EAA_rapa <- cor(d[d$diet=="EAA",3:ncol(d)], method="spearman") - cor(d[d$diet=="rapa",3:ncol(d)], method="spearman")

reord <- order(as.numeric(substr(rownames(EAA_SY), 3, 5)))
EAA_SY_obs <- EAA_SY[reord,reord]
EAA_rapa_obs <- EAA_rapa[reord, reord]

myMat <- EAA_SY
myMat[row(myMat) > col(myMat)] <- EAA_rapa[row(myMat) > col(myMat)]

diag(myMat) <- NA

heatmap.2(myMat, trace="n", col=colorRampPalette(c("blue", "white", "red"))(n=1000), Rowv=F, Colv=F, na.color="grey", density="n", cellnote=round(myMat, 2), notecex=0.5, notecol="purple")

heatmap.2(EAA_SY_obs, trace="n", col=colorRampPalette(c("blue", "white", "red"))(n=1000), Rowv=F, Colv=F, na.color="grey", density="n", cellnote=round(EAA_SY_obs, 2), notecex=0.5, notecol="purple", dendrogram="none", key=F)

heatmap.2(EAA_rapa_obs, trace="n", col=colorRampPalette(c("blue", "white", "red"))(n=1000), Rowv=F, Colv=F, na.color="grey", density="n", cellnote=round(EAA_rapa_obs, 2), notecex=0.5, notecol="purple", dendrogram="none", key=F)

	#function to make unique permutations
rperm <- function(m, size=2) { # Obtain m unique permutations of 1:size
    max.failures <- 10

    # Function to index into the upper-level cache.
    prefix <- function(p, k) {    # p is a permutation, k is the prefix size
        sum((p[1:k] - 1) * (size ^ ((1:k)-1))) + 1
    } # Returns a value from 1 through size^k

    # Function to obtain a new permutation.
    newperm <- function() {
        # References cache, k.head, and failures in parent context.
        # Modifies cache and failures.        

        count <- 0                # Protects against infinite loops
        repeat {
            # Generate a permutation and check against previous ones.
            p <- sample(1:size)
            k <- prefix(p, k.head)
            ip <- cache[[k]]
            hash.p <- paste(tail(p,-k.head), collapse="")
            if (is.null(ip[[hash.p]])) break

            # Prepare to try again.
            n.failures <<- n.failures + 1
            count <- count+1
            if (count > max.failures) {  
                p <- NA           # NA indicates a new permutation wasn't found
                hash.p <- ""
                break
            }
        }
        if (count <= max.failures) {
            ip[[hash.p]] <- TRUE      # Update the list of permutations found
            cache[[k]] <<- ip
        }
        p                         # Return this (new) permutation
    }

    # Initialize the cache.
    k.head <- min(size-1, max(1, floor(log(m / log(m)) / log(size))))
    cache <- as.list(1:(size^k.head))
    for (i in 1:(size^k.head)) cache[[i]] <- list()

    # Count failures (for benchmarking and error checking).
    n.failures <- 0

    # Obtain (up to) m unique permutations.
    s <- replicate(m, newperm())
    s[is.na(s)] <- NULL
    list(failures=n.failures, sample=matrix(unlist(s), ncol=size))
} # Returns an m by size matrix; each row is a permutation of 1:size.

	#make matrices of 10000 unique permutations
set.seed(1984)
nPerm <- 10000
myPermsEAA <- rperm(nPerm, size=nrow(subset(d, diet=="EAA")))
myPermsSY <- rperm(nPerm, size=nrow(subset(d, diet=="SY")))
myPermsRapa <- rperm(nPerm, size=nrow(subset(d, diet=="rapa")))

	#save the data at this point
dd <- d

theNewIndex <- d[,1:2]
theData <- d[,3:ncol(d)]

	#remove "ME" from eigengene names
colnames(theData) <- gsub("ME", "", colnames(theData))
theData <- theData[,order(as.numeric(colnames(theData)))]

	#split the data by conditions
EAA <- theData[theNewIndex$diet=="EAA",] 
RAPA <- theData[theNewIndex$diet=="rapa",]
SY <- theData[theNewIndex$diet=="SY",]

	#list to recieve permutation results
permResults <- list(sy=list(), rapa=list())

	#permute the Eigengenes 
for( i in 1:ncol(theData)){
	
	modname <- colnames(theData)[i]
	
	permed_sy <- matrix(ncol=length(colnames(theData)), nrow=nrow(myPermsSY$sample))
	permed_rapa <- matrix(ncol=length(colnames(theData)), nrow=nrow(myPermsRapa$sample))
	
for(xx in 1:nPerm){
	
	eaa <- cbind(EAA[myPermsEAA$sample[xx,],i], EAA[,-i])
	rapa <- cbind(RAPA[myPermsRapa$sample[xx,],i], RAPA[,-i])
	sy <- cbind(SY[myPermsSY$sample[xx,],i], SY[,-i])
	
	colnames(eaa)[1] <- colnames(rapa)[1] <- colnames(sy)[1] <- modname
	
	rapaDiff <- cor(eaa, method="spearman") - cor(rapa, method="spearman")
	syDiff <- cor(eaa, method="spearman") - cor(sy, method="spearman")

	
	permed_sy[xx,] <- syDiff[1,]
	permed_rapa[xx,] <- rapaDiff[1,]
}
colnames(permed_rapa) <- colnames(permed_sy) <- c(modname, colnames(theData)[colnames(theData)!=modname])
permResults$sy[i] <- list(permed_sy)
permResults$rapa[i] <- list(permed_rapa)
names(permResults$sy)[i] <- names(permResults$rapa)[i] <- modname
}

	#determine significance (changes that occur <5% of the time)
EAA_rapa_sig <- matrix(nrow=nrow(EAA_rapa_obs), ncol=ncol(EAA_rapa_obs), dimnames=dimnames(EAA_rapa_obs))

for(i in 1:nrow(EAA_rapa_sig)){
	for(j in 1:ncol(EAA_rapa_sig)){
		obsVal <- EAA_rapa_obs[i,j]
		if(!is.na(obsVal)){
			testMod <- gsub("ME", "", colnames(EAA_rapa_sig))[j]
			permRes <- permResults$rapa[i][[1]][,colnames(permResults$rapa[i][[1]])==testMod]
			EAA_rapa_sig[i,j] <- obsVal >= quantile(permRes, 0.975, na.rm=T) | obsVal <= quantile(permRes, 0.025, na.rm=T)
		}
	}
}
EAA_rapa_sig

EAA_SY_sig <- matrix(nrow=nrow(EAA_SY_obs), ncol=ncol(EAA_SY_obs), dimnames=dimnames(EAA_SY_obs))

for(i in 1:nrow(EAA_SY_sig)){
	for(j in 1:ncol(EAA_SY_sig)){
		obsVal <- EAA_SY_obs[i,j]
		if(!is.na(obsVal)){
			testMod <- gsub("ME", "", colnames(EAA_SY_sig))[j]
			permRes <- permResults$sy[i][[1]][,colnames(permResults$sy[i][[1]])==testMod]
			EAA_SY_sig[i,j] <- obsVal >= quantile(permRes, 0.975, na.rm=T) | obsVal <= quantile(permRes, 0.025, na.rm=T)
		}
	}
}
EAA_SY_sig
EAA_rapa_sig

diag(EAA_SY_sig) <- NA
diag(EAA_rapa_sig) <- NA

EAA_SY_sig[upper.tri(EAA_SY_sig)] == EAA_SY_sig[lower.tri(EAA_SY_sig)]

rowSums(EAA_SY_sig, na.rm=T)
rowSums(EAA_rapa_sig, na.rm=T)
rowSums(EAA_SY_sig == T & EAA_rapa_sig == T, na.rm=T)

	#re-make the heatmaps
myMatNoteLogic <- 	EAA_SY_sig
myMatNoteLogic[row(myMatNoteLogic) > col(myMatNoteLogic)] <- EAA_rapa_sig[row(myMatNoteLogic) > col(myMatNoteLogic)]

myMat <- t(myMat[reord, reord])

diag(myMat) <- diag(EAA_SY_obs) <- diag(EAA_rapa_obs) <- NA

dimnames(EAA_SY_obs) <- dimnames(EAA_rapa_obs) <- dimnames(myMat) <- lapply(dimnames(myMat), function(x){substr(x, 3, 4)})

colInd <- {a <- ifelse(myMatNoteLogic==T, round(myMat,2), "");
	ifelse(as.numeric(a[which(a!="")])<0, "red", "blue")}
keyPar <- list(cex=0.45, las=1)

	#make the final plot
pdf("./CorrelationChangeMatrixHeatmap.pdf")
heatmap.2(myMat, trace="n", col=colorRampPalette(c("steelblue1", "white", "hotpink"))(n=1000), Rowv=F, Colv=F, na.color="lightgrey", density="n", cellnote=ifelse(myMatNoteLogic==T, round(myMat,2), ""), notecex=0.68, notecol="black", rowsep=1:nrow(myMat), colsep=1:ncol(myMat), key.par=keyPar, key.xlab="Change in correlation", key.title="")
dev.off()

sessionInfo()