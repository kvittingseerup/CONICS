colMeans = function (expmat) {
  n = apply(expmat,2,mean)
  return(n)
}

rowMeans = function (expmat) {
  n = apply(expmat,1,mean)
  return(n)
}


#' Calculate the BIC (information-theoretic criteria) for a mixture model.
#'
#' This function calculate the BIC (information-theoretic criteria) for a mixture model.
#' @param l1 Log likelihood 1
#' @param l2 Log likelihood 1
#' @param n Number of total tests for Benfferroni correction
#' @keywords Detect varibale genes
#' @export
#' @examples
#' detectVarGenes(suva_exp,500)

BIC.mix = function(out){
	bc=-2*out$loglik+log(length(out$x))*(length(unlist(out[2:4]))-1)
	return(bc)
}

#' Perform likelihood ratio test given two log likelihoods and a integer numbe rindicating the number of tests performed for multiple testing correction
#'
#' This function performs a likelihood ratio test given two log likelihoods and a integer numbe rindicating the number of tests performed for multiple testing correction
#' @param l1 Log likelihood 1
#' @param l2 Log likelihood 1
#' @param n Number of total tests for Benfferroni correction
#' @keywords Detect varibale genes
#' @export
#' @examples
#' detectVarGenes(suva_exp,500)

likelihoodRatioTest = function(l1,l2,n){
	lrt= -2*l1 -(-2*l2)
	res=min(pchisq(lrt, df=2, lower.tail=FALSE)*n,1)
	#print(n) # KVS outcomment
	return(res)
}


#' Test difference of the z-score distribution in the estimated mixed gausian model
#'
#' Test difference of the z-score distribution in the estimated mixed gausian model
#' both using mean and median as well as the summary statistics. Not meant to be used by the end user
#' @param set1 the z-score values of the first estimated gausian model.
#' @param set1 the z-score values of the second estimated gausian model.
#' @keywords NA
#' @export
#' @examples
#' detectVarGenes(suva_exp,500)

testWrapper <- function(set1, set2) {
    data.frame(
        n1 = length(set1),
        n2 = length(set2),
        mean1 = mean(set1),
        mean2 = mean(set2),
        median1 = median(set1),
        median2 = median(set2),
        mean_shift = mean(set2) - mean(set1),
        median_shift = median(set2) - median(set1),
        wilcox_p_value = suppressWarnings(wilcox.test(set1, set2)$p.value),
        ttest_p_value  = t.test(set1, set2)$p.value,
        stringsAsFactors=F
    )
}


#' Extract coordinats of a sliding window across expressed genes
#'
#' This function extract coordinats of a sliding window across expressed genes. It uses the matrix supplied to \'mat\'
#' to extract expressed genes and then generates the regions corresponding a sliding window across the genome
#' (defined by fracChrGenesInWindow and windowSlideFraction) which can directly be used with the \code{plotAll} function.
#' @param mat A genes X samples expression matrix of log2(CPM/10+1) scaled (single cell) RNA-seq counts.
#' @param gene_pos A matrix of positions for every gene is the expression matrix, calculated via getGenePositions()
#' @param fracChrGenesInWindow A fraction indicating the fraction of genes included in each window
#' @param windowSlideFraction A fraction indicating the number of genes to move the sliding window as a fraction of the window size (defined via fracChrGenesInWindow)
#' @param minWindowSize A minimum number of genes in each region.
#' @keywords Get regions
#' @export
#' @examples
#' detectVarGenes(suva_exp,500)

extractSlidingRegionsFromGenes <- function(mat, gene_pos, fracChrGenesInWindow = 0.1, windowSlideFraction = 0.25, minWindowSize = 100 ) {
    expressedGenes <- gene_pos[which(
        gene_pos$ensembl_gene_id %in% rownames(mat) |
            gene_pos$hgnc_symbol %in% rownames(mat)
    ),]
    chrGeneList <- split(expressedGenes, expressedGenes$chromosome_name)

    chrRegionList <- plyr::llply(chrGeneList, .inform = TRUE, function(aChr) { # aChr <- chrGeneList[[18]]
        if( nrow(aChr) < minWindowSize ) {
            return(NULL)
        }
        aChr <- aChr[order(aChr$start_position, aChr$end_position, decreasing = FALSE),]

        nTotal <- nrow(aChr)
        windowSize <- max( c(round(nTotal * fracChrGenesInWindow), minWindowSize))
        windowSizeHalf <- round(windowSize / 2)

        windowSizeSize <- round( windowSize * windowSlideFraction )

        indexSeq <- seq(
            from = 1 + windowSizeHalf,
            to = nTotal - windowSizeHalf,
            by = windowSizeSize
        )

        indexDf <- data.frame(
            Chrom = aChr$chromosome_name[1],
            fromIndex = indexSeq - windowSizeHalf,
            toIndex = indexSeq + windowSizeHalf
        )
        indexDf$Start <- aChr$start_position[ indexDf$fromIndex ] - 1
        indexDf$End   <- aChr$end_position  [ indexDf$toIndex ] + 1

        indexDf$Length <- indexDf$End - indexDf$Start + 1
        indexDf$nGenes <- indexDf$toIndex - indexDf$fromIndex + 1

        return(indexDf[,c('Chrom','Start','End','Length','nGenes')])
    })

    chrRegionDf <- do.call(rbind, chrRegionList)
    rownames(chrRegionDf) <- paste0('region_', 1:nrow(chrRegionDf))
    return(chrRegionDf)
}



#' Estimate a two-component Gaussian Mixture Model based on a log2(CPM/10+1) expression matrix
#'
#' This function estimates a two component Gaussian Mixture Model on log2(CPM/10+1) expression matrix.
#' @param expmat A genes X samples expression matrix of log2(CPM/10+1) scaled (single cell) RNA-seq counts.
#' @param chr Number of the chromosome to be considered.
#' @param normFactor A vector of normalization factors for each cell calculated with calcNormFactors()
#' @param groups1 A vector of indices for normal cells
#' @param groups2 A vector of indices for tumor cells
#' @param gene_positions A matrix of positions for every gene in the expression matrix, calculated via getGenePositions()
#' @param start Optional: The start position of the region to be considered
#' @param end Optional: The end position of the region to be considered
#' @keywords Detect varibale genes
#' @export
#' @examples
#' detectVarGenes(suva_exp,500)

plotChrEnichment = function(expmat,chr,normFactor,gene_positions,n=1,groups1=NULL,groups2=NULL,start=NULL,end=NULL,k=2,vis=T,postProb=0.95,repetitions=5){
    if( vis ) {
        par(mfrow=c(2,2))
    }
    if (!is.null(groups1)){
        cellcolor=rep("black",(length(groups1)+length(groups2)));cellcolor[groups2]="red"
    } else{
        cellcolor=NULL
    }
    if (!is.null(start)){
        chr_genes=gene_positions[which(gene_positions[,3]==chr & gene_positions[,4]>start &  gene_positions[,5]<end),2]
    } else{
        chr_genes=gene_positions[which(gene_positions[,3]==chr),2]
    }
    if(length(chr_genes)>100){
        ### Scale data
        chr_exp=scale(colMeans(expmat[intersect(chr_genes,row.names(expmat)),])-normFactor)

        ### Initialize vectors to store best result
        bestlog=(-Inf)
        bestmix=NULL
        loglik=NULL

        ### Fit the bimodular distribution "repetitions" times retain the best
        for (i in 1:repetitions){
            #print(paste("Fitting GMM for chr",chr," ",start,":",end," iteration ",i,sep="")) # outcommented KVS

            ### Add capture.output to prevent printing
            tmp <- capture.output(
                mixmdl <- tryCatch(
                    mixtools::normalmixEM(chr_exp,k=k,maxit = 1000,maxrestarts=10),
                    error=function(e) {
                        print(paste("EM algorithm did not converge for region",chr," ",start," ",end))
                        mixmdl=NULL
                    }
                )
            )

            if(!is.null(mixmdl)){
                if (mixmdl$loglik>bestlog){
                    bestlog=mixmdl$loglik
                    bestmix=mixmdl
                }
            }
        }

        if(is.null(bestmix)){
            hist(chr_exp,breaks=50,main=paste("Chr: ",chr,":",start,":",end,"\n","Unable to fit 2 component mixture model",sep=""),lwd2=3,xlab2="Expression z-score")
            plot( runif(length(chr_exp), 0,100),chr_exp,pch=16,ylab="Expression z-score",ylim=c(min(chr_exp),(max(chr_exp)+2)),xlab="Cells")
            hist(chr_exp,breaks=50,main=paste("Chr: ",chr,":",start,":",end,"\n","Unable to fit 2 component mixture model",sep=""),lwd2=3,xlab2="Expression z-score")
            plot( runif(length(chr_exp), 0,100),chr_exp,pch=16,ylab="Expression z-score",ylim=c(min(chr_exp),(max(chr_exp)+2)),xlab="Cells")
        }
        else{
            out1 = list(
                x=chr_exp,
                mu=mean(chr_exp),
                sigma=sd(chr_exp),
                lambda=1,
                loglik=sum(dnorm(chr_exp,mean(chr_exp),sd(chr_exp),log=TRUE))
            )
            bics = c(max(BIC.mix(out1),1),max(BIC.mix(bestmix),1))
            lrt = round( likelihoodRatioTest(out1$loglik,bestmix$loglik,n),6)
            bestmix$BIC=bics
            bestmix$lrt=lrt

            ### Visualize
            if (vis==T){
                plot(bestmix,which=2,breaks=50,col1=c("red","green"),main2=paste("Chr: ",chr,":",start,":",end,"\n","Log likelihood ",round(bestmix$loglik,1),sep=""),lwd2=3,xlab2="Expression z-score")
            }
            if (length(cellcolor)>1 & vis==T){
                g1=length(which(bestmix$posterior[groups1,1]>postProb))/length(groups1)*100
                g2=length(which(bestmix$posterior[groups1,2]>postProb))/length(groups1)*100
                g3=length(which(bestmix$posterior[groups1,2]<postProb & bestmix$posterior[groups1,1]<postProb))/length(groups1)*100
                g4=length(which(bestmix$posterior[groups2,1]>postProb))/length(groups2)*100
                g5=length(which(bestmix$posterior[groups2,2]>postProb))/length(groups2)*100
                g6=length(which(bestmix$posterior[groups2,2]<postProb & bestmix$posterior[groups2,1]<postProb))/(length(groups2)+length(groups1))*100
                barplot(rbind(c(g1,g2,g3),c(g4,g5,g6)),ylim=c(0,100),beside=T,ylab="Percentage of cells",names=c("Cluster","Cluster","Ambigu"),legend = c("Non-malignant", "Malignant"),args.legend = list(title = "Pred. via transcript.", x = "topright", cex = .65),xlab="Predicted via transcriptomics")
                axis(1, at=c(0.5,1,2,3,3.3), line=2, tick=T, labels=rep("",5), lwd=3, lwd.ticks=0,col="red")
                axis(1, at=c(3.5,4,5,6,6.5), line=2, tick=T, labels=rep("",5), lwd=3, lwd.ticks=0,col="green")
                barplot(bics,names=c("1","2"),ylab="BIC",pch=16,xlab="Number of components",log="y")
                plot( runif(length(chr_exp), 0,100),chr_exp,pch=16,col=cellcolor,ylab="Expression z-score",ylim=c(min(chr_exp),(max(chr_exp)+2)),xlab="Cells")
                legend("topright", col=c("black","red"), c("Non-malignant","Malignant"), bty="o",  box.col="darkgreen", cex=.65,pch=16,title="Pred. via transcript.")
            } else{
                if (vis==T){
                    plot( runif(length(chr_exp), 0,100),chr_exp,pch=16,ylab="Expression z-score",ylim=c(min(chr_exp),(max(chr_exp)+2)),xlab="Cells")
                    barplot(bics,names=c("1","2"),ylab="BIC",pch=16,xlab="Number of components",log="y")
                    hist(bestmix$posterior[,1],main="Posterior probablility distribution\n component 1",xlab="Posterior probability",breaks=20,xlim=c(0,1))
                }
            }
            return(bestmix)
        }
    }
}

#' Estimate and visualize a two-component Gaussian Mixture Model based on a log2(CPM/10+1) expression matrix
#'
#' This function estimates a two component Gaussian Mixture Model on log2(CPM/10+1) expression matrix and visualizes the results.
#' @param mat A genes X samples expression matrix of log2(CPM/10+1) scaled (single cell) RNA-seq counts.
#' @param normFactor A vector of normalization factors for each cell calculated with calcNormFactors()
#' @param regions A matrix of regions to be evaluated
#' @param gene_pos A matrix of positions for every gene is the expression matrix, calculated via getGenePositions()
#' @param fname A string indicating the name of the output pdf
#' @param normal Optional: A vector of indices for normal cells
#' @param tumor Optional: A vector of indices for tumor cells
#' @param vis A logici indicating whether to produce visuals or not
#' @param returnInternalData A logic indicating whether to return more of the data internally generated instead of the default (posterior data). If true a list with a 1) The scaled avg exp for each region of interest matrix (CNV estimate), 2) The result of biomodual test (extention of orignal table) and 3) the classification of each sample into group1 or group2 (which are compared in the in the result of biomodual test).
#' @keywords Detect varibale genes
#' @export
#' @examples
#' detectVarGenes(suva_exp,500)

plotAll = function (mat,normFactor,regions,gene_pos,fname,normal=NULL,tumor=NULL,postProb=0.8,repetitions=4, vis=T, showProgress=TRUE, returnInternalData=FALSE){
    if( vis ) {
        pdf(paste(fname,"_CNVs.pdf",sep=""))
    }

    loglik=c()
    bic=c()
    lrt=c()
    l=c()
    testDiff <- list() # KVS added
    zScores <- list()  # KVS added
    cellClass <- list()

    if( showProgress ) {
        pb <- txtProgressBar(min = 1, max = nrow(regions), style = 3)
    }
    for(i in 1:nrow(regions)){ # i <- 2
        mixmdl=plotChrEnichment(mat,regions[i,1],normFactor,gene_pos,nrow(regions),normal,tumor,regions[i,2],regions[i,3],postProb=postProb,repetitions=repetitions,vis=vis)
        if (!is.null(mixmdl)){
            loglik=c(loglik,mixmdl$loglik)
            bic=c(bic,mixmdl$BIC)
            lrt=c(lrt,mixmdl$lrt)
            names(loglik)[length(loglik)]=rownames(regions)[i]
            names(bic)[length(bic)]=paste(rownames(regions)[i],"1_comp",sep="_")
            names(bic)[length(bic)-1]=paste(rownames(regions)[i],"2_comp",sep="_")

            if( ! returnInternalData ) { # KVS added If (else it is to timeconsuming)
                if (mixmdl$mu[1]>mixmdl$mu[2]){
                    r=mixmdl$posterior[,1]
                } else{
                    r=mixmdl$posterior[,2]
                }
                l=cbind(l,r)
                colnames(l)[ncol(l)]=rownames(regions)[i]
            }

            ### KVS added
            if( returnInternalData ) {
                # testing
                zList <- split(
                    mixmdl$x,
                    apply(mixmdl$posterior, 1, which.max)
                )
                if(
                    length(zList) == 2 &
                    all(sapply(zList, length) >= 5)
                ) {
                    testDiff[[ rownames(regions)[i] ]] <- testWrapper(
                        zList[[1]],
                        zList[[2]]
                    )
                } else {
                    testDiff[[ rownames(regions)[i] ]] <- data.frame(
                        n1 = length(mixmdl$x),
                        n2 = 0,
                        mean1 = mean(mixmdl$x),
                        mean2 = NA,
                        median1 = median(length(mixmdl$x)),
                        median2 = NA,
                        mean_shift = NA,
                        median_shift = NA,
                        wilcox_p_value = NA,
                        ttest_p_value  = NA,
                        stringsAsFactors=F
                    )
                }

                # z-scores
                localZ <- data.frame(zscore=mixmdl$x)
                colnames(localZ) <- rownames(regions)[i]
                zScores[[ rownames(regions)[i] ]] <-localZ

                # cell class
                localClass <- apply(mixmdl$posterior, 1, which.max)
                localCount <- table(localClass)
                localLarge <- as.integer(names(localCount)[which.max(localCount)])
                cellClass[[ rownames(regions)[i] ]] <- as.integer(localClass == localLarge)

            }

        }

        if( showProgress ) {
            setTxtProgressBar(pb, i)
        }
    }

    if( showProgress ) {
        close(pb)
    }
    if( vis ) {
        par(mfrow=c(1,1))
        barplot(sort(loglik),names=names(sort(loglik)),cex.axis=0.8,cex.names=0.7,las=2,ylab="log-likelihood")
        dev.off()
    }
    bicLRmat=matrix(ncol=4,nrow=length(loglik))
    bicLRmat[,1]=bic[seq(1,(length(bic)-1),2)]
    bicLRmat[,2]=bic[seq(2,length(bic),2)]
    bicLRmat[,3]=bicLRmat[,1]-bicLRmat[,2]
    bicLRmat[,ncol(bicLRmat)]=lrt
    colnames(bicLRmat)=c("BIC 1 component","BIC 2 components","BIC difference","LRT adj. p-val")
    if( ! returnInternalData ) {
        rownames(l)=colnames(mat)
    }

    ### KVS added
    if( returnInternalData ) {
        ### Summary statistics
        bicLRdf <- as.data.frame(bicLRmat)
        bicLRdf$region <- names(loglik)

        testDiffDf <- do.call(rbind, testDiff)
        testDiffDf$region <- names(testDiff)

        testDiffDf$ttest_q_value  <- p.adjust(testDiffDf$ttest_p_value, method = 'fdr')
        testDiffDf$wilcox_p_value <- p.adjust(testDiffDf$wilcox_p_value, method = 'fdr')

        combinedRes <- merge(
            bicLRdf,
            testDiffDf,
            by = 'region',
            all.x = T,
            all.y = T
        )

        ### z-scores
        zScoresDf <- as.data.frame(t(do.call(cbind, zScores)))
        colnames(zScoresDf) <- colnames(mat)

        ### Class
        cellClassDf <- as.data.frame(t(do.call(cbind, cellClass)))
        colnames(cellClassDf) <- colnames(mat)


    }

    if( returnInternalData ) {
        return(
            list(
                summaryStatistics = combinedRes,
                zScore = zScoresDf,
                cellClassification = cellClassDf
            )
        )
    } else {
        rownames(bicLRmat)=names(loglik)
        write.table(bicLRmat,paste(fname,"BIC_LR.txt",sep="_"),sep="\t")
        return(l)
    }
}

#' Visualize a heatmap of posterior probabilities of cells for component 1 of Gaussian Mixture Models fit for each user-defined region
#'
#' This function allows to plot a heatmap of posterior probabilities of cells for component 1 of Gaussian Mixture Models fit for each user-defined region.
#' @param pmat Matrix of posterior probabilites, as returned by plotAll().
#' @param expmat A genes X samples expression matrix of log2(CPM/10+1) scaled (single cell) RNA-seq counts.
#' @param clusters The number of clusters to cut the dendrogram into
#' @param zscoreThreshold A threshold on the row-wise z-scores for visualiztion. Default: 4
#' @param gene_pos A matrix of positions for every gene is the expression matrix, calculated via getGenePositions()
#' @param patients Optional: A vector naming the patient for each cell
#' @param celltypes Optional: A vector naming a putative celltype for each cell
#' @keywords Heatmap of posteriors
#' @export
#' @examples
#' plotHistogram(pmat,suva_expr,2,zscoreThreshold=4)

plotHistogram = function (pmat,expmat,clusters,zscoreThreshold=4,patients=NULL,celltypes=NULL){
  t=zscoreThreshold
  pmat=scale(pmat)
  if (max(pmat)>t){
    pmat[which(pmat>t)]=t
    pmat[which(pmat<(-t))]=(-t)
  }
  else {
    mx=min(max(pmat),abs(min(pmat)))
    sc=t/mx
    pmat=pmat*sc
    pmat[which(pmat>t)]=t
    pmat[which(pmat<(-t))]=(-t)
  }
  if(is.null(patients) & is.null(celltypes)){
    p=pheatmap::pheatmap(t(pmat),cluster_rows=F, cutree_cols = clusters, col=squash::bluered(100),gaps_col=50,show_colnames = F,clustering_distance_cols="euclidean")
  }
  else if (is.null(celltypes)){
    patientcolors =data.frame(patients)
    rownames(patientcolors)=colnames(expmat)
    pmat=pmat[colnames(expmat),]
    p=pheatmap::pheatmap(t(pmat),cluster_rows=F,silent=F, cutree_cols = clusters, col=squash::bluered(100),gaps_col=50,annotation=patientcolors,show_colnames = F,clustering_distance_cols="euclidean")
  }
  else if(is.null(patients)){
    patientcolors =data.frame(celltypes)
    rownames(patientcolors)=colnames(expmat)
    pmat=pmat[colnames(expmat),]
    p=pheatmap::pheatmap(t(pmat),cluster_rows=F, cutree_cols = clusters, col=squash::bluered(100),gaps_col=50,annotation=patientcolors,show_colnames = F,clustering_distance_cols="euclidean")
  }
  else {
    patientcolors =data.frame(celltypes)
    patientcolors=cbind(patientcolors,patients)
    rownames(patientcolors)=colnames(expmat)
	pmat=pmat[colnames(expmat),]
    p=pheatmap::pheatmap(t(pmat),cluster_rows=F, cutree_cols = clusters, col=squash::bluered(100),gaps_col=50,annotation=patientcolors,show_colnames = F,clustering_distance_cols="euclidean")
  }
 ord=unique(cutree(p$tree_col, k = clusters)[p$tree_col[["order"]]])
 numb=table(cutree(p$tree_col, k = clusters))[ord]
 n=length(numb)
 grid::grid.text(expression(bold("Cluster ID \n(left to right)")),x=rep(0.92),y=c(n*0.03+0.03),gp=grid::gpar(fontsize=8, col="grey"))
 grid::grid.text(ord,x=rep(0.92,length(numb)),y=seq(n*0.03, 0.03, -0.03),gp=grid::gpar(fontsize=8, col="grey"))
 return(cutree(p$tree_col, k = clusters))
}

