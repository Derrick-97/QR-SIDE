#Factor Analysis for Princomp
factpc<-function(X, m=2,rotation="none",scores="regression")
{  
  options(digits=4)
  S=cor(X); 
  p<-nrow(S); diag_S<-diag(S); sum_rank<-sum(diag_S)
  #rowname<-paste("X", 1:p, sep="")
  rowname = names(X)
  colname<-paste("Factor", 1:p, sep="")
  A<-matrix(0, nrow=p, ncol=p, dimnames=list(rowname, colname))
  eig<-eigen(S)
  for (i in 1:p){
    if (min(eig$values[i]) < 0){
      A[,i]<-sqrt(eig$values[i] - min(eig$values[i]))*eig$vectors[,i] ## make positive
    }else{
      A[,i]<-sqrt(eig$values[i])*eig$vectors[,i] ## make positive
    }
  }

  for (i in 1:p) { if(sum(A[,i])<0) A[,i] = -A[,i] }
  h<-diag(A%*%t(A))
  rowname<-c("SS loadings", "Proportion Var", "Cumulative Var")
  B<-matrix(0, nrow=3, ncol=p, dimnames=list(rowname, colname))
  for (i in 1:p){
    B[1,i]<-sum(A[,i]^2)
    B[2,i]<-B[1,i]/sum_rank
    B[3,i]<-sum(B[1,1:i])/sum_rank
  }
  #cat("\nFactor Analysis for Princomp: \n\n");
  #cat("\n"); print(B[,1:m]);
  W=B[2,1:m]*100; 
  Vars=cbind('Vars'=B[1,],'Vars.Prop'=B[2,],'Vars.Cum'=B[3,]*100)
  #cat("\n"); print(Vars[1:m,])
  #cat("\n"); print(A[,1:m]);
  A=A[,1:m] 
  #{ cat("\n common and specific \n"); print(cbind(common=h, spcific=diag_S-h)); }
  if(rotation=="varimax")
  {   
    #stop(" factor number >= 2 !")
    cat("\n Factor Analysis for Princomp in Varimax: \n\n");
    VA=varimax(A); A=VA$loadings; 
    s2=apply(A^2,2,sum); 
    k=rank(-s2); s2=s2[k]; 
    W=s2/sum(B[1,])*100; 
    Vars=cbind('Vars'=s2,'Vars.Prop'=W,'Vars.Cum'=cumsum(W))
    rownames(Vars) <- paste("Factor", 1:m, sep="")
    #print(Vars[1:m,])
    A=A[,k]
    for (i in 1:m) { if(sum(A[,i])<0) A[,i] = -A[,i] }
    A=A[,1:m]; 
    colnames(A) <- paste("Factor", 1:m, sep="")
    #cat("\n"); print(A) 
  }
  fit<-NULL
  fit$Vars<-Vars[1:m,]
  fit$loadings <- A
  X=as.matrix(scale(X));
  PCs=X%*%solve(S + diag(0.000001, p, p))%*%A
  #if(scores) cat("\n"); print(PCs)
  fit$scores <- PCs
  #if(rank)
  { 
    W=W/sum(W);
    PC=PCs%*%W;
    #cat("\n"); print(cbind(PCs,'PC'=PC[,1],'rank'=rank(-PC[,1])))
    Ri=data.frame('F'=PC,'Ri'=rank(-PC))
    fit$Rank <- Ri
  }
  #if(plot)
  #{ plot(PCs);abline(h=0,v=0,lty=3); text(PCs,label=rownames(PCs),pos=1.1,adj=0.5,cex=0.85) }
  #if(biplot)
  #{ biplot(PCs,A) } #text(PCs,label=rownames(PCs),pos=1.1,adj=0.5,cex=0.85) 
  common=apply(A^2,1,sum);
  fit$common <- common
  fit
  #list(Vars=B[,1:m],loadings=A,scores=PCs,Ri=Ri,common=common)
} #fa=factpc(X,2)


Createsce <- function(X, n){
  library(SingleCellExperiment)
  # -------------------------------------------------
  # make BayesSpace metadata used in BayesSpace
  counts <- t(X)
  ## Make array coordinates - filled rectangle
  cdata <- list()
  cdata$row <- rep(1, n)
  cdata$col <- rep(1, n)
  cdata <- as.data.frame(do.call(cbind, cdata))
  ## Scale and jitter image coordinates
  #scale.factor <- rnorm(1, 8);  n_spots <- n
  #cdata$imagerow <- scale.factor * cdata$row + rnorm(n_spots)
  #cdata$imagecol <- scale.factor * cdata$col + rnorm(n_spots)
  cdata$imagerow <- cdata$row
  cdata$imagecol <- cdata$col 
  ## Make SCE
  ## note: scater::runPCA throws warning on our small sim data, so use prcomp
  sce <- SingleCellExperiment(assays=list(counts=counts), colData=cdata)
  #reducedDim(sce, "PCA") <- fit_y
  sce$spatial.cluster <- floor(runif(ncol(sce), 1, 3))
  
  metadata(sce)$BayesSpace.data <- list()
  metadata(sce)$BayesSpace.data$platform <- "Visium"
  metadata(sce)$BayesSpace.data$is.enhanced <- FALSE
  sce
}


spatialPreprocess <- function(sce, platform=c("Visium", "ST"),
                              n.PCs=15, n.HVGs=2000, skip.PCA=FALSE,
                              log.normalize=TRUE, assay.type="logcounts",
                              BSPARAM=ExactParam()) {
  
  ## Set BayesSpace metadata
  metadata(sce)$BayesSpace.data <- list()
  metadata(sce)$BayesSpace.data$platform <- match.arg(platform)
  metadata(sce)$BayesSpace.data$is.enhanced <- FALSE
  # metadata(sce)$BayesSpace.data$use_dimred <- use.dimred
  # metadata(sce)$BayesSpace.data$d <- n.PCs
  
  ## Run PCA on HVGs, log-normalizing if necessary
  if (!skip.PCA) {
    if (log.normalize)
      sce <- logNormCounts(sce)
    
    dec <- modelGeneVar(sce, assay.type=assay.type)
    top <- getTopHVGs(dec, n=n.HVGs)
    sce <- runPCA(sce, subset_row=top, ncomponents=n.PCs, 
                  exprs_values=assay.type, BSPARAM=BSPARAM)
    rowData(sce)[["is.HVG"]] <- (rownames(sce) %in% top)
  }
  
  sce
}



## calculate RMSE
calRMSE = function(groundtruth, deconProp, permutation = TRUE){
  library(combinat)
  K = dim(groundtruth)[2]
  n = length(permn(K))
  if (permutation == TRUE){
    ## permutation
    perm = matrix(0, n, K)
    perm_tmp=permn(K)
    for (i in 1:n){
      perm[i,] = perm_tmp[[i]]
    }
    out = matrix(0, n, 1)
    for (i in 1:n){
      out[i] = sqrt(norm(groundtruth - deconProp[,perm[i,]], "2")^2/K)
    }
    idx = which.min(out)
    output = list()
    output[[1]] = sqrt(norm(groundtruth - deconProp[,perm[idx,]], "2")^2/K)
    output[[2]] = perm[idx,]
    return(output)
  }else{
    return(sqrt(norm(groundtruth - deconProp, "2")^2/K))
  }
}



generate_groundtruth <- function(){
  n_s = c(50, 100, 100, 100, 100, 100, 100, 100, 100, 50)
  true_label = c(rep(1,50),rep(2:9, each = 100), rep(10,50))
  cum_n_s = cumsum(c(50, 100, 100, 100, 100, 100, 100, 100, 100, 50))
  n = sum(n_s)
  beta = matrix(c(1,0,0,
                  0.75,0.25,0,
                  5/8,3/8,0,
                  0.5,0.5,0,
                  0.25,0.75,0,
                  0,0.75,0.25,
                  0,0.5,0.5,
                  0,3/8,5/8,
                  0,0.25,0.75,
                  0,0,1),10,3, byrow = TRUE)
  
  groundtruth = matrix(0,900,3)
  for (i in 1:900){
    if (i < 51){
      groundtruth[i,] = beta[1,]
    }
    if (i > 50 & i < 151){
      groundtruth[i,] = beta[2,]
    }
    if (i > 150 & i < 251){
      groundtruth[i,] = beta[3,]
    }
    if (i > 250 & i < 351){
      groundtruth[i,] = beta[4,]
    }
    if (i > 350 & i < 451){
      groundtruth[i,] = beta[5,]
    }
    if (i > 450 & i < 551){
      groundtruth[i,] = beta[6,]
    }
    if (i > 550 & i < 651){
      groundtruth[i,] = beta[7,]
    }
    if (i > 650 & i < 751){
      groundtruth[i,] = beta[8,]
    }
    if (i > 750 & i < 851){
      groundtruth[i,] = beta[9,]
    }
    if (i > 850){
      groundtruth[i,] = beta[10,]
    }
  }
  return(groundtruth)
}

