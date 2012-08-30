#source("simctest/R/simctest.r");

# mmctest -- Multiple Monte-Carlo Tests





# export: mmctest, hBH, hPC, hIndep
# some FDR controlling procedures

# FDR control by Benjamini-Hochberg
hBH <- function(p, threshold) {

    return(rank(p) <= max( c(which(sort(p)<=(1:length(p))*threshold/length(p)),-1) ));
}

# FDR control by Pounds&Cheng
hPC <- function(p, threshold) {

  m <- length(p);
  pi0_hat <- min(1, 2/m*sum(p));
  return( rank(p) <= max( c(which(sort(p)<=(1:length(p))*threshold/ (pi0_hat*m) ),-1) ));
}

# independent testing via Bonferroni
hBonferroni <- function(p, threshold) {

  m <- length(p);
  return(p<=threshold/m);
}





# class mmctSamplerGeneric
# derive a class from mmctSamplerGeneric to implement the interface used to draw new samples
# or use class mmctSampler to directly pass a function and the number of hypotheses
setClass("mmctSamplerGeneric")
setGeneric("getSamples", def=function(obj, ind, n){standardGeneric("getSamples")})
setGeneric("getNumber", def=function(obj){standardGeneric("getNumber")})





# class mmctSampler
# directly pass a function "f" and the number of hypotheses "n", class has a slot for additional data
setClass("mmctSampler", contains="mmctSamplerGeneric", representation=representation(f="function", num="numeric", data="numeric"))

# implemented method
# for each hypotheses in vector ind[i], request n[i] new samples
setMethod("getSamples", signature(obj="mmctSampler"), function(obj, ind, n) {
    return(obj@f(ind, n, obj@data));
  }
)

# implemented method
setMethod("getNumber", signature(obj="mmctSampler"), function(obj) {
    return(obj@num);
  }
)

# pseudo constructor
mmctSampler <- function(f, num, data=NULL) {
  obj <- new("mmctSampler");
  obj@f <- f;
  obj@num <- num;
  obj@data <- data;
  return(obj);
}





# class mmctestres
# exportMethods: cont, show, pEstimate
setClass("mmctestres", contains="sampalgPrecomp", representation=representation(internal="environment", epsilon="numeric", threshold="numeric", h="function",
gensample="mmctSamplerGeneric", g="numeric", num="numeric", A="numeric", B="numeric", C="numeric"))

setGeneric("mainalg", def=function(obj, stopcrit){standardGeneric("mainalg")})
setMethod("mainalg", signature(obj="mmctestres"), function(obj, stopcrit) {

    m <- getNumber(obj@gensample);
    pu <- rep(0, m);
    pl <- rep(0, m);
    k <- 1000;
    it <- 0;
    timer <- proc.time()[[3]];

    g <- obj@g;
    num <- obj@num;
    batchN <- obj@internal$batchN;
    currentBatch <- max(batchN);

    A <- obj@A;
    B <- obj@B;
    C <- obj@C;
    copying <- 0;

    tryCatch({
      while(length(B)>stopcrit$undecided) {

	# sample all \hat{p}_i in B
	currentBatch <- floor(1.25*currentBatch);
	batchN[B] <- currentBatch;
	g[B] <- g[B] + getSamples(obj@gensample, B, currentBatch);
	num[B] <- num[B] + batchN[B];

	# compute upper "exact" confidence level (Clopper-Pearson)
	a <- (num/(num+k) - (num-batchN)/(num-batchN+k)) * (1-(1-obj@epsilon)**(1/m));
	qindex <- (g>0) & (g<num);
	pu[qindex] <- 1 - qbeta(a[qindex]/2, num[qindex]-g[qindex], g[qindex]+1);
	pu[g==0] <- 1-(a[g==0]/2)**(1/num[g==0]);
	pu[g==num] <- 1;

	# ...and lower "exact" confidence level (Clopper-Pearson)
	pl[qindex] <- 1 - qbeta(1-a[qindex]/2, num[qindex]+1-g[qindex], g[qindex]);
	pl[g==0] <- 0;
	pl[g==num] <- (a[g==num]/2)**(1/num[g==num]);

	A <- which(obj@h(pu, obj@threshold));
	C <- which(obj@h(pl, obj@threshold));
	B <- setdiff(C, A);

	copying <- 1;
	obj@g <- g;
	obj@num <- num;
	obj@internal$batchN <- batchN;

	obj@A <- A;
	obj@B <- B;
	obj@C <- C;
	copying <- 0;

	it <- it+1;
	if((stopcrit$maxnum>0) && (sum(num)+length(B)*floor(1.25*currentBatch)>=stopcrit$maxnum)) { break; }
	if((stopcrit$maxit>0) && (it>=stopcrit$maxit)) { break; }
	if((stopcrit$elapsedsec>0) && (proc.time()[[3]]-timer>=stopcrit$elapsedsec)) { break; }
      } # end while
    }, interrupt = function(interrupt){
      # finish copying if appropriate
      if(copying) {
	obj@g <- g;
	obj@num <- num;
	obj@internal$batchN <- batchN;

	obj@A <- A;
	obj@B <- B;
	obj@C <- C;
	return(obj);
      }
    }) # end tryCatch

    return(obj);
  }
)

# cont offers four stopping criteria given as list "steps"
# steps=list(maxit=0, maxnum=0, undecided=0, elapsedsec=0)
# corresonding to a maximal number of iterations, maximal number of samples drawn, number of yet undecided hypotheses, elapsed time
setMethod("cont", signature(data="mmctestres"), function(data, steps=list(maxit=0, maxnum=0, undecided=0, elapsedsec=0)) {

    if(length(setdiff(ls(steps),c("maxit","maxnum","undecided","elapsedsec")))>0) { stop("parameter steps contains invalid data"); }

    if(length(steps)==0) { stopcrit <- list(maxit=0, maxnum=0, undecided=0, elapsedsec=0); }
    else {
      stopcrit <- list();
      if("maxit" %in% ls(steps)) { stopcrit<-c(stopcrit,maxit=steps$maxit); } else { stopcrit<-c(stopcrit,maxit=0); }
      if("maxnum" %in% ls(steps)) { stopcrit<-c(stopcrit,maxnum=steps$maxnum); } else { stopcrit<-c(stopcrit,maxnum=0); }
      if("undecided" %in% ls(steps)) { stopcrit<-c(stopcrit,undecided=steps$undecided); } else { stopcrit<-c(stopcrit,undecided=0); }
      if("elapsedsec" %in% ls(steps)) { stopcrit<-c(stopcrit,elapsedsec=steps$elapsedsec); } else { stopcrit<-c(stopcrit,elapsedsec=0); }
    }

    return(mainalg(data, stopcrit));
  }
)

setMethod("show", signature(object="mmctestres"), function(object) {

    m <- getNumber(object@gensample);
    cat(paste("Number of rejected hypotheses: ",length(object@A),"\n",sep=""));
    cat(paste("Number of non-rejected hypotheses: ",length(setdiff(1:m, object@C)),"\n",sep=""));
    cat(paste("Number of unclassified hypotheses: ",length(object@B),"\n",sep=""));
    cat(paste("Total number of samples: ",sum(object@num),"\n",sep=""));
  }
)

setGeneric("pEstimate", def=function(obj){standardGeneric("pEstimate")})
setMethod("pEstimate", signature(obj="mmctestres"), function(obj) {

    return(obj@g/obj@num);
  }
)

setGeneric("confidenceLimits", def=function(obj){standardGeneric("confidenceLimits")})
setMethod("confidenceLimits", signature(obj="mmctestres"), function(obj) {

    m <- getNumber(obj@gensample);
    k <- 1000;
    g <- obj@g;
    num <- obj@num;
    batchN <- obj@internal$batchN;
    pu <- rep(0, m);
    pl <- rep(0, m);

    # compute upper "exact" confidence level (Clopper-Pearson)
    a <- (num/(num+k) - (num-batchN)/(num-batchN+k)) * (1-(1-obj@epsilon)**(1/m));
    qindex <- (g>0) & (g<num);
    pu[qindex] <- 1 - qbeta(a[qindex]/2, num[qindex]-g[qindex], g[qindex]+1);
    pu[g==0] <- 1-(a[g==0]/2)**(1/num[g==0]);
    pu[g==num] <- 1;

    # ...and lower "exact" confidence level (Clopper-Pearson)
    pl[qindex] <- 1 - qbeta(1-a[qindex]/2, num[qindex]+1-g[qindex], g[qindex]);
    pl[g==0] <- 0;
    pl[g==num] <- (a[g==num]/2)**(1/num[g==num]);

    l <- list();
    l$lowerLimits <- pl;
    l$upperLimits <- pu;
    return(l);
  }
)

setGeneric("testResult", def=function(obj){standardGeneric("testResult")})
setMethod("testResult", signature(obj="mmctestres"), function(obj) {

    m <- getNumber(obj@gensample);
    l <- list();
    l$rejected <- obj@A;
    l$nonrejected <- setdiff(1:m, obj@C);
    l$undecided <- obj@B;
    return(l);
  }
)

setGeneric("summary.mmctestres", def=function(object,...){standardGeneric("summary.mmctestres")})
setMethod("summary.mmctestres", signature(object="mmctestres"), function(object) {

    cat(paste("Number of hypotheses: ",getNumber(object@gensample),sep=""));
    cat(strwrap(paste("Indices of rejected hypotheses:", paste(object@A,collapse=" ")),prefix="\n"));
    cat(strwrap(paste("Indices of unclassified hypotheses:", paste(object@B,collapse=" ")),prefix="\n"));
    cat("\nAll hypotheses not listed are classified as not rejected.\n");
  }
)





# class mmctest
# exportMethods: run
setClass("mmctest", contains="mmctestres", representation=representation(internal="environment"))

setMethod("run", signature(alg="mmctest", gensample="mmctSamplerGeneric"), function(alg, gensample, maxsteps=list(maxit=0, maxnum=0, undecided=0, elapsedsec=0)) {

    obj <- new("mmctestres");

    obj@epsilon <- alg@internal$epsilon;
    obj@threshold <- alg@internal$threshold;
    obj@h <- alg@internal$h;
    obj@gensample <- gensample;

    m <- getNumber(obj@gensample);
    obj@g <- rep(0, m);
    obj@num <- rep(0, m);
    obj@internal$batchN <- rep(10, m);

    obj@A <- 0;
    obj@B <- 1:m;
    obj@C <- 0;

    return(cont(obj, steps=maxsteps));
  }
)

# pseudo constructor
mmctest <- function(epsilon=0.01, threshold=0.1, h) {

  obj <- new("mmctest");
  obj@internal$epsilon=epsilon;
  obj@internal$threshold=threshold;
  obj@internal$h=h;
  return(obj);
}
