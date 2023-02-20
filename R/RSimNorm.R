#' @title A scaling normalization method based on rank correlation.
#'
#' @description A method for efficient normalization of RNA-seq data. This method is designed to eliminate the bias from sequencing depth and reveal the true biological variations. It can be directly applied to PCoA, two-sample test and differential abundant test. 
#'
#' @details We present a new normalization method on OTUs/ASVs which is robust to zero prevalence
#' and has a theoretical guarantee of the accuracy of normalization.
#'
#' @importFrom matrixStats colMedians
#' @importFrom microbiome abundances
#'
#' @examples
#'
#'  data(gut_cn)
#'      Datanorm <- RSimNorm(count_table = A)$P
#'
#' @param physeq A \code{phyloseq} object which consists of a count table, a sample metadata. The row
#' names of the meta data must match the sample names of the count table.
#' @param count_table A numeric matrix or data frame. Row represents taxa. Column represents samples.
#' @param eta The classification error control level. We recommend to choose it between 0 and 0.07.
#' @param lib_cut A numerical threshold for filtering samples based on library sizes. Samples with
#' library sizes less than \code{lib_cut} will be excluded in the analysis. Default is 0, i.e. do not
#' discard any sample.
#' @param bootstrap_num Bootstrap times. Default setting is 3 if the number of taxa is large.
#'
#' return a \code{list} with components:
#'\itemize{
#'\item{ \code{P}, a matrix. Normalized data. Row represents taxa. Column represents samples.}
#'\item{ \code{I0}, a vector. Selected reference set.}
#'\item{ \code{pi0}, a number. The estimated proportion of reference set.}
#'}
RSimNorm = function(physeq = NULL, count_table = NULL, eta=0.01, lib_cut = 0, bootstrap_num = 3){
    # 1. data preprocessing
    if(is.null(count_table)){
        if(!is.null(physeq)){
            count_table = microbiome::abundances(physeq)
            X = data_preprocess(count_table, lib_cut)
        }else{
            stop('Must have count data.')
        }
    }else{
      X = data_preprocess(count_table, lib_cut)
    }
    if(any(is.na(X))) {
        stop('The OTU/ASV table contains NAs! Please remove!\n')
      }
    d = nrow(X)
    if(d < 1000){
      warn_txt = sprintf(paste("Results might be invalid when the number of taxa is too small. Please conduct the normalization on a OTUs/ASVs level."))
      warning(warn_txt, call. = FALSE)
    }
    
    # 2. Calculate M statistic
    v = CStat(X)
    I0.1 = which(v>0.8)
    X0 = X[I0.1,]
    if(bootstrap_num*0.5*nrow(X0)<100){
        bootstrap_num = ceiling(100/(0.5*nrow(X0)))
    }
    v0 = replicate(bootstrap_num,CStat(X0[sample(1:nrow(X0),0.5*nrow(X0)),]))
    w = v[v>0.8]
    if(is.null(w)){
        stop('Please conduct the normalization on a OTUs/ASVs level without prevalence filtering.')
    }
    f1 = sapply(w,function(x)mean(v>x))
    f0 = sapply(w,function(x)mean(v0>x))
    pi = sum(f1*f0)/sum(f0^2)
    vord = order(v,decreasing = T)
    res = sapply(1:length(vord),function(x)(1-pi*length(vord)*mean(v0>v[x])/(which(vord==x))))
    lowerx = max(which(res[vord]<eta))
    ref = vord[1:lowerx]
    if(length(ref)<5){
        stop('Please choose a larger eta for a valid normalization results.')
    }
    tc.cn = apply(X,2,function(x)sum(x[ref]))
    f.cn = tc.cn/(mean(tc.cn))
    f.cn = ifelse(f.cn==0,1,f.cn)
    cn.res = scale(X,center=FALSE,scale=f.cn)
    return(list('P' = cn.res, 'I0' = ref, 'pi0'= pi))
}

# calculate the median of spearman correlation between one taxon with others.
CStat = function(X){
  d <- nrow(X)
  R <- X
  S1 <- apply(R,1,order)
  S1 <- S1 - colMeans(S1);
  S1 <- S1 / sqrt(colSums(S1^2));
  corr_s <- crossprod(S1)
  med <- as.data.frame(matrixStats::colMedians(corr_s))
  return(as.numeric(med[,1]))
}
