#' @title Two sample t-test with data normalized by RSimNorm.
#'
#' @description This function conduct differential abundant test by t-test with data normalized by RSimNorm.
#'
#' @importFrom phyloseq tax_table otu_table sample_data phyloseq
#' @importFrom microbiome abundances aggregate_taxa
#'
#' @param physeq A \code{phyloseq} object which consists of a count table, a sample metadata, a taxonomy table.
#' The row names of the meta data must match the sample names of the count table.
#' @param count_table A numeric matrix or data frame. Row represents taxa. Column represents samples.
#' @param taxa_level Taxonomy level to conduct the test.
#' @param meta Meta table.
#' @param main_var The name of the two-level variable of interest.
#' @param method Correction method for multiple testing effect.
#' @param alpha Significance level.
#' @param eta The classification error control level. We recommend to choose it between 0 and 0.07.
#' @param lib_cut A numerical threshold for filtering samples based on library sizes. Samples with
#' library sizes less than \code{lib_cut} will be excluded in the analysis. Default is 0, i.e. do not
#' discard any sample.
#' @param bootstrap_num Bootstrap times. Default setting is 3 if the number of taxa is large.
#'
#' @returns a \code{dataframe} with columns:
#' \itemize{
#' \item{ \code{taxa}, a vector. Names of taxa.}
#' \item{ \code{p-value}, a vector. p-value.}
#' \item{ \code{adjusted p-value}, a vector. adjusted p-value.}
#' \item{ \code{is.DA}, a vector. T is differential abundant, F is not differential abundant.}
#' }
#'
#' @export
t_normalized = function(physeq = NULL, count_table = NULL, tax_level = NULL, meta = NULL, main_var,
                        method = "BH", alpha = 0.05, eta = 0, lib_cut = 0, bootstrap_num = 3){
    # 1. data preprocessing
    if(is.null(count_table)){
        if(!is.null(physeq)){
            count_table = microbiome::abundances(physeq)
            meta = phyloseq::sample_data(physeq)
            X = data_preprocess(count_table, lib_cut)
        }else{
            stop('Must have count data.')
        }
    }else{
      X = data_preprocess(count_table, lib_cut)
      meta = meta
    }
    meta = data.frame(meta)
    if(any(is.na(X))) {
        stop('The OTU/ASV table contains NAs! Please remove!\n')
      }
    X.cn = RSimNorm(count_table = X, eta, lib_cut, bootstrap_num)$P
    if(!is.null(tax_level)){
        TAXA = phyloseq::tax_table(physeq)
        META = phyloseq::sample_data(physeq)
        OTU = phyloseq::otu_table(X.cn,taxa_are_rows = T)
        physeq1 = phyloseq::phyloseq(OTU,TAXA,META)
        physeq.agg = microbiome::aggregate_taxa(physeq1,tax_level)
        X = microbiome::abundances(physeq.agg)
    }else{
        X = X.cn
    }
    Y = meta[,main_var]
    if(length(unique(Y)) > 2) {
        stop('Main varaible contains too many categories. Only two are allowed for t-test.')
    }
    if(any(is.na(Y))){
        stop('Main varaible contains NAs! Please remove!\n')
    }
    if(ncol(X)!=length(Y)){
        stop('Main varaible has different size with samples.')
    }
    d = nrow(X)
    data = data.frame(cbind(t(X),Y))
    colnames(data)[ncol(data)] = 'Y'
    t_p <- sapply(1:d, function(i)t.test(as.numeric(data[, i]) ~ Y, data = data)$p.value)
    p <- p.adjust(t_p, method=method)
    res = data.frame('taxa' = rownames(X), 'p-value' = t_p, 'adjusted p-value' = p, 'is.DA' = (p<alpha))
    return(res)
}
