data_preprocess = function(feature_table, lib_cut){
    # Discard samples with library size < lib_cut
    lib_size = colSums(feature_table, na.rm = T)
    samp_keep = which(lib_size >= lib_cut)
    feature_table = feature_table[,samp_keep]
    
    # Discard taxa with total count = 0
    feature_table = feature_table[rowSums(feature_table)>0,]
    
    core_data = as.matrix(feature_table)
    return(core_data)
}
