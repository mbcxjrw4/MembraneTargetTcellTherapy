source("config/tissue_risk_categories.R")

# function to calculate percentage above hard cutoff
percentage_above_hard_cutoff <- function(hard_cutoff = 10, expression_data){
    # membrane target candidate list
    membrane <- colnames(expression_data)[!(colnames(expression_data) %in% c("Sample.ID", "type", "tissue.cancer"))]
    
    for(i in membrane){
        tmp <- expression_data[, c("Sample.ID", "type", "tissue.cancer", i), with = F]
        data.table::setnames(tmp, old = i, new = "TPM")
        tmp[, TPM := (2^TPM - 0.001)]
        freq <- tmp[, list(Frq_above_cutoff = (sum(TPM>=hard_cutoff)/length(TPM)*100)), by=list(type, tissue.cancer)]
        data.table::setnames(freq, old = "Frq_above_cutoff", new = i)
    
        if(exists("freq.above.cutoff")){
            freq.above.cutoff <- merge(freq.above.cutoff, freq, by = c("type", "tissue.cancer"), all = T)
        }else{
            freq.above.cutoff <- freq
        }
    
        rm(tmp, freq)
    }

    sample.size <- expression_data[, list(SampleSize = length(Sample.ID)), by = list(type, tissue.cancer)]
    sample.size[, Risk := fifelse(tissue.cancer%in%high, "high", fifelse(tissue.cancer%in%medium, "medium", fifelse(tissue.cancer%in%low, "low", fifelse(tissue.cancer%in%immunoprivileged, "Immunoprivileged", "NA"))))]

    freq.above.cutoff <- merge(sample.size, freq.above.cutoff, by = c("type", "tissue.cancer"))
    return(freq.above.cutoff)
}

# function to calculate abs max and relevant max
relevant_max <- function(y, Risk, hi_f, med_f, lo_f){
    if(IQR(y)==0){
        background <- quantile(y, probs = 0.75)
        y <- y[y>background]
    }
    # tissue risk factor
    Risk <- unique(Risk)
    risk.factor <- data.table::fifelse(Risk=="high", hi_f, data.table::fifelse(Risk=="medium", med_f, lo_f))
    
    if(length(y)>1){
        return((quantile(y, probs = c(0.75)) + risk.factor*IQR(y)))
    }else if(length(y)==1){
        return(y)
    }else{
        return(background)
    }
}

absolute_max <- function(y, Risk, hi_f, med_f, lo_f){
    # tissue risk factor
    Risk <- unique(Risk)
    return(data.table::fifelse(Risk=="high", quantile(y, probs = hi_f), data.table::fifelse(Risk=="medium", quantile(y, probs = med_f), quantile(y, probs = lo_f))))

}

# function to calculate percentage above target specific cutoff
percentage_above_target_specific_cutoff <- function(expression_data){
    # membrane target candidate list
    membrane <- colnames(expression_data)[!(colnames(expression_data) %in% c("Sample.ID", "type", "tissue.cancer"))]
    
    # normal tissue
    normal <- expression_data[type=="normal"]
    normal[, Risk := fifelse(tissue.cancer%in%high, "high", fifelse(tissue.cancer%in%medium, "medium", fifelse(tissue.cancer%in%low, "low", fifelse(tissue.cancer%in%immunoprivileged, "Immunoprivileged", "TBD"))))]

    # calculate target specific cut-off
    for(i in membrane){
        tmp <- normal[, c("Sample.ID", "type", "tissue.cancer", "Risk", i), with = F]
        data.table::setnames(tmp, old = i, new = "TPM")
        tmp[, TPM := (2^TPM - 0.001)]
        threshold <- tmp[, list(abs_max=absolute_max(y=TPM, Risk=Risk, hi_f = 0.975, med_f = 0.9, lo_f = 0.75), rel_max=relevant_max(y=TPM, Risk=Risk, hi_f = 1, med_f = 0.9, lo_f = 0.8)), by=list(tissue.cancer)]
        cutoff.adap <- min(quantile(threshold$rel_max, probs = 0.95), quantile(threshold$abs_max, probs = 0.95))
        # set minmal as 10 TPM
        if(cutoff.adap < 10){
            cutoff.adap <- 10
        }
    
        res <- data.frame(Gene = i, Cutoff = cutoff.adap)
        if(exists("membrane.cutoff")){
            membrane.cutoff <- rbind(membrane.cutoff, res)
        }else{
            membrane.cutoff <- res
        }
        rm(tmp, threshold, cutoff.adap, res)
    }

    # Percentage above cut-off
    for(i in membrane){
        tmp <- expression_data[, c("Sample.ID", "type", "tissue.cancer", i), with = F]
        data.table::setnames(tmp, old = i, new = "TPM")
        tmp[, TPM := (2^TPM - 0.001)]
        cutoff <- membrane.cutoff$Cutoff[membrane.cutoff$Gene == i]
        freq <- tmp[, list(Frq_above_cutoff = (sum(TPM>=cutoff)/length(TPM)*100)), by=list(type, tissue.cancer)]
        data.table::setnames(freq, old = "Frq_above_cutoff", new = i)
    
        if(exists("freq.above.cutoff")){
            freq.above.cutoff <- merge(freq.above.cutoff, freq, by = c("type", "tissue.cancer"), all = T)
        }else{
            freq.above.cutoff <- freq
        }
        rm(tmp, cutoff, freq)
    }
    
    sample.size <- expression_data[, list(SampleSize = length(Sample.ID)), by = list(type, tissue.cancer)]
    sample.size[, Risk := fifelse(tissue.cancer%in%high, "high", fifelse(tissue.cancer%in%medium, "medium", fifelse(tissue.cancer%in%low, "low", fifelse(tissue.cancer%in%immunoprivileged, "Immunoprivileged", "NA"))))]

    freq.above.cutoff <- merge(sample.size, freq.above.cutoff, by = c("type", "tissue.cancer"))

    res <- data.table::melt(freq.above.cutoff, id.var = c("tissue.cancer"), variable.name = "Gene", value.name = "VALUE")
    res <- data.table::dcast(res, Gene~tissue.cancer, value.var = "VALUE")
    res <- merge(membrane.cutoff, res, by = "Gene", all = T)
    
    res <- data.table::as.data.table(res)
    res[, Gene := factor(Gene, levels = colnames(freq.above.cutoff))]
    data.table::setorder(res, Gene)
    return(res)
}

# Percentage above max of normal
percentage_above_normal_max <- function(expression_data){
    # membrane target candidate list
    membrane <- colnames(expression_data)[!(colnames(expression_data) %in% c("Sample.ID", "type", "tissue.cancer"))]
    
    # normal tissue
    normal <- expression_data[type=="normal"]
    normal[, Risk := fifelse(tissue.cancer%in%high, "high", fifelse(tissue.cancer%in%medium, "medium", fifelse(tissue.cancer%in%low, "low", fifelse(tissue.cancer%in%immunoprivileged, "Immunoprivileged", "TBD"))))]

    # calculate max expression of normal
    for(i in membrane){
        tmp <- normal[, c("Sample.ID", "type", "tissue.cancer", "Risk", i), with = F]
        data.table::setnames(tmp, old = i, new = "TPM")
        tmp[, TPM := (2^TPM - 0.001)]
        threshold <- tmp[, list(abs_max=absolute_max(y=TPM, Risk=Risk, hi_f = 0.975, med_f = 0.9, lo_f = 0.75), rel_max=relevant_max(y=TPM, Risk=Risk, hi_f = 1, med_f = 0.9, lo_f = 0.8)), by=list(tissue.cancer)]
        cutoff.adap <- min(quantile(threshold$rel_max, probs = 0.95), quantile(threshold$abs_max, probs = 0.95))
    
        res <- data.frame(Gene = i, Cutoff = cutoff.adap)
        if(exists("membrane.cutoff")){
            membrane.cutoff <- rbind(membrane.cutoff, res)
        }else{
            membrane.cutoff <- res
        }
        rm(tmp, threshold, cutoff.adap, res)
    }

    # Percentage above normal max
    for(i in membrane){
        tmp <- expression_data[, c("Sample.ID", "type", "tissue.cancer", i), with = F]
        data.table::setnames(tmp, old = i, new = "TPM")
        tmp[, TPM := (2^TPM - 0.001)]
        cutoff <- membrane.cutoff$Cutoff[membrane.cutoff$Gene == i]
        freq <- tmp[, list(Frq_above_cutoff = (sum(TPM>=cutoff)/length(TPM)*100)), by=list(type, tissue.cancer)]
        data.table::setnames(freq, old = "Frq_above_cutoff", new = i)
    
        if(exists("freq.above.normal.max")){
            freq.above.normal.max <- merge(freq.above.normal.max, freq, by = c("type", "tissue.cancer"), all = T)
        }else{
            freq.above.normal.max <- freq
        }
        rm(tmp, cutoff, freq)
    }
    
    sample.size <- expression_data[, list(SampleSize = length(Sample.ID)), by = list(type, tissue.cancer)]
    sample.size[, Risk := fifelse(tissue.cancer%in%high, "high", fifelse(tissue.cancer%in%medium, "medium", fifelse(tissue.cancer%in%low, "low", fifelse(tissue.cancer%in%immunoprivileged, "Immunoprivileged", "NA"))))]

    freq.above.normal.max <- merge(sample.size, freq.above.normal.max, by = c("type", "tissue.cancer"))

    res <- data.table::melt(freq.above.normal.max, id.var = c("tissue.cancer"), variable.name = "Gene", value.name = "VALUE")
    res <- data.table::dcast(res, Gene~tissue.cancer, value.var = "VALUE")
    res <- merge(membrane.cutoff, res, by = "Gene", all = T)

    res <- data.table::as.data.table(tmp2)
    res[, Gene := factor(Gene, levels = colnames(freq.above.normal.max))]
    data.table::setorder(res, Gene)

    return(res)
}
