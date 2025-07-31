# define hard cutoff
hard_cutoff <- 10

# 1. load data
# TSV file of integrated RNAseq data are brought into R.
rna <- data.table::fread(file="data/processed/TCGA_GTEX_integrated_membrane_target.txt", check.names=FALSE, stringsAsFactors=F)

# 2. Data analysis
# 2.1. Data preprocessing
# membrane target candidate list
membrane <- colnames(rna)[!(colnames(rna) %in% c("Sample.ID", "type", "tissue.cancer"))]

# Normal tissue
# normal tissue risk categories

high <- c("Adipose Tissue", "Blood", "Blood Vessel", "Soft tissue/Bone", "Brain", "Esophagus", "Eye", "Heart", "Liver", "Lung", "Lymphatic tissue", "Muscle", "Nerve", "Paraganglia", "Small Intestine", "White blood cell")
 
medium <- c("Adrenal Gland", "Bladder", "Colon", "Kidney", "Pancreas", "Pituitary", "Stomach", "Skin")

low <- c("Breast", "Cervix", "Fallopian Tube", "Ovary", "Salivary Gland", "Prostate", "Rectum", "Spleen", "Thymus", "Thyroid", "Uterus", "Endometrium", "Vagina")

immunoprivileged <- c("Testis")

tbd <- c("Lining of body cavities", "Bile duct", "Head and Neck region") 

normal <- rna[type=="normal"]

normal[, Risk := fifelse(tissue.cancer%in%high, "high", fifelse(tissue.cancer%in%medium, "medium", fifelse(tissue.cancer%in%low, "low", fifelse(tissue.cancer%in%immunoprivileged, "Immunoprivileged", "TBD"))))]

# 2.2. Percentage above hard cut-off
for(i in membrane){
    tmp <- rna[, c("Sample.ID", "type", "tissue.cancer", i), with = F]
    data.table::setnames(tmp, old = i, new = "TPM")
    tmp[, TPM := (2^TPM - 0.001)]
    freq <- tmp[, list(Frq_above_cutoff = (sum(TPM>=hard_cutoff)/length(TPM)*100)), by=list(type, tissue.cancer)]
    data.table::setnames(freq, old = "Frq_above_cutoff", new = i)
    
    if(exists("freq.above.cutoff1")){
        freq.above.cutoff1 <- merge(freq.above.cutoff1, freq, by = c("type", "tissue.cancer"), all = T)
    }else{
        freq.above.cutoff1 <- freq
    }
    
    rm(tmp, freq)
}

sample.size <- rna[, list(SampleSize = length(Sample.ID)), by = list(type, tissue.cancer)]
sample.size[, Risk := fifelse(tissue.cancer%in%high, "high", fifelse(tissue.cancer%in%medium, "medium", fifelse(tissue.cancer%in%low, "low", fifelse(tissue.cancer%in%immunoprivileged, "Immunoprivileged", "NA"))))]

freq.above.cutoff1 <- merge(sample.size, freq.above.cutoff1, by = c("type", "tissue.cancer"))

# 2.3. Percentage above target specific cut-off
# 2.3.1. Cut-off calculation
# function to calculate cut-off
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

# 2.3.2. Percentage above cut-off
for(i in membrane){
    tmp <- rna[, c("Sample.ID", "type", "tissue.cancer", i), with = F]
    data.table::setnames(tmp, old = i, new = "TPM")
    tmp[, TPM := (2^TPM - 0.001)]
    cutoff <- membrane.cutoff$Cutoff[membrane.cutoff$Gene == i]
    freq <- tmp[, list(Frq_above_cutoff = (sum(TPM>=cutoff)/length(TPM)*100)), by=list(type, tissue.cancer)]
    data.table::setnames(freq, old = "Frq_above_cutoff", new = i)
    
    if(exists("freq.above.cutoff2")){
        freq.above.cutoff2 <- merge(freq.above.cutoff2, freq, by = c("type", "tissue.cancer"), all = T)
    }else{
        freq.above.cutoff2 <- freq
    }
    
    rm(tmp, cutoff, freq)
}

freq.above.cutoff2 <- merge(sample.size, freq.above.cutoff2, by = c("type", "tissue.cancer"))

tmp <- data.table::melt(freq.above.cutoff2, id.var = c("tissue.cancer"), variable.name = "Gene", value.name = "VALUE")
tmp <- data.table::dcast(tmp, Gene~tissue.cancer, value.var = "VALUE")
tmp <- merge(membrane.cutoff, tmp, by = "Gene", all = T)

tmp <- data.table::as.data.table(tmp)
tmp[, Gene := factor(Gene, levels = colnames(freq.above.cutoff))]
data.table::setorder(tmp, Gene)

# 3. Result
filename <- paste0("processed/freq_above_hard_cutoff.csv")
if(!exists(filename)){
    data.table::fwrite(freq.above.cutoff1, file=filename)
}

filename <- paste0("freq_above_target_specific_cutoff.csv")
if(!exists(filename)){
    data.table::fwrite(tmp, file=filename)
}
