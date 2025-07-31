# define cutoff
cutoff <- 10

# TSV file of integrated RNAseq data are brought into R.
rna <- data.table::fread(file="data/processed/TCGA_GTEX_integrated_membrane_target.txt", check.names=FALSE, stringsAsFactors=F)

# Data analysis
# Data preprocessing
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

# Percentage above cut-off
for(i in membrane){
    tmp <- rna[, c("Sample.ID", "type", "tissue.cancer", i), with = F]
    data.table::setnames(tmp, old = i, new = "TPM")
    tmp[, TPM := (2^TPM - 0.001)]
    freq <- tmp[, list(Frq_above_cutoff = (sum(TPM>=cutoff)/length(TPM)*100)), by=list(type, tissue.cancer)]
    data.table::setnames(freq, old = "Frq_above_cutoff", new = i)
    
    if(exists("freq.above.cutoff")){
        freq.above.cutoff <- merge(freq.above.cutoff, freq, by = c("type", "tissue.cancer"), all = T)
    }else{
        freq.above.cutoff <- freq
    }
    
    rm(tmp, freq)
}

sample.size <- rna[, list(SampleSize = length(Sample.ID)), by = list(type, tissue.cancer)]
sample.size[, Risk := fifelse(tissue.cancer%in%high, "high", fifelse(tissue.cancer%in%medium, "medium", fifelse(tissue.cancer%in%low, "low", fifelse(tissue.cancer%in%immunoprivileged, "Immunoprivileged", "NA"))))]

freq.above.cutoff <- merge(sample.size, freq.above.cutoff, by = c("type", "tissue.cancer"))

# Result
filename <- paste0("processed/freq_above_cutoff.csv")
if(!exists(filename)){
    data.table::fwrite(freq.above.cutoff, file=filename)
}
