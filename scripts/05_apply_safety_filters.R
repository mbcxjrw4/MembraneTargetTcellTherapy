# To show potential membrane target safety bucket
# Various sub-categories of good safety profile
# Expression restricted to gastro-intestinal tract (GIT) / epithelial tissues – scientist should further validate tissue polarity including the literature
# Various sub-categories of therapeutic window – manually check >10 TPM and target threshold graphs to deduce likely therapeutic window
# Various sub-categories for predicted unsafe, often after calculation of FC between normal and tumour

# 1. load data
# CSV file of % prevalence of potential target expression above 10 TPM in tumours and normal tissues are brought into R.
prevalence1 <- data.table::fread(file="data/processed/freq_above_hard_cutoff.csv", check.names=FALSE, stringsAsFactors=F)

# CSV file of % prevalence of potential target expression above target-specific threshold in tumours and normal tissues are brought into R.
prevalence2 <- data.table::fread(file="data/processed/freq_above_target_specific_cutoff.csv", check.names=FALSE, stringsAsFactors=F)

# TSV file of integrated RNAseq data are brought into R.
rna <- data.table::fread(file="data/processed/TCGA_GTEX_integrated_membrane_target.txt", check.names=FALSE, stringsAsFactors=F)

# 2. Data analysis 
# 2.1. Data preprocessing
membrane.target <- colnames(prevalence1)[5:5117]
prevalence2[, Cutoff := NULL]
prevalence2 <- data.table::melt(prevalence2, id.vars = c("Gene"), variable.name = "tissue.cancer", value.name = "value")
prevalence2 <- data.table::dcast(prevalence2, tissue.cancer~Gene, value.var = "value")

# exclude Immunoprivileged tissue
prevalence1[is.na(Risk), Risk := ""]
prevalence1 <- prevalence1[Risk != "Immunoprivileged"]
prevalence2 <- prevalence2[Risk != "Immunoprivileged"]

# normal tissue risk categories
high <- c("Adipose Tissue", "Blood", "Blood Vessel", "Soft tissue/Bone", "Brain", "Esophagus", "Eye", "Heart", "Liver", "Lung", "Lymphatic tissue", "Muscle", "Nerve", "Paraganglia", "Small Intestine", "White blood cell")
 
medium <- c("Adrenal Gland", "Bladder", "Colon", "Kidney", "Pancreas", "Pituitary", "Stomach", "Skin")

low <- c("Breast", "Cervix", "Fallopian Tube", "Ovary", "Salivary Gland", "Prostate", "Rectum", "Spleen", "Thymus", "Thyroid", "Uterus", "Endometrium", "Vagina")

tbd <- c("Lining of body cavities", "Bile duct", "Head and Neck region") 

immunoprivileged <- c("Testis")

# 2.2.Safety bucket
# 2.2.1. Likely safe / good therapeutic window
# Filtering for tumour prevalence ≥20% in at least one indication 
tmp <- prevalence1[type=="cancer"]
for(i in membrane.target){
    if(sum(tmp[, i, with = F] >= 20) > 1){
        if(exists("morethan20perc")){
            morethan20perc <- c(morethan20perc, i)
        }else{
            morethan20perc <- i
        }
    }
}
rm(tmp)

# Are all normal tissues <5% using a 10 TPM threshold?
tmp <- prevalence1[type == "normal"]
for(i in morethan20perc){
    if(sum(tmp[, i, with = F] > 5) == 0)
        if(exists("cleantarget")){
            cleantarget <- c(cleantarget, i)
        }else{
            cleantarget <- i
        }
}
rm(tmp)

# Are all normal tissues <5% using a target-specific threshold?
tmp <- prevalence2[type == "normal"]
for(i in morethan20perc[!morethan20perc %in% cleantarget]){
    if(sum(tmp[, i, with = F] > 5) == 0){
        if(exists("likelysafe")){
            likelysafe <- c(likelysafe, i)
        }else{
            likelysafe <- i
        }
    }else{
        if(exists("potentiallysafe")){
            potentiallysafe <- c(potentiallysafe, i)
        }else{
            potentiallysafe <- i
        }
    }
}
rm(tmp)

# Are all tissues >5% using a 10 TPM threshold ‘low-risk’?
tmp <- prevalence1[type == "normal"]
for(i in likelysafe){
    if( length(unique(tmp$Risk[tmp[, i, with = F] > 5])) == 1 ){
        if(unique(tmp$Risk[tmp[, i, with = F] > 5]) == "low"){
            if(exists("tw.lowrisk")){
                tw.lowrisk <- c(tw.lowrisk, i)
            }else{
                tw.lowrisk <- i
            }
        }
    }
}
rm(tmp)

# Are the only tissues* >5% at 10 TPM GIT (gastro-intestinal tract: Small Intestine, Colon, Stomach, Esophagus)?
tmp <- prevalence1[type == "normal"]
git.tissue <- c("Small Intestine", "Colon", "Stomach", "Esophagus")
for (i in likelysafe[!(likelysafe %in% tw.lowrisk)]){
    if(sum(!tmp$tissue.cancer[tmp[, i, with = F] > 5] %in% git.tissue)==0){
        if(exists("tw.git")){
            tw.git <- c(tw.git, i)
        }else{
            tw.git <- i
        }
    }
}
rm(tmp)

# Are the only tissues* >5% at 10 TPM blood or spleen?
immune.tissue <- c("Blood", "Spleen")
tmp <- prevalence1[type == "normal"]
for (i in likelysafe[!(likelysafe %in% c(tw.lowrisk, tw.git))]){
    if(sum(!tmp$tissue.cancer[tmp[, i, with = F] > 5] %in% immune.tissue)==0){
        if(exists("tw.immune")){
            tw.immune <- c(tw.immune, i)
        }else{
            tw.immune <- i
        }
    }else{
        if(exists("tw.else")){
            tw.else <- c(tw.else, i)
        }else{
            tw.else <- i
        }
    }
}
rm(tmp)

# 2.2.2.Potentially safe / manually assess therapeutic window
# Are all tissues >5% using a target threshold ‘low-risk’?
tmp <- prevalence2[type == "normal"]
for(i in potentiallysafe){
    if(length(unique(tmp$Risk[tmp[, i, with = F] > 5])) == 1 ){
        if(unique(tmp$Risk[tmp[, i, with = F] > 5]) == "low"){
            if(exists("he.lowrisk")){
                he.lowrisk <- c(he.lowrisk, i)
            }else{
                he.lowrisk <- i
            }
        }
    }
}
rm(tmp)

# Are all tissues >5% using a target threshold GIT?
tmp <- prevalence2[type == "normal"]
for(i in potentiallysafe[!(potentiallysafe %in% he.lowrisk)]){
    if(sum(!tmp$tissue.cancer[tmp[, i, with = F] > 5] %in% git.tissue)==0){
        if(exists("tmp.git")){
            tmp.git <- c(tmp.git, i)
        }else{
            tmp.git <- i
        }
    }
}
rm(tmp)

# Are all tissues* >5% at 10 TPM GIT?
tmp <- prevalence1[type == "normal"]
for(i in tmp.git){
    if(sum(!tmp$tissue.cancer[tmp[, i, with = F] > 5] %in% git.tissue)==0){
        tw.git <- c(tw.git, i)
    }else{
        if(exists("he.git")){
            he.git <- c(he.git, i)
        }else{
            he.git <- i
        }
    }
}
rm(tmp)

# Are the only tissues* >5% using a target threshold blood or spleen?
tmp <- prevalence2[type == "normal"]
for(i in potentiallysafe[!(potentiallysafe %in% c(he.lowrisk, tmp.git))]){
    if(sum(!tmp$tissue.cancer[tmp[, i, with = F] > 5] %in% immune.tissue)==0){
        if(exists("tmp.immune")){
            tmp.immune <- c(tmp.immune, i)
        }else{
            tmp.immune <- i
        }
    }
}
rm(tmp)

# Are the only tissues* >5% at 10 TPM blood or spleen?
tmp <- prevalence1[type == "normal"]
for(i in tmp.immune){
    if(sum(!tmp$tissue.cancer[tmp[, i, with = F] > 5] %in% immune.tissue)==0){
        if(exists("immune.only")){
            immune.only <- c(immune.only, i)
        }else{
            immune.only <- i
        }
    }
}
rm(tmp)

# Is the FC between tumour and normal (median blood/spleen) ≥2X
fc_immune <- function(target){
    dt <- rna[, c("Sample.ID", "type", "tissue.cancer", target), with = F]
    data.table::setnames(dt, old = target, new = "TPM")
    dt[, TPM := (2^TPM - 0.001)]
    normal <- dt[tissue.cancer %in% immune.tissue]
    normal.median <- median(normal$TPM)
    tumour <- dt[type=="cancer"]
    tumour <- tumour[, list(Median = median(TPM)), by = "tissue.cancer"]
    tumour[, FC := Median/normal.median]
    tumour <- tumour[FC>=2]
    
    if(nrow(tumour)>=1){
        return(TRUE)
    }else{
        return(FALSE)
    }
}

for(i in tmp.immune){ # tmp.immune[!tmp.immune %in% immune.only]
    if(fc_immune(target = i)){
        if(exists("he.immune")){
            he.immune <- c(he.immune, i)
        }else{
            he.immune <- i
        }
    }else{
        if(exists("nw.immune")){
            nw.immune <- c(nw.immune, i)
        }else{
            nw.immune <- i
        }
    }
}

# Are there ≤2 med/high-risk tissues >5% using a target threshold
tmp <- prevalence2[type == "normal"]
for(i in potentiallysafe[!(potentiallysafe %in% c(tmp.git, tmp.immune, he.lowrisk))]){
    if(sum(tmp$Risk[tmp[, i, with = F] > 5] %in% c("high", "medium")) <= 2){
        if(exists("tmp.restricted")){
            tmp.restricted <- c(tmp.restricted, i)
        }else{
            tmp.restricted <- i
        }
    }else{
        if(exists("highrisk")){
            highrisk <- c(highrisk, i)
        }else{
            highrisk <- i
        }
    }
}
rm(tmp)

# Is the FC between tumour and all high/med risk normal ≥2X
fc_hmrisk <- function(target){
    dt <- rna[, c("Sample.ID", "type", "tissue.cancer", target), with = F]
    data.table::setnames(dt, old = target, new = "TPM")
    dt[, TPM := (2^TPM - 0.001)]
    normal <- dt[tissue.cancer %in% c(high, medium)]
    normal.median <- median(normal$TPM)
    tumour <- dt[type=="cancer"]
    tumour <- tumour[, list(Median = median(TPM)), by = "tissue.cancer"]
    tumour[, FC := Median/normal.median]
    tumour <- tumour[FC>=2]
    
    if(nrow(tumour)>=1){
        return(TRUE)
    }else{
        return(FALSE)
    }
}

for(i in tmp.restricted){
    if(fc_hmrisk(target = i)){
        if(exists("he.restricted")){
            he.restricted <- c(he.restricted, i)
        }else{
            he.restricted <- i
        }
    }else{
        if(exists("nw.risk")){
            nw.risk <- c(nw.risk, i)
        }else{
            nw.risk <- i
        }
    }
}

# 2.3.Final table
# clean target
res <- data.frame(Gene = cleantarget, Recommendations = rep("Clean target", length(cleantarget)))

# Therapeutic window: low risk only
res <- rbind(res, data.frame(Gene = tw.lowrisk, Recommendations = rep("Therapeutic window: low risk only", length(tw.lowrisk))))

# GIT only
res <- rbind(res, data.frame(Gene = tw.git, Recommendations = rep("GIT only", length(tw.git))))

# Therapeutic window: immune only
res <- rbind(res, data.frame(Gene = tw.immune, Recommendations = rep("Therapeutic window: immune only", length(tw.immune))))

# Therapeutic window
res <- rbind(res, data.frame(Gene = tw.else, Recommendations = rep("Therapeutic window", length(tw.else))))

# High expression: low risk only
res <- rbind(res, data.frame(Gene = he.lowrisk, Recommendations = rep("High expression: low risk only", length(he.lowrisk))))

# High expression: GIT only
res <- rbind(res, data.frame(Gene = he.git, Recommendations = rep("High expression: GIT only", length(he.git))))

# Immune only
# res <- rbind(res, data.frame(Gene = immune.only, Recommendations = rep("Immune only", length(immune.only))))

# High expression: Immune only
res <- rbind(res, data.frame(Gene = he.immune, Recommendations = rep("High expression: Immune only", length(he.immune))))

# High expression: restricted
res <- rbind(res, data.frame(Gene = he.restricted, Recommendations = rep("High expression: restricted", length(he.restricted))))

# Widespread high risk expression
res <- rbind(res, data.frame(Gene = highrisk, Recommendations = rep("Widespread high risk expression", length(highrisk))))

# No window with high/med risk tissues
res <- rbind(res, data.frame(Gene = nw.risk, Recommendations = rep("No window with high/med risk tissues", length(nw.risk))))

# No window with immune system
res <- rbind(res, data.frame(Gene = nw.immune, Recommendations = rep("No window with immune system", length(nw.immune))))

# 2.Result 
filename <- paste0("results/candidate_target_list.csv")
if(!exists(filename)){
    data.table::fwrite(res, file=filename)
}
