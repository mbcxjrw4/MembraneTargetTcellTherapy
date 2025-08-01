# Fold change to target specific threshold
fc_max_normal <- function(target){
    dt <- rna[, c("Sample.ID", "type", "tissue.cancer", target), with = F]
    data.table::setnames(dt, old = target, new = "TPM")
    dt[, TPM := (2^TPM - 0.001)]
    # zerotpm <- min(dt$TPM[dt$TPM>0])
    # dt[TPM<0, TPM := zerotpm]
    target.cutoff <- max.normal$Cutoff[max.normal$Gene == target]
    tumour <- dt[type=="cancer"]
    tumour <- tumour[TPM >= target.cutoff, ]
    tumour <- tumour[, list(Median = median(TPM)), by = "tissue.cancer"]
    tumour[, FC := Median/target.cutoff]
    
    tumour <- tumour[, c("tissue.cancer", "FC"), with = F]
    data.table::setnames(tumour, old = "FC", new = target)
    return(tumour)
}

# overview plot
plot_overview <- function(gene_list, expression_data){
    # % prevalence of potential target expression above target-specific threshold in tumours and normal tissues
    prevalence2 <- percentage_above_target_specific_cutoff(expression_data = expression_data)

    # CSV file of max expression in normal tissue
    max.normal <- percentage_above_normal_max(expression_data = expression_data)

    # max expression in normal tissue
    max.normal <- max.normal[, c("Gene", "Cutoff"), with = F]
    max.normal <- max.normal[-c(1,2,3), ]

    # Prevalence on tumours
    prevalence2[, Cutoff := NULL]
    prevalence2 <- data.table::melt(prevalence2, id.vars = c("Gene"), variable.name = "tissue.cancer", value.name = "value")
    prevalence2 <- data.table::dcast(prevalence2, tissue.cancer~Gene, value.var = "value")
    prevalence2 <- prevalence2[type=="cancer"]
    prevalence2[, c("type", "Risk") := NULL]    
    prevalence2 <- prevalence2[, c("tissue.cancer", "SampleSize", gene_list), with = F]
    prevalence2 <- data.table::melt(prevalence2, id.vars = c("tissue.cancer", "SampleSize"), variable.name = "Gene", value.name = "prevalence2")
    prevalence2[, prevalence2 := as.numeric(prevalence2)]

    # Fold change to Normal tissue
    expression_data <- expression_data[, c("Sample.ID", "type", "tissue.cancer", gene_list), with = F]

    # exclude immunoprivileged tissue
    expression_data <- expression_data[tissue.cancer!="Testis"]
    
    for(i in gene_list){
        tmp <- fc_max_normal(target = i)    
        if(exists("fc2maxnormal")){
            fc2maxnormal <- merge(fc2maxnormal, tmp, by = "tissue.cancer", all = T)
        }else{
            fc2maxnormal <- tmp
        }
    }
    rm(tmp)

    fc2maxnormal <- data.table::melt(fc2maxnormal, id.vars = c("tissue.cancer"), variable.name = "Gene", value.name = "FC2maxnormal")
    
    sheet <- merge(prevalence2, fc2maxnormal, by = c("tissue.cancer", "Gene"), all = T)

    # indication
    sheet[tissue.cancer == "acute myeloid leukemia", tissue.cancer := "LAML"]
    sheet[tissue.cancer == "adrenocortical cancer", tissue.cancer := "ACC"]
    sheet[tissue.cancer == "bladder urothelial carcinoma", tissue.cancer := "BLCA"]
    sheet[tissue.cancer == "brain lower grade glioma", tissue.cancer := "LGG"]
    sheet[tissue.cancer == "breast invasive carcinoma", tissue.cancer := "BRCA"]
    sheet[tissue.cancer == "cervical & endocervical cancer", tissue.cancer := "CESC"]
    sheet[tissue.cancer == "cholangiocarcinoma", tissue.cancer := "CHOL"]
    sheet[tissue.cancer == "chronic myelogenous leukemia", tissue.cancer := "LCML"]
    sheet[tissue.cancer == "colon adenocarcinoma", tissue.cancer := "COAD"]
    sheet[tissue.cancer == "controls", tissue.cancer := "CNTL"]
    sheet[tissue.cancer == "esophageal carcinoma", tissue.cancer := "ESCA"]
    sheet[tissue.cancer == "glioblastoma multiforme", tissue.cancer := "GBM"]
    sheet[tissue.cancer == "head & neck squamous cell carcinoma", tissue.cancer := "HNSC"]
    sheet[tissue.cancer == "kidney chromophobe", tissue.cancer := "KICH"]
    sheet[tissue.cancer == "kidney clear cell carcinoma", tissue.cancer := "KIRC"]
    sheet[tissue.cancer == "kidney papillary cell carcinoma", tissue.cancer := "KIRP"]
    sheet[tissue.cancer == "liver hepatocellular carcinoma", tissue.cancer := "LIHC"]
    sheet[tissue.cancer == "lung adenocarcinoma", tissue.cancer := "LUAD"]
    sheet[tissue.cancer == "lung squamous cell carcinoma", tissue.cancer := "LUSC"]
    sheet[tissue.cancer == "diffuse large B-cell lymphoma", tissue.cancer := "DLBC"]
    sheet[tissue.cancer == "mesothelioma", tissue.cancer := "MESO"]
    sheet[tissue.cancer == "miscellaneous", tissue.cancer := "MISC"]
    sheet[tissue.cancer == "ovarian serous cystadenocarcinoma", tissue.cancer := "OV"]
    sheet[tissue.cancer == "pancreatic adenocarcinoma", tissue.cancer := "PAAD"]
    sheet[tissue.cancer == "pheochromocytoma & paraganglioma", tissue.cancer := "PCPG"]
    sheet[tissue.cancer == "prostate adenocarcinoma", tissue.cancer := "PRAD"]
    sheet[tissue.cancer == "rectum adenocarcinoma", tissue.cancer := "READ"]
    sheet[tissue.cancer == "sarcoma", tissue.cancer := "SARC"]
    sheet[tissue.cancer == "skin cutaneous melanoma", tissue.cancer := "SKCM"]
    sheet[tissue.cancer == "stomach adenocarcinoma", tissue.cancer := "STAD"]
    sheet[tissue.cancer == "testicular germ cell tumor", tissue.cancer := "TGCT"]
    sheet[tissue.cancer == "thymoma", tissue.cancer := "THYM"]
    sheet[tissue.cancer == "thyroid carcinoma", tissue.cancer := "THCA"]
    sheet[tissue.cancer == "uterine carcinosarcoma", tissue.cancer := "UCS"]
    sheet[tissue.cancer == "uterine corpus endometrioid carcinoma", tissue.cancer := "UCEC"]
    sheet[tissue.cancer == "uveal melanoma", tissue.cancer := "UVM"]

    sheet[, tissue.cancer := paste0(tissue.cancer, " (N=", SampleSize, ")")]
    sheet[, SampleSize := NULL]
    sheet[, prevalence2 := prevalence2/100]

    fig <- ggplot(data = sheet[prevalence2 > 0.05], aes(x = Gene, y = tissue.cancer))+
           geom_point(aes(fill = FC2maxnormal, size = prevalence2), shape = 21)+
           scale_size_continuous(name = "% above target specific threshold", breaks = seq(0.2, 1, 0.2), labels = c("20%", "40%", "60%", "80%", "100"))+
           scale_fill_gradientn(name="Fold change to max normal", colours=colorRampPalette(c("white", "firebrick"))(n=20), trans="log10", na.value ="grey")+
           # scale_fill_gradient2(name="Fold change to max normal", low="blue", mid="white", high="firebrick", midpoint=0, na.value ="lightgrey", trans="log10", labels = abs)+
           theme_bw()+
           theme(axis.text.x =  element_blank(), # element_text(size=10, angle=90, hjust=1, vjust=0.5), #
                 axis.text.y = element_text(size=14),
                 axis.title = element_blank(),
                 plot.title = element_blank(),
                 legend.text = element_text(size=10))
    return(fig)
}
