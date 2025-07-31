# 1. Load data
# TSV file of integrated RNAseq data are brought into R.
rna <- data.table::fread(file="data/processed/TCGA_GTEX_integrated_membrane_target.txt", check.names=FALSE, stringsAsFactors=F)

# CSV file of candidate target list is brought into R.
target.list <- data.table::fread(file="results/candidate_target_list.csv", check.names=FALSE, stringsAsFactors=F)

# CSV file of % prevalence of potential target expression above target-specific threshold in tumours and normal tissues are brought into R.
prevalence2 <- data.table::fread(file="data/processed/freq_above_cutoff.csv", check.names=FALSE, stringsAsFactors=F)

# CSV file of max expression in normal tissue are brought into R.
max.normal <- data.table::fread(file="data/processed/freq_above_normal_max.csv", check.names=FALSE, stringsAsFactors=F)

# 2. Data analysis {.tabset}
# 2.1 Data preprocessing
# target candidate list
clean.targets <- c("CLDN6", "CD70", "TRPM1", "DLL3", "PCDHB10", "ULBP2", "RAB42", "LRRC55", "WNT7A")
git.only <- c("GUCY2C", "CLCA1", "LYPD8", "MUC2", "NOX1", "GPA33")
therapeutic.window <- c("MSLN", "CD47", "EPCAM", "CEACAM5", "MUC16", "PSCA")
hi.expression.lowrisk <- c("FOLH1", "TRPM8")
hi.expression.restricted <- c("GPC3", "CD276", "CLDN1", "FOLR1")
widespread.normal.tissue.expression <- c("PVR", "RASD1")
no.window.with.hi.med.risk <- c("RHOQ", "RASL12", "VPS4A", "GLG1")
no.window.with.immune <- c("GRK6", "ITGAM", "ADAM8")

target.list <- c(clean.targets, git.only, therapeutic.window, hi.expression.lowrisk, hi.expression.restricted, widespread.normal.tissue.expression, no.window.with.hi.med.risk, no.window.with.immune)

# max expression in normal tissue
max.normal <- max.normal[, c("Gene", "Cutoff"), with = F]
max.normal <- max.normal[-c(1,2,3), ]

# 2.2. Prevalence on tumours
prevalence2[, Cutoff := NULL]

prevalence2 <- data.table::melt(prevalence2, id.vars = c("Gene"), variable.name = "tissue.cancer", value.name = "value")

prevalence2 <- data.table::dcast(prevalence2, tissue.cancer~Gene, value.var = "value")

prevalence2 <- prevalence2[type=="cancer"]

prevalence2[, c("type", "Risk") := NULL]

prevalence2 <- prevalence2[, c("tissue.cancer", "SampleSize", target.list), with = F]

prevalence2 <- data.table::melt(prevalence2, id.vars = c("tissue.cancer", "SampleSize"), variable.name = "Gene", value.name = "prevalence2")

prevalence2[, prevalence2 := as.numeric(prevalence2)]

# 2.3. Fold change to Normal tissue
rna <- rna[, c("Sample.ID", "type", "tissue.cancer", target.list), with = F]

# exclude immunoprivileged tissue
rna <- rna[tissue.cancer!="Testis"]

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

for(i in target.list){
    tmp <- fc_max_normal(target = i)
    
    if(exists("fc2maxnormal")){
        fc2maxnormal <- merge(fc2maxnormal, tmp, by = "tissue.cancer", all = T)
    }else{
        fc2maxnormal <- tmp
    }
}
rm(tmp)

fc2maxnormal <- data.table::melt(fc2maxnormal, id.vars = c("tissue.cancer"), variable.name = "Gene", value.name = "FC2maxnormal")
rm(rna)

sheet <- merge(prevalence2, fc2maxnormal, by = c("tissue.cancer", "Gene"), all = T)
rm(prevalence2, fc2maxnormal)

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

# 2.4. Data presentation
fig4 <- ggplot(data = sheet[prevalence2 > 0.05], aes(x = Gene, y = tissue.cancer))+
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

target.plot <- data.frame(Gene = target.list, Group = c(rep("Clean Targets", length(clean.targets)), rep("GIT only", length(git.only)), rep("Therapeutic window", length(therapeutic.window)), rep("High expression in low-risk only", length(hi.expression.lowrisk)), rep("High expression restricted", length(hi.expression.restricted)), rep("Widespread normal tissue expression", length(widespread.normal.tissue.expression)), rep("No window with hi/med risk tissues", length(no.window.with.hi.med.risk)), rep("No window with immune cells", length(no.window.with.immune))))
target.plot$Gene <- factor(target.plot$Gene, levels = target.list)

fig5 <- ggplot(data = target.plot, aes(x = Gene, y = 1))+
        geom_raster(aes(fill = Group))+
        scale_fill_manual(name = "", values = c("Clean Targets" = "darkgreen", "GIT only" = "purple", "Therapeutic window" = "light blue", "High expression in low-risk only" = "cornsilk", "High expression restricted" = "yellow", "Widespread normal tissue expression" = "orange", "No window with hi/med risk tissues" = "brown4", "No window with immune cells" = "red"))+
        ylim(0.5,1.5)+
        theme_bw()+
        theme(axis.text.x = element_text(size=14, angle=45, hjust=1, vjust=1),
              axis.text.y = element_blank(),
              axis.title = element_blank(),
              plot.title = element_blank(),
              legend.position = c(1.07, 3))

p <- cowplot::plot_grid(fig4, fig5, ncol = 1, align = "v", rel_heights = c(9,1))

# 3. results
filename <- paste0(figure.stamp, "percent_above_TargetCutoff_vs_fc2maxnormal.png")
if(!file.exists(filename)){
    png(file=filename, width=1700, height=920, units = "px")
    print(p)
    dev.off()
}
