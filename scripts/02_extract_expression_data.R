extract_expression_data <- function(membrane_targets){
    # 1. load data
    # tsv.gz file of gene id and name is brought into R.
    UCSC_Path <- "data/UCSC/"
    file4 <- "TOIL-GTEX_TARGET_TCGA-ROW_DATA.tsv.gz"
    gene <- data.table::fread(file = paste0(UCSC_Path, file4), check.names=FALSE, stringsAsFactors=F)

    # TSV file of meta data is brought into R.
    file5 <-  "TOIL-GTEX_TARGET_TCGA-COLUMN_DATA.tsv"
    meta <- data.table::fread(file=paste0(UCSC_Path, file5), check.names=FALSE, stringsAsFactors=F)

    # TSV file of integrated RNAseq data are brought into R.
    file6 <- "TOIL-GTEX_TARGET_TCGA-ASSAY_DATA-LOG2TPM.tsv"
    rna <- data.table::fread(file=paste0(UCSC_Path, file6), check.names=FALSE, stringsAsFactors=F)

    # 2.Data analysis 
    # Membrane antigen
    membrane_targets <- membrane_targets[membrane_targets %in% gene$gene]

    # UCSC Xena
    # remove data from TARGET
    meta <- meta[Study!="TARGET"]
    meta <- meta[Primary.site!=""]
    meta <- meta[Diagnosis!=""]

    # Thyroid
    meta[Primary.site=="Thyroid Gland", Primary.site := "Thyroid"] 

    # tissue type
    meta[, type := fifelse(Diagnosis=="Normal", "normal", "cancer")]

    # tissue.cancer
    meta[, tissue.cancer := fifelse(Diagnosis=="Normal", Primary.site, Diagnosis)]

    # remove immunoprivileged tissue
    # meta <- meta[Primary.site!="Testis"]

    # remove column
    meta[, c("Primary.site", "Diagnosis", "Gender", "Study") := NULL]

    # rna data
    rna <- data.table::transpose(rna, keep.names = "Sample.ID")

    data.table::setnames(rna, old=colnames(rna), new=c("Sample.ID", gene$gene))

    # select samples
    rna <- rna[Sample.ID%in%meta$Sample.ID]

    # select features
    feature <- c("Sample.ID", membrane_targets)
    rna <- rna[, ..feature]
    rna <- merge(meta, rna, by="Sample.ID")

    # 3. results
    return(rna)
}
