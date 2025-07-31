# Import data
# TSV files of data from protein subcellular localisation database COMPARTMENTS [https://compartments.jensenlab.org/Search] is brought into R.
filePath1 <- "data/COMPARTMENT/"
file1 <-  "human_compartment_experiments_full.tsv"
experiments <- data.table::fread(file=paste0(filePath1, file1), check.names=FALSE, stringsAsFactors=F)

file2 <-  "human_compartment_knowledge_full.tsv"
knowledge <- data.table::fread(file=paste0(filePath1, file2), check.names=FALSE, stringsAsFactors=F)

file3 <-  "human_compartment_predictions_full.tsv"
predictions <- data.table::fread(file=paste0(filePath1, file3), check.names=FALSE, stringsAsFactors=F)

file4 <-  "human_compartment_textmining_full.tsv"
textmining <- data.table::fread(file=paste0(filePath1, file4), check.names=FALSE, stringsAsFactors=F)

# Antigen searching space {.tabset}
# Data selection
# COMPARTMENTS database - In knowledge and predication channels, for each gene, data from different resource, get the highest confidence 
# knowledge channel - select condidence index >= 3
knowledge <- knowledge[(V4 %in% c("Plasma membrane"))]
knowledge <- knowledge[V7 >= 3]
knowledge <- knowledge[, c("V1", "V2", "V7"), with = F]
data.table::setnames(knowledge, old = c("V1", "V2", "V7"), new = c("Gene ID", "Gene", "Knowledge"))
knowledge <- knowledge[, list(Knowledge = max(Knowledge)), by = list(`Gene ID`, Gene)]

# experiments channel 
experiments <- experiments[(V1 %in% knowledge$`Gene ID`) & (V4 %in% c("Plasma membrane"))]
experiments <- experiments[, c("V1", "V2", "V7"), with = F]
data.table::setnames(experiments, old = c("V1", "V2", "V7"), new = c("Gene ID", "Gene", "Experiments"))

# predictions channel
predictions <- predictions[(V1 %in% knowledge$`Gene ID`) & (V4 %in% c("Plasma membrane"))]
predictions <- predictions[, c("V1", "V2", "V7"), with = F]
data.table::setnames(predictions, old = c("V1", "V2", "V7"), new = c("Gene ID", "Gene", "Predictions"))
predictions <- predictions[, list(Predictions = max(Predictions)), by = list(`Gene ID`, Gene)]

# textmining channel
textmining <- textmining[(V1 %in% knowledge$`Gene ID`) & (V4 %in% c("Plasma membrane"))]
textmining <- textmining[, c("V1", "V2", "V6"), with = F]
data.table::setnames(textmining, old = c("V1", "V2", "V6"), new = c("Gene ID", "Gene", "Textmining"))

# Combine the dataset
sheet <- merge(knowledge, experiments, by = c("Gene ID", "Gene"), all.x = T)
sheet <- merge(sheet, predictions, by = c("Gene ID", "Gene"), all.x = T)
sheet <- merge(sheet, textmining, by = c("Gene ID", "Gene"), all.x = T)

# Result output
filename <- paste0("data/processed/membrane_target_candidates.csv")
if(!exists(filename)){
    data.table::fwrite(sheet, file=filename)
}
