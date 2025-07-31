select_membrane_targets <- function(score_threshold = 3) {
  # Load your membrane protein annotations
  knowledge <- data.table::fread(file="data/COMPARTMENT/human_compartment_knowledge_full.tsv", check.names=FALSE, stringsAsFactors=F)

  # Filter based on plasma membrane AND knowledge score
  knowledge <- knowledge[(V4 %in% c("Plasma membrane"))]
  knowledge <- knowledge[V7 >= score_threshold]

  return(unique(knowledge$V2))
}
