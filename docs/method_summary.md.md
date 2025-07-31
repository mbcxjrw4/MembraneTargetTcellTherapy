# Introduction
With access to both TRuC and HiT platforms, the identification and prioritization of membrane targets stand as pivotal components within T-cell therapy's future pipeline. Previous methodologies tailored for interrogating gene expression primarily focused on pHLA targets, rendering them often inadequate for precise assessment of membrane targets. Consequently, we've found ourselves reliant on competitors and existing literature for insights. However, through the implementation of a novel workflow, we've integrated subcellular localization predictions, relative risk profiles of normal tissues, and distinctive features of membrane proteins in cancer. This innovative approach enables us to conduct automated and unbiased identification and triage of membrane targets, offering a more precise estimation of their prevalence in tumours and furnishing a quantitative prediction of the therapeutic window.

# Method
- Antigen search space: To assemble the list of potential membrane targets we initiated the process by utilizing a list of putative transmembrane genes then filtered by localization of the plasma membrane as annotated in COMPARTMENTS database with high confidence (level 3 or higher in the knowledge channel).

- Gene expression data: RNA-seq data across 9795 samples taken from TCGA and 7425 samples from GTEx were extracted in the TPM format from the UCSC database, which had been uniformly processed by the Toil pipeline

- Target-specific thresholding is determined by considering the expression level of the target in healthy tissues, utilizing a scaling factor weighted by the risk categories of normal tissues. The minimal is set at 10 TPM

- Safety bucket: Target candidates are first filtered by tumour prevalence (≥20% in at least one indication) and assessed for good safety profiles. Expression is restricted to gastrointestinal tract (GIT) or epithelial tissue, validated with tissue polarity and literature. Targets are further examined for the therapeutic window (≥10 TPM) and for predicted unsafe instances based on fold change between normal and tumour tissues

<img width="395" height="221" alt="image" src="https://github.com/user-attachments/assets/d6785f9c-72be-47ad-becb-08d367983e4c" />
