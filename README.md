# 🧬 Membrane Target Discovery Pipeline for Cell Therapy

This project provides an end-to-end R pipeline — now with an interactive **Shiny web application** — to identify membrane-bound antigens suitable for targeted cell therapies (e.g., CAR-T, TCR). The workflow integrates membrane protein annotation, RNA expression analysis from tumor and normal tissues, and safety profiling based on tissue specificity and therapeutic window.

---

## 🚀 What's New?

✅ Now includes a fully interactive **RShiny app** to:
- Adjust **TPM expression threshold**
- Set **tumor prevalence** criteria
- Tune **COMPARTMENTS knowledge score threshold**
- Explore and download filtered targets and plots

---

## 📁 Project Structure
membrane-target-discovery/ \
├── app.R # Shiny app interface \
├── scripts/ # Modular R scripts for pipeline steps \
├── data/ # Input data (e.g. COMPARTMENTS, expression) \
├── results/ # Auto-generated output files \
├── config/ # Risk tier definitions, thresholds \
├── docs/ # Method summary and figures \
├── LICENSE \
└── README.md # This file
