# Membrane Target Discovery Pipeline for Cell Therapy

This repository contains a bioinformatics pipeline for discovering membrane-bound antigen candidates suitable for targeted cell therapy. The workflow integrates membrane protein annotation, RNA expression analysis from tumor and normal tissues, and safety profiling based on tissue specificity and therapeutic window.

## 📁 Project Structure

membrane-target-discovery/

├── data/ # Input data or symbolic links

├── results/ # Pipeline outputs (figures, filtered candidates)

├── scripts/ # Core R scripts for each step of the pipeline

├── config/ # Configuration files (e.g., risk definitions)

├── docs/ # Method description and figures

├── run_pipeline.R # Master script to execute the full workflow

├── README.md # This file

├── LICENSE # License (e.g., MIT, GPL)

└── .gitignore # Files/directories to exclude from Git
