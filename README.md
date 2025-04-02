# soybean-ai
A Novel Framework for Soybean Phenotype Prediction and Salient Loci Mining via Machine Learning and Interpretability Analysis

# Soybean Genomic Analysis and Prediction Project

This project aims to predict soybean oil content, protein content, and water-soluble protein content based on Single Nucleotide Polymorphism (SNP) data using various machine learning (ML) and deep learning (DL) models. It also includes SHAP (SHapley Additive exPlanations) for model interpretability analysis, and performs GWAS (Genome-Wide Association Studies) analysis on the SNP data.

## Project Structure

The project contains the following folders and files:


## Main Features

### 1. **ML Folder**
- This folder contains code that uses machine learning models (e.g., Support Vector Machines, Random Forests, etc.) to predict soybean oil content, protein content, and water-soluble protein content.
- The code processes SNP data and trains and validates models for prediction.

### 2. **DL Folder**
- This folder contains code for predicting using deep learning models (e.g., neural networks).
- Similar to the ML models, the deep learning models train on SNP data to predict soybean characteristics.

### 3. **SHAP Folder**
- This folder uses SHAP (SHapley Additive exPlanations) to perform model interpretability analysis on the best-performing Support Vector Regression (SVR) model.
- The SHAP analysis helps to understand how different features (such as SNPs) influence the model's predictions.

### 4. **data_processing Folder**
- This folder contains code for performing GWAS (Genome-Wide Association Studies) analysis on SNP data.
- The data is pre-processed, and statistical methods are applied to perform association analysis, producing results for further research.




