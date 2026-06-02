 
# Geodetic Network Adjustment with Python

This project implements a comprehensive geodetic network adjustment using the **Least Squares Method**. It processes observational geodetic data, performs statistical analysis, and computes optimized coordinates for network points with high precision.

## 🎯 Purpose
The main objective of this project is to automate the rigorous adjustment of local/regional geodetic networks. It eliminates random errors from field measurements (such as distances and angles) to ensure spatial data integrity, which is fundamental for high-precision surveying and mapping engineering.

## 🛠️ Tech Stack & Libraries
- **Language:** Python
- **Data Manipulation:** NumPy, Pandas
- **Visualization (Optional):** Matplotlib (for error ellipses or network geometry)

## 📊 Key Features & Methodology
- **Mathematical Framework:** Built upon Gauss-Markov and Gauss-Helmert adjustment models.
- **Least Squares Estimation:** Solves the matrix equation $A^T P A \hat{x} = A^T P l$ to find optimal coordinate corrections.
- **Statistical Testing:** Implements Global Test (Chi-Square test) for outlier detection and variance component estimation.
- **Precision Analysis:** Computes cofactor matrices, standard deviations, and error ellipses for adjusted points.

## 📁 Repository Structure
- `geodetic_ajustment.py`: Core Python code/notebook executing the matrix operations.
- `Baseline Processing Report.pdf and referans_istasyonlarin_turef_koordinatlari`: Input files containing raw geodetic observations (angles, distances, coordinates).

## 🏃 How to Run
