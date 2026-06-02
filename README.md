 
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
- `Baseline Processing Report.pdf and referans_istasyonlarin_turef_koordinatlari.pdf`: Input files containing raw geodetic observations (angles, distances, coordinates).

## 🏃 How to Run 
1. Reference Coordinates (referans_istasyonlarin_turef_koordinatlari.pdf)
The adjustment requires a minimum of two fixed stations (Control Points) in the TUREF (Turkish National Reference Frame) system.
Format: Point ID, Y (Easting), X (Northing), Z (Ellipsoidal Height).
Source: Official TUREF coordinates from your project area.

2. GNSS Baseline Observations (Baseline Processing Report.pdf)
The static GNSS data must be exported from Trimble Business Center (TBC).
Extraction: Open your TBC project and generate the "Baseline Processing Report."Data Points: Look for the "Adjusted Grid Deltas" or "Vector Components" section in the report.
What to copy: You need to extract the Delta Y, Delta X, Delta H (or Delta N, Delta E, Delta h) values and their corresponding variance-covariance values (Standard Deviations) for each baseline.

3. Execution
Place your exported TBC data into the data/ directory.
Ensure your fixed TUREF station IDs match the IDs in your baseline file.
Run the main script:
                      geodetic_adjustment.py
The script will output the Adjusted Coordinates, Precision Criteria, and Chi-Square Test results.
