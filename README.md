# Compressed-Sensing-to-Pooling-and-Deconvolution-of-Pd-Catalyzed-Cross-Coupling-Reactions-Experiments
This repository contains Python implementations of the deconvolution algorithms presented in the paper titled "Compressed Sensing Approach to Pooling and Deconvolution of Pd-Catalyzed Cross-Coupling Reactions." The code uses compressed sensing techniques and group testing principles to efficiently screen catalysts in high-throughput experiments. Additionally, this repository includes a guide on integrating these methods as TIBCO Spotfire data functions.

## Table of Contents
- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Getting Started](#getting-started)
- [Integration with TIBCO Spotfire](#integration-with-tibco-spotfire)
- [References](#references)

## Overview
This repository demonstrates how to use combinatorial designs and compressed sensing methods to identify optimal reaction conditions for Pd-catalyzed cross-coupling reactions. The approach allows for efficient deconvolution of pooled experiments, greatly expanding the number of catalysts that can be screened in a given experimental setup.

## Repository Structure
```
/CompressedSensingPooling
│
├── notebooks
│   ├── Deconvolution_Algorithm.ipynb   # Main notebook implementing the deconvolution algorithm
│   ├── Experiment_Analysis.ipynb       # Notebook for analyzing experimental results
│   └── Linearity_Test.ipynb            # Notebook testing the linearity assumption
│
├── data
│   └── sample_data.csv                 # Example data for running the notebooks
│
├── spotfire
│   ├── Spotfire_Integration_Guide.md   # Detailed guide on integrating with Spotfire
│   └── data_functions                  # Folder containing Spotfire data functions (.sfd files)
│
├── requirements.txt                    # Required Python packages
└── README.md                           # This README file
```

## Getting Started
### Prerequisites
- Python 3.8 or later
- Jupyter Notebook

### Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/CompressedSensingPooling.git
   cd CompressedSensingPooling
   ````
2. Install the additional packages:
   ```bash
   pip install -r requirements.txt
   ````
### Required Python Packages
The `requirement.txt` includes the following essential packages:

  ```txt
 ... add stuff
  ```

## Integration with TIBCO Spotfire
### Overview
You can integrate the deconvolution algorithm as a data function in TIBCO Spotfire to enable interactive data analysis. This integration allows the algorithm to process data directly from Spotfire and return results within the Spotfire environment.

## Required Packages
Make sure the following Python packages are installed in the environment connected to Spotfire:

 - `numpy` (our verison is 1.24.1)
 - `pandas` (our verison is 1.5.2)
 - `cvxpy` (our verison is 1.5.3)
- `scipy` (our verison is 1.14.1)

You can search the package already installed and install new ones in `Tools > Python Tools > Package Managment`

## Setting Up Spotfire Integration
1. Enable Python in Spotfire:

 - Go to `Tools > Options > Data Functions` and set up the Python environment. Ensure the path to your Python executable matches the one where you installed the required packages.

2. Import Data Functions:

 - Go to `Tools > Register Data Functions` and select New.
 - Click Import and select the .sfd files located in the spotfire/data_functions folder.

3. Configure Data Function Parameters:

 - Map the input parameters to your data table columns, ensuring that the catalyst pools, experimental results, and solvent/base combinations are correctly linked.
 - Map the output parameters to Spotfire visualization elements.

4. Run the Data Function:

 - Once configured, run the data function to perform the deconvolution on your data directly within Spotfire.

## References
If you find this repository helpful, please cite the following paper:

 - **Compressed Sensing Approach to Pooling and Deconvolution of Pd-Catalyzed Cross-Coupling Reactions** by Lorenzo Tarricone, Jesse Ahlbrecht, Vera Jost, and Georg Wuitschik
