# Compressed-Sensing-to-Pooling-and-Deconvolution-of-Pd-Catalyzed-Cross-Coupling-Reactions-Experiments
This repository contains Python implementations of the deconvolution algorithms presented in the paper titled "Compressed Sensing Approach to Pooling and Deconvolution of Pd-Catalysed Cross-Coupling Reactions." An updated draft of the paper can be found [here](https://drive.google.com/file/d/10rbDDRqr5-LMHrAWEbvfnymOKy3wxHzL/view?usp=drive_link). The code uses compressed sensing techniques and group testing principles to efficiently screen catalysts in high-throughput experiments. Additionally, this repository includes a guide on integrating these methods as TIBCO Spotfire data functions.

## Table of Contents
- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Getting Started](#getting-started)
- [Integration with TIBCO Spotfire](#integration-with-tibco-spotfire)
- [References](#references)

## Overview
This repository demonstrates how to use combinatorial designs and compressed sensing methods to identify optimal reaction conditions for Pd-catalysed cross-coupling reactions. The approach allows for efficient deconvolution of pooled experiments, greatly expanding the number of catalysts that can be screened in a given experimental setup.

## Repository Structure
```
/CompressedSensingPooling
│
├── Python Notebooks
│   ├── 2_0_Algorithm.ipynb                   # Object-oriented Implementation of some of the most promising algorithms tested. Assessment of the performance of different algorithms and different averaging strategies was made here with Buchwald-Hartwig reactions
│   ├── 3_0_Algorithm.ipynb                   # Latest version of the algorithm applied to Arylation of Meldrum’s acid and of Barbituric acid. With some intermediate visualisation too
│   ├── Algorithm_playground.ipynb            # In this notebook, there are some early experiments with a  custom implementation of an optimiser and the use of a greedy approach for the optimisation (that didn't give rigorous results
│   ├── Data_analysis_playground.ipynb        # Here, some Exploratory Data Analysis is done on the data coming from Spotfire, and towards the end, there is the script to build a distance matrix (that might be wrong because not a Gaussian score, just distance)
│   ├── Data_preprocessing.ipynb              # This notebook contains some useful functions to visualise and process the plates data that are imported from Spotfire as a .csv file
│   ├── Genetic_Algorithm.ipynb               # Codebase for a genetic algorithm implementation of the solver. The result shows that there is a lot of variability in the parameter and in the definition of exploration and exploitation, and therefore not a robust optimisation technique, probably for this scope
│   ├── Kirkman_matrix_designs.ipynb          # Code for producing the full Kirkman matrices mentioned in the paper
│   ├── Plate_constructor.ipynb               # Plate that from a list of catalysts (here as .csv file, but in the Spotfire documentation it takes directly from a GSheet) gives back the instruction for building a mixed experiment
│   ├── SpotfireProcessing.ipynb              # Implementation of the final version of the code that will be used as a Python data function in Spotfire. The most recent version (with the implementation of the chemical prior) is at the bottom.
|   └── Embeddings_to_D_matrix.ipynb          # Script that allows to pass directly from some embeddings of catalysts to the matrix D to use in the Lasso loss
│ 
├── Data
│   ├── cats_names_35.csv                     # List of 35 catalysts we use for our standard plate design
│   └── embedding_df.csv                      # Dataframe containing the embedding coordinates (2D) that we obtain using our empirical-hand-made prior
│
├── Spotfire
│   ├── Mixing_experiments_user_guide.md       # Detailed guide on integrating with Spotfire our method. 
│   ├── Lorenzo_M_E_deconvolution.sfd         # Spotfire data function file to import to implement the deconvolution algorithm
│   └── Lorenzo_M_E_pooling.sfd         # Spotfire data function file to import to implement the pooling algorithm (a.k.a. designing the experiment)
│
├── requirements.txt                          # Required Python packages
└── README.md                                 # This README file
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
    numpy==1.24.1
    pandas==1.5.2
    cvxpy==1.5.3
    scipy==1.14.1
  ```

## Integration with TIBCO Spotfire
### Overview
You can integrate the deconvolution algorithm as a data function in TIBCO Spotfire to enable interactive data analysis. This integration allows the algorithm to process data directly from Spotfire and return results within the Spotfire environment.

## Required Packages
Make sure the following Python packages are installed in the environment connected to Spotfire:

 - `numpy` (our version is 1.24.1)
 - `pandas` (our version is 1.5.2)
 - `cvxpy` (our version is 1.5.3)
- `scipy` (our version is 1.14.1)

You can search the packages already installed and install new ones in `Tools > Python Tools > Package Management`

## Setting Up Spotfire Integration
1. Enable Python in Spotfire:

 - Go to `Tools > Options > Data Functions` and set up the Python environment. Ensure the path to your Python executable matches the one where you installed the required packages.

2. Import Data Functions:

 - Go to `Tools > Register Data Functions` and select New.
 - Click Import and select the .sfd files located in the Spotfire/data_functions folder.

3. Configure Data Function Parameters:

 - Map the input parameters to your data table columns, ensuring that the catalyst pools, experimental results, and solvent/base combinations are correctly linked.
 - Map the output parameters to Spotfire visualisation elements.

4. Allow for the debugging tool to show printing statements:
 - In order to be able to see the printing messages or assertion errors that are set to pop out if the input is not in the correct format, you should allow the debugging tool to show its output. This setting can be found at `Tools > Options... > Data functions > Enable data function debugging`

5. Run the Data Function:

 - Once configured, run the data function to perform the deconvolution on your data directly within Spotfire.

The analysis of the data can be performed following the instructions of the .md file in the Spotfire folder of this repository

## Usage and suggested workflow
A complete explanation of all the required inputs and the suggested workflow is specified [here](https://github.com/LorenzoTarricone/Compressed-Sensing-to-Pooling-and-Deconvolution-of-Pd-Catalyzed-Cross-Coupling-Reactions-Experiments/blob/main/Spotfire/Mixing_experiments_user_guide.md)

## References
If you find this repository helpful, please cite the following paper:

 - **Compressed Sensing Approach to Pooling and Deconvolution of Pd-Catalyzed Cross-Coupling Reactions** by Lorenzo Tarricone, Jesse Ahlbrecht, Vera Jost, and Georg Wuitschik
