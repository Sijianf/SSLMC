# **Spike and Slab Lasso Matrix Completion (SSLMC)**

## 📖 Table of Contents
- [Introduction](#introduction)
- [Tutorial of SSLMC](#tutorial-of-sslmc)
- [Competing Methods](#competing-methods)
  - [NRLMF](#nrlmf)

---

## Introduction
The **Spike and Slab Lasso Matrix Completion (SSLMC)** algorithm is designed for efficient matrix completion using a Bayesian framework. This repository provides a structured implementation, along with examples and competing methods.

---

## Tutorial of SSLMC
The [`main.R`](https://github.com/Sijianf/SSLMC/blob/main/codes/SSLMC/main.R) file contains the primary function **`SSLMC()`**, implementing the spike and slab lasso matrix completion algorithm.

Additionally, the [`functions.cpp`](https://github.com/Sijianf/SSLMC/blob/main/codes/SSLMC/functions.cpp) file includes the Rcpp functions used in the main implementation.

### 🚀 **Example Usage**
To run the SSLMC algorithm in **R**, use the following example:

```r
library(Rcpp)

# Load the SSLMC functions
source("main.R")
sourceCpp("functions.cpp")

# Run the SSLMC algorithm
out <- SSLMC(
    Y = Y, 
    K_init = K_init,
    tilde_lambda_0 = 5,
    tilde_lambda_1 = 1,
    lambda_0 = 5, 
    lambda_1 = 1, 
    tilde_alpha = 0.1, 
    tilde_beta = 1,
    alpha = 0.1, 
    beta = 1,
    max_iter = 5000, 
    tol = 1e-5, 
    IBP = 1, 
    show_plot = FALSE,
    eta = 0.001,
    xi = 2
)
```

### 📌 **Parameter Descriptions**
- **`Y`**: Binary dataset to be completed.
- **`K_init`**: Number of latent space columns (default: **≥ 20**).
- **`max_iter`**: Maximum iterations (adjustable to **200** or **500**).
- **`show_plot`**: Used for simulations; set to `FALSE` for real data analysis.
- **`eta`**: Learning rate.
- **`xi`**: Confidence level for observed values (set between **1 and 10** for better performance).

---

## Competing Methods

### NRLMF
This method is proposed by:

> **Liu, Y., Wu, M., Miao, C., Zhao, P., & Li, X. L. (2016).** *Neighborhood regularized logistic matrix factorization for drug-target interaction prediction.* PLoS Computational Biology, **12(2)**, e1004760.  
> 🔗 [Original Code Repository - PyDTI](https://github.com/stephenliu0423/PyDTI)

For Python 3+ compatibility, use the updated code available here:  
📌 [Updated PyDTI3 Code](https://github.com/Sijianf/SSLMC/tree/main/codes/PyDTI3)

### 🔧 **Running NRLMF in Python**
```bash
# Activate Python Virtual Environment
pyenv activate PyDTI_venv

# Run NRLMF with PyDTI
python PyDTI.py --method="nrlmf" --dataset="simulation" --predict-num=1 --data-dir="./datasets" --output-dir="./outputs"

# Deactivate Environment
pyenv deactivate
```

---

### 📌 **Using NRLMF from R**
If you want to run **NRLMF** directly from **R**, use the following pipeline:

```r
# Load required packages
library(glue)
library(reticulate)

# Set Python environment
Sys.setenv(RETICULATE_PYTHON = "~/.pyenv/versions/PyDTI_venv/bin/python")
use_virtualenv("~/.pyenv/versions/3.13.0/envs/PyDTI_venv", required = TRUE)
py_config() # Check Python setup

# Define paths
python_dir <- glue("{local_path}/codes/PyDTI3/PyDTI.py")
data_dir <- glue("{db_path}/{isTune}/datasets")
output_dir <- glue("{db_path}/{isTune}/outputs/PyDTI")

# Create directories if they don’t exist
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Prepare input data
write.table(Y, file = glue("{data_dir}/{db}_admat_dgc.txt"), sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(simD, file = glue("{data_dir}/{db}_simmat_dc.txt"), sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(simT, file = glue("{data_dir}/{db}_simmat_dg.txt"), sep = "\t", row.names = FALSE, col.names = FALSE)

# Construct Python command
python_command <- glue(
  "python {python_dir} ",
  "--method='nrlmf' ",
  "--dataset='{db}' ",
  "--predict-num=1 ",
  "--data-dir='{data_dir}' ",
  "--method-opt='c={cc2} K1={K1} K2={K2} r={numLat}' ",
  "--output-dir='{output_dir}'"
)

# Execute Python script from R
system(python_command)

# Load NRLMF results into R
A_out <- as.matrix(read.table(glue("{output_dir}/U.txt")))
B_out <- as.matrix(read.table(glue("{output_dir}/V.txt")))
K_out <- ncol(B_out)
```
