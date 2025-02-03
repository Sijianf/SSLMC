# Spike and Slab Lasso Matrix Completion (SSLMC)

## 1. Tutorial
The [`main.R`](https://github.com/Sijianf/SSLMC/blob/main/codes/main.R) contains the main function `SSLMC()` of our spike and slab lasso matrix completion algorithm.    
The [`functions.cpp`](https://github.com/Sijianf/SSLMC/blob/main/codes/functions.cpp) contains the Rcpp codes that used in the main function.     

This is one example to use the algorithm:

```r
library(Rcpp)
#library(RcppArmadillo)

source("main.R")
sourceCpp("functions.cpp")

out = SSLMC(Y = Y, 
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

- `Y`: The binary dataset to be completed. 
- `K_init`: The column number of the latent space, use $20$ or a larger value for general usage. 
- `max_iter`: The maximum iteration number, you can change as $200$ or $500$ as needed. 
- `show_plot`: This is only for simulation study, set as `FALSE` in real data analysis. 
- `eta`: The learning rate of the algorithm.  
- `xi`: The confidence level for the observed values, use positive integers like $1,2,3,\cdots,10$ for a better performance.  


## 2. Competing methods

### 2.1 NRLMF

This method is proposed by: 

Liu, Y., Wu, M., Miao, C., Zhao, P., & Li, X. L. (2016). Neighborhood regularized logistic matrix factorization for drug-target interaction prediction. PLoS computational biology, 12(2), e1004760.

Their orginal codes are available at: [PyDTI](https://github.com/stephenliu0423/PyDTI). If you are using python 3 or more recent python versions, you will need to modify these codes or you can directly use my updated codes [here](). 

```bash
#--------------------------------------#
#--------------- PyDTI3 ---------------#
#--------------------------------------#
# This is the compatible codes for python 3.13.0

# I used pyenv to manage the environment: 
pyenv activate PyDTI_venv
pyenv deactivate
pyenv uninstall PyDTI_env

# Running codes are similar, for example: 
python PyDTI.py --method="nrlmf" --dataset="simulation" --predict-num=1 --data-dir="./datasets" --output-dir="./outputs"

```


Once you setup the python environment, you can directly call this function without leaving R. Below is my pipeline to directly run the codes from R console:

```r
# Load in the python environment
library(glue) # useful package to combine strings
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "~/.pyenv/versions/PyDTI_venv/bin/python") # your python path
use_virtualenv("~/.pyenv/versions/3.13.0/envs/PyDTI_venv", required = TRUE) # your virtual environment path
py_config() # configure the above settings


# Define path of your working space:
python_dir <- glue("{local_path}/codes/PyDTI3/PyDTI.py")
data_dir <- glue("{db_path}/{isTune}/datasets")
output_dir <- glue("{db_path}/{isTune}/outputs/PyDTI")
# Create the directories if they don't exist
if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
  message("Created directory: ", data_dir)
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("Created directory: ", output_dir)
}


# Construct the python command
python_command <- glue(
  "python {python_dir} ",
  "--method='nrlmf' ",
  "--dataset='{db}' ",
  "--predict-num=1 ",
  "--data-dir='{data_dir}' ",
  "--method-opt='c={cc2} K1={K1} K2={K2} r={numLat}' ",
  "--output-dir='{output_dir}'"
)


# Modify the file names just for this algorithm: (matched with the name rules in PyDTI.py)
write.table(Y, file = glue("{data_dir}/{db}_admat_dgc.txt"), sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(simD, file = glue("{data_dir}/{db}_simmat_dc.txt"), sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(simT, file = glue("{data_dir}/{db}_simmat_dg.txt"), sep = "\t", row.names = FALSE, col.names = FALSE)


# Execute the python code: 
system(python_command)


# Extract what you will need for downstream analysis: 
A_out <- as.matrix(read.table(glue("{output_dir}/U.txt")))
B_out <- as.matrix(read.table(glue("{output_dir}/V.txt")))
K_out <- ncol(B_out)
```






