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


## 2. ## Competing methods

### 2.1 NRLMF

This method is proposed by: 

```bibtex
@article{liu2016neighborhood,
  title={Neighborhood regularized logistic matrix factorization for drug-target interaction prediction},
  author={Liu, Yong and Wu, Min and Miao, Chunyan and Zhao, Peilin and Li, Xiao-Li},
  journal={PLoS computational biology},
  volume={12},
  number={2},
  pages={e1004760},
  year={2016},
  publisher={Public Library of Science San Francisco, CA USA}
}
```

Their orginal codes are available at: [PyDTI](https://github.com/stephenliu0423/PyDTI). If you are using python 3 or more recent python versions, you will need to modify these codes or you can directly use my updated codes [here](). 

```bash
# Install python for a specific project 

# brew install python

# which python3

# echo $PATH | grep --color=auto "$(pyenv root)/shims"

# pyenv virtualenv 3.13.0 my_env


#--------------------------------------#
#--------------- PyDTI ----------------#
#--------------------------------------#
# This is the orginal codes written in python 2.7.9

# I used conda to manage the environment (because pyenv cannot fully install the old python in Macbook M chip): 
# conda env create -n PyDTI_conda python=2.7.13
# conda env create -f environment.yaml
# conda activate PyDTI_conda
# conda deactivate
# conda env remove --n PyDTI_conda
# conda install -c conda-forge numpy scipy

python PyDTI.py --method="nrlmf" --dataset="nr" --predict-num=1 --data-dir="./datasets" --output-dir="./outputs"

python sat_analysis.py

python PyDTI.py --method="nrlmf" --dataset="nr" --cvs=1 --specify-arg=0 --data-dir='./datasets'

python PyDTI.py --method="nrlmf" --dataset="nr" --cvs=1 --specify-arg=1 --method-opt="c=5 K1=5 K2=5 r=100 lambda_d=0.125 lambda_t=0.125 alpha=0.25 beta=0.125 theta=0.5" --data-dir="./datasets"

#--------------------------------------#
#--------------- PyDTI3 ---------------#
#--------------------------------------#
# This is the compatible codes for python 3.13.0

# I used pyenv to manage the environment: 
# pyenv activate PyDTI_venv
# pyenv deactivate
# pyenv uninstall PyDTI_env

# running codes are similar


python PyDTI.py --method="nrlmf" --dataset="simulation" --predict-num=1 --data-dir="./datasets" --output-dir="./outputs"

```









