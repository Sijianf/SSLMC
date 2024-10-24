# SSLMC
 Spike and Slab Lasso Matrix Completion

## 1. Tutorial
The [`main.R`](https://github.com/Sijianf/SSLMC/blob/main/codes/main.R) contains the main function `SSLMC()` of our spike and slab lasso matrix completion algorithm.    
The [`functions.cpp`](https://github.com/Sijianf/SSLMC/blob/main/codes/functions.cpp) contains the Rcpp codes that used in the main function.     

This is one example to use the algorithm:

```
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
            xi = 1
            )
```

- `Y`: The binary dataset to be completed. 
- `K_init`: The column number of the latent space, use $20$ or a larger value for general usage. 
- `max_iter`: The maximum iteration number, you can change as $200$ or $500$ as needed. 
- `show_plot`: This is only for simulation study, set as `FALSE` in real data analysis. 
- `eta`: The learning rate of the algorithm.  
- `xi`: The confidence level for the observed values, use positive integers like $1,2,3,...,10$ for a better performance.  