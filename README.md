# SSLMC
 Spike and Slab Lasso Matrix Completion

## 1. Tutorial
The [`main.R`](https://github.com/Sijianf/SSLMC/blob/main/codes/main.R) contains the main function `SSLMC()` of our spike and slab lasso matrix completion algorithm.    
The [`functions.cpp`](https://github.com/Sijianf/SSLMC/blob/main/codes/functions.cpp) contains the Rcpp codes that used in the main function.     

This is one example to use the algorithm:
```
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
            # rescale = FALSE
            )
```