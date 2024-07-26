# ResamplingBRPortfolios

Reproducible code for the paper by Oliveira, Trucíos, and Valls (2024). The script `01_Empirical_Application_IBRX.R` replicates the tables in the paper and `02_out_of_sample_tests.R` performs the bootstrap test. Additional results (for other values of $\lambda$ are in the folder `Results`. 

The function in the `covShrinkage` paste were obtained from https://github.com/MikeWolf007/covShrinkage.

> Due to Economatica license, the data provided here cannot be used elsewhere and is made available only for reproducibility purposes.

## References

- Oliveira, A., Trucíos, C., and Valls, P. (2024). _Portfolio resampling in the Brazilian stock market: Can it outperforms Markowitz optimization?._ (Submitted).


**Obs:**
After conducting the empirical application, we decided to name the methods Linear Shrinkage 1, 2, 3, and 4 based on the dates they were published. Thus:

| Text Files in `Results` | Paper       |
|:----------:|:-----------:|
| ls_0       | Sample covariance  |
| ls_1       | Linear Shrinkage 4 |
| ls_2       | Linear Shrinkage 1 |
| ls_3       | Linear Shrinkage 3 |
| ls_4       | Linear Shrinkage 2 |

