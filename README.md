# ResamplingBRPortfolios

Reproducible codes for the Paper of Oliveira, Trucíos and Valls (2022).

- `Empirical_Application.R` is the main code and perform the results in the empirical application. It calls functions in `Resampling_Techniques.R`, `Auxiliary_Functions.R` and `Performance_Measures.R`
- `Resampling_Techniques` contains the portfolio resampling techniques used in the paper.
- `Auxiliary_Functions.R` contains the auxiliary function `calculate_portfolio_weights` which runs all portfolio allocation methods (Markovitz, Michaud and Factor-Based). The functions receives three arguments: the returns, the constrains and the number of bootstrap replications.
- `Performance_Measures.R` contains functions for the out-of-sample comparison (AV, SD, SR, ASR, SO, TP, SSPW).

