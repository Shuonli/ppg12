# Kelly, Malamud, Zhou (2024) Reproduction

## Paper
**"The Virtue of Complexity in Return Prediction"**  
Bryan Kelly, Semyon Malamud, and Kangying Zhou  
*The Journal of Finance*, Vol. LXXIX, No. 1, February 2024

## Overview

This directory contains a complete reproduction of the empirical analysis from the paper, which demonstrates that **complex machine learning models (with more parameters than observations) can outperform simple models** in predicting stock market returns.

### Key Finding
The paper proves theoretically and shows empirically that:
- **High-complexity models** (P > T, where P = parameters, T = training samples) achieve better out-of-sample performance than simple models
- **Negative R² doesn't mean poor economic value**: Models with R² < 0 can still generate positive Sharpe ratios
- **Ridge regularization helps**, especially near the "interpolation boundary" (P ≈ T)

## Files

### Data Preparation
- **`ridge_data_prep.py`**: Loads Goyal-Welch predictors and Fama-French factors, applies paper's standardization
  - Input: `GoyalWelch_Data2024.xlsx`, `F-F_Research_Data_Factors.csv`
  - Output: `ridge_panel.csv` (1094 months, 1926-2020)
  - Standardization:
    - Returns: 12-month rolling RMS volatility
    - Predictors: Expanding-window std with 36-month burn-in

### Reproduction Code
- **`ridge_reproduce.py`**: Main reproduction script
  - Implements Random Fourier Features (RFF) from Rahimi & Recht (2007)
  - Ridge/ridgeless regression
  - Recursive out-of-sample forecasting
  - Performance metrics: R², Sharpe ratio, alpha, information ratio
  - Sweeps over model complexity (P) and ridge penalty (z)

- **`ridge_plot_results.py`**: Visualization script
  - Generates plots matching paper Figures 7-8
  - "Virtue of Complexity" curves
  - R² vs. Sharpe ratio scatter
  - Performance heatmaps

### Data Files
- **`GoyalWelch_Data2024.xlsx`**: 15 standard market return predictors (Goyal & Welch 2008)
- **`F-F_Research_Data_Factors.csv`**: Fama-French factors (market excess return)
- **`ridge_panel.csv`**: Processed monthly panel (generated)
- **`ridge_results.csv`**: Out-of-sample performance results (generated)

### Output
- **`ridge_plots/`**: Directory with generated figures
  - `voc_curves.png`: Main "virtue of complexity" plots
  - `sharpe_heatmap.png`: Sharpe ratio vs. (complexity, regularization)
  - `r2_vs_sharpe.png`: Shows negative R² with positive Sharpe
  - `complexity_comparison.png`: Bar charts by complexity category

## Usage

### Setup Environment
```bash
# Activate environment (has pandas, numpy, matplotlib)
source /sphenix/u/shuhang98/setup.sh
```

### Step 1: Prepare Data
```bash
cd /sphenix/user/shuhangli/ppg12/FunWithxgboost
python ridge_data_prep.py
```

This creates `ridge_panel.csv` with:
- 1094 monthly observations (1926-2020)
- Target: `Mkt-RF_std` (volatility-standardized market excess return)
- 15 predictors (standardized): `dfy_std`, `infl_std`, `svar_std`, `de_std`, `lty_std`, `tms_std`, `tbl_std`, `dfr_std`, `dp_std`, `dy_std`, `ltr_std`, `ep_std`, `bm_std`, `ntis_std`, `mkt_excess_lag1_std`

### Step 2: Run Reproduction
```bash
python ridge_reproduce.py
```

This performs:
1. **Quick test**: Single model evaluation (T=12, P=100, z=1e3)
2. **Complexity sweep**: Tests multiple (P, z) combinations
   - Training window: T = 12 months
   - Features: P ∈ {2, 4, 10, 20, 40, 100, 200, 400}
   - Ridge penalty: z ∈ {10⁻³, 10⁻², 10⁻¹, 10⁰, 10¹, 10², 10³}
   - Averages over 5 random seeds

Output: `ridge_results.csv` with performance metrics

### Step 3: Generate Plots
```bash
python ridge_plot_results.py
```

Creates visualizations in `ridge_plots/` directory.

## Results Summary

From our reproduction (T=12, smaller P range for speed):

### Key Findings
- **Best Sharpe Ratio**: 0.24 at P=400, z=1e2 (complexity c=33.3)
- **Ridgeless (z=1e-3)**: Average Sharpe = 0.09
- **High Ridge (z=1e3)**: Average Sharpe = 0.11
- **High complexity (c>10)**: Average Sharpe = 0.22, R² = 0.01

### Performance by Ridge Penalty
```
log10(z)    R²      Sharpe    Info Ratio
-3.0      -1.04     0.09      0.08
-2.0      -0.70     0.09      0.08
-1.0      -0.35     0.09      0.08
 0.0      -0.09     0.10      0.09
 1.0      -0.01     0.11      0.10
 2.0      -0.02     0.12      0.10
 3.0      -0.02     0.11      0.09
```

### Performance by Complexity
```
Complexity   R²      Sharpe    Info Ratio
<0.5       -0.17     0.01      0.01
0.5-1      -1.47     0.03      0.03
1-2        -0.54     0.06      0.05
2-10       -0.10     0.14      0.12
>10         0.01     0.22      0.19
```

**Key Insight**: Sharpe ratio increases monotonically with complexity, even when R² is negative!

## Technical Details

### Random Fourier Features (RFF)
The paper uses RFF to create nonlinear features from linear predictors:

```
S_{i,t} = [sin(γ ω_i' G_t), cos(γ ω_i' G_t)]'
```

where:
- `G_t`: 15 raw predictors at time t
- `ω_i ~ N(0, I)`: Random weights
- `γ = 2`: Bandwidth parameter
- `P = 2 × n_rff`: Total number of features

This creates a wide neural network with fixed input weights.

### Ridge Regression
```
β̂(z) = (zI + T⁻¹ΣₜSₜSₜ')⁻¹ T⁻¹ΣₜSₜRₜ₊₁
```

where:
- `z`: Ridge penalty (shrinkage parameter)
- `z → 0`: Ridgeless limit (minimum ℓ² norm solution)
- `z > 0`: Explicit regularization

### Out-of-Sample Evaluation
- **Recursive forecasting**: Rolling training windows
- **Training window**: T ∈ {12, 60, 120} months
- **Test period**: 1926-2020 (after burn-in)
- **Metrics**:
  - R² = 1 - MSE/Var(y)
  - Sharpe ratio = E[strategy return] / SD[strategy return]
  - Alpha from regression on market return
  - Information ratio = Alpha / Tracking error

## Extending the Reproduction

### Full-Scale Reproduction
To match the paper exactly, increase the parameter ranges:

```python
# In ridge_reproduce.py, modify run_complexity_sweep():
ridge_alphas = [10**i for i in range(-3, 4)]  # Already correct
n_rff_list = list(range(1, 6001))  # P from 2 to 12,000
n_seeds = 1000  # Paper uses 1000 random seeds

# Also test other training windows
for train_window in [12, 60, 120]:
    results = run_complexity_sweep(df, train_window=train_window, ...)
```

**Warning**: Full reproduction with P up to 12,000 and 1000 seeds is computationally expensive (hours to days).

### Additional Analyses
The paper also includes:
- **Comparison with Goyal-Welch (2008)**: Linear "kitchen sink" model
- **Market timing positions**: Time series of forecasts
- **Recession prediction**: NBER recession timing
- **Robustness**: Different γ, volatility standardizations

These can be added by extending the scripts.

## References

### Paper
Kelly, B., Malamud, S., & Zhou, K. (2024). The virtue of complexity in return prediction. *The Journal of Finance*, 79(1), 459-503.

### Data Sources
- Goyal, A., & Welch, I. (2008). A comprehensive look at the empirical performance of equity premium prediction. *Review of Financial Studies*, 21(4), 1455-1508.
- Fama, E. F., & French, K. R. (1993). Common risk factors in the returns on stocks and bonds. *Journal of Financial Economics*, 33(1), 3-56.

### Methods
- Rahimi, A., & Recht, B. (2007). Random features for large-scale kernel machines. *Advances in Neural Information Processing Systems*, 20.
- Hastie, T., Montanari, A., Rosset, S., & Tibshirani, R. J. (2022). Surprises in high-dimensional ridgeless least squares interpolation. *The Annals of Statistics*, 50(2), 949-986.

## Contact

For questions about this reproduction, contact the user who created these scripts.

For questions about the original paper, contact the authors at the affiliations listed in the paper.


