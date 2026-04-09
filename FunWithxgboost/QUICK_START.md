# Quick Start: Reproduce Kelly-Malamud-Zhou (2024)

## TL;DR
```bash
# Setup
source /sphenix/u/shuhang98/setup.sh
cd /sphenix/user/shuhangli/ppg12/FunWithxgboost

# Run (data already prepared)
python ridge_reproduce.py    # ~2 minutes, generates ridge_results.csv
python ridge_plot_results.py # ~10 seconds, generates plots in ridge_plots/
```

## What You Get

### Results (`ridge_results.csv`)
- 56 rows: different (complexity, ridge_penalty) combinations
- Columns: R², Sharpe ratio, alpha, information ratio, etc.
- Shows **high complexity → high Sharpe**, even with **negative R²**!

### Plots (`ridge_plots/`)
1. **`voc_curves.png`**: Main result - 6 subplots showing how R², returns, volatility, Sharpe, etc. vary with model complexity
2. **`sharpe_heatmap.png`**: Sharpe ratio as function of (complexity, regularization)
3. **`r2_vs_sharpe.png`**: Scatter plot showing negative R² can yield positive Sharpe
4. **`complexity_comparison.png`**: Bar charts comparing low/medium/high complexity

## Key Result
```
Complexity Category    Sharpe Ratio    R²
Low (c < 0.5)          0.01           -0.17
Medium (0.5-1.5)       0.03           -1.47
High (c > 1.5)         0.22            0.01    ← BEST!
```

**Virtue of Complexity**: More parameters = better performance!

## Customize

Edit `ridge_reproduce.py` to change:
- `train_window`: 12, 60, or 120 months
- `n_rff_list`: Range of feature counts (P/2)
- `ridge_alphas`: Range of regularization strengths
- `n_seeds`: Number of random seeds to average

For full paper reproduction, use P up to 6000 and 1000 seeds (takes hours).

## Files
- `ridge_data_prep.py`: Data preparation (already run)
- `ridge_reproduce.py`: Main analysis
- `ridge_plot_results.py`: Visualization
- `RIDGE_REPRODUCTION_README.md`: Full documentation
