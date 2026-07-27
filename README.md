# Want It Like You Mean It: Revisiting Goal Self-Concordance Through the Dissociation of Autonomous and Controlled Motivation in Relation to Meaning in Life

Analysis code and data for the paper published in the *International Journal of Applied Positive Psychology*.

## Citation

Le Guellaff Pallin, L. (2026). Want it like you mean it: Revisiting goal self-concordance through the dissociation of autonomous and controlled motivation in relation to meaning in life. *International Journal of Applied Positive Psychology*. https://doi.org/10.1007/s41042-026-00321-w

## Data Provenance

The dataset (`dataset.csv`) is derived from Sangeorzan (2023), archived at Harvard Dataverse:
[doi:10.7910/DVN/RBZUED](https://doi.org/10.7910/DVN/RBZUED), licensed under CC0 1.0.

Processing applied to the original deposit:

- Aggregation of two independent samples (Study 1 and Study 2; combined N = 437).
- Correction of a data-entry error in CESD-10 item 8 (out-of-range value in the original deposit).
- Computation of composite scores: Autonomous, Controlled, MILjudgements, Depression, Anxiety, GADscore.

All scripts expect `dataset.csv` in the repository root.

## Script Execution Order

| Order | Script | Description |
|-------|--------|-------------|
| 1 | `script1-10.R` | Core analyses: H1 comparative regression, path decomposition, Bayesian ordinal models, H2 controlled motivation, H3 SEM DAG comparison, H4 robustness, depression models, item-level analyses |
| 2 | `script11.R` | Frequentist latent SEM (CFA + DAG comparison + path decomposition) |
| 3 | `script11b.R` | Bayesian latent SEM (blavaan CFA + DAG comparison + path decomposition) |
| 4 | `script13.R` | IEA framework distress-mediation model |
| 5 | `script14.R` | Bidirectional causality test |
| 6 | `script15.R` | Autonomous × Controlled interaction test |
| 7 | `power_analysis.R` | A priori power analysis (standalone; can be run independently) |
| 8 | `make_figures.R` | Publication figures (requires model outputs from scripts above) |

## Requirements

- **R version:** 4.5.2
- **Required packages:** brms, cmdstanr, lavaan, blavaan, dplyr, tidyr, ggplot2, ppcor, pwr, showtext, sysfonts, gridExtra

## Reproducibility

All Bayesian models use `seed = 42` for exact reproducibility of MCMC sampling.

**Note on reported values:** The published article reports the Model 3 goal self-concordance p-value as .215, from the original analysis run; this pipeline reproduces it as .220. All coefficients, credible intervals, and conclusions are identical.

## License

MIT License. See [LICENSE](LICENSE).
