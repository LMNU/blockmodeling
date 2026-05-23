# Multiplex Blockmodeling Explorer

A Shiny app to explore how weighting three task interdependence rationales shifts the task grouping produced by generalised blockmodeling.
Produced for a Masters’ thesis submitted to obtain the degree of Master in Sociology at KU Leuven.

**Live app:** https://lmnu-blockmodeling.share.connect.posit.cloud

## What it does

The app simulates a production network of tasks (products × machines) and runs binary blockmodeling across three network layers:

- **M1 — Machine-sharing:** tasks that share the same machine (functional dependency)
- **M2 — Product-sharing:** tasks that belong to the same product (workflow dependency)
- **M3 — Goal conflict:** last-step tasks vs. earlier steps within the same product (i.e. control function)

You adjust the relative weight of each layer and the algorithm selects the optimal number of clusters *k* and partition. The main panel shows the error-by-k table, cluster assignments, and reordered block matrices.

## Requirements

```r
install.packages(c("shiny", "blockmodeling"))
```

## Run locally

```r
shiny::runApp("path/to/app")
```

Or open `app.R` in RStudio and click **Run App**.

## Files

| File | Description |
|------|-------------|
| `app.R` | Single-file Shiny app (UI + server + helpers) |
| `manifest.json` | Deployment manifest for Posit Connect Cloud |

## Notes

- Block type prespecification: `comnul` for M1 and M2, `nuldnc` for M3.
- Restart counts are reduced from the thesis run for interactive speed; the simulated data is well-structured and reliably finds the optimum at 10–20 restarts.
- Helper functions (`sort_partition`, `per_layer_err`, etc.) are copied verbatim from the companion analysis script `sim_unimode.qmd`.
