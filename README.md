# raschkidmap
# RaschKIDMAP for Baby Growth Online

RaschKIDMAP for Baby Growth Online is a lightweight Shiny application for multidimensional assessment of infant growth from 0 to 24 months. It compares an input infant's height, weight, and head circumference with 50 sex- and age-matched reference infants, then presents Rasch-based and conventional percentile results.

## Features

- Accepts sex, month age, height, weight, and head-circumference inputs.
- Matches each input with 50 reference infants of the same sex and month of age.
- Supports three analysis options:
  - Probability-based Continuous Response Model (CRM)
  - TAM Rating Scale Model (linear 0–4)
  - TAM Partial Credit Model (linear 0–4)
- Produces individual and group forest plots, KIDMAP, Wright Map, item/person tables, input-item residual Z-score details, and traditional percentile plots.
- Includes a homepage language menu for English, Simplified Chinese, Traditional Chinese, French, Japanese, Korean, Spanish, German, and Arabic.
- Saves the selected homepage language in a browser cookie for one year.
- Runs locally without a database backend.

## Requirements

- R 4.0 or later
- R packages:

```r
install.packages("shiny")
install.packages("TAM") # Required only for the TAM models
```

`shiny` is required. The Continuous Response Model does not require `TAM`; select a TAM model only after installing that package.

## Run locally

1. Clone or download this repository.
2. Keep `app.R` and `baby.csv` in the same folder.
3. Open R or RStudio in that folder.
4. Run:

```r
shiny::runApp()
```

Alternatively:

```r
shiny::runApp("path/to/jmirbaby")
```

## Input data

The application reads the bundled `baby.csv` file. It must include these columns:

```text
sample_id, gender, month, height, weight, head_circ
```

For each sex-by-month combination used in the application, the reference data must contain exactly 50 records. The input infant is appended as the 51st record during analysis.

## Interpretation

- **Person measure:** Rasch-based estimate of the infant's relative multidimensional growth position.
- **Item difficulty (Delta):** Relative location of height, weight, or head circumference on the model scale.
- **Infit/Outfit MNSQ:** Model-fit statistics for persons and items.
- **KIDMAP residual Z-score:** Difference between observed and model-expected item values, scaled by its residual standard deviation.
- **Traditional percentile:** Empirical percentile of each original measurement within the 50 matched reference infants.

The results are designed to support exploration and communication of multidimensional growth patterns. They do not replace clinical evaluation or validated growth standards.

## Project structure

```text
.
├── app.R        # Shiny application
├── baby.csv     # Bundled reference data
└── README.md
```

## License

No license has been specified for this repository. Add a license file before distributing or reusing the code.
