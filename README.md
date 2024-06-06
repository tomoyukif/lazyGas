# lazyGas

`lazyGas` is an R package designed for conducting genetic association studies. It provides tools for loading, preprocessing, and analyzing genomic and phenotype data, as well as for visualizing the results.

## Installation

### Prerequisites

Before installing `lazyGas`, ensure you have the following R packages installed:

- `GBScleanR`
- `ggplot2`

You can install these packages using the following commands in R:
```R
install.packages("ggplot2")
```

To install `GBScleanR`, please visit the [GitHub repository](https://github.com/tomoyukif/GBScleanR)

### Installing lazyGas

To install `lazyGas`, clone the repository and install the package using the following commands in R:
```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
    
devtools::install_github("tomoyukif/lazyGas", build_vignettes = TRUE)
```
## Vignette
For more information, run the following code on a R console to see a vignette.

browseVignettes(package = "lazyGas")

## Usage

Below is a brief overview of how to use `lazyGas` to conduct a genetic association study. For more detailed instructions, refer to the vignette.

### Loading Input Data

Load an input GDS file and build a `LazyGas` class object:
```R
library(lazyGas)
library(GBScleanR)

gds_fn <- "path/to/your/input.gds"
lg <- buildLazyGas(gds_fn = gds_fn, load_filter = TRUE, overwrite = FALSE)
```

### Loading and Preprocessing Phenotype Data

Load phenotype data and prepare it for analysis:
```R
pheno_fn <- "path/to/your/phenotype_data.csv"
pheno <- read.csv(file = pheno_fn)
sample_id <- getSamID(object = lg)
pheno$ID <- paste0("Sample_", sprintf("%03d", pheno$ID))

# Calculate average scores and handle invalid data
pheno$score <- rowMeans(pheno[, -1], na.rm = TRUE)
pheno$score[pheno$rep4 <= 3 & !is.na(pheno$rep4)] <- NA

# Assign phenotype data to LazyGas object
lg <- assignPheno(object = lg,
                  pheno = pheno,
                  rename = c("Trait1", "Trait2", "Trait3"))
```

### Conducting Association Analysis

Run the association scan:
```R
makeDF_FUN <- function(h){
    out <- cbind(colSums(h == 1), colSums(h == 3))
    colnames(out) <- c("hap1", "hap3")
    return(out)
}

scanAssoc(object = lg,
          formula = "hap1 + hap3",
          makeDF_FUN = makeDF_FUN,
          geno_format = "haplotype",
          kruskal = NULL)
```

### Visualizing Results

Generate and save Manhattan plots:
```R
pheno <- getPheno(object = lg)
pdf("manhattan.pdf")
for(i in seq_along(pheno$pheno_names)){
    p <- plotManhattan(object = lg, pheno = i)
    p <- p + labs(title = pheno$pheno_names[i])
    print(p)
}
dev.off()
```

### Extracting Data

Extract data for further analysis:
```R
scan_pheno_1 <- lazyData(object = lg, dataset = "scan", pheno = pheno$pheno_names[1])
peakcall_pheno_5 <- lazyData(object = lg, dataset = "peakcall", pheno = pheno$pheno_names[5])
recalc_pheno_2 <- lazyData(object = lg, dataset = "recalc", pheno = pheno$pheno_names[2])
```

## Contributing

We welcome contributions to `lazyGas`. Please submit issues and pull requests on GitHub.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
```
