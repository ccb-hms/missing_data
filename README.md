This repo has two git-backed components, slides/ and workshop/ . If you update the R dependencies in either, make sure to re-run `rsconnect::writeManifest()` in the respective directory.

Rendered output here: https://ccb.connect.hms.harvard.edu/missing_data/

Rendered workshop here: https://ccb.connect.hms.harvard.edu/missing_data_workshop/

# Description

Missingness in data is a ubiquitous and important issue in biomedical research, whether it’s an assay that fails, a subject that declines to answer a survey question, or a lost sample. When it comes time to analyze the data, the choices on handling these missing values can impact the downstream scientific conclusions.

This workshop will provide an introduction to important concepts, strategies, and tools in missing data analysis. We will also discuss several real-world examples of missing data analysis in biomedical research, including examples from epidemiological surveys, genomics, and single-cell multi-omics data. The first hour will be a seminar, the second will be a hands-on workshop where attendees run code.

# Learning objectives

How to:
* Assess the character of missingness in data
* Assess feasible modeling strategies for missing data
* Visualize missing data with `ggmice`
* Impute and model missing data with `mice` and `brms`

# Who should attend

HMS graduate students, postdocs, or faculty who are interested in analyzing data with missing values.

# Prerequisites

The hands-on workshop requires an installation of R and several packages which can be installed with the following commands:

```{r}
pkgs <- c("ggplot2", "dplyr", "mice", "ggmice", "brms")
install.packages(pkgs, Ncpus = 4)
```

Contact andrew_ghazi@hms.harvard.edu with installation questions (or if you have a good missing data analysis problem/example you’d like to share). 

