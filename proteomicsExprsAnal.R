library(plyr)
library(tidyverse)
library(gplots)
library(pheatmap)
library(DEP)

raw <- read.csv("raw.csv", header = TRUE)

colnames(raw)

glimpse(raw)

raw1 <- raw %>% select(-c(Protein.Name, Protein.Code, Molecular.Weight:dummy_HEK_Pool_20201025035036.mzML))

data_unique <- make_unique(raw1, "Accession.Number", "Protein.ID")

# Make SummarizedExperiment
expDesign <- data.frame(label = c("dummy_PS25.mzML", "dummy_PS26.mzML", "dummy_PS27.mzML",
                                  "dummy_PS28.mzML", "dummy_PS29.mzML", "dummy_PS30.mzML",
                                  "dummy_PS31.mzML", "dummy_PS32.mzML", "dummy_PS33.mzML"),
                        condition = c("A1", "A1", "A1", "S1", "S1", "S1", "S2", "S2", "S2"),
                        replicate = c(1,2,3,1,2,3,1,2,3))

columns <- 3:11

se <- make_se(data_unique, columns, expDesign)

# Filter, normalize and impute missing values
rawNorm <- se %>% normalize_vsn()

rawImputed <- impute(rawNorm, fun = "MinProb", q = 0.01) 
  
# Test for differentially expressed proteins
diff <- test_diff(rawImputed, "all")
dep <- add_rejections(diff, alpha = 0.05, lfc = 1)
