package.skeleton(name="beast", code_files="fullMCMC_complexityPrior.R")
R CMD build beast
R CMD check --as-cran beast_version.tar.gz
R CMD INSTALL beast


