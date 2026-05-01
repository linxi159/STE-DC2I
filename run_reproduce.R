set.seed(1234)

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  install.packages("rmarkdown")
}

rmarkdown::render("reproduce_GSE146771_C14.Rmd")