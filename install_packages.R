if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("graph", "RBGL"))
list_of_packages <- c("readr","pcalg")

new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]

if(length(new_packages)){
  install.packages(new_packages)
} 

