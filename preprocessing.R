require(readr)
path <- getwd()
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  cat("\nOne argument must be supplied: T105641_p1@U2AF1.csv")
  cat("\n")
  stop("", call.=FALSE)
} else {
  expansion_list_dir <- c(paste(path, paste("/", args[1], sep = "", collapse = NULL) , sep = "", collapse = NULL))
}

header <- read.table(expansion_list_dir, header = F, sep = ",", nrows = 1)
target_gene <- header[1, "V4"]
target_isoform <- header[1, "V5"]

expansion_list <- read_csv(expansion_list_dir, skip = 1, progress=FALSE, show_col_types = FALSE)
#nrow(expansion_list)
#head(expansion_list)

expansion_list <- expansion_list[expansion_list$type == "gene with protein product",]
n <- round(nrow(expansion_list)/4)
#print(n)
#print(expansion_list[n,]$rank)
expansion_list <- head(expansion_list, n)
#print(expansion_list)
nrow(expansion_list)
#print(header)


nodes_list_name <- paste(target_isoform,"_",n, "_explist.csv", sep = "", collapse = NULL)
nodes_dir <- paste(path, nodes_list_name , sep = "/", collapse = NULL)
write.csv(expansion_list, nodes_dir, row.names=FALSE, quote=TRUE)
