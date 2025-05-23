require(readr)
require(pcalg)
source("functions.R")
options(max.print = .Machine$integer.max)
#options(warn=-1)
path <- getwd()
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0 | length(args)==1 | length(args)==2 | length(args)==3) {
  cat("\nAt least 4 arguments must be supplied: \nexplist.csv organism_type n pcVer")
  cat("\n")
  stop("", call.=FALSE)
} else if (length(args)==4 & (args[4] == "mrf" | args[4] == "cf" | args[4] == "o" | args[4] == "s" | args[4] == "sf" | args[4] == "mr" | args[4] == "c")) {
  expansion_list_dir <- c(paste(path, paste("/", args[1], sep = "", collapse = NULL) , sep = "", collapse = NULL))
} else {
  cat("\nThe available pc versions are: mrf (maj.rule fast), cf (conservative fast), o (original), s (stable), sf (stable.fast), mr (maj.rule), c (conservative)")
  cat("\n")
  stop("", call.=FALSE)
}


# EXPANSION LIST LOADING
expansion_list <- read_csv(expansion_list_dir, skip = 1, progress=FALSE, show_col_types = FALSE)
#head(expansion_list)


if(args[2] == "Hs"){

  # FANTOM MATRIX LOADING
  cat("\nLoading fantom matrix...")
  matrix_dir <- paste(path, "/fantom_mat.csv", sep = "", collapse = NULL)
  input_matrix <- load_matrix(matrix_dir, "Hs")
  cat("\nLoaded fantom matrix.")
  
  # Hs PC INPUT SETTING
  header <- read.table(expansion_list_dir, header = F, sep = ",", nrows = 1)
  target_gene <- header[1, "V4"]
  target_isoform <- header[1, "V5"]
  #print(target_isoform)
  target_tid <- strsplit(target_gene, "-")[[1]][1]
  #print(target_tid)
  pc_input_matrix <- set_pc_input(input_matrix, target_gene, expansion_list, args[3], args[2])
  
} else if (args[2] == "Vv") {
  
  # Vv EXPANSION LIST PROCESSING
  colnames(expansion_list)[2] <- "TID"
  colnames(expansion_list)[11] <- "gene_name_synonyms"
  expansion_list[ , 16:17] <- list(NULL)
  expansion_list$TID <- gsub('.br.','<br>',expansion_list$TID)
  #print(expansion_list)
  
  # VESPUCCI MATRIX LOADING
  cat("\nLoading vespucci matrix...")
  matrix_dir <- paste(path, "/vespucci_mat.csv", sep = "", collapse = NULL)
  input_matrix <- load_matrix(matrix_dir, "Vv")
  cat("\nLoaded vespucci matrix.")
  
  # Vv PC INPUT SETTING
  target_gene <- read.table(expansion_list_dir, header = F, sep = ",", nrows = 1)
  #print(target_gene)
  target_vit <- target_gene[1, "V4"]
  #target_vit <- strsplit(target_vit, "_")[[1]][3] # new expansion list
  #target_vit = paste("VIT_", target_vit, sep="")
  #print(target_vit)
  pc_input_matrix <- set_pc_input(input_matrix, target_vit, expansion_list, args[3], args[2])

} else {
  cat("\nThe second argument can be either Hs (Homo sapiens) or Vv (Vitis vinifera)")
  cat("\n")
  stop("", call.=FALSE)
}

# PC APPLICATION
pc_fit <- apply_pc(pc_input_matrix)
print(pc_fit)

# PC OUTPUT ADJACENCY MATRIX (function as amat)
adj_mat <- as(pc_fit, "amat")
pc_output_df <- edgeListFromPCamat(adj_mat)

if (nrow(pc_output_df) == 0){
  cat("\nThe pc() function did not find any edges, terminating execution...")
  cat("\n")
  stop("", call.=FALSE)
} else {
  
  if(args[2] == "Hs"){
    
    # Hs OUTPUT ANNOTATION
    cat("\nLoading Hs annotation information")
    fantom_anno_table <- read.csv("anno-hsf5.csv", header=FALSE, sep = ",", dec = ".")
    cat("\nLoaded Hs annotation information")
    
    gene_of_interest <- data.frame(
      rank = "0",
      TID = fantom_anno_table$V1[fantom_anno_table$V1 == target_tid],
      Fabs = NA, 
      Frel = 1,
      coords = "?",
      transcript = fantom_anno_table$V2[fantom_anno_table$V1 == target_tid],
      entrezgene_id = paste0("e:",fantom_anno_table$V3[fantom_anno_table$V1 == target_tid]), 
      hgnc_id = paste0("h:",fantom_anno_table$V4[fantom_anno_table$V1 == target_tid]),
      uniprot_id = paste0("u:", fantom_anno_table$V5[fantom_anno_table$V1 == target_tid]),
      gene_name = fantom_anno_table$V6[fantom_anno_table$V1 == target_tid],
      description = fantom_anno_table$V7[fantom_anno_table$V1 == target_tid], 
      type = fantom_anno_table$V8[fantom_anno_table$V1 == target_tid]
    )
    
    expansion_list <- rbind(gene_of_interest, expansion_list)
    nodes_list <- get_nodes(pc_output_df, expansion_list)
    check_deleted_nodes(nodes_list, colnames(pc_input_matrix), expansion_list, args[2])
    nodes_list <- format_Hs(nodes_list)
    nodes_list <- nodes_list[, c(6,10,5,7,8,9,11,12,1,4,2)]
    nodes_list[is.na(nodes_list)] <- ""
    nodes_list$rank <- as.numeric(nodes_list$rank) 
    nodes_list <- nodes_list[order(nodes_list$rank),] 
    cpdag <- get_pearson_corr (pc_output_df, pc_input_matrix) # INITIAL PEARSON CORRELATION (EDGES)
    
    # OUTPUT NODES
    if(grepl("L", args[3])){
      n <- substr(args[3], 1, nchar(args[3])-1)
      n <- as.integer(n)+1
    } else {
      n <- args[3]
    }
    nodes_list_name <- paste(target_isoform,"_",n,"_",args[4], "_nodes.csv", sep = "", collapse = NULL)
    nodes_dir <- paste(path, nodes_list_name , sep = "/Hs/", collapse = NULL)
    write.csv(nodes_list, nodes_dir, row.names=FALSE, quote=FALSE)
    
    # OUTPUT EDGES
    edges_list_name <- paste(target_isoform,"_",n,"_",args[4], "_edges.csv", sep = "", collapse = NULL)
    cpdag_dir <- paste(path, edges_list_name , sep = "/Hs/", collapse = NULL)
    write.csv(cpdag, cpdag_dir, row.names=FALSE, quote=FALSE)
    
  } else if (args[2] == "Vv") {
    
    # Vv OUTPUT ANNOTATION
    cat("\nLoading Vv annotation information")
    vespucci_anno_table_dir <- paste(path, "/anno-vvv2.csv", sep = "", collapse = NULL)
    vespucci_anno_table <- read.csv(vespucci_anno_table_dir, header=FALSE, sep = ",", dec = ".", skip=1)
    cat("\nLoaded Vv annotation information")
    
    gene_of_interest <- data.frame(
      rank = "0",
      TID = vespucci_anno_table$V2[vespucci_anno_table$V2 == target_vit],
      Fabs = NA,
      Frel = 1,
      gene_symbol = vespucci_anno_table$V3[vespucci_anno_table$V2 == target_vit],
      v3 = vespucci_anno_table$V4[vespucci_anno_table$V2 == target_vit],
      group = vespucci_anno_table$V5[vespucci_anno_table$V2 == target_vit],
      subgroup = vespucci_anno_table$V6[vespucci_anno_table$V2 == target_vit],
      pathway_or_family = vespucci_anno_table$V7[vespucci_anno_table$V2 == target_vit],
      full_name = vespucci_anno_table$V8[vespucci_anno_table$V2 == target_vit],
      gene_name_synonyms = vespucci_anno_table$V9[vespucci_anno_table$V2 == target_vit],
      type = vespucci_anno_table$V10[vespucci_anno_table$V2 == target_vit],
      reaction_EC = vespucci_anno_table$V11[vespucci_anno_table$V2 == target_vit],
      Vitisnet_Network1 = vespucci_anno_table$V12[vespucci_anno_table$V2 == target_vit],
      Vitisnet_Network2 = vespucci_anno_table$V13[vespucci_anno_table$V2 == target_vit]
    )
    
    #print(gene_of_interest)
    #print(expansion_list)
    expansion_list <- rbind(gene_of_interest, expansion_list)
    nodes_list <-  get_nodes(pc_output_df,expansion_list)
    check_deleted_nodes(nodes_list, colnames(pc_input_matrix), expansion_list, args[2])
    #nodes_list <- nodes_list[, c(2,5,6,7,8,9,10,11,12,13,14,15,1,4)]
    nodes_list[is.na(nodes_list)] <- ""
    nodes_list$rank <- as.numeric(nodes_list$rank) 
    nodes_list <- nodes_list[order(nodes_list$rank),] 
    cpdag <- get_pearson_corr (pc_output_df, pc_input_matrix) # INITIAL PEARSON CORRELATION (EDGES)
    
    # OUTPUT NODES
    if(grepl("L", args[3])){
      n <- substr(args[3], 1, nchar(args[3])-1)
      n <- as.integer(n)+1
    } else {
      n <- args[3]
    }
    nodes_list_name <- paste(target_vit,"_",n,"_",args[4], "_nodes.csv", sep = "", collapse = NULL)
    nodes_dir <- paste(path, nodes_list_name , sep = "/Vv/", collapse = NULL)
    write.csv(nodes_list, nodes_dir, row.names=FALSE, quote=FALSE)
    
    # OUTPUT EDGES
    edges_list_name <- paste(target_vit,"_",n,"_",args[4], "_edges.csv", sep = "", collapse = NULL)
    cpdag_dir <- paste(path, edges_list_name , sep = "/Vv/", collapse = NULL)
    write.csv(cpdag, cpdag_dir, row.names=FALSE, quote=FALSE)
    
  }
  
}


