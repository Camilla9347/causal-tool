require(readr)
source("functions.R")
options(max.print = .Machine$integer.max)
#options(warn=-1)
path <- getwd()
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0 | length(args)==1 | length(args)==2) {
  cat("\nAt least 3 argument must be supplied: \nexplist.csv organism_type n")
  cat("\n")
  stop("", call.=FALSE)
} else if (length(args)==3) {
  expansion_list_dir <- c(paste(path, paste("/", args[1], sep = "", collapse = NULL) , sep = "", collapse = NULL))
}


# EXPANSION LIST LOADING
expansion_list <- read_csv(expansion_list_dir, skip = 1, progress=FALSE, show_col_types = FALSE)

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
  target_tid <- strsplit(target_gene, "-")[[1]][1]
  pc_input_matrix <- set_pc_input(input_matrix, target_gene, expansion_list, args[3], args[2])
  
} else if (args[2] == "Vv") {
  
  # Vv EXPANSION LIST PROCESSING
  colnames(expansion_list)[2] <- "ID"
  colnames(expansion_list)[11] <- "gene_name_synonyms"
  expansion_list[ , 16:17] <- list(NULL)
  expansion_list$ID <- gsub('.br.','<br>',expansion_list$ID)
  
  # VESPUCCI MATRIX LOADING
  cat("\nLoading vespucci matrix...")
  matrix_dir <- paste(path, "/vespucci_mat.csv", sep = "", collapse = NULL)
  input_matrix <- load_matrix(matrix_dir, "Vv")
  cat("\nLoaded vespucci matrix.")
  
  # Vv PC INPUT SETTING
  target_gene <- read.table(expansion_list_dir, header = F, sep = ",", nrows = 1)
  target_vit <- target_gene[1, "V4"]
  pc_input_matrix <- set_pc_input(input_matrix, target_vit, expansion_list, args[3], args[2])

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
      ID = fantom_anno_table$V1[fantom_anno_table$V1 == target_tid],
      Fabs = NA, 
      Frel = 1,
      association_with_transcript = fantom_anno_table$V2[fantom_anno_table$V1 == target_tid],
      entrezgene_id = fantom_anno_table$V3[fantom_anno_table$V1 == target_tid], 
      hgnc_id = fantom_anno_table$V4[fantom_anno_table$V1 == target_tid],
      uniprot_id = fantom_anno_table$V5[fantom_anno_table$V1 == target_tid],
      gene_name = fantom_anno_table$V6[fantom_anno_table$V1 == target_tid],
      description = fantom_anno_table$V7[fantom_anno_table$V1 == target_tid], 
      type = fantom_anno_table$V8[fantom_anno_table$V1 == target_tid]
    )
    
    expansion_list <- rbind(gene_of_interest, expansion_list)
    nodes_list <- get_nodes(pc_output_df, expansion_list)
    check_deleted_nodes(nodes_list, colnames(pc_input_matrix), expansion_list, args[2])
    nodes_list <- format_Hs(nodes_list)
    nodes_list <- nodes_list[, c(2,5,9,6,7,8,10,1,4,11)]
    nodes_list[is.na(nodes_list)] <- ""
    cpdag <- get_pearson_corr (pc_output_df, pc_input_matrix) # INITIAL PEARSON CORRELATION (EDGES)
    
    # OUTPUT NODES
    nodes_dir <- paste(path, paste(target_isoform, "_nodes.csv", sep = "", collapse = NULL) , sep = "/Hs/", collapse = NULL)
    write.csv(nodes_list, nodes_dir, row.names=FALSE, quote=FALSE)
    
    # OUTPUT EDGES
    cpdag_dir <- paste(path, paste(target_isoform, "_edges.csv", sep = "", collapse = NULL) , sep = "/Hs/", collapse = NULL)
    write.csv(cpdag, cpdag_dir, row.names=FALSE, quote=FALSE)
    
  } else if (args[2] == "Vv") {
    
    # Vv OUTPUT ANNOTATION
    cat("\nLoading Vv annotation information")
    vespucci_anno_table_dir <- paste(path, "/anno-vvv2.csv", sep = "", collapse = NULL)
    vespucci_anno_table <- read.csv(vespucci_anno_table_dir, header=FALSE, sep = ",", dec = ".", skip=1)
    cat("\nLoaded Vv annotation information")
    
    gene_of_interest <- data.frame(
      rank = "0",
      ID = vespucci_anno_table$V2[vespucci_anno_table$V2 == target_vit],
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
    
    expansion_list <- rbind(gene_of_interest, expansion_list)
    nodes_list <-  get_nodes(pc_output_df,expansion_list)
    check_deleted_nodes(nodes_list, colnames(pc_input_matrix), expansion_list, args[2])
    nodes_list <- nodes_list[, c(2,5,6,7,8,9,10,11,12,13,14,15,1,4)]
    nodes_list[is.na(nodes_list)] <- ""
    cpdag <- get_pearson_corr (pc_output_df, pc_input_matrix) # INITIAL PEARSON CORRELATION (EDGES)
    
    # OUTPUT NODES
    nodes_dir <- paste(path, paste(target_vit, "_nodes.csv", sep = "", collapse = NULL) , sep = "/Vv/", collapse = NULL)
    write.csv(nodes_list, nodes_dir, row.names=FALSE, quote=FALSE)
    
    # OUTPUT EDGES
    cpdag_dir <- paste(path, paste(target_vit, "_edges.csv", sep = "", collapse = NULL) , sep = "/Vv/", collapse = NULL)
    write.csv(cpdag, cpdag_dir, row.names=FALSE, quote=FALSE)
    
  }
  
}


