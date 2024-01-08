load_matrix <- function(matrix_dir, organism){
  input_matrix <-  as.data.frame(read_csv(matrix_dir, show_col_types = FALSE))
  rownames(input_matrix) <- input_matrix[,1]
  input_matrix <- input_matrix[, -1]
  if(organism == "Hs"){
    cat("\nDimension of Fantom matrix (row x col): ", dim(input_matrix))
  }else if (organism == "Vv") {
    cat("\nDimension of Vespucci matrix (row x col): ", dim(input_matrix))
  }
  return(input_matrix)
}

set_pc_input <- function(input_matrix, target_gene, exp_list, cut, organism){
  if(grepl("L", cut)){
    cut <- gsub("L", "", cut)
    cut <- strtoi(cut, base = 0L)
    #cut the list according to number of genes
    tid_list <- exp_list$ID[as.numeric(exp_list$rank) <= as.numeric(cut)]
    rank_list <- exp_list$rank[as.numeric(exp_list$rank) <= as.numeric(cut)]
  } else {
    #cut the list according to relative frequency
    tid_list <- exp_list$ID[as.numeric(exp_list$Frel) >= as.numeric(cut)]
    rank_list <- exp_list$rank[as.numeric(exp_list$Frel) >= as.numeric(cut)]
  }
  if(organism == "Hs"){
    target_tid <- strsplit(target_gene, "-")[[1]][1]
    tid_list <- append(tid_list, target_tid,after=0) #, after=0
  } else if (organism == "Vv") {
    tid_list <- append(tid_list, target_gene,after=0) #, after=0
  }
  rank_list <- append(rank_list, "0", after=0)
  #check_rank = data.frame(unlist(tid_list),unlist(rank_list))
  #print(check_rank)
  cat("\nNumber of pc() input nodes: ", length(tid_list))
  subset_matrix <- data.frame()
  for (i in tid_list) {
    for (j in row.names(input_matrix)) {
      if (i == j) {
        subset_matrix <- rbind(subset_matrix, input_matrix[j, ])
      }
    }
  }
  t_subset_matrix <- t(subset_matrix)
  cat("\nDimension of pc() input matrix (row x col): ", dim( t_subset_matrix))
  cat("\n")
  return(t_subset_matrix)
  
}

apply_pc <- function (pc_input_matrix){
  
  n <- nrow(pc_input_matrix)
  v <- colnames(pc_input_matrix)
  
  start_time <- Sys.time()
  cat("\npc() started at", format(start_time, "%a %b %d %X %Y"))
  
  if(args[4]== "mrf"){
    
    pc_fit <-pc(
      suffStat = list(C = cor(pc_input_matrix), n = n),
      indepTest = gaussCItest, 
      alpha = 0.05,
      labels = v,
      u2pd = "relaxed",
      skel.method = "stable.fast",
      conservative = FALSE,
      maj.rule = TRUE,
      solve.confl = TRUE,
      verbose = FALSE
      
    )
    
  } else if(args[4]== "cf") {
    
    pc_fit <-pc(
      suffStat = list(C = cor(pc_input_matrix), n = n),
      indepTest = gaussCItest, 
      alpha = 0.05,
      labels = v,
      u2pd = "relaxed",
      skel.method = "stable.fast",
      conservative = TRUE,
      maj.rule = FALSE,
      solve.confl = TRUE,
      verbose = FALSE
      
    )
  }
  
  end_time <- Sys.time()
  cat("\npc() ended at", format(end_time, "%a %b %d %X %Y"))
  cat("\npc() took", round(end_time -  start_time, 2))
  cat("\n")
  
  return(pc_fit)
    
}


edgeListFromPCamat <- function (adj_mat){
  
  twos <- which( adj_mat==2, arr.in=TRUE)
  ones <- which( adj_mat==1, arr.in=TRUE)
  
  # BI-DIRECTED (UNDIRECTED/AMBIGUOUS) EDGES (2)
  twos_df = cbind.data.frame(source = colnames(adj_mat)[twos[, 2]],
                             interaction = adj_mat[twos],
                             target = rownames( adj_mat)[twos[, 1]] 
  )
  # there are no bidirected edges
  if (nrow(twos_df)!=0){
    twos_df$source <- gsub('.br.','<br>',twos_df$source)
    twos_df$target <- gsub('.br.','<br>',twos_df$target)
    twos_df <- twos_df[!duplicated(apply(twos_df,1,function(x) paste(sort(x),collapse=' '))),]
    twos_df$interaction <- "<->"
  }
  # DIRECTED EDGES (1) : if gene 1 in column Y and gene 2 in row X have 1, but NOT the other way around
  # colIDs -> source
  # rowIDS -> target
  ones_df = cbind.data.frame(source = colnames( adj_mat)[ones[, 2]],
                             interaction =  adj_mat[ones],
                             target = rownames( adj_mat)[ones[, 1]] 
  )
  ones_df$target <- gsub('.br.','<br>',ones_df$target)
  ones_df$source <- gsub('.br.','<br>',ones_df$source)
  temp <- merge(ones_df, ones_df, by.x = c("source", "target"), by.y = c("target", "source"), all.x = TRUE)
  ones_df <- temp[is.na(temp$interaction.y), ]
  ones_df <- subset(ones_df, select = c(-interaction.y))
  colnames(ones_df)[3] <- "interaction"
  ones_df <- ones_df[, c("source", "interaction", "target")]
  # there are no directed edges
  if(nrow(ones_df!=0)){
    ones_df$interaction <- "-->"
    rownames(ones_df) <- 1:nrow(ones_df)
    ones_df <- ones_df[order(ones_df$source), ]
  }
  # UNDIRECTED EDGES (double 1s) : if gene 1 in column Y and gene 2 in row X have 1, and gene 1 in row Y and gene 2 in column X have 1.
  double_ones_df <- temp[!is.na(temp$interaction.y), ]
  double_ones_df <- subset(double_ones_df, select = c(-interaction.y, -interaction.x))
  double_ones_df <- t(apply(double_ones_df, 1, sort))
  double_ones_df = double_ones_df[!duplicated(double_ones_df),]
  double_ones_df = as.data.frame(double_ones_df)
  if(ncol(double_ones_df)==0){
    # there are no undirected edges
    double_ones_df = cbind.data.frame(source = NULL ,
                                      interaction =  NULL,
                                      target = NULL
    )
  } else if(ncol(double_ones_df) == 1 && nrow(double_ones_df) ==2){
    # only a pair is present
    double_ones_df = cbind.data.frame(source = double_ones_df[1,1] ,
                                      interaction =  "---",
                                      target = double_ones_df[2,1] 
    )
  } else {
    colnames(double_ones_df)[1] <- "source"
    colnames(double_ones_df)[2] <- "target"
    double_ones_df$interaction = "---"
    double_ones_df <- double_ones_df[, c("source", "interaction", "target")]
    double_ones_df$target <- gsub('.br.','<br>',double_ones_df$target)
    double_ones_df$source <- gsub('.br.','<br>',double_ones_df$source)
    rownames(double_ones_df) <- 1:nrow(double_ones_df)
  }
  
  # WHOLE PC OUTPUT
  if(nrow(double_ones_df)==0 && nrow(ones_df)==0 && nrow(twos_df)==0){
    pc_output_df = cbind.data.frame(source = NULL ,
                                    interaction =  NULL,
                                    target = NULL
    )} else {
      pc_output_df <- rbind(double_ones_df, ones_df,twos_df)
      rownames(pc_output_df) <- 1:nrow(pc_output_df)
      return(pc_output_df)
    }
  
}

get_nodes <- function (cpdag, expansion_list){
  lazy_nodes_list <-c(cpdag[["source"]], cpdag[["target"]])
  unique_nodes_list <- unique(lazy_nodes_list)
  nodes_list <- data.frame()
  for (node in unique_nodes_list) {
    nodes_list <- rbind(nodes_list, expansion_list[expansion_list$ID == node,])
  }
  row.names(nodes_list) <- 1:nrow(nodes_list)
  return(nodes_list)
}

check_deleted_nodes <- function(nodes_list, input_nodes, expansion_list,organism) {
  if (length(setdiff(input_nodes, nodes_list[,"ID"])) == 0){ 
    cat("\nAll input nodes belong to the output graph")
  } else {
    #print(setdiff(input_nodes, nodes_list[,"ID"]))
    del_nodes <- setdiff(input_nodes, nodes_list[,"ID"])
    del_nodes_list <- list()
    if(organism == "Hs"){
      for (i in del_nodes) {
        for (j in expansion_list$ID) {
          if(i == j){
            del_nodes_list <- append(del_nodes_list, expansion_list$association_with_transcript[expansion_list$ID==i])
          }
        }
      }
    }else if (organism == "Vv") {
      for (i in del_nodes) {
        for (j in expansion_list$ID) {
          if(i == j){
            del_nodes_list <- append(del_nodes_list, expansion_list$gene_name[expansion_list$ID==i])
          }
        }
      }
    }
    cat("\nThe following input nodes are isolated in the output graph: ", unlist(del_nodes_list))
  }
}


format_Hs <- function(nodes_list){
  i <- 0
  for (i in 1:nrow(nodes_list)){
    if(grepl(",", nodes_list$association_with_transcript[i])){
      nodes_list$association_with_transcript[i]<- gsub(",",";", nodes_list$association_with_transcript[i])
    }
    if(grepl(",", nodes_list$entrezgene_id[i])){
      nodes_list$entrezgene_id[i]<- gsub(",",";", nodes_list$entrezgene_id[i])
    }
    if(grepl(",", nodes_list$hgnc_id[i])){
      nodes_list$hgnc_id[i]<- gsub(",",";", nodes_list$hgnc_id[i])
    }
    if(grepl(",", nodes_list$uniprot_id[i])){
      nodes_list$uniprot_id[i]<- gsub(",",";", nodes_list$uniprot_id[i])
    }
    if(grepl(",", nodes_list$description[i])){
      nodes_list$description[i]<- gsub(",",";", nodes_list$description[i])
    }
    if(grepl(",", nodes_list$type[i])){
      nodes_list$type[i]<- gsub(",",";", nodes_list$type[i])
    }
  }
  return(nodes_list)
}

get_pearson_corr <- function (cpdag, pc_input_matrix){
  for (i in 1:nrow(cpdag)) {
    x <- cpdag[i, "source"]
    x <- gsub(".br.", "<br>", x )
    y <- cpdag[i, "target"]
    y <- gsub(".br.", "<br>", y )
    
    pearson_corr <- cor(pc_input_matrix[,x],pc_input_matrix[,y])
    cpdag[i, "cor"] <- pearson_corr
  }
  cpdag[cpdag$cor > 0, "cor_sign"] <- "+"
  cpdag[cpdag$cor <= 0, "cor_sign"] <- "-"
  return(cpdag)
}
