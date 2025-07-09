library(tidyverse)
library(readxl)
library(purrr)

###############################################################
# This function combine_df will be used to combine
# Multiple sheets from a excel files that have same columns

combine_df<-function(path,sheet_name){
# path is the directory of the excel file in your computer
# sheet_name is a vector of name the sheets that you want to combine
  
  combined_df=data.frame() # Creating empty data frame
  for (sheets in sheet_name) {
    df=read_excel(path,sheet=sheets)
    combined_df<-rbind(combined_df,df)
  }
  return(combined_df)
}
####################################################################
# cleaning_columns function clean/remove unnecessary stuff from a 
# column provided by the user

cleaning_columns<-function(df,named_list){

# df is the data frame that needs to be cleaned
# named_list is where each name is the column name of the data frame and 
# the value is the sub-string that you want to remove from that column
  
  for (name in names(named_list)) {
    for (pattern in named_list[[name]]) {
      df<-df %>% 
        mutate(!!sym(name):=gsub(pattern,"",!!sym(name)))
    }
  }
  return(df) # return the cleaned columns
}

#########################################################################
# This function is used to merge any number of data frames and retains
# all the unique columns and 1 copy of the each common columns
# It takes list of data frames, column name by which all the data frame
# will be joined and join type: either full or inner

merge_df<-function(list_of_df,by="Transcript_ID",join="full",...){
  
  args=list(...)
  if (join=='inner') {
    merged = reduce(list_of_df,inner_join,by=by)
  }
  else{
    merged = reduce(list_of_df, full_join,by=by)
  }
  
  common_cols = Reduce(intersect,lapply(list_of_df, colnames))
  common_cols = common_cols[!common_cols%in%by]

  if (length(args)>0) {
    merged<-merged %>% 
      filter((!!rlang::sym(args$col1) == !!rlang::sym(args$col2)) |
               is.na(!!rlang::sym(args$col1)) |
               is.na(!!rlang::sym(args$col2)))
  } 
  
  
  for (col in common_cols) {
    coalesce_list=names(merged)[grepl(paste0("^",col),names(merged))]
    if (length(coalesce_list) > 1) {  # Only coalesce if duplicates exist
      merged <- merged %>%
        mutate(!!col := coalesce(!!!syms(coalesce_list))) %>%
        select(-all_of(coalesce_list[coalesce_list != col]))  # Keep one version
    }
  }
  return (merged)
}

##########################################################################
# This function is used to merge any number of data frames and retains
# all the unique columns and 1 copy of the each common columns
# It takes list of data frames, column name by which all the data frame
# will be joined and join type: either full or inner

merge_df_v2 <- function(list_of_df, by = "Transcript_ID", join = "full") {
  
  # Start merging from the first data frame
  merged <- list_of_df[[1]]
  
  # Iterate through remaining data frames sequentially
  for (i in 2:length(list_of_df)) {
    next_df <- list_of_df[[i]]
    
    # Find common columns (excluding the join column)
    common_cols <- intersect(colnames(merged), colnames(next_df))
    common_cols <- setdiff(common_cols, by)
    
    # Perform the appropriate join (inner or full)
    if (join == "inner") {
      merged <- inner_join(merged, next_df, by = by)
    } else {
      merged <- full_join(merged, next_df, by = by)
    }
    
    
    
    # Coalesce common columns and remove duplicates
    for (col in common_cols) {
      coalesce_list <- names(merged)[grepl(paste0("^", col), names(merged))]
      if (length(coalesce_list) > 1) {
        merged <- merged %>%
          mutate(!!col := coalesce(!!!syms(coalesce_list))) %>%
          select(-all_of(coalesce_list[coalesce_list != col]))  # Keep only one version
      }
    }
  }
  
  return(merged)
}

##########################################################################
# Function for pca analysis
# Takes a data frame and retun a pca result as a data frame

pca<-function(df,col_to_keep,rowname_column_number=1){
  rn<-df[[rowname_column_number]]
  
  pca_matrix<-df %>%
    select(all_of(col_to_keep)) %>% 
    mutate(across(everything(), ~ replace_na(.x, 0))) %>% 
    as.matrix()
  
  rownames(pca_matrix)<-rn
  
  # Log-transform the data (optional for transcriptomic data)
  log_data <- log2(pca_matrix + 1)
  
  # Transpose the data to have samples as rows
  transposed_data <- t(log_data)
  
  # Perform PCA
  pca_results <- prcomp(transposed_data, scale. = TRUE)
  
  pca_scores <- as.data.frame(pca_results$x)
  return(pca_scores)
}

##########################################################################
# Pathway Assignment
# Function to match and fetch the Pathway-Group
# Pathway_map is a named vector, nmes are pathway, values are group

pathway_assignment <- function(df, pathway_map_named_list) {
  df<-df %>% 
    mutate(Pathway_groups=NA)
  
  for (name in names(pathway_map_named_list)){
    df<-df %>% 
      mutate(Pathway_groups=if_else(is.na(Pathway_groups),
                                    if_else(Pathway%in%pathway_map_named_list[[name]],name,NA),
                                    Pathway_groups))
  } 
  return(df)
}

##########################################################################
# Pathway enrichment

pathway_enrichment <- function(df, gene_col, pathway_col, regulation_col, regulation_type = "UP") {
  # Filter favored and not_favored data
  favored <- df %>% 
    filter(!!sym(regulation_col) == regulation_type) %>% 
    group_by(!!sym(pathway_col)) %>% 
    summarise(gene_count = n_distinct(!!sym(gene_col))) %>% 
    arrange(desc(gene_count)) %>% 
    na.omit() %>% 
    ungroup()
  
  not_favored <- df %>% 
    filter(!!sym(regulation_col) != regulation_type) %>% 
    group_by(!!sym(pathway_col)) %>% 
    summarise(gene_count = n_distinct(!!sym(gene_col))) %>% 
    arrange(desc(gene_count)) %>% 
    na.omit() %>% 
    ungroup()
  
  # Initialize vectors for results
  pathway = character(nrow(favored))
  gene = numeric(nrow(favored))
  a = numeric(nrow(favored))
  b = numeric(nrow(favored))
  c = numeric(nrow(favored))
  d = numeric(nrow(favored))
  p_value = numeric(nrow(favored))
  adjusted_p_value = numeric(nrow(favored))
  enrichment_score = numeric(nrow(favored))
  
  for (i in 1:nrow(favored)) {
    # Assign pathway and gene count from favored data
    pathway[i] = favored[[pathway_col]][i]
    gene[i] = favored$gene_count[i]
    
    # Calculate a and b for Fisher's test (favored group)
    a[i] = gene[i]
    b[i] = sum(favored$gene_count) - a[i]
    
    # Find matching pathways in not_favored and calculate c and d
    matched_index <- which(not_favored[[pathway_col]] == pathway[i])
    if (length(matched_index) > 0) {
      c[i] = not_favored$gene_count[matched_index]
    } else {
      c[i] = 0
    }
    d[i] = sum(not_favored$gene_count) - c[i]
    
    # Ensure all values in the contingency table are non-negative and finite
    if (any(c(a[i], b[i], c[i], d[i]) < 0) || any(is.na(c(a[i], b[i], c[i], d[i])))) {
      p_value[i] = NA  # Set to NA if any value is negative or NA
    } else {
      # Perform Fisher's test
      p_value[i] = fisher.test(matrix(c(a[i], b[i], c[i], d[i]), nrow = 2, byrow = TRUE))$p.value
    }
    
    # Calculate expected count and enrichment score
    total_gene_count = sum(favored$gene_count) + sum(not_favored$gene_count)
    expected_count = total_gene_count / length(unique(df[[pathway_col]]))
    enrichment_score[i] = a[i] / expected_count
  }
  
  # Create final data frame with results
  final_df <- data.frame(Pathway = pathway, Gene = gene, A = a, B = b, C = c, D = d, P_Value = p_value, 
                         Enrichment_Score = enrichment_score)
  final_df <- final_df %>%
    mutate(p_adjust = p.adjust(P_Value, method = "BH")) %>% 
    mutate(Regulation=regulation_type)
  
  return(final_df)
}

##########################################################################
# GO data frame processing function that splits any
# GO column, separate GO ID and Description

go_data_processing<-function(df,type="Biological_Process",separated_by=";"){
  processed_df<-df %>% 
    filter(!is.na(!!sym(type))) %>% 
    separate_rows(!!sym(type),sep=separated_by) %>% 
    extract(!!sym(type), into = c("GO_ID", "Description"), regex = "([^~]+)~(.+)", remove = FALSE)
  
  return(processed_df)
}

#############################################################################
# GO group assignment

go_assignment <- function(df, go_map_named_list,go_type='Description') {
  df<-df %>% 
    mutate(GO_groups=NA)
  
  for (name in names(go_map_named_list)){
    df<-df %>% 
      mutate(GO_groups=if_else(is.na(GO_groups),
                                    if_else(!!sym(go_type)%in%go_map_named_list[[name]],name,NA),
                                    GO_groups))
  } 
  return(df)
}
