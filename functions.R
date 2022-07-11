#---------------------------------------------------------------------------------------
# Generate protein descriptors
#---------------------------------------------------------------------------------------

generateProteinsDescriptors <- function(){
  
  library(protr)
  
  options(stringsAsFactors = FALSE)
  
  sequences <- list()
  
  for(file in list.files("work/fasta/")){
    sequences[[substr(file,1,6)]] <- readFASTA(paste0("work/fasta/",file))
  }
  
  sequences_df <- data.frame(matrix(unlist(sequences), ncol=1, byrow=FALSE))
  sequences_df <- cbind(names(sequences),sequences_df)
  
  colnames(sequences_df) <- c("uniprot","sequence")
  rownames(sequences_df) <- sequences_df[,1]
  
  sequences_df <- cbind(sequences_df,t(sapply(sequences_df[,2], extractSOCN)))
  
  write.table(sequences_df, file="work/data/protein_descriptors.tsv", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  
  return(sequences_df[,-c(1,2)])
}

#---------------------------------------------------------------------------------------
# Prepare features (without balancing)
#---------------------------------------------------------------------------------------

prepareFeaturesUnbalanced <- function(pdb){
  
  options(stringsAsFactors = FALSE)
  
  pdb_ba <- read.table(paste0("work/data/",pdb,"/bindingAffinity-official-",pdb,".tsv"))
  pdb_ba[,5] <- paste0(pdb_ba[,2],"-",pdb_ba[,3])
  colnames(pdb_ba)[5] <- "ID"  
  
  mut_features <- read.table(paste0("work/data/",pdb,"/pdbbind_pocket_variants_features_",pdb,".tsv"), sep="\t", header=TRUE)
  
  mut_cols_to_keep <- which(colnames(mut_features) %in% c("uniprot",
                                                          "FoldXname",
                                                          "secStruct",
                                                          "Entropy_Complex"))
  mut_cols_to_keep <- c(mut_cols_to_keep[1],mut_cols_to_keep[2],
                        mut_cols_to_keep[3]:mut_cols_to_keep[4])
  
  
  mut_features_reduced <- mut_features[,mut_cols_to_keep]

  lig_features <- read.table(paste0("work/data/",pdb,"/chembl_ligands_features_",pdb,".tsv"), sep="\t", header=TRUE)
  lig_features <- lig_features[!is.na(lig_features[,5]),]
  
  lig_features_reduced <- lig_features[,-grep("Circular",colnames(lig_features))]
  lig_features_reduced <- lig_features_reduced[,-grep("WHIM",colnames(lig_features_reduced))]
  
  features_df <- merge(mut_features_reduced, lig_features_reduced, by = NULL)
  
  features_df[,ncol(features_df)+1] <- paste0(features_df[,2],
                                              "-",
                                              features_df[,ncol(mut_features_reduced)+2])
  colnames(features_df)[ncol(features_df)] <- "ID"  
  
  features_df_joined <- left_join(features_df, pdb_ba[,4:5], by=c("ID" = "ID"))
  colnames(features_df_joined)[ncol(features_df_joined)] <- "ba"
  
 
  features_df_clean <- features_df_joined[,-c(ncol(features_df_joined)-1)]
  
  features_df_clean <- features_df_clean %>% select(uniprot, pdb, FoldXname, ligand_file, chembl_id, tanimoto_index, everything()) 
  
  features_df_clean <- features_df_clean[!is.na(features_df_clean[,ncol(features_df_clean)]),]
  

  protein_descriptors <- read.table("work/data/protein_descriptors.tsv", sep="\t", header=TRUE)
  
  features_df_clean <- left_join(features_df_clean, protein_descriptors[,-2], by=c("uniprot" = "uniprot"))
  
  features_df_clean <- features_df_clean %>% select(-ba,ba) 
    
  
  pocket_descriptors <- read.table(paste0("work/data/",pdb,"/pdbbind_pocket_features_",pdb,".tsv"), sep="\t", header=TRUE)
  
  features_df_clean <- left_join(features_df_clean, pocket_descriptors, by=c("pdb" = "pdb"))
  
  features_df_clean <- features_df_clean %>% select(-ba,ba) 
  
  return(features_df_clean)
}

#---------------------------------------------------------------------------------------
# Prepare features (with balancing)
#---------------------------------------------------------------------------------------

prepareFeaturesBalanced <- function(pdb){
  
  options(stringsAsFactors = FALSE)
  
  pdb_ba <- read.table(paste0("work/data/",pdb,"/bindingAffinity-official-",pdb,".tsv"))
  pdb_ba[,5] <- paste0(pdb_ba[,2],"-",pdb_ba[,3])
  colnames(pdb_ba)[5] <- "ID"  
  
  mut_features <- read.table(paste0("work/data/",pdb,"/pdbbind_pocket_variants_features_",pdb,".tsv"), sep="\t", header=TRUE)
  
  mut_cols_to_keep <- which(colnames(mut_features) %in% c("uniprot",
                                                          "FoldXname",
                                                          "secStruct",
                                                          "Entropy_Complex"))
  
  mut_cols_to_keep <- c(mut_cols_to_keep[1],mut_cols_to_keep[2],
                        mut_cols_to_keep[3]:mut_cols_to_keep[4])
  
  mut_features_reduced <- mut_features[,mut_cols_to_keep]

  lig_features <- read.table(paste0("work/data/",pdb,"/chembl_ligands_features_",pdb,".tsv"), sep="\t", header=TRUE)
  lig_features <- lig_features[!is.na(lig_features[,5]),]
  
  #print(paste("number of ligand rows: ",nrow(lig_features)))
  
  
  if(nrow(lig_features) > 350){
    
    split.index <- createDataPartition(lig_features$tanimoto_index, p = round(350/nrow(lig_features),3), list = FALSE)
    
    lig_features <- lig_features[split.index,]    
  }
    
  lig_features_reduced <- lig_features[,-grep("Circular",colnames(lig_features))]
  lig_features_reduced <- lig_features_reduced[,-grep("WHIM",colnames(lig_features_reduced))]
  
  features_df <- merge(mut_features_reduced, lig_features_reduced, by = NULL)
  
  #print(paste(pdb,ncol(features_df), nrow(features_df)))
  
  features_df[,ncol(features_df)+1] <- paste0(features_df[,2],
                                              "-",
                                              features_df[,ncol(mut_features_reduced)+2])
  colnames(features_df)[ncol(features_df)] <- "ID"  
  
  features_df_joined <- left_join(features_df, pdb_ba[,4:5], by=c("ID" = "ID"))
  colnames(features_df_joined)[ncol(features_df_joined)] <- "ba"
  
  features_df_clean <- features_df_joined[,-c(ncol(features_df_joined)-1)]
  
  features_df_clean <- features_df_clean %>% select(uniprot, pdb, FoldXname, ligand_file, chembl_id, tanimoto_index, everything()) 
  
  features_df_clean <- features_df_clean[!is.na(features_df_clean[,ncol(features_df_clean)]),]
  

  protein_descriptors <- read.table("work/data/protein_descriptors.tsv", sep="\t", header=TRUE)
    
  features_df_clean <- left_join(features_df_clean, protein_descriptors[,-2], by=c("uniprot" = "uniprot"))
    
  features_df_clean <- features_df_clean %>% select(-ba,ba) 
  

  pocket_descriptors <- read.table(paste0("work/data/",pdb,"/pdbbind_pocket_features_",pdb,".tsv"), sep="\t", header=TRUE)
  
  features_df_clean <- left_join(features_df_clean, pocket_descriptors, by=c("pdb" = "pdb"))
    
  features_df_clean <- features_df_clean %>% select(-ba,ba) 
  
  return(features_df_clean)
}

#---------------------------------------------------------------------------------------
# Prepare protein-ligand features (The second ML approach)
#---------------------------------------------------------------------------------------

prepareFeaturesProteinLigand <- function(pdb){
  
  options(stringsAsFactors = FALSE)
  
  pdb_ba <- read.table(paste0("work/data/",pdb,"/bindingAffinity-official-",pdb,".tsv"))
  
  pdb_ba <- pdb_ba[which(grepl("WT",pdb_ba[,2])),]
  pdb_ba <- pdb_ba[,-2]
  
  mut_features <- read.table(paste0("work/data/",pdb,"/pdbbind_pocket_variants_features_",pdb,".tsv"), sep="\t", header=TRUE)
  
  uniprot <- as.data.frame(rep(mut_features$uniprot[1],nrow(pdb_ba)))
  
  pdb_ba <- cbind(uniprot, pdb_ba)
  
  colnames(pdb_ba) <- c("uniprot", "pdb", "ligand_file", "ba")
  
  
  lig_features <- read.table(paste0("work/data/",pdb,"/chembl_ligands_features_",pdb,".tsv"), sep="\t", header=TRUE)
  lig_features <- lig_features[!is.na(lig_features[,5]),]
  
  if(nrow(lig_features) > 350){
    split.index <- createDataPartition(lig_features$tanimoto_index, p = round(350/nrow(lig_features),3), list = FALSE)
    
    lig_features <- lig_features[split.index,]    
  }

  lig_features_reduced <- lig_features[,-grep("Circular",colnames(lig_features))]
  lig_features_reduced <- lig_features_reduced[,-grep("WHIM",colnames(lig_features_reduced))]
  
  features_df <- lig_features_reduced
  
  features_df_joined <- left_join(features_df, pdb_ba[,-2], by=c("ligand_file" = "ligand_file"))
  
  features_df_clean <- features_df_joined[!is.na(features_df_joined[,ncol(features_df_joined)]),]
  

  protein_descriptors <- read.table("work/data/protein_descriptors.tsv", sep="\t", header=TRUE)
  
  features_df_clean <- left_join(features_df_clean, protein_descriptors[,-2], by=c("uniprot" = "uniprot"))
  
  features_df_clean <- features_df_clean %>% select(-ba,ba) 
    
    
  pocket_descriptors <- read.table(paste0("work/data/",pdb,"/pdbbind_pocket_features_",pdb,".tsv"), sep="\t", header=TRUE)
  
  features_df_clean <- left_join(features_df_clean, pocket_descriptors, by=c("pdb" = "pdb"))
  
  features_df_clean <- features_df_clean %>% select(-ba,ba) 
  
  features_cols_to_rm <- which(colnames(features_df_clean) %in% c("uniprot"))
	
  features_df_clean <- features_df_clean[, -features_cols_to_rm]
  
  return(features_df_clean)
}

#---------------------------------------------------------------------------------------
# Prepare mutation features (The second ML approach)
#---------------------------------------------------------------------------------------

prepareFeaturesMutations <- function(pdb, selected_ligands, outputClasses=FALSE){
  
  options(stringsAsFactors = FALSE)
  
  pdb_ba <- read.table(paste0("work/data/",pdb,"/bindingAffinity-official-",pdb,".tsv"))
  pdb_ba <- pdb_ba[,-1]
  colnames(pdb_ba) <- c("FoldXname", "ligand_file","ba")
  
  pdb_ba <- pdb_ba[pdb_ba$ligand_file %in% selected_ligands, ]
  
  pdb_ba_wt <- pdb_ba[which(grepl("WT",pdb_ba[,1])),]
  
  colnames(pdb_ba_wt) <- c("FoldXname", "ligand_file","ba")
  
  pdb_ba_rest <- pdb_ba[which(!grepl("WT",pdb_ba[,1])),]
  
  colnames(pdb_ba_rest) <- c("FoldXname", "ligand_file","ba")
  
  mut_features <- read.table(paste0("work/data/",pdb,"/pdbbind_pocket_variants_features_",pdb,".tsv"), sep="\t", header=TRUE)
  
  features_df <- left_join(pdb_ba_rest, mut_features, by=c("FoldXname"="FoldXname"))
  
  colnames(features_df)[3] <- "mutation_ba"
  
  features_df_joined <- left_join(features_df, pdb_ba_wt, by=c("ligand_file" = "ligand_file"))

  
  features_df_clean <- features_df_joined[!is.na(features_df_joined[,ncol(features_df_joined)]),]
  features_df_clean <- features_df_joined[!is.na(features_df_joined[,3]),]
  
  features_df_clean <- features_df_clean[,-c(4:24)]
  
  colnames(features_df_clean)[ncol(features_df_clean)] <- "wt_ba"
  colnames(features_df_clean)[3] <- "ba"
  
  features_df_clean <- features_df_clean[,-c(ncol(features_df_clean)-1)]
  
  features_df_clean <- cbind(features_df_clean[,-3],features_df_clean[,3])

  colnames(features_df_clean)[ncol(features_df_clean)] <- "ba"
  
  if(outputClasses){
  
	#features_df_clean <- features_df_clean %>% mutate(ba_class = ifelse(ba > wt_ba, "DECREASED", ifelse(ba < wt_ba, "INCREASED", "NEUTRAL")))
	features_df_clean <- features_df_clean %>% mutate(ba_class = ifelse(ba > wt_ba, "DECREASED", "INCREASED"))
	features_df_clean <- features_df_clean[,-(ncol(features_df_clean)-1)]
  }
  
  return(features_df_clean)
}

#---------------------------------------------------------------------------------------
# Split by ligand (weight, TPSA or Volume)
#---------------------------------------------------------------------------------------

getTrainSetByLigand <- function(pdb, full_dataset, by="weight"){
  
  pdb1 <- pdb[1]
  pdb <- pdb[-1]
  
  lig_features <- read.table(paste0("work/data/",pdb1,"/chembl_ligands_features_",pdb1,".tsv"), sep="\t", header=TRUE)
  lig_features <- lig_features[!is.na(lig_features[,5]),]
  
  for(i in pdb){
    
    lig_features_tmp <- read.table(paste0("work/data/",i,"/chembl_ligands_features_",i,".tsv"), sep="\t", header=TRUE)
    lig_features_tmp <- lig_features_tmp[!is.na(lig_features_tmp[,5]),]
    
    lig_features <- rbind(lig_features,lig_features_tmp)
  }
  
  set.seed(123) 
  
  if(by=="weight"){
	train.index <- splitWithCategoricalCheck(lig_features, lig_features$WeightDescriptor)
  }else if(by=="tpsa"){
	train.index <- splitWithCategoricalCheck(lig_features, lig_features$TPSADescriptor)
  }else if(by=="volume"){
	train.index <- splitWithCategoricalCheck(lig_features, lig_features$VABCDescriptor)
  }else{
	train.index <- splitWithCategoricalCheck(lig_features, lig_features$WeightDescriptor)
  }
  
  lig_features_train <- lig_features[train.index,]
  lig_features_test <- lig_features[-train.index,]
  
  train.index <- which(full_dataset$ligand_file %in% lig_features_train$ligand_file)
  
  return(train.index)
}

splitWithCategoricalCheck <- function(results, results_col){

	train.index <- createDataPartition(results_col, p = .8, list = FALSE)

	training_data <- results[train.index,]
	test_data <- results[-train.index,]

	rof <- unique(results$RuleOfFiveDescriptor)
	train_rof <- unique(training_data$RuleOfFiveDescriptor)
	test_rof <- unique(test_data$RuleOfFiveDescriptor)
	
	if(length(intersect(train_rof, rof)) == length(rof) && 
	   length(intersect(test_rof, rof)) == length(rof)){
		return(train.index)
	}else{
		splitWithCategoricalCheck(results, results_col)
	}
}
