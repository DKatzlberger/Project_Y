make.dir <- function(fp, overwrite = FALSE) {
  # Creates a directory 
  # can be specified if you want to overwrite.
  # fb: Path to directory.
  # overwrite: If the existing directory should be overwritten.
  if(!file.exists(fp)) {  
    make.dir(dirname(fp))
    dir.create(fp)
  } else {   if(isTRUE(overwrite)){
    unlink(fp, recursive = TRUE)
    dir.create(fp)
    print('Overwritting existing directory')
  } else {
    print('File exists!')
    } 
  }
} 

normalize_log <- function(X, e = 0.01){
  # Log normalization.
  # X: Vector or Matrix to be transformed.
  # e: threshold to correct for infinity values.
  X = log2(X + e)
}

create_design <- function(output_column, meta){
  # Creates design matrix 
  # Currently supports means model.
  # output_column: Condition to be calculated.
  
  # Design matrix
  design <- model.matrix(~0 + eval(as.name(output_column)), data = meta) 
  # Rename the columns of design matrix
  repl_str <- 'eval(as.name(output_column))'
  colnames(design) <- gsub(repl_str, '', colnames(design), fixed=TRUE)
  # Return
  design
}

create_contrast <- function(coef){
  # Creates contrast matrix.
  # Compares the average of one condition to the average of all other.
  # coef: Coefficients of the design matrix.
  
  # Contrast matrix
  coef1 <- coef
  contrast.matrix <- matrix(0, ncol = length(coef1), nrow=length(coef1))
  row.names(contrast.matrix) <- coef1
  colnames(contrast.matrix) <- coef1
  
  cx <- coef1[1]
  for(cx in coef1){
    colx <- cx
    # colx <- str_replace(cx, 'eval(as.name(setup$classification$output_column))', '')
    contrast.matrix[coef1, colx] <- 1
    contrast.matrix[row.names(contrast.matrix) != cx, colx] <- -1 * (1/(length(coef1)-1))
  }
  # Return
  contrast.matrix
}

extract_results <- function(limmaFit){
  # Extracts results from a topTable for all coefficients.
  # limmaFit: Results from limma fitted model.
  
  limmaRes <- list()
  for (x in colnames(coef(limmaFit))){
    # Extract for each coefficient all genes
    limmaRes[[x]] <- topTable(limmaFit, coef=x, number = Inf) |> 
      rownames_to_column('Gene')
  }
  limmaRes <- bind_rows(limmaRes, .id = "coef") 
  limmaRes <- filter(limmaRes, coef != "(Intercept)") 
  # Return
  limmaRes
}


