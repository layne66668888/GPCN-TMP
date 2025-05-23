# Load necessary packages and functions
library(data.table)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)

# Setting the working directory
setwd("C:/Users/51195/Desktop/exposure/SNP")

# Create Error Folder
dir.create("C:/Users/51195/Desktop/MRresult/SNPERROR", showWarnings = FALSE)

# Get all .snplocs objects
snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38

# Get all .tsv and .tsv.gz files
files <- list.files(pattern = "\\.(tsv|tsv\\.gz)$")

# Process each file
for (file in files) {
  tryCatch({
    # Reading a file
    df1 <- fread(file)
    setnames(df1, c("chromosome", "base_pair_location", "effect_allele", "other_allele"), c("chr", "pos", "effect_allele", "other_allele"))
    
    # Generate an empty rsid column for the entire dataset
    df1[, rsid := NA_character_]
    
    # Grouping by chromosome
    df1 <- df1[, .SD, by = chr]  # Here we use .SD to represent all columns
    
    for (i in 1:22) {
      chr_snps <- snpsBySeqname(snps, as.character(i))
      idx <- df1[chr == i, match(pos, pos(chr_snps))]
      df1[chr == i, rsid := mcols(chr_snps)$RefSNP_id[idx]]
      print(paste(as.character(i), "processing completed"))
    }
    
    # Generate output file name
    out_file <- gsub("\\.(tsv|tsv\\.gz)$", ".txt", file)
    
    # Writing to a file
    fwrite(df1, out_file, row.names = F, sep = "\t", quote = F)
    
    # Copy the file to the destination folder "C:/Users/51195/Desktop/outcome"
    file.copy(out_file, paste0("C:/Users/51195/Desktop/outcome/", out_file), overwrite = TRUE)
    
    # At the same time copy the files to the target folder "C:/Users/51195/Desktop/exposure"
    file.copy(out_file, paste0("C:/Users/51195/Desktop/exposure/", out_file), overwrite = TRUE)
    
    # Delete the original file
    file.remove(file)
  }, error = function(e) {
    # If an error occurs, write the error information to SNPERROR.txt
    write(file, file = "C:/Users/51195/Desktop/MRresult/SNPERROR/SNPERROR.txt", append = TRUE)
  })
}
