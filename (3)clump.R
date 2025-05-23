# Load necessary packages and functions
library(data.table)
r_files <- c("ld_clump.R", "ld_matrix.R", "afl2.r", "api.R", "backwards.R", "query.R", "utils-pipe.R", "variants.R", "zzz.R")
lapply(paste0("C:/Users/51195/Desktop/RCODE/", r_files), source)

# Setting the working directory
setwd("C:/Users/51195/Desktop/exposure")

# Get all files
tsv_files <- list.files(pattern = "\\.(gz|tsv|txt|tsv\\.gz|h\\.tsv\\.gz)$")

# Define possible column names
snp_cols <- c("variant_id", "rsid", "rs_id", "SNP", "variant_ID", "RSID", "rsID", "marker_id", "ID", "marker", "rsids")
pval_cols <- c("pval", "p_value", "p")

# Create Error Folder
dir.create("C:/Users/51195/Desktop/MRresult/CLUMPERROR", showWarnings = FALSE)

# Process each file
for (file in tsv_files) {
  tryCatch({
    # Read the first row of data
    first_line <- readLines(gzfile(file), n = 1)
    # Get column names
    col_names <- strsplit(first_line, "\t")[[1]]
    
    # Adjustment of logic for determining SNP column names
    if("rsid" %in% col_names) {
      snp_col <- "rsid"
    } else if("rs_id" %in% col_names) {
      snp_col <- "rs_id"
    } else {
      snp_col <- col_names[which(col_names %in% snp_cols)[1]] # Modify the logic here
    }
    
    pval_col <- col_names[col_names %in% pval_cols][1]
    
    # Reading Data
    expo_rt <- fread(file, header = T)
    
    # Select pval threshold based on file name and format
    if (grepl("^finn.*\\.gz$", file)) {
      pval_threshold <- 5e-6
    } else {
      pval_threshold <- 5e-6
    }
    
    # Select rows where pval is less than the threshold
    expo_rt <- expo_rt[get(pval_col) < pval_threshold,]
    
    # Select the required columns and rename them
    expo_rt2 <- expo_rt[,c(snp_col, pval_col), with = F]
    colnames(expo_rt2) <- c("rsid", "pval")
    
    # Perform LD clumping
    suppressWarnings(clumdf <- ld_clump_local(dat = expo_rt2, clump_kb = 10000, clump_r2 = 0.001, clump_p = 1,
                                              bfile = "C:/Users/51195/Desktop/RCODE/data_maf0/data_maf0.01_rs_ref", 
                                              plink_bin = "C:/Users/51195/Desktop/RCODE/plink_win64_20231018/plink.exe"))
    
    # Results after selecting LD clumping
    expo_rt3 <- expo_rt[which(get(snp_col) %in% clumdf$rsid),]
    
    # Generate a new file name
    new_file <- sub("\\.(gz|tsv|txt|tsv\\.gz|h\\.tsv\\.gz)$", ".txt", file)
    
    # Save the results to a new file
    write.table(expo_rt3, paste0("C:/Users/51195/Desktop/clump/expo_", new_file), row.names = F, sep = "\t", quote = F)
  }, error = function(e) {
    # If an error occurs, write the error message to CLUMPERROR.txt
    write(file, file = "C:/Users/51195/Desktop/MRresult/CLUMPERROR/CLUMPERROR.txt", append = TRUE)
  })
}