# Load necessary libraries
library(TwoSampleMR)
library(writexl)
library(readxl)

# Create Error Folder
dir.create("C:/Users/51195/Desktop/MRresult/NORESULT", showWarnings = FALSE)

# Define possible column names
snp_cols <- c("rsids", "rsid", "rs_id", "variant_id", "SNP", "variant_ID", "RSID", "rsID", "marker_id", "ID", "marker")
pval_cols <- c("p_value", "pval", "p")
beta_cols <- c("beta", "BETA", "b")
se_cols <- c("standard_error", "se", "SE", "ste", "sebeta")
effect_allele_cols <- c("effect_allele", "EA", "alt", "A1")
other_allele_cols <- c("other_allele", "non_effect_allele", "OA", "ref", "A2")
eaf_cols <- c("effect_allele_frequency", "eaf", "effect_allele_Freq", "EAF", "freq_effect_allele", "af_alt", "freq")

# Define function to check which columns exist in the file
check_cols <- function(filename, cols) {
  data <- read.table(filename, header = TRUE, sep = "\t", nrows = 5)
  colnames <- colnames(data)
  for (col in cols) {
    if (col %in% colnames) {
      return(col)
    }
  }
  stop(paste("The required column was not found in the file", filename))
}

# List all exposure and outcome data files
expo_files <- list.files("C:\\Users\\51195\\Desktop\\clump", pattern = "\\.txt$", full.names = TRUE)
outc_files <- list.files("C:\\Users\\51195\\Desktop\\outcome", pattern = "\\.txt$|\\.tsv$|\\.tsv\\.gz$|\\.h\\.tsv\\.gz$|\\.gz$", full.names = TRUE)

# Process each exposure and outcome data file
for (i in 1:length(expo_files)) {
  # Create an empty data frame to store all the results
  all_results <- data.frame()
  for (j in 1:length(outc_files)) {
    tryCatch({
      # Check and get column names
      if (grepl("expo_finngen_R", expo_files[i])) {
        snp_col <- "rsids"
        beta_col <- "beta"
        se_col <- "sebeta"
        effect_allele_col <- "alt"
        other_allele_col <- "ref"
        eaf_col <- "af_alt"
        pval_col <- "pval"
      } else {
        col_data <- read.table(expo_files[i], header = TRUE, sep = "\t", nrows = 5)
        col_names <- colnames(col_data)
        # Adjustment of logic for determining SNP column names
        if("rsid" %in% col_names) {
          snp_col <- "rsid"
        } else if("rs_id" %in% col_names) {
          snp_col <- "rs_id"
        } else {
          snp_col <- col_names[which(col_names %in% snp_cols)[1]] # Modify the logic here
        }
        pval_col <- check_cols(expo_files[i], pval_cols)
        beta_col <- check_cols(expo_files[i], beta_cols)
        se_col <- check_cols(expo_files[i], se_cols)
        effect_allele_col <- check_cols(expo_files[i], effect_allele_cols)
        other_allele_col <- check_cols(expo_files[i], other_allele_cols)
        eaf_col <- check_cols(expo_files[i], eaf_cols)
      }
      
      # Reading exposed data
      expo_rt <- read_exposure_data(
        filename = expo_files[i],
        sep = "\t",
        snp_col = snp_col,
        beta_col = beta_col,
        se_col = se_col,
        effect_allele_col = effect_allele_col,
        other_allele_col = other_allele_col,
        eaf_col = eaf_col,
        pval_col = pval_col
      )
      
      # Check and get the column name of the result data
      if (grepl("finn", outc_files[j])) {
        snp_col <- "rsids"
        beta_col <- "beta"
        se_col <- "sebeta"
        effect_allele_col <- "alt"
        other_allele_col <- "ref"
        eaf_col <- "af_alt"
        pval_col <- "pval"
      } else {
        col_data <- read.table(outc_files[j], header = TRUE, sep = "\t", nrows = 5)
        col_names <- colnames(col_data)
        # Adjustment of logic for determining SNP column names
        if("rsid" %in% col_names) {
          snp_col <- "rsid"
        } else if("rs_id" %in% col_names) {
          snp_col <- "rs_id"
        } else {
          snp_col <- col_names[which(col_names %in% snp_cols)[1]] # Modify the logic here
        }
        pval_col <- check_cols(outc_files[j], pval_cols)
        beta_col <- check_cols(outc_files[j], beta_cols)
        se_col <- check_cols(outc_files[j], se_cols)
        effect_allele_col <- check_cols(outc_files[j], effect_allele_cols)
        other_allele_col <- check_cols(outc_files[j], other_allele_cols)
        eaf_col <- check_cols(outc_files[j], eaf_cols)
      }
      
      # Read result data
      outc_rt <- read_outcome_data(
        snps = expo_rt$SNP,
        filename = outc_files[j],
        sep = "\t",
        snp_col = snp_col,
        beta_col = beta_col,
        se_col = se_col,
        effect_allele_col = effect_allele_col,
        other_allele_col = other_allele_col,
        eaf_col = eaf_col,
        pval_col = pval_col
      )
      
      # If the number of rows in the result data is 0, skip the subsequent operations and record the file name
      if (nrow(outc_rt) == 0) {
        write(paste(basename(expo_files[i]), "&", basename(outc_files[j])), file = "C:\\Users\\51195\\Desktop\\MRresult\\NORESULT\\NORESULT.txt", append = TRUE)
        next
      }
      
      # Integrate and merge data
      harm_rt <- harmonise_data(
        exposure_dat =  expo_rt, 
        outcome_dat = outc_rt,
        action = 2
      )
      
      # Mendelian randomization analysis was performed, limited to MR Egger and Inverse variance weighted methods
      mr_result <- mr(harm_rt, method_list = c("mr_egger_regression", "mr_ivw"))
      
      # Generation ratio
      OR = generate_odds_ratios(mr_result)
      mr_result <- cbind(mr_result, OR[,5:ncol(OR)])
      
      # Remove Duplicate Columns
      mr_result <- mr_result[ , !duplicated(colnames(mr_result))]
      
      # Remove the "outcome" and "exposure" columns
      mr_result$outcome <- NULL
      mr_result$exposure <- NULL
      
      # Generate output file name
      expo_filename <- basename(expo_files[i])
      outc_filename <- basename(outc_files[j])
      expo_filename <- sub("\\.txt$", "", expo_filename)
      outc_filename <- sub("\\.txt$", "", outc_filename)
      outc_filename <- sub("\\.h\\.tsv\\.gz$", "", outc_filename)
      outc_filename <- sub("\\.gz$", "", outc_filename)
      outc_filename <- sub("\\.tsv\\.gz$", "", outc_filename)
      outc_filename <- sub("\\.tsv$", "", outc_filename)
      output_filename <- paste0("C:\\Users\\51195\\Desktop\\MRresult\\", 
                                expo_filename, 
                                "&", 
                                outc_filename, 
                                ".xlsx")
      # Create a dedicated results folder
      result_folder <- paste0("C:\\Users\\51195\\Desktop\\MRresult\\", expo_filename)
      dir.create(result_folder, showWarnings = FALSE)
      
      output_filename <- paste0(result_folder, "\\", expo_filename, "&", outc_filename, ".xlsx")
      
      # Modify the contents of the id.exposure and id.outcome columns
      mr_result$id.exposure <- gsub("^expo_", "", expo_filename)
      mr_result$id.outcome <- outc_filename
      
      # Write the results to the output file
      write_xlsx(mr_result, output_filename)
      
      # Add the result to all results
      all_results <- rbind(all_results, mr_result)
    }, error = function(e) {
      # If an error occurs, write the error message to NORESULT.txt
      write(paste(basename(expo_files[i]), "&", basename(outc_files[j])), file = "C:\\Users\\51195\\Desktop\\MRresult\\NORESULT\\NORESULT.txt", append = TRUE)
    })
  }
  # Create a general results folder
  dir.create("C:/Users/51195/Desktop/MRresult/ALLRESULT", showWarnings = FALSE)
  
  # Generate output file name
  expo_filename <- basename(expo_files[i])
  expo_filename <- sub("\\.txt$", "", expo_filename)
  expo_filename <- gsub("^expo_", "", expo_filename)
  all_output_filename <- paste0("C:\\Users\\51195\\Desktop\\MRresult\\ALLRESULT\\", expo_filename, ".xlsx")
  
  # Write all the results to the total table
  write_xlsx(all_results, all_output_filename)
}
