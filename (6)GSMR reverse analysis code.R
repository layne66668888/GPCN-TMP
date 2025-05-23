# Load the required libraries
library(gsmr2)
library(TwoSampleMR)
library(ggplot2)
library(foreach)
library(readxl)
library(openxlsx)  # For writing xlsx files
library(data.table)  # For fast reading of large files

# Read N.xlsx file
N_df <- read_excel("C:/Users/51195/Desktop/N/N.xlsx")

# Get all files in the exposure_data and outcome_data folders
exposure_files <- list.files("C:/Users/51195/Desktop/clump", pattern = "\\.txt$", full.names = TRUE)
outcome_files <- list.files(
  "C:/Users/51195/Desktop/outcome",
  pattern = "\\.txt$|\\.tsv$|\\.tsv\\.gz$|\\.h\\.tsv\\.gz$|\\.gz$",
  full.names = TRUE
)

# Create a results folder if it does not exist
mrresult_dir <- "C:/Users/51195/Desktop/MRresult"
if (!dir.exists(mrresult_dir)) {
  dir.create(mrresult_dir)
}

# Create the ALLRESULT folder if it does not exist
allresult_dir <- file.path(mrresult_dir, "ALLRESULT")
if (!dir.exists(allresult_dir)) {
  dir.create(allresult_dir)
}

# Create the NORESULT folder if it does not exist
noresult_dir <- file.path(mrresult_dir, "NORESULT")
if (!dir.exists(noresult_dir)) {
  dir.create(noresult_dir)
}

# Setting the Threshold
thresh <- 5e-6

# Define possible column names
snp_cols <- c("rsids", "rsid", "rs_id", "variant_id", "SNP", "variant_ID", "RSID", "rsID", "marker_id", "ID", "marker")
pval_cols <- c("p_value", "pval", "p")
beta_cols <- c("beta", "BETA", "b")
se_cols <- c("standard_error", "se", "SE", "ste", "sebeta")
effect_allele_cols <- c("effect_allele", "EA", "alt", "A1")
other_allele_cols <- c("other_allele", "non_effect_allele", "OA", "ref", "A2")
eaf_cols <- c("effect_allele_frequency", "eaf", "effect_allele_Freq", "EAF", "freq_effect_allele", "af_alt", "freq")
samplesize_cols <- c("N", "samplesize", "n", "N_samples", "sample_size", "samples", "sample_count", "total_n", "total_samples", "exposure_N", "outcome_N", "control_N", "case_N", "study_N", "population_N", "cohort_size", "participant_count", "meta_N", "met")

# Define function to check which columns exist in the file
check_cols <- function(filename, cols) {
  data <- read.table(filename, header = TRUE, sep = "\t", nrows = 5)
  colnames <- colnames(data)
  for (col in cols) {
    if (col %in% colnames) {
      return(col)
    }
  }
  stop(paste("文件", filename, "中未找到所需的列"))
}

# Loop through each exposed file
for (expo_file in exposure_files) {
  tryCatch({
    # Extract the exposure ID (assuming the file name is in the format "expo_ID.txt")
    expo_filename <- basename(expo_file)
    
    # Extract GCST ID or other ID (based on file name)
    gcst_match <- regexpr("GCST[0-9]+", expo_filename)
    
    if (gcst_match != -1) {
      # If GCST is included, extract the GCST ID
      expo_id <- regmatches(expo_filename, gcst_match)
    } else {
      # If GCST is not included, the file name after removing the extension is used as the ID
      expo_id <- sub("^expo_", "", sub("\\.txt$|\\.tsv$|\\.tsv\\.gz$|\\.h\\.tsv\\.gz$|\\.gz$", "", expo_filename))
    }
    
    # Find the corresponding sample size in N.xlsx
    samplesize <- N_df$N[N_df$ID == expo_id]
    
    # If the corresponding sample size is not found, a prompt is given and the current loop is skipped.
    if (length(samplesize) == 0) {
      warning(paste("未在 N.xlsx 中找到 ID 为", expo_id, "的样本量。"))
      next
    }
    
    # Create a folder to store all the results of this exposure
    expo_result_dir <- file.path("C:/Users/51195/Desktop/MRresult", expo_id)
    if (!dir.exists(expo_result_dir)) dir.create(expo_result_dir, recursive = TRUE)
    
    # Check and get the column names of the exposed files
    if (grepl("expo_finngen_R", expo_file)) {
      snp_col <- "rsids"
      beta_col <- "beta"
      se_col <- "sebeta"
      effect_allele_col <- "alt"
      other_allele_col <- "ref"
      eaf_col <- "af_alt"
      pval_col <- "pval"
    } else {
      col_data <- read.table(expo_file, header = TRUE, sep = "\t", nrows = 5)
      col_names <- colnames(col_data)
      # Adjustment of logic for determining SNP column names
      if ("rsid" %in% col_names) {
        snp_col <- "rsid"
      } else if ("rs_id" %in% col_names) {
        snp_col <- "rs_id"
      } else {
        snp_col <- col_names[which(col_names %in% snp_cols)[1]] # Modify the logic
      }
      pval_col <- check_cols(expo_file, pval_cols)
      beta_col <- check_cols(expo_file, beta_cols)
      se_col <- check_cols(expo_file, se_cols)
      effect_allele_col <- check_cols(expo_file, effect_allele_cols)
      other_allele_col <- check_cols(expo_file, other_allele_cols)
      eaf_col <- check_cols(expo_file, eaf_cols)
    }
    
    # Reading exposed data
    expo_rt <- read_exposure_data(
      filename = expo_file,
      sep = "\t",
      snp_col = snp_col,
      beta_col = beta_col,
      se_col = se_col,
      effect_allele_col = effect_allele_col,
      other_allele_col = other_allele_col,
      eaf_col = eaf_col,
      pval_col = pval_col
      # ,samplesize_col = samplesize_col
    )
    
    if (nrow(expo_rt) == 0) {
      cat(paste("暴露数据", expo_id, "为空。\n"))
      next
    }
    
    # Setting the sample size
    expo_rt$samplesize.exposure <- samplesize
    
    # Initialize the result data frame for this exposure
    result <- data.frame()
    
    # Traversing the ending files
    for (outcome_file in outcome_files) {
      tryCatch({
        # Get the file name
        filename_base <- basename(outcome_file)
        
        # The file name after removing the extension is used as the ID
        outcome_id <- sub("\\.txt$|\\.tsv$|\\.tsv\\.gz$|\\.h\\.tsv\\.gz$|\\.gz$", "", filename_base)
        
        # Check and get the column names of the final file
        if (grepl("finn", outcome_file)) {
          snp_col_outcome <- "rsids"
          beta_col_outcome <- "beta"
          se_col_outcome <- "sebeta"
          effect_allele_col_outcome <- "alt"
          other_allele_col_outcome <- "ref"
          eaf_col_outcome <- "af_alt"
          pval_col_outcome <- "pval"
        } else {
          col_data_outcome <- read.table(outcome_file, header = TRUE, sep = "\t", nrows = 5)
          col_names_outcome <- colnames(col_data_outcome)
          # Adjustment of logic for determining SNP column names
          if ("rsid" %in% col_names_outcome) {
            snp_col_outcome <- "rsid"
          } else if ("rs_id" %in% col_names_outcome) {
            snp_col_outcome <- "rs_id"
          } else {
            snp_col_outcome <- col_names_outcome[which(col_names_outcome %in% snp_cols)[1]] # Modify the logic here
          }
          pval_col_outcome <- check_cols(outcome_file, pval_cols)
          beta_col_outcome <- check_cols(outcome_file, beta_cols)
          se_col_outcome <- check_cols(outcome_file, se_cols)
          effect_allele_col_outcome <- check_cols(outcome_file, effect_allele_cols)
          other_allele_col_outcome <- check_cols(outcome_file, other_allele_cols)
          eaf_col_outcome <- check_cols(outcome_file, eaf_cols)
          samplesize_col_outcome <- check_cols(outcome_file, samplesize_cols)
        }
        
        # Read the outcome data
        outc_rt <- read_outcome_data(
          snps = expo_rt$SNP,
          filename = outcome_file,
          sep = "\t",
          snp_col = snp_col_outcome,
          beta_col = beta_col_outcome,
          se_col = se_col_outcome,
          effect_allele_col = effect_allele_col_outcome,
          other_allele_col = other_allele_col_outcome,
          eaf_col = eaf_col_outcome,
          pval_col = pval_col_outcome,
          samplesize_col = samplesize_col_outcome
        )
        
        if (nrow(outc_rt) == 0) {
          cat(paste("结局数据", outcome_id, "为空或没有匹配的 SNP。\n"))
          next
        }
        
        # Data coordination
        harm_rt <- harmonise_data(
          exposure_dat = expo_rt,
          outcome_dat = outc_rt,
          action = 2
        )
        
        if (is.null(harm_rt) || nrow(harm_rt) == 0) {
          cat(paste("暴露", expo_id, "和结局", outcome_id, "没有协调后的数据。\n"))
          next
        }
        
        # Calculate R2 and F values
        harm_rt$R2 <- (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure)) / 
          (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) + 
             2 * harm_rt$samplesize.exposure * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) * harm_rt$se.exposure^2)
        harm_rt$f <- harm_rt$R2 * (harm_rt$samplesize.exposure - 2) / (1 - harm_rt$R2)
        harm_rt$meanf <- mean(harm_rt$f, na.rm = TRUE)
        harm_rt <- harm_rt[harm_rt$f > 10, ]
        
        if (nrow(harm_rt) < 3) {
          cat(paste("在过滤后，暴露", expo_id, "和结局", outcome_id, "的 SNP 数量不足。\n"))
          # Write error messages to NORESULT.txt
          write(paste("在过滤后，暴露", expo_id, "和结局", outcome_id, "的 SNP 数量不足。"), file = file.path(noresult_dir, "NORESULT.txt"), append = TRUE)
          next
        }
        
        # Preparing GSMR data
        gsmr_data <- data.frame(
          SNP = harm_rt$SNP,
          a1 = harm_rt$effect_allele.exposure,
          a2 = harm_rt$other_allele.exposure,
          a1_freq = harm_rt$eaf.exposure,
          bzx_n = harm_rt$samplesize.exposure,
          bzx = harm_rt$beta.exposure,
          bzx_se = harm_rt$se.exposure,
          bzx_pval = harm_rt$pval.exposure,
          bzy = harm_rt$beta.outcome,
          bzy_n = harm_rt$samplesize.outcome,
          bzy_se = harm_rt$se.outcome,
          bzy_pval = harm_rt$pval.outcome
        )
        
        # Preparation of LD matrix
        ldrho <- diag(nrow(gsmr_data))
        colnames(ldrho) <- rownames(ldrho) <- snp_coeff_id <- as.character(gsmr_data$SNP)
        
        # GSMR parameters
        bzx <- gsmr_data$bzx
        bzx_se <- gsmr_data$bzx_se
        bzx_pval <- gsmr_data$bzx_pval
        bzy <- gsmr_data$bzy
        bzy_se <- gsmr_data$bzy_se
        bzy_pval <- gsmr_data$bzy_pval
        n_ref <- 7703
        gwas_thresh <- thresh
        multi_snps_heidi_thresh <- 0.01
        nsnps_thresh <- 3
        heidi_outlier_flag <- TRUE
        ld_r2_thresh <- 0.05
        ld_fdr_thresh <- 0.05
        gsmr2_beta <- 1
        
        # Run GSMR analysis
        gsmr_results <- gsmr(
          bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval,
          ldrho, snp_coeff_id, n_ref, heidi_outlier_flag, gwas_thresh,
          multi_snps_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr2_beta
        )
        
        if (is.null(gsmr_results)) {
          next
        }
        
        # Extract results
        beta <- gsmr_results[["bxy"]]
        beta_se <- gsmr_results[["bxy_se"]]
        beta_lci <- beta - 1.96 * beta_se
        beta_uci <- beta + 1.96 * beta_se
        or <- exp(beta)
        or_lci <- exp(beta_lci)
        or_uci <- exp(beta_uci)
        pvalue <- gsmr_results[["bxy_pval"]]
        
        # Save the results
        result_df <- data.frame(
          exposure = expo_id,
          outcome = outcome_id,
          beta = beta,
          beta_lci = beta_lci,
          beta_uci = beta_uci,
          or = or,
          or_lci = or_lci,
          or_uci = or_uci,
          pvalue = pvalue
        )
        
        # Add the results to the exposed results dataframe
        result <- rbind(result, result_df)
        
        # Save each result as a separate file
        output_file <- file.path(expo_result_dir, paste0(expo_id, "&", outcome_id, ".xlsx"))
        write.xlsx(result_df, file = output_file, rowNames = FALSE)
        
        cat(paste("GSMR 分析完成，暴露", expo_id, "和结局", outcome_id, "\n"))
      }, error = function(e) {
        cat(paste("在处理暴露", expo_id, "和结局", outcome_id, "时出错：", e$message, "\n"))
        write(paste("在处理暴露", expo_id, "和结局", outcome_id, "时出错：", e$message), file = file.path(noresult_dir, "NORESULT.txt"), append = TRUE)
      })
    }
    
    # After processing all outcomes for that exposure, merge into a single results file
    all_results_files <- list.files(expo_result_dir, pattern = "*.xlsx", full.names = TRUE)
    if (length(all_results_files) > 0) {
      combined_result <- do.call(rbind, lapply(all_results_files, function(file) {
        df <- read.xlsx(file)
        df
      }))
      
      # Save the merged result
      combined_output_file <- file.path(allresult_dir, paste0(expo_id, ".xlsx"))
      write.xlsx(combined_result, file = combined_output_file, rowNames = FALSE)
      
      cat(paste("所有结果已合并并保存为", expo_id, ".xlsx\n"))
    } else {
      cat(paste("暴露", expo_id, "没有可合并的结果。\n"))
    }
  }, error = function(e) {
    cat(paste("在处理暴露文件", expo_file, "时出错：", e$message, "\n"))
    write(paste("在处理暴露文件", expo_file, "时出错：", e$message), file = file.path(noresult_dir, "NORESULT.txt"), append = TRUE)
  })
}

cat("所有 GSMR 分析和结果保存已完成！\n")