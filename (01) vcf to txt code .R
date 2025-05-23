# Load the required packages
library(MungeSumstats) 
library(data.table)
library(tidyverse)

# Setting the working directory
setwd("C:/Users/51195/Desktop/exposure/VCF")

# Get all .vcf.gz files
vcf_files <- list.files(pattern = "\\.vcf.gz$")

# Create Error Folder
dir.create("C:/Users/51195/Desktop/MRresult/VCFERROR", showWarnings = FALSE)

# Process each .vcf.gz file
for (file in vcf_files) {
  # Generate .tbi file name and .txt file name
  tbi_file <- gsub("\\.vcf.gz$", ".vcf.gz.tbi", file)
  txt_file <- gsub("\\.vcf.gz$", ".txt", file)
  
  tryCatch({
    # Formatting summary statistics
    mydata <- MungeSumstats::format_sumstats(file, ref_genome="GRCh37")
    
    # Reading Data
    df <- fread(mydata)
    
    # Rename the columns and select the required columns
    df <- df %>% 
      rename(SNP=SNP, A1=A2, A2=A1, freq=FRQ, b=BETA, se=SE, p=P) %>% 
      select(SNP, A1, A2, freq, b, se, p)
    
    # Save data to a .txt file
    fwrite(df, txt_file, row.names = F, sep = "\t", quote = F)
    
    # Delete the .vcf.gz.tbi file and the .vcf.gz file
    file.remove(tbi_file)
    file.remove(file)
    
    # Copy the .txt file to the "C:/Users/51195/Desktop/outcome" folder
    file.copy(txt_file, "C:/Users/51195/Desktop/outcome", overwrite = TRUE)
    
    # Copy the .txt file to the "C:/Users/51195/Desktop/exposure" folder
    file.copy(txt_file, "C:/Users/51195/Desktop/exposure", overwrite = TRUE)
  }, error = function(e) {
    # If an error occurs, write the error information to VCFERROR.txt
    write(file, file = "C:/Users/51195/Desktop/MRresult/VCFERROR/VCFERROR.txt", append = TRUE)
    
    # Delete the erroneous .vcf.gz.tbi file if it exists
    if (file.exists(tbi_file)) {
      file.remove(tbi_file)
    }
  })
}  