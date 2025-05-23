library(TwoSampleMR)
library(dplyr)
library(openxlsx)
library(writexl)
library(ggplot2)
library(flextable)
library(officer)
library(igraph)
library(ggraph)
library(grid)  # For arrow()

run_mr_validation <- function(
    input_dir,
    output_dir = "MR_OUTPUT"
) {
  # 1. Set the output directory
  validated_dir <- file.path(output_dir, "VALIDATED")
  plot_dir <- file.path(output_dir, "PLOTS")
  report_dir <- file.path(output_dir, "REPORTS")
  total_summary_file <- file.path(output_dir, "All_Validated_Summary.xlsx")
  network_plot_file <- file.path(output_dir, "Causal_Network.png")
  
  # 2. Create Output Directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(validated_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 3. Get all .xlsx files
  files <- list.files(input_dir, pattern = "\\.xlsx$", full.names = TRUE)
  
  total_summary <- data.frame()
  
  for (file in files) {
    dat <- read.xlsx(file)
    if (!("id.exposure" %in% colnames(dat))) {
      cat("File", basename(file), "does not contain 'id.exposure'. Skipping!\n")
      next
    }
    
    cat("Processing file:", basename(file), "\n")
    cat("Initial rows:", nrow(dat), "\n")
    
    # 4. Data Filtering
    dat <- dat %>%
      filter(method %in% c("MR Egger", "Inverse variance weighted")) %>%
      filter(!is.na(pval))
    
    cat("After method & NA filtering, rows:", nrow(dat), "\n")
    
    # Calculate FDR (reference only, no filtering)
    dat <- dat %>%
      group_by(id.exposure, id.outcome) %>%
      mutate(pval.FDR = p.adjust(pval, method = "fdr")) %>%
      ungroup()
    
    # Core screening
    validated <- dat %>%
      group_by(id.exposure, id.outcome) %>%
      filter(n() >= 1) %>%
      filter(any(pval < 0.05)) %>%
      filter(all(nsnp >= 10)) %>%
      filter(if(n() == 1) TRUE else (all(b < 0) | all(b > 0))) %>%
      ungroup()
    
    cat("After filtering, rows:", nrow(validated), "\n")
    
    # 5. If you pass the screening
    if (nrow(validated) > 0) {
      # Save the results
      write.xlsx(validated, file.path(validated_dir, basename(file)))
      
      # Generate a volcano plot
      p <- ggplot(validated, aes(x = b, y = -log10(pval), color = method)) +
        geom_point(alpha = 0.7) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
        labs(title = paste0("Volcano Plot: ", unique(validated$id.exposure)),
             x = "Effect Size (beta)", y = "-log10(P-value)") +
        theme_bw()
      
      plot_path <- file.path(plot_dir, paste0(gsub("\\.xlsx$", "", basename(file)), "_volcano.png"))
      ggsave(plot_path, p, width = 7, height = 5)
      
      # Word Report
      doc <- read_docx()
      exp_str <- paste(unique(validated$id.exposure), collapse = ", ")
      out_str <- paste(unique(validated$id.outcome), collapse = ", ")
      
      doc <- doc %>%
        body_add_par(paste("MR Validation Report -", basename(file)), style = "heading 1") %>%
        body_add_par(paste("Exposure:", exp_str), style = "Normal") %>%
        body_add_par(paste("Outcome:", out_str), style = "Normal") %>%
        body_add_par("Volcano Plot:", style = "heading 2") %>%
        body_add_img(src = plot_path, width = 5, height = 4) %>%
        body_add_par("Validation Results Table:", style = "heading 2") %>%
        body_add_flextable(flextable(validated))
      
      save_docx <- file.path(report_dir, paste0(gsub("\\.xlsx$", "", basename(file)), "_report.docx"))
      print(doc, target = save_docx)
      
      total_summary <- bind_rows(total_summary, validated)
    } else {
      cat("File", basename(file), "has no records that meet filtering criteria.\n")
    }
  }
  
  # 6. Save Summary
  write.xlsx(total_summary, total_summary_file)
  
  # 7. Generate causal network diagram
  if (nrow(total_summary) > 0) {
    edges <- total_summary %>%
      select(from = id.exposure, to = id.outcome, b) %>%
      mutate(weight = abs(b)) %>%
      filter(!is.na(from) & !is.na(to))
    
    g <- graph_from_data_frame(edges, directed = TRUE)
    png(network_plot_file, width = 900, height = 700)
    
    # If the number of nodes is <= 2, use circle, otherwise use fr
    layout_choice <- if (length(V(g)) <= 2) "circle" else "fr"
    
    net_plot <- ggraph(g, layout = layout_choice) +
      # Keep the original color gradient with effect size + edge width changes with weight
      geom_edge_link(
        aes(width = weight, color = weight), 
        alpha = 0.6,
        arrow = arrow(length = unit(4, "mm"), type = "closed")  # Add arrows only
      ) +
      scale_edge_color_gradient2(
        low = "blue", high = "red", mid = "gray",
        midpoint = median(edges$weight)
      ) +
      geom_node_point(size = 6, color = "steelblue") +
      geom_node_text(aes(label = name), repel = TRUE, size = 5) +
      theme_void()
    
    print(net_plot)
    dev.off()
  }
  
  cat("MR validation process completed! Output directory:", output_dir, "\n")
}

# Run the function
run_mr_validation(
  input_dir = "C:\\Users\\51195\\Desktop\\test",
  output_dir = "C:\\Users\\51195\\Desktop\\result"
)
