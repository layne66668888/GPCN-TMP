library(openxlsx)
library(dplyr)
library(ggplot2)
library(writexl)
library(officer)
library(flextable)
library(igraph)
library(ggraph)
library(grid)  # For arrow() function

gsmr_verify_pipeline <- function(result_dir, output_dir) {
  validated_dir <- file.path(output_dir, "VALIDATED")
  plot_dir      <- file.path(output_dir, "PLOTS")
  report_dir    <- file.path(output_dir, "REPORTS")
  summary_file  <- file.path(output_dir, "Summary_Validated.xlsx")
  network_file  <- file.path(output_dir, "Causal_Network.png")
  
  dir.create(validated_dir, showWarnings = FALSE)
  dir.create(plot_dir,      showWarnings = FALSE)
  dir.create(report_dir,    showWarnings = FALSE)
  
  result_files   <- list.files(result_dir, pattern = "\\.xlsx$", full.names = TRUE)
  total_validated <- data.frame()
  
  validate_gsmr <- function(df) {
    df %>%
      mutate(
        logOR        = log(or),
        logOR_lci    = log(or_lci),
        logOR_uci    = log(or_uci),
        significance = case_when(
          pvalue < 5e-8 ~ "Genome-wide significant",
          pvalue < 0.05 ~ "Suggestive",
          TRUE          ~ "Non-significant"
        ),
        stable_CI = (beta_lci * beta_uci > 0),
        direction = ifelse(beta > 0, "Positive", "Negative")
      ) %>%
      filter(pvalue < 0.05, stable_CI)
  }
  
  for (file in result_files) {
    dat <- read.xlsx(file)
    if (!"id.exposure" %in% colnames(dat)) next
    
    validated <- validate_gsmr(dat)
    if (nrow(validated) == 0) next
    
    # Save verification results
    write.xlsx(validated, file.path(validated_dir, basename(file)))
    
    # Volcano image (white background)
    volcano_plot <- ggplot(validated, aes(x = beta, y = -log10(pvalue))) +
      geom_point(color = "tomato", alpha = 0.7, size = 2) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      geom_vline(xintercept = 0,       linetype = "dashed", color = "gray") +
      labs(
        title = paste("GSMR Volcano -", basename(file)),
        x     = "Beta",
        y     = "-log10(P-value)"
      ) +
      theme_bw()
    
    volcano_file <- file.path(
      plot_dir,
      paste0(gsub("\\.xlsx$", "", basename(file)), "_volcano.png")
    )
    ggsave(volcano_file, volcano_plot, width = 7, height = 5)
    
    # Word report (removed forest plot part)
    doc <- read_docx() %>%
      body_add_par(paste("GSMR结果验证报告 -", basename(file)), style = "heading 1") %>%
      body_add_par("火山图：", style = "heading 2") %>%
      body_add_img(volcano_file, width = 5, height = 4) %>%
      body_add_par("验证通过结果表：", style = "heading 2") %>%
      body_add_flextable(flextable(validated))
    
    report_file <- file.path(
      report_dir,
      paste0(gsub("\\.xlsx$", "", basename(file)), "_report.docx")
    )
    print(doc, target = report_file)
    
    total_validated <- bind_rows(total_validated, validated)
  }
  
  # Summary and Network Diagram
  if (nrow(total_validated) > 0) {
    write.xlsx(total_validated, summary_file)
    
    edges <- total_validated %>%
      select(from = id.exposure, to = id.outcome, weight = beta)
    
    g <- graph_from_data_frame(edges, directed = TRUE)
    
    png(network_file, width = 1000, height = 800)
    p <- ggraph(g, layout = "fr", weights = abs(E(g)$weight)) +
      geom_edge_link(
        aes(width = abs(weight), color = weight),
        alpha = 0.7,
        arrow = arrow(length = unit(4, "mm"), type = "closed")
      ) +
      geom_node_point(color = "steelblue", size = 6) +
      geom_node_text(aes(label = name), repel = TRUE, size = 5) +
      scale_edge_color_gradient2(
        low      = "darkblue",
        mid      = "grey",
        high     = "darkred",
        midpoint = 0
      ) +
      theme_void()
    print(p)
    dev.off()
    
    cat("✅ 验证完成：结果保存在 VALIDATED，图保存在 PLOTS，网络图已生成，报告完成。\n")
  } else {
    cat("⚠️ 无结果通过验证。\n")
  }
}

# Example call
gsmr_verify_pipeline(
  result_dir = "C:/Users/51195/Desktop/test",
  output_dir = "C:/Users/51195/Desktop/result"
)
