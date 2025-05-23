setwd("C:\\Users\\51195\\Desktop")

# Try to specify the file encoding as UTF-8
mydata <- read.table("input.txt", header = TRUE, sep = "\t", fileEncoding = "UTF-8")

# Processing Data
mydata$pvalue <- ifelse(mydata$pvalue < 0.001, "<0.001", sprintf("%.4f", mydata$pvalue))
mydata$` ` <- paste(rep(" ", 30), collapse = " ")
mydata$'OR(95%CI)' <- ifelse(is.na(mydata$or), "", sprintf('%.4f(%.4f to %.4f)', mydata$or, mydata$or_lci, mydata$or_uci))
mydata[is.na(mydata)] <- " "

# Make sure the data is numeric
mydata$or <- as.numeric(mydata$or)
mydata$or_lci <- as.numeric(mydata$or_lci)
mydata$or_uci <- as.numeric(mydata$or_uci)

mydata$"P value" <- mydata$pvalue

# Set the theme of the forest plot
tm <- forest_theme(base_size = 10,
                   ci_pch = 20,
                   ci_col = "#4575b4",
                   ci_lty = 1,
                   ci_lwd = 2.3,
                   ci_Theight = 0.2, 
                   refline_gp = gpar(lwd = 1.5, lty = "dashed", col = "red"),
                   summary_fill = "#4575b4",
                   summary_col = "#4575b4",
                   footnote_gp = gpar(cex = 1.1, fontface = "italic", col = "blue"))

# Draw a forest plot
forest(mydata[, c(1:2, 7:9)],
       est = mydata$or,
       lower = mydata$or_lci,
       upper = mydata$or_uci,
       sizes = 0.6,
       ci_column = 3,
       ref_line = 1,
       xlim = c(0, 2),
       ticks_at = c(0, 1, 2),
       xlab = "GSMR effect size, OR",
       theme = tm)