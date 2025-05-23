install.packages("forestploter")
install.packages("grid")

library(grid)
library(forestploter)

setwd("C:\\Users\\51195\\Desktop\\MR28SCI\\3_forest")

mydata=read.table("input.txt",header = T,sep = "\t")
mydata$pval=ifelse(mydata$pval<0.001, "<0.001", sprintf("%.4f", mydata$pval))
mydata$` ` <- paste(rep(" ", 20), collapse = " ")
mydata$`OR (95% CI)` <- ifelse(is.na(mydata$or), "",sprintf("%.4f (%.4f - %.4f)",
                                                            mydata$or, mydata$or_lci95, 
                                                            mydata$or_uci95))
mydata$'P-value'=mydata$pval
mydata$Used_SNPS=ifelse(is.na(mydata$Used_SNPS), "", mydata$Used_SNPS)
mydata$pval=ifelse(is.na(mydata$pval), "", mydata$pval)
mydata$or=ifelse(is.na(mydata$or), "", mydata$or)
mydata$or_lci95=ifelse(is.na(mydata$or_lci95), "", mydata$or_lci95)
mydata$or_uci95=ifelse(is.na(mydata$or_uci95), "", mydata$or_uci95)
mydata$'P-value'=ifelse(is.na(mydata$'P-value'), "", mydata$'P-value')

tm1 <- forest_theme(core=list(fg_params=list(hjust = 0.9, x = 0.9),
                              bg_params=list(fill = c("#edf8e9", "#c7e9c0", "#a1d99b"))),
                    colhead=list(fg_params=list(hjust=0.5, x=0.5)))
tm2 <- forest_theme(core=list(fg_params=list(hjust=c(1, 0, 0, 0.5),
                                             x=c(0.9, 0.1, 0, 0.5)),
                              bg_params=list(fill = c("#f6eff7", "#d0d1e6", "#a6bddb", "#67a9cf"))),
                    colhead=list(fg_params=list(hjust=c(1, 0, 0, 0, 0.5),
                                                x=c(0.9, 0.1, 0, 0, 0.5))))
tm3=forest_theme(ci_Theight = 0.2)


forest(mydata[,c(1:2,7:9)],
       est = as.numeric(mydata$or)
       ,
       lower =as.numeric(mydata$or_lci95)
       , 
       upper = as.numeric(mydata$or_uci95)
       ,
       sizes =0.3,
       ci_column =3 ,
       ref_line = 1,
       xlim = c(0.05, 2),
       theme = tm3)

