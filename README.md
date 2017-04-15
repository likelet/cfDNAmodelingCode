# cfDNAmodelingCode
Rscripts for analysis cfDNA methylation level in multi cancers.
Currently, only code for HCC analysis was released. All analysis started with well-imputed methlyation datamatrix from sequencing platform.
# SessionInfor in R env
```r 
sessionInfo()
R version 3.2.3 (2015-12-10)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

locale:
 [1] LC_CTYPE=zh_CN.UTF-8          LC_NUMERIC=C                  LC_TIME=zh_CN.UTF-8          
 [4] LC_COLLATE=zh_CN.UTF-8        LC_MONETARY=zh_CN.UTF-8       LC_MESSAGES=zh_CN.UTF-8      
 [7] LC_PAPER=zh_CN.UTF-8          LC_NAME=zh_CN.UTF-8           LC_ADDRESS=zh_CN.UTF-8       
[10] LC_TELEPHONE=zh_CN.UTF-8      LC_MEASUREMENT=zh_CN.UTF-8    LC_IDENTIFICATION=zh_CN.UTF-8

attached base packages:
[1] grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] R2HTML_2.3.2           limma_3.26.9           risksetROC_1.0.4      
 [4] MASS_7.3-45            RColorBrewer_1.1-2     scales_0.4.1          
 [7] easyGgplot2_1.0.0.9000 devtools_1.12.0        ggthemes_3.3.0        
[10] ggsci_2.0              pheatmap_1.0.8         varSelRF_0.7-5        
[13] randomForest_4.6-12    papeR_1.0-1            xtable_1.8-2          
[16] car_2.1-3              caret_6.0-73           glmnet_2.0-5          
[19] foreach_1.4.3          pROC_1.8               mice_2.25             
[22] Rcpp_0.12.8            rms_5.0-0              SparseM_1.74          
[25] Hmisc_4.0-0            Formula_1.2-1          lattice_0.20-34       
[28] survival_2.40-1        ggplot2_2.2.0.9000     ROCR_1.0-7            
[31] gplots_3.0.1           Matrix_1.2-7.1        

loaded via a namespace (and not attached):
 [1] splines_3.2.3       RWeka_0.4-29        gtools_3.5.0        assertthat_0.1     
 [5] stats4_3.2.3        latticeExtra_0.6-28 quantreg_5.29       chron_2.3-47       
 [9] digest_0.6.10       minqa_1.2.4         colorspace_1.3-1    sandwich_2.3-4     
[13] htmltools_0.3.5     plyr_1.8.4          gmodels_2.16.2      mvtnorm_1.0-5      
[17] gdata_2.17.0        lme4_1.1-12         MatrixModels_0.4-1  htmlTable_1.7      
[21] tibble_1.2          mgcv_1.8-16         FSelector_0.21      TH.data_1.0-7      
[25] withr_1.0.2         nnet_7.3-12         lazyeval_0.2.0      pbkrtest_0.4-6     
[29] magrittr_1.5        memoise_1.0.0       polspline_1.1.12    nlme_3.1-128       
[33] foreign_0.8-67      RWekajars_3.9.0-1   tools_3.2.3         data.table_1.9.6   
[37] multcomp_1.4-6      stringr_1.1.0       munsell_0.4.3       cluster_2.0.5      
[41] entropy_1.2.1       caTools_1.17.1      nloptr_1.0.4        iterators_1.0.8    
[45] bitops_1.0-6        gtable_0.2.0        ModelMetrics_1.1.0  codetools_0.2-15   
[49] reshape2_1.4.2      gridExtra_2.2.1     zoo_1.7-13          knitr_1.15.1       
[53] KernSmooth_2.23-15  rJava_0.9-8         stringi_1.1.2       rpart_4.1-10       
[57] acepack_1.4.1      `
```
