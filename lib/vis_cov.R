
argv <- commandArgs(TRUE)

dir <- argv[1]
covdata <- argv[2]
phenodata <- argv[3]
out <- argv[4]


source(paste0(dir, "/vis_utils.R"))
if (!require("APAVplot", quietly = TRUE))
    stop("Please install R package \"APAVplot\" first")
library(APAVplot)

cov_data <- read.table(covdata, header = T, comment.char = "#", sep = "\t", check.names = F)

myregion <- cov_data[, 1:4]
rownames(myregion) <- cov_data$Annotation
mycov <- cov_data[, 6:ncol(cov_data)]
rownames(mycov) <- cov_data$Annotation

mycov <- check_region(mycov)

if(phenodata == "NULL"){
  my_cov <- get_cov_obj(mycov, 
                        region_info = myregion)
}else{
  mypheno <- read.table(phenodata, header = T, comment.char = "#", sep = "\t", quote = "")
  rownames(mypheno) <- mypheno[, 1]
  mypheno <- mypheno[, -1, drop = F]
  
  my_cov <- get_cov_obj(mycov, 
                        region_info = myregion, 
                        pheno_info = mypheno)
}

pdf(file = paste0(out, ".pdf"), width = as.numeric(argv[65]), height = as.numeric(argv[66]))
suppressMessages(cov_heatmap(my_cov, 
            add_pheno_info = names(my_cov@sample),
	    #add_region_info = "Chr",
            
            cov_colors = labelsnull(argv[5]),
            region_info_color_list = colorlist(argv[6]),
            pheno_info_color_list = colorlist(argv[7]),
            
            border = tobool(argv[8]),
            
            cluster_rows = tobool(argv[9]),
            clustering_distance_rows = argv[10], 
            clustering_method_rows = argv[11],
            row_dend_side = argv[12],
            row_dend_width = grid::unit(as.numeric(argv[13]), "mm"),
            row_sorted = toc(argv[14]),
            
            show_row_names = tobool(argv[15]),
            row_names_side = argv[16],
            row_names_size = as.numeric(argv[17]),
            row_names_rot = as.numeric(argv[18]),
            
            cluster_columns = tobool(argv[19]),
            clustering_distance_columns = argv[20],
            clustering_method_columns = argv[21],
            column_dend_side = argv[22],
            column_dend_height = grid::unit(as.numeric(argv[23]), "mm"),
            column_sorted = toc(argv[24]),
            
            show_column_names = tobool(argv[25]),
            column_names_side = argv[26],
            column_names_size = as.numeric(argv[27]),
            column_names_rot = as.numeric(argv[28]),
            
            anno_param_row_pheno = list(show = TF(argv[29], T), 
                                       width = as.numeric(argv[30]),
                                       border = TF(argv[31], F),
                                       name_size = numnull(argv[32]),
                                       name_rot = as.numeric(argv[33]),
                                       name_side = argv[34]),
            anno_param_column_region = list(show = TF(argv[35], T), 
                                            height = as.numeric(argv[36]),
                                            border = TF(argv[37], F),
                                            name_size = numnull(argv[38]),
                                            name_rot = as.numeric(argv[39]),
                                            name_side = argv[40]),
            anno_param_row_stat = list(show = TF(argv[41], F), 
                                       width = as.numeric(argv[42]),
                                       border = TF(argv[43], F),
                                       title = strnull(argv[44]),
                                       title_size = as.numeric(argv[45]),
                                       title_side = argv[46],
                                       title_rot = as.numeric(argv[47]),
                                       axis_side = argv[48],
                                       axis_labels_size = as.numeric(argv[49])),
            anno_param_column_stat = list(show = TF(argv[50], F), 
                                          hight = as.numeric(argv[51]),
                                          border = TF(argv[52], F),
                                          title = strnull(argv[53]),
                                          title_size = as.numeric(argv[54]),
                                          title_side = argv[55],
                                          title_rot = as.numeric(argv[56]),
                                          axis_side = argv[57],
                                          axis_labels_size = as.numeric(argv[58])),
            
            legend_side = argv[59],
            legend_title = strnull(argv[60]),
            legend_title_size = numnull(argv[61]),
            legend_text_size = numnull(argv[62]),
            legend_grid_size = grid::unit(as.numeric(argv[63]), "mm"),
            
            use_raster = tobool(argv[64])
            
            ))
dev.off()






