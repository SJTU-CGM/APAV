
################### Get arguments

argv <- commandArgs(TRUE)

dir <- argv[1]
pavdata <- argv[2]
phenodata <- argv[3]
out <- argv[4]
command <- argv[5]


################### Import functions

source(paste0(dir, "/vis_utils.R"))
if (!require("APAVplot", quietly = TRUE))
    stop("Please install R package \"APAVplot\" first")
library(APAVplot)

################### Construct a PAV object

pav_data <- read.table(pavdata, header = T, comment.char = "#", sep = "\t", check.names = F)

myregion <- pav_data[, 1:4]
rownames(myregion) <- pav_data$Annotation

if(command == "manhattan"){
  if(is.numeric(myregion$Start)){
    myregion$start_n <- myregion$Start
  }else{
    myregion$start_n <- as.numeric(unlist(strsplit(myregion$Start, ','))[1])
  }
}

mypav <- pav_data[, 6:ncol(pav_data)]
rownames(mypav) <- pav_data$Annotation

mypav <- check_region(mypav)

if(phenodata == "NULL"){
  my_pav <- get_pav_obj(mypav,
                        region_info = myregion)
}else{
  mypheno <- read.table(phenodata, header=T, comment.char = "#", sep = "\t", quote = "" )
  rownames(mypheno) <- mypheno[, 1]
  mypheno <- mypheno[, -1, drop = F]

  my_pav <- get_pav_obj(mypav,
                        pheno_info = mypheno,
                        region_info = myregion)
}

################### 

if(command == "stat"){
  
  mypheno <- pheno_stat(my_pav, 
                      pheno_name = names(my_pav@sample$pheno),
                      p_adjust_method = argv[6],
                      pallel = ifelse(argv[7] > 0, T, F),
                      parallel_n = argv[8])
  mypheno <- mypheno[order(mypheno$p_adjusted, mypheno$p_value),]
  
  write.table(mypheno, out, row.names = F, sep = "\t", quote = F)
  
}else if (command == "heatmap"){

  pheno_res <- read.table(argv[6], header = T, sep = "\t")

  pdf(file = paste0(out, ".pdf"), width = as.numeric(argv[7]), height = as.numeric(argv[8]))
  suppressMessages(pheno_heatmap(my_pav, pheno_res,
		add_region_info = "Chr",
		
		p_threshold = as.numeric(argv[9]),
                adjust_p = tobool(argv[10]),
                only_show_significant = tobool(argv[11]),
                flip = tobool(argv[12]),

		p_colors = labelsnull(argv[13]),
                na_col = argv[14],
                cell_border_color = argv[15],
                region_info_color_list = colorlist(argv[16]),

		cluster_rows = tobool(argv[17]),
            	clustering_distance_rows = argv[18],
            	clustering_method_rows = argv[19],
            	row_dend_side = argv[20],
            	row_dend_width = grid::unit(as.numeric(argv[21]), "mm"),
            	row_sorted = toc(argv[22]),

            	show_row_names = !tobool(argv[23]),
            	row_names_side = argv[24],
            	row_names_size = as.numeric(argv[25]),
            	row_names_rot = as.numeric(argv[26]),

            	cluster_columns = tobool(argv[27]),
            	clustering_distance_columns = argv[28],
            	clustering_method_columns = argv[29],
            	column_dend_side = argv[30],
            	column_dend_height = grid::unit(as.numeric(argv[31]), "mm"),
            	column_sorted = toc(argv[32]),

            	show_column_names = !tobool(argv[33]),
            	column_names_side = argv[34],
            	column_names_size = as.numeric(argv[35]),
            	column_names_rot = as.numeric(argv[36]),

            	anno_param_region = list(show = TF(argv[37], T),
                                       width = as.numeric(argv[38]),
                                       border = TF(argv[39], F),
                                       name_size = numnull(argv[40]),
                                       name_rot = as.numeric(argv[41]),
                                       name_side = argv[42]),

		legend_side = argv[43],
            	legend_title = strnull(argv[44]),
            	legend_title_size = numnull(argv[45]),
            	legend_text_size = numnull(argv[46]),
            	legend_grid_size = grid::unit(as.numeric(argv[47]), "mm")

		))
  dev.off()

}else if(command == "block"){
  pheno_res <- read.table(argv[6], header = T, sep = "\t")

  pdf(file = paste0(out, ".pdf"), width = as.numeric(argv[7]), height = as.numeric(argv[8]))
  suppressMessages(pheno_block(my_pav, pheno_res,
		pheno_name = argv[9],

                add_region_info = ifelse(tobool(argv[11]), "p_adjusted", "p_value"),

                p_threshold = as.numeric(argv[10]),
                adjust_p = tobool(argv[11]),
                only_show_significant = tobool(argv[12]),
                flip = tobool(argv[13]),

                per_colors = labelsnull(argv[14]),
                na_col = argv[15],
                cell_border_color = argv[16],
                region_info_color_list = colorlist(argv[17]),

                cluster_rows = tobool(argv[18]),
                clustering_distance_rows = argv[19],
                clustering_method_rows = argv[20],
                row_dend_side = argv[21],
                row_dend_width = grid::unit(as.numeric(argv[22]), "mm"),
                row_sorted = toc(argv[23]),

                show_row_names = !tobool(argv[24]),
                row_names_side = argv[25],
                row_names_size = as.numeric(argv[26]),
                row_names_rot = as.numeric(argv[27]),

                cluster_columns = tobool(argv[28]),
                clustering_distance_columns = argv[29],
                clustering_method_columns = argv[30],
                column_dend_side = argv[31],
                column_dend_height = grid::unit(as.numeric(argv[32]), "mm"),
                column_sorted = toc(argv[33]),

                show_column_names = !tobool(argv[34]),
                column_names_side = argv[35],
                column_names_size = as.numeric(argv[36]),
                column_names_rot = as.numeric(argv[37]),

                anno_param_region = list(show = TF(argv[38], T),
                                       width = as.numeric(argv[39]),
                                       border = TF(argv[40], F),
                                       name_size = numnull(argv[41]),
                                       name_rot = as.numeric(argv[42]),
                                       name_side = argv[43]),

                legend_side = argv[44],
                legend_title = strnull(argv[45]),
                legend_title_size = numnull(argv[46]),
                legend_text_size = numnull(argv[47]),
                legend_grid_size = grid::unit(as.numeric(argv[48]), "mm")
		))
  dev.off()

}else if(command == "bar"){

  pdf(file = paste0(out, ".pdf"), width = as.numeric(argv[8]), height = as.numeric(argv[9]))

  suppressWarnings(print(pheno_bar(my_pav, 
	    pheno_name = argv[6],
	    region_name = argv[7],

	    pav_colors = unlist(strsplit(argv[10], ",")),
            bar_width = as.numeric(argv[11]),
            x_text_size = numnull(argv[12]),
            x_title_size = numnull(argv[13]),
            y_text_size = numnull(argv[14]),
            y_title_size = numnull(argv[15]),
            legend_side = argv[16],
            legend_title_size = numnull(argv[17]),
            legend_text_size = numnull(argv[18])
	    )))

  dev.off()
}else if(command == "violin"){

  pdf(file = paste0(out, ".pdf"), width = as.numeric(argv[8]), height = as.numeric(argv[9]))

  suppressWarnings(print(pheno_violin(my_pav,
            pheno_name = argv[6],
            region_name = argv[7],

            pav_colors = unlist(strsplit(argv[10], ",")),
            
            x_text_size = numnull(argv[11]),
            x_title_size = numnull(argv[12]),
            y_text_size = numnull(argv[13]),
            y_title_size = numnull(argv[14]),
            legend_side = argv[15],
            legend_title_size = numnull(argv[16]),
            legend_text_size = numnull(argv[17])
            )))

  dev.off()

}else if(command == "manhattan"){

  pheno_res <- read.table(argv[6], header = T, sep = "\t")

  pdf(file = paste0(out, ".pdf"), width = as.numeric(argv[8]), height = as.numeric(argv[9]))
  
  print(pheno_manhattan(my_pav,
                  pheno_stat_res = pheno_res,
                  pheno_name = argv[7],
                  chr = "Chr",
                  bp = "start_n",
                  adjust_p = tobool(argv[10]),
                  
                  highlight_top_n = as.numeric(argv[11]),
                  highlight_text_size = as.numeric(argv[12]),
                  point_size = as.numeric(argv[13]),

                  x_text_size = numnull(argv[14]),
		  x_text_angle = as.numeric(argv[15]),
                  x_title_size = numnull(argv[16]),
                  y_text_size = numnull(argv[17]),
                  y_title_size = numnull(argv[18])))


  dev.off()

}




