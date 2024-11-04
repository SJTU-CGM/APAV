
################### Get arguments
argv <- commandArgs(TRUE)

dir <- argv[1]
pavdata <- argv[2]
phenodata <- argv[3]
out <- argv[4]
command <- argv[5]

fig_width <- as.numeric(argv[11])
fig_height <- as.numeric(argv[12])

################### Functions
source(paste0(dir, "/vis_utils.R"))
if (!require("APAVplot", quietly = TRUE))
    stop("Please install R package \"APAVplot\" first")
library(APAVplot)

################### Construct a PAV object
pav_data <- read.table(pavdata, header = T, comment.char = "#", sep = "\t", check.names = F)

myregion <- pav_data[, 1:4]
rownames(myregion) <- pav_data$Annotation
mypav <- pav_data[, 6:ncol(pav_data)]
rownames(mypav) <- pav_data$Annotation

mypav <- check_region(mypav)

if(phenodata == "NULL"){
  my_pav <- get_pav_obj(mypav, 
                        region_info = myregion, 
                        add_softcore = !tobool(argv[6]),
                        add_private = !tobool(argv[7]),
                        softcore_loss_rate = as.numeric(argv[8]),
			use_binomial = tobool(argv[9]),
                        softcore_p_value = as.numeric(argv[10]))
}else{
  mypheno <- read.table(phenodata, header=T, comment.char = "#", sep = "\t", quote = "")
  rownames(mypheno) <- mypheno[,1]
  mypheno <- mypheno[, -1, drop = F]
  
  my_pav <- get_pav_obj(mypav, 
                        pheno_info = mypheno, 
                        region_info = myregion, 
                        add_softcore = !tobool(argv[6]),
                        add_private = !tobool(argv[7]),
                        softcore_loss_rate = as.numeric(argv[8]),
			use_binomial = tobool(argv[9]),
                        softcore_p_value = as.numeric(argv[10]))
}

################### Plotting

pdf(file = paste0(out, ".pdf"), width = fig_width, height = fig_height)

if(command == "heatmap"){
  
  ## Region types
  if(argv[13] == "NULL"){
    mytype <- intersect(unique(my_pav@region$type), c("Core", "Softcore", "Distributed", "Private"))
  }else{
    mytype <- intersect(unique(my_pav@region$type), unlist(strsplit(argv[13], ",")))
  }
  
  ## P and A colors
  mycols <- c(presence = "steelblue", absence = "gray70")
  incols <- unlist(strsplit(argv[14], ","))
  if(!is.na(incols[1])) mycols[1] <- incols[1]
  if(!is.na(incols[2])) mycols[2] <- incols[2]
  
  ## Legend title
  mylgtitle <- list(pav = "PAV", type = "Region")
  inlgtitle <- unlist(strsplit(argv[72], ","))
  if(!is.na(inlgtitle[1])) mylgtitle[[1]] <- inlgtitle[1]
  if(!is.na(inlgtitle[2])) mylgtitle[[2]] <- inlgtitle[2]

   suppressMessages(pav_heatmap(my_pav,
              add_pheno_info = names(my_pav@sample$pheno),
	      #add_region_info = "Chr",

              region_type = mytype,

              pav_colors = mycols,
              type_colors = labelsnull(argv[15]),
              region_info_color_list = colorlist(argv[16]),
              pheno_info_color_list = colorlist(argv[17]),

              border = !tobool(argv[18]),
              split_block = T,
              block_name_size = numnull(argv[19]),
              block_name_rot = as.numeric(argv[20]),

              cluster_rows = tobool(argv[21]),
              clustering_distance_rows = argv[22],
              clustering_method_rows = argv[23],
              row_dend_side = argv[24],
              row_dend_width = grid::unit(as.numeric(argv[25]), "mm"),
              row_sorted = toc(argv[26]),

              show_row_names = tobool(argv[27]),
              row_names_side = argv[28],
              row_names_size = as.numeric(argv[29]),
              row_names_rot = as.numeric(argv[30]),

              cluster_columns = tobool(argv[31]),
              clustering_distance_columns = argv[32],
              clustering_method_columns = argv[33],
              column_dend_side = argv[34],
              column_dend_height = grid::unit(as.numeric(argv[35]), "mm"),
              column_sorted = toc(argv[36]),

              show_column_names = tobool(argv[37]),
              column_names_side = argv[38],
              column_names_size = as.numeric(argv[39]),
              column_names_rot = as.numeric(argv[40]),

              anno_param_row_pheno = list(show = TF(argv[41], T),
                                         width = as.numeric(argv[42]),
                                         border = TF(argv[43], F),
                                         name_size = numnull(argv[44]),
                                         name_rot = as.numeric(argv[45]),
                                         name_side = argv[46]),
              anno_param_column_region = list(show = TF(argv[47], T),
                                              height = as.numeric(argv[48]),
                                              border = TF(argv[49], F),
                                              name_size = numnull(argv[50]),
                                              name_rot = as.numeric(argv[51]),
                                              name_side = argv[52]),
              anno_param_row_stat = list(show = TF(argv[53], T),
                                         width = as.numeric(argv[54]),
                                         border = TF(argv[55], F),
                                         title = strnull(argv[56]),
                                         title_size = as.numeric(argv[57]),
                                         title_side = argv[58],
                                         title_rot = as.numeric(argv[59]),
                                         axis_side = argv[60],
                                         axis_at = NULL,
                                         axis_labels = NULL,
                                         axis_labels_size = as.numeric(argv[61])),
              anno_param_column_stat = list(show = TF(argv[62], T),
                                            hight = as.numeric(argv[63]),
                                            border = TF(argv[64], F),
                                            title = strnull(argv[65]),
                                            title_size = as.numeric(argv[66]),
                                            title_side = argv[67],
                                            title_rot = as.numeric(argv[68]),
                                            axis_side = argv[69],
                                            axis_at = NULL,
                                            axis_labels = NULL,
                                            axis_labels_size = as.numeric(argv[70])),

              legend_side = argv[71],
              legend_title = mylgtitle,
              legend_title_size = numnull(argv[73]),
              legend_text_size = numnull(argv[74]),
              legend_grid_size = grid::unit(as.numeric(argv[75]), "mm"),

              use_raster = tobool(argv[76])
              ))
}else if(command == "hist"){
  
  pav_hist(my_pav,
           show_ring = !tobool(argv[13]),
           ring_pos_x = as.numeric(argv[14]),
           ring_pos_y = as.numeric(argv[15]),
           ring_r = as.numeric(argv[16]),
           ring_label_size = numtona(argv[17]),
           type_colors = labelsnull(argv[18]),
           x_title = strnull(argv[19]),
           x_title_size = numnull(argv[20]),
           y_title = strnull(argv[21]),
           y_title_size = numnull(argv[22]),
           x_breaks = breaksnull(argv[23]),
           x_text_size = numnull(argv[24]),
           y_text_size = numnull(argv[25]))
  
}else if(command == "halfviolin"){
  
  mycolor <- NULL
  if(!is.null(argv[14]) && !is.null(argv[15])) mycolor[[argv[14]]] = labelsnull(argv[15])
  
  pav_halfviolin(my_pav,
                 violin_color = argv[13],
                 add_pheno_info = argv[14],
                 pheno_info_color_list = mycolor,
                 x_text_size  = numnull(argv[16]),
                 y_text_size = numnull(argv[17]),
                 x_title_size = numnull(argv[18]),
                 y_title_size = numnull(argv[19]))
  
}else if(command == "stackbar"){
  
  mycolor <- NULL
  if(!is.null(argv[14]) && !is.null(argv[16])) mycolor[[argv[14]]] = labelsnull(argv[16])
  
  pav_stackbar(my_pav,
               show_relative = tobool(argv[13]),
               add_pheno_info = argv[14],
               
               type_colors = labelsnull(argv[15]),
               pheno_info_color_list = mycolor,
               
               clustering_distance = argv[17],
               clustering_method = argv[18],
               
               bar_width = as.numeric(argv[19]),
               sample_name_size = as.numeric(argv[20]),
               
               legend_side = argv[21],
               legend_title = strnull(argv[22]),
               legend_title_size = numnull(argv[23]),
               legend_text_size = numnull(argv[24]),
               
               dend_width = as.numeric(argv[25]),
               name_width = as.numeric(argv[26]))
  
}else if(command == "cluster"){
  
  mycolor <- NULL
  if(!is.null(argv[15]) && !is.null(argv[16])) mycolor[[argv[15]]] = labelsnull(argv[16])
  
  pav_cluster(my_pav,
              clustering_distance = argv[13],
              clustering_method = argv[14],
              
              add_pheno_info = strnull(argv[15]),
              pheno_info_color_list = mycolor,
              
              sample_name_size = as.numeric(argv[17]),
              mult = as.numeric(argv[18]),
              
              legend_side = argv[19],
              legend_title_size = numnull(argv[20]),
              legend_text_size =numnull(argv[21]))
  
}else if(command == "pca"){
  
  mycolor <- NULL
  if(!is.null(argv[13]) && !is.null(argv[14])) mycolor[[argv[13]]] = labelsnull(argv[14])
  
  pav_pca(my_pav,
          
          add_pheno_info = strnull(argv[13]),
          pheno_info_color_list = mycolor,
          
          axis_text_size = numnull(argv[15]),
          axis_title_size = numnull(argv[16]),
          
          legend_side = argv[17],
          legend_text_size = numnull(argv[18]),
          legend_title_size = numnull(argv[19]))
}

dev.off()
