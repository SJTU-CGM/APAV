
argv <- commandArgs(TRUE)

dir <- argv[1]
data <- argv[2]
ele_region <- argv[3]
out <- argv[4]

source(paste0(dir, "/vis_utils.R"))
if (!require("APAVplot", quietly = TRUE))
    stop("Please install R package \"APAVplot\" first")
library(APAVplot)

ele_region <- strnull(ele_region)

eledata <- read.table(data, header = T, sep = "\t")

if(argv[5] == "NULL"){
	phenodata <- NULL
}else{
	phenodata <- read.table(argv[5], header = T)
	rownames(phenodata) <- phenodata[, 1]
	phenodata <- phenodata[, -1]
}

if(argv[6] == "NULL"){
	gffdata <- NULL
}else{
	gffdata <- read.table(argv[6], header=F)
}

if(argv[7] == "depth"){
	depth_table <- read.table(argv[8], header = T, sep = "\t")
	pdf(file = paste0(out, ".pdf"), width = as.numeric(argv[9]), height = as.numeric(argv[10]))
	print(plot_ele_depth(depth_data = depth_table,
		       ele_data = eledata,
		       ele_region = ele_region,
		       gff_data = gffdata,
		       pheno_data = phenodata,

		       cluster_samples = tobool(argv[11]),
			clustering_distance = argv[12],
                        clustering_method = argv[13],

			log10 = tobool(argv[14]),

			top_anno_height = as.numeric(argv[15]),
                        left_anno_width = as.numeric(argv[16]),

			ele_color = argv[17],
			depth_colors = labelsnull(argv[18]),
			lowlight_colors = labelsnull(argv[19]),
			gene_colors = labelsnull(argv[20]),
                        pheno_info_color_list = colorlist(argv[21]),

			seq_name_size = numnull(argv[22]),
                        loc_name_size = numnull(argv[23]),

                        show_sample_name = !tobool(argv[24]),
                        sample_name_size = numnull(argv[25]),

			legend_title = argv[26],
                        legend_title_size = numnull(argv[27]),
                        legend_text_size = numnull(argv[28])			       
	))
	dev.off()

}else{
	pdf(file = paste0(out, ".pdf"), width = as.numeric(argv[8]), height = as.numeric(argv[9]))
	if(argv[7] == "cov"){
		print(plot_ele_cov(
			eledata,
			ele_region = ele_region,

			gff_data = gffdata,
	 		pheno_data = phenodata,

			cluster_samples = tobool(argv[10]),
	 		clustering_distance = argv[11],
	 		clustering_method = argv[12],
	 
	 		top_anno_height = as.numeric(argv[13]),
	 		left_anno_width = as.numeric(argv[14]),

	 		ele_color = argv[15],
	 		ele_line_color = argv[16],
	 		cov_colors = labelsnull(argv[17]),
	 		gene_colors = labelsnull(argv[18]),
	 		pheno_info_color_list = colorlist(argv[19]),

	 		seq_name_size = numnull(argv[20]),
	 		loc_name_size = numnull(argv[21]),

	 		show_ele_name = !tobool(argv[22]),
	 		ele_name_size = numnull(argv[23]),
	 		ele_name_rot = as.numeric(argv[24]),

	 		show_sample_name = !tobool(argv[25]),
	 		sample_name_size = numnull(argv[26]),

	 		legend_title_size = numnull(argv[27]),
	 		legend_text_size = numnull(argv[28])
	 	))
	}else{
		print(plot_ele_pav(
			eledata,
			ele_region = ele_region,

         		gff_data = gffdata,
         		pheno_data = phenodata,

			cluster_samples = tobool(argv[10]),
         		clustering_distance = argv[11],
         		clustering_method = argv[12],

         		top_anno_height = as.numeric(argv[13]),
         		left_anno_width = as.numeric(argv[14]),

         		ele_color = argv[15],
         		ele_line_color = argv[16],
         		pav_colors = labelsnull(argv[17]),
         		gene_colors = labelsnull(argv[18]),
         		pheno_info_color_list = colorlist(argv[19]),

         		seq_name_size = numnull(argv[20]),
         		loc_name_size = numnull(argv[21]),

         		show_ele_name = !tobool(argv[22]),
         		ele_name_size = numnull(argv[23]),
         		ele_name_rot = as.numeric(argv[24]),

         		show_sample_name = !tobool(argv[25]),
         		sample_name_size = numnull(argv[26]),

         		legend_title_size = numnull(argv[27]),
         		legend_text_size = numnull(argv[28])
         	))
	}
	dev.off()
}



