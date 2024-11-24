
argv <- commandArgs(TRUE)

dir <- argv[1]
data <- argv[2]
out <- argv[3]

source(paste0(dir, "/vis_utils.R"))
if (!require("APAVplot", quietly = TRUE))
    stop("Please install R package \"APAVplot\" first")
library(APAVplot)

sim_data <- read.table(data, header = T, sep = "\t")

pdf(file = paste0(out, ".pdf"), width = as.numeric(argv[35]), height = as.numeric(argv[36]))
plot_sim(sim_data,
         chart_type = argv[4],
         data_type = argv[5],
         
         x_title = strnull(argv[6]),
         y_title = strnull(argv[7]),
         x_title_size = numnull(argv[8]),
         y_title_size = numnull(argv[9]),
         x_breaks = breaksnull(argv[10]),
         y_breaks = breaksnull(argv[11]),
         x_text_size = numnull(argv[12]),
         y_text_size = numnull(argv[13]),
         
         legend_show = !tobool(argv[14]),
         legend_side = argv[15],
         legend_title = strnull(argv[16]),
         legend_title_size = numnull(argv[17]),
         legend_text_size = numnull(argv[18]),
         
         errorbar_width = as.numeric(argv[19]),
         errorbar_size = as.numeric(argv[20]),
         errorbar_color = strnull(argv[21]),
         errorbar_alpha = as.numeric(argv[22]),
         errorbar_point_size = as.numeric(argv[23]),
         errorbar_point_color = strnull(argv[24]),
         
         jitter_width = as.numeric(argv[25]),
         jitter_size = as.numeric(argv[26]),
         jitter_color = strnull(argv[27]),
         jitter_alpha = as.numeric(argv[28]),
         jitter_point_size = as.numeric(argv[279]),
         jitter_point_color = strnull(argv[30]),
         
         path_size = as.numeric(argv[31]),
         path_color = labelsnull(argv[32]),
         ribbon_fill = labelsnull(argv[33]),
         ribbon_alpha = as.numeric(argv[34]))
dev.off()

