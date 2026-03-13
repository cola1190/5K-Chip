# R

library(circlize)

output_dir <- choose.dir(caption = "")
setwd(output_dir) 

input_file <- file.choose()

snp_data <- read.table(input_file, header = FALSE,
                       col.names = c("chr", "start", "end", "pvalue"))

custom_genome <- data.frame(
  chr = c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07",
          "Chr08", "Chr09", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14",
          "Chr15", "Chr16", "Chr17"),
  start = rep(1, 17),
  end = c(25410864, 28589932, 29843988, 24992597, 34903675, 24994454,
          29252404, 25159700, 27821969, 32276254, 32930255, 25426305,
          31173752, 26917991, 41619619, 30522836, 28188000)
)

snp_data$pvalue_log <- -log10(snp_data$pvalue)
snp_data <- snp_data[order(snp_data$chr, snp_data$start), ]

base_name <- tools::file_path_sans_ext(basename(input_file))

plot_circos <- function(file_path, file_type = c("png", "pdf")) {
  file_type <- match.arg(file_type)
  if (file_type == "png") {
    png(file_path, width = 8, height = 8, units = "in", res = 600)
  } else if (file_type == "pdf") {
    pdf(file_path, width = 8, height = 8)
  }
  
  circos.clear()
  circos.par(
    "start.degree" = 90,    
    "clock.wise" = TRUE,    
    "track.height" = 0.15,  
    cell.padding = c(0, 0, 0, 0),
    gap.degree = 1        
  )
  circos.genomicInitialize(custom_genome)
  
  circos.genomicDensity(
    snp_data,
    col = "skyblue",
    track.height = 0.2,
    bg.border = NA
  )
  
  circos.genomicTrackPlotRegion(
    snp_data,
    panel.fun = function(region, value, ...) {
      circos.genomicPoints(region, value, col = "pink", pch = 16, cex = 0.3)
    },
    numeric.column = 5,
    ylim = c(0, max(snp_data$pvalue_log)),
    track.height = 0.3,
    bg.border = NA
  )
  
  circos.clear()
  dev.off()
}

png_path <- file.path(output_dir, paste0(base_name, "_circos.png"))
pdf_path <- file.path(output_dir, paste0(base_name, "_circos.pdf"))

plot_circos(png_path, "png")
plot_circos(pdf_path, "pdf")

cat("Finish！\n")
cat("PNG : ", png_path, "\n")
cat("PDF : ", pdf_path, "\n")
