#!/usr/bin/env Rscript
library(circlize)
library(Cairo)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
coverage_file <- args[1]
cytoband_file <- args[2]
output_pdf <- args[3]


# Read cytoband data
cytoband.df <- read.csv(cytoband_file, sep = "\t")
cytoband.df$start <- as.numeric(cytoband.df$start)
cytoband.df$end <- as.numeric(cytoband.df$end)
cytoband.df <- cytoband.df[, c('chrom', 'start', 'end', 'Name', 'gieStain')]
cytoband.df <- mutate(cytoband.df, start = as.integer(start) +2193827)
cytoband.df <- mutate(cytoband.df, end = as.integer(end) + 2193827)

# Hard fail if the file is missing or empty
if (!file.exists(coverage_file) || file.info(coverage_file)$size == 0) {
  message("plotPKS.R: coverage file is missing or empty: ", coverage_file)
  quit(status = 0)  # 0 = pipeline doesn't fail; use 1 if you want it to fail
}

# Also catch "looks empty" (whitespace only)
first_line <- readLines(coverage_file, n = 1, warn = FALSE)
if (length(first_line) == 0 || grepl("^\\s*$", first_line)) {
  message("plotPKS.R: coverage file has no data lines: ", coverage_file)
  quit(status = 0)
}

# Filter and adjust coverage
coverage <- read.csv(coverage_file, sep = "\t", header = FALSE)
colnames(coverage) <- c('chr', 'start', 'end', 'value')
coverage <- coverage[, c('chr', 'start', 'end', 'value')]

# Subset to region of interest
coverage <- coverage[
  (coverage$chr == 'NC_017628.1') &
  (coverage$start >= min(cytoband.df$start)) &
  (coverage$end <= max(cytoband.df$end)),
]

 

if (nrow(coverage) < 1) {
  coverage <- data.frame(
    chr   = cytoband.df$chrom,
    start = as.integer(cytoband.df$start),
    end   = as.integer(cytoband.df$end),
    value = rep(0, nrow(cytoband.df)),
    stringsAsFactors = FALSE
  )
}

# Extract sample name from coverage filename (optional)
sample_name <- basename(coverage_file)
sample_name <- sub("\\.coverage\\.bedgraph$", "", sample_name)

# Open PDF device
CairoPDF(output_pdf, width = 5, height = 5)

circos.clear()
circos.par(gap.after=3, gap.before=3, ADD = TRUE)
circos.initializeWithIdeogram(cytoband.df)
#title(sample_name,  cex = 1.5, side = "top", adj = 0.5)
#circos.text(0, 1.1, sample_name, facing = "clockwise", cex = 1.2)

circos.labels(cytoband.df$chrom, x = as.numeric(cytoband.df$start), 
              labels = cytoband.df$Name, side = "inside", facing="clockwise", niceFacing=TRUE)

if(length(unique(coverage$value)) == 1){
  circos.genomicTrack(coverage, ylim = c(0, 1), numeric.column = 4,
    panel.fun = function(region, value, ...) {
      circos.genomicLines(region, value, type = "h", numeric.column = 1, col = "#008000", track.height = 0.3)
    })  
} else {
  circos.genomicTrack(coverage, numeric.column = 4,
    panel.fun = function(region, value, ...) {
      circos.genomicLines(region, value, type = "h", numeric.column = 1, col = "#008000", track.height = 0.3)
    })
}

max_y_value <- max(coverage$value) + 100
y_axis_breaks <- seq(0, max_y_value, by = 1e5)

circos.yaxis(
  side = "left",
  at = y_axis_breaks,
  sector.index = unique(coverage$chr),
  labels.cex = 0.1,
  labels.niceFacing = TRUE,
  labels.col = "black"
)
title(sample_name,  cex = 1.5, side = "top", adj = 0.5)
circos.clear()
dev.off()
