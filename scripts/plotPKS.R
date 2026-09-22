#!/usr/bin/env Rscript
library(circlize)
library(Cairo)
library(dplyr)

# F12: the contig and the island offset used to be written into this script, while
# the pipeline held them as parameters. They are arguments now, so there is one
# source of truth. The interval is 0-based half-open, as the bedgraph is.
args <- commandArgs(trailingOnly = TRUE)
coverage_file <- args[1]
cytoband_file <- args[2]
output_pdf    <- args[3]
contig        <- if (length(args) >= 4) args[4] else "NC_017628.1"
island_start  <- if (length(args) >= 5) as.numeric(args[5]) else 2193827
island_end    <- if (length(args) >= 6) as.numeric(args[6]) else 2244594


# Read cytoband data
cytoband.df <- read.csv(cytoband_file, sep = "\t")
cytoband.df$start <- as.numeric(cytoband.df$start)
cytoband.df$end <- as.numeric(cytoband.df$end)
cytoband.df <- cytoband.df[, c('chrom', 'start', 'end', 'Name', 'gieStain')]
# Gene coordinates in the cytoband file are island-relative; shift them onto the
# reference using the same island start the rest of the pipeline uses.
cytoband.df <- mutate(cytoband.df, start = as.integer(start) + island_start)
cytoband.df <- mutate(cytoband.df, end   = as.integer(end)   + island_start)

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

# F12: keep every bin that INTERSECTS the island and clip it to the island, rather
# than keeping only bins wholly inside it. A bedgraph bin that straddles either edge
# -- starting before the island and ending inside, or the reverse -- was dropped
# entirely, so coverage at the first and last bins never reached the plot. Those are
# exactly the boundaries the mobility context cares about.
coverage <- coverage[
  (coverage$chr == contig) &
  (coverage$end   > island_start) &
  (coverage$start < island_end),
]
if (nrow(coverage) > 0) {
  coverage$start <- pmax(coverage$start, island_start)
  coverage$end   <- pmin(coverage$end,   island_end)
  coverage <- coverage[coverage$end > coverage$start, ]
}
message(sprintf("plotPKS.R: %d coverage bins in %s:%d-%d",
                nrow(coverage), contig, island_start, island_end))

 

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

# F13: the track is raw depth now, so fixed 1e5 steps -- sized for RPKM values --
# would leave a single tick at zero. Let the breaks follow the data.
max_y_value <- max(coverage$value)
y_axis_breaks <- pretty(c(0, max(max_y_value, 1)), n = 4)

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
