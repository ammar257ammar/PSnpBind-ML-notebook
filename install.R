requiredPackages <- c("dplyr", # Data manipulation library
                      "data.table", # Fast processing of large data
                      "protr", # Generating numerical representation of protein sequences
                      "ggplot2", # Plotting Data package
                      "gplots",  # Plotting Data package
                      "RColorBrewer",  # Ready-to-use color palettes for graphics
                      "ggpubr",
                      "caret",
                      "gridExtra",
                      "BiocManager",
                      "ggfortify")

for (pkg in requiredPackages) { 
    if(! pkg %in% row.names(installed.packages())) install.packages(pkg)
}

if(! "pcaMethods" %in% row.names(installed.packages()))
    BiocManager::install("pcaMethods", warning=stop)