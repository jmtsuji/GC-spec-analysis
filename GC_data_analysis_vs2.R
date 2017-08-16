# Gas chromatograph data analysis template
# Created Aug 15, 2017 by Jackson Tsuji (Neufeld Lab PhD student)
# Description: Imports and processes GC data. Now handles multiple days of data
# Version 2.0.0
# Last updated: Aug 15, 2017

# Notes:
#   - Imported data must match template. Script can print template on request - set print_data_template to TRUE.
#     - "Type" should be either 'Sample' or 'Standard'
#   - This script requires three extra packages to be installed: plyr, dplyr, and ggplot2. You'll get an error if you try to run without these installed.
#   - To install these packages, run:
#       install.packages(c('plyr', 'dplyr', 'ggplot2'), dependencies = TRUE) # Do for each package name

#####################################################
## Script settings: #################################
setwd("/Users/JTsuji/Documents/Research_General/PhD/04f_microcosm/02_methane_oxidation/03_round2a_J2017/script_vs2/") # your working directory where files are stored
GC_data_filename <- "data_test_2_noExtraStd.csv"
print_std_curve <- TRUE # print a PDF of the standard curve?
print_processed_data <- TRUE # print data tables for standards and samples?
force_zero <- TRUE # force the standard curve plots to go through (0,0)?

print_data_template <- FALSE # prints a template file to fill in with your GC data, then makes the script exit.
#####################################################


#####################################################
## Load required packages: ##########################
library(plyr)
library(dplyr)
library(ggplot2)
library(grid)
#####################################################


#####################################################
### Part A: Input data and clean up #################
#####################################################
# Print data template if requested
if (print_data_template == TRUE) {
  template_table <- data.frame('Peak_time_min' = as.character(), 'Timepoint' = as.character(), 'Day' = as.character(), 'Septum_number' = as.character(), 'Type' = as.character(), 'Std_conc_ppm' = as.character(), 'Sample_name' = as.character(), 'Computer_peak_time_min' = as.character(), 'Area' = as.character(), 'Height' = as.character(), 'Notes' = as.character())
  write.table(template_table, file = "GC_data_template.csv", sep = ",", col.names = TRUE, row.names = FALSE)
  
  print("Template data file printed, and saved as GC_data_template.csv. Script exiting... change print_data_template to FALSE and then run rest of script with your filled out template file.")
  stop("(ending script early)")
}

# Import raw plate data
GC_data <- read.table(GC_data_filename, sep = ",", header = TRUE, comment.char = "", stringsAsFactors = FALSE)

# Check that imported data in 'Type' column is okay, or exit.
if (("Sample" %in% unique(GC_data$Type) && "Standard" %in% unique(GC_data$Type)) == FALSE) {
  stop("column 'Type' has inappropriate values. Should only be 'Sample' or 'Standard'. Script is exiting.")
}

# In the future: consider writing some code to compare observed peak time to computational one...


###############################################
### Part B: Process standards #################
###############################################

### First, split the raw data into multiple sub-tables by standard curve
### This is determined by day AND septum number
# Make a new vector for splitting by standards
GC_data$std_code <- paste(GC_data$Day, "_", GC_data$Septum_number, sep = "")

# Split table into sub-tables in a list
GC_data_sep <- lapply(unique(GC_data$std_code), function(code) {filter(GC_data, std_code == code)})
names(GC_data_sep) <- unique(GC_data$std_code)

# Make a function for converting raw area units to ppm based on standards
raw_to_ppm <- function(data_table) {
  # test variable
  # data_table <- GC_data_sep[[4]]
  
  ### Summarize standard curve data
  std_raw <- filter(data_table, Type == "Standard")
  std_grouped <- group_by(std_raw, Day, Std_conc_ppm)
  std_summ <- summarise(std_grouped, Ave_area = mean(Area), StdDev_area = sd(Area))
  
  # Make linear trendline
  if (force_zero == TRUE) {
    trendline <- lm(Ave_area ~ 0 + Std_conc_ppm, data = std_summ) # Forcing through origin: https://stackoverflow.com/a/18947967, accessed 170815
  } else {
    trendline <- lm(Ave_area ~ Std_conc_ppm, data = std_summ)
  }
  
  # Summarize some trendline values
  if (force_zero == TRUE) {
    trendline_coeff <- c(0, coefficients(trendline)) # to make this match what the coefficients look like when not fixed to origin
  } else {
    trendline_coeff <- coefficients(trendline)
  }
  trendline_Rsquared <- summary(trendline)$r.squared
  trendline_summ <- data.frame("Intercept" = unname(trendline_coeff[1]), "Slope" = unname(trendline_coeff[2]), "R_squared" = trendline_Rsquared)
  
  
  ### Summarize unknowns data to ppm
  unk_raw <- filter(data_table, Type == "Sample")
  unk_grouped <- group_by(unk_raw, Day, Sample_name)
  unk_summ <- summarise(unk_grouped, Ave_area = mean(Area), StdDev_area = sd(Area))
  
  # Convert to concentration using standard curve (linear model)
  unk_summ$Ave_concentration_ppm <- (unk_summ$Ave_area - trendline_coeff[1])/trendline_coeff[2]
  
  # Get standard devations
  unk_summ$StdDev_Concentration_ppm <- (unk_summ$StdDev_area / trendline_coeff[2])
  # For calculating std dev like this, see http://www.psychstat.missouristate.edu/introbook/sbk15.htm, accessed ~Jan. 2017
  
  
  # Figure out where to put the trend line equation text on the plot
  min_x <- min(std_summ$Std_conc_ppm)
  max_x <- max(std_summ$Std_conc_ppm)
  min_y <- min(std_summ$Ave_area)
  max_y <- max(std_summ$Ave_area)
  x_coord <- (max_x - min_x) * 0.05 # arbitrary position 5% to the right of the y axis
  y_coord <- (max_y - min_y) * 0.8 # arbitrary position 80% above the x axis
  
  ## Make standards plot without samples
  std_plot <- ggplot() +
    geom_smooth(data = std_summ, aes(x = Std_conc_ppm, y = Ave_area), method = "lm", se = T, colour = "purple") +
    geom_errorbar(data = std_summ, aes(x = Std_conc_ppm, ymin = Ave_area-StdDev_area, ymax = Ave_area+StdDev_area), width = 0, size = 0.8, alpha = 0.8) +
    geom_point(data = std_summ, aes(Std_conc_ppm, Ave_area), alpha = 0.9, size = 3) +
    annotate("text", x = x_coord, y = y_coord, label = paste("y = ", round(trendline_summ[1,2], digits = 5), "x + ", round(trendline_summ[1,1], digits = 5), "\n R^2 = ", round(trendline_summ[1,3], digits = 5), sep = ""), size = 3) +
    theme_bw() +
    # Add theme elements to make the plot look nice
    theme(panel.grid = element_blank(), title = element_text(size = 10), axis.title = element_text(size = 12), 
          strip.text = element_text(size = 10), strip.background = element_rect(fill = "#e6e6e6"),
          panel.border = element_rect(colour = "black", size = 1), panel.spacing.y = unit(3, "mm"), panel.spacing.x = unit(0, "mm"),
          axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 10, colour = "black"),
          axis.ticks = element_line(size = 0.5), axis.line = element_line(colour = "black", size = 0.5),
          legend.text = element_text(size = 7), legend.title = element_blank(),
          legend.key = element_rect(colour = "grey", size = 0.3), legend.key.size = unit(3, "mm"),
          legend.spacing = unit(1, "mm"), legend.box.just = "left", plot.margin = unit(c(2,2,2,2), "mm")) +
    #         See http://stackoverflow.com/questions/1330989/rotating-and-spacing-axis-labels-in-ggplot2 (approx. Jan 27, 2016) for vertical x-axis labels
    # scale_x_continuous(limits = c(0, NA)) +
    # scale_y_continuous(limits = c(0, NA)) +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks() +
    xlab("Methane concentration (ppm)") +
    ylab("Peak area") +
    ggtitle(paste("Day ", unique(data_table$Day), ", Septum ", unique(data_table$Septum_number), sep = ""))
  
  # # Print plot to the screen
  # print(std_plot)
  
  # if (print_std_curve == TRUE) {
  #   # Modify input file name as the output name for the table
  #   std_plot_name <- paste("Std_plot__Day", unique(data_table$Day), "_Septum_", unique(data_table$Septum_number), ".pdf", sep = "")
  #   
  #   # Save the plot
  #   ggsave(std_plot_name, width = 130, height = 130, units = c("mm"))
  # }
  
  # Add on samples to the standard plot
  std_plot_2 <- ggplot() +
    geom_smooth(data = std_summ, aes(x = Std_conc_ppm, y = Ave_area), method = "lm", se = T, colour = "purple") +
    geom_errorbar(data = std_summ, aes(x = Std_conc_ppm, ymin = Ave_area-StdDev_area, ymax = Ave_area+StdDev_area), width = 0, size = 0.8, alpha = 0.8) +
    geom_point(data = std_summ, aes(Std_conc_ppm, Ave_area), alpha = 0.9, size = 3) +
    geom_errorbar(data = unk_summ, aes(x = Ave_concentration_ppm, ymin = Ave_area-StdDev_area, ymax = Ave_area+StdDev_area), width = 0, size = 0.8, alpha = 0.8, colour = "darkcyan") +
    geom_point(data = unk_summ, aes(Ave_concentration_ppm, Ave_area), alpha = 0.9, size = 3, colour = "darkcyan") +
    annotate("text", x = x_coord, y = y_coord, label = paste("y = ", round(trendline_summ[1,2], digits = 5), "x + ", round(trendline_summ[1,1], digits = 5), "\n R^2 = ", round(trendline_summ[1,3], digits = 5), sep = ""), size = 3) +
    theme_bw() +
    # Add theme elements to make the plot look nice
    theme(panel.grid = element_blank(), title = element_text(size = 10), axis.title = element_text(size = 12), 
          strip.text = element_text(size = 10), strip.background = element_rect(fill = "#e6e6e6"),
          panel.border = element_rect(colour = "black", size = 1), panel.spacing.y = unit(3, "mm"), panel.spacing.x = unit(0, "mm"),
          axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 10, colour = "black"),
          axis.ticks = element_line(size = 0.5), axis.line = element_line(colour = "black", size = 0.5),
          legend.text = element_text(size = 7), legend.title = element_blank(),
          legend.key = element_rect(colour = "grey", size = 0.3), legend.key.size = unit(3, "mm"),
          legend.spacing = unit(1, "mm"), legend.box.just = "left", plot.margin = unit(c(2,2,2,2), "mm")) +
    #         See http://stackoverflow.com/questions/1330989/rotating-and-spacing-axis-labels-in-ggplot2 (approx. Jan 27, 2016) for vertical x-axis labels
    # scale_x_continuous(limits = c(0, NA)) +
    # scale_y_continuous(limits = c(0, NA)) +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks() +
    xlab("Methane concentration (ppm)") +
    ylab("Peak area") +
    ggtitle(paste("Day ", unique(data_table$Day), ", Septum ", unique(data_table$Septum_number), sep = ""))
  
  # # Print plot to the screen
  # print(std_plot_2)
    
  # if (print_std_curve == TRUE) {
  #   # Modify input file name as the output name for the table
  #   std_plot_name <- paste("Std_plot__Day", unique(data_table$Day), "_Septum_", unique(data_table$Septum_number), "withSamples.pdf", sep = "")
  #   
  #   # Save the plot
  #   ggsave(std_plot_name, width = 130, height = 130, units = c("mm"))
  # }
  
  output <- list(std_summ, unk_summ, trendline_summ, std_plot, std_plot_2)
  names(output) <- c("std_summ", "unk_summ", "trendline_summ", "std_plot", "std_plot_with_samples")
  return(output)
}
  
# Run the function for all standard sets
GC_data_calc <- lapply(names(GC_data_sep), function(x) {raw_to_ppm(GC_data_sep[[x]])})
names(GC_data_calc) <- names(GC_data_sep)


####Summarizing output
# Summarize the unknowns (samples) data
GC_data_unknowns <- dplyr::bind_rows(lapply(names(GC_data_calc), function(x) {GC_data_calc[[x]][["unk_summ"]]}))

# Write summary table of biological replicates if desired
if (print_processed_data == TRUE) {
  # Modify input file name as the output name for the table
  unk_table_name <- paste(substr(GC_data_filename, 1, nchar(GC_data_filename)-4), "_samples.csv", sep = "")
  
  # Write table
  write.table(GC_data_unknowns, file = unk_table_name, sep = ",", col.names = TRUE, row.names = FALSE)
}

######## CODE NOT ADDED YET #######
# Summarize the trendline equations into a data table
# Columns: Day, Septum_number, Number_of_standards, Trendline_slope, Trendline_y_int, Trendline_R2
# do this as a new feature on another day...

# Summarize standards like with the unknowns above.
##################################

# Make multi-panel standard curve plots, if desired
std_plots_list <- lapply(names(GC_data_calc), function(x) {GC_data_calc[[x]][["std_plot"]]})
std_plot_name <- paste(substr(GC_data_filename, 1, nchar(GC_data_filename)-4), "_std_curves.pdf", sep = "")
if (print_std_curve == TRUE) {
  pdf(file = std_plot_name)
  print(std_plots_list) # See https://stackoverflow.com/a/29834646, accessed 170815
  dev.off()
}

# Do again for the standard curves with samples overlaid
std_plots_list2 <- lapply(names(GC_data_calc), function(x) {GC_data_calc[[x]][["std_plot_with_samples"]]})
std_plot_name2 <- paste(substr(GC_data_filename, 1, nchar(GC_data_filename)-4), "_std_curves_with_samples.pdf", sep = "")
if (print_std_curve == TRUE) {
  pdf(file = std_plot_name2)
  print(std_plots_list2)
  dev.off()
}

#############################################
### Part C: Plot time series ######
#############################################

ts_data <- as.data.frame(GC_data_unknowns) # to work with old code

# Split the sample names by "-" and then make into a data frame
sample_names <- as.data.frame(matrix(unlist(strsplit(ts_data$Sample_name, split = "-")), ncol = 6, byrow = T), stringsAsFactors = F)
colnames(sample_names) <- c("Term", "Experiment_abbrev", "Treatment", "Lake", "Depth", "Replicate")

# Add those sample names back to the main data frame
ts_data$Depth <- sample_names$Depth
ts_data$Treatment <- sample_names$Treatment
ts_data$Replicate <- sample_names$Replicate

# Add units to end of depths
ts_data$Depth <- paste(ts_data$Depth, " m", sep = "")

# Assign levels and nice names to the treatments
# Old and new names in the desired order (top to bottom)
treatment_names_old <- c("O2", "L", "NO3", "Fe", "SO4", "BAC", "AC", "MC")
treatment_names_new <- c("O2", "Light", "Nitrate", "Iron", "Sulfate", "Baseline activity control", "Killed control", "Methanogenesis control")

# Make the treatment into a factor and make the levels in the order of the vectors above
ts_data$Treatment <- factor(ts_data$Treatment, levels = treatment_names_old, ordered = T)

# Change the old names into the human-readable names above
ts_data$Treatment <- mapvalues(ts_data$Treatment, treatment_names_old, treatment_names_new)


#### Plot the data

# Get averages and StdDev
unk_grouped <- dplyr::group_by(ts_data, Depth, Treatment, Day)
unk_summ <- summarise(unk_grouped, Ave = mean(Ave_concentration_ppm), StdDev = sd(Ave_concentration_ppm))

## Write output data to table if desired:
if (print_processed_data == TRUE) {
  # Modify input file name as the output name for the table
  samples_table_name <- paste(substr(GC_data_filename, 1, nchar(GC_data_filename)-4), "_samples_avrg.csv", sep = "")

  # Write table
  write.table(unk_summ, file = samples_table_name, sep = ",", col.names = TRUE, row.names = FALSE)
}
  
## Make a plot
# Load colour-blind palette, from http://jfly.iam.u-tokyo.ac.jp/color/ via http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/ (accessed Jan. 11th, 2016)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

unk_plot <- ggplot(unk_summ, aes(Day, Ave)) +
  geom_errorbar(aes(x = Day, ymin = Ave-StdDev, ymax = Ave+StdDev, colour = Treatment), alpha = 0.8, size = 0.5, width = 0.2) + # colour = "#0D0D0D", 
  geom_line(aes(colour = Treatment), size = 1, alpha = 0.8) +
  geom_point(aes(x = Day, y = Ave, fill = Treatment), shape = 21, alpha = 0.85, size = 4) +
  facet_grid(Depth ~ ., scales = "free") +
  theme_bw() +
  # Add theme elements to make the plot look nice
  theme(panel.grid = element_blank(), axis.title = element_text(size = 14), 
        strip.text = element_text(size = 12), strip.background = element_rect(fill = "#e6e6e6"),
        panel.border = element_rect(colour = "black", size = 2),
        axis.text = element_text(size = 9, colour = "black"),
        axis.ticks = element_line(size = 0.5), axis.line = element_line(colour = "black", size = 0.5),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = "bold"),
        legend.key = element_rect(colour = "transparent", size = 0.5), legend.key.size = unit(5, "mm"),
        legend.spacing = unit(1, "mm"), legend.box.just = "left", plot.margin = unit(c(2,2,2,2), "mm")) +
  scale_fill_manual(values = rev(cbbPalette)) +
  scale_color_manual(values = rev(cbbPalette)) +
  xlab("Time (days)") +
  ylab("Methane concentration (ppm)")

# Print plot to the screen
print(unk_plot)

# Modify input file name as the output name for the table
unk_plot_name <- paste(substr(GC_data_filename, 1, nchar(GC_data_filename)-4), "_ts_plot.pdf", sep = "")

# Save the plot
ggsave(unk_plot_name, width = 160, height = 160, units = c("mm"))



#############################################
### Part D: Keep biological replicates separate and plot again ######
#############################################

## Make a plot using non-summarized data as input
# Load colour-blind palette, from http://jfly.iam.u-tokyo.ac.jp/color/ via http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/ (accessed Jan. 11th, 2016)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

unk_plot_2 <- ggplot(ts_data, aes(Day, Ave_concentration_ppm)) +
  geom_errorbar(aes(x = Day, ymin = Ave_concentration_ppm - StdDev_Concentration_ppm, ymax = Ave_concentration_ppm + StdDev_Concentration_ppm, colour = Treatment), alpha = 0.8, size = 0.5, width = 0.2) + # colour = "#0D0D0D", 
  geom_line(aes(colour = Treatment), size = 1, alpha = 0.85) +
  geom_point(aes(x = Day, y = Ave_concentration_ppm, fill = Treatment), shape = 21, alpha = 0.8, size = 4) +
  facet_grid(Depth ~ Replicate, scales = "free") +
  theme_bw() +
  # Add theme elements to make the plot look nice
  theme(panel.grid = element_blank(), axis.title = element_text(size = 14), 
        strip.text = element_text(size = 12), strip.background = element_rect(fill = "#e6e6e6"),
        panel.border = element_rect(colour = "black", size = 2),
        axis.text = element_text(size = 9, colour = "black"),
        axis.ticks = element_line(size = 0.5), axis.line = element_line(colour = "black", size = 0.5),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = "bold"),
        legend.key = element_rect(colour = "transparent", size = 0.5), legend.key.size = unit(5, "mm"),
        legend.spacing = unit(1, "mm"), legend.box.just = "left", plot.margin = unit(c(2,2,2,2), "mm")) +
  scale_fill_manual(values = rev(cbbPalette)) +
  scale_color_manual(values = rev(cbbPalette)) +
  xlab("Time (days)") +
  ylab("Methane concentration (ppm)")

# Print plot to the screen
print(unk_plot_2)

# Modify input file name as the output name for the table
unk_plot_name <- paste(substr(GC_data_filename, 1, nchar(GC_data_filename)-4), "_ts_plot_ALT.pdf", sep = "")

# Save the plot
ggsave(unk_plot_name, width = 300, height = 160, units = c("mm"))


