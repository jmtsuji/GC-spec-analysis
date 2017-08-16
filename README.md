# GC-spec-analysis
A collection of scripts for working with GC/spectrophotometer data. Includes converting from standard curves, and so on.

## GC_data_analysis_vs2.R
A script for working with GC data. Run in a console (e.g., RStudio) and modify the variables in the first few lines of the code.

### Usage:
1. To start, set `print_data_template <- TRUE` and run. Will print a template to add your data into.
**NOTE: sample names MUST correspond to the following syntax, or else the script will break: Term-Experiment_abbrev-Treatment-Lake-Depth-Replicate (i.e., six terms separated by dashes; this is something that should be better generalized for the future).**
2. Set the other variables as you'd like, and run...
