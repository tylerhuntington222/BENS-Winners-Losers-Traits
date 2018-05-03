###############################################################################
# WLT_main.R

# PURPOSE: Manages control flow for the entire WLT analysis pipeline from
#          raw data to final deliverables (plots, data products etc.).

# INPUT: Raw data files located in `Data\` directory and source code in 
#        various subdirectories of `Rcode`.

# OUTPUTS: All intermediate data products generated throughout the analysis
#          pipeline, and final plots/data products.

# AUTHOR: Tyler Huntington
###############################################################################

# function for setting/resetting working directory to `Rcode/` dir
resetWorkingDir <- function () {
  tryCatch (
    {
      this.dir <- dirname(parent.frame(2)$ofile)
    },
    error = function(x) {
      this.dir <<- dirname(rstudioapi::getActiveDocumentContext()$path)
    }
  )
  setwd(this.dir)
}
resetWorkingDir()

# generate vector of source code scripts to be run in sequential order
src <- c(
  "avian_scripts/data_cleaning/a_df_construct.R",
  "beetle_scripts/data_cleaning/b_df_construct.R",
  "avian+beetle_scripts.R/a+b_trait_response_models_at_spp_specific_spatial_scales.R"
)

for (s in src) {
  source(s)
  resetWorkingDir()
}

