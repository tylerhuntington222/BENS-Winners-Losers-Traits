################################################################################
# make_all_plots.R

# PURPOSE: Runs all plotting scripts.


# OUTPUTS: All plots exported to `Figures/` directory.
            

# AUTHOR: Tyler Huntington

################################################################################

# Set WLT Analysis directory as working directory for this script.
# Filepaths are relative, so running the below code
# without modification on any machine should set the wd properly.
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
  setwd("..")
}
resetWorkingDir()

scripts <- c(
  "avian_scripts/determine_spatial_scale/a_responses_spatial_compare_top_mod.R",
  "avian_scripts/determine_spatial_scale/a_responses_spatial_compare_top_mod.R",
  "beetle_scripts/determine_spatial_scale/b_responses_spatial_compare_top_mod.R",
  "beetle_scripts/determine_spatial_scale/b_responses_spatial_compare_cand_set.R"
)

for (s in scripts) {
  source(s)
  resetWorkingDir()
}

