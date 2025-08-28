library(rmarkdown)
library(config)

cfg <- config::get(pediatric)  # or config::get("test")

# Create dynamic output name
output_file <- paste0("report_", cfg$output_suffix, ".html")

# Render the R Markdown file
rmarkdown::render(
  input = "01_cohort_analysis.Rmd",
  output_file = output_file,
  params = cfg
)