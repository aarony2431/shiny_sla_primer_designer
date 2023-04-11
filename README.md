# shiny_sla_primer_designer
The **shiny_sla_primer_designer** is a [ShinyApp](https://shiny.rstudio.com/) hosted for free on https://shinyapps.io used to design oligonucleotide primers for stem-loop assay PCR (SLA), as previously described by [Kramer et al. 2011](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3152947/).

### Usage
The application is located at https://aaronyu.shinyapps.io/shiny-sla-primer-designer/. The app may require several minutes to load.

The user will enter all of the entire targeting (anti-sense, guide) strand sequences without modifications into the suggested input area, separated by spaces or newline characters. The user can press the `Calculate Primers` button to generate a list output of the recommended SLA primer and probe sequences for each inputted sequence. This output may be manually copy and pasted or downloaded as a `.csv` file.

### Notes
- This tool picks the primers which yield the highest number of probes and the overall longest in 50 iterations of optimization.
- The Stem-Loop Assay is a qPCR method designed to detect the Guide Strand from the 3' end.
- **Please note that it is possible to get different primers using the same sequence.** This is due to to a pseudo-random generator for GC content and melting temperature.

## Shiny Reticulate
The primer calculation logic utilizes Python packages without R equivalents. As such, the [reticulate package](https://rstudio.github.io/reticulate/) is required in R to parse native Python code and interface between R and Python objects. A ShinyApp was created and hosted on the web based on [this example](https://github.com/ranikay/shiny-reticulate-app) for making a Shiny app using reticulate.

When the user loads the application, a temporary virtual version of Python is created and installed with all the relevant packages, allowing the application to be self-contained and independent of the user's system.

## Stem-Loop Assay qPCR (SLA)

## Contributions
This tool was originally developed in Python by Sam Beppler (2020), after which it was edited and transferred to https://shinyapps.io as an online tool by Aaron Yu (2022).

The Python code utilizes the C-based oligo tool package [Primer3](https://pypi.org/project/primer3-py/).

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[Apache License 2.0](http://www.apache.org/licenses/)
