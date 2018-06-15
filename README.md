# Rscript clone for NCBI remap API
R script for accessing NCBI remap API
## Purpose [Introduction]
NCBI remap service API is provided in perl script. It is possible to use perl to access the service or simply go to remap web service. 

However, for life scientists who wants to include the process of remapping into their R pipeline, it's a hassle since wrapping perl in R script does not work properly (which I have tried to do previously). So, I made an R script that clones the behavior of perl script provided by NCBI.

## Usage
For now, there is no clean way of executing the function (which is a thing to do for me).

You would need to call the [main] function with variables manually set. If done properly, the function would continue returning "Refreshing..." message until the result is returned from the NCBI remap service.

And, again with limited functionality, the R script should return two different files automatically (remapped report and annotation file). Again, I would need to clean up ways to output files.

## To-do
- Pretty main function to execute with
- Ways of outputting the file
- Functionality that were lost in translation
