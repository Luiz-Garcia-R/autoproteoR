.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "\n",
    crayon::green("autoproteoR "), "loaded successfully!\n",
    "--------------------------------------------------\n",
    "A package for streamlined proteomics data analysis.\n",
    "GitHub: https://github.com/Luiz-Garcia-R/autoproteoR\n"
  )
}
