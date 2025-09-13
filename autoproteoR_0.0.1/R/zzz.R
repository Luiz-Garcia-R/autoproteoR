.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "\n",
    crayon::green("autoproteoR "), "loaded successfully!\n",
    "--------------------------------------------------\n",
    "A package for streamlined proteomics data analysis.\n",
    "Use ", crayon::green("?autoproteoR"), " for general help.\n",
    "--------------------------------------------------\n"
  )
}
