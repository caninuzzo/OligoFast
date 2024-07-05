<a name="readme-top"></a>


<!-- PROJECT LOGO -->


# OligoFast R package

<img src="img/logo_oligofast_resized_400px.png" alt="Logo" width="200" height="200" align="left">

Your R companion to design primers couples for various projects!


The main goal of `OligoFast` package is to provide an easy way to design primers couples in R for various projects.  
It also enables to test the primers couples on datasets and to conduct *in silico* investigations before to move into the laboratory.


## Installation

To install the package `OligoFast` you need `devtools` package on R, then use:

``` r
devtools::install_github("https://github.com/caninuzzo/OligoFast")
```

If the installation failed with the following error: "HTTP error 403. API rate limit exceeded (...)", then follow the
instruction on the console:
``` r
usethis::create_github_token()
# for that you need a github account
usethis::edit_r_environ()
# in the file .Renviron (hidden) add the following line: GITHUB_PAT=token (where token is the token given in github using the previous line)
```
After running that, save .Renviron, restart R, and use devtools::install_github() function again.


## Tutorial

The complete package description and tutorial is available <a href="https://caninuzzo.github.io/OligoFast/"><strong>here!</strong></a>


## Citation

If you use `OligoFast` please cite CANINO, A. (2024). R package OligoFast: presentation and pipeline tutorial (v.0.1.0). Zenodo. [doi:10.5281/zenodo.12664517](https://caninuzzo.github.io/OligoFast/)

----

## Current version: v.0.1.0

**Modifications:**  
*(v.0.0.1->v.0.1.0)*  
After running the pipeline several times with different datasets and confronting it to a large dataset, some issues were revealed and improvements were thus made on the following functions:  

- **OligoFindR()**: the processing time of this function was too long with large datasets. The algorithm has been improved to be more time-efficient.   
- **OligoCheckR()**: redesign of the function. In the previous version, the oligonucleotides with degenerated bases appeared to be poorly considered when counting the occurrences. The argument 'nOcc' has been replaced by 'pOcc' which is now a proportion. The user can choose the maximum number of degenerated base allowed with the new argument 'maxDeg'.ces. The current version corrects it and allows any tolerance (degenerated bases) around the resulting oligonucleotides. A new 'score' has replaced the field 'warnings' in the output dataframe and the graphical output of the function has been improved. More details in the tutorial!  
- **OligoTestR()**: a small issue when testing a given couple of primers has been corrected.

**Previous versions:**

- *v.0.1.0* (current) 05.07.2024
- *v.0.0.1* first release, 22.04.2024


<!-- DOI -->
<!-- old one [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11001268.svg)](https://doi.org/10.5281/zenodo.11001268)-->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12664517.svg)](https://doi.org/10.5281/zenodo.12664517)



