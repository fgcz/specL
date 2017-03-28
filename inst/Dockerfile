FROM bioconductor/release_core2
MAINTAINER Christian Panse <cp@fgcz.ethz.ch>
LABEL version="1.0"
RUN R -e "install.packages(pkgs=c('devtools', 'knitr'), quiet=FALSE, repos='https://stat.ethz.ch/CRAN/', type='source'); \
   library(devtools); \
   source('https://bioconductor.org/biocLite.R'); \
   biocLite('specL'); \ 
   install_git('https://github.com/protViz/prozor', build_vignettes = TRUE, quiet = FALSE)" 
