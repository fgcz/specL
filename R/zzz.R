#R

.onAttach <- function(lib, pkg){
        if(interactive()){
                version <- packageVersion('specL')
                packageStartupMessage("Package 'specL' version ", version)
          invisible()
        }

}

