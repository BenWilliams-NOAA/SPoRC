library(renv) # load in library to test pkg installation in new env
renv::init(bare = TRUE)  # create isolated lib
Sys.unsetenv("GITHUB_PAT") # dont use authenticaiton token
devtools::install_github("chengmatt/SPoRC") # try to install SPoRC from git
library(SPoRC) # load in SPoRC

# clean up renv
renv::deactivate()
unlink("renv", recursive = TRUE)
unlink("renv.lock")
unlink(".Rprofile")
