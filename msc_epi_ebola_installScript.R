#~ Modified from script by T Jombart (E. Volz e.volz@imperial.ac.uk )
#~ http://adegenet.r-forge.r-project.org/files/hackLib/hackLib.R
## HACK TO INSTALL R PACKAGES ON WINDOWS SYSTEMS
## WITHOUT ADMIN RIGHTS
##
## Thibaut Jombart, March 2013
## tjombart@imperial.ac.uk
##
##
myPath="C:/Users/Public/Rlibs"
if(file.exists(myPath)) myPath <- paste(myPath,"Rlibs", sep="/")
myPath <- gsub("/+", "/", myPath)
if(!file.exists(myPath)){
	if(!dir.create(myPath)) stop(paste(myPath, "could not be created."))
}

.libPaths(myPath)
cat(paste("\nR library set to:",myPath,"\n"))
if(as.numeric(file.info(myPath)$mode)<500) warning(paste(myPath, "may not be 'writable'."))
return(invisible())


install.packages('ape')
install.packages('devtools')
require(devtools)
install_github( 'mrc-ide/skygrowth')
install_github( 'emvolz/treedater')

