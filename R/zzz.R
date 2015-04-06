# ##
# Anything needed for the functioning of the namespace should 
# be handled at load/unload times by the .onLoad and .onUnload 
# hooks. For example, DLLs can be loaded (unless done by a 
#                                         useDynLib directive in 
#                                         the NAMESPACE file) 
# and initialized in .onLoad and unloaded in .onUnload. 
# Use .onAttach only for actions that are needed only when the 
# package becomes visible to the user (for example a start-up message) 
# or need to be run after the package environment has been created.

#.onUnload <- function(libpath){
#  library.dynam.unload("deform", libpath)
#}

##
# executed after .onLoad is executed, once the namespace is visible to user
.onAttach <- function(...) {
  pack <- 'deform'
  ##
  packageStartupMessage(
    sprintf("Loaded %s (%s) -- crustal deformation modeling tools", pack, utils::packageVersion(pack))
    )
  ## options which can be reset by the user
  options(deform.ops=list())
}
