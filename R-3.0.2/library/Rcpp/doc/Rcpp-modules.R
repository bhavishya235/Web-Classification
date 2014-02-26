### R code from vignette source 'Rcpp-modules.Rnw'

###################################################
### code chunk number 1: Rcpp-modules.Rnw:30-33
###################################################
suppressMessages(require(Rcpp))
prettyVersion <- packageDescription("Rcpp")$Version
prettyDate <- format(Sys.Date(), "%B %e, %Y")


###################################################
### code chunk number 7: Rcpp-modules.Rnw:204-206 (eval = FALSE)
###################################################
## f1 <- cxxfunction( , "", includes = unifModCode, plugin = "Rcpp" )
## getDynLib(f1)  ## will display information about 'f1' just created


###################################################
### code chunk number 8: Rcpp-modules.Rnw:212-231 (eval = FALSE)
###################################################
## setClass( "Uniform", representation( pointer = "externalptr" ) )
## 
## # helper
## Uniform_method <- function(name) {
##     paste( "Uniform", name, sep = "__" )
## }
## 
## # syntactic sugar to allow object$method( ... )
## setMethod( "$", "Uniform", function(x, name ) {
##     function(...) .Call( Uniform_method(name) , x@pointer, ... )
## } )
## # syntactic sugar to allow new( "Uniform", ... )
## setMethod( "initialize", "Uniform", function(.Object, ...) {
##     .Object@pointer <- .Call( Uniform_method("new"), ... )
##     .Object
## } )
## 
## u <- new( "Uniform", 0, 10 )
## u$draw( 10L )


###################################################
### code chunk number 10: Rcpp-modules.Rnw:282-296 (eval = FALSE)
###################################################
## inc <- '
## using namespace Rcpp;
## 
## double norm( double x, double y ) {
##     return sqrt( x*x + y*y );
## }
## 
## RCPP_MODULE(mod) {
##     function( "norm", &norm );
## }
## '
## 
## fx <- cxxfunction(signature(), plugin="Rcpp", include=inc)
## mod <- Module("mod", getDynLib(fx))


###################################################
### code chunk number 11: Rcpp-modules.Rnw:306-309 (eval = FALSE)
###################################################
## require(nameOfMyModulePackage)
## mod <- new( mod )
## mod$norm( 3, 4 )


###################################################
### code chunk number 14: Rcpp-modules.Rnw:361-370 (eval = FALSE)
###################################################
## require( Rcpp )
## 
## yd <- Module("yada", getDynLib(fx))
## yd$bar(2L)
## yd$foo(2L, 10.0)
## yd$hello()
## yd$bla()
## yd$bla1(2L)
## yd$bla2(2L, 5.0)


###################################################
### code chunk number 15: Rcpp-modules.Rnw:377-385 (eval = FALSE)
###################################################
## require(myModulePackage)    ## or whichever name was chose
## 
## bar(2L)
## foo(2L, 10.0)
## hello()
## bla()
## bla1(2L)
## bla2(2L, 5.0)


###################################################
### code chunk number 17: Rcpp-modules.Rnw:422-424 (eval = FALSE)
###################################################
## mod <- Module("mod", getDynLib(fx))
## show(mod$norm)


###################################################
### code chunk number 19: Rcpp-modules.Rnw:450-455 (eval = FALSE)
###################################################
## norm <- mod$norm
## norm()
## norm( y = 2 )
## norm( x = 2, y = 3 )
## args( norm )


###################################################
### code chunk number 21: Rcpp-modules.Rnw:476-478 (eval = FALSE)
###################################################
## norm <- mod$norm
## args(norm)


###################################################
### code chunk number 23: Rcpp-modules.Rnw:497-499 (eval = FALSE)
###################################################
## norm <- mod$norm
## args( norm )


###################################################
### code chunk number 25: Rcpp-modules.Rnw:547-556 (eval = FALSE)
###################################################
## ## assumes   fx_unif <- cxxfunction(...)   has ben run
## unif_module <- Module( "unif_module", getDynLib(fx_unif ) )
## Uniform <- unif_module$Uniform
## u <- new( Uniform, 0, 10 )
## u$draw( 10L )
## u$range()
## u$max <- 1
## u$range()
## u$draw( 10 )


###################################################
### code chunk number 35: Rcpp-modules.Rnw:721-727 (eval = FALSE)
###################################################
## Bar <- mod_bar$Bar
## b <- new( Bar, 10 )
## b$x + b$x
## b$stats()
## b$x <- 10
## b$stats()


###################################################
### code chunk number 38: Rcpp-modules.Rnw:804-808 (eval = FALSE)
###################################################
## setMethod( "show", yada$World , function(object) {
##     msg <- paste( "World object with message : ", object$greet() )
##     writeLines( msg )
## } )


###################################################
### code chunk number 40: Rcpp-modules.Rnw:871-884 (eval = FALSE)
###################################################
## ## for code compiled on the fly using cxxfunction() into 'fx_vec', we use
## mod_vec <- Module( "mod_vec", getDynLib(fx_vec), mustStart = TRUE )
## vec <- mod_vec$vec
## ## and that is not needed in a package setup as e.g. one created
## ## via Rcpp.package.skeleton(..., module=TRUE)
## v <- new( vec )
## v$reserve( 50L )
## v$assign( 1:10 )
## v$push_back( 10 )
## v$size()
## v$capacity()
## v[[ 0L ]]
## v$as.vector()


###################################################
### code chunk number 41: Rcpp-modules.Rnw:898-899
###################################################
options( prompt = " ", continue = " " )


###################################################
### code chunk number 42: Rcpp-modules.Rnw:902-903 (eval = FALSE)
###################################################
## import( Rcpp )


###################################################
### code chunk number 43: Rcpp-modules.Rnw:910-913 (eval = FALSE)
###################################################
## .onLoad <- function(libname, pkgname) {
##     loadRcppModules()
## }


###################################################
### code chunk number 44: Rcpp-modules.Rnw:939-942 (eval = FALSE)
###################################################
## loadModule("yada")
## loadModule("stdVector")
## loadModule("NumEx")


###################################################
### code chunk number 45: Rcpp-modules.Rnw:954-959 (eval = FALSE)
###################################################
## yada <- Module( "yada" )
## 
## .onLoad <- function(libname, pkgname) {
##     # placeholder
## }


###################################################
### code chunk number 46: Rcpp-modules.Rnw:962-963
###################################################
options( prompt = "> ", continue = "+ " )


###################################################
### code chunk number 47: Rcpp-modules.Rnw:973-974 (eval = FALSE)
###################################################
## Rcpp.package.skeleton( "testmod", module = TRUE )


###################################################
### code chunk number 48: Rcpp-modules.Rnw:987-989 (eval = FALSE)
###################################################
## yada <- Module( "yada" )
## prompt( yada, "yada-module.Rd" )


