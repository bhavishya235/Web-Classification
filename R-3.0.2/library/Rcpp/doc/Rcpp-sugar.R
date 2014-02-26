### R code from vignette source 'Rcpp-sugar.Rnw'

###################################################
### code chunk number 1: Rcpp-sugar.Rnw:32-34
###################################################
prettyVersion <- packageDescription("Rcpp")$Version
prettyDate <- format(Sys.Date(), "%B %e, %Y")


###################################################
### code chunk number 2: Rcpp-sugar.Rnw:41-59
###################################################
link <- function( f, package, text = f, root = "http://finzi.psych.upenn.edu/R/library/" ){
	h <- if( missing(package) ) {
		as.character( help( f ) )
	} else {
		as.character( help( f, package = paste( package, sep = "" ) ) )
	}
	if( ! length(h) ){
		sprintf( "\\\\textbf{%s}", f )
	} else {
		rx <- "^.*/([^/]*?)/help/(.*?)$"
		package <- sub( rx, "\\1", h, perl = TRUE )
		page <- sub( rx, "\\2", h, perl = TRUE )
		sprintf( "\\\\href{%s%s/html/%s.html}{\\\\texttt{%s}}", root, package, page, text )
	}
}
linkS4class <- function( cl, package, text = cl, root = "http://finzi.psych.upenn.edu/R/library/" ){
	link( sprintf("%s-class", cl), package, text, root )
}


###################################################
### code chunk number 4: Rcpp-sugar.Rnw:116-119 (eval = FALSE)
###################################################
## foo <- function(x, y){
## 	ifelse( x < y, x*x, -(y*y) )
## }


