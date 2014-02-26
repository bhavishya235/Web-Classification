### R code from vignette source 'Rcpp-extending.Rnw'

###################################################
### code chunk number 1: Rcpp-extending.Rnw:36-38
###################################################
prettyVersion <- packageDescription("Rcpp")$Version
prettyDate <- format(Sys.Date(), "%B %e, %Y")


###################################################
### code chunk number 2: Rcpp-extending.Rnw:45-67
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

require(inline)
require(highlight)
require(Rcpp)


###################################################
### code chunk number 4: Rcpp-extending.Rnw:105-121
###################################################
code <- '
// we get a list from R
List input(input_) ;

// pull std::vector<double> from R list
// this is achieved through an implicit call to Rcpp::as
std::vector<double> x = input["x"] ;

// return an R list
// this is achieved through implicit call to Rcpp::wrap
return List::create(
    _["front"] = x.front(),
    _["back"]  = x.back()
    ) ;
'
writeLines( code, "code.cpp" )


###################################################
### code chunk number 5: Rcpp-extending.Rnw:123-124
###################################################
external_highlight( "code.cpp", type = "LATEX", doc = FALSE )


###################################################
### code chunk number 6: Rcpp-extending.Rnw:127-133
###################################################
fx <- cxxfunction( signature( input_ = "list"),
	paste( readLines( "code.cpp" ), collapse = "\n" ),
	plugin = "Rcpp"
	)
input <- list( x = seq(1, 10, by = 0.5) )
fx( input )


###################################################
### code chunk number 14: Rcpp-extending.Rnw:339-340
###################################################
unlink( "code.cpp" )


