### R code from vignette source 'Rcpp-FAQ.Rnw'

###################################################
### code chunk number 1: Rcpp-FAQ.Rnw:37-39
###################################################
prettyVersion <- packageDescription("Rcpp")$Version
prettyDate <- format(Sys.Date(), "%B %e, %Y")


###################################################
### code chunk number 2: Rcpp-FAQ.Rnw:46-66
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


###################################################
### code chunk number 3: Rcpp-FAQ.Rnw:85-86 (eval = FALSE)
###################################################
## vignette( "Rcpp-introduction" )


###################################################
### code chunk number 4: Rcpp-FAQ.Rnw:176-177 (eval = FALSE)
###################################################
## vignette( "Rcpp-package" )


###################################################
### code chunk number 5: Rcpp-FAQ.Rnw:188-194
###################################################
fx <- cxxfunction(signature( x = "numeric" ),
    'NumericVector xx(x);
     return wrap( std::accumulate( xx.begin(), xx.end(), 0.0));',
    plugin = "Rcpp")
res <- fx( seq( 1, 10, by = 0.5 ) )
res


###################################################
### code chunk number 6: Rcpp-FAQ.Rnw:195-196
###################################################
stopifnot( identical( res, sum( seq( 1, 10, by = 0.5 ) ) ) )


###################################################
### code chunk number 7: Rcpp-FAQ.Rnw:202-203 (eval = FALSE)
###################################################
## list.files( system.file( "unitTests", package = "Rcpp" ), pattern = "^runit[.]" )


###################################################
### code chunk number 8: Rcpp-FAQ.Rnw:210-213 (eval = FALSE)
###################################################
## fx <- cxxfunction( signature(),
## 	paste( readLines( "myfile.cpp"), collapse = "\n" ),
## 	plugin = "Rcpp" )


###################################################
### code chunk number 13: Rcpp-FAQ.Rnw:409-428
###################################################
inc <- 'template <typename T>
        class square : public std::unary_function<T,T> {
        public:
            T operator()( T t) const { return t*t ;}
        };
       '

src <- '
       double x = Rcpp::as<double>(xs);
       int i = Rcpp::as<int>(is);
       square<double> sqdbl;
       square<int> sqint;
       return Rcpp::DataFrame::create(Rcpp::Named("x", sqdbl(x)),
                                      Rcpp::Named("i", sqint(i)));
       '
fun <- cxxfunction(signature(xs="numeric", is="integer"),
                   body=src, include=inc, plugin="Rcpp")

fun(2.2, 3L)


###################################################
### code chunk number 14: Rcpp-FAQ.Rnw:462-476
###################################################
writeLines( '
    // copy the data to armadillo structures
    arma::colvec x = Rcpp::as<arma::colvec> (x_);
    arma::mat Y = Rcpp::as<arma::mat>( Y_ ) ;
    arma::colvec z = Rcpp::as<arma::colvec>( z_ ) ;

    // calculate the result
    double result = arma::as_scalar(
        arma::trans(x) * arma::inv(Y) * z
        );

    // return it to R
    return Rcpp::wrap( result );
', "myfile.cpp" )


###################################################
### code chunk number 15: Rcpp-FAQ.Rnw:478-479
###################################################
external_highlight( "myfile.cpp", type = "LATEX", doc = FALSE )


###################################################
### code chunk number 16: Rcpp-FAQ.Rnw:484-489 (eval = FALSE)
###################################################
## fx <- cxxfunction(
## 	signature(x_ = "numeric", Y_ = "matrix", z_ = "numeric" ),
## 	paste( readLines( "myfile.cpp" ), collapse = "\n" ),
## 	plugin = "RcppArmadillo" )
## fx( 1:4, diag( 4 ), 1:4 )


###################################################
### code chunk number 17: Rcpp-FAQ.Rnw:491-492
###################################################
unlink( "myfile.cpp" )


###################################################
### code chunk number 18: Rcpp-FAQ.Rnw:520-526
###################################################
fx <- cxxfunction(signature(),
                  'RNGScope();
                   return rnorm(5, 0, 100);',
                  plugin="Rcpp")
fx()
fx()


###################################################
### code chunk number 19: Rcpp-FAQ.Rnw:537-545
###################################################
src <- 'Rcpp::NumericVector v(4);
        v[0] = R_NegInf;  // -Inf
        v[1] = NA_REAL;   // NA
        v[2] = R_PosInf;  // Inf
        v[3] = 42;        // see the Hitchhiker Guide
        return Rcpp::wrap(v);'
fun <- cxxfunction(signature(), src, plugin="Rcpp")
fun()


###################################################
### code chunk number 20: Rcpp-FAQ.Rnw:558-566 (eval = FALSE)
###################################################
## txt <- 'arma::mat Am = Rcpp::as< arma::mat >(A);
##         arma::mat Bm = Rcpp::as< arma::mat >(B);
##         return Rcpp::wrap( Am * Bm );'
## mmult <- cxxfunction(signature(A="numeric", B="numeric"),
##                      body=txt, plugin="RcppArmadillo")
## A <- matrix(1:9, 3, 3)
## B <- matrix(9:1, 3, 3)
## C <- mmult(A, B)


###################################################
### code chunk number 21: Rcpp-FAQ.Rnw:589-607 (eval = FALSE)
###################################################
## ## simple example of seeding RNG and drawing one random number
## gslrng <- '
## int seed = Rcpp::as<int>(par) ;
## gsl_rng_env_setup();
## gsl_rng *r = gsl_rng_alloc (gsl_rng_default);
## gsl_rng_set (r, (unsigned long) seed);
## double v = gsl_rng_get (r);
## gsl_rng_free(r);
## return Rcpp::wrap(v);'
## 
## plug <- Rcpp:::Rcpp.plugin.maker(
##     include.before = "#include <gsl/gsl_rng.h>",
##     libs = paste("-L/usr/local/lib/R/site-library/Rcpp/lib -lRcpp",
##                  "-Wl,-rpath,/usr/local/lib/R/site-library/Rcpp/lib",
##                  "-L/usr/lib -lgsl -lgslcblas -lm"))
## registerPlugin("gslDemo", plug )
## fun <- cxxfunction(signature(par="numeric"), gslrng, plugin="gslDemo")
## fun(0)


###################################################
### code chunk number 22: Rcpp-FAQ.Rnw:627-634 (eval = FALSE)
###################################################
## myplugin <- getPlugin("Rcpp")
## myplugin$env$PKG_CXXFLAGS <- "-std=c++0x"
## f <- cxxfunction(signature(), settings=myplugin, body='
## +    std::vector<double> x = { 1.0, 2.0, 3.0 };  // fails without -std=c++0x
## +    return Rcpp::wrap(x);
## + ')
## f()


###################################################
### code chunk number 23: Rcpp-FAQ.Rnw:647-659 (eval = FALSE)
###################################################
## src <- '
##   Rcpp::NumericMatrix x(2,2);
##   x.fill(42);                           // or more interesting values
##   Rcpp::List dimnms =                   // two vec. with static names
##       Rcpp::List::create(Rcpp::CharacterVector::create("cc", "dd"),
##                          Rcpp::CharacterVector::create("ee", "ff"));
##   // and assign it
##   x.attr("dimnames") = dimnms;
##   return(x);
## '
## fun <- cxxfunction(signature(), body=src, plugin="Rcpp")
## fun()


###################################################
### code chunk number 24: Rcpp-FAQ.Rnw:671-684 (eval = FALSE)
###################################################
## BigInts <- cxxfunction(signature(),
##   'std::vector<long> bigints;
##    bigints.push_back(12345678901234567LL);
##    bigints.push_back(12345678901234568LL);
##    Rprintf("Difference of %ld\\n", 12345678901234568LL - 12345678901234567LL);
##    return wrap(bigints);', plugin="Rcpp", includes="#include <vector>")
## 
## retval<-BigInts()
## 
## # Unique 64-bit integers were cast to identical lower precision numerics
## # behind my back with no warnings or errors whatsoever.  Error.
## 
## stopifnot(length(unique(retval)) == 2)


