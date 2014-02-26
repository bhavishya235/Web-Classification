#  File src/library/base/R/Bessel.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

besselI <- function(x, nu, expon.scaled = FALSE)
{
    ## Oddly, this is passed as a double to fit Math3 semantics
    .Internal(besselI(x,nu, 1 + as.logical(expon.scaled)))
}
besselK <- function(x, nu, expon.scaled = FALSE)
{
    .Internal(besselK(x,nu, 1 + as.logical(expon.scaled)))
}
besselJ <- function(x, nu) .Internal(besselJ(x,nu))
besselY <- function(x, nu) .Internal(besselY(x,nu))
#  File src/library/base/R/Defunct.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

.Defunct <- function(new, package=NULL, msg) {
    if (missing(msg)) {
	msg <- gettextf("'%s' is defunct.\n",
			as.character(sys.call(sys.parent())[[1L]]))
	if(!missing(new))
	    msg <- c(msg, gettextf("Use '%s' instead.\n", new))
	msg <- c(msg,
		 if(!is.null(package))
		 gettextf("See help(\"Defunct\") and help(\"%s-defunct\").",
			  package)
		 else gettext("See help(\"Defunct\")"))
    }
    else msg <- as.character(msg)

    stop(paste(msg, collapse=""), call. = FALSE, domain = NA)
}

## Version <- function() .Defunct("R.Version")
## provide <- function(package) .Defunct()

## <entry>
## Deprecated in 1.2.0
## Defunct in 1.3.0
# getenv <- function(...) .Defunct("Sys.getenv")
## </entry>

## <entry>
## Deprecated in 1.2.3
## Defunct in 1.3.0
## Removed in 1.4.0: conflicts with lattice
## dotplot <- function(...) .Defunct()
## stripplot <- function(...) .Defunct()
## </entry>

## <entry>
## Deprecated in 1.3.0
## Defunct in 1.4.0
## read.table.url <- function(url, method, ...) .Defunct("read.table(url())")
## scan.url <- function(url, file = tempfile(), method, ...)
##     .Defunct("scan(url())")
## source.url <- function(url, file = tempfile(), method, ...)
##     .Defunct("source(url())")
## httpclient <- function(url, port=80, error.is.fatal=TRUE, check.MIME.type=TRUE,
##                        file=tempfile(), drop.ctrl.z=TRUE)
##     .Defunct()
## parse.dcf <- function(text = NULL, file = "", fields = NULL,
##                       versionfix = FALSE) .Defunct("read.dcf")
## </entry>

## <entry>
## Deprecated in 1.4.0
## Defunct in 1.5.0
# .Alias <- function(expr) .Defunct()
## </entry>

## <entry>
## Deprecated in 1.6.0
## Defunct in 1.7.0
## machine <- function() .Defunct()
## Machine <- function() .Defunct(".Machine")
## Platform <- function() .Defunct(".Platform")
## restart <- function() .Defunct("try")
## </entry>

## <entry>
## Deprecated in 1.7.0
## Defunct in 1.8.0
## printNoClass <- function(x, digits = NULL, quote = TRUE, na.print = NULL,
##                          print.gap = NULL, right = FALSE, ...)
##     .Defunct()
## </entry>

## <entry>
## Deprecated in 1.8.0
## Defunct in 1.9.0
## codes <- function(x, ...) .Defunct()
## codes.factor <- function(x, ...) .Defunct("unclass")
## codes.ordered <- function(x, ...) .Defunct("unclass")
## `codes<-` <- function(x, ..., value) .Defunct()
# removed in 2.9.1, as it caused confusion for an S4 class union of this name.
#print.atomic <- function(x, quote = TRUE, ...) .Defunct("print.default")
## </entry>

## <entry>
## Deprecated in 1.9.0
## Defunct in 2.0.0
## La.eigen <- function(x, symmetric, only.values = FALSE,
##                      method = c("dsyevr", "dsyev")) .Defunct("eigen")
## tetragamma <- function(x) .Defunct("psigamma")
## pentagamma <- function(x) .Defunct("psigamma")
## package.description <- function(pkg, lib.loc = NULL, fields = NULL)
##     .Defunct("packageDescription")
## </entry>

## <entry>
## Deprecated in 2.1.0
## Defunct in 2.2.0
## delay <- function(x, env=.GlobalEnv) .Defunct("delayedAssign")
## loadURL <- function (url, envir = parent.frame(), quiet = TRUE, ...)
##     .Defunct("load(url())")
## </entry>

## Defunct in 2.3.0
## write.table0 <-
## function (x, file = "", append = FALSE, quote = TRUE, sep = " ",
##           eol = "\n", na = "NA", dec = ".", row.names = TRUE,
##           col.names = TRUE, qmethod = c("escape", "double"))
##     .Defunct("write.table")
## format.char <- function(x, width = NULL, flag = "-")
##     .Defunct("format.default")
## </entry>

## <entry>
## Deprecated in 2.3.0
## Defunct in 2.4.0
# La.chol <- function(x) .Defunct("chol")
# La.chol2inv <- function(x, size = ncol(x)) .Defunct("chol2inv")
## </entry>

## <entry>
## Deprecated in 2.4.0
## Defunct in 2.5.0
## symbol.C <- function(name)
## {
##     warning("'symbol.C' is not needed: please remove it", immediate.=TRUE)
##     name
## }
## symbol.For <- function(name)
## {
##     warning("'symbol.For' is not needed: please remove it", immediate.=TRUE)
##     name
## }
## </entry>

## <entry>
## Deprecated in 1999
## Defunct in 2.5.0
# unix <- function(call, intern = FALSE) .Defunct("system")
## </entry>

## <entry>
## Deprecated in 2.7.0
## Defunct in 2.8.0
## gammaCody <- function(x) .Defunct("gamma")
## </entry>

## <entry>
## Deprecated inter alia in 2.8.1
## Defunct in 2.9.0
## manglePackageName <- function (pkgName, pkgVersion) .Defunct()
## </entry>

## <entry>
## Deprecated in 2.12.2 (and only ever experimental)
## Defunct in 2.13.0
## .Import <- function(...)
##     .Defunct(msg = "namespaces should be specified via the 'NAMESPACE' file")
## .ImportFrom <- function(name, ...)
##     .Defunct(msg = "namespaces should be specified via the 'NAMESPACE' file")
## .Export <- function(...)
##     .Defunct(msg = "namespaces should be specified via the 'NAMESPACE' file")
## .S3method <- function(generic, class, method)
##     .Defunct(msg = "namespaces should be specified via the 'NAMESPACE' file")
## </entry>

## <entry>
## Deprecated in 2.14.0
## Defunct in 2.15.0
mem.limits <- function(nsize=NA, vsize=NA) .Defunct("gc")
## </entry>

## <entry>
## Deprecated in 2.13.1
## Defunct in 2.15.0
.readRDS <- function(...) .Defunct("readRDS")
.saveRDS <- function(...) .Defunct("saveRDS")
## </entry>

## <entry>
## Deprecated in 2.5.0
## Removed in 2.15.0
# Sys.putenv <- function(...) .Defunct("Sys.setenv")
## </entry>
#  File src/library/base/R/Deprecated.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2013 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

###----- NOTE:	../man/base-deprecated.Rd   must be synchronized with this file!
###		-------------------------
.Deprecated <- function(new, package = NULL, msg,
			old = as.character(sys.call(sys.parent()))[1L])
{
    msg <- if( missing(msg) ) {
	msg <- gettextf("'%s' is deprecated.\n", old)
	if(!missing(new))
	    msg <- c(msg, gettextf("Use '%s' instead.\n", new))
	c(msg,
	  if(!is.null(package))
	  gettextf("See help(\"Deprecated\") and help(\"%s-deprecated\").",
		   package)
	  else gettext("See help(\"Deprecated\")"))
    }
    else as.character(msg)
    warning(paste(msg, collapse=""), call. = FALSE, domain = NA)
}

## consider keeping one (commented) entry here, for easier additions

## <entry>
## Deprecated in 3.0.0
.find.package <- function(...)
{
    .Deprecated("find.package")
    find.package(...)
}

.path.package <- function(...)
{
    .Deprecated("path.package")
    path.package(...)
}
## </entry>
#  File src/library/base/R/LAPACK.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

La.svd <- function(x, nu = min(n, p), nv = min(n, p))
{
    if(!is.logical(x) && !is.numeric(x) && !is.complex(x))
	stop("argument to 'La.svd' must be numeric or complex")
    if (any(!is.finite(x))) stop("infinite or missing values in 'x'")
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    if(!n || !p) stop("a dimension is zero")

    if(is.complex(x)) {
        if(nu == 0L) {
            jobu <- "N"
            u <- matrix(0+0i, 1L, 1L)  # dim is checked
        }
        else if(nu == n) {
            jobu <- ifelse(n > p, "A", "S")
            u <- matrix(0+0i, n, n)
        }
        else if(nu == p) {
            jobu <- ifelse(n > p, "S", "A")
            u <- matrix(0+0i, n, p)
        }
        else
            stop("'nu' must be 0, nrow(x) or ncol(x)")

        if (nv == 0L) {
            jobv <- "N"
            v <- matrix(0+0i, 1L, 1L) # dim is checked
        }
        else if (nv == n) {
            jobv <- ifelse(n > p, "A", "S")
            v <- matrix(0+0i, min(n, p), p)
        }
        else if (nv == p) {
            jobv <- ifelse(n > p, "S", "A")
            v <- matrix(0+0i, p, p)
        }
        else
            stop("'nv' must be 0, nrow(x) or ncol(x)")
        res <- .Internal(La_svd_cmplx(jobu, jobv, x, double(min(n, p)), u, v))
        return(res[c("d", if(nu) "u", if(nv) "vt")])
    } else {
        if(nu || nv) {
            np <- min(n, p)
            if(nu <= np && nv <= np) {
                jobu <- "S"
                u <- matrix(0, n, np)
                vt <- matrix(0, np, p)
                nu0 <- nv0 <- np
            } else {
                jobu <- "A"
                u <- matrix(0, n, n)
                vt <- matrix(0, p, p)
                nu0 <- n; nv0 <- p
            }
        } else {
            jobu <- "N"
            # these dimensions _are_ checked, but unused
            u <- matrix(0, 1L, 1L)
            vt <- matrix(0, 1L, 1L)
        }
        ## type is now coerced internally
        res <- .Internal(La_svd(jobu, x, double(min(n,p)), u, vt))
        res <- res[c("d", if(nu) "u", if(nv) "vt")]
        if(nu && nu < nu0) res$u <- res$u[, 1L:min(n, nu), drop = FALSE]
        if(nv && nv < nv0) res$vt <- res$vt[1L:min(p, nv), , drop = FALSE]
        return(res)
    }
    ## not reached
}
#  File src/library/base/R/New-Internal.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

geterrmessage <- function() .Internal(geterrmessage())

try <- function(expr, silent = FALSE) {
    tryCatch(expr, error = function(e) {
        call <- conditionCall(e)
        if (! is.null(call)) {
            ## Patch up the call to produce nicer result for testing as
            ## try(stop(...)).  This will need adjusting if the
            ## implementation of tryCatch changes.
            ## Use identical() since call[[1L]] can be non-atomic.
            if (identical(call[[1L]], quote(doTryCatch)))
                call <- sys.call(-4L)
            dcall <- deparse(call)[1L]
            prefix <- paste("Error in", dcall, ": ")
            LONG <- 75L # to match value in errors.c
            msg <- conditionMessage(e)
            sm <- strsplit(msg, "\n")[[1L]]
            w <- 14L + nchar(dcall, type="w") + nchar(sm[1L], type="w")
            ## this could be NA if any of this is invalid in a MBCS
            if(is.na(w))
                w <-  14L + nchar(dcall, type="b") + nchar(sm[1L], type="b")
            if (w > LONG)
                prefix <- paste0(prefix, "\n  ")
        }
        else prefix <- "Error : "
        msg <- paste0(prefix, conditionMessage(e), "\n")
        ## Store the error message for legacy uses of try() with
        ## geterrmessage().
        .Internal(seterrmessage(msg[1L]))
        if (! silent && identical(getOption("show.error.messages"), TRUE)) {
            cat(msg, file = stderr())
            .Internal(printDeferredWarnings())
        }
        invisible(structure(msg, class = "try-error", condition = e))
    })
}

comment <- function(x) .Internal(comment(x))
`comment<-` <- function(x, value) .Internal("comment<-"(x, value))

logb <- function(x, base=exp(1)) if(missing(base)) log(x) else log(x, base)

atan2 <- function(y, x) .Internal(atan2(y, x))

beta <- function(a, b) .Internal( beta(a, b))
lbeta <- function(a, b) .Internal(lbeta(a, b))

psigamma <- function(x, deriv = 0L) .Internal(psigamma(x, deriv))

factorial <- function(x) gamma(x + 1)
lfactorial <- function(x) lgamma(x + 1)

choose <- function(n, k) .Internal(choose(n, k))
lchoose <- function(n, k) .Internal(lchoose(n, k))

##-- 2nd part --
R.Version <- function() .Internal(Version())

commandArgs <- function(trailingOnly = FALSE) {
    args <- .Internal(commandArgs())
    if(trailingOnly) {
        m <- match("--args", args, 0L)
        if(m) args[-seq_len(m)] else character()
    } else args
}

args <- function(name) .Internal(args(name))

cbind <- function(..., deparse.level = 1)
    .Internal(cbind(deparse.level, ...))

rbind <- function(..., deparse.level = 1)
    .Internal(rbind(deparse.level, ...))

## for methods:::bind_activation
.__H__.cbind <- cbind
.__H__.rbind <- rbind


# convert deparsing options to bitmapped integer

.deparseOpts <- function(control) {
    opts <- pmatch(as.character(control),
                   ## the exact order of these is determined by the integer codes in
                   ## ../../../include/Defn.h
                   c("all",
                     "keepInteger", "quoteExpressions", "showAttributes",
                     "useSource", "warnIncomplete", "delayPromises",
                     "keepNA", "S_compatible"))
    if (any(is.na(opts)))
        stop(sprintf(ngettext(as.integer(sum(is.na(opts))),
                              "deparse option %s is not recognized",
                              "deparse options %s are not recognized"),
                     paste(sQuote(control[is.na(opts)]), collapse=", ")),
             call. = FALSE, domain = NA)
    if (any(opts == 1L))
        opts <- unique(c(opts[opts != 1L], 2L,3L,4L,5L,6L,8L)) # not (7,9)
    return(sum(2^(opts-2)))
}

deparse <-
    function(expr, width.cutoff = 60L,
	     backtick = mode(expr) %in% c("call", "expression", "(", "function"),
	     control = c("keepInteger", "showAttributes", "keepNA"),
             nlines = -1L)
    .Internal(deparse(expr, width.cutoff, backtick,
                      .deparseOpts(control), nlines))

do.call <- function(what, args, quote = FALSE, envir = parent.frame())
{
    if (!is.list(args))
	stop("second argument must be a list")
    if (quote)
	args <- lapply(args, enquote)
    .Internal(do.call(what, args, envir))
}

drop <- function(x) .Internal(drop(x))

format.info <- function(x, digits = NULL, nsmall = 0L)
    .Internal(format.info(x, digits, nsmall))

gc <- function(verbose = getOption("verbose"),	reset=FALSE)
{
    res <- .Internal(gc(verbose, reset))
    res <- matrix(res, 2L, 7L,
		  dimnames = list(c("Ncells","Vcells"),
		  c("used", "(Mb)", "gc trigger", "(Mb)",
		    "limit (Mb)", "max used", "(Mb)")))
    if(all(is.na(res[, 5L]))) res[, -5L] else res
}
gcinfo <- function(verbose) .Internal(gcinfo(verbose))
gctorture <- function(on = TRUE) .Internal(gctorture(on))
gctorture2 <- function(step, wait = step, inhibit_release = FALSE)
    .Internal(gctorture2(step, wait, inhibit_release))

is.unsorted <- function(x, na.rm = FALSE, strictly = FALSE)
{
    if(length(x) <= 1L) return(FALSE)
    if(!na.rm && any(is.na(x)))## "FIXME" is.na(<large>) is "too slow"
	return(NA)
    ## else
    if(na.rm && any(ii <- is.na(x)))
	x <- x[!ii]
    .Internal(is.unsorted(x, strictly))
}

nchar <- function(x, type = "chars", allowNA = FALSE)
    .Internal(nchar(x, type, allowNA))

polyroot <- function(z) .Internal(polyroot(z))

readline <- function(prompt = "") .Internal(readline(prompt))
search <- function() .Internal(search())
searchpaths <- function()
{
    s <- search()
    paths <-
        lapply(seq_along(s), function(i) attr(as.environment(i), "path"))
    paths[[length(s)]] <- system.file()
    m <- grep("^package:", s)
    if(length(m)) paths[-m] <- as.list(s[-m])
    unlist(paths)
}

sprintf <- function(fmt, ...) .Internal(sprintf(fmt, ...))

##-- DANGER ! ---   substitute(list(...))  inside functions !!!
##substitute <- function(expr, env=baseenv()) .Internal(substitute(expr, env))

t.default <- function(x) .Internal(t.default(x))
typeof <- function(x) .Internal(typeof(x))


memory.profile <- function() .Internal(memory.profile())

capabilities <- function(what = NULL)
{
    z  <- .Internal(capabilities())
    if(!is.null(what))
        z <- z[match(what, names(z), 0L)]
    if(.Platform$OS.type == "windows") return(z)
    ## Now we need to deal with any NA entries if X11 is unknown.
    nas <- names(z[is.na(z)])
    if(any(nas %in% c("X11", "jpeg", "png", "tiff"))) {
        ## This might throw an X11 error
         z[nas] <- tryCatch(.Internal(capabilitiesX11()),
                            error = function(e) FALSE)
    }
    z
}

inherits <- function(x, what, which = FALSE)
	.Internal(inherits(x, what, which))

NextMethod <- function(generic=NULL, object=NULL, ...)
    .Internal(NextMethod(generic, object,...))

data.class <- function(x) {
    if (length(cl <- oldClass(x)))
	cl[1L]
    else {
	l <- length(dim(x))
        if (l == 2L) "matrix" else if(l) "array" else mode(x)
    }
}

encodeString <- function(x, width = 0L, quote = "", na.encode = TRUE,
                         justify = c("left", "right", "centre", "none"))
{
    at <- attributes(x)
    x <- as.character(x) # we want e.g. NULL to work
    attributes(x) <- at  # preserve names, dim etc
    oldClass(x) <- NULL  # but not class
    justify <- match(match.arg(justify),
                     c("left", "right", "centre", "none")) - 1L
    .Internal(encodeString(x, width, quote, justify, na.encode))
}

l10n_info <- function() .Internal(l10n_info())

iconv <- function(x, from = "", to = "", sub = NA, mark = TRUE, toRaw = FALSE)
{
    if(! (is.character(x) || (is.list(x) && is.null(oldClass(x)))))
        x <- as.character(x)
    .Internal(iconv(x, from, to, as.character(sub), mark, toRaw))
}

iconvlist <- function()
{
    int <- .Internal(iconv(NULL, "", "", "", TRUE, FALSE))
    if(length(int)) return(sort.int(int))
    icfile <- system.file("iconvlist", package="utils")
    if(!nchar(icfile, type="bytes"))
        stop("'iconvlist' is not available on this system")
    ext <- readLines(icfile)
    if(!length(ext)) stop("'iconvlist' is not available on this system")
    ## glibc has lines ending //, some versions with a header and some without.
    ## libiconv has lines with multiple entries separated by spaces
    cnt <- grep("//$", ext)
    if(length(cnt)/length(ext) > 0.5) {
        ext <- grep("//$", ext, value = TRUE)
        ext <- sub("//$", "", ext)
    }
    sort.int(unlist(strsplit(ext, "[[:space:]]")))
}

Cstack_info <- function() .Internal(Cstack_info())

reg.finalizer <- function(e, f, onexit = FALSE)
    .Internal(reg.finalizer(e, f, onexit))

Encoding <- function(x) .Internal(Encoding(x))
`Encoding<-` <- function(x, value) .Internal(setEncoding(x, value))

setTimeLimit <- function(cpu = Inf, elapsed = Inf, transient = FALSE)
    .Internal(setTimeLimit(cpu, elapsed, transient))
setSessionTimeLimit <- function(cpu = Inf, elapsed = Inf)
    .Internal(setSessionTimeLimit(cpu, elapsed))

icuSetCollate <- function(...) .Internal(icuSetCollate(...))

## base has no S4 generics
.noGenerics <- TRUE
#  File src/library/base/R/RNG.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

## Random Number Generator

## The available kinds are in
## ../../../include/Random.h  and ../../../main/RNG.c [RNG_Table]
##
RNGkind <- function(kind = NULL, normal.kind = NULL)
{
    kinds <- c("Wichmann-Hill", "Marsaglia-Multicarry", "Super-Duper",
               "Mersenne-Twister", "Knuth-TAOCP", "user-supplied",
               "Knuth-TAOCP-2002", "L'Ecuyer-CMRG", "default")
    n.kinds <- c("Buggy Kinderman-Ramage", "Ahrens-Dieter", "Box-Muller",
                 "user-supplied", "Inversion", "Kinderman-Ramage",
		 "default")
    do.set <- length(kind) > 0L
    if(do.set) {
	if(!is.character(kind) || length(kind) > 1L)
	    stop("'kind' must be a character string of length 1 (RNG to be used).")
	if(is.na(i.knd <- pmatch(kind, kinds) - 1L))
	    stop(gettextf("'%s' is not a valid abbreviation of an RNG", kind),
                 domain = NA)
        if(i.knd == length(kinds) - 1L) i.knd <- -1L
    } else i.knd <- NULL

    if(!is.null(normal.kind)) {
	if(!is.character(normal.kind) || length(normal.kind) != 1L)
	    stop("'normal.kind' must be a character string of length 1")
        normal.kind <- pmatch(normal.kind, n.kinds) - 1L
        if(is.na(normal.kind))
	    stop(gettextf("'%s' is not a valid choice", normal.kind),
                 domain = NA)
	if (normal.kind == 0L)
            warning("buggy version of Kinderman-Ramage generator used",
                    domain = NA)
         if(normal.kind == length(n.kinds) - 1L) normal.kind <- -1L
    }
    r <- 1L + .Internal(RNGkind(i.knd, normal.kind))
    r <- c(kinds[r[1L]], n.kinds[r[2L]])
    if(do.set || !is.null(normal.kind)) invisible(r) else r
}

set.seed <- function(seed, kind = NULL, normal.kind = NULL)
{
    kinds <- c("Wichmann-Hill", "Marsaglia-Multicarry", "Super-Duper",
               "Mersenne-Twister", "Knuth-TAOCP", "user-supplied",
               "Knuth-TAOCP-2002", "L'Ecuyer-CMRG", "default")
    n.kinds <- c("Buggy Kinderman-Ramage", "Ahrens-Dieter", "Box-Muller",
                 "user-supplied", "Inversion", "Kinderman-Ramage",
		 "default")
    if(length(kind) ) {
	if(!is.character(kind) || length(kind) > 1L)
	    stop("'kind' must be a character string of length 1 (RNG to be used).")
	if(is.na(i.knd <- pmatch(kind, kinds) - 1L))
	    stop(gettextf("'%s' is not a valid abbreviation of an RNG", kind),
                 domain = NA)
        if(i.knd == length(kinds) - 1L) i.knd <- -1L
    } else i.knd <- NULL

    if(!is.null(normal.kind)) {
	if(!is.character(normal.kind) || length(normal.kind) != 1L)
	    stop("'normal.kind' must be a character string of length 1")
        normal.kind <- pmatch(normal.kind, n.kinds) - 1L
        if(is.na(normal.kind))
	    stop(gettextf("'%s' is not a valid choice", normal.kind),
                 domain = NA)
	if (normal.kind == 0L)
            stop("buggy version of Kinderman-Ramage generator is not allowed",
                 domain = NA)
         if(normal.kind == length(n.kinds) - 1L) normal.kind <- -1L
    }
    .Internal(set.seed(seed, i.knd, normal.kind))
}

# Compatibility function to set RNGkind as in a given R version

RNGversion <- function(vstr)
{
    vnum <- as.numeric(strsplit(vstr,".", fixed=TRUE)[[1L]])
    if (length(vnum) < 2L)
	stop("malformed version string")
    if (vnum[1L] == 0 && vnum[2L] < 99)
        RNGkind("Wichmann-Hill", "Buggy Kinderman-Ramage")
    else if (vnum[1L] == 0 || vnum[1L] == 1 && vnum[2L] <= 6)
	RNGkind("Marsaglia-Multicarry", "Buggy Kinderman-Ramage")
    else
	RNGkind("Mersenne-Twister", "Inversion")
}
#  File src/library/base/R/Scripts.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

.Script <-
function(interpreter, script, args, ...)
{
    if(.Platform$OS.type == "windows") {
        cmd <- paste(file.path(R.home("bin"), "Rcmd"),
                     file.path("..", "share", interpreter, script),
                     args)
        system(cmd, invisible = TRUE)
    }
    else
        system(paste(shQuote(file.path(R.home("bin"), "Rcmd")),
                     interpreter,
                     shQuote(file.path(R.home("share"),
                                       interpreter, script)),
                     args),
               ...)
}
#  File src/library/base/R/files.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 2007 R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

.TAOCP1997init <- function(seed)
{
    KK <- 100L; LL <- 37L; MM <- as.integer(2^30)
    KKK <- KK + KK - 1L; KKL <- KK - LL
    ss <- seed - (seed %% 2L) + 2L
    X <- integer(KKK)
    for(j in 1L:KK) {
        X[j] <- ss
        ss <- ss+ss
        if(ss >= MM) ss <- ss - MM + 2L
    }
    X[2L] <- X[2L] + 1L
    ss <- seed
    T <- 69L
    while(T > 0) {
        for(j in KK:2L) X[j + j - 1L] <- X[j]
        for(j in seq(KKK, KKL + 1L, -2L))
            X[KKK - j + 2L] <- X[j] - (X[j] %% 2L)
        for(j in KKK:(KK+1L))
            if(X[j] %% 2L == 1L) {
                X[j - KKL] <- (X[j - KKL] - X[j]) %% MM
                X[j - KK] <- (X[j - KK] - X[j]) %% MM
            }
        if(ss %% 2L == 1L) {
            for(j in KK:1L) X[j + 1L] <- X[j]
            X[1L] <- X[KK + 1L]
            if(X[KK + 1L] %% 2L == 1L)
                X[LL + 1L] <- (X[LL + 1L] - X[KK + 1L]) %% MM
        }
        if(ss) ss <- ss %/% 2L else T <- T - 1L
    }
    rs <- c(X[(LL+1L):KK], X[1L:LL])
    invisible(rs)
}
#  File src/library/base/R/all.equal.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2013 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

all.equal <- function(target, current, ...) UseMethod("all.equal")

all.equal.default <-
    function(target, current, ...)
{
    ## Really a dispatcher given mode() of args :
    ## use data.class as unlike class it does not give "integer"
    if(is.language(target) || is.function(target) || is.environment(target))
	return(all.equal.language(target, current, ...))
    if(is.recursive(target))
	return(all.equal.list(target, current, ...))
    msg <- switch (mode(target),
                   integer = ,
                   complex = ,
                   numeric = all.equal.numeric(target, current, ...),
                   character = all.equal.character(target, current, ...),
                   logical = ,
                   raw = all.equal.raw(target, current, ...),
		   ## assumes that slots are implemented as attributes :
		   S4 = attr.all.equal(target, current, ...),
                   if(data.class(target) != data.class(current)) {
                       gettextf("target is %s, current is %s",
                                data.class(target), data.class(current))
                   } else NULL)
    if(is.null(msg)) TRUE else msg
}

all.equal.numeric <-
    function(target, current, tolerance = .Machine$double.eps ^ .5,
             scale = NULL, check.attributes = TRUE, ...)
{
    msg <- if(check.attributes)
	attr.all.equal(target, current, tolerance=tolerance, scale=scale, ...)
    if(data.class(target) != data.class(current)) {
	msg <- c(msg, paste0("target is ", data.class(target), ", current is ",
                             data.class(current)))
	return(msg)
    }

    lt <- length(target)
    lc <- length(current)
    cplx <- is.complex(target) # and so current must be too.
    if(lt != lc) {
	## *replace* the 'Lengths' msg[] from attr.all.equal():
	if(!is.null(msg)) msg <- msg[- grep("\\bLengths\\b", msg)]
	msg <- c(msg, paste0(if(cplx) "Complex" else "Numeric",
                             ": lengths (", lt, ", ", lc, ") differ"))
	return(msg)
    }
    ## remove atttributes (remember these are both numeric or complex vectors)
    ## one place this is needed is to unclass Surv objects in the rpart test suite.
    target <- as.vector(target)
    current <- as.vector(current)
    out <- is.na(target)
    if(any(out != is.na(current))) {
	msg <- c(msg, paste("'is.NA' value mismatch:", sum(is.na(current)),
			    "in current", sum(out), "in target"))
	return(msg)
    }
    out <- out | target == current
    if(all(out)) { if (is.null(msg)) return(TRUE) else return(msg) }

    target <- target[!out]
    current <- current[!out]
    if(is.integer(target) && is.integer(current)) target <- as.double(target)
    xy <- mean((if(cplx) Mod else abs)(target - current))
    what <-
	if(is.null(scale)) {
	    xn <- mean(abs(target))
	    if(is.finite(xn) && xn > tolerance) {
		xy <- xy/xn
		"relative"
	    } else "absolute"
	} else {
	    xy <- xy/scale
	    "scaled"
	}

    if (cplx) what <- paste(what, "Mod") # PR#10575
    if(is.na(xy) || xy > tolerance)
        msg <- c(msg, paste("Mean", what, "difference:", format(xy)))

    if(is.null(msg)) TRUE else msg
}

all.equal.character <-
    function(target, current, check.attributes = TRUE, ...)
{
    msg <-  if(check.attributes) attr.all.equal(target, current, ...)
    if(data.class(target) != data.class(current)) {
	msg <- c(msg, paste0("target is ", data.class(target), ", current is ",
                             data.class(current)))
	return(msg)
    }
    lt <- length(target)
    lc <- length(current)
    if(lt != lc) {
	if(!is.null(msg)) msg <- msg[- grep("\\bLengths\\b", msg)]
	msg <- c(msg,
                 paste0("Lengths (", lt, ", ", lc,
                        ") differ (string compare on first ",
                        ll <- min(lt, lc), ")"))
	ll <- seq_len(ll)
	target <- target[ll]
	current <- current[ll]
    }
    nas <- is.na(target); nasc <- is.na(current)
    if (any(nas != nasc)) {
	msg <- c(msg, paste("'is.NA' value mismatch:", sum(nasc),
                            "in current", sum(nas), "in target"))
	return(msg)
    }
    ne <- !nas & (target != current)
    if(!any(ne) && is.null(msg)) TRUE
    else if(sum(ne) == 1L) c(msg, paste("1 string mismatch"))
    else if(sum(ne) > 1L) c(msg, paste(sum(ne), "string mismatches"))
    else msg
}

## visible, so need to test both args
all.equal.factor <- function(target, current, check.attributes = TRUE, ...)
{
    if(!inherits(target, "factor"))
	return("'target' is not a factor")
    if(!inherits(current, "factor"))
	return("'current' is not a factor")
    msg <-  if(check.attributes) attr.all.equal(target, current, ...)
    n <- all.equal(as.character(target), as.character(current),
                   check.attributes = check.attributes, ...)
    if(is.character(n)) msg <- c(msg, n)
    if(is.null(msg)) TRUE else msg
}

all.equal.formula <- function(target, current, ...)
{
    ## NB: this assumes the default method for class formula, not
    ## the misquided one in package Formula
    if(length(target) != length(current))
	return(paste("target, current differ in having response: ",
		     length(target) == 3L, ", ",
                     length(current) == 3L, sep=""))
    ## <NOTE>
    ## This takes same-length formulas as all equal if they deparse
    ## identically.  As of 2010-02-24, deparsing strips attributes; if
    ## this is changed, the all equal behavior will change unless the
    ## test is changed.
    ## </NOTE>
    if(!identical(deparse(target), deparse(current)))
	"formulas differ in contents"
    else TRUE
}

all.equal.language <- function(target, current, ...)
{
    mt <- mode(target)
    mc <- mode(current)
    if(mt == "expression" && mc == "expression")
	return(all.equal.list(target, current, ...))
    ttxt <- paste(deparse(target), collapse = "\n")
    ctxt <- paste(deparse(current), collapse = "\n")
    msg <- c(if(mt != mc)
	     paste0("Modes of target, current: ", mt, ", ", mc),
	     if(ttxt != ctxt) {
		 if(pmatch(ttxt, ctxt, 0L))
		     "target is a subset of current"
		 else if(pmatch(ctxt, ttxt, 0L))
		     "current is a subset of target"
		 else "target, current do not match when deparsed"
	     })
    if(is.null(msg)) TRUE else msg
}

all.equal.list <- function(target, current, check.attributes = TRUE, ...)
{
    msg <- if(check.attributes) attr.all.equal(target, current, ...)
    ## Unclass to ensure we get the low-level components
    target <- unclass(target)
    current <- unclass(current)
    iseq <-
	if(length(target) == length(current)) {
	    seq_along(target)
	} else {
            if(!is.null(msg)) msg <- msg[- grep("\\bLengths\\b", msg)]
	    nc <- min(length(target), length(current))
	    msg <- c(msg, paste("Length mismatch: comparison on first",
				nc, "components"))
	    seq_len(nc)
	}
    for(i in iseq) {
	mi <- all.equal(target[[i]], current[[i]], check.attributes = check.attributes, ...)
	if(is.character(mi))
	    msg <- c(msg, paste0("Component ", i, ": ", mi))
    }
    if(is.null(msg)) TRUE else msg
}

## also used for logical
all.equal.raw <-
    function(target, current, check.attributes = TRUE, ...)
{
    msg <-  if(check.attributes) attr.all.equal(target, current, ...)
    if(data.class(target) != data.class(current)) {
	msg <- c(msg, paste0("target is ", data.class(target), ", current is ",
                             data.class(current)))
	return(msg)
    }
    lt <- length(target)
    lc <- length(current)
    if(lt != lc) {
	if(!is.null(msg)) msg <- msg[- grep("\\bLengths\\b", msg)]
	msg <- c(msg, paste0("Lengths (", lt, ", ", lc,
                             ") differ (comparison on first ",
                             ll <- min(lt, lc), " components)"))
	ll <- seq_len(ll)
	target <- target[ll]
	current <- current[ll]
    }
    # raws do not have NAs, but logicals do
    nas <- is.na(target); nasc <- is.na(current)
    if (any(nas != nasc)) {
	msg <- c(msg, paste("'is.NA' value mismatch:", sum(nasc),
                            "in current", sum(nas), "in target"))
	return(msg)
    }
    ne <- !nas & (target != current)
    if(!any(ne) && is.null(msg)) TRUE
    else if(sum(ne) == 1L) c(msg, paste("1 element mismatch"))
    else if(sum(ne) > 1L) c(msg, paste(sum(ne), "element mismatches"))
    else msg
}


## attributes are a pairlist, so never 'long'
attr.all.equal <- function(target, current,
                           check.attributes = TRUE, check.names = TRUE, ...)
{
    ##--- "all.equal(.)" for attributes ---
    ##---  Auxiliary in all.equal(.) methods --- return NULL or character()
    msg <- NULL
    if(mode(target) != mode(current))
	msg <- paste0("Modes: ", mode(target), ", ", mode(current))
    if(length(target) != length(current))
	msg <- c(msg,
                 paste0("Lengths: ", length(target), ", ", length(current)))
    ax <- attributes(target)
    ay <- attributes(current)
    if(check.names) {
        nx <- names(target)
        ny <- names(current)
        if((lx <- length(nx)) | (ly <- length(ny))) {
            ## names() treated now; hence NOT with attributes()
            ax$names <- ay$names <- NULL
            if(lx && ly) {
                if(is.character(m <- all.equal.character(nx, ny, check.attributes = check.attributes)))
                    msg <- c(msg, paste("Names:", m))
            } else if(lx)
                msg <- c(msg, "names for target but not for current")
            else msg <- c(msg, "names for current but not for target")
        }
    } else {
	ax[["names"]] <- NULL
	ay[["names"]] <- NULL
    }
	
    if(check.attributes && (length(ax) || length(ay))) {# some (more) attributes
	## order by names before comparison:
	nx <- names(ax)
	ny <- names(ay)
	if(length(nx)) ax <- ax[order(nx)]
	if(length(ny)) ay <- ay[order(ny)]
	tt <- all.equal(ax, ay, check.attributes = check.attributes, ...)
	if(is.character(tt)) msg <- c(msg, paste("Attributes: <", tt, ">"))
    }
    msg # NULL or character
}
#  File src/library/base/R/allnames.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

all.names <- function(expr, functions = TRUE, max.names = -1L, unique = FALSE)
    .Internal(all.names(expr, functions, max.names, unique))

all.vars <- function(expr, functions = FALSE, max.names = -1L, unique = TRUE)
    .Internal(all.names(expr, functions, max.names, unique))
#  File src/library/base/R/aperm.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

aperm <- function(a, perm, ...) UseMethod("aperm")

aperm.default <- function (a, perm = NULL, resize = TRUE, ...)
     .Internal(aperm(a, perm, resize))

aperm.table <- function(a, perm = NULL, resize = TRUE, keep.class = TRUE, ...)
{
     r <- aperm.default(a, perm, resize=resize)
     if(keep.class) class(r) <- class(a)
     r
}
#  File src/library/base/R/append.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

append <- function (x, values, after = length(x))
{
    lengx <- length(x)
    if (!after) c(values, x)
    else if (after >= lengx) c(x, values)
    else c(x[1L:after], values, x[(after + 1L):lengx])
}
#  File src/library/base/R/apply.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

apply <- function(X, MARGIN, FUN, ...)
{
    FUN <- match.fun(FUN)

    ## Ensure that X is an array object
    dl <- length(dim(X))
    if(!dl) stop("dim(X) must have a positive length")
    if(is.object(X))
	X <- if(dl == 2L) as.matrix(X) else as.array(X)
    ## now record dim as coercion can change it
    ## (e.g. when a data frame contains a matrix).
    d <- dim(X)
    dn <- dimnames(X)
    ds <- seq_len(dl)

    ## Extract the margins and associated dimnames

    if (is.character(MARGIN)) {
        if(is.null(dnn <- names(dn))) # names(NULL) is NULL
           stop("'X' must have named dimnames")
        MARGIN <- match(MARGIN, dnn)
        if (any(is.na(MARGIN)))
            stop("not all elements of 'MARGIN' are names of dimensions")
    }
    s.call <- ds[-MARGIN]
    s.ans  <- ds[MARGIN]
    d.call <- d[-MARGIN]
    d.ans <- d[MARGIN]
    dn.call <- dn[-MARGIN]
    dn.ans <- dn[MARGIN]
    ## dimnames(X) <- NULL

    ## do the calls

    d2 <- prod(d.ans)
    if(d2 == 0L) {
        ## arrays with some 0 extents: return ``empty result'' trying
        ## to use proper mode and dimension:
        ## The following is still a bit `hackish': use non-empty X
        newX <- array(vector(typeof(X), 1L), dim = c(prod(d.call), 1L))
        ans <- FUN(if(length(d.call) < 2L) newX[,1] else
                   array(newX[, 1L], d.call, dn.call), ...)
        return(if(is.null(ans)) ans else if(length(d.ans) < 2L) ans[1L][-1L]
               else array(ans, d.ans, dn.ans))
    }
    ## else
    newX <- aperm(X, c(s.call, s.ans))
    dim(newX) <- c(prod(d.call), d2)
    ans <- vector("list", d2)
    if(length(d.call) < 2L) {# vector
        if (length(dn.call)) dimnames(newX) <- c(dn.call, list(NULL))
        for(i in 1L:d2) {
            tmp <- FUN(newX[,i], ...)
            if(!is.null(tmp)) ans[[i]] <- tmp
        }
    } else
       for(i in 1L:d2) {
           tmp <- FUN(array(newX[,i], d.call, dn.call), ...)
           if(!is.null(tmp)) ans[[i]] <- tmp
        }

    ## answer dims and dimnames

    ans.list <- is.recursive(ans[[1L]])
    l.ans <- length(ans[[1L]])

    ans.names <- names(ans[[1L]])
    if(!ans.list)
	ans.list <- any(unlist(lapply(ans, length)) != l.ans)
    if(!ans.list && length(ans.names)) {
        all.same <- vapply(ans, function(x) identical(names(x), ans.names), NA)
        if (!all(all.same)) ans.names <- NULL
    }
    len.a <- if(ans.list) d2 else length(ans <- unlist(ans, recursive = FALSE))
    if(length(MARGIN) == 1L && len.a == d2) {
	names(ans) <- if(length(dn.ans[[1L]])) dn.ans[[1L]] # else NULL
	return(ans)
    }
    if(len.a == d2)
	return(array(ans, d.ans, dn.ans))
    if(len.a && len.a %% d2 == 0L) {
        if(is.null(dn.ans)) dn.ans <- vector(mode="list", length(d.ans))
        dn.ans <- c(list(ans.names), dn.ans)
	return(array(ans, c(len.a %/% d2, d.ans),
		     if(!all(vapply(dn.ans, is.null, NA))) dn.ans))
    }
    return(ans)
}
#  File src/library/base/R/array.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

array <-
function(data = NA, dim = length(data), dimnames = NULL)
{
    ## allow for as.vector.factor (converts to character)
    if(is.atomic(data) && !is.object(data))
        return(.Internal(array(data, dim, dimnames)))
    data <- as.vector(data)
    ## package rv has an as.vector() method which leave this as a classed list
    if(is.object(data)) {
        dim <- as.integer(dim)
        if (!length(dim)) stop("'dims' cannot be of length 0")
        vl <- prod(dim)
        if(length(data) != vl) {
            ## C code allows long vectors, but rep() does not.
            if(vl > .Machine$integer.max)
                stop("'dim' specifies too large an array")
            data <- rep_len(data, vl)
        }
        if(length(dim)) dim(data) <- dim
        if(is.list(dimnames) && length(dimnames)) dimnames(data) <- dimnames
        data
    } else .Internal(array(data, dim, dimnames))
}

slice.index <-
function(x, MARGIN)
{
    d <- dim(x)
    if(is.null(d))
        d <- length(x)
    n <- length(d)

    if((length(MARGIN) > 1L) || (MARGIN < 1L) || (MARGIN > n))
        stop("incorrect value for 'MARGIN'")

    if(any(d == 0L)) return(array(integer(), d))

    y <- rep.int(rep.int(1L:d[MARGIN],
			 prod(d[seq_len(MARGIN - 1L)]) * rep.int(1L, d[MARGIN])),
		 prod(d[seq.int(from = MARGIN + 1L, length.out = n - MARGIN)]))
    dim(y) <- d
    y
}

provideDimnames <- function(x, sep="", base = list(LETTERS)) {
    ## provide dimnames where missing - not copying x unnecessarily
    dx <- dim(x)
    dnx <- dimnames(x)
    if(new <- is.null(dnx))
	dnx <- vector("list", length(dx))
    k <- length(M <- vapply(base, length, 1L))
    for(i in which(vapply(dnx, is.null, NA))) {
	ii <- 1L+(i-1L) %% k # recycling
	dnx[[i]] <-
	    make.unique(base[[ii]][1L+ 0:(dx[i]-1L) %% M[ii]],
			sep = sep)
	new <- TRUE
    }
    if(new) dimnames(x) <- dnx
    x
}
#  File src/library/base/R/as.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

as.single <- function(x,...) UseMethod("as.single")
as.single.default <- function(x,...)
    structure(.Internal(as.vector(x,"double")), Csingle=TRUE)

# as.character is now internal.  The default method remains here to
# preserve the semantics that for a call with an object argument
# dispatching is done first on as.character and then on as.vector.
as.character.default <- function(x,...) .Internal(as.vector(x, "character"))

as.expression <- function(x,...) UseMethod("as.expression")
as.expression.default <- function(x,...) .Internal(as.vector(x, "expression"))

as.list <- function(x,...) UseMethod("as.list")
## This if() avoid dispatch on methods for as.vector.
as.list.default <- function (x, ...)
    if (typeof(x) == "list") x else .Internal(as.vector(x, "list"))

as.list.function <- function (x, ...) c(formals(x), list(body(x)))

## FIXME:  Really the above  as.vector(x, "list")  should work for data.frames!
as.list.data.frame <- function(x,...) {
    x <- unclass(x)
    attr(x,"row.names") <- NULL
    x
}

as.list.environment <- function(x, all.names=FALSE, ...)
    .Internal(env2list(x, all.names))

## NB: as.vector is used for several other as.xxxx, including
## as.expression, as.list, as.pairlist, as.single, as.symbol.
## as.vector dispatches internally so no need for a generic
as.vector <- function(x, mode = "any") .Internal(as.vector(x, mode))

as.matrix <- function(x, ...) UseMethod("as.matrix")
as.matrix.default <- function(x, ...) {
    if (is.matrix(x)) x
    else
	array(x, c(length(x), 1L),
	      if(!is.null(names(x))) list(names(x), NULL) else NULL)
}
as.null <- function(x,...) UseMethod("as.null")
as.null.default <- function(x,...) NULL

as.function <- function(x,...) UseMethod("as.function")
as.function.default <- function (x, envir = parent.frame(), ...)
    if (is.function(x)) x else .Internal(as.function.default(x, envir))

as.array <- function(x, ...) UseMethod("as.array")
as.array.default <- function(x, ...)
{
    if(is.array(x)) return(x)
    n <- names(x)
    dim(x) <- length(x)
    if(length(n)) dimnames(x) <- list(n)
    return(x)
}

as.symbol <- function(x) .Internal(as.vector(x, "symbol"))
as.name <- as.symbol
## would work too: as.name <- function(x) .Internal(as.vector(x, "name"))

as.qr <- function(x) stop("you cannot be serious", domain = NA)
#  File src/library/base/R/assign.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

assign <-
    function (x, value, pos = -1, envir = as.environment(pos),
              inherits = FALSE, immediate = TRUE)
    .Internal(assign(x, value, envir, inherits))

## do_list2env in ../../../main/envir.c
list2env <- function(x, envir = NULL, parent = parent.frame(),
		     hash = (length(x) > 100), size = max(29L, length(x)))
{
    if (is.null(envir)) envir <- new.env(hash=hash, parent=parent, size=size)
    .Internal(list2env(x, envir))
}
#  File src/library/base/R/attach.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

attach <- function(what, pos = 2L, name = deparse(substitute(what)),
                   warn.conflicts = TRUE)
{
    checkConflicts <- function(env)
    {
        dont.mind <- c("last.dump", "last.warning", ".Last.value",
                       ".Random.seed", ".Last.lib", ".onDetach",
                       ".packageName", ".noGenerics", ".required",
                       ".no_S3_generics", ".requireCachedGenerics")
        sp <- search()
        for (i in seq_along(sp)) {
            if (identical(env, as.environment(i))) {
                db.pos <- i
                break
            }
        }
        ob <- objects(db.pos, all.names = TRUE)
        if(.isMethodsDispatchOn()) {
            ## <FIXME>: this is wrong-headed: see library().
            these <- objects(db.pos, all.names = TRUE)
            these <- these[substr(these, 1L, 6L) == ".__M__"]
            gen <- gsub(".__M__(.*):([^:]+)", "\\1", these)
            from <- gsub(".__M__(.*):([^:]+)", "\\2", these)
            gen <- gen[from != ".GlobalEnv"]
            ob <- ob[!(ob %in% gen)]
        }
        ipos <- seq_along(sp)[-c(db.pos, match(c("Autoloads", "CheckExEnv"), sp, 0L))]
        for (i in ipos) {
            obj.same <- match(objects(i, all.names = TRUE), ob, nomatch = 0L)
            if (any(obj.same > 0L)) {
                same <- ob[obj.same]
                same <- same[!(same %in% dont.mind)]
                Classobjs <- grep("^\\.__", same)
                if(length(Classobjs)) same <- same[-Classobjs]
                ## report only objects which are both functions or
                ## both non-functions.
                is_fn1 <- sapply(same, function(x)
                                 exists(x, where = i, mode = "function",
                                        inherits = FALSE))
                is_fn2 <- sapply(same, function(x)
                                 exists(x, where = db.pos, mode = "function",
                                        inherits = FALSE))
                same <- same[is_fn1 == is_fn2]
                if(length(same)) {
                    objs <- strwrap(paste(same, collapse=", "),
                                    indent = 4L, exdent = 4L)
                    pkg <-
                        if (sum(sp == sp[i]) > 1L) {
                            sprintf("%s (position %d)", sp[i], i)
                        } else {
                            sp[i]
                        }
                    msg <- sprintf(ngettext(length(same),
                                            "The following object is masked %s %s:\n\n%s\n",
                                            "The following objects are masked %s %s:\n\n%s\n"),
                                   if (i < db.pos) "_by_" else "from",
                                   pkg, paste(objs, collapse="\n"))
                    cat(msg)
                }
            }
        }
    }

    if(pos == 1L) {
        warning("*** 'pos=1' is not possible; setting 'pos=2' for now.\n",
                "*** Note that 'pos=1' will give an error in the future")
        pos <- 2L
    }
    if (is.character(what) && (length(what) == 1L)){
        if (!file.exists(what))
            stop(gettextf("file '%s' not found", what), domain = NA)
        if(missing(name)) name <- paste0("file:", what)
        value <- .Internal(attach(NULL, pos, name))
        load(what, envir = as.environment(pos))
    }
    else
        value <- .Internal(attach(what, pos, name))
    if(warn.conflicts &&
       !exists(".conflicts.OK", envir = value, inherits = FALSE)) {
        checkConflicts(value)
    }
    if( length(ls(envir = value, all.names = TRUE)) && .isMethodsDispatchOn() )
        methods:::cacheMetaData(value, TRUE)
    invisible(value)
}

detach <- function(name, pos = 2L, unload = FALSE, character.only = FALSE,
                   force = FALSE)
{
    if(!missing(name)) {
	if(!character.only) name <- substitute(name)
	pos <-
	    if(is.numeric(name)) name
	    else {
                if (!is.character(name)) name <- deparse(name)
                match(name, search())
            }
	if(is.na(pos)) stop("invalid 'name' argument")
    }

    packageName <- search()[[pos]]

    ## we need to treat packages differently from other objects, so get those
    ## out of the way now
    if (! grepl("^package:", packageName) )
        return(invisible(.Internal(detach(pos))))

    ## From here down we are detaching a package.
    pkgname <- sub("^package:", "", packageName)
    for(pkg in search()[-1L]) {
        if(grepl("^package:", pkg) &&
           exists(".Depends", pkg, inherits = FALSE) &&
           pkgname %in% get(".Depends", pkg, inherits = FALSE))
            if(force)
                warning(gettextf("package %s is required by %s, which may no longer work correctly",
                                 sQuote(pkgname), sQuote(sub("^package:", "", pkg))),
                     call. = FALSE, domain = NA)
            else
                stop(gettextf("package %s is required by %s so will not be detached",
                              sQuote(pkgname), sQuote(sub("^package:", "", pkg))),
                     call. = FALSE, domain = NA)
    }
    env <- as.environment(pos)
    libpath <- attr(env, "path")
    hook <- getHook(packageEvent(pkgname, "detach")) # might be a list
    for(fun in rev(hook)) try(fun(pkgname, libpath))
    ## some people, e.g. package g.data, have faked pakages without namespaces
    ns <- .getNamespace(pkgname)
    if(!is.null(ns) &&
       exists(".onDetach", mode = "function", where = ns, inherits = FALSE)) {
        .onDetach <- get(".onDetach",  mode = "function", pos = ns,
                         inherits = FALSE)
        if(!is.null(libpath)) {
            res <- tryCatch(.onDetach(libpath), error = identity)
            if (inherits(res, "error")) {
                warning(gettextf("%s failed in %s() for '%s', details:\n  call: %s\n  error: %s",
                                 ".onDetach", "detach", pkgname,
                                 deparse(conditionCall(res))[1L],
                                 conditionMessage(res)),
                        call. = FALSE, domain = NA)
            }
        }
    }
    else if(exists(".Last.lib", mode = "function", where = pos, inherits = FALSE)) {
        .Last.lib <- get(".Last.lib",  mode = "function", pos = pos,
                         inherits = FALSE)
        if(!is.null(libpath)) {
            res <- tryCatch(.Last.lib(libpath), error = identity)
            if (inherits(res, "error")) {
                warning(gettextf("%s failed in %s() for '%s', details:\n  call: %s\n  error: %s",
                                 ".Last.lib", "detach", pkgname,
                                 deparse(conditionCall(res))[1L],
                                 conditionMessage(res)),
                        call. = FALSE, domain = NA)
            }
        }
    }
    .Internal(detach(pos))

    if(pkgname %in% loadedNamespaces()) {
        ## the lazyload DB is flushed when the namespace is unloaded
        if(unload) {
            tryCatch(unloadNamespace(pkgname),
                     error = function(e)
                     warning(gettextf("%s namespace cannot be unloaded:\n  ",
                                      sQuote(pkgname)),
                             conditionMessage(e),
                             call. = FALSE, domain = NA))
        }
    } else {
        if(.isMethodsDispatchOn() && methods:::.hasS4MetaData(env))
            methods:::cacheMetaData(env, FALSE)
        .Internal(lazyLoadDBflush(paste0(libpath, "/R/", pkgname, ".rdb")))
    }
    invisible()
}

.detach <- function(pos) .Internal(detach(pos))

ls <- objects <-
    function (name, pos = -1L, envir = as.environment(pos), all.names = FALSE,
              pattern)
{
    if (!missing(name)) {
        nameValue <- try(name, silent = TRUE)
        if(identical(class(nameValue), "try-error")) {
            name <- substitute(name)
            if (!is.character(name))
                name <- deparse(name)
            warning(gettextf("%s converted to character string", sQuote(name)),
                    domain = NA)
            pos <- name
        }
        else
            pos <- nameValue
    }
    all.names <- .Internal(ls(envir, all.names))
    if (!missing(pattern)) {
        if ((ll <- length(grep("[", pattern, fixed = TRUE))) &&
            ll != length(grep("]", pattern, fixed = TRUE))) {
            if (pattern == "[") {
                pattern <- "\\["
                warning("replaced regular expression pattern '[' by  '\\\\['")
            }
            else if (length(grep("[^\\\\]\\[<-", pattern))) {
                pattern <- sub("\\[<-", "\\\\\\[<-", pattern)
                warning("replaced '[<-' by '\\\\[<-' in regular expression pattern")
            }
        }
        grep(pattern, all.names, value = TRUE)
    }
    else all.names
}
#  File src/library/base/R/attr.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

`mostattributes<-` <- function(obj, value)
{
    if(length(value)) {
	if(!is.list(value)) stop("'value' must be a list")
	if(h.nam <- !is.na(inam <- match("names", names(value)))) {
	    n1 <- value[[inam]];	value <- value[-inam] }
	if(h.dim <- !is.na(idin <- match("dim", names(value)))) {
	    d1 <- value[[idin]];	value <- value[-idin] }
	if(h.dmn <- !is.na(idmn <- match("dimnames", names(value)))) {
	    dn1 <- value[[idmn]];	value <- value[-idmn] }
	attributes(obj) <- value
        dm <- attr(obj, "dim")
	## for list-like objects with a length() method, e.g. POSIXlt
	L <- length(if(is.list(obj)) unclass(obj) else obj)
        ## Be careful to set dim before dimnames.
	if(h.dim && L == prod(d1)) attr(obj, "dim") <- dm <- d1
	if(h.dmn && !is.null(dm)) {
            ddn <- sapply(dn1, length)
            if( all((dm == ddn)[ddn > 0]) ) attr(obj, "dimnames") <- dn1
        }
        ## don't set if it has 'dim' now
	if(h.nam && is.null(dm) && L == length(n1)) attr(obj, "names") <- n1
    }
    obj
}
#  File src/library/base/R/autoload.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

autoload <- function(name, package, reset=FALSE, ...)
{
    if (!reset && exists(name, envir = .GlobalEnv, inherits = FALSE))
	stop("an object with that name already exists")
    m <- match.call()
    m[[1L]] <- as.name("list")
    newcall <- eval(m, parent.frame())
    newcall <- as.call(c(as.name("autoloader"), newcall))
    newcall$reset <- NULL
    if (is.na(match(package, .Autoloaded)))
	assign(".Autoloaded", c(package, .Autoloaded), envir =.AutoloadEnv)
    do.call("delayedAssign", list(name, newcall, .GlobalEnv, .AutoloadEnv))
    ## no longer return the result, which is a promise
    invisible()
}

autoloader <- function (name, package, ...)
{
    name <- paste0(name, "")
    rm(list = name, envir = .AutoloadEnv, inherits = FALSE)
    m <- match.call()
    m$name <- NULL
    m[[1L]] <- as.name("library")
    ## load the package
    eval(m, .GlobalEnv)
    ## reset the autoloader
    autoload(name, package, reset = TRUE, ...)
    ## reevaluate the object
    where <- match(paste("package", package, sep = ":"), search())
    if (exists(name, where = where, inherits = FALSE))
	eval(as.name(name), as.environment(where))
    else
	stop(gettextf("autoloader did not find '%s' in '%s'", name, package),
             domain = NA)
}
#  File src/library/base/R/backquote.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

## quote() is .Primitive

### PR15077: need to substitute in a length-one pairlist, so
### handle pairlists first
bquote <- function(expr, where=parent.frame())
{
    unquote <- function(e)
        if (is.pairlist(e)) as.pairlist(lapply(e, unquote))
        else if (length(e) <= 1L) e
        else if (e[[1L]] == as.name(".")) eval(e[[2L]], where)
        else as.call(lapply(e, unquote))

    unquote(substitute(expr))
}

## utility we've used ourselves
enquote <- function(cl) as.call(list(as.name("quote"), cl))
#  File src/library/base/R/backsolve.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

forwardsolve <-
    function(l, x, k = ncol(l), upper.tri = FALSE, transpose = FALSE)
{
    l <- as.matrix(l)
    x.mat <- is.matrix(x)
    if(!x.mat) x <- as.matrix(x)
    z <- .Internal(backsolve(l, x, k, upper.tri, transpose))
    if(x.mat) z else drop(z)
}

backsolve <- function(r, x, k  = ncol(r), upper.tri = TRUE, transpose = FALSE)
{
    r <- as.matrix(r) # so ncol(r) works
    x.mat <- is.matrix(x)
    if(!x.mat) x <- as.matrix(x)
    z <- .Internal(backsolve(r, x, k, upper.tri, transpose))
    if(x.mat) z else drop(z)
}
#  File src/library/base/R/bindenv.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

lockEnvironment <- function(env, bindings = FALSE)
    .Internal(lockEnvironment(env, bindings))

environmentIsLocked <- function(env)
    .Internal(environmentIsLocked(env))

lockBinding <- function(sym, env) {
    if (is.character(sym)) sym <- as.name(sym)
    .Internal(lockBinding(sym, env))
}

bindingIsLocked <- function(sym, env) {
    if (is.character(sym)) sym <- as.name(sym)
    .Internal(bindingIsLocked(sym, env))
}

makeActiveBinding <- function(sym, fun, env) {
    if (is.character(sym)) sym <- as.name(sym)
    .Internal(makeActiveBinding(sym, fun, env))
}

bindingIsActive <- function(sym, env) {
    if (is.character(sym)) sym <- as.name(sym)
    .Internal(bindingIsActive(sym, env))
}

unlockBinding <- function(sym, env) {
    if (is.character(sym)) sym <- as.name(sym)
    .Internal(unlockBinding(sym, env))
}
#  File src/library/base/R/octhex.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

bitwNot <- function(a) .Internal(bitwiseNot(a))

bitwAnd <- function(a, b) .Internal(bitwiseAnd(a, b))

bitwOr <- function(a, b) .Internal(bitwiseOr(a, b))

bitwXor <- function(a, b) .Internal(bitwiseXor(a, b))

bitwShiftL <- function(a, n) .Internal(bitwiseShiftL(a, n))

bitwShiftR <- function(a, n) .Internal(bitwiseShiftR(a, n))
#  File src/library/base/R/builtins.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

builtins <- function(internal=FALSE)
    .Internal(builtins(internal))
#  File src/library/base/R/by.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

by <- function(data, INDICES, FUN, ..., simplify = TRUE) UseMethod("by")

## prior to 2.7.0 this promoted vectors to data frames, but
## the data frame method dropped to a single column.
by.default <- function(data, INDICES, FUN, ..., simplify = TRUE)
{
    dd <- as.data.frame(data)
    if(length(dim(data)))
        by(dd, INDICES, FUN, ..., simplify = simplify)
    else {
        if(!is.list(INDICES)) {        # record the names for print.by
            IND <- vector("list", 1L)
            IND[[1L]] <- INDICES
            names(IND) <- deparse(substitute(INDICES))[1L]
        } else IND <- INDICES
        FUNx <- function(x) FUN(dd[x,], ...)
        nd <- nrow(dd)
	structure(eval(substitute(tapply(seq_len(nd), IND, FUNx,
				      simplify = simplify)), dd),
		  call = match.call(),
		  class = "by")
    }
}

by.data.frame <- function(data, INDICES, FUN, ..., simplify = TRUE)
{
    if(!is.list(INDICES)) { # record the names for print.by
        IND <- vector("list", 1L)
        IND[[1L]] <- INDICES
        names(IND) <- deparse(substitute(INDICES))[1L]
    } else IND <- INDICES
    FUNx <- function(x) FUN(data[x,, drop=FALSE], ...) # (PR#10506)
    nd <- nrow(data)
    structure(eval(substitute(tapply(seq_len(nd), IND, FUNx,
				     simplify = simplify)), data),
	      call = match.call(),
	      class = "by")
}

print.by <- function(x, ..., vsep)
{
    d <- dim(x)
    dn <- dimnames(x)
    dnn <- names(dn)
    if(missing(vsep))
        vsep <- paste(rep.int("-", 0.75*getOption("width")), collapse = "")
    lapply(X = seq_along(x), FUN = function(i, x, vsep, ...) {
        if(i != 1L && !is.null(vsep)) cat(vsep, "\n")
        ii <- i - 1L
        for(j in seq_along(dn)) {
            iii <- ii %% d[j] + 1L; ii <- ii %/% d[j]
            cat(dnn[j], ": ", dn[[j]][iii], "\n", sep = "")
        }
        print(x[[i]], ...)
    } , x, vsep, ...)
    invisible(x)
}
#  File src/library/base/R/callCC.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

callCC <- function(fun) {
    value <- NULL
    delayedAssign("throw", return(value))
    fun(function(v) { value <<- v; throw })
}
#  File src/library/base/R/cat.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

cat <- function(..., file = "", sep = " ", fill = FALSE,
                labels = NULL, append = FALSE)
{
    if(is.character(file))
        if(file == "") file <- stdout()
        else if(substring(file, 1L, 1L) == "|") {
            file <- pipe(substring(file, 2L), "w")
            on.exit(close(file))
        } else {
            file <- file(file, ifelse(append, "a", "w"))
            on.exit(close(file))
        }
    .Internal(cat(list(...), file, sep, fill, labels, append))
}
#  File src/library/base/R/character.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

substr <- function(x, start, stop)
{
    if(!is.character(x)) x <- as.character(x)
    .Internal(substr(x, as.integer(start), as.integer(stop)))
}

substring <- function(text, first, last=1000000L)
{
    if(!is.character(text)) text <- as.character(text)
    n <- max(lt <- length(text), length(first), length(last))
    if(lt && lt < n) text <- rep_len(text, length.out = n)
    .Internal(substr(text, as.integer(first), as.integer(last)))
}

`substr<-` <- function(x, start, stop, value)
    .Internal(`substr<-`(x, as.integer(start), as.integer(stop), value))

`substring<-` <- function(text, first, last=1000000L, value)
    .Internal(`substr<-`(text, as.integer(first), as.integer(last), value))

abbreviate <-
    function(names.arg, minlength = 4L, use.classes = TRUE, dot = FALSE,
             strict = FALSE, method = c("left.kept", "both.sides"))
{
    ## we just ignore use.classes
    if(minlength <= 0L)
	return(rep.int("", length(names.arg)))
    ## need to remove leading/trailing spaces before we check for dups
    ## This is inefficient but easier than modifying do_abbrev (=> FIXME !)
    names.arg <- sub("^ +", "", sub(" +$", "", as.character(names.arg)))
    dups <- duplicated(names.arg)
    old <- names.arg
    if(any(dups))
	names.arg <- names.arg[!dups]
    x <- names.arg
    if(strict) {
	x[] <- .Internal(abbreviate(x, minlength, use.classes))
    } else {
	method <- match.arg(method)
	if(method == "both.sides")
	    ## string reversion: FIXME reverse .Internal(abbreviate(.))
	    chRev <- function(x)
		sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
	dup2 <- rep.int(TRUE, length(names.arg))
	these <- names.arg
	repeat {
	    ans <- .Internal(abbreviate(these, minlength, use.classes))
	    ## NB: fulfills   max(nchar(ans)) <= minlength
	    x[dup2] <- ans
	    if(!any(dup2 <- duplicated(x))) break
	    if(method == "both.sides") { ## abbreviate the dupl. ones from the other side:
		x[dup2] <- chRev(.Internal(abbreviate(chRev(names.arg[dup2]),
						      minlength, use.classes)))
		if(!any(dup2 <- duplicated(x))) break
	    }
	    minlength <- minlength+1
	    dup2 <- dup2 | match(x, x[dup2], 0L)
	    these <- names.arg[dup2]
	}
    }
    if(any(dups))
	x <- x[match(old,names.arg)]
    if(dot) {			    # add "." where we did abbreviate:
	chgd <- x != old
	x[chgd] <- paste0(x[chgd],".")
    }
    names(x) <- old
    x
}

make.names <- function(names, unique = FALSE, allow_ = TRUE)
{
    names <- .Internal(make.names(as.character(names), allow_))
    if(unique) names <- make.unique(names)
    names
}

make.unique <- function (names, sep = ".") .Internal(make.unique(names, sep))

chartr <- function(old, new, x)
{
    if(!is.character(x)) x <- as.character(x)
    .Internal(chartr(old, new, x))
}
tolower <- function(x)
{
    if(!is.character(x)) x <- as.character(x)
    .Internal(tolower(x))
}
toupper <- function(x)
{
    if(!is.character(x)) x <- as.character(x)
    .Internal(toupper(x))
}

casefold <- function(x, upper = FALSE)
    if(upper) toupper(x) else tolower(x)

sQuote <- function(x)
{
    if (!length(x)) return(character())
    before <- after <- "'"
    q <- getOption("useFancyQuotes")
    if(!is.null(q)) {
        if(identical(q, TRUE)) {
            li <- l10n_info()
            if(li$"UTF-8") q <- "UTF-8"
            if(!is.null(li$codepage) && li$codepage > 0L) {
                ## we can't just use iconv, as that seems to think
                ## it is in latin1 in CP1252
                if(li$codepage >= 1250L && li$codepage <= 1258L
                   || li$codepage == 874L) {
                    before <- "\x91"; after <- "\x92"
                } else {
                    z <- iconv(c("\xe2\x80\x98", "\xe2\x80\x99"), "UTF-8", "")
                    before <- z[1L]; after <- z[2L]
                }
            }
        }
        if(identical(q, "TeX")) {
            before <- "`"; after <- "'"
        }
        if(identical(q, "UTF-8")) {
            before <- "\xe2\x80\x98"; after <- "\xe2\x80\x99"
        }
        if(is.character(q) && length(q) >= 4L) {
            before <- q[1L]; after <- q[2L]
        }
        ## we do not want these strings marked as in the encoding
        ## R was built under
        Encoding(before) <- Encoding(after) <- "unknown"
    }
    paste0(before, x, after)
}

dQuote <- function(x)
{
    if (!length(x)) return(character())
    before <- after <- "\""
    q <- getOption("useFancyQuotes")
    if(!is.null(q)) {
        if(identical(q, TRUE)) {
            li <- l10n_info()
            if(li$"UTF-8") q <- "UTF-8"
            if(!is.null(li$codepage) && li$codepage > 0L) {
                if(li$codepage >= 1250L && li$codepage <= 1258L
                    || li$codepage == 874L) {
                    before <- "\x93"; after <- "\x94"
                } else {
                    z <- iconv(c("\xe2\x80\x9c", "\xe2\x80\x9d"), "UTF-8", "")
                    before <- z[1L]; after <- z[2L]
                }
            }
        }
        if(identical(q, "TeX")) {
            before <- "``"; after <- "''"
        }
        if(identical(q, "UTF-8")) {
            before <- "\xe2\x80\x9c"; after <- "\xe2\x80\x9d"
        }
        if(is.character(q) && length(q) >= 4L) {
            before <- q[3L]; after <- q[4L]
        }
        Encoding(before) <- Encoding(after) <- "unknown"
    }
    paste0(before, x, after)
}

strtoi <-
function(x, base = 0L)
    .Internal(strtoi(as.character(x), as.integer(base)))
#  File src/library/base/R/chol.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2013 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

chol <- function(x, ...) UseMethod("chol")

chol.default <- function(x, pivot = FALSE, LINPACK = FALSE, tol = -1, ...)
{
    if (is.complex(x))
        stop("complex matrices not permitted at present")

    if (!pivot && LINPACK) {
        warning("LINPACK = TRUE, pivot = FALSE is defunct and will be ignored",
                domain = NA)
        LINPACK <- FALSE
    }

    if(!LINPACK) return(.Internal(La_chol(as.matrix(x), pivot, tol)))

    ## deprecated in R 3.0.0
    warning("LINPACK = pivot = TRUE is deprecated", domain = NA)

    if(is.matrix(x)) {
	if(nrow(x) != ncol(x)) stop("non-square matrix in 'chol'")
	n <- nrow(x)
    } else {
	if(length(x) != 1L) stop("non-matrix argument to 'chol'")
	n <- 1L
    }

    ## sanity checks
    n <- as.integer(n)
    if(is.na(n)) stop("invalid nrow(x)")
    if(n > 46340) stop("too large a matrix for LINPACK")

    storage.mode(x) <- "double"

    xx <- x
    xx[lower.tri(xx)] <- 0
    z <- .Fortran(.F_dchdc, x = xx, n, n, double(n), piv = integer(n),
                  as.integer(pivot), rank = integer(1L), DUP = FALSE)
    if (z$rank < n)
        if(!pivot) stop("matrix not positive definite")
        else warning("matrix not positive definite")
    robj <- z$x
    if (pivot) {
        attr(robj, "pivot") <- z$piv
        attr(robj, "rank") <- z$rank
        if(!is.null(cn <- colnames(x)))
            colnames(robj) <- cn[z$piv]
    }
    robj
}

chol2inv <- function(x, size = NCOL(x), LINPACK = FALSE)
{
    if(LINPACK) warning("LINPACK = TRUE is defunct", domain = NA)
    .Internal(La_chol2inv(x, size))
}
#  File src/library/base/R/colSums.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


## NB: we now have  implicitGeneric() on these,
##     in ../../methods/R/makeBasicFunsList.R

colSums <- function(x, na.rm = FALSE, dims = 1L)
{
    if(is.data.frame(x)) x <- as.matrix(x)
    if(!is.array(x) || length(dn <- dim(x)) < 2L)
        stop("'x' must be an array of at least two dimensions")
    if(dims < 1L || dims > length(dn) - 1L)
        stop("invalid 'dims'")
    n <- prod(dn[1L:dims])
    dn <- dn[-(1L:dims)]
    z <- if(is.complex(x))
        .Internal(colSums(Re(x), n, prod(dn), na.rm)) +
            1i * .Internal(colSums(Im(x), n, prod(dn), na.rm))
    else .Internal(colSums(x, n, prod(dn), na.rm))
    if(length(dn) > 1L) {
        dim(z) <- dn
        dimnames(z) <- dimnames(x)[-(1L:dims)]
    } else names(z) <- dimnames(x)[[dims+1]]
    z
}

colMeans <- function(x, na.rm = FALSE, dims = 1L)
{
    if(is.data.frame(x)) x <- as.matrix(x)
    if(!is.array(x) || length(dn <- dim(x)) < 2L)
        stop("'x' must be an array of at least two dimensions")
    if(dims < 1L || dims > length(dn) - 1L)
        stop("invalid 'dims'")
    n <- prod(dn[1L:dims])
    dn <- dn[-(1L:dims)]
    z <- if(is.complex(x))
        .Internal(colMeans(Re(x), n, prod(dn), na.rm)) +
            1i * .Internal(colMeans(Im(x), n, prod(dn), na.rm))
    else .Internal(colMeans(x, n, prod(dn), na.rm))
    if(length(dn) > 1L) {
        dim(z) <- dn
        dimnames(z) <- dimnames(x)[-(1L:dims)]
    } else names(z) <- dimnames(x)[[dims+1]]
    z
}

rowSums <- function(x, na.rm = FALSE, dims = 1L)
{
    if(is.data.frame(x)) x <- as.matrix(x)
    if(!is.array(x) || length(dn <- dim(x)) < 2L)
        stop("'x' must be an array of at least two dimensions")
    if(dims < 1L || dims > length(dn) - 1L)
        stop("invalid 'dims'")
    p <- prod(dn[-(1L:dims)])
    dn <- dn[1L:dims]
    z <- if(is.complex(x))
        .Internal(rowSums(Re(x), prod(dn), p, na.rm)) +
            1i * .Internal(rowSums(Im(x), prod(dn), p, na.rm))
    else .Internal(rowSums(x, prod(dn), p, na.rm))
    if(length(dn) > 1L) {
        dim(z) <- dn
        dimnames(z) <- dimnames(x)[1L:dims]
    } else  names(z) <- dimnames(x)[[1L]]
    z
}

rowMeans <- function(x, na.rm = FALSE, dims = 1L)
{
    if(is.data.frame(x)) x <- as.matrix(x)
    if(!is.array(x) || length(dn <- dim(x)) < 2L)
        stop("'x' must be an array of at least two dimensions")
    if(dims < 1L || dims > length(dn) - 1L)
        stop("invalid 'dims'")
    p <- prod(dn[-(1L:dims)])
    dn <- dn[1L:dims]
    z <- if(is.complex(x))
        .Internal(rowMeans(Re(x), prod(dn), p, na.rm)) +
            1i * .Internal(rowMeans(Im(x), prod(dn), p, na.rm))
    else .Internal(rowMeans(x, prod(dn), p, na.rm))
    if(length(dn) > 1L) {
        dim(z) <- dn
        dimnames(z) <- dimnames(x)[1L:dims]
    } else  names(z) <- dimnames(x)[[1L]]
    z
}

.colSums <- function(X, m, n, na.rm = FALSE)
    .Internal(colSums(X, m, n, na.rm))
.colMeans <- function(X, m, n, na.rm = FALSE)
    .Internal(colMeans(X, m, n, na.rm))

.rowSums <- function(X, m, n, na.rm = FALSE)
    .Internal(rowSums(X, m, n, na.rm))
.rowMeans <- function(X, m, n, na.rm = FALSE)
    .Internal(rowMeans(X, m, n, na.rm))
#  File src/library/base/R/conditions.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

##
## Handling Conditions
##

## CARE:  try() in ./New-Internal.R  depends on *internal* coding of tryCatch()!
## ----   If you change this, be sure to adapt  try().
tryCatch <- function(expr, ..., finally) {
    tryCatchList <- function(expr, names, parentenv, handlers) {
	nh <- length(names)
	if (nh > 1L)
	    tryCatchOne(tryCatchList(expr, names[-nh], parentenv,
                                     handlers[-nh]),
			names[nh], parentenv, handlers[[nh]])
	else if (nh == 1L)
	    tryCatchOne(expr, names, parentenv, handlers[[1L]])
	else expr
    }
    tryCatchOne <- function(expr, name, parentenv, handler) {
	doTryCatch <- function(expr, name, parentenv, handler) {
	    .Internal(.addCondHands(name, list(handler), parentenv,
				    environment(), FALSE))
	    expr
	}
	value <- doTryCatch(return(expr), name, parentenv, handler)
	# The return in the call above will exit withOneRestart unless
	# the handler is invoked; we only get to this point if the handler
	# is invoked.  If we get here then the handler will have been
	# popped off the internal handler stack.
	if (is.null(value[[1L]])) {
	    # a simple error; message is stored internally
	    # and call is in result; this defers all allocs until
	    # after the jump
	    msg <- .Internal(geterrmessage())
	    call <- value[[2L]]
	    cond <- simpleError(msg, call)
	}
	else cond <- value[[1L]]
	value[[3L]](cond)
    }
    if (! missing(finally))
        on.exit(finally)
    handlers <- list(...)
    classes <- names(handlers)
    parentenv <- parent.frame()
    if (length(classes) != length(handlers))
        stop("bad handler specification")
    tryCatchList(expr, classes, parentenv, handlers)
}

withCallingHandlers <- function(expr, ...) {
    handlers <- list(...)
    classes <- names(handlers)
    parentenv <- parent.frame()
    if (length(classes) != length(handlers))
        stop("bad handler specification")
    .Internal(.addCondHands(classes, handlers, parentenv, NULL, TRUE))
    expr
}

suppressWarnings <- function(expr) {
    ops <- options(warn = -1) ## FIXME: temporary hack until R_tryEval
    on.exit(options(ops))     ## calls are removed from methods code
    withCallingHandlers(expr,
                        warning=function(w)
                            invokeRestart("muffleWarning"))
}


##
## Conditions and Condition Signaling
##

simpleCondition <- function(message, call = NULL) {
    class <- c("simpleCondition", "condition")
    structure(list(message=as.character(message), call = call), class=class)
}

simpleError <- function(message, call = NULL) {
    class <- c("simpleError", "error", "condition")
    structure(list(message=as.character(message), call = call), class=class)
}

simpleWarning <- function(message, call = NULL) {
    class <- c("simpleWarning", "warning", "condition")
    structure(list(message=as.character(message), call = call), class=class)
}

conditionMessage <- function(c) UseMethod("conditionMessage")
conditionCall <- function(c) UseMethod("conditionCall")

conditionMessage.condition <- function(c) c$message
conditionCall.condition <- function(c) c$call

print.condition <- function(x, ...) {
    msg <- conditionMessage(x)
    call <- conditionCall(x)
    cl <- class(x)[1L]
    if (! is.null(call))
        cat("<", cl, " in ", deparse(call), ": ", msg, ">\n", sep="")
    else
        cat("<", cl, ": ", msg, ">\n", sep="")
    invisible(x)
}

as.character.condition <- function(x, ...) {
    msg <- conditionMessage(x)
    call <- conditionCall(x)
    cl <- class(x)[1L]
    if (! is.null(call))
        paste0(cl, " in ", deparse(call)[1L], ": ", msg, "\n")
    else
        paste0(cl, ": ", msg, "\n")
}

as.character.error <- function(x, ...) {
    msg <- conditionMessage(x)
    call <- conditionCall(x)
    if (! is.null(call))
        paste0("Error in ", deparse(call)[1L], ": ", msg, "\n")
    else
        paste0("Error: ", msg, "\n")
}

signalCondition <- function(cond) {
    if (! inherits(cond, "condition"))
        cond <- simpleCondition(cond)
    msg <- conditionMessage(cond)
    call <- conditionCall(cond)
    .Internal(.signalCondition(cond, msg, call))
}


##
##  Restarts
##

restartDescription <- function(r) r$description
restartFormals <- function(r) formals(r$handler)

print.restart <- function(x, ...) {
    cat(paste("<restart:", x[[1L]], ">\n"))
    invisible(x)
}

isRestart <- function(x) inherits(x, "restart")

findRestart <- function(name, cond = NULL) {
    i <- 1L
    repeat {
        r <- .Internal(.getRestart(i))
        if (is.null(r))
            return(NULL)
        else if (name == r[[1L]] &&
                 (is.null(cond) || is.null(r$test) || r$test(cond)))
            return(r)
        else i <- i + 1L
    }
}

computeRestarts <- function(cond = NULL) {
    val <- NULL
    i <- 1L
    repeat {
        r <- .Internal(.getRestart(i))
        if (is.null(r))
            return(val)
        else if (is.null(cond) || is.null(r$test) || r$test(cond))
            val <- c(val, list(r))
        i <- i + 1L
    }
}

invokeRestart <- function(r, ...) {
    if (! isRestart(r)) {
        res <- findRestart(r)
        if (is.null(res))
            stop(gettextf("no 'restart' '%s' found", as.character(r)),
                 domain = NA)
        r <- res
    }
    .Internal(.invokeRestart(r, list(...)))
}

invokeRestartInteractively <- function(r) {
    if (! interactive())
        stop("not an interactive session")
    if (! isRestart(r)) {
        res <- findRestart(r)
        if (is.null(res))
            stop(gettextf("no 'restart' '%s' found", as.character(r)),
                 domain = NA)
        r <- res
    }
    if (is.null(r$interactive)) {
        pars <- names(restartFormals(r))
        args <- NULL
        if (length(pars)) {
            cat("Enter values for restart arguments:\n\n")
            for (p in pars) {
            if (p == "...") {
		    prompt <- "... (a list): "
		    args <- c(args, eval(parse(prompt = prompt)))
		}
		else {
		    prompt <- paste0(p, ": ")
		    args <- c(args, list(eval(parse(prompt = prompt))))
		}
	    }
	}
    }
    else args <- r$interactive()
    .Internal(.invokeRestart(r, args))
}

withRestarts <- function(expr, ...) {
    docall <- function(fun, args) {
	if ((is.character(fun) && length(fun) == 1L) || is.name(fun))
	    fun <- get(as.character(fun), envir = parent.frame(),
                       mode = "function")
	do.call("fun", lapply(args, enquote))
    }
    makeRestart <- function(name = "",
			   handler = function(...) NULL,
			   description = "",
			   test = function(c) TRUE,
			   interactive = NULL) {
	structure(list(name = name, exit = NULL, handler = handler,
		       description = description, test = test,
		       interactive = interactive),
		  class = "restart")
    }
    makeRestartList <- function(...) {
        specs <- list(...)
        names <- names(specs)
        restarts <- vector("list", length(specs))
        for (i in seq_along(specs)) {
            spec <- specs[[i]]
            name <- names[i]
            if (is.function(spec))
                restarts[[i]] <- makeRestart(handler = spec)
            else if (is.character(spec))
                restarts[[i]] <- makeRestart(description = spec)
            else if (is.list(spec))
                restarts[[i]] <- docall("makeRestart", spec)
            else
               stop("not a valid restart specification")
            restarts[[i]]$name <- name
        }
        restarts
    }
    withOneRestart <- function(expr, restart) {
	doWithOneRestart <- function(expr, restart) {
	    restart$exit <- environment()
	    .Internal(.addRestart(restart))
	    expr
	}
	restartArgs <- doWithOneRestart(return(expr), restart)
	# The return in the call above will exit withOneRestart unless
	# the restart is invoked; we only get to this point if the restart
	# is invoked.  If we get here then the restart will have been
	# popped off the internal restart stack.
	docall(restart$handler, restartArgs)
    }
    withRestartList <- function(expr, restarts) {
	nr <- length(restarts)
	if (nr > 1L)
	    withOneRestart(withRestartList(expr, restarts[-nr]),
                           restarts[[nr]])
	else if (nr == 1L)
	    withOneRestart(expr, restarts[[1L]])
	else expr
    }
    restarts <- makeRestartList(...)
    if (length(restarts) == 0L)
        expr
    else if (length(restarts) == 1L)
        withOneRestart(expr, restarts[[1L]])
    else withRestartList(expr, restarts)
}


##
## Callbacks
##

.signalSimpleWarning <- function(msg, call)
    withRestarts({
           .Internal(.signalCondition(simpleWarning(msg, call), msg, call))
           .Internal(.dfltWarn(msg, call))
        }, muffleWarning = function() NULL)

.handleSimpleError <- function(h, msg, call)
    h(simpleError(msg, call))
#  File src/library/base/R/conflicts.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1998 B. D. Ripley
#  Copyright (C) 2005-2011 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

conflicts <- function(where = search(), detail = FALSE)
{
    if(length(where) < 1L) stop("argument 'where' of length 0")
    z <- vector(length(where), mode="list")
    names(z) <- where
    for(i in seq_along(where)) z[[i]] <- objects(pos = where[i])
    all <- unlist(z, use.names=FALSE)
    dups <- duplicated(all)
    dups <- all[dups]
    if(detail) {
	for(i in where) z[[i]] <- z[[i]][match(dups, z[[i]], 0L)]
	z[vapply(z, function(x) length(x) == 0L, NA)] <- NULL
	z
    } else dups
}
#  File src/library/base/R/connections.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

stdin <- function() .Internal(stdin())
stdout <- function() .Internal(stdout())
stderr <- function() .Internal(stderr())

isatty <- function(con) {
    if (!inherits(con, "terminal")) FALSE
    else .Internal(isatty(con))
}

readLines <- function(con = stdin(), n = -1L, ok = TRUE, warn = TRUE,
                      encoding = "unknown")
{
    if(is.character(con)) {
        con <- file(con, "r")
        on.exit(close(con))
    }
    .Internal(readLines(con, n, ok, warn, encoding))
}


writeLines <- function(text, con = stdout(), sep = "\n", useBytes = FALSE)
{
    if(is.character(con)) {
        con <- file(con, "w")
        on.exit(close(con))
    }
    .Internal(writeLines(text, con, sep, useBytes))
}

open <- function(con, ...)
    UseMethod("open")

open.connection <- function(con, open = "r", blocking = TRUE, ...)
    .Internal(open(con, open, blocking))

isOpen <- function(con, rw = "")
{
    rw <- pmatch(rw, c("read", "write"), 0L)
    .Internal(isOpen(con, rw))
}

isIncomplete <- function(con)
    .Internal(isIncomplete(con))

isSeekable <- function(con)
    .Internal(isSeekable(con))

close <- function(con, ...)
    UseMethod("close")

close.connection <- function (con, type = "rw", ...)
    .Internal(close(con, type))

flush <- function(con) UseMethod("flush")

flush.connection <- function (con)
    .Internal(flush(con))

file <- function(description = "", open = "", blocking = TRUE,
                 encoding = getOption("encoding"), raw = FALSE)
    .Internal(file(description, open, blocking, encoding, raw))

pipe <- function(description, open = "", encoding = getOption("encoding"))
    .Internal(pipe(description, open, encoding))

fifo <- function(description, open = "", blocking = FALSE,
                 encoding = getOption("encoding"))
    .Internal(fifo(description, open, blocking, encoding))

url <- function(description, open = "", blocking = TRUE,
                encoding = getOption("encoding"))
    .Internal(url(description, open, blocking, encoding))

gzfile <- function(description, open = "",
                   encoding = getOption("encoding"), compression = 6)
    .Internal(gzfile(description, open, encoding, compression))

unz <- function(description, filename, open = "",
                encoding = getOption("encoding"))
    .Internal(unz(paste(description, filename, sep=":"), open, encoding))

bzfile <- function(description, open = "", encoding = getOption("encoding"),
                   compression = 9)
    .Internal(bzfile(description, open, encoding, compression))

xzfile <- function(description, open = "", encoding = getOption("encoding"),
                   compression = 6)
    .Internal(xzfile(description, open, encoding, compression))

socketConnection <- function(host = "localhost", port, server = FALSE,
                             blocking = FALSE, open = "a+",
                             encoding = getOption("encoding"),
                             timeout = getOption("timeout"))
    .Internal(socketConnection(host, port, server, blocking, open, encoding,
                               timeout))

rawConnection <- function(object, open = "r") {
    .Internal(rawConnection(deparse(substitute(object)), object, open))
}

rawConnectionValue <- function(con) .Internal(rawConnectionValue(con))

textConnection <- function(object, open = "r", local = FALSE,
                           encoding = c("", "bytes", "UTF-8"))
{
    env <- if (local) parent.frame() else .GlobalEnv
    type <- match(match.arg(encoding), c("", "bytes", "UTF-8"))
    nm <- deparse(substitute(object))
    if(length(nm) != 1)
        stop("argument 'object' must deparse to a single character string")
    .Internal(textConnection(nm, object, open, env, type))
}

textConnectionValue <- function(con) .Internal(textConnectionValue(con))

seek <- function(con, ...)
    UseMethod("seek")

seek.connection <- function(con, where = NA, origin = "start", rw = "", ...)
{
    origin <- pmatch(origin, c("start", "current", "end"))
    rw <- pmatch(rw, c("read", "write"), 0L)
    if(is.na(origin))
        stop("'origin' must be one of 'start', 'current' or 'end'")
    .Internal(seek(con, as.double(where), origin, rw))
}

truncate <- function(con, ...)
    UseMethod("truncate")

truncate.connection <- function(con, ...)
{
    if(!isOpen(con)) stop("can only truncate an open connection")
    .Internal(truncate(con))
}

pushBack <- function(data, connection, newLine = TRUE)
    .Internal(pushBack(data, connection, newLine))

pushBackLength <- function(connection)
    .Internal(pushBackLength(connection))

clearPushBack <- function(connection)
    .Internal(clearPushBack(connection))

print.connection <- function(x, ...)
{
    print(unlist(summary(x)))
    invisible(x)
}

summary.connection <- function(object, ...)
    .Internal(summary.connection(object))

showConnections <- function(all = FALSE)
{
    set <- getAllConnections()
    if(!all) set <- set[set > 2L]
    ans <- matrix("", length(set), 7L)
    for(i in seq_along(set)) ans[i, ] <- unlist(summary.connection(set[i]))
    rownames(ans) <- set
    colnames(ans) <- c("description", "class", "mode", "text", "isopen",
                       "can read", "can write")
    if(!all) ans[ans[, 5L] == "opened", , drop = FALSE]
    else ans[, , drop = FALSE]
}

getAllConnections <- function()
    .Internal(getAllConnections())

getConnection <- function(what) .Internal(getConnection(what))

closeAllConnections <- function()
{
    # first re-divert any diversion of stderr.
    i <- sink.number(type = "message")
    if(i > 0L) sink(stderr(), type = "message")
    # now unwind the sink diversion stack.
    n <- sink.number()
    if(n > 0L) for(i in seq_len(n)) sink()
    # get all the open connections.
    set <- getAllConnections()
    set <- set[set > 2L]
    # and close all user connections.
    for(i in seq_along(set)) close(getConnection(set[i]))
    invisible()
}

readBin <- function(con, what, n = 1L, size = NA_integer_, signed = TRUE,
                    endian = .Platform$endian)
{
    if(is.character(con)) {
        con <- file(con, "rb")
        on.exit(close(con))
    }
    swap <- endian != .Platform$endian
    if(!is.character(what) || is.na(what) ||
       length(what) != 1L || ## hence length(what) == 1:
       !any(what == c("numeric", "double", "integer", "int", "logical",
	    "complex", "character", "raw")))
	what <- typeof(what)
    .Internal(readBin(con, what, n, size, signed, swap))
}

writeBin <-
    function(object, con, size = NA_integer_, endian = .Platform$endian,
             useBytes = FALSE)
{
    swap <- endian != .Platform$endian
    if(!is.vector(object) || mode(object) == "list")
        stop("can only write vector objects")
    if(is.character(con)) {
        con <- file(con, "wb")
        on.exit(close(con))
    }
    .Internal(writeBin(object, con, size, swap, useBytes))
}

readChar <- function(con, nchars, useBytes = FALSE)
{
    if(is.character(con)) {
        con <- file(con, "rb")
        on.exit(close(con))
    }
    .Internal(readChar(con, as.integer(nchars), useBytes))
}

writeChar <- function(object, con, nchars = nchar(object, type="chars"),
                      eos = "", useBytes = FALSE)
{
    if(!is.character(object))
        stop("can only write character objects")
    if(is.character(con)) {
        con <- file(con, "wb")
        on.exit(close(con))
    }
    .Internal(writeChar(object, con, as.integer(nchars), eos, useBytes))
}

gzcon <- function(con, level = 6, allowNonCompressed = TRUE)
    .Internal(gzcon(con, level, allowNonCompressed))

socketSelect <- function(socklist, write = FALSE, timeout = NULL) {
    if (is.null(timeout))
        timeout <- -1
    else if (timeout < 0)
        stop("'timeout' must be NULL or a non-negative number")
    if (length(write) < length(socklist))
        write <- rep_len(write, length(socklist))
    .Internal(sockSelect(socklist, write, timeout))
}

memCompress <-
    function(from, type = c("gzip", "bzip2", "xz", "none"))
{
    if(is.character(from))
        from <- charToRaw(paste(from, collapse = "\n"))
    else if(!is.raw(from)) stop("'from' must be raw or character")
    type <- match(match.arg(type), c("none", "gzip", "bzip2", "xz"))
    .Internal(memCompress(from, type))
}

memDecompress <-
    function(from,
             type = c("unknown", "gzip", "bzip2", "xz", "none"),
             asChar = FALSE)
{
    type <- match(match.arg(type),
                  c("none", "gzip", "bzip2", "xz", "unknown"))
    ans <- .Internal(memDecompress(from, type))
    if(asChar) rawToChar(ans) else ans
}
#  File src/library/base/R/constants.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

pi <- 4*atan(1)

letters <- c("a","b","c","d","e","f","g","h","i","j","k","l", "m",
	     "n","o","p","q","r","s","t","u","v","w","x","y","z")

LETTERS <- c("A","B","C","D","E","F","G","H","I","J","K","L", "M",
	     "N","O","P","Q","R","S","T","U","V","W","X","Y","Z")

month.name <-
    c("January", "February", "March", "April", "May", "June",
      "July", "August", "September", "October", "November", "December")

month.abb <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
	       "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
#  File src/library/base/R/contributors.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

contributors <- function()
{
    outFile <- tempfile()
    outConn <- file(outFile, open = "w")
    writeLines(paste0("R is a project which is attempting to provide a ",
                      "modern piece of\nstatistical software for the ",
                      "GNU suite of software.\n\n",
                      "The current R is the result of a collaborative ",
                      "effort with\ncontributions from all over the ",
                      "world.\n\n"), outConn)
    writeLines(readLines(file.path(R.home("doc"), "AUTHORS")), outConn)
    writeLines("", outConn)
    writeLines(readLines(file.path(R.home("doc"), "THANKS")), outConn)
    close(outConn)
    file.show(outFile, delete.file = TRUE)
}
#  File src/library/base/R/cut.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

cut <- function(x, ...) UseMethod("cut")

cut.default <-
    function (x, breaks, labels = NULL, include.lowest = FALSE,
              right = TRUE, dig.lab = 3L, ordered_result = FALSE, ...)
{
    if (!is.numeric(x)) stop("'x' must be numeric")
    if (length(breaks) == 1L) {
	if (is.na(breaks) || breaks < 2L)
	    stop("invalid number of intervals")
	nb <- as.integer(breaks + 1) # one more than #{intervals}
	dx <- diff(rx <- range(x, na.rm = TRUE))
	if(dx == 0) dx <- abs(rx[1L])
	breaks <- seq.int(rx[1L] - dx/1000,
                          rx[2L] + dx/1000, length.out = nb)
    } else nb <- length(breaks <- sort.int(as.double(breaks)))
    if (anyDuplicated(breaks)) stop("'breaks' are not unique")
    codes.only <- FALSE
    if (is.null(labels)) {#- try to construct nice ones ..
	for(dig in dig.lab:max(12L, dig.lab)) {
	    ch.br <- formatC(breaks, digits = dig, width = 1L)
	    if(ok <- all(ch.br[-1L] != ch.br[-nb])) break
	}
	labels <-
	    if(ok) paste0(if(right)"(" else "[",
                          ch.br[-nb], ",", ch.br[-1L],
                          if(right)"]" else ")")
	    else paste("Range", seq_len(nb - 1L), sep="_")
        if (ok && include.lowest) {
            if (right)
                substr(labels[1L], 1L, 1L) <- "[" # was "("
            else
                substring(labels[nb-1L],
                          nchar(labels[nb-1L], "c")) <- "]" # was ")"
        }
    } else if (is.logical(labels) && !labels)
        codes.only <- TRUE
    else if (length(labels) != nb - 1L)
        stop("lengths of 'breaks' and 'labels' differ")
    code <- .bincode(x, breaks, right, include.lowest)
    if(codes.only) code
    else factor(code, seq_along(labels), labels, ordered = ordered_result)
}

## called from image.default and for use in packages.
.bincode <- function(x, breaks, right = TRUE, include.lowest = FALSE)
    .Internal(bincode(x, breaks, right, include.lowest))

#  File src/library/base/R/data.matrix.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

data.matrix <- function(frame, rownames.force = NA)
{
    if(!is.data.frame(frame)) return(as.matrix(frame))

    d <- dim(frame)
    rn <- if(rownames.force %in% FALSE) NULL
    else if(rownames.force %in% TRUE) row.names(frame)
    else {if(.row_names_info(frame) <= 0L) NULL else row.names(frame)}

    for(i in seq_len(d[2L])) {
        xi <- frame[[i]]
        ## at present is.numeric suffices, but let's be cautious
        if(is.integer(xi) || is.numeric(xi)) next
        if(is.logical(xi) || is.factor(xi)) {
            frame[[i]] <- as.integer(xi)
            next
        }
        frame[[i]] <- if(isS4(xi)) methods::as(xi, "numeric") else as.numeric(xi)
    }

    ## it makes sense to find the type needed first.
    intOK <- all(unlist(lapply(frame, is.integer)))
    x <- matrix(if(intOK) NA_integer_ else NA_real_,
                nrow = d[1L], ncol = d[2L],
		dimnames = list(rn, names(frame)) )
    for(i in seq_len(d[2L])) x[, i] <- frame[[i]]
    x
}
#  File src/library/base/R/dataframe.R
#  Part of the R package, http://www.R-project.org
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

# Statlib code by John Chambers, Bell Labs, 1994
# Changes Copyright (C) 1998-2013 The R Core Team


## As from R 2.4.0, row.names can be either character or integer.
## row.names() will always return character.
## attr(, "row.names") will return either character or integer.
##
## Do not assume that the internal representation is either, since
## 1L:n is stored as the integer vector c(NA, n) to save space (and
## the C-level code to get/set the attribute makes the appropriate
## translations.
##
## As from 2.5.0 c(NA, n > 0) indicates deliberately assigned row names,
## and c(NA, n < 0) automatic row names.

## We cannot allow long vectors as elements until we can handle
## duplication of row names.

.row_names_info <- function(x, type = 1L)
    .Internal(shortRowNames(x, type))

row.names <- function(x) UseMethod("row.names")
row.names.data.frame <- function(x) as.character(attr(x, "row.names"))
row.names.default <- function(x) if(!is.null(dim(x))) rownames(x)# else NULL

.set_row_names <- function(n)
    if(n > 0) c(NA_integer_, -n) else integer()

`row.names<-` <- function(x, value) UseMethod("row.names<-")

`row.names<-.data.frame` <- function(x, value)
{
    if (!is.data.frame(x)) x <- as.data.frame(x)
    n <- .row_names_info(x, 2L)
    if(is.null(value)) { # set automatic row.names
        attr(x, "row.names") <- .set_row_names(n)
        return(x)
    }
    ## do this here, as e.g. POSIXlt changes length when coerced.
    if( is.object(value) || !is.integer(value) )
        value <- as.character(value)
    if(n == 0L) {
        ## we have to be careful here.  This could be a
        ## 0-row data frame or an invalid one being constructed.
        if(!is.null(attr(x, "row.names")) && length(value) > 0L)
           stop("invalid 'row.names' length")
    }
    else if (length(value) != n)
	stop("invalid 'row.names' length")
    if (anyDuplicated(value)) {
        nonuniq <- sort(unique(value[duplicated(value)]))
        warning(ngettext(length(nonuniq),
                         sprintf("non-unique value when setting 'row.names': %s",
                                 sQuote(nonuniq[1L])),
                         sprintf("non-unique values when setting 'row.names': %s",
                                 paste(sQuote(nonuniq), collapse = ", "))),
                domain = NA, call. = FALSE)
	stop("duplicate 'row.names' are not allowed")
    }
    if (any(is.na(value)))
	stop("missing values in 'row.names' are not allowed")
    attr(x, "row.names") <- value
    x
}

`row.names<-.default` <- function(x, value) `rownames<-`(x, value)

is.na.data.frame <- function (x)
{
    ## need to special-case no columns
    y <- if (length(x)) {
        do.call("cbind", lapply(x, "is.na")) # gives a matrix
    } else matrix(FALSE, length(row.names(x)), 0)
    if(.row_names_info(x) > 0L) rownames(y) <- row.names(x)
    y
}

is.data.frame <- function(x) inherits(x, "data.frame")

I <- function(x) { structure(x, class = unique(c("AsIs", oldClass(x)))) }

print.AsIs <- function (x, ...)
{
    cl <- oldClass(x)
    oldClass(x) <- cl[cl != "AsIs"]
    NextMethod("print")
    invisible(x)
}


t.data.frame <- function(x)
{
    x <- as.matrix(x)
    NextMethod("t")
}

dim.data.frame <- function(x) c(.row_names_info(x, 2L), length(x))

dimnames.data.frame <- function(x) list(row.names(x), names(x))

`dimnames<-.data.frame` <- function(x, value)
{
    d <- dim(x)
    if(!is.list(value) || length(value) != 2L)
	stop("invalid 'dimnames' given for data frame")
    ## do the coercion first, as might change length
    value[[1L]] <- as.character(value[[1L]])
    value[[2L]] <- as.character(value[[2L]])
    if(d[[1L]] != length(value[[1L]]) || d[[2L]] != length(value[[2L]]))
	stop("invalid 'dimnames' given for data frame")
    row.names(x) <- value[[1L]] # checks validity
    names(x) <- value[[2L]]
    x
}

as.data.frame <- function(x, row.names = NULL, optional = FALSE, ...)
{
    if(is.null(x))			# can't assign class to NULL
	return(as.data.frame(list()))
    UseMethod("as.data.frame")
}

as.data.frame.default <- function(x, ...)
    stop(gettextf("cannot coerce class \"%s\" to a data.frame",
                  deparse(class(x))),
         domain = NA)

###  Here are methods ensuring that the arguments to "data.frame"
###  are in a form suitable for combining into a data frame.

as.data.frame.data.frame <- function(x, row.names = NULL, ...)
{
    cl <- oldClass(x)
    i <- match("data.frame", cl)
    if(i > 1L)
	class(x) <- cl[ - (1L:(i-1L))]
    if(!is.null(row.names)){
        nr <- .row_names_info(x, 2L)
	if(length(row.names) == nr)
	    attr(x, "row.names") <- row.names
	else
            stop(sprintf(ngettext(nr,
                                  "invalid 'row.names', length %d for a data frame with %d row",
                                  "invalid 'row.names', length %d for a data frame with %d rows"),
                         length(row.names), nr), domain = NA)
    }
    x
}

## prior to 1.8.0 this coerced names - PR#3280
as.data.frame.list <-
    function(x, row.names = NULL, optional = FALSE, ...,
             stringsAsFactors = default.stringsAsFactors())
{
    ## need to protect names in x.
    cn <- names(x)
    m <- match(c("row.names", "check.rows", "check.names", "stringsAsFactors"),
               cn, 0L)
    if(any(m)) {
        cn[m] <- paste0("..adfl.", cn[m])
        names(x) <- cn
    }
    x <- eval(as.call(c(expression(data.frame), x, check.names = !optional,
                        stringsAsFactors = stringsAsFactors)))
    if(any(m)) names(x) <- sub("^\\.\\.adfl\\.", "", names(x))
    if(!is.null(row.names)) {
	# row.names <- as.character(row.names)
	if(length(row.names) != dim(x)[[1L]])
            stop(sprintf(ngettext(length(row.names),
                                  "supplied %d row name for %d rows",
                                  "supplied %d row names for %d rows"),
                         length(row.names), dim(x)[[1L]]), domain = NA)
	attr(x, "row.names") <- row.names
    }
    x
}

as.data.frame.vector <- function(x, row.names = NULL, optional = FALSE, ...,
                                 nm = paste(deparse(substitute(x),
                                 width.cutoff = 500L), collapse=" ")  )
{
    force(nm)
    nrows <- length(x)
    if(is.null(row.names)) {
	if (nrows == 0L)
	    row.names <- character()
	else if(length(row.names <- names(x)) == nrows &&
		!anyDuplicated(row.names)) {}
	else row.names <- .set_row_names(nrows)
    }
    if(!is.null(names(x))) names(x) <- NULL # remove names as from 2.0.0
    value <- list(x)
    if(!optional) names(value) <- nm
    attr(value, "row.names") <- row.names
    class(value) <- "data.frame"
    value
}

as.data.frame.ts <- function(x, ...)
{
    if(is.matrix(x))
	as.data.frame.matrix(x, ...)
    else
	as.data.frame.vector(x, ...)
}

as.data.frame.raw  <- as.data.frame.vector
as.data.frame.factor  <- as.data.frame.vector
as.data.frame.ordered <- as.data.frame.vector
as.data.frame.integer <- as.data.frame.vector
as.data.frame.numeric <- as.data.frame.vector
as.data.frame.complex <- as.data.frame.vector

default.stringsAsFactors <- function()
{
    val <- getOption("stringsAsFactors")
    if(is.null(val)) val <- TRUE
    if(!is.logical(val) || is.na(val) || length(val) != 1L)
        stop('options("stringsAsFactors") not set to TRUE or FALSE')
    val
}

## in case someone passes 'nm'
as.data.frame.character <-
    function(x, ..., stringsAsFactors = default.stringsAsFactors())
{
    nm <- deparse(substitute(x), width.cutoff=500L)
    if(stringsAsFactors) x <- factor(x)
    if(!"nm" %in% names(list(...)))
        as.data.frame.vector(x, ..., nm = nm)
    else as.data.frame.vector(x, ...)
}

as.data.frame.logical <- as.data.frame.vector

as.data.frame.matrix <- function(x, row.names = NULL, optional = FALSE, ...,
                                 stringsAsFactors = default.stringsAsFactors())
{
    d <- dim(x)
    nrows <- d[1L]; ir <- seq_len(nrows)
    ncols <- d[2L]; ic <- seq_len(ncols)
    dn <- dimnames(x)
    ## surely it cannot be right to override the supplied row.names?
    ## changed in 1.8.0
    if(is.null(row.names)) row.names <- dn[[1L]]
    collabs <- dn[[2L]]
    if(any(empty <- !nzchar(collabs)))
	collabs[empty] <- paste0("V", ic)[empty]
    value <- vector("list", ncols)
    if(mode(x) == "character" && stringsAsFactors) {
	for(i in ic)
	    value[[i]] <- as.factor(x[,i])
    } else {
	for(i in ic)
	    value[[i]] <- as.vector(x[,i])
    }
    ## Explicitly check for NULL in case nrows==0
    if(is.null(row.names) || length(row.names) != nrows)
	row.names <- .set_row_names(nrows)
    if(length(collabs) == ncols)
	names(value) <- collabs
    else if(!optional)
	names(value) <- paste0("V", ic)
    attr(value, "row.names") <- row.names
    class(value) <- "data.frame"
    value
}

as.data.frame.model.matrix <-
    function(x, row.names = NULL, optional = FALSE, ...)
{
    d <- dim(x)
    nrows <- d[1L]
    dn <- dimnames(x)
    row.names <- dn[[1L]]
    value <- list(x)
    if(!is.null(row.names)) {
	row.names <- as.character(row.names)
	if(length(row.names) != nrows)
            stop(sprintf(ngettext(length(row.names),
                                  "supplied %d row name for %d rows",
                                  "supplied %d row names for %d rows"),
                          length(row.names), nrows), domain = NA)
    }
    else row.names <- .set_row_names(nrows)
    if(!optional) names(value) <- deparse(substitute(x))[[1L]]
    attr(value, "row.names") <- row.names
    class(value) <- "data.frame"
    value
}

as.data.frame.array <- function(x, row.names = NULL, optional = FALSE, ...)
{
    d <- dim(x)
    if(length(d) == 1L) { ## same as as.data.frame.vector, but deparsed here
        value <- as.data.frame.vector(drop(x), row.names, optional, ...)
        if(!optional) names(value) <- deparse(substitute(x))[[1L]]
        value
    } else if (length(d) == 2L) {
        as.data.frame.matrix(x, row.names, optional, ...)
    } else {
        dn <- dimnames(x)
        dim(x) <- c(d[1L], prod(d[-1L]))
        if(!is.null(dn)) {
            if(length(dn[[1L]])) rownames(x) <- dn[[1L]]
            for(i in 2L:length(d))
                if(is.null(dn[[i]])) dn[[i]] <- seq_len(d[i])
            colnames(x) <- interaction(expand.grid(dn[-1L]))
        }
        as.data.frame.matrix(x, row.names, optional, ...)
    }
}

## Allow extraction method to have changed the underlying class,
## so re-assign the class based on the result.
`[.AsIs` <- function(x, i, ...) I(NextMethod("["))


as.data.frame.AsIs <- function(x, row.names = NULL, optional = FALSE, ...)
{
    ## why not remove class and NextMethod here?
    if(length(dim(x)) == 2L)
	as.data.frame.model.matrix(x, row.names, optional)
    else { # as.data.frame.vector without removing names
        nrows <- length(x)
        nm <- paste(deparse(substitute(x), width.cutoff=500L), collapse=" ")
        if(is.null(row.names)) {
            if (nrows == 0L)
                row.names <- character()
            else if(length(row.names <- names(x)) == nrows &&
                    !anyDuplicated(row.names)) {}
            else row.names <- .set_row_names(nrows)
        }
        value <- list(x)
        if(!optional) names(value) <- nm
        attr(value, "row.names") <- row.names
        class(value) <- "data.frame"
        value
    }

}

###  This is the real "data.frame".
###  It does everything by calling the methods presented above.

data.frame <-
    function(..., row.names = NULL, check.rows = FALSE, check.names = TRUE,
             stringsAsFactors = default.stringsAsFactors())
{
    data.row.names <-
	if(check.rows && is.null(row.names))
	    function(current, new, i) {
		if(is.character(current)) new <- as.character(new)
		if(is.character(new)) current <- as.character(current)
		if(anyDuplicated(new))
		    return(current)
		if(is.null(current))
		    return(new)
		if(all(current == new) || all(current == ""))
		    return(new)
		stop(gettextf("mismatch of row names in arguments of 'data.frame\', item %d", i), domain = NA)
	    }
	else function(current, new, i) {
	    if(is.null(current)) {
		if(anyDuplicated(new)) {
		    warning(gettextf("some row.names duplicated: %s --> row.names NOT used",
                                     paste(which(duplicated(new)), collapse=",")),
                            domain = NA)
		    current
		} else new
	    } else current
	}
    object <- as.list(substitute(list(...)))[-1L]
    mirn <- missing(row.names) # record before possibly changing
    mrn <- is.null(row.names) # missing or NULL
    x <- list(...)
    n <- length(x)
    if(n < 1L) {
        if(!mrn) {
            if(is.object(row.names) || !is.integer(row.names))
                row.names <- as.character(row.names)
            if(any(is.na(row.names)))
                stop("row names contain missing values")
            if(anyDuplicated(row.names))
                stop(gettextf("duplicate row.names: %s",
                              paste(unique(row.names[duplicated(row.names)]),
                                    collapse = ", ")),
                     domain = NA)
        } else row.names <- integer()
	return(structure(list(), names = character(),
                         row.names = row.names,
			 class = "data.frame"))
    }
    vnames <- names(x)
    if(length(vnames) != n)
	vnames <- character(n)
    no.vn <- !nzchar(vnames)
    vlist <- vnames <- as.list(vnames)
    nrows <- ncols <- integer(n)
    for(i in seq_len(n)) {
        ## do it this way until all as.data.frame methods have been updated
	xi <- if(is.character(x[[i]]) || is.list(x[[i]]))
                 as.data.frame(x[[i]], optional = TRUE,
                               stringsAsFactors = stringsAsFactors)
        else as.data.frame(x[[i]], optional = TRUE)

        nrows[i] <- .row_names_info(xi) # signed for now
	ncols[i] <- length(xi)
	namesi <- names(xi)
	if(ncols[i] > 1L) {
	    if(length(namesi) == 0L) namesi <- seq_len(ncols[i])
	    if(no.vn[i]) vnames[[i]] <- namesi
	    else vnames[[i]] <- paste(vnames[[i]], namesi, sep=".")
	}
	else {
            if(length(namesi)) vnames[[i]] <- namesi
            else if (no.vn[[i]]) {
                tmpname <- deparse(object[[i]])[1L]
                if( substr(tmpname, 1L, 2L) == "I(" ) {
                    ntmpn <- nchar(tmpname, "c")
                    if(substr(tmpname, ntmpn, ntmpn) == ")")
                        tmpname <- substr(tmpname, 3L, ntmpn - 1L)
                }
                vnames[[i]] <- tmpname
            }
        } # end of ncols[i] <= 1
	if(mirn && nrows[i] > 0L) {
            rowsi <- attr(xi, "row.names")
            ## Avoid all-blank names
            nc <- nchar(rowsi, allowNA = FALSE)
            nc <- nc[!is.na(nc)]
            if(length(nc) && any(nc))
                row.names <- data.row.names(row.names, rowsi, i)
        }
        nrows[i] <- abs(nrows[i])
	vlist[[i]] <- xi
    }
    nr <- max(nrows)
    for(i in seq_len(n)[nrows < nr]) {
	xi <- vlist[[i]]
	if(nrows[i] > 0L && (nr %% nrows[i] == 0L)) {
            ## make some attempt to recycle column i
            xi <- unclass(xi) # avoid data-frame methods
            fixed <- TRUE
            for(j in seq_along(xi)) {
                xi1 <- xi[[j]]
                if(is.vector(xi1) || is.factor(xi1))
                    xi[[j]] <- rep(xi1, length.out = nr)
		else if(is.character(xi1) && inherits(xi1, "AsIs"))
                    xi[[j]] <- structure(rep(xi1, length.out = nr),
                                         class = class(xi1))
                else if(inherits(xi1, "Date") || inherits(xi1, "POSIXct"))
                    xi[[j]] <- rep(xi1, length.out = nr)
                else {
                    fixed <- FALSE
                    break
                }
            }
            if (fixed) {
                vlist[[i]] <- xi
                next
            }
        }
        stop(gettextf("arguments imply differing number of rows: %s",
                      paste(unique(nrows), collapse = ", ")),
             domain = NA)
    }
    value <- unlist(vlist, recursive=FALSE, use.names=FALSE)
    ## unlist() drops i-th component if it has 0 columns
    vnames <- unlist(vnames[ncols > 0L])
    noname <- !nzchar(vnames)
    if(any(noname))
	vnames[noname] <- paste("Var", seq_along(vnames), sep = ".")[noname]
    if(check.names)
	vnames <- make.names(vnames, unique=TRUE)
    names(value) <- vnames
    if(!mrn) { # non-null row.names arg was supplied
        if(length(row.names) == 1L && nr != 1L) {  # one of the variables
            if(is.character(row.names))
                row.names <- match(row.names, vnames, 0L)
            if(length(row.names) != 1L ||
               row.names < 1L || row.names > length(vnames))
                stop("'row.names' should specify one of the variables")
            i <- row.names
            row.names <- value[[i]]
            value <- value[ - i]
        } else if ( !is.null(row.names) && length(row.names) != nr )
            stop("row names supplied are of the wrong length")
    } else if( !is.null(row.names) && length(row.names) != nr ) {
        warning("row names were found from a short variable and have been discarded")
        row.names <- NULL
    }
    if(is.null(row.names))
        row.names <- .set_row_names(nr) #seq_len(nr)
    else {
        if(is.object(row.names) || !is.integer(row.names))
            row.names <- as.character(row.names)
        if(any(is.na(row.names)))
            stop("row names contain missing values")
        if(anyDuplicated(row.names))
            stop(gettextf("duplicate row.names: %s",
                          paste(unique(row.names[duplicated(row.names)]),
                                collapse = ", ")),
                 domain = NA)
    }
    attr(value, "row.names") <- row.names
    attr(value, "class") <- "data.frame"
    value
}


###  Subsetting and mutation methods
###  These are a little less general than S

`[.data.frame` <-
    function(x, i, j, drop = if(missing(i)) TRUE else length(cols) == 1)
{
    mdrop <- missing(drop)
    Narg <- nargs() - !mdrop  # number of arg from x,i,j that were specified
    has.j <- !missing(j)
    if(!all(names(sys.call()) %in% c("", "drop"))
       && !isS4(x)) # at least don't warn for callNextMethod!
        warning("named arguments other than 'drop' are discouraged")

    if(Narg < 3L) {  # list-like indexing or matrix indexing
        if(!mdrop) warning("'drop' argument will be ignored")
	if(missing(i)) return(x)
	if(is.matrix(i))
	    return(as.matrix(x)[i])  # desperate measures
        ## zero-column data frames prior to 2.4.0 had no names.
        nm <- names(x); if(is.null(nm)) nm <- character()
        ## if we have NA names, character indexing should always fail
        ## (for positive index length)
        if(!is.character(i) && any(is.na(nm))) { # less efficient version
            names(nm) <- names(x) <- seq_along(x)
            y <- NextMethod("[")
            cols <- names(y)
            if(any(is.na(cols))) stop("undefined columns selected")
            cols <- names(y) <- nm[cols]
        } else {
            y <- NextMethod("[")
            cols <- names(y)
            if(!is.null(cols) && any(is.na(cols)))
                stop("undefined columns selected")
        }
        ## added in 1.8.0
        if(anyDuplicated(cols)) names(y) <- make.unique(cols)
        ## since we have not touched the rows, copy over the raw row.names
	return(structure(y, class = oldClass(x),
                         row.names = .row_names_info(x, 0L)))
    }

    if(missing(i)) { # df[, j] or df[ , ]
        ## not quite the same as the 1/2-arg case, as 'drop' is used.
        if(drop && !has.j && length(x) == 1L) return(.subset2(x, 1L))
        nm <- names(x); if(is.null(nm)) nm <- character()
        if(has.j && !is.character(j) && any(is.na(nm))) {
            ## less efficient version
            names(nm) <- names(x) <- seq_along(x)
            y <- .subset(x, j)
            cols <- names(y)
            if(any(is.na(cols))) stop("undefined columns selected")
            cols <- names(y) <- nm[cols]
        } else {
            y <- if(has.j) .subset(x, j) else x
            cols <- names(y)
            if(any(is.na(cols))) stop("undefined columns selected")
        }
        if(drop && length(y) == 1L) return(.subset2(y, 1L))
        if(anyDuplicated(cols)) names(y) <- make.unique(cols)
        nrow <- .row_names_info(x, 2L)
        if(drop && !mdrop && nrow == 1L)
            return(structure(y, class = NULL, row.names = NULL))
        else
            return(structure(y, class = oldClass(x),
                             row.names = .row_names_info(x, 0L)))
    }

    ### df[i, j] or df[i , ]
    ## rewritten for R 2.5.0 to avoid duplicating x.
    xx <- x
    cols <- names(xx)  # needed for computation of 'drop' arg
    ## make a shallow copy
    x <- vector("list", length(x))
    ## attributes(x) <- attributes(xx) expands row names
    x <- .Internal(copyDFattr(xx, x))
    oldClass(x) <- attr(x, "row.names") <- NULL

    if(has.j) { # df[i, j]
        nm <- names(x); if(is.null(nm)) nm <- character()
        if(!is.character(j) && any(is.na(nm)))
            names(nm) <- names(x) <- seq_along(x)
        x <- x[j]
        cols <- names(x)  # needed for 'drop'
        if(drop && length(x) == 1L) {
            ## for consistency with [, <length-1>]
            if(is.character(i)) {
                rows <- attr(xx, "row.names")
                i <- pmatch(i, rows, duplicates.ok = TRUE)
            }
            ## need to figure which col was selected:
            ## cannot use .subset2 directly as that may
            ## use recursive selection for a logical index.
            xj <- .subset2(.subset(xx, j), 1L)
            return(if(length(dim(xj)) != 2L) xj[i] else xj[i, , drop = FALSE])
        }
        if(any(is.na(cols))) stop("undefined columns selected")
        ## fix up names if we altered them.
        if(!is.null(names(nm))) cols <- names(x) <- nm[cols]
        ## sxx <- match(cols, names(xx)) fails with duplicate names
        nxx <- structure(seq_along(xx), names=names(xx))
        sxx <- match(nxx[j], seq_along(xx))
    } else sxx <- seq_along(x)

    rows <- NULL # placeholder: only create row names when needed
                 # as this can be expensive.
    if(is.character(i)) {
        rows <- attr(xx, "row.names")
        i <- pmatch(i, rows, duplicates.ok = TRUE)
    }
    for(j in seq_along(x)) {
        xj <- xx[[ sxx[j] ]]
        ## had drop = drop prior to 1.8.0
        x[[j]] <- if(length(dim(xj)) != 2L) xj[i] else xj[i, , drop = FALSE]
    }

    if(drop) {
	n <- length(x)
	if(n == 1L) return(x[[1L]]) # drops attributes
	if(n > 1L) {
	    xj <- x[[1L]]
	    nrow <- if(length(dim(xj)) == 2L) dim(xj)[1L] else length(xj)
            ## for consistency with S: don't drop (to a list)
            ## if only one row, unless explicitly asked for
            drop <- !mdrop && nrow == 1L
	} else drop <- FALSE ## for n == 0
    }

    if(!drop) { # not else as previous section might reset drop
        ## row names might have NAs.
        if(is.null(rows)) rows <- attr(xx, "row.names")
        rows <- rows[i]
	if((ina <- any(is.na(rows))) | (dup <- anyDuplicated(rows))) {
	    ## both will coerce integer 'rows' to character:
	    if (!dup && is.character(rows)) dup <- "NA" %in% rows
	    if(ina)
		rows[is.na(rows)] <- "NA"
	    if(dup)
		rows <- make.unique(as.character(rows))
	}
        ## new in 1.8.0  -- might have duplicate columns
	if(has.j && anyDuplicated(nm <- names(x)))
            names(x) <- make.unique(nm)
        if(is.null(rows)) rows <- attr(xx, "row.names")[i]
	attr(x, "row.names") <- rows
	oldClass(x) <- oldClass(xx)
    }
    x
}

`[[.data.frame` <- function(x, ..., exact=TRUE)
{
    ## use in-line functions to refer to the 1st and 2nd ... arguments
    ## explicitly. Also will check for wrong number or empty args
    na <- nargs() - !missing(exact)
    if(!all(names(sys.call()) %in% c("", "exact")))
        warning("named arguments other than 'exact' are discouraged")

    if(na < 3L)
	(function(x, i, exact)
	  if(is.matrix(i)) as.matrix(x)[[i]]
 	  else .subset2(x, i, exact=exact))(x, ..., exact=exact)
    else {
        col <- .subset2(x, ..2, exact=exact)
        i <- if(is.character(..1))
            pmatch(..1, row.names(x), duplicates.ok = TRUE)
        else ..1
        ## we do want to dispatch on methods for a column.
        ## .subset2(col, i, exact=exact)
        col[[i, exact = exact]]
    }
}

`[<-.data.frame` <- function(x, i, j, value)
{
    if(!all(names(sys.call()) %in% c("", "value")))
        warning("named arguments are discouraged")

    nA <- nargs() # 'value' is never missing, so 3 or 4.
    if(nA == 4L) { ## df[,] or df[i,] or df[, j] or df[i,j]
	has.i <- !missing(i)
	has.j <- !missing(j)
    }
    else if(nA == 3L) {
        ## this collects both df[] and df[ind]
        if (is.atomic(value) && !is.null(names(value)))
            names(value) <- NULL
        if(missing(i) && missing(j)) { # case df[]
            i <- j <- NULL
            has.i <- has.j <- FALSE
            ## added in 1.8.0
            if(is.null(value)) return(x[logical()])
        } else { # case df[ind]
            ## really ambiguous, but follow common use as if list
            ## except for two column numeric matrix or full-sized logical matrix
            if(is.numeric(i) && is.matrix(i) && ncol(i) == 2) {
                # Rewrite i as a logical index
                index <- rep.int(FALSE, prod(dim(x)))
                dim(index) <- dim(x)
                tryCatch(index[i] <- TRUE,
                         error = function(e) stop(conditionMessage(e), call.=FALSE))
                # Put values in the right order
                o <- order(i[,2], i[,1])
                N <- length(value)
                if (length(o) %% N != 0L)
                    warning("number of items to replace is not a multiple of replacement length")
                if (N < length(o))
                    value <- rep(value, length.out=length(o))
                value <- value[o]
                i <- index
            }
            if(is.logical(i) && is.matrix(i) && all(dim(i) == dim(x))) {
                nreplace <- sum(i, na.rm=TRUE)
                if(!nreplace) return(x) # nothing to replace
                ## allow replication of length(value) > 1 in 1.8.0
                N <- length(value)
                if(N > 1L && N < nreplace && (nreplace %% N) == 0L)
                    value <- rep(value, length.out = nreplace)
                if(N > 1L && (length(value) != nreplace))
                    stop("'value' is the wrong length")
                n <- 0L
                nv <- nrow(x)
                for(v in seq_len(dim(i)[2L])) {
                    thisvar <- i[, v, drop = TRUE]
                    nv <- sum(thisvar, na.rm = TRUE)
                    if(nv) {
                        if(is.matrix(x[[v]]))
                            x[[v]][thisvar, ] <- if(N > 1L) value[n+seq_len(nv)] else value
                        else
                            x[[v]][thisvar] <- if(N > 1L) value[n+seq_len(nv)] else value
                    }
                    n <- n+nv
                }
                return(x)
            }  # end of logical matrix
            if(is.matrix(i))
                stop("unsupported matrix index in replacement")
            j <- i
            i <- NULL
            has.i <- FALSE
            has.j <- TRUE
        }
    }
    else {
	stop("need 0, 1, or 2 subscripts")
    }
    ## no columns specified
    if(has.j && length(j) == 0L) return(x)

    cl <- oldClass(x)
    ## delete class: S3 idiom to avoid any special methods for [[, etc
    class(x) <- NULL
    new.cols <- NULL
    nvars <- length(x)
    nrows <- .row_names_info(x, 2L)
    if(has.i) { # df[i, ] or df[i, j]
        rows <- NULL  # indicator that it is not yet set
        if(any(is.na(i)))
            stop("missing values are not allowed in subscripted assignments of data frames")
	if(char.i <- is.character(i)) {
            rows <- attr(x, "row.names")
	    ii <- match(i, rows)
	    nextra <- sum(new.rows <- is.na(ii))
	    if(nextra > 0L) {
		ii[new.rows] <- seq.int(from = nrows + 1L, length.out = nextra)
		new.rows <- i[new.rows]
	    }
	    i <- ii
	}
	if(all(i >= 0L) && (nn <- max(i)) > nrows) {
	    ## expand
            if(is.null(rows)) rows <- attr(x, "row.names")
	    if(!char.i) {
		nrr <- (nrows + 1L):nn
		if(inherits(value, "data.frame") &&
		   (dim(value)[1L]) >= length(nrr)) {
		    new.rows <- attr(value, "row.names")[seq_along(nrr)]
		    repl <- duplicated(new.rows) | match(new.rows, rows, 0L)
		    if(any(repl)) new.rows[repl] <- nrr[repl]
		}
		else new.rows <- nrr
	    }
	    x <- xpdrows.data.frame(x, rows, new.rows)
	    rows <- attr(x, "row.names")
	    nrows <- length(rows)
	}
	iseq <- seq_len(nrows)[i]
	if(any(is.na(iseq))) stop("non-existent rows not allowed")
    }
    else iseq <- NULL

    if(has.j) {
        if(any(is.na(j)))
            stop("missing values are not allowed in subscripted assignments of data frames")
	if(is.character(j)) {
            if("" %in% j) stop("column name \"\" cannot match any column")
	    jj <- match(j, names(x))
	    nnew <- sum(is.na(jj))
	    if(nnew > 0L) {
		n <- is.na(jj)
		jj[n] <- nvars + seq_len(nnew)
		new.cols <- j[n]
	    }
	    jseq <- jj
	}
	else if(is.logical(j) || min(j) < 0L)
	    jseq <- seq_along(x)[j]
	else {
	    jseq <- j
	    if(max(jseq) > nvars) {
		new.cols <- paste0("V",
                                   seq.int(from = nvars + 1L, to = max(jseq)))
		if(length(new.cols)  != sum(jseq > nvars))
		    stop("new columns would leave holes after existing columns")
                ## try to use the names of a list `value'
                if(is.list(value) && !is.null(vnm <- names(value))) {
                    p <- length(jseq)
                    if(length(vnm) < p) vnm <- rep_len(vnm, p)
                    new.cols <- vnm[jseq > nvars]
                }
	    }
	}
    }
    else jseq <- seq_along(x)

    ## addition in 1.8.0
    if(anyDuplicated(jseq))
        stop("duplicate subscripts for columns")
    n <- length(iseq)
    if(n == 0L) n <- nrows
    p <- length(jseq)
    m <- length(value)
    if(!is.list(value)) {
        if(p == 1L) {
            N <- NROW(value)
            if(N > n)
                stop(sprintf(ngettext(N,
                                      "replacement has %d row, data has %d",
                                      "replacement has %d rows, data has %d"),
                             N, n), domain = NA)
            if(N < n && N > 0L)
                if(n %% N == 0L && length(dim(value)) <= 1L)
                    value <- rep(value, length.out = n)
                else
                    stop(sprintf(ngettext(N,
                                          "replacement has %d row, data has %d",
                                          "replacement has %d rows, data has %d"),
                                 N, nrows), domain = NA)
            if (!is.null(names(value))) names(value) <- NULL
            value <- list(value)
         } else {
            if(m < n*p && (m == 0L || (n*p) %% m))
                stop(sprintf(ngettext(m,
                                      "replacement has %d item, need %d",
                                      "replacement has %d items, need %d"),
                             m, n*p), domain = NA)
            value <- matrix(value, n, p)  ## will recycle
            value <- split(value, col(value))
        }
	dimv <- c(n, p)
    } else { # a list
        ## careful, as.data.frame turns things into factors.
	## value <- as.data.frame(value)
        value <- unclass(value) # to avoid data frame indexing
        lens <- vapply(value, NROW, 1L)
        for(k in seq_along(lens)) {
            N <- lens[k]
            if(n != N && length(dim(value[[k]])) == 2L)
                stop(sprintf(ngettext(N,
                                      "replacement element %d is a matrix/data frame of %d row, need %d",
                                      "replacement element %d is a matrix/data frame of %d rows, need %d"),
                             k, N, n),
                     domain = NA)
            if(N > 0L && N < n && n %% N)
                stop(sprintf(ngettext(N,
                                      "replacement element %d has %d row, need %d",
                                      "replacement element %d has %d rows, need %d"),
                             k, N, n), domain = NA)
            ## these fixing-ups will not work for matrices
            if(N > 0L && N < n) value[[k]] <- rep(value[[k]], length.out = n)
            if(N > n) {
                warning(sprintf(ngettext(N,
                                         "replacement element %d has %d row to replace %d rows",
                                         "replacement element %d has %d rows to replace %d rows"),
                                k, N, n), domain = NA)
                value[[k]] <- value[[k]][seq_len(n)]
            }
        }
	dimv <- c(n, length(value))
    }
    nrowv <- dimv[1L]
    if(nrowv < n && nrowv > 0L) {
	if(n %% nrowv == 0L)
	    value <- value[rep_len(seq_len(nrowv), n),,drop = FALSE]
	else
            stop(sprintf(ngettext(nrowv,
                                  "%d row in value to replace %d rows",
                                  "%d rows in value to replace %d rows"),
                         nrowv, n), domain = NA)
    }
    else if(nrowv > n)
        warning(sprintf(ngettext(nrowv,
                                 "replacement data has %d row to replace %d rows",
                                 "replacement data has %d rows to replace %d rows"),
                        nrowv, n), domain = NA)
    ncolv <- dimv[2L]
    jvseq <- seq_len(p)
    if(ncolv < p) jvseq <- rep_len(seq_len(ncolv), p)
    else if(ncolv > p) {
        warning(sprintf(ngettext(ncolv,
                                 "provided %d variable to replace %d variables",
                                 "provided %d variables to replace %d variables"),
                        ncolv, p), domain = NA)
        new.cols <- new.cols[seq_len(p)]
    }
    if(length(new.cols)) {
        ## extend and name now, as assignment of NULL may delete cols later.
        nm <- names(x)
        rows <- .row_names_info(x, 0L)
        a <- attributes(x); a["names"] <- NULL
        x <- c(x, vector("list", length(new.cols)))
        attributes(x) <- a
        names(x) <- c(nm, new.cols)
        attr(x, "row.names") <- rows
    }
    if(has.i)
	for(jjj in seq_len(p)) {
	    jj <- jseq[jjj]
	    vjj <- value[[ jvseq[[jjj]] ]]
            if(jj <= nvars) {
                ## if a column exists, preserve its attributes
                if(length(dim(x[[jj]])) != 2L) x[[jj]][iseq] <- vjj
                else x[[jj]][iseq, ] <- vjj
            } else {
                ## try to make a new column match in length: may be an error
                x[[jj]] <- vjj[FALSE]
                if(length(dim(vjj)) == 2L) {
                    length(x[[j]]) <- nrows * ncol(vjj)
                    dim(x[[j]]) <- c(nrows, ncol(vjj))
                    x[[jj]][iseq, ] <- vjj
                } else {
                    length(x[[j]]) <- nrows
                    x[[jj]][iseq] <- vjj
                }
            }
	}
    else if(p > 0L)
      for(jjj in p:1L) { # we might delete columns with NULL
        ## ... and for that reason, we'd better ensure that jseq is increasing!
        o <- order(jseq)
        jseq <- jseq[o]
        jvseq <- jvseq[o]

        jj <- jseq[jjj]
        v <- value[[ jvseq[[jjj]] ]]
        ## This is consistent with the have.i case rather than with
        ## [[<- and $<- (which throw an error).  But both are plausible.
        if (nrows > 0L && !length(v)) length(v) <- nrows
	x[[jj]] <- v
        if (!is.null(v) && is.atomic(x[[jj]]) && !is.null(names(x[[jj]])))
            names(x[[jj]]) <- NULL
    }
    if(length(new.cols) > 0L) {
        new.cols <- names(x) # we might delete columns with NULL
        ## added in 1.8.0
        if(anyDuplicated(new.cols)) names(x) <- make.unique(new.cols)
    }
    class(x) <- cl
    x
}

`[[<-.data.frame` <- function(x, i, j, value)
{
    if(!all(names(sys.call()) %in% c("", "value")))
        warning("named arguments are discouraged")

    cl <- oldClass(x)
    ## delete class: Version 3 idiom
    ## to avoid any special methods for [[<-
    class(x) <- NULL
    nrows <- .row_names_info(x, 2L)
    if(is.atomic(value) && !is.null(names(value))) names(value) <- NULL
    if(nargs() < 4L) {
	## really ambiguous, but follow common use as if list
        nc <- length(x)
	if(!is.null(value)) {
            N <- NROW(value)
            if(N > nrows)
                stop(sprintf(ngettext(N,
                                      "replacement has %d row, data has %d",
                                      "replacement has %d rows, data has %d"),
                             N, nrows), domain = NA)
            if(N < nrows)
                if(N > 0L && (nrows %% N == 0L) && length(dim(value)) <= 1L)
                    value <- rep(value, length.out = nrows)
                else
                    stop(sprintf(ngettext(N,
                                          "replacement has %d row, data has %d",
                                          "replacement has %d rows, data has %d"),
                                 N, nrows), domain = NA)
	}
	x[[i]] <- value
        ## added in 1.8.0 -- make sure there is a name
        if(length(x) > nc) {
            nc <- length(x)
            if(names(x)[nc] == "") names(x)[nc] <- paste0("V", nc)
            names(x) <- make.unique(names(x))
        }
	class(x) <- cl
	return(x)
    }
    if(missing(i) || missing(j))
	stop("only valid calls are x[[j]] <- value or x[[i,j]] <- value")
    rows <- attr(x, "row.names")
    nvars <- length(x)
    if(n <- is.character(i)) {
	ii <- match(i, rows)
	n <- sum(new.rows <- is.na(ii))
	if(n > 0L) {
	    ii[new.rows] <- seq.int(from = nrows + 1L, length.out = n)
	    new.rows <- i[new.rows]
	}
	i <- ii
    }
    if(all(i >= 0L) && (nn <- max(i)) > nrows) {
	## expand
	if(n == 0L) {
	    nrr <- (nrows + 1L):nn
	    if(inherits(value, "data.frame") &&
	       (dim(value)[1L]) >= length(nrr)) {
		new.rows <- attr(value, "row.names")[seq_len(nrr)]
		repl <- duplicated(new.rows) | match(new.rows, rows, 0L)
		if(any(repl)) new.rows[repl] <- nrr[repl]
	    }
	    else new.rows <- nrr
	}
	x <- xpdrows.data.frame(x, rows, new.rows)
	rows <- attr(x, "row.names")
	nrows <- length(rows)
    }

    ## FIXME: this is wasteful and probably unnecessary
    iseq <- seq_len(nrows)[i]
    if(any(is.na(iseq)))
	stop("non-existent rows not allowed")

    if(is.character(j)) {
        if("" %in% j) stop("column name \"\" cannot match any column")
	jseq <- match(j, names(x))
	if(any(is.na(jseq)))
            stop(gettextf("replacing element in non-existent column: %s",
                          j[is.na(jseq)]), domain = NA)
    }
    else if(is.logical(j) || min(j) < 0L)
	jseq <- seq_along(x)[j]
    else {
	jseq <- j
	if(max(jseq) > nvars)
            stop(gettextf("replacing element in non-existent column: %s",
                          jseq[jseq > nvars]), domain = NA)
    }
    if(length(iseq) > 1L || length(jseq) > 1L)
	stop("only a single element should be replaced")
    x[[jseq]][[iseq]] <- value
    class(x) <- cl
    x
}

## added in 1.8.0
`$<-.data.frame` <- function(x, name, value)
{
    cl <- oldClass(x)
    ## delete class: Version 3 idiom
    ## to avoid any special methods for [[<-
    ## This forces a copy, but we are going to need one anyway
    ## and NAMED=1 prevents any further copying.
    class(x) <- NULL
    nrows <- .row_names_info(x, 2L)
    if(!is.null(value)) {
        N <- NROW(value)
        if(N > nrows)
            stop(sprintf(ngettext(N,
                                  "replacement has %d row, data has %d",
                                  "replacement has %d rows, data has %d"),
                         N, nrows), domain = NA)
        if (N < nrows)
            if (N > 0L && (nrows %% N == 0L) && length(dim(value)) <= 1L)
                value <- rep(value, length.out = nrows)
            else
                stop(sprintf(ngettext(N,
                                      "replacement has %d row, data has %d",
                                      "replacement has %d rows, data has %d"),
                             N, nrows), domain = NA)
        if(is.atomic(value) && !is.null(names(value))) names(value) <- NULL
    }
    x[[name]] <- value
    class(x) <- cl
    return(x)
}

xpdrows.data.frame <- function(x, old.rows, new.rows)
{
    nc <- length(x)
    nro <- length(old.rows)
    nrn <- length(new.rows)
    nr <- nro + nrn
    for (i in seq_len(nc)) {
	y <- x[[i]]
	dy <- dim(y)
	cy <- oldClass(y)
	class(y) <- NULL
	if (length(dy) == 2L) {
	    dny <- dimnames(y)
	    if (length(dny[[1L]]) > 0L)
		dny[[1L]] <- c(dny[[1L]], new.rows)
	    z <- array(y[1L], dim = c(nr, nc), dimnames = dny)
	    z[seq_len(nro), ] <- y
	    class(z) <- cy
	    x[[i]] <- z
	}
	else {
	    ay <- attributes(y)
	    if (length(names(y)) > 0L)
		ay$names <- c(ay$names, new.rows)
	    length(y) <- nr
	    attributes(y) <- ay
	    class(y) <- cy
	    x[[i]] <- y
	}
    }
    attr(x, "row.names") <- c(old.rows, new.rows)
    x
}


### Here are the methods for rbind and cbind.

cbind.data.frame <- function(..., deparse.level = 1)
    data.frame(..., check.names = FALSE)

rbind.data.frame <- function(..., deparse.level = 1)
{
    match.names <- function(clabs, nmi)
    {
	if(identical(clabs, nmi)) NULL
	else if(length(nmi) == length(clabs) && all(nmi %in% clabs)) {
            ## we need 1-1 matches here
	    m <- pmatch(nmi, clabs, 0L)
            if(any(m == 0L))
                stop("names do not match previous names")
            m
	} else stop("names do not match previous names")
    }
    Make.row.names <- function(nmi, ri, ni, nrow)
    {
	if(nzchar(nmi)) {
            if(ni == 0L) character()  # PR8506
	    else if(ni > 1L) paste(nmi, ri, sep = ".")
	    else nmi
	}
	else if(nrow > 0L && identical(ri, seq_len(ni)))
	    as.integer(seq.int(from = nrow + 1L, length.out = ni))
	else ri
    }
    allargs <- list(...)
    allargs <- allargs[vapply(allargs, length, 1L) > 0L]
    if(length(allargs)) {
    ## drop any zero-row data frames, as they may not have proper column
    ## types (e.g. NULL).
        nr <- vapply(allargs, function(x)
                     if(is.data.frame(x)) .row_names_info(x, 2L)
                     else if(is.list(x)) length(x[[1L]]) # mismatched lists are checked later
                     else length(x), 1L)
        if(any(nr > 0L)) allargs <- allargs[nr > 0L]
        else return(allargs[[1L]]) # pretty arbitrary
    }
    n <- length(allargs)
    if(n == 0L)
	return(structure(list(),
			 class = "data.frame",
			 row.names = integer()))
    nms <- names(allargs)
    if(is.null(nms))
	nms <- character(n)
    cl <- NULL
    perm <- rows <- rlabs <- vector("list", n)
    nrow <- 0L
    value <- clabs <- NULL
    all.levs <- list()
    for(i in seq_len(n)) {
	## check the arguments, develop row and column labels
	xi <- allargs[[i]]
	nmi <- nms[i]
        ## coerce matrix to data frame
        if(is.matrix(xi)) allargs[[i]] <- xi <- as.data.frame(xi)
	if(inherits(xi, "data.frame")) {
	    if(is.null(cl))
		cl <- oldClass(xi)
	    ri <- attr(xi, "row.names")
	    ni <- length(ri)
	    if(is.null(clabs))
		clabs <- names(xi)
	    else {
                if(length(xi) != length(clabs))
                    stop("numbers of columns of arguments do not match")
		pi <- match.names(clabs, names(xi))
		if( !is.null(pi) ) perm[[i]] <- pi
	    }
	    rows[[i]] <- seq.int(from = nrow + 1L, length.out = ni)
	    rlabs[[i]] <- Make.row.names(nmi, ri, ni, nrow)
	    nrow <- nrow + ni
	    if(is.null(value)) {
		value <- unclass(xi)
		nvar <- length(value)
		all.levs <- vector("list", nvar)
		has.dim <- logical(nvar)
                facCol <- logical(nvar)
                ordCol <- logical(nvar)
		for(j in seq_len(nvar)) {
		    xj <- value[[j]]
		    if( !is.null(levels(xj)) ) {
			all.levs[[j]] <- levels(xj)
                        facCol[j] <- TRUE # turn categories into factors
                    } else facCol[j] <- is.factor(xj)
                    ordCol[j] <- is.ordered(xj)
		    has.dim[j] <- length(dim(xj)) == 2L
		}
	    }
	    else for(j in seq_len(nvar)) {
                xij <- xi[[j]]
                if(is.null(pi) || is.na(jj <- pi[[j]])) jj <- j
                if(facCol[jj]) {
                    if(length(lij <- levels(xij))) {
                        all.levs[[jj]] <- unique(c(all.levs[[jj]], lij))
                        ordCol[jj] <- ordCol[jj] & is.ordered(xij)
                    } else if(is.character(xij))
                        all.levs[[jj]] <- unique(c(all.levs[[jj]], xij))
                }
            }
	}
	else if(is.list(xi)) {
	    ni <- range(vapply(xi, length, 1L))
	    if(ni[1L] == ni[2L])
		ni <- ni[1L]
	    else stop("invalid list argument: all variables should have the same length")
	    rows[[i]] <- ri <-
                as.integer(seq.int(from = nrow + 1L, length.out = ni))
	    nrow <- nrow + ni
	    rlabs[[i]] <- Make.row.names(nmi, ri, ni, nrow)
	    if(length(nmi <- names(xi)) > 0L) {
		if(is.null(clabs))
		    clabs <- nmi
		else {
                    if(length(xi) != length(clabs))
                        stop("numbers of columns of arguments do not match")
		    pi <- match.names(clabs, nmi)
		    if( !is.null(pi) ) perm[[i]] <- pi
		}
	    }
	}
	else if(length(xi)) {
	    rows[[i]] <- nrow <- nrow + 1L
	    rlabs[[i]] <- if(nzchar(nmi)) nmi else as.integer(nrow)
	}
    }
    nvar <- length(clabs)
    if(nvar == 0L)
	nvar <- max(vapply(allargs, length, 1L)) # only vector args
    if(nvar == 0L)
	return(structure(list(), class = "data.frame",
			 row.names = integer()))
    pseq <- seq_len(nvar)
    if(is.null(value)) { # this happens if there has been no data frame
	value <- list()
	value[pseq] <- list(logical(nrow)) # OK for coercion except to raw.
        all.levs <- vector("list", nvar)
        has.dim <- logical(nvar)
        facCol <- logical(nvar)
        ordCol <- logical(nvar)
    }
    names(value) <- clabs
    for(j in pseq)
	if(length(lij <- all.levs[[j]]))
            value[[j]] <-
                factor(as.vector(value[[j]]), lij, ordered = ordCol[j])
    if(any(has.dim)) {
	rmax <- max(unlist(rows))
	for(i in pseq[has.dim])
	    if(!inherits(xi <- value[[i]], "data.frame")) {
		dn <- dimnames(xi)
		rn <- dn[[1L]]
		if(length(rn) > 0L) length(rn) <- rmax
		pi <- dim(xi)[2L]
		length(xi) <- rmax * pi
		value[[i]] <- array(xi, c(rmax, pi), list(rn, dn[[2L]]))
	    }
    }
    for(i in seq_len(n)) {
	xi <- unclass(allargs[[i]])
	if(!is.list(xi))
	    if(length(xi) != nvar)
		xi <- rep(xi, length.out = nvar)
	ri <- rows[[i]]
	pi <- perm[[i]]
	if(is.null(pi)) pi <- pseq
	for(j in pseq) {
	    jj <- pi[j]
            xij <- xi[[j]]
	    if(has.dim[jj]) {
		value[[jj]][ri,	 ] <- xij
                ## copy rownames
                rownames(value[[jj]])[ri] <- rownames(xij)
	    } else {
                ## coerce factors to vectors, in case lhs is character or
                ## level set has changed
                value[[jj]][ri] <- if(is.factor(xij)) as.vector(xij) else xij
                ## copy names if any
                if(!is.null(nm <- names(xij))) names(value[[jj]])[ri] <- nm
            }
	}
    }
    rlabs <- unlist(rlabs)
    if(anyDuplicated(rlabs))
        rlabs <- make.unique(as.character(unlist(rlabs)), sep = "")
    if(is.null(cl)) {
	as.data.frame(value, row.names = rlabs)
    } else {
	class(value) <- cl
	attr(value, "row.names") <- rlabs
	value
    }
}


### coercion and print methods

print.data.frame <-
    function(x, ..., digits = NULL, quote = FALSE, right = TRUE,
	     row.names = TRUE)
{
    n <- length(row.names(x))
    if(length(x) == 0L) {
        cat(gettextf("data frame with 0 columns and %d rows\n", n))
    } else if(n == 0L) {
        ## FIXME: header format is inconsistent here
	print.default(names(x), quote = FALSE)
	cat(gettext("<0 rows> (or 0-length row.names)\n"))
    } else {
	## format.<*>() : avoiding picking up e.g. format.AsIs
	m <- as.matrix(format.data.frame(x, digits = digits, na.encode = FALSE))
	if(!isTRUE(row.names))
	    dimnames(m)[[1L]] <- if(identical(row.names,FALSE))
		rep.int("", n) else row.names
	print(m, ..., quote = quote, right = right)
    }
    invisible(x)
}

as.matrix.data.frame <- function (x, rownames.force = NA, ...)
{
    dm <- dim(x)
    rn <- if(rownames.force %in% FALSE) NULL
    else if(rownames.force %in% TRUE) row.names(x)
    else {if(.row_names_info(x) <= 0L) NULL else row.names(x)}
    dn <- list(rn, names(x))
    if(any(dm == 0L))
	return(array(NA, dim = dm, dimnames = dn))
    p <- dm[2L]
    pseq <- seq_len(p)
    n <- dm[1L]
    X <- x # will contain the result;
    ## the "big question" is if we return a numeric or a character matrix
    class(X) <- NULL
    non.numeric <- non.atomic <- FALSE
    all.logical <- TRUE
    for (j in pseq) {
        if(inherits(X[[j]], "data.frame") && ncol(xj) > 1L)
            X[[j]] <- as.matrix(X[[j]])
        xj <- X[[j]]
        j.logic <- is.logical(xj)
        if(all.logical && !j.logic) all.logical <- FALSE
	if(length(levels(xj)) > 0L || !(j.logic || is.numeric(xj) || is.complex(xj))
	   || (!is.null(cl <- attr(xj, "class")) && # numeric classed objects to format:
	       any(cl %in% c("Date", "POSIXct", "POSIXlt"))))
	    non.numeric <- TRUE
	if(!is.atomic(xj))
	    non.atomic <- TRUE
    }
    if(non.atomic) {
	for (j in pseq) {
	    xj <- X[[j]]
	    if(!is.recursive(xj))
		X[[j]] <- as.list(as.vector(xj))
	}
    } else if(all.logical) {
        ## do nothing for logical columns if a logical matrix will result.
    } else if(non.numeric) {
	for (j in pseq) {
	    if (is.character(X[[j]]))
		next
	    xj <- X[[j]]
            miss <- is.na(xj)
	    xj <- if(length(levels(xj))) as.vector(xj) else format(xj)
            is.na(xj) <- miss
            X[[j]] <- xj
	}
    }
    ## These coercions could have changed the number of columns
    ## (e.g. class "Surv" coerced to character),
    ## so only now can we compute collabs.
    collabs <- as.list(dn[[2L]])
    for (j in pseq) {
        xj <- X[[j]]
        dj <- dim(xj)
	if(length(dj) == 2L && dj[2L] > 1L) { # matrix with >=2 col
	    dnj <- colnames(xj)
	    collabs[[j]] <- paste(collabs[[j]],
				  if(length(dnj)) dnj else seq_len(dj[2L]),
				  sep = ".")
	}
     }
    X <- unlist(X, recursive = FALSE, use.names = FALSE)
    dim(X) <- c(n, length(X)/n)
    dimnames(X) <- list(dn[[1L]], unlist(collabs, use.names = FALSE))
    ##NO! don't copy buggy S-plus!  either all matrices have class or none!!
    ##NO class(X) <- "matrix"
    X
}

Math.data.frame <- function (x, ...)
{
    mode.ok <- vapply(x, function(x) is.numeric(x) || is.complex(x), NA)
    if (all(mode.ok)) {
	x[] <- lapply(X = x, FUN = .Generic, ...)
	return(x)
    } else {
	vnames <- names(x)
	if (is.null(vnames)) vnames <- seq_along(x)
	stop("non-numeric variable in data frame: ", vnames[!mode.ok])
    }
}

Ops.data.frame <- function(e1, e2 = NULL)
{
    isList <- function(x) !is.null(x) && is.list(x)
    unary <- nargs() == 1L
    lclass <- nzchar(.Method[1L])
    rclass <- !unary && (nzchar(.Method[2L]))
    value <- list()
    rn <- NULL
    ## set up call as op(left, right)
    FUN <- get(.Generic, envir = parent.frame(), mode = "function")
    f <- if (unary) quote(FUN(left)) else quote(FUN(left, right))
    lscalar <- rscalar <- FALSE
    if(lclass && rclass) {
        nr <- .row_names_info(e1, 2L)
	if(.row_names_info(e1) > 0L) rn <- attr(e1, "row.names")
	cn <- names(e1)
	if(any(dim(e2) != dim(e1)))
	    stop(.Generic, " only defined for equally-sized data frames")
    } else if(lclass) {
	## e2 is not a data frame, but e1 is.
        nr <- .row_names_info(e1, 2L)
	if(.row_names_info(e1) > 0L) rn <- attr(e1, "row.names")
	cn <- names(e1)
	rscalar <- length(e2) <= 1L # e2 might be null
	if(isList(e2)) {
	    if(rscalar) e2 <- e2[[1L]]
	    else if(length(e2) != ncol(e1))
		stop(gettextf("list of length %d not meaningful", length(e2)),
                     domain = NA)
	} else {
	    if(!rscalar)
		e2 <- split(rep_len(as.vector(e2), prod(dim(e1))),
			    rep.int(seq_len(ncol(e1)),
                                    rep.int(nrow(e1), ncol(e1))))
	}
    } else {
	## e1 is not a data frame, but e2 is.
        nr <- .row_names_info(e2, 2L)
	if(.row_names_info(e2) > 0L) rn <- attr(e2, "row.names")
	cn <- names(e2)
	lscalar <- length(e1) <= 1L
	if(isList(e1)) {
	    if(lscalar) e1 <- e1[[1L]]
	    else if(length(e1) != ncol(e2))
		stop(gettextf("list of length %d not meaningful", length(e1)),
                     domain = NA)
	} else {
	    if(!lscalar)
		e1 <- split(rep_len(as.vector(e1), prod(dim(e2))),
			    rep.int(seq_len(ncol(e2)),
                                    rep.int(nrow(e2), ncol(e2))))
	}
    }
    for(j in seq_along(cn)) {
	left <- if(!lscalar) e1[[j]] else e1
	right <- if(!rscalar) e2[[j]] else e2
	value[[j]] <- eval(f)
    }
    if(.Generic %in% c("+","-","*","/","%%","%/%") ) {
	names(value) <- cn
	data.frame(value, row.names = rn, check.names = FALSE,
                   check.rows = FALSE)
    }
    else matrix(unlist(value, recursive = FALSE, use.names = FALSE),
		nrow = nr, dimnames = list(rn,cn))
}

Summary.data.frame <- function(..., na.rm)
{
    args <- list(...)
    args <- lapply(args, function(x) {
        x <- as.matrix(x)
        if(!is.numeric(x) && !is.complex(x))
            stop("only defined on a data frame with all numeric variables")
        x
    })
    do.call(.Generic, c(args, na.rm=na.rm))
}
#  File src/library/base/R/dates.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

## First shot at adding a "Date" class to base R.
## Representation is the number of whole days since 1970-01-01.

## The difftime class already covers time differences in days.

## Need to take timezone into account here
Sys.Date <- function() as.Date(as.POSIXlt(Sys.time()))

as.Date <- function(x, ...) UseMethod("as.Date")

as.Date.POSIXct <- function(x, tz = "UTC", ...)
{
    if(tz == "UTC") {
        z <- floor(unclass(x)/86400)
        attr(z, "tzone") <- NULL
        structure(z, class = "Date")
    } else
        as.Date(as.POSIXlt(x, tz = tz))
}

as.Date.POSIXlt <- function(x, ...) .Internal(POSIXlt2Date(x))

as.Date.factor <- function(x, ...) as.Date(as.character(x), ...)


as.Date.character <- function(x, format="", ...)
{
    charToDate <- function(x) {
	xx <- x[1L]
        if(is.na(xx)) {
            j <- 1L
            while(is.na(xx) && (j <- j+1L) <= length(x)) xx <- x[j]
            if(is.na(xx)) f <- "%Y-%m-%d" # all NAs
        }
	if(is.na(xx) ||
	   !is.na(strptime(xx, f <- "%Y-%m-%d", tz="GMT")) ||
	   !is.na(strptime(xx, f <- "%Y/%m/%d", tz="GMT"))
           ) return(strptime(x, f))
	stop("character string is not in a standard unambiguous format")
    }
    res <- if(missing(format)) charToDate(x) else strptime(x, format, tz="GMT")
    as.Date(res)
}

as.Date.numeric <- function(x, origin, ...)
{
    if(missing(origin)) stop("'origin' must be supplied")
    as.Date(origin, ...) + x
}

as.Date.default <- function(x, ...)
{
    if(inherits(x, "Date")) return(x)
    if(is.logical(x) && all(is.na(x)))
        return(structure(as.numeric(x), class = "Date"))
    stop(gettextf("do not know how to convert '%s' to class %s",
                  deparse(substitute(x)),
                  dQuote("Date")),
         domain = NA)
}

## convert from package date
as.Date.date <- function(x, ...)
{
    if(inherits(x, "date")) {
        x <- (x - 3653) # origin 1960-01-01
        return(structure(x, class = "Date"))
    } else stop(gettextf("'%s' is not a \"date\" object",
                         deparse(substitute(x)) ))
}

## convert from package chron
as.Date.dates <- function(x, ...)
{
    if(inherits(x, "dates")) {
        z <- attr(x, "origin")
        x <- trunc(as.numeric(x))
        if(length(z) == 3L && is.numeric(z))
            x  <- x + as.numeric(as.Date(paste(z[3L], z[1L], z[2L], sep="/")))
        return(structure(x, class = "Date"))
    } else stop(gettextf("'%s' is not a \"dates\" object",
                         deparse(substitute(x)) ))
}

format.Date <- function(x, ...)
{
    xx <- format(as.POSIXlt(x), ...)
    names(xx) <- names(x)
    xx
}

## could handle arrays for max.print
print.Date <- function(x, max = NULL, ...)
{
    if(is.null(max)) max <- getOption("max.print", 9999L)
    if(max < length(x)) {
	print(format(x[seq_len(max)]), max=max, ...)
	cat(' [ reached getOption("max.print") -- omitted',
	    length(x) - max, 'entries ]\n')
    } else print(format(x), max=max, ...)
    invisible(x)
}

summary.Date <- function(object, digits = 12L, ...)
{
    x <- summary.default(unclass(object), digits = digits, ...)
    if(m <- match("NA's", names(x), 0)) {
        NAs <- as.integer(x[m])
        x <- x[-m]
        attr(x, "NAs") <- NAs
    }
    class(x) <- c("summaryDefault", "table", oldClass(object))
    x
}

`+.Date` <- function(e1, e2)
{
    ## need to drop "units" attribute here
    coerceTimeUnit <- function(x)
        as.vector(round(switch(attr(x,"units"),
                               secs = x/86400, mins = x/1440, hours = x/24,
                               days = x, weeks = 7*x)))

    if (nargs() == 1) return(e1)
    # only valid if one of e1 and e2 is a scalar.
    if(inherits(e1, "Date") && inherits(e2, "Date"))
        stop("binary + is not defined for \"Date\" objects")
    if (inherits(e1, "difftime")) e1 <- coerceTimeUnit(e1)
    if (inherits(e2, "difftime")) e2 <- coerceTimeUnit(e2)
    structure(unclass(e1) + unclass(e2), class = "Date")
}

`-.Date` <- function(e1, e2)
{
    coerceTimeUnit <- function(x)
        as.vector(round(switch(attr(x,"units"),
                               secs = x/86400, mins = x/1440, hours = x/24,
                               days = x, weeks = 7*x)))
    if(!inherits(e1, "Date"))
        stop("can only subtract from \"Date\" objects")
    if (nargs() == 1) stop("unary - is not defined for \"Date\" objects")
    if(inherits(e2, "Date")) return(difftime(e1, e2, units="days"))
    if (inherits(e2, "difftime")) e2 <- coerceTimeUnit(e2)
    if(!is.null(attr(e2, "class")))
        stop("can only subtract numbers from \"Date\" objects")
    structure(unclass(as.Date(e1)) - e2, class = "Date")
}

Ops.Date <- function(e1, e2)
{
    if (nargs() == 1)
        stop(gettextf("unary %s not defined for \"Date\" objects", .Generic),
             domain = NA)
    boolean <- switch(.Generic, "<" =, ">" =, "==" =,
                      "!=" =, "<=" =, ">=" = TRUE,
                      FALSE)
    if (!boolean)
        stop(gettextf("%s not defined for \"Date\" objects", .Generic),
             domain = NA)
    ## allow character args to be coerced to dates
    if (is.character(e1)) e1 <- as.Date(e1)
    if (is.character(e2)) e2 <- as.Date(e2)
    NextMethod(.Generic)
}

Math.Date <- function (x, ...)
    stop(gettextf("%s not defined for \"Date\" objects", .Generic),
         domain = NA)

Summary.Date <- function (..., na.rm)
{
    ok <- switch(.Generic, max = , min = , range = TRUE, FALSE)
    if (!ok) stop(gettextf("%s not defined for \"Date\" objects", .Generic),
                  domain = NA)
   val <- NextMethod(.Generic)
    class(val) <- oldClass(list(...)[[1L]])
    val
}

`[.Date` <- function(x, ..., drop = TRUE)
{
    cl <- oldClass(x)
    class(x) <- NULL
    val <- NextMethod("[")
    class(val) <- cl
    val
}

`[[.Date` <- function(x, ..., drop = TRUE)
{
    cl <- oldClass(x)
    class(x) <- NULL
    val <- NextMethod("[[")
    class(val) <- cl
    val
}

`[<-.Date` <- function(x, ..., value)
{
    if(!length(value)) return(x)
    value <- unclass(as.Date(value))
    cl <- oldClass(x)
    class(x) <- NULL
    x <- NextMethod(.Generic)
    class(x) <- cl
    x
}

as.character.Date <- function(x, ...) format(x, ...)

as.data.frame.Date <- as.data.frame.vector

as.list.Date <- function(x, ...)
    lapply(seq_along(x), function(i) x[i])

c.Date <- function(..., recursive=FALSE)
    structure(c(unlist(lapply(list(...), unclass))), class="Date")

mean.Date <- function (x, ...)
    structure(mean(unclass(x), ...), class = "Date")

seq.Date <- function(from, to, by, length.out=NULL, along.with=NULL, ...)
{
    if (missing(from)) stop("'from' must be specified")
    if (!inherits(from, "Date")) stop("'from' must be a \"Date\" object")
        if(length(as.Date(from)) != 1L) stop("'from' must be of length 1")
    if (!missing(to)) {
        if (!inherits(to, "Date")) stop("'to' must be a \"Date\" object")
        if (length(as.Date(to)) != 1L) stop("'to' must be of length 1")
    }
    if (!missing(along.with)) {
        length.out <- length(along.with)
    }  else if (!is.null(length.out)) {
        if (length(length.out) != 1L) stop("'length.out' must be of length 1")
        length.out <- ceiling(length.out)
    }
    status <- c(!missing(to), !missing(by), !is.null(length.out))
    if(sum(status) != 2)
        stop("exactly two of 'to', 'by' and 'length.out' / 'along.with' must be specified")
    if (missing(by)) {
        from <- unclass(as.Date(from))
        to <- unclass(as.Date(to))
        res <- seq.int(from, to, length.out = length.out)
        return(structure(res, class = "Date"))
    }

    if (length(by) != 1L) stop("'by' must be of length 1")
    valid <- 0L
    if (inherits(by, "difftime")) {
        by <- switch(attr(by,"units"), secs = 1/86400, mins = 1/1440,
                     hours = 1/24, days = 1, weeks = 7) * unclass(by)
    } else if(is.character(by)) {
        by2 <- strsplit(by, " ", fixed=TRUE)[[1L]]
        if(length(by2) > 2L || length(by2) < 1L)
            stop("invalid 'by' string")
        valid <- pmatch(by2[length(by2)],
                        c("days", "weeks", "months", "years"))
        if(is.na(valid)) stop("invalid string for 'by'")
        if(valid <= 2) {
            by <- c(1, 7)[valid]
            if (length(by2) == 2L) by <- by * as.integer(by2[1L])
        } else
            by <- if(length(by2) == 2L) as.integer(by2[1L]) else 1
    } else if(!is.numeric(by)) stop("invalid mode for 'by'")
    if(is.na(by)) stop("'by' is NA")

    if(valid <= 2L) { # days or weeks
        from <- unclass(as.Date(from))
        if(!is.null(length.out))
            res <- seq.int(from, by = by, length.out = length.out)
        else {
            to0 <- unclass(as.Date(to))
            ## defeat test in seq.default
            res <- seq.int(0, to0 - from, by) + from
        }
        res <- structure(res, class="Date")
    } else {  # months or years or DSTdays
        r1 <- as.POSIXlt(from)
        if(valid == 4L) {
            if(missing(to)) { # years
                yr <- seq.int(r1$year, by = by, length.out = length.out)
            } else {
                to0 <- as.POSIXlt(to)
                yr <- seq.int(r1$year, to0$year, by)
            }
            r1$year <- yr
            res <- as.Date(r1)
        } else if(valid == 3L) { # months
            if(missing(to)) {
                mon <- seq.int(r1$mon, by = by, length.out = length.out)
            } else {
                to0 <- as.POSIXlt(to)
                mon <- seq.int(r1$mon, 12*(to0$year - r1$year) + to0$mon, by)
            }
            r1$mon <- mon
            res <- as.Date(r1)
        }
    }
    ## can overshoot
    if (!missing(to)) {
        to <- as.Date(to)
        res <- if (by > 0) res[res <= to] else res[res >= to]
    }
    res
}

## *very* similar to cut.POSIXt [ ./datetime.R ] -- keep in sync!
cut.Date <-
    function (x, breaks, labels = NULL, start.on.monday = TRUE,
              right = FALSE, ...)
{
    if(!inherits(x, "Date")) stop("'x' must be a date-time object")
    x <- as.Date(x)

    if (inherits(breaks, "Date")) {
	breaks <- sort(as.Date(breaks))
    } else if(is.numeric(breaks) && length(breaks) == 1L) {
	## specified number of breaks
    } else if(is.character(breaks) && length(breaks) == 1L) {
	by2 <- strsplit(breaks, " ", fixed=TRUE)[[1L]]
	if(length(by2) > 2L || length(by2) < 1L)
	    stop("invalid specification of 'breaks'")
	valid <-
	    pmatch(by2[length(by2)],
		   c("days", "weeks", "months", "years", "quarters"))
	if(is.na(valid)) stop("invalid specification of 'breaks'")
	start <- as.POSIXlt(min(x, na.rm=TRUE))
	if(valid == 1L) incr <- 1L
	if(valid == 2L) {		# weeks
	    start$mday <- start$mday - start$wday
	    if(start.on.monday)
		start$mday <- start$mday + ifelse(start$wday > 0L, 1L, -6L)
            start$isdst <- -1L
	    incr <- 7L
	}
	if(valid == 3L) {		# months
	    start$mday <- 1L
            start$isdst <- -1L
	    end <- as.POSIXlt(max(x, na.rm = TRUE))
	    step <- ifelse(length(by2) == 2L, as.integer(by2[1L]), 1L)
	    end <- as.POSIXlt(end + (31 * step * 86400))
	    end$mday <- 1L
            end$isdst <- -1L
	    breaks <- as.Date(seq(start, end, breaks))
	} else if(valid == 4L) {	# years
	    start$mon <- 0L
	    start$mday <- 1L
            start$isdst <- -1L
	    end <- as.POSIXlt(max(x, na.rm = TRUE))
	    step <- ifelse(length(by2) == 2L, as.integer(by2[1L]), 1L)
	    end <- as.POSIXlt(end + (366 * step * 86400))
	    end$mon <- 0L
	    end$mday <- 1L
            end$isdst <- -1L
	    breaks <- as.Date(seq(start, end, breaks))
	} else if(valid == 5L) {	# quarters
	    qtr <- rep(c(0L, 3L, 6L, 9L), each = 3L)
	    start$mon <- qtr[start$mon + 1L]
	    start$mday <- 1L
            start$isdst <- -1L
	    maxx <- max(x, na.rm = TRUE)
	    end <- as.POSIXlt(maxx)
	    step <- ifelse(length(by2) == 2L, as.integer(by2[1L]), 1L)
	    end <- as.POSIXlt(end + (93 * step * 86400))
	    end$mon <- qtr[end$mon + 1L]
	    end$mday <- 1L
            end$isdst <- -1L
	    breaks <- as.Date(seq(start, end, paste(step * 3L, "months")))
	    ## 93 days ahead could give an empty level, so
	    lb <- length(breaks)
	    if(maxx < breaks[lb-1]) breaks <- breaks[-lb]
	} else {
	    start <- as.Date(start)
	    if (length(by2) == 2L) incr <- incr * as.integer(by2[1L])
	    maxx <- max(x, na.rm = TRUE)
	    breaks <- seq(start, maxx + incr, breaks)
	    breaks <- breaks[seq_len(1L+max(which(breaks <= maxx)))]
	}
    } else stop("invalid specification of 'breaks'")
    res <- cut(unclass(x), unclass(breaks), labels = labels,
	       right = right, ...)
    if(is.null(labels)) {
	levels(res) <-
	    as.character(if (is.numeric(breaks)) x[!duplicated(res)]
			 else breaks[-length(breaks)])
    }
    res
}

julian.Date <- function(x, origin = as.Date("1970-01-01"), ...)
{
    if(length(origin) != 1L) stop("'origin' must be of length one")
    structure(unclass(x) - unclass(origin), "origin" = origin)
}

weekdays.Date <- function(x, abbreviate = FALSE)
    format(x, ifelse(abbreviate, "%a", "%A"))

months.Date <- function(x, abbreviate = FALSE)
    format(x, ifelse(abbreviate, "%b", "%B"))

quarters.Date <- function(x, ...)
{
    x <- (as.POSIXlt(x)$mon) %/% 3L
    paste0("Q", x+1L)
}

## These only make sense for negative digits, but still ...
round.Date <- function(x, ...)
{
    cl <- oldClass(x)
    class(x) <- NULL
    val <- NextMethod()
    class(val) <- cl
    val
}

## must avoid truncating forwards dates prior to 1970-01-01.
trunc.Date <- function(x, ...) round(x - 0.4999999)

rep.Date <- function(x, ...)
{
    y <- NextMethod()
    structure(y, class="Date")
}

diff.Date <- function (x, lag = 1L, differences = 1L, ...)
{
    ismat <- is.matrix(x)
    xlen <- if (ismat) dim(x)[1L] else length(x)
    if (length(lag) > 1L || length(differences) > 1L || lag < 1L || differences < 1L)
        stop("'lag' and 'differences' must be integers >= 1")
    if (lag * differences >= xlen)
        return(structure(numeric(), class="difftime", units="days"))
    r <- x
    i1 <- -seq_len(lag)
    if (ismat)
        for (i in seq_len(differences)) r <- r[i1, , drop = FALSE] -
            r[-nrow(r):-(nrow(r) - lag + 1L), , drop = FALSE]
    else for (i in seq_len(differences))
        r <- r[i1] - r[-length(r):-(length(r) - lag + 1L)]
    r
}

## ---- additions in 2.6.0 -----

is.numeric.Date <- function(x) FALSE

## ---- additions in 2.8.0 -----

split.Date <-
function(x, f, drop = FALSE, ...)
{
    y <- split.default(as.integer(x), f, drop = drop)
    for(i in seq_along(y)) class(y[[i]]) <- "Date"
    y
}

xtfrm.Date <- function(x) as.numeric(x)
#  File src/library/base/R/datetime.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

Sys.time <- function() .POSIXct(.Internal(Sys.time()))

Sys.timezone <- function() Sys.getenv("TZ", names = FALSE)

as.POSIXlt <- function(x, tz = "", ...) UseMethod("as.POSIXlt")

as.POSIXlt.Date <- function(x, ...)
{
    ## <FIXME>
    ## Move names handling to C code eventually ...
    y <- .Internal(Date2POSIXlt(x))
    names(y$year) <- names(x)
    y
    ## </FIXME>
}

as.POSIXlt.date <- as.POSIXlt.dates <- function(x, ...)
    as.POSIXlt(as.POSIXct(x), ...)

as.POSIXlt.POSIXct <- function(x, tz = "", ...)
{
    tzone <- attr(x, "tzone")
    if((missing(tz) || is.null(tz)) && !is.null(tzone)) tz <- tzone[1L]
    ## <FIXME>
    ## Move names handling to C code eventually ...
    y <- .Internal(as.POSIXlt(x, tz))
    names(y$year) <- names(x)
    y
    ## </FIXME>
}

as.POSIXlt.factor <- function(x, ...)
{
    y <- as.POSIXlt(as.character(x), ...)
    names(y$year) <- names(x)
    y
}

as.POSIXlt.character <- function(x, tz = "", format, ...)
{
    x <- unclass(x) # precaution PR7826
    if(!missing(format)) {
        res <- strptime(x, format, tz=tz)
        if(nzchar(tz)) attr(res, "tzone") <- tz
        return(res)
    }
    xx <- x[!is.na(x)]
    if (!length(xx)) {
        res <- strptime(x, "%Y/%m/%d")
        if(nzchar(tz)) attr(res, "tzone") <- tz
        return(res)
    } else if(all(!is.na(strptime(xx, f <- "%Y-%m-%d %H:%M:%OS", tz=tz))) ||
            all(!is.na(strptime(xx, f <- "%Y/%m/%d %H:%M:%OS", tz=tz))) ||
            all(!is.na(strptime(xx, f <- "%Y-%m-%d %H:%M", tz=tz))) ||
            all(!is.na(strptime(xx, f <- "%Y/%m/%d %H:%M", tz=tz))) ||
            all(!is.na(strptime(xx, f <- "%Y-%m-%d", tz=tz))) ||
            all(!is.na(strptime(xx, f <- "%Y/%m/%d", tz=tz)))
            ) {
        res <- strptime(x, f, tz=tz)
        if(nzchar(tz)) attr(res, "tzone") <- tz
        return(res)
    }
    stop("character string is not in a standard unambiguous format")
}

as.POSIXlt.numeric <- function(x, tz = "", origin, ...)
{
    if(missing(origin)) stop("'origin' must be supplied")
    as.POSIXlt(as.POSIXct(origin, tz = "UTC", ...) + x, tz = tz)
}

as.POSIXlt.default <- function(x, tz = "", ...)
{

    if(inherits(x, "POSIXlt")) return(x)
    if(is.logical(x) && all(is.na(x)))
        return(as.POSIXlt(as.POSIXct.default(x), tz=tz))
    stop(gettextf("do not know how to convert '%s' to class %s",
                  deparse(substitute(x)),
                  dQuote("POSIXlt")),
         domain = NA)
}

as.POSIXct <- function(x, tz = "", ...) UseMethod("as.POSIXct")

as.POSIXct.Date <- function(x, ...) .POSIXct(unclass(x)*86400)


## convert from package date
as.POSIXct.date <- function(x, ...)
{
    if(inherits(x, "date")) {
        x <- (x - 3653) * 86400 # origin 1960-01-01
        return(.POSIXct(x))
    } else stop(gettextf("'%s' is not a \"date\" object",
                         deparse(substitute(x)) ))
}

## convert from package chron
as.POSIXct.dates <- function(x, ...)
{
    if(inherits(x, "dates")) {
        z <- attr(x, "origin")
        x <- as.numeric(x) * 86400
        if(length(z) == 3L && is.numeric(z))
            x  <- x + as.numeric(ISOdate(z[3L], z[1L], z[2L], 0))
        return(.POSIXct(x))
    } else stop(gettextf("'%s' is not a \"dates\" object",
                         deparse(substitute(x)) ))
}

as.POSIXct.POSIXlt <- function(x, tz = "", ...)
{
    tzone <- attr(x, "tzone")
    if(missing(tz) && !is.null(tzone)) tz <- tzone[1L]
    ## <FIXME>
    ## Move names handling to C code eventually ...
    y <- .Internal(as.POSIXct(x, tz))
    names(y) <- names(x$year)
    .POSIXct(y, tz)
    ## </FIXME>
}

as.POSIXct.numeric <- function(x, tz = "", origin, ...)
{
    if(missing(origin)) stop("'origin' must be supplied")
    .POSIXct(as.POSIXct(origin, tz="GMT", ...) + x, tz)
}

as.POSIXct.default <- function(x, tz = "", ...)
{
    if(inherits(x, "POSIXct")) return(x)
    if(is.character(x) || is.factor(x))
	return(as.POSIXct(as.POSIXlt(x, tz, ...), tz, ...))
    if(is.logical(x) && all(is.na(x)))
        return(.POSIXct(as.numeric(x)))
    stop(gettextf("do not know how to convert '%s' to class %s",
                  deparse(substitute(x)),
                  dQuote("POSIXct")),
         domain = NA)
}

as.double.POSIXlt <- function(x, ...) as.double(as.POSIXct(x))

## POSIXlt is not primarily a list, but primarily an abstract vector of
## time stamps:
length.POSIXlt <- function(x) length(x[[1L]])

format.POSIXlt <- function(x, format = "", usetz = FALSE, ...)
{
    if(!inherits(x, "POSIXlt")) stop("wrong class")
    if(format == "") {
        ## need list [ method here.
        times <- unlist(unclass(x)[1L:3L])
        secs <- x$sec; secs <- secs[!is.na(secs)]
        np <- getOption("digits.secs")
        if(is.null(np)) np <- 0L else np <- min(6L, np)
        if(np >= 1L)
            for (i in seq_len(np)- 1L)
                if(all( abs(secs - round(secs, i)) < 1e-6 )) {
                    np <- i
                    break
                }
        format <- if(all(times[!is.na(times)] == 0)) "%Y-%m-%d"
        else if(np == 0L) "%Y-%m-%d %H:%M:%S"
        else paste0("%Y-%m-%d %H:%M:%OS", np)
    }
    ## <FIXME>
    ## Move names handling to C code eventually ...
    y <- .Internal(format.POSIXlt(x, format, usetz))
    names(y) <- names(x$year)
    y
    ## </FIXME>
}

## prior to 2.9.0 the same as format.POSIXlt.
## now more or less the same as format.POSIXct but also works for Dates.
strftime <- function(x, format = "", tz = "", usetz = FALSE, ...)
    format(as.POSIXlt(x, tz = tz), format = format, usetz = usetz, ...)

strptime <- function(x, format, tz = "")
{
    ## <FIXME>
    ## Move names handling to C code eventually ...
    y <- .Internal(strptime(as.character(x), format, tz))
    ## Assuming we can rely on the names of x ...
    names(y$year) <- names(x)
    y
    ## </FIXME>
}

format.POSIXct <- function(x, format = "", tz = "", usetz = FALSE, ...)
{
    if(!inherits(x, "POSIXct")) stop("wrong class")
    if(missing(tz) && !is.null(tzone <- attr(x, "tzone"))) tz <- tzone
    structure(format.POSIXlt(as.POSIXlt(x, tz), format, usetz, ...),
              names=names(x))
}

## could handle arrays for max.print
print.POSIXct <- function(x, ...)
{
    max.print <- getOption("max.print", 9999L)
    if(max.print < length(x)) {
        print(format(x[seq_len(max.print)], usetz=TRUE), ...)
        cat(' [ reached getOption("max.print") -- omitted',
            length(x) - max.print, 'entries ]\n')
    } else print(format(x, usetz=TRUE), ...)
    invisible(x)
}

print.POSIXlt <- function(x, ...)
{
    max.print <- getOption("max.print", 9999L)
    if(max.print < length(x)) {
        print(format(x[seq_len(max.print)], usetz=TRUE), ...)
        cat(' [ reached getOption("max.print") -- omitted',
            length(x) - max.print, 'entries ]\n')
   } else print(format(x, usetz=TRUE), ...)
    invisible(x)
}

summary.POSIXct <- function(object, digits = 15L, ...)
{
    x <- summary.default(unclass(object), digits = digits, ...)
    if(m <- match("NA's", names(x), 0)) {
        NAs <- as.integer(x[m])
        x <- x[-m]
        attr(x, "NAs") <- NAs
    }
    class(x) <- c("summaryDefault", "table", oldClass(object))
    attr(x, "tzone") <- attr(object, "tzone")
    x
}

summary.POSIXlt <- function(object, digits = 15, ...)
    summary(as.POSIXct(object), digits = digits, ...)


`+.POSIXt` <- function(e1, e2)
{
    ## need to drop "units" attribute here
    coerceTimeUnit <- function(x)
        as.vector(switch(attr(x,"units"),
                         secs = x, mins = 60*x, hours = 60*60*x,
                         days = 60*60*24*x, weeks = 60*60*24*7*x))

    if (nargs() == 1) return(e1)
    # only valid if one of e1 and e2 is a scalar/difftime
    if(inherits(e1, "POSIXt") && inherits(e2, "POSIXt"))
        stop("binary '+' is not defined for \"POSIXt\" objects")
    if(inherits(e1, "POSIXlt")) e1 <- as.POSIXct(e1)
    if(inherits(e2, "POSIXlt")) e2 <- as.POSIXct(e2)
    if (inherits(e1, "difftime")) e1 <- coerceTimeUnit(e1)
    if (inherits(e2, "difftime")) e2 <- coerceTimeUnit(e2)
    .POSIXct(unclass(e1) + unclass(e2), check_tzones(e1, e2))
}

`-.POSIXt` <- function(e1, e2)
{
    ## need to drop "units" attribute here
    coerceTimeUnit <- function(x)
        as.vector(switch(attr(x,"units"),
                         secs = x, mins = 60*x, hours = 60*60*x,
                         days = 60*60*24*x, weeks = 60*60*24*7*x))
    if(!inherits(e1, "POSIXt"))
        stop("can only subtract from \"POSIXt\" objects")
    if (nargs() == 1) stop("unary '-' is not defined for \"POSIXt\" objects")
    if(inherits(e2, "POSIXt")) return(difftime(e1, e2))
    if (inherits(e2, "difftime")) e2 <- coerceTimeUnit(e2)
    if(!is.null(attr(e2, "class")))
        stop("can only subtract numbers from \"POSIXt\" objects")
    e1 <- as.POSIXct(e1)
    .POSIXct(unclass(e1) - e2, attr(e1, "tzone"))
}

Ops.POSIXt <- function(e1, e2)
{
    if (nargs() == 1)
        stop(gettextf("unary '%s' not defined for \"POSIXt\" objects",
                      .Generic), domain = NA)
    boolean <- switch(.Generic, "<" = , ">" = , "==" = ,
                      "!=" = , "<=" = , ">=" = TRUE, FALSE)
    if (!boolean)
        stop(gettextf("'%s' not defined for \"POSIXt\" objects", .Generic),
             domain = NA)
    if(inherits(e1, "POSIXlt") || is.character(e1)) e1 <- as.POSIXct(e1)
    if(inherits(e2, "POSIXlt") || is.character(e1)) e2 <- as.POSIXct(e2)
    check_tzones(e1, e2)
    NextMethod(.Generic)
}

Math.POSIXt <- function (x, ...)
{
    stop(gettextf("'%s' not defined for \"POSIXt\" objects", .Generic),
         domain = NA)
}

check_tzones <- function(...)
{
    tzs <- unique(sapply(list(...), function(x) {
        y <- attr(x, "tzone")
        if(is.null(y)) "" else y[1L]
    }))
    tzs <- tzs[tzs != ""]
    if(length(tzs) > 1L)
        warning("'tzone' attributes are inconsistent")
    if(length(tzs)) tzs[1L] else NULL
}

Summary.POSIXct <- function (..., na.rm)
{
    ok <- switch(.Generic, max = , min = , range = TRUE, FALSE)
    if (!ok)
        stop(gettextf("'%s' not defined for \"POSIXt\" objects", .Generic),
             domain = NA)
    args <- list(...)
    tz <- do.call("check_tzones", args)
    val <- NextMethod(.Generic)
    class(val) <- oldClass(args[[1L]])
    attr(val, "tzone") <- tz
    val
}

Summary.POSIXlt <- function (..., na.rm)
{
    ok <- switch(.Generic, max = , min = , range = TRUE, FALSE)
    if (!ok)
        stop(gettextf("'%s' not defined for \"POSIXt\" objects", .Generic),
             domain = NA)
    args <- list(...)
    tz <- do.call("check_tzones", args)
    args <- lapply(args, as.POSIXct)
    val <- do.call(.Generic, c(args, na.rm = na.rm))
    as.POSIXlt(.POSIXct(val, tz))
}

`[.POSIXct` <-
function(x, ..., drop = TRUE)
{
    cl <- oldClass(x)
    class(x) <- NULL
    val <- NextMethod("[")
    class(val) <- cl
    attr(val, "tzone") <- attr(x, "tzone")
    val
}

`[[.POSIXct` <-
function(x, ..., drop = TRUE)
{
    cl <- oldClass(x)
    class(x) <- NULL
    val <- NextMethod("[[")
    class(val) <- cl
    attr(val, "tzone") <- attr(x, "tzone")
    val
}

`[<-.POSIXct` <-
function(x, ..., value) {
    if(!length(value)) return(x)
    value <- unclass(as.POSIXct(value))
    cl <- oldClass(x)
    tz <- attr(x, "tzone")
    class(x) <- NULL
    x <- NextMethod(.Generic)
    class(x) <- cl
    attr(x, "tzone") <- tz
    x
}

as.character.POSIXt <- function(x, ...) format(x, ...)

as.data.frame.POSIXct <- as.data.frame.vector

as.list.POSIXct <- function(x, ...)
{
    nms <- names(x)
    names(x) <- NULL
    y <- lapply(seq_along(x), function(i) x[i])
    names(y) <- nms
    y
}

is.na.POSIXlt <- function(x) is.na(as.POSIXct(x))

## <FIXME> check the argument validity
## This is documented to remove the timezone
c.POSIXct <- function(..., recursive=FALSE)
    .POSIXct(c(unlist(lapply(list(...), unclass))))

## we need conversion to POSIXct as POSIXlt objects can be in different tz.
c.POSIXlt <- function(..., recursive=FALSE)
    as.POSIXlt(do.call("c", lapply(list(...), as.POSIXct)))

## force absolute comparisons
all.equal.POSIXct <- function(target, current, ..., tolerance = 1e-3, scale=1)
{
    check_tzones(target, current)
    NextMethod("all.equal", tolerance = tolerance, scale = scale)
}


ISOdatetime <- function(year, month, day, hour, min, sec, tz="")
{
    if(min(sapply(list(year, month, day, hour, min, sec), length)) == 0L)
        .POSIXct(numeric(), tz=tz)
    else {
        x <- paste(year, month, day, hour, min, sec, sep = "-")
        as.POSIXct(strptime(x, "%Y-%m-%d-%H-%M-%OS", tz=tz), tz=tz)
    }
}

ISOdate <- function(year, month, day, hour=12, min=0, sec=0, tz="GMT")
    ISOdatetime(year, month, day, hour, min, sec, tz)

as.matrix.POSIXlt <- function(x, ...)
{
    as.matrix(as.data.frame(unclass(x)), ...)
}

mean.POSIXct <- function (x, ...)
    .POSIXct(mean(unclass(x), ...), attr(x, "tzone"))

mean.POSIXlt <- function (x, ...)
    as.POSIXlt(mean(as.POSIXct(x), ...))

## ----- difftime -----

difftime <-
    function(time1, time2, tz,
             units = c("auto", "secs", "mins", "hours", "days", "weeks"))
{
    if (missing(tz)) {
        time1 <- as.POSIXct(time1)
        time2 <- as.POSIXct(time2)
    } else {
        ## Wishlist PR#14182
        time1 <- as.POSIXct(time1, tz = tz)
        time2 <- as.POSIXct(time2, tz = tz)
    }
    z <- unclass(time1) - unclass(time2)
    units <- match.arg(units)
    if(units == "auto") {
        if(all(is.na(z))) units <- "secs"
        else {
            zz <- min(abs(z),na.rm=TRUE)
            if(is.na(zz) || zz < 60) units <- "secs"
            else if(zz < 3600) units <- "mins"
            else if(zz < 86400) units <- "hours"
            else units <- "days"
        }
    }
    switch(units,
           "secs" = .difftime(z, units="secs"),
           "mins" = .difftime(z/60, units="mins"),
           "hours" = .difftime(z/3600, units="hours"),
           "days" = .difftime(z/86400, units="days"),
           "weeks" = .difftime(z/(7*86400), units="weeks")
           )
}

## "difftime" constructor
## Martin Maechler, Date: 16 Sep 2002
## Numeric input version Peter Dalgaard, December 2006
as.difftime <- function(tim, format="%X", units="auto")
{
    if (inherits(tim, "difftime")) return(tim)
    if (is.character(tim)){
        difftime(strptime(tim, format=format),
             strptime("0:0:0", format="%X"), units=units)
    } else {
        if (!is.numeric(tim)) stop("'tim' is not character or numeric")
	if (units=="auto") stop("need explicit units for numeric conversion")
        if (!(units %in% c("secs", "mins", "hours", "days", "weeks")))
	    stop("invalid units specified")
        structure(tim, units=units, class="difftime")
    }
}

### For now, these have only difftime methods, but you never know...
units <- function(x) UseMethod("units")

`units<-` <- function(x, value) UseMethod("units<-")

units.difftime <- function(x) attr(x, "units")

`units<-.difftime` <- function(x, value)
{
    from <- units(x)
    if (from == value) return(x)
    if (!(value %in% c("secs", "mins", "hours", "days", "weeks")))
        stop("invalid units specified")
    sc <- cumprod(c(secs=1, mins=60, hours=60, days=24, weeks=7))
    newx <- unclass(x) * as.vector(sc[from]/sc[value])
    .difftime(newx, value)
}

as.double.difftime <- function(x, units="auto", ...)
{
    if (units != "auto") units(x) <- units
    as.vector(x, "double")
}

as.data.frame.difftime <- as.data.frame.vector

format.difftime <- function(x,...) paste(format(unclass(x),...), units(x))



print.difftime <- function(x, digits = getOption("digits"), ...)
{
    if(is.array(x) || length(x) > 1L) {
        cat("Time differences in ", attr(x, "units"), "\n", sep = "")
        y <- unclass(x); attr(y, "units") <- NULL
        print(y)
    }
    else
        cat("Time difference of ", format(unclass(x), digits = digits), " ",
            attr(x, "units"), "\n", sep = "")

    invisible(x)
}

`[.difftime` <- function(x, ..., drop = TRUE)
{
    cl <- oldClass(x)
    class(x) <- NULL
    val <- NextMethod("[")
    class(val) <- cl
    attr(val, "units") <- attr(x, "units")
    val
}

Ops.difftime <- function(e1, e2)
{
    coerceTimeUnit <- function(x)
    {
        switch(attr(x, "units"),
               secs = x, mins = 60*x, hours = 60*60*x,
               days = 60*60*24*x, weeks = 60*60*24*7*x)
    }
    if (nargs() == 1) {
        switch(.Generic, "+" = {}, "-" = {e1[] <- -unclass(e1)},
               stop(gettextf("unary '%s' not defined for \"difftime\" objects",
                             .Generic), domain = NA, call. = FALSE)
               )
        return(e1)
    }
    boolean <- switch(.Generic, "<" = , ">" = , "==" = ,
                      "!=" = , "<=" = , ">=" = TRUE, FALSE)
    if (boolean) {
        ## assume user knows what he/she is doing if not both difftime
        if(inherits(e1, "difftime") && inherits(e2, "difftime")) {
            e1 <- coerceTimeUnit(e1)
            e2 <- coerceTimeUnit(e2)
        }
        NextMethod(.Generic)
    } else if(.Generic == "+" || .Generic == "-") {
        if(inherits(e1, "difftime") && !inherits(e2, "difftime"))
            return(structure(NextMethod(.Generic),
                             units = attr(e1, "units"), class = "difftime"))
        if(!inherits(e1, "difftime") && inherits(e2, "difftime"))
            return(structure(NextMethod(.Generic),
                             units = attr(e2, "units"), class = "difftime"))
        u1 <- attr(e1, "units")
        if(attr(e2, "units") == u1) {
            structure(NextMethod(.Generic), units=u1, class="difftime")
        } else {
            e1 <- coerceTimeUnit(e1)
            e2 <- coerceTimeUnit(e2)
            structure(NextMethod(.Generic), units="secs", class="difftime")
        }
    } else {
        ## '*' is covered by a specific method
        stop(gettextf("'%s' not defined for \"difftime\" objects", .Generic),
             domain = NA)
    }
}

`*.difftime` <- function (e1, e2)
{
    ## need one scalar, one difftime.
    if(inherits(e1, "difftime") && inherits(e2, "difftime"))
        stop("both arguments of * cannot be \"difftime\" objects")
    if(inherits(e2, "difftime")) {tmp <- e1; e1 <- e2; e2 <- tmp}
    .difftime(e2 * unclass(e1), attr(e1, "units"))
}

`/.difftime` <- function (e1, e2)
{
    ## need one scalar, one difftime.
    if(inherits(e2, "difftime"))
        stop("second argument of / cannot be a \"difftime\" object")
    .difftime(unclass(e1) / e2, attr(e1, "units"))
}

## "Math": some methods should work; the other ones are meaningless :
Math.difftime <- function (x, ...)
{
    switch(.Generic,
           "abs" =, "sign" =, "floor" =, "ceiling" =, "trunc" =,
           "round" =, "signif" = {
               units <- attr(x, "units")
               .difftime(NextMethod(), units)
           },
           ### otherwise :
           stop(gettextf("'%s' not defined for \"difftime\" objects", .Generic),
                domain = NA))
}


mean.difftime <- function (x, ...)
    .difftime(mean(unclass(x), ...), attr(x, "units"))

Summary.difftime <- function (..., na.rm)
{
    ## FIXME: this could return in the smallest of the units of the inputs.
    coerceTimeUnit <- function(x)
    {
        as.vector(switch(attr(x,"units"),
                         secs = x, mins = 60*x, hours = 60*60*x,
                         days = 60*60*24*x, weeks = 60*60*24*7*x))
    }
    ok <- switch(.Generic, max = , min = , sum=, range = TRUE, FALSE)
    if (!ok)
        stop(gettextf("'%s' not defined for \"difftime\" objects", .Generic),
             domain = NA)
    x <- list(...)
    Nargs <- length(x)
    if(Nargs == 0) {
        .difftime(do.call(.Generic), "secs")
    } else {
        units <- sapply(x, function(x) attr(x, "units"))
        if(all(units == units[1L])) {
            args <- c(lapply(x, as.vector), na.rm = na.rm)
        } else {
            args <- c(lapply(x, coerceTimeUnit), na.rm = na.rm)
            units <- "secs"
        }
        .difftime(do.call(.Generic, args), units[[1L]])
    }
}


## ----- convenience functions -----

seq.POSIXt <-
    function(from, to, by, length.out = NULL, along.with = NULL, ...)
{
    if (missing(from)) stop("'from' must be specified")
    if (!inherits(from, "POSIXt")) stop("'from' must be a \"POSIXt\" object")
    cfrom <- as.POSIXct(from)
    if(length(cfrom) != 1L) stop("'from' must be of length 1")
    tz <- attr(cfrom , "tzone")
    if (!missing(to)) {
        if (!inherits(to, "POSIXt")) stop("'to' must be a \"POSIXt\" object")
        if (length(as.POSIXct(to)) != 1) stop("'to' must be of length 1")
    }
    if (!missing(along.with)) {
        length.out <- length(along.with)
    }  else if (!is.null(length.out)) {
        if (length(length.out) != 1L) stop("'length.out' must be of length 1")
        length.out <- ceiling(length.out)
    }
    status <- c(!missing(to), !missing(by), !is.null(length.out))
    if(sum(status) != 2L)
        stop("exactly two of 'to', 'by' and 'length.out' / 'along.with' must be specified")
    if (missing(by)) {
        from <- unclass(cfrom)
        to <- unclass(as.POSIXct(to))
        ## Till (and incl.) 1.6.0 :
        ##- incr <- (to - from)/length.out
        ##- res <- seq.default(from, to, incr)
        res <- seq.int(from, to, length.out = length.out)
        return(.POSIXct(res, tz))
    }

    if (length(by) != 1L) stop("'by' must be of length 1")
    valid <- 0L
    if (inherits(by, "difftime")) {
        by <- switch(attr(by,"units"), secs = 1, mins = 60, hours = 3600,
                     days = 86400, weeks = 7*86400) * unclass(by)
    } else if(is.character(by)) {
        by2 <- strsplit(by, " ", fixed=TRUE)[[1L]]
        if(length(by2) > 2L || length(by2) < 1L)
            stop("invalid 'by' string")
        valid <- pmatch(by2[length(by2)],
                        c("secs", "mins", "hours", "days", "weeks",
                          "months", "years", "DSTdays"))
        if(is.na(valid)) stop("invalid string for 'by'")
        if(valid <= 5L) {
            by <- c(1, 60, 3600, 86400, 7*86400)[valid]
            if (length(by2) == 2L) by <- by * as.integer(by2[1L])
        } else
            by <- if(length(by2) == 2L) as.integer(by2[1L]) else 1
    } else if(!is.numeric(by)) stop("invalid mode for 'by'")
    if(is.na(by)) stop("'by' is NA")

    if(valid <= 5L) { # secs, mins, hours, days, weeks
        from <- unclass(as.POSIXct(from))
        if(!is.null(length.out))
            res <- seq.int(from, by = by, length.out = length.out)
        else {
            to0 <- unclass(as.POSIXct(to))
            ## defeat test in seq.default
            res <- seq.int(0, to0 - from, by) + from
        }
        return(.POSIXct(res, tz))
    } else {  # months or years or DSTdays
        r1 <- as.POSIXlt(from)
        if(valid == 7L) { # years
            if(missing(to)) { # years
                yr <- seq.int(r1$year, by = by, length.out = length.out)
            } else {
                to <- as.POSIXlt(to)
                yr <- seq.int(r1$year, to$year, by)
            }
            r1$year <- yr
        } else if(valid == 6L) { # months
            if(missing(to)) {
                mon <- seq.int(r1$mon, by = by, length.out = length.out)
            } else {
                to0 <- as.POSIXlt(to)
                mon <- seq.int(r1$mon, 12*(to0$year - r1$year) + to0$mon, by)
            }
            r1$mon <- mon
        } else if(valid == 8L) { # DSTdays
            if(!missing(to)) {
                ## We might have a short day, so need to over-estimate.
                length.out <- 2L + floor((unclass(as.POSIXct(to)) -
                                          unclass(as.POSIXct(from)))/86400)
            }
            r1$mday <- seq.int(r1$mday, by = by, length.out = length.out)
        }
	r1$isdst <- -1L
	res <- as.POSIXct(r1)
	## now shorten if necessary.
	if(!missing(to)) {
	    to <- as.POSIXct(to)
	    res <- if(by > 0) res[res <= to] else res[res >= to]
	}
	res
    }
}

## *very* similar to cut.Date [ ./dates.R ] -- keep in sync!
cut.POSIXt <-
    function (x, breaks, labels = NULL, start.on.monday = TRUE,
              right = FALSE, ...)
{
    if(!inherits(x, "POSIXt")) stop("'x' must be a date-time object")
    x <- as.POSIXct(x)

    if (inherits(breaks, "POSIXt")) {
	breaks <- sort(as.POSIXct(breaks))
    } else if(is.numeric(breaks) && length(breaks) == 1L) {
	## specified number of breaks
    } else if(is.character(breaks) && length(breaks) == 1L) {
        by2 <- strsplit(breaks, " ", fixed=TRUE)[[1L]]
        if(length(by2) > 2L || length(by2) < 1L)
            stop("invalid specification of 'breaks'")
	valid <-
	    pmatch(by2[length(by2)],
		   c("secs", "mins", "hours", "days", "weeks",
		     "months", "years", "DSTdays", "quarters"))
	if(is.na(valid)) stop("invalid specification of 'breaks'")
	start <- as.POSIXlt(min(x, na.rm=TRUE))
	incr <- 1
	if(valid > 1L) { start$sec <- 0L; incr <- 60 }
	if(valid > 2L) { start$min <- 0L; incr <- 3600 }
        ## start of day need not be on the same DST, PR#14208
	if(valid > 3L) { start$hour <- 0L; start$isdst <- -1L; incr <- 86400 }
	if(valid == 5L) {               # weeks
	    start$mday <- start$mday - start$wday
	    if(start.on.monday)
		start$mday <- start$mday + ifelse(start$wday > 0L, 1L, -6L)
	    incr <- 7*86400
	}
        if(valid == 8L) incr <- 25*3600 # DSTdays
        if(valid == 6L) {               # months
            start$mday <- 1L
            end <- as.POSIXlt(max(x, na.rm = TRUE))
            step <- ifelse(length(by2) == 2L, as.integer(by2[1L]), 1L)
            end <- as.POSIXlt(end + (31 * step * 86400))
            end$mday <- 1L
            end$isdst <- -1L
            breaks <- seq(start, end, breaks)
        } else if(valid == 7L) {        # years
            start$mon <- 0L
            start$mday <- 1L
            end <- as.POSIXlt(max(x, na.rm = TRUE))
            step <- ifelse(length(by2) == 2L, as.integer(by2[1L]), 1L)
            end <- as.POSIXlt(end + (366 * step* 86400))
            end$mon <- 0L
            end$mday <- 1L
            end$isdst <- -1L
            breaks <- seq(start, end, breaks)
        } else if(valid == 9L) {        # quarters
            qtr <- rep(c(0L, 3L, 6L, 9L), each = 3L)
            start$mon <- qtr[start$mon + 1L]
            start$mday <- 1L
            maxx <- max(x, na.rm = TRUE)
            end <- as.POSIXlt(maxx)
            step <- ifelse(length(by2) == 2L, as.integer(by2[1L]), 1L)
            end <- as.POSIXlt(end + (93 * step * 86400))
            end$mon <- qtr[end$mon + 1L]
            end$mday <- 1L
            end$isdst <- -1L
            breaks <- seq(start, end, paste(step * 3, "months"))
            ## 93 days ahead could give an empty level, so
            lb <- length(breaks)
            if(maxx < breaks[lb-1]) breaks <- breaks[-lb]
        } else {                        # weeks or shorter
            if (length(by2) == 2L) incr <- incr * as.integer(by2[1L])
            maxx <- max(x, na.rm = TRUE)
            breaks <- seq(start, maxx + incr, breaks)
            breaks <- breaks[seq_len(1+max(which(breaks <= maxx)))]
        }
    } else stop("invalid specification of 'breaks'")
    res <- cut(unclass(x), unclass(breaks), labels = labels,
               right = right, ...)
    if(is.null(labels)) {
	levels(res) <-
	    as.character(if (is.numeric(breaks)) x[!duplicated(res)]
			 else breaks[-length(breaks)])
    }
    res
}

julian <- function(x, ...) UseMethod("julian")

julian.POSIXt <- function(x, origin = as.POSIXct("1970-01-01", tz="GMT"), ...)
{
    origin <- as.POSIXct(origin)
    if(length(origin) != 1L) stop("'origin' must be of length one")
    res <- difftime(as.POSIXct(x), origin, units = "days")
    structure(res, "origin" = origin)
}

weekdays <- function(x, abbreviate) UseMethod("weekdays")
weekdays.POSIXt <- function(x, abbreviate = FALSE)
{
    format(x, ifelse(abbreviate, "%a", "%A"))
}

months <- function(x, abbreviate) UseMethod("months")
months.POSIXt <- function(x, abbreviate = FALSE)
{
    format(x, ifelse(abbreviate, "%b", "%B"))
}

quarters <- function(x, abbreviate) UseMethod("quarters")
quarters.POSIXt <- function(x, ...)
{
    x <- (as.POSIXlt(x)$mon)%/%3
    paste0("Q", x+1)
}

trunc.POSIXt <- function(x, units=c("secs", "mins", "hours", "days"), ...)
{
    units <- match.arg(units)
    x <- as.POSIXlt(x)
    if(length(x$sec))
	switch(units,
	       "secs" = {x$sec <- trunc(x$sec)},
	       "mins" = {x$sec[] <- 0},
	       "hours" = {x$sec[] <- 0; x$min[] <- 0L},
               ## start of day need not be on the same DST.
	       "days" = {x$sec[] <- 0; x$min[] <- 0L; x$hour[] <- 0L; x$isdst[] <- -1L}
	       )
    x
}

round.POSIXt <- function(x, units=c("secs", "mins", "hours", "days"))
{
    ## this gets the default from the generic, as that has two args.
    if(is.numeric(units) && units == 0.0) units <-"secs"
    units <- match.arg(units)
    x <- as.POSIXct(x)
    x <- x + switch(units,
                    "secs" = 0.5, "mins" = 30, "hours" = 1800, "days" = 43200)
    trunc.POSIXt(x, units = units)
}

## ---- additions in 1.5.0 -----

`[.POSIXlt` <- function(x, ..., drop = TRUE)
{
    val <- lapply(X = x, FUN = "[", ..., drop = drop)
    attributes(val) <- attributes(x) # need to preserve timezones
    val
}

`[<-.POSIXlt` <- function(x, i, value)
{
    if(!length(value)) return(x)
    value <- unclass(as.POSIXlt(value))
    cl <- oldClass(x)
    class(x) <- NULL
    for(n in names(x)) x[[n]][i] <- value[[n]]
    class(x) <- cl
    x
}

as.data.frame.POSIXlt <- function(x, row.names = NULL, optional = FALSE, ...)
{
    value <- as.data.frame.POSIXct(as.POSIXct(x), row.names, optional, ...)
    if (!optional)
        names(value) <- deparse(substitute(x))[[1L]]
    value
}

## ---- additions in 1.8.0 -----

rep.POSIXct <- function(x, ...)
{
    y <- NextMethod()
    .POSIXct(y, attr(x, "tzone"))
}

rep.POSIXlt <- function(x, ...)
{
    y <- lapply(X = x, FUN = rep, ...)
    attributes(y) <- attributes(x)
    y
}

diff.POSIXt <- function (x, lag = 1L, differences = 1L, ...)
{
    ismat <- is.matrix(x)
    r <- if(inherits(x, "POSIXlt")) as.POSIXct(x) else x
    xlen <- if (ismat) dim(x)[1L] else length(r)
    if (length(lag) > 1L || length(differences) > 1L || lag < 1L || differences < 1L)
        stop("'lag' and 'differences' must be integers >= 1")
    if (lag * differences >= xlen) return(.difftime(numeric(), "secs"))
    i1 <- -seq_len(lag)
    if (ismat) for (i in seq_len(differences)) r <- r[i1, , drop = FALSE] -
            r[-nrow(r):-(nrow(r) - lag + 1), , drop = FALSE]
    else for (i in seq_len(differences))
        r <- r[i1] - r[-length(r):-(length(r) - lag + 1L)]
    r
}

## ---- additions in 2.2.0 -----

duplicated.POSIXlt <- function(x, incomparables = FALSE, ...)
{
    x <- as.POSIXct(x)
    NextMethod("duplicated", x)
}

unique.POSIXlt <- function(x, incomparables = FALSE, ...)
    x[!duplicated(x, incomparables, ...)]

## ---- additions in 2.4.0 -----

sort.POSIXlt <- function(x, decreasing = FALSE, na.last = NA, ...)
    x[order(as.POSIXct(x), na.last = na.last, decreasing = decreasing)]


## ---- additions in 2.6.0 -----

is.numeric.POSIXt <- function(x) FALSE

## ---- additions in 2.8.0 -----

split.POSIXct <-
function(x, f, drop = FALSE, ...)
    lapply(split.default(as.double(x), f, drop = drop), .POSIXct,
           tz = attr(x, "tzone"))

xtfrm.POSIXct <- function(x) as.numeric(x)
xtfrm.POSIXlt <- function(x) as.double(x)  # has POSIXlt method
xtfrm.difftime <- function(x) as.numeric(x)
is.numeric.difftime <- function(x) FALSE


# class generators added in 2.11.0, class order changed in 2.12.0
.POSIXct <- function(xx, tz = NULL)
    structure(xx, class = c("POSIXct", "POSIXt"), tzone = tz)

.POSIXlt <- function(xx, tz = NULL)
    structure(xx, class = c("POSIXlt", "POSIXt"), tzone = tz)

.difftime <- function(xx, units)
    structure(xx, units = units, class = "difftime")

## ---- additions in 2.13.0 -----

names.POSIXlt <-
function(x)
    names(x$year)

`names<-.POSIXlt` <-
function(x, value)
{
    names(x$year) <- value
    x
}
#  File src/library/base/R/dcf.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

read.dcf <-
function(file, fields = NULL, all = FALSE, keep.white = NULL)
{
    if(is.character(file)){
        file <- gzfile(file)
        on.exit(close(file))
    }
    if(!inherits(file, "connection"))
        stop("'file' must be a character string or connection")

    ## For historical reasons, the default is not to accumulate repeated
    ## fields in a record (in fact picking the *last* field occurrence).
    ## Use the internal code for performance reasons, but note that we
    ## could of course as well use
    ##   do.call("cbind",
    ##           lapply(out,
    ##                  function(s)
    ##                  if(is.atomic(s)) s
    ##                  else mapply("[[", s, sapply(s, length))))
    if(!all) return(.Internal(readDCF(file, fields, keep.white)))

    .assemble_things_into_a_data_frame <- function(tags, vals, nums) {
        tf <- factor(tags, levels = unique(tags))

        cnts <- table(nums, tf)
        out <- array(NA_character_, dim = dim(cnts),
                     dimnames = list(NULL, levels(tf)))
        if(all(cnts <= 1L)) {
            ## No repeated tags ...
            out[cbind(nums, tf)] <- vals
            out <- as.data.frame(out, stringsAsFactors = FALSE)
        }
        else {
            levs <- colSums(cnts > 1L) == 0L
            if(any(levs)) {
                inds <- tf %in% levels(tf)[levs]
                out[cbind(nums[inds], tf[inds])] <- vals[inds]
            }
            out <- as.data.frame(out, stringsAsFactors = FALSE)
            for(l in levels(tf)[!levs]) {
                out[[l]] <- rep.int(list(NA_character_), nrow(cnts))
                i <- tf == l
                out[[l]][unique(nums[i])] <- split(vals[i], nums[i])
            }
        }

        out
    }

    on.exit(Sys.setlocale("LC_CTYPE", Sys.getlocale("LC_CTYPE")), add = TRUE)
    Sys.setlocale("LC_CTYPE", "C")

    lines <- readLines(file)

    ## Try to find out about invalid things: mostly, lines which do not
    ## start with blanks but have no ':' ...
    ind <- grep("^[^[:blank:]][^:]*$", lines)
    if(length(ind)) {
        lines <- strtrim(lines[ind], 0.7 * getOption("width"))
        stop(gettextf("Invalid DCF format.\nRegular lines must have a tag.\nOffending lines start with:\n%s",
                      paste0("  ", lines, collapse = "\n")),
             domain = NA)
    }

    line_is_not_empty <- !grepl("^[[:space:]]*$", lines)
    nums <- cumsum(diff(c(FALSE, line_is_not_empty) > 0L) > 0L)
    ## Remove the empty ones so that nums knows which record each line
    ## belongs to.
    nums <- nums[line_is_not_empty]
    lines <- lines[line_is_not_empty]

    ## Deal with escaped blank lines (used by Debian at least for the
    ## Description: values, see man 5 deb-control):
    line_is_escaped_blank <- grepl("^[[:space:]]+\\.[[:space:]]*$", lines)
    if(any(line_is_escaped_blank))
        lines[line_is_escaped_blank] <- ""

    line_has_tag <- grepl("^[^[:blank:]][^:]*:", lines)

    ## Check that records start with tag lines.
    ind <- which(!line_has_tag[which(diff(nums) > 0L) + 1L])
    if(length(ind)) {
        lines <- strtrim(lines[ind], 0.7 * getOption("width"))
        stop(gettextf("Invalid DCF format.\nContinuation lines must not start a record.\nOffending lines start with:\n%s",
                      paste0("  ", lines, collapse = "\n")),
             domain = NA)
    }

    lengths <- rle(cumsum(line_has_tag))$lengths
    ## End positions of field entries.
    pos <- cumsum(lengths)

    tags <- sub(":.*", "", lines[line_has_tag])
    lines[line_has_tag] <-
        sub("[^:]*:[[:space:]]*", "", lines[line_has_tag])
    foldable <- rep.int(is.na(match(tags, keep.white)), lengths)
    lines[foldable] <- sub("^[[:space:]]*", "", lines[foldable])
    lines[foldable] <- sub("[[:space:]]*$", "", lines[foldable])

    vals <- mapply(function(from, to) paste(lines[from:to],
                                            collapse = "\n"),
                   c(1L, pos[-length(pos)] + 1L), pos)

    out <- .assemble_things_into_a_data_frame(tags, vals, nums[pos])

    if(!is.null(fields))
        out <- out[fields]

    out
}

write.dcf <-
function(x, file = "", append = FALSE,
         indent = 0.1 * getOption("width"),
         width = 0.9 * getOption("width"),
         keep.white = NULL)
{
    if(file == "")
        file <- stdout()
    else if(is.character(file)) {
        file <- file(file, ifelse(append, "a", "w"))
        on.exit(close(file))
    }
    if(!inherits(file, "connection"))
        stop("'file' must be a character string or connection")

    ## We need to take care of two things:
    ## * We really should not write out NA entries.
    ## * We have to handle multiple fields per record.

    ## do not assume that the input is valid in this locale
    escape_paragraphs <- function(s)
	gsub("\n \\.([^\n])","\n  .\\1",
	     gsub("\n[ \t]*\n", "\n .\n ", s, perl = TRUE, useBytes = TRUE),
             perl = TRUE, useBytes = TRUE)
    fmt <- function(tag, val, fold = TRUE) {
        s <- if(fold)
            formatDL(rep.int(tag, length(val)), val, style = "list",
                     width = width, indent = indent)
        else {
            ## Need to ensure a leading whitespace for continuation
            ## lines.
            sprintf("%s: %s", tag,
                    gsub("\n([^[:blank:]])", "\n \\1", val))
        }
        escape_paragraphs(s)
    }


    if(!is.data.frame(x))
        x <- as.data.frame(x, stringsAsFactors = FALSE)
    nmx <- names(x)
    out <- matrix("", nrow(x), ncol(x))

    foldable <- is.na(match(nmx, keep.white))

    for(j in seq_along(x)) {
        xj <- x[[j]]
        if(is.atomic(xj)) {
            ## For atomic ("character") columns, things are simple ...
            i <- !is.na(xj)
            out[i, j] <- fmt(nmx[j], xj[i], foldable[j])
        }
        else {
            ## Should be a list ...
            nmxj <- nmx[j]
            fold <- foldable[j]
            i <- !vapply(xj, function(s) (length(s) == 1L) && is.na(s), NA)
            out[i, j] <-
		vapply(xj[i],
                       function(s) {
                           paste(fmt(nmxj, s, fold), collapse = "\n")
                       }, "")
        }
    }
    out <- t(out)
    is_not_empty <- c(out != "")
    eor <- character(sum(is_not_empty))
    if(length(eor)) {
        ## Newline for end of record.
        ## Note that we do not write a trailing blank line.
        eor[ diff(c(col(out))[is_not_empty]) >= 1L ] <- "\n"
    }
    writeLines(paste0(c(out[is_not_empty]), eor), file)
}
#  File src/library/base/R/debug.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

debug <- function(fun, text="", condition=NULL)
    .Internal(debug(fun, text, condition))
debugonce <- function(fun, text="", condition=NULL)
    .Internal(debugonce(fun, text, condition))
undebug <- function(fun) .Internal(undebug(fun))
isdebugged <- function(fun) .Internal(isdebugged(fun))

browserText <- function(n=1L) .Internal(browserText(n))
browserCondition <- function(n=1L) .Internal(browserCondition(n))
browserSetDebug <- function(n=1L) .Internal(browserSetDebug(n))
#  File src/library/base/R/delay.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

delayedAssign <-
    function(x, value, eval.env=parent.frame(1), assign.env=parent.frame(1))
    .Internal(delayedAssign(x, substitute(value), eval.env, assign.env))
#  File src/library/base/R/det.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

## det now uses Lapack and an LU decomposition.  The method argument is
##     no longer used.
## S-plus' Matrix pkg has arg. "logarithm = TRUE" and returns list
##        (which is necessary for keeping the sign when taking log ..)
## S-plus v 6.x has incorporated the Matrix pkg det as determinant

det <- function(x, ...)
{
    z <- determinant(x, logarithm = TRUE, ...)
    c(z$sign * exp(z$modulus))
}

determinant <- function(x, logarithm = TRUE, ...) UseMethod("determinant")

determinant.matrix <- function(x, logarithm = TRUE, ...)
{
    if ((n <- ncol(x)) != nrow(x))
        stop("'x' must be a square matrix")
    if (n < 1L)
	return(structure(list(modulus =
			      structure(if(logarithm) 0 else 1,
					logarithm = logarithm),
			      sign = 1L),
			 class = "det"))
    if (is.complex(x))
        stop("'determinant' not currently defined for complex matrices")
    ## FIXME: should not be so hard to implement; see
    ##      moddet_ge_real() in ../../../modules/lapack/Lapack.c
    ## the 'sign' would have to be complex z, with |z|=1
    .Internal(det_ge_real(x, logarithm))
}
#  File src/library/base/R/diag.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

diag <- function(x = 1, nrow, ncol)
{
    if (is.matrix(x)) {
        if (nargs() > 1L)
            stop("'nrow' or 'ncol' cannot be specified when 'x' is a matrix")

        if((m <- min(dim(x))) == 0L) return(vector(typeof(x), 0L))
        ## NB: need double index to avoid overflows.
        y <- c(x)[1 + 0L:(m - 1L) * (dim(x)[1L] + 1)]
        nms <- dimnames(x)
        if (is.list(nms) && !any(sapply(nms, is.null)) &&
            identical((nm <- nms[[1L]][seq_len(m)]), nms[[2L]][seq_len(m)]))
            names(y) <- nm
        return(y)
    }
    if (is.array(x) && length(dim(x)) != 1L)
        stop("'x' is an array, but not one-dimensional.")

    if (missing(x)) n <- nrow
    else if (length(x) == 1L && nargs() == 1L) {
	n <- as.integer(x)
	x <- 1
    } else n <- length(x)
    if (!missing(nrow)) n <- nrow
    if (missing(ncol)) ncol <- n
    ## some people worry about speed
    .Internal(diag(x, n, ncol))
}

`diag<-` <- function(x, value)
{
    dx <- dim(x)
    if (length(dx) != 2L)
	## no further check, to also work with 'Matrix'
	stop("only matrix diagonals can be replaced")
    len.i <- min(dx)
    len.v <- length(value)
    if (len.v != 1L && len.v != len.i)
	stop("replacement diagonal has wrong length")
    if (len.i) {
	i <- seq_len(len.i)
	x[cbind(i, i)] <- value
    }
    x
}
#  File src/library/base/R/diff.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

diff <- function(x, ...) UseMethod("diff")

diff.default <- function(x, lag = 1L, differences = 1L, ...)
{
    ismat <- is.matrix(x)
    xlen <- if(ismat) dim(x)[1L] else length(x)
    if (length(lag) > 1L || length(differences) > 1L ||
        lag < 1L || differences < 1L)
	stop("'lag' and 'differences' must be integers >= 1")
    if (lag * differences >= xlen)
	return(x[0L]) # empty, but of proper mode
    r <- unclass(x)  # don't want class-specific subset methods
    i1 <- -seq_len(lag)
    if (ismat)
	for (i in seq_len(differences))
	    r <- r[i1, , drop = FALSE] -
                r[-nrow(r):-(nrow(r)-lag+1L), , drop = FALSE]
    else
        for (i in seq_len(differences))
            r <- r[i1] - r[-length(r):-(length(r)-lag+1L)]
    class(r) <- oldClass(x)
    r
}
#  File src/library/base/R/dput.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

dput <-
    function(x, file = "",
             control = c("keepNA", "keepInteger", "showAttributes"))
{
    if(is.character(file))
        if(nzchar(file)) {
            file <- file(file, "wt")
            on.exit(close(file))
        } else file <- stdout()
    opts <- .deparseOpts(control)
    ## FIXME: this should happen in C {deparse2() in ../../../main/deparse.c}
    ##        but we are missing a C-level slotNames()
    ## Fails e.g. if an S3 list-like object has S4 components
    if(isS4(x)) {
        clx <- class(x)
	cat('new("', clx,'"\n', file = file, sep = "")
	for(n in .slotNames(clx)) {
	    cat("    ,", n, "= ", file = file)
	    dput(slot(x, n), file = file, control = control)
	}
	cat(")\n", file = file)
	invisible()
    }
    else .Internal(dput(x, file, opts))
}

dget <- function(file)
    eval(parse(file = file))
#  File src/library/base/R/dump.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

dump <- function (list, file = "dumpdata.R", append = FALSE,
		  control = "all", envir = parent.frame(),
		  evaluate = TRUE)
{
    if(is.character(file)) {
	## avoid opening a file if there is nothing to dump
	ex <- sapply(list, exists, envir=envir)
	if(!any(ex)) return(invisible(character()))
	if(nzchar(file)) {
	    file <- file(file, ifelse(append, "a", "w"))
	    on.exit(close(file), add = TRUE)
	} else file <- stdout()
    }
    opts <- .deparseOpts(control)
    .Internal(dump(list, file, envir, opts, evaluate))
}

#  File src/library/base/R/duplicated.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

duplicated <- function(x, incomparables = FALSE, ...) UseMethod("duplicated")

duplicated.default <-
    function(x, incomparables = FALSE, fromLast = FALSE, nmax = NA, ...)
    .Internal(duplicated(x, incomparables, fromLast,
                         if(is.factor(x)) min(length(x), nlevels(x) + 1L) else nmax))

duplicated.data.frame <-
    function(x, incomparables = FALSE, fromLast = FALSE, ...)
{
    if(!identical(incomparables, FALSE))
	.NotYetUsed("incomparables != FALSE")
    if(length(x) != 1L)
        duplicated(do.call("paste", c(x, sep="\r")), fromLast = fromLast)
    else duplicated(x[[1L]], fromLast = fromLast, ...)
}

duplicated.matrix <- duplicated.array <-
    function(x, incomparables = FALSE, MARGIN = 1L, fromLast = FALSE, ...)
{
    if(!identical(incomparables, FALSE))
	.NotYetUsed("incomparables != FALSE")
    dx <- dim(x)
    ndim <- length(dx)
    if (length(MARGIN) > ndim || any(MARGIN > ndim))
        stop(gettextf("MARGIN = %d is invalid for dim = %d", MARGIN, dx),
             domain = NA)
    collapse <- (ndim > 1L) && (prod(dx[-MARGIN]) > 1L)
    temp <- if(collapse) apply(x, MARGIN, function(x) paste(x, collapse = "\r")) else x
    res <- duplicated.default(temp, fromLast = fromLast, ...)
    dim(res) <- dim(temp)
    dimnames(res) <- dimnames(temp)
    res
}

anyDuplicated <- function(x, incomparables = FALSE, ...)
    UseMethod("anyDuplicated")

anyDuplicated.default <-
    function(x, incomparables = FALSE, fromLast = FALSE, ...)
    .Internal(anyDuplicated(x, incomparables, fromLast))


anyDuplicated.data.frame <-
    function(x, incomparables = FALSE, fromLast = FALSE, ...)
{
    if(!identical(incomparables, FALSE))
	.NotYetUsed("incomparables != FALSE")
    anyDuplicated(do.call("paste", c(x, sep="\r")), fromLast = fromLast)
}

anyDuplicated.matrix <- anyDuplicated.array <-
    function(x, incomparables = FALSE, MARGIN = 1L, fromLast = FALSE, ...)
{
    if(!identical(incomparables, FALSE))
	.NotYetUsed("incomparables != FALSE")
    dx <- dim(x)
    ndim <- length(dx)
    if (length(MARGIN) > ndim || any(MARGIN > ndim))
        stop(gettextf("MARGIN = %d is invalid for dim = %d", MARGIN, dx),
             domain = NA)
    collapse <- (ndim > 1L) && (prod(dx[-MARGIN]) > 1L)
    temp <- if(collapse) apply(x, MARGIN, function(x) paste(x, collapse = "\r")) else x
    anyDuplicated.default(temp, fromLast = fromLast)
}

unique <- function(x, incomparables = FALSE, ...) UseMethod("unique")


## NB unique.default is used by factor to avoid unique.matrix,
## so it needs to handle some other cases.
unique.default <-
    function(x, incomparables = FALSE, fromLast = FALSE, nmax = NA, ...)
{
    if(is.factor(x)) {
        z <- .Internal(unique(x, incomparables, fromLast,
                              min(length(x), nlevels(x) + 1L)))
 	return(factor(z, levels = seq_len(nlevels(x)), labels = levels(x),
               ordered = is.ordered(x)))
    }
    z <- .Internal(unique(x, incomparables, fromLast, nmax))
    if(inherits(x, "POSIXct"))
        structure(z, class = class(x), tzone = attr(x, "tzone"))
    else if(inherits(x, "Date"))
        structure(z, class = class(x))
    else z
}

unique.data.frame <- function(x, incomparables = FALSE, fromLast = FALSE, ...)
{
    if(!identical(incomparables, FALSE))
	.NotYetUsed("incomparables != FALSE")
    x[!duplicated(x, fromLast = fromLast, ...),  , drop = FALSE]
}

unique.matrix <- unique.array <-
    function(x, incomparables = FALSE , MARGIN = 1, fromLast = FALSE, ...)
{
    if(!identical(incomparables, FALSE))
	.NotYetUsed("incomparables != FALSE")
    dx <- dim(x)
    ndim <- length(dx)
    if (length(MARGIN) > ndim || any(MARGIN > ndim))
        stop(gettextf("MARGIN = %d is invalid for dim = %d", MARGIN, dx),
             domain = NA)
    collapse <- (ndim > 1L) && (prod(dx[-MARGIN]) > 1L)
    temp <- if(collapse) apply(x, MARGIN, function(x) paste(x, collapse = "\r")) else x
    args <- rep(alist(a=), ndim)
    names(args) <- NULL
    args[[MARGIN]] <- !duplicated.default(temp, fromLast = fromLast, ...)
    do.call("[", c(list(x), args, list(drop = FALSE)))
}
#  File src/library/base/R/dynload.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

if(.Platform$OS.type == "windows") {
    dyn.load <- function(x, local = TRUE, now = TRUE, ...) {
        inDL <- function(x, local, now, ..., DLLpath = "")
            .Internal(dyn.load(x, local, now, DLLpath))
        inDL(x, as.logical(local), as.logical(now), ...)
    }
} else {
    dyn.load <- function(x, local = TRUE, now = TRUE, ...)
        .Internal(dyn.load(x, as.logical(local), as.logical(now), ""))
}

dyn.unload <- function(x)
    .Internal(dyn.unload(x))

is.loaded <- function(symbol, PACKAGE = "", type = "")
    .Internal(is.loaded(symbol, PACKAGE, type))

getNativeSymbolInfo <- function(name, PACKAGE, unlist = TRUE,
                                 withRegistrationInfo = FALSE)
{
    if(missing(PACKAGE)) PACKAGE <- ""

    if(is.character(PACKAGE))
        pkgName <- PACKAGE
    else if(inherits(PACKAGE, "DLLInfo")) {
        pkgName <- PACKAGE[["path"]]
        PACKAGE <- PACKAGE[["info"]]
    } else if(inherits(PACKAGE, "DLLInfoReference")) {
        pkgName <- character()
    } else
        stop(gettextf("must pass a package name, %s or %s object",
                      dQuote("DLLInfo"),
                      dQuote("DllInfoReference")),
             domain = NA)

    syms <- lapply(name, function(id) {
	v <- .Internal(getSymbolInfo(as.character(id), PACKAGE,
                                     as.logical(withRegistrationInfo)))
	if(is.null(v)) {
	    msg <- paste("no such symbol", id)
	    if(length(pkgName) && nzchar(pkgName))
		msg <- paste(msg, "in package", pkgName)
	    stop(msg, domain = NA)
	}
	names(v) <- c("name", "address", "package", "numParameters")[seq_along(v)]
	v
    })

   if(length(name) == 1L && unlist)
     syms <- syms[[1L]]
   else
     names(syms) <- name

   syms
}

getLoadedDLLs <- function() .Internal(getLoadedDLLs())


getDLLRegisteredRoutines <- function(dll, addNames = TRUE)
    UseMethod("getDLLRegisteredRoutines")


getDLLRegisteredRoutines.character <- function(dll, addNames = TRUE)
{
    dlls <- getLoadedDLLs()
    w <- vapply(dlls, function(x) x[["name"]] == dll || x[["path"]] == dll, NA)

    if(!any(w))
        stop(gettextf("No DLL currently loaded with name or path %s", sQuote(dll)),
             domain = NA)

    dll <- which(w)[1L]
    if(sum(w) > 1L)
        warning(gettextf("multiple DLLs match '%s'. Using '%s'",
                         dll, dll[["path"]]), domain = NA)

    getDLLRegisteredRoutines(dlls[[dll]], addNames)
}


getDLLRegisteredRoutines.DLLInfo <- function(dll, addNames = TRUE)
{
    ## Provide methods for the different types.
    if(!inherits(dll, "DLLInfo"))
        stop(gettextf("must specify DLL via a %s object. See getLoadedDLLs()",
                      dQuote("DLLInfo")),
             domain = NA)

    info <- dll[["info"]]
    els <- .Internal(getRegisteredRoutines(info))
    ## Put names on the elements by getting the names from each element.
    if(addNames) {
      els <- lapply(els, function(x) {
                              if(length(x))
                                 names(x) <- vapply(x, function(z) z$name, "")
                              x
                         })
    }
    class(els) <- "DLLRegisteredRoutines"
    els
}


print.NativeRoutineList <-
function(x, ...)
{
    m <- data.frame(numParameters = sapply(x, function(x) x$numParameters),
                    row.names = sapply(x, function(x) x$name))
    print(m, ...)
    invisible(x)
}

### This is arranged as a ragged data frame.  It may be confusing
### if one reads it row-wise as the columns are related in pairs
### but not across pairs.  We might leave it as  a list of lists
### but that spans a great deal of vertical space and involves
### a lot of scrolling for the user.
print.DLLRegisteredRoutines <-
function(x, ...)
{
    ## Create a data frame with as many rows as the maximum number
    ## of routines in any category. Then fill the column with ""
    ## and then the actual entries.

    n <- vapply(x, length, 1L)
    x <- x[n > 0]
    n <- max(n)
    d <- list()
    sapply(names(x),
             function(id) {
		d[[id]] <<- rep.int("", n)
		names <- vapply(x[[id]], function(x) x$name, "")
                if(length(names)) d[[id]][seq_along(names)] <<- names
                d[[paste(id, "numParameters")]] <<- rep.int("", n)
                names <- sapply(x[[id]], function(x) x$numParameters)
                if(length(names))
                    d[[paste(id, "numParameters")]][seq_along(names)] <<- names
             })
    print(as.data.frame(d), ...)
    invisible(x)
}

getCallingDLLe <- function(e)
{
    if (is.null(env <- e$".__NAMESPACE__.")) env <- baseenv()
    if(exists("DLLs", envir = env) &&
       length(env$DLLs))
        return(env$DLLs[[1L]])
    NULL
}

getCallingDLL <-
function(f = sys.function(-1), doStop = FALSE)
{
    e <- environment(f)

    if(!isNamespace(e)) {
        if(doStop)
            stop("function is not in a namespace, so cannot locate associated DLL")
        else
            return(NULL)
    }

    ## Please feel free to replace with a more encapsulated way to do this.
    if (is.null(env <- e$".__NAMESPACE__.")) env <- baseenv()
    if(exists("DLLs", envir = env) && length(env$DLLs))
        return(env$DLLs[[1L]])
    else {
        if(doStop)
            stop("looking for DLL for native routine call, but no DLLs in namespace of call")
        else
            NULL
    }
    NULL
}

print.DLLInfo <- function(x, ...)
{
    tmp <- as.data.frame.list(x[c("name", "path", "dynamicLookup")])
    names(tmp) <- c("DLL name", "Filename", "Dynamic lookup")
    write.dcf(tmp, ...)
    invisible(x)
}

print.DLLInfoList <- function(x, ...)
{
    if(length(x)) {
        m <- data.frame(Filename = sapply(x, function(x) x[["path"]]),
                        "Dynamic Lookup" =
                        sapply(x, function(x) x[["dynamicLookup"]]))
        print(m, ...)
    }
    invisible(x)
}


`$.DLLInfo` <- function(x, name)
    getNativeSymbolInfo(as.character(name), PACKAGE = x)




#  File src/library/base/R/eapply.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

eapply <- function (env, FUN, ..., all.names = FALSE, USE.NAMES = TRUE)
{
    FUN <- match.fun(FUN)
    .Internal(eapply(env, FUN, all.names, USE.NAMES))
}
#  File src/library/base/R/eigen.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2013 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


isSymmetric <- function(object, ...) UseMethod("isSymmetric")

isSymmetric.matrix <- function(object, tol = 100*.Machine$double.eps, ...)
{
    if(!is.matrix(object)) return(FALSE) ## we test for  symmetric *matrix*
    ## cheap pretest: is it square?
    d <- dim(object)
    if(d[1L] != d[2L]) return(FALSE)
    test <-
        if(is.complex(object))
            all.equal.numeric(object, Conj(t(object)), tolerance = tol, ...)
        else # numeric, character, ..
            all.equal(object, t(object), tolerance = tol, ...)
    isTRUE(test)
}

eigen <- function(x, symmetric, only.values = FALSE, EISPACK = FALSE)
{
    if (EISPACK)
        warning("EISPACK = TRUE is defunct and will be ignored", domain = NA)

    x <- as.matrix(x)
    n <- nrow(x)
    if (!n) stop("0 x 0 matrix")
    if (n != ncol(x)) stop("non-square matrix in 'eigen'")
    n <- as.integer(n)
    if(is.na(n)) stop("invalid nrow(x)")

    complex.x <- is.complex(x)
    if (!all(is.finite(x))) stop("infinite or missing values in 'x'")

    if(missing(symmetric)) symmetric <- isSymmetric.matrix(x)

    if (symmetric) {
        z <- if(!complex.x) .Internal(La_rs(x, only.values))
        else .Internal(La_rs_cmplx(x, only.values))
        ord <- rev(seq_along(z$values))
    } else {
        z <- if(!complex.x) .Internal(La_rg(x, only.values))
        else .Internal(La_rg_cmplx(x, only.values))
        ord <- sort.list(Mod(z$values), decreasing = TRUE)
    }
    return(list(values = z$values[ord],
                vectors = if (!only.values) z$vectors[, ord, drop = FALSE]))
}
#  File src/library/base/R/environment.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

environment <- function(fun=NULL) .Internal(environment(fun))

environmentName <- function(env) .Internal(environmentName(env))

env.profile <- function(env) .Internal(env.profile(env))
#  File src/library/base/R/eval.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

.GlobalEnv <- environment()
parent.frame <- function(n = 1) .Internal(parent.frame(n))

eval <-
    function(expr, envir = parent.frame(),
	     enclos = if(is.list(envir) || is.pairlist(envir))
                       parent.frame() else baseenv())
    .Internal(eval(expr, envir, enclos))

eval.parent <- function(expr, n = 1) {
    p <- parent.frame(n + 1)
    eval(expr, p)
}

evalq <-
    function (expr, envir = parent.frame(), enclos = if (is.list(envir) ||
    is.pairlist(envir)) parent.frame() else baseenv())
      .Internal(eval(substitute(expr), envir, enclos))

new.env <- function (hash = TRUE, parent = parent.frame(), size = 29L)
    .Internal(new.env(hash, parent, size))

parent.env <- function(env)
    .Internal(parent.env(env))

`parent.env<-` <- function(env, value)
    .Internal("parent.env<-"(env, value))

local <-
    function (expr, envir = new.env())
    eval.parent(substitute(eval(quote(expr), envir)))

Recall <- function(...) .Internal(Recall(...))

with <- function(data, expr, ...) UseMethod("with")
within <- function(data, expr, ...) UseMethod("within")

with.default <- function(data, expr, ...)
    eval(substitute(expr), data, enclos=parent.frame())

within.data.frame <- function(data, expr, ...)
{
    parent <- parent.frame()
    e <- evalq(environment(), data, parent)
    eval(substitute(expr), e)
    l <- as.list(e)
    l <- l[!sapply(l, is.null)]
    ## del: variables to *del*ete from data[]
    nD <- length(del <- setdiff(names(data), (nl <- names(l))))
    data[nl] <- l
    if(nD)
	data[del] <- if(nD == 1) NULL else vector("list", nD)
    data
}

within.list <- within.data.frame



force <- function(x) x
#  File src/library/base/R/exists.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

exists <-
    function (x, where = -1,
              envir = if(missing(frame)) as.environment(where) else sys.frame(frame),
              frame, mode = "any", inherits = TRUE)
    .Internal(exists(x, envir, mode, inherits))
#  File src/library/base/R/expand.grid.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

expand.grid <- function(..., KEEP.OUT.ATTRS = TRUE, stringsAsFactors = TRUE)
{
    ## x should either be a list or a set of vectors or factors
    nargs <- length(args <- list(...))
    if(!nargs) return(as.data.frame(list()))
    if(nargs == 1L && is.list(a1 <- args[[1L]]))
	nargs <- length(args <- a1)
    if(nargs == 0L) return(as.data.frame(list()))
    ## avoid classed args such as data frames: cargs <- args
    cargs <- vector("list", nargs)
    iArgs <- seq_len(nargs)
    nmc <- paste0("Var", iArgs)
    nm <- names(args)
    if(is.null(nm))
	nm <- nmc
    else if(any(ng0 <- nzchar(nm)))
	nmc[ng0] <- nm[ng0]
    names(cargs) <- nmc
    rep.fac <- 1L
    d <- sapply(args, length)
    if(KEEP.OUT.ATTRS) {
	dn <- vector("list", nargs)
	names(dn) <- nmc
    }
    orep <- prod(d)
    if(orep == 0L) {
        for(i in iArgs) cargs[[i]] <- args[[i]][FALSE]
    } else {
        for(i in iArgs) {
            x <- args[[i]]
            if(KEEP.OUT.ATTRS)
                dn[[i]] <-
                    paste0(nmc[i], "=", if(is.numeric(x)) format(x) else x)
            nx <- length(x)
            orep <- orep/nx
            x <- x[rep.int(rep.int(seq_len(nx),
                                   rep.int(rep.fac, nx)), orep)]
	    ## avoid sorting the levels of character variates
	    if(stringsAsFactors && !is.factor(x) && is.character(x))
		x <- factor(x, levels = unique(x))
            cargs[[i]] <- x
            rep.fac <- rep.fac * nx
        }
    }
    if(KEEP.OUT.ATTRS)
	attr(cargs, "out.attrs") <- list(dim=d, dimnames=dn)
    rn <- .set_row_names( as.integer(prod(d)) )
    structure(cargs, class = "data.frame", row.names = rn)
}
#  File src/library/base/R/factor.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2013 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

factor <- function(x = character(), levels, labels = levels,
                   exclude = NA, ordered = is.ordered(x), nmax = NA)
{
    if(is.null(x)) x <- character()
    nx <- names(x)
    if (missing(levels)) {
	y <- unique(x, nmax = nmax)
	ind <- sort.list(y) # or possibly order(x) which is more (too ?) tolerant
	y <- as.character(y)
	levels <- unique(y[ind])
    }
    force(ordered) # check if original x is an ordered factor
    exclude <- as.vector(exclude, typeof(x)) # may result in NA
    x <- as.character(x)
    ## levels could be a long vectors, but match will not handle that.
    levels <- levels[is.na(match(levels, exclude))]
    f <- match(x, levels)
    if(!is.null(nx))
	names(f) <- nx
    nl <- length(labels)
    nL <- length(levels)
    if(!any(nl == c(1L, nL)))
	stop(gettextf("invalid 'labels'; length %d should be 1 or %d", nl, nL),
	     domain = NA)
    levels(f) <- ## nl == nL or 1
	if (nl == nL) as.character(labels)
	else paste0(labels, seq_along(levels))
    class(f) <- c(if(ordered) "ordered", "factor")
    f
}

is.factor <- function(x) inherits(x, "factor")
as.factor <- function(x) if (is.factor(x)) x else factor(x)

levels <- function(x) UseMethod("levels")
levels.default <- function(x) attr(x, "levels")
nlevels <- function(x) length(levels(x))

`levels<-.factor` <- function(x, value)
{
    xlevs <- levels(x)
    if (is.list(value)) {
        nlevs <- rep.int(names(value), lapply(value, length))
        value <- unlist(value)
        m <- match(value, xlevs, nomatch = 0L)
        xlevs[m] <- nlevs[m > 0L]
    } else {
        if (length(xlevs) > length(value))
            stop("number of levels differs")
        nlevs <- xlevs <- as.character(value)
        nlevs <- nlevs[!is.na(nlevs)]
    }
    ## take care here not to drop attributes, including class.
    ## factor(xlevs[x], levels = unique(nlevs))
    nlevs <- unique(nlevs)
    at <- attributes(x)
    at$levels <- nlevs
    y <- match(xlevs[x], nlevs)
    attributes(y) <- at
    y
}

droplevels <- function(x, ...) UseMethod("droplevels")
droplevels.factor <- function(x, ...) factor(x)
droplevels.data.frame <- function(x, except = NULL, ...)
  {
    ix <- vapply(x, is.factor, NA)
    if (!is.null(except)) ix[except] <- FALSE
    x[ix] <- lapply(x[ix], factor)
    x
  }

as.vector.factor <- function(x, mode="any")
{
    if(mode=="list") as.list(x)
    else if(mode== "any" || mode== "character" || mode== "logical")
	as.vector(levels(x)[x], mode)
    else
	as.vector(unclass(x), mode)
}

as.character.factor <- function(x,...) levels(x)[x]

as.logical.factor <- function(x,...) as.logical(levels(x))[x]

as.list.factor <- function(x,...)
{
    res <- vector("list", length(x))
    for(i in seq_along(x)) res[[i]] <- x[i]
    res
}

## for `factor' *and* `ordered' :
print.factor <- function (x, quote = FALSE, max.levels = NULL,
                          width = getOption("width"), ...)
{
    ord <- is.ordered(x)
    if (length(x) == 0L)
        cat(if(ord)"ordered" else "factor", "(0)\n", sep = "")
    else {
        ## The idea here is to preserve all relevant attributes such as
        ## names and dims
        xx <- x
        class(xx) <- NULL
        levels(xx) <- NULL
        xx[] <- as.character(x)
        print(xx, quote = quote, ...)
    }
    maxl <- if(is.null(max.levels)) TRUE else max.levels
    if (maxl) {
        n <- length(lev <- encodeString(levels(x), quote=ifelse(quote, '"', '')))
        colsep <- if(ord) " < " else " "
        T0 <- "Levels: "
        if(is.logical(maxl))
            maxl <- { ## smart default
                width <- width - (nchar(T0, "w") + 3L + 1L + 3L)
                                        # 3='...', 3=#lev, 1=extra
                lenl <- cumsum(nchar(lev, "w") + nchar(colsep, "w"))
                if(n <= 1L || lenl[n] <= width) n
		else max(1L, which.max(lenl > width) - 1L)
            }
        drop <- n > maxl
        cat(if(drop) paste(format(n), ""), T0,
            paste(if(drop)c(lev[1L:max(1,maxl-1)],"...",if(maxl > 1) lev[n])
                      else lev, collapse = colsep),
            "\n", sep = "")
    }
    invisible(x)
}


Math.factor <- function(x, ...) {
    stop(.Generic, " not meaningful for factors")
}
## The next two have an .ordered method:
Summary.factor <- function(..., na.rm)
    stop(.Generic, " not meaningful for factors")
Ops.factor <- function(e1, e2)
{
    ok <- switch(.Generic, "=="=, "!="=TRUE, FALSE)
    if(!ok) {
	warning(.Generic, " not meaningful for factors")
	return(rep.int(NA, max(length(e1), if(!missing(e2)) length(e2))))
    }
    nas <- is.na(e1) | is.na(e2)
    ## Need this for NA *levels* as opposed to missing
    noNA.levels <- function(f) {
	r <- levels(f)
	if(any(ina <- is.na(r))) {
	    n <- "  NA "
	    while(n %in% r) n <- paste(n, ".")
	    r[ina] <- n
	}
	r
    }
    if (nzchar(.Method[1L])) { # e1 *is* a factor
	l1 <- noNA.levels(e1)
	e1 <- l1[e1]
    }
    if (nzchar(.Method[2L])) { # e2 *is* a factor
	l2 <- noNA.levels(e2)
	e2 <- l2[e2]
    }
    if (all(nzchar(.Method)) &&
	(length(l1) != length(l2) || !all(sort.int(l2) == sort.int(l1))))
	stop("level sets of factors are different")
    value <- NextMethod(.Generic)
    value[nas] <- NA
    value
}

## NB for next four:
## a factor has levels before class in attribute list (PR#6799)
`[.factor` <- function(x, ..., drop = FALSE)
{
    y <- NextMethod("[")
    attr(y,"contrasts") <- attr(x,"contrasts")
    attr(y,"levels") <- attr(x,"levels")
    class(y) <- oldClass(x)
    lev <- levels(x)
    if (drop)
        factor(y, exclude = if(any(is.na(levels(x)))) NULL else NA ) else y
}

`[<-.factor` <- function(x, ..., value)
{
    lx <- levels(x)
    cx <- oldClass(x)
    if (is.factor(value)) value <- levels(value)[value]
    m <- match(value, lx)
    if (any(is.na(m) & !is.na(value)))
	warning("invalid factor level, NA generated")
    class(x) <- NULL
    x[...] <- m
    attr(x,"levels") <- lx
    class(x) <- cx
    x
}

`[[.factor` <- function(x, ...)
{
    y <- NextMethod("[[")
    attr(y,"contrasts") <- attr(x,"contrasts")
    attr(y,"levels") <- attr(x,"levels")
    class(y) <- oldClass(x)
    y
}

## added for 2.12.0
`[[<-.factor` <- function(x, ..., value)
{
    lx <- levels(x)
    cx <- oldClass(x)
    if (is.factor(value)) value <- levels(value)[value]
    m <- match(value, lx)
    if (any(is.na(m) & !is.na(value)))
	warning("invalid factor level, NA generated")
    class(x) <- NULL
    x[[...]] <- m
    attr(x,"levels") <- lx
    class(x) <- cx
    x
}


## ordered factors ...

ordered <- function(x, ...) factor(x, ..., ordered=TRUE)

is.ordered <- function(x) inherits(x, "ordered")
as.ordered <- function(x) if(is.ordered(x)) x else ordered(x)

Ops.ordered <- function (e1, e2)
{
    ok <- switch(.Generic,
		 "<" = , ">" = , "<=" = , ">=" = ,"=="=, "!=" =TRUE,
		 FALSE)
    if(!ok) {
	warning(sprintf("'%s' is not meaningful for ordered factors",
                        .Generic))
	return(rep.int(NA, max(length(e1), if(!missing(e2)) length(e2))))
    }
    if (.Generic %in% c("==", "!="))
      return(NextMethod(.Generic))  ##not S-PLUS compatible, but saner
    nas <- is.na(e1) | is.na(e2)
    ord1 <- FALSE
    ord2 <- FALSE
    if (nzchar(.Method[1L])) {
	l1 <- levels(e1)
	ord1 <- TRUE
    }
    if (nzchar(.Method[2L])) {
	l2 <- levels(e2)
	ord2 <- TRUE
    }
    if (all(nzchar(.Method)) &&
        (length(l1) != length(l2) || !all(l2 == l1)))
	stop("level sets of factors are different")
    if (ord1 && ord2) {
	e1 <- as.integer(e1) # was codes, but same thing for ordered factor.
	e2 <- as.integer(e2)
    }
    else if (!ord1) {
	e1 <- match(e1, l2)
	e2 <- as.integer(e2)
    }
    else if (!ord2) {
	e2 <- match(e2, l1)
	e1 <- as.integer(e1)
    }
    value <- get(.Generic, mode = "function")(e1, e2)
    value[nas] <- NA
    value
}

Summary.ordered <- function(..., na.rm)
{
    ok <- switch(.Generic, max = , min = , range = TRUE,
		 FALSE)
    if (!ok)
	stop(gettextf("'%s' not defined for ordered factors", .Generic),
	     domain = NA)
    args <- list(...)
    levl <- lapply(args, levels)
    levset <- levl[[1]]
    if (!all(vapply(args, is.ordered, NA)) ||
	!all(sapply(levl, identical, levset)))
	stop(gettextf("'%s' is only meaningful for ordered factors if all arguments have the same level sets",
		      .Generic))
    codes <- lapply(args, as.integer)
    ind <- do.call(.Generic, c(codes, na.rm = na.rm))
    ordered(levset[ind], levels = levset)
}

`is.na<-.factor` <- function(x, value)
{
    lx <- levels(x)
    cx <- oldClass(x)
    class(x) <- NULL
    x[value] <- NA
    structure(x, levels = lx, class = cx)
}

`length<-.factor` <- function(x, value)
{
    cl <- class(x)
    levs <- levels(x)
    x <- NextMethod()
    structure(x, levels=levs, class=cl)
}

addNA <- function(x, ifany=FALSE)
{
    if (!is.factor(x)) x <- factor(x)
    if (ifany & !any(is.na(x))) return(x)
    ll <- levels(x)
    if (!any(is.na(ll))) ll <- c(ll, NA)
    factor(x, levels=ll, exclude=NULL)
}
#  File src/library/base/R/files.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

R.home <- function(component="home")
{
    rh <- .Internal(R.home())
    switch(component,
           "home" = rh,
           "bin" = if(.Platform$OS.type == "windows" &&
                      nzchar(p <- .Platform$r_arch)) file.path(rh, component, p)
           else file.path(rh, component),
           "share" = if(nzchar(p <- Sys.getenv("R_SHARE_DIR"))) p
           else file.path(rh, component),
	   "doc" = if(nzchar(p <- Sys.getenv("R_DOC_DIR"))) p
           else file.path(rh, component),
           "include" = if(nzchar(p <- Sys.getenv("R_INCLUDE_DIR"))) p
           else file.path(rh, component),
           "modules" = if(nzchar(p <- .Platform$r_arch)) file.path(rh, component, p)
           else file.path(rh, component),
           file.path(rh, component))
}

file.show <-
    function (..., header = rep("", nfiles), title = "R Information",
              delete.file = FALSE, pager = getOption("pager"), encoding = "")
{
    files <- path.expand(c(...))
    nfiles <- length(files)
    if(nfiles == 0L)
        return(invisible(NULL))
    ## avoid re-encoding files to the current encoding.
    if(l10n_info()[["UTF-8"]] && encoding == "UTF-8") encoding <- ""
    if(l10n_info()[["Latin-1"]] && encoding == "latin1") encoding <- ""
    if(!is.na(encoding) && encoding != "") {
        for(i in seq_along(files)) {
            f <- files[i]
            tf <- tempfile()
            tmp <- readLines(f, warn = FALSE)
            tmp2 <- try(iconv(tmp, encoding, "", "byte"))
            if(inherits(tmp2, "try-error")) file.copy(f, tf)
            else writeLines(tmp2, tf)
            files[i] <- tf
            if(delete.file) unlink(f)
        }
        delete.file <- TRUE
    }
    if(is.function(pager))
	pager(files, header = header, title = title, delete.file = delete.file)
    else
        .Internal(file.show(files, header, title, delete.file, pager))
}

file.append <- function(file1, file2)
    .Internal(file.append(file1, file2))

file.remove <- function(...)
    .Internal(file.remove(c(...)))

file.rename <- function(from, to)
    .Internal(file.rename(from, to))

list.files <-
    function(path = ".", pattern = NULL, all.files = FALSE,
             full.names = FALSE, recursive = FALSE,
             ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
    .Internal(list.files(path, pattern, all.files, full.names,
			 recursive, ignore.case, include.dirs, no..))

dir <- list.files

list.dirs <- function(path = ".", full.names = TRUE, recursive = TRUE)
    .Internal(list.dirs(path, full.names, recursive))


file.path <-
function(..., fsep=.Platform$file.sep)
    .Internal(file.path(list(...), fsep))


file.exists <- function(...) .Internal(file.exists(c(...)))

file.create <- function(..., showWarnings =  TRUE)
    .Internal(file.create(c(...), showWarnings))

file.choose <- function(new=FALSE) .Internal(file.choose(new))

file.copy <- function(from, to,
                      overwrite = recursive, recursive = FALSE,
                      copy.mode = TRUE)
{
    if (!(nf <- length(from))) return(logical())
    if (!(nt <- length(to)))   stop("no files to copy to")
    ## we don't use file_test as that is in utils.
    if (nt == 1 && isTRUE(file.info(to)$isdir)) {
        if (recursive && to %in% from)
            stop("attempt to copy a directory to itself")
        ## on Windows we need \ for the compiled code (e.g. mkdir).
        if(.Platform$OS.type == "windows") {
            from <- gsub("/", "\\", from, fixed = TRUE)
            to <- gsub("/", "\\", to, fixed = TRUE)
        }
        return(.Internal(file.copy(from, to, overwrite, recursive, copy.mode)))
    } else if (nf > nt) stop("more 'from' files than 'to' files")
    else if (recursive)
        warning("'recursive' will be ignored as 'to' is not a single existing directory")
    if(nt > nf) from <- rep_len(from, length.out = nt)
    okay <- file.exists(from)
    if (!overwrite) okay[file.exists(to)] <- FALSE
    if (any(from[okay] %in% to[okay]))
        stop("file can not be copied both 'from' and 'to'")
    if (any(okay)) { # care: file.create could fail but file.append work.
    	okay[okay] <- file.create(to[okay])
    	if(any(okay)) {
            okay[okay] <- file.append(to[okay], from[okay])
            if(copy.mode)
                Sys.chmod(to[okay], file.info(from[okay])$mode, TRUE)
        }
    }
    okay
}

file.symlink <- function(from, to) {
    if (!(length(from))) stop("no files to link from")
    if (!(nt <- length(to)))   stop("no files/directory to link to")
    if (nt == 1 && file.exists(to) && file.info(to)$isdir)
        to <- file.path(to, basename(from))
    .Internal(file.symlink(from, to))
}

file.link <- function(from, to) {
    if (!(length(from))) stop("no files to link from")
    if (!(nt <- length(to)))   stop("no files to link to")
    .Internal(file.link(from, to))
}

file.info <- function(...)
{
    res <- .Internal(file.info(fn <- c(...)))
    res$mtime <- .POSIXct(res$mtime)
    res$ctime <- .POSIXct(res$ctime)
    res$atime <- .POSIXct(res$atime)
    class(res) <- "data.frame"
    attr(res, "row.names") <- fn # not row.names<- as that does a length check
    res
}

file.access <- function(names, mode = 0)
{
    res <- .Internal(file.access(names, mode))
    names(res) <- names
    res
}

dir.create <- function(path, showWarnings = TRUE, recursive = FALSE,
                       mode = "0777")
    .Internal(dir.create(path, showWarnings, recursive, as.octmode(mode)))

system.file <- function(..., package = "base", lib.loc = NULL, mustWork = FALSE)
{
    if(nargs() == 0L)
        return(file.path(.Library, "base"))
    if(length(package) != 1L)
        stop("'package' must be of length 1")
    packagePath <- find.package(package, lib.loc, quiet = TRUE)
    ans <- if(length(packagePath)) {
        FILES <- file.path(packagePath, ...)
        present <- file.exists(FILES)
        if(any(present)) FILES[present] else ""
    } else ""
    if (mustWork && identical(ans, "")) stop("no file found")
    ans
}

getwd <- function()
    .Internal(getwd())
setwd <- function(dir)
    .Internal(setwd(dir))
basename <- function(path)
    .Internal(basename(path))
dirname <- function(path)
    .Internal(dirname(path))

Sys.info <- function()
    .Internal(Sys.info())

Sys.sleep <- function(time)
    .Internal(Sys.sleep(time))

path.expand <- function(path)
    .Internal(path.expand(path))

Sys.glob <- function(paths, dirmark = FALSE)
    .Internal(Sys.glob(path.expand(paths), dirmark))

unlink <- function(x, recursive = FALSE, force = FALSE)
    .Internal(unlink(as.character(x), recursive, force))

Sys.chmod <- function(paths, mode = "0777", use_umask = TRUE)
    .Internal(Sys.chmod(paths, as.octmode(mode), use_umask))

Sys.umask <- function(mode = NA)
    .Internal(Sys.umask(if(is.na(mode)) NA_integer_ else as.octmode(mode)))

Sys.readlink <- function(paths)
    .Internal(Sys.readlink(paths))

readRenviron <- function(path)
    .Internal(readRenviron(path))

normalizePath <- function(path, winslash = "\\", mustWork = NA)
    .Internal(normalizePath(path.expand(path), winslash, mustWork))

Sys.setFileTime <- function(path, time)
{
    if (!is.character(path) || length(path) != 1L)
        stop("invalid 'path' argument")
    time <- as.POSIXct(time)
    if (is.na(time))  stop("invalid 'time' argument")
    .Internal(setFileTime(path, time))
}
#  File src/library/base/R/findInt.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

### This is a `variant' of  approx( method = "constant" ) :
findInterval <- function(x, vec, rightmost.closed = FALSE, all.inside = FALSE)
{
    ## Purpose: returns the indices of  x in vec;  vec[] sorted
    ## ---------------------------------------------------------
    ## Author: Martin Maechler, Date:  4 Jan 2002, 10:16 (of very different .C version)

    if(any(is.na(vec)))
	stop("'vec' contains NAs")
    if(is.unsorted(vec))
	stop("'vec' must be sorted non-decreasingly")
    .Internal(findInterval(as.double(vec), as.double(x),
                           rightmost.closed, all.inside))
}
#  File src/library/base/R/formals.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

formals <- function(fun = sys.function(sys.parent())) {
    if(is.character(fun))
	fun <- get(fun, mode = "function", envir = parent.frame())
    .Internal(formals(fun))
}

body <- function(fun = sys.function(sys.parent())) {
    if(is.character(fun))
	fun <- get(fun, mode = "function", envir = parent.frame())
    .Internal(body(fun))
}

alist <- function (...) as.list(sys.call())[-1L]

`body<-` <- function (fun, envir = environment(fun), value) {
    if (is.expression(value)) {
        if (length(value) > 1L)
            warning("using the first element of 'value' of type \"expression\"")
        value <- value[[1L]]
    }
    as.function(c(as.list(formals(fun)), list(value)), envir)
}

`formals<-` <- function (fun, envir = environment(fun), value)
{
    bd <- body(fun)
    as.function(c(value,
                  if(is.null(bd) || is.list(bd)) list(bd) else bd),
                envir)
}
#  File src/library/base/R/format.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

format <- function(x, ...) UseMethod("format")

format.default <-
    function(x, trim = FALSE, digits = NULL, nsmall = 0L,
	     justify = c("left", "right", "centre", "none"),
	     width = NULL, na.encode = TRUE, scientific = NA,
	     big.mark = "", big.interval = 3L,
	     small.mark = "", small.interval = 5L, decimal.mark = ".",
	     zero.print = NULL, drop0trailing = FALSE, ...)
{
    justify <- match.arg(justify)
    adj <- match(justify, c("left", "right", "centre", "none")) - 1L
    if(is.list(x)) {
	## do it this way to force evaluation of args
	if(missing(trim)) trim <- TRUE
	if(missing(justify)) justify <- "none"
	res <- lapply(X = x,
                      FUN = function(xx, ...) format.default(unlist(xx),...),
		      trim = trim, digits = digits, nsmall = nsmall,
		      justify = justify, width = width, na.encode = na.encode,
		      scientific = scientific,
		      big.mark = big.mark, big.interval = big.interval,
		      small.mark = small.mark, small.interval = small.interval,
		      decimal.mark = decimal.mark, zero.print = zero.print,
		      drop0trailing = drop0trailing, ...)
	sapply(res, paste, collapse = ", ")
    } else {
	switch(mode(x),
	       NULL = "NULL",
	       character = .Internal(format(x, trim, digits, nsmall, width,
					    adj, na.encode, scientific)),
	       call=, expression=, "function"=, "(" = deparse(x),
	       raw = as.character(x),
           {
	       ## else: logical, numeric, complex, .. :
	       prettyNum(.Internal(format(x, trim, digits, nsmall, width,
					  3L, na.encode, scientific)),
			 big.mark = big.mark, big.interval = big.interval,
			 small.mark = small.mark,
			 small.interval = small.interval,
			 decimal.mark = decimal.mark,
			 zero.print = zero.print, drop0trailing = drop0trailing,
			 is.cmplx = is.complex(x),
			 preserve.width = if (trim) "individual" else "common")
           })
    }
}

format.pval <- function(pv, digits = max(1L, getOption("digits") - 2L),
			eps = .Machine$double.eps, na.form = "NA", ...)
{
    ## Format  P values; auxiliary for print.summary.[g]lm(.)

    if((has.na <- any(ina <- is.na(pv)))) pv <- pv[!ina]
    ## Better than '0.0' for very small values `is0':
    r <- character(length(is0 <- pv < eps))
    if(any(!is0)) {
	rr <- pv <- pv[!is0]
	## be smart -- differ for fixp. and expon. display:
	expo <- floor(log10(ifelse(pv > 0, pv, 1e-50)))
	fixp <- expo >= -3 | (expo == -4 & digits>1)
	if(any( fixp)) rr[ fixp] <- format(pv[ fixp], digits = digits, ...)
	if(any(!fixp)) rr[!fixp] <- format(pv[!fixp], digits = digits, ...)
	r[!is0] <- rr
    }
    if(any(is0)) {
	digits <- max(1L, digits - 2L)
	if(any(!is0)) {
	    nc <- max(nchar(rr, type="w"))
	    if(digits > 1L && digits + 6L > nc)
		digits <- max(1L, nc - 7L)
	    sep <- if(digits == 1L && nc <= 6L) "" else " "
	} else sep <- if(digits == 1) "" else " "
	r[is0] <- paste("<", format(eps, digits = digits, ...), sep = sep)
    }
    if(has.na) { ## rarely
	rok <- r
	r <- character(length(ina))
	r[!ina] <- rok
	r[ina] <- na.form
    }
    r
}

## Martin Maechler <maechler@stat.math.ethz.ch> , 1994-1998,
## many corrections by R-core.
formatC <- function (x, digits = NULL, width = NULL,
		     format = NULL, flag = "", mode = NULL,
		     big.mark = "", big.interval = 3L,
		     small.mark = "", small.interval = 5L,
		     decimal.mark = ".", preserve.width = "individual",
                     zero.print = NULL, drop0trailing = FALSE)
{
    if(is.object(x)) {
        x <- unclass(x)
        warning("class of 'x' was discarded")
    }

    format.char <- function (x, width, flag)
    {
	if(is.null(width)) width <- 0L
	else if(width < 0L) { flag <- "-"; width <- -width }
	format.default(x, width=width,
		       justify = if(flag=="-") "left" else "right")
    }
     blank.chars <- function(no)
 	vapply(no+1L, function(n) paste(character(n), collapse=" "), "")

    if (!(n <- length(x))) return("")
    if (is.null(mode))	  mode <- storage.mode(x)
    else if (any(mode == c("double", "real", "integer")))  {
      ## for .C call later on
	if(mode=="real") mode <- "double"
	storage.mode(x) <- mode
    }
    else if (mode != "character") stop("'mode' must be \"double\" (\"real\"), \"integer\" or \"character\"")
    if (mode == "character" || (!is.null(format) && format == "s")) {
	if (mode != "character") {
	    warning('coercing argument to "character" for format="s"')
	    x <- as.character(x)
	}
	return(format.char(x, width=width, flag=flag))
    }
    if (missing(format) || is.null(format))
	format <- if (mode == "integer") "d" else "g"
    else {
	if (any(format == c("f", "e", "E", "g", "G", "fg"))) {
	    if (mode == "integer") mode <- storage.mode(x) <- "double"
	}
	else if (format == "d") {
	    if (mode != "integer") mode <- storage.mode(x) <- "integer"
	}
	else stop('\'format\' must be one of {"f","e","E","g","G", "fg", "s"}')
    }
    some.special <- !all(Ok <- is.finite(x))
    if (some.special) {
	rQ <- as.character(x[!Ok])
	rQ[is.na(rQ)] <- "NA"
	x[!Ok] <- as.vector(0, mode = mode)
    }
    if(is.null(width) && is.null(digits))
	width <- 1L
    if (is.null(digits))
	digits <- if (mode == "integer") 2L else 4L
    else if(digits < 0L)
	digits <- 6L
    else {
	maxDigits <- if(format != "f") 50L else ceiling(-(.Machine$double.neg.ulp.digits + .Machine$double.min.exp) / log2(10))
	if (digits > maxDigits) {
            warning(gettextf("'digits' reduced to %d", maxDigits),
                    domain = NA)
	    digits <- maxDigits
	}
    }
    if(is.null(width))	width <- digits + 1L
    else if (width == 0L) width <- digits
    i.strlen <-
	pmax(abs(as.integer(width)),
	     if(format == "fg" || format == "f") {
		 xEx <- as.integer(floor(log10(abs(x+ifelse(x==0,1,0)))))
		 as.integer(x < 0 | flag!="") + digits +
		     if(format == "f") {
			 2L + pmax(xEx, 0L)
		     } else {# format == "fg"
			 pmax(xEx, digits, digits + (-xEx) + 1L) +
			     ifelse(flag != "", nchar(flag, "b"), 0L) + 1L
		     }
	     } else # format == "g" or "e":
	     rep.int(digits + 8L, n)
	     )
    ## sanity check for flags added 2.1.0
    flag <- as.character(flag)
    nf <- strsplit(flag, "")[[1L]]
    if(!all(nf %in% c("0", "+", "-", " ", "#")))
	stop("'flag' can contain only '0+- #'")
    if(digits > 0 && any(nf == "#"))
	digits <- -digits # C-code will notice "do not drop trailing zeros"

    attr(x, "Csingle") <- NULL	# avoid interpreting as.single
    r <- .Internal(formatC(x, as.character(mode), width, digits,
                           as.character(format), as.character(flag),
                           i.strlen))
    if (some.special) r[!Ok] <- format.char(rQ, width = width, flag = flag)

    if(big.mark != "" || small.mark != "" || decimal.mark != "." ||
       !is.null(zero.print) || drop0trailing)
	r <- prettyNum(r, big.mark = big.mark, big.interval = big.interval,
		       small.mark = small.mark, small.interval = small.interval,
		       decimal.mark = decimal.mark, preserve.width = preserve.width,
		       zero.print = zero.print, drop0trailing = drop0trailing,
		       is.cmplx = FALSE)

    if (!is.null(x.atr <- attributes(x)))
	attributes(r) <- x.atr
    r
}

format.factor <- function (x, ...)
    format(structure(as.character(x), names=names(x),
                     dim=dim(x), dimnames=dimnames(x)), ...)


format.data.frame <- function(x, ..., justify = "none")
{
    nr <- .row_names_info(x, 2L)
    nc <- length(x)
    rval <- vector("list", nc)
    for(i in 1L:nc)
	rval[[i]] <- format(x[[i]], ..., justify = justify)
    lens <- sapply(rval, NROW)
    if(any(lens != nr)) { # corrupt data frame, must have at least one column
	warning("corrupt data frame: columns will be truncated or padded with NAs")
	for(i in 1L:nc) {
	    len <- NROW(rval[[i]])
	    if(len == nr) next
	    if(length(dim(rval[[i]])) == 2L) {
		rval[[i]] <- if(len < nr)
		    rbind(rval[[i]], matrix(NA, nr-len, ncol(rval[[i]])))
		else rval[[i]][1L:nr,]
	    } else {
		rval[[i]] <- if(len < nr) c(rval[[i]], rep.int(NA, nr-len))
		else rval[[i]][1L:nr]
	    }
	}
    }
    for(i in 1L:nc) {
	if(is.character(rval[[i]]) && inherits(rval[[i]], "character"))
	    oldClass(rval[[i]]) <- "AsIs"
    }
    cn <- names(x)
    m <- match(c("row.names", "check.rows", "check.names", ""), cn, 0L)
    if(any(m)) cn[m] <- paste0("..dfd.", cn[m])
    ## This requires valid symbols for the columns, so we need to
    ## truncate any of more than 256 bytes.
    long <- nchar(cn, "bytes") > 256L
    cn[long] <- paste(substr(cn[long], 1L, 250L), "...")
    names(rval) <- cn
    rval$check.names <- FALSE
    rval$row.names <- row.names(x)
    x <- do.call("data.frame", rval)
    ## x will have more cols than rval if there are matrix/data.frame cols
    if(any(m)) names(x) <- sub("^..dfd.", "", names(x))
    x
}

format.AsIs <- function(x, width = 12, ...)
{
    if(is.character(x)) return(format.default(x, ...))
    if(is.null(width)) width = 12L
    n <- length(x)
    rvec <- rep.int(NA_character_, n)
    for(i in 1L:n) {
        y <- x[[i]]
        ## need to remove class AsIs to avoid an infinite loop.
        cl <- oldClass(y)
        if(m <- match("AsIs", cl, 0L)) oldClass(y) <- cl[-m]
        rvec[i] <- toString(y, width = width, ...)
    }
    ## AsIs might be around a matrix, which is not a class.
    dim(rvec) <- dim(x)
    dimnames(rvec) <- dimnames(x)
    format.default(rvec, justify = "right")
}

prettyNum <-
    function(x,
	     big.mark = "", big.interval = 3L,
	     small.mark = "", small.interval = 5L,
	     decimal.mark = ".",
	     preserve.width = c("common", "individual", "none"),
	     zero.print = NULL, drop0trailing = FALSE, is.cmplx = NA, ...)
{
    if(!is.character(x)) {
        is.cmplx <- is.complex(x)
	x <- sapply(X = x, FUN = format, ...)
    }
    ## be fast in trivial case (when all options have their default):
    nMark <- big.mark== "" && small.mark== "" && decimal.mark== "."
    nZero <- is.null(zero.print) && !drop0trailing
    if(nMark && nZero)
	return(x)

    ## else
    if(!is.null(zero.print) && any(i0 <- as.numeric(x) == 0)) {
	## print zeros according to 'zero.print' (logical or string):
	if(length(zero.print) > 1L) stop("'zero.print' has length > 1")
	if(is.logical(zero.print))
	    zero.print <- if(zero.print) "0" else " "
	if(!is.character(zero.print))
	    stop("'zero.print' must be character, logical or NULL")
	blank.chars <- function(no) # as in formatC()
	    vapply(no+1L, function(n) paste(character(n), collapse=" "), "")
	nz <- nchar(zero.print, "c")
	nc <- nchar(x[i0], "c")
	ind0 <- regexpr("0", x[i0], fixed = TRUE)# first '0' in string
	substr(x[i0],ind0, (i1 <- ind0+nz-1L)) <- zero.print
	substr(x[i0],ind0+nz, nc) <- blank.chars(nc - i1)
    }
    if(nMark && !drop0trailing)# zero.print was only non-default
	return(x)

    ## else
    if(is.na(is.cmplx)) { ## find if 'x' is format from a *complex*
	ina <- is.na(x) | x == "NA"
	is.cmplx <-
	    if(all(ina)) FALSE
	    else length(grep("[0-9].*[-+][0-9].*i$", x)) > 0
    }
    if(is.cmplx) {
	## should be rare .. taking an easy route
	z.sp <- strsplit(sub("([0-9] *)([-+])( *[0-9])",
			     "\\1::\\2::\\3", x), "::", fixed=TRUE)
	## be careful, if x had an  "	NA":
	i3 <- vapply(z.sp, length, 0L) == 3L # those are re + im *i
	if(any(i3)) {
	    z.sp <- z.sp[i3]
	    z.im <- sapply(z.sp, `[[`, 3L)
	    ## drop ending 'i' (and later re-add it)
	    has.i <- grep("i$", z.im)
	    z.im[has.i] <- sub("i$", '', z.im[has.i])
	    r <- lapply(list(sapply(z.sp, `[[`, 1L), z.im),
			function(.)
			prettyNum(.,
				  big.mark=big.mark, big.interval=big.interval,
				  small.mark=small.mark, small.interval=small.interval,
				  decimal.mark=decimal.mark, preserve.width=preserve.width,
				  zero.print=zero.print, drop0trailing=drop0trailing,
				  is.cmplx=FALSE, ...))
	    r[[2]][has.i] <- paste0(r[[2]][has.i], "i")
	    x[i3] <- paste0(r[[1]], sapply(z.sp, `[[`, 2L), r[[2]])
	}
	return(x)
    }
    preserve.width <- match.arg(preserve.width)
    x.sp <- strsplit(x, ".", fixed=TRUE)
    revStr <- function(cc)
	sapply(lapply(strsplit(cc,NULL), rev), paste, collapse="")
    B. <- sapply(x.sp, `[`, 1L)	    # Before "."
    A. <- sapply(x.sp, `[`, 2)	    # After  "." ; empty == NA
    if(any(iN <- is.na(A.))) A.[iN] <- ""

    if(nzchar(big.mark) &&
       length(i.big <- grep(paste0("[0-9]{", big.interval + 1L,",}"), B.))
       ) { ## add 'big.mark' in decimals before "." :
	B.[i.big] <-
	    revStr(gsub(paste0("([0-9]{",big.interval,"})\\B"),
			paste0("\\1",revStr(big.mark)), revStr(B.[i.big])))
    }
    if(nzchar(small.mark) &&
       length(i.sml <- grep(paste0("[0-9]{", small.interval + 1L,",}"), A.))
       ) { ## add 'small.mark' in decimals after "."  -- but *not* trailing
	A.[i.sml] <- gsub(paste0("([0-9]{",small.interval,"}\\B)"),
			  paste0("\\1",small.mark), A.[i.sml])
    }
    if(drop0trailing) {
	a <- A.[!iN]
	if(length(hasE <- grep("e", a, fixed=TRUE))) {
	    a[ hasE] <- sub("e[+-]0+$", '', a[ hasE]) # also drop "e+00"
	    a[-hasE] <- sub("0+$",	'', a[-hasE])
	} else a <- sub("0+$", '', a)
	A.[!iN] <- a
	## iN := TRUE for those A.[]  which are ""
	iN <- !nzchar(A.)
    }
    ## extraneous trailing dec.marks: paste(B., A., sep = decimal.mark)
    A. <- paste0(B., c(decimal.mark, "")[iN+ 1L], A.)
    if(preserve.width != "none") {
	nnc <- nchar(A., "c")
	d.len <- nnc - nchar(x, "c") # extra space added by 'marks' above
	if(any(ii <- d.len > 0L)) {
	    switch(preserve.width,
		   "individual" = {
		       ## drop initial blanks preserving original width
		       ## where possible:
		       A.[ii] <- sapply(which(ii), function(i)
					sub(sprintf("^ {1,%d}", d.len[i]), "",
					    A.[i]))
		   },
		   "common" = {
		       A. <- format(A., justify = "right")
		   })
	}
    }
    attributes(A.) <- attributes(x)
    class(A.) <- NULL
    A.
}
#  File src/library/base/R/frametools.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2013 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

subset.data.frame <- function (x, subset, select, drop = FALSE, ...)
{
    r <- if(missing(subset))
	rep_len(TRUE, nrow(x))
    else {
	e <- substitute(subset)
	r <- eval(e, x, parent.frame())
        if(!is.logical(r)) stop("'subset' must be logical")
	r & !is.na(r)
    }
    vars <- if(missing(select))
	TRUE
    else {
	nl <- as.list(seq_along(x))
	names(nl) <- names(x)
	eval(substitute(select), nl, parent.frame())
    }
    x[r, vars, drop = drop]
}

subset <- function(x, ...) UseMethod("subset")

subset.default <- function(x, subset, ...) {
    if(!is.logical(subset)) stop("'subset' must be logical")
    x[subset & !is.na(subset)]
}

subset.matrix <- function(x, subset, select, drop = FALSE, ...)
{
    if(missing(select))
	vars <- TRUE
    else {
	nl <- as.list(1L:ncol(x))
	names(nl) <- colnames(x)
	vars <- eval(substitute(select), nl, parent.frame())
    }
    if(missing(subset)) subset <- TRUE
    else if(!is.logical(subset)) stop("'subset' must be logical")
    x[subset & !is.na(subset), vars, drop = drop]
}

### Notice use of non-syntactic variable name for the first argument
### This used to be "x", but then you couldn't create a variable
### called "x"...

transform.data.frame <- function (`_data`, ...)
{
    e <- eval(substitute(list(...)), `_data`, parent.frame())
    tags <- names(e)
    inx <- match(tags, names(`_data`))
    matched <- !is.na(inx)
    if (any(matched)) {
	`_data`[inx[matched]] <- e[matched]
	`_data` <- data.frame(`_data`)
    }
    if (!all(matched))  # add as separate arguments to get replication
	do.call("data.frame", c(list(`_data`), e[!matched]))
    else `_data`
}

transform <- function(`_data`,...) UseMethod("transform")

## Actually, I have no idea what to transform(), except dataframes.
## The default converts its argument to a dataframe and transforms
## that. This is probably marginally useful at best. --pd
transform.default <- function(`_data`,...)
    transform.data.frame(data.frame(`_data`),...)
#  File src/library/base/R/funprog.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

Reduce <-
function(f, x, init, right = FALSE, accumulate = FALSE)
{
    mis <- missing(init)
    len <- length(x)

    if(len == 0L) return(if(mis) NULL else init)

    f <- match.fun(f)

    ## Try to avoid the "obvious"
    ##   if(!mis) x <- if(right) c(x, init) else c(init, x)
    ## to be more efficient ...

    if(!is.vector(x) || is.object(x))
        x <- as.list(x)

    ind <- seq_len(len)

    if(mis) {
        if(right) {
            init <- x[[len]]
            ind <- ind[-len]
        }
        else {
            init <- x[[1L]]
            ind <- ind[-1L]
        }
    }

    if(!accumulate) {
        if(right) {
            for(i in rev(ind))
                init <- f(x[[i]], init)
        }
        else {
            for(i in ind)
                init <- f(init, x[[i]])
        }
        init
    }
    else {
        len <- length(ind) + 1L
        ## We need a list to accumulate the results as these do not
        ## necessarily all have length one (e.g., reducing with c()).
        out <- vector("list", len)
        if(mis) {
            if(right) {
                out[[len]] <- init
                for(i in rev(ind)) {
                    init <- f(x[[i]], init)
                    out[[i]] <- init
                }
            } else {
                out[[1L]] <- init
                for(i in ind) {
                    init <- f(init, x[[i]])
                    out[[i]] <- init
                }
            }
        } else {
            if(right) {
                out[[len]] <- init
                for(i in rev(ind)) {
                    init <- f(x[[i]], init)
                    out[[i]] <- init
                }
            }
            else {
                for(i in ind) {
                    out[[i]] <- init
                    init <- f(init, x[[i]])
                }
                out[[len]] <- init
            }
        }
        ## If all results have length one, we can simplify.
        ## (Note that we do not simplify to arrays in case all results
        ## have a common length > 1.)
	if(all(vapply(out, length, 1.) == 1L))
            out <- unlist(out, recursive = FALSE)
        out
    }
}

Filter <-
function(f, x)
{
    ind <- as.logical(unlist(lapply(x, f)))
    x[!is.na(ind) & ind]
}


Map <-
function(f, ...)
{
    f <- match.fun(f)
    mapply(FUN = f, ..., SIMPLIFY = FALSE)
}

Negate <-
function(f)
{
    f <- match.fun(f) # effectively force f, avoid lazy eval.
    function(...) ! f(...)
}

Position <-
function(f, x, right = FALSE, nomatch = NA_integer_)
{
    ind <- if(right) rev(seq_along(x)) else seq_along(x)

    for(i in ind)
        if(f(x[[i]]))
            return(i)

    nomatch
}

Find <-
function(f, x, right = FALSE, nomatch = NULL)
{
    f <- match.fun(f)
    if((pos <- Position(f, x, right, nomatch = 0L)) > 0L)
        x[[pos]]
    else
        nomatch
}

identity <-
function(x)
    x
#  File src/library/base/R/get.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

get <-
    function (x, pos = -1L, envir = as.environment(pos), mode = "any",
              inherits = TRUE)
    .Internal(get(x, envir, mode, inherits))

mget <- function(x, envir = as.environment(-1L), mode = "any",
                 ifnotfound, inherits = FALSE)
    .Internal(mget(x, envir, mode,
                   if(missing(ifnotfound))
                   list(function(x) stop(gettextf("value for %s not found", sQuote(x)), call. = FALSE)) else ifnotfound,
                   inherits))

## DB's proposed name "getSlotOrComponent" is more precise but harder to type
getElement <- function(object, name) {
    if(isS4(object)) slot(object, name) else object[[name, exact=TRUE]]
}
#  File src/library/base/R/getenv.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

Sys.getenv <- function(x = NULL, unset = "", names = NA)
{
    if (is.null(x)) {
        ## This presumes that '=' does not appear as part of the name
        ## of an environment variable.  That used to happen on Windows.
	x <- strsplit(.Internal(Sys.getenv(character(), "")), "=", fixed=TRUE)
	v <- n <- character(LEN <- length(x))
	for (i in 1L:LEN) {
	    n[i] <- x[[i]][1L]
	    v[i] <- paste(x[[i]][-1L], collapse = "=")
	}
        if (!identical(names, FALSE)) v <- structure(v, names = n)
	v[sort.list(n)]
    } else {
        v <- .Internal(Sys.getenv(as.character(x), as.character(unset)))
	if (isTRUE(names) || (length(x) > 1L && !identical(names, FALSE)))
            structure(v, names = x)
        else v
    }
}

Sys.setenv <- function(...)
{
    x <- list(...)
    nm <- names(x)
    if(is.null(nm) || "" %in% nm)
        stop("all arguments must be named")
    .Internal(Sys.setenv(nm, as.character(unlist(x))))
}

Sys.unsetenv <- function(x) .Internal(Sys.unsetenv(as.character(x)))

Sys.getpid <- function() .Internal(Sys.getpid())
#  File src/library/base/R/gl.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

## gl function of GLIM

gl <- function (n, k, length = n*k, labels=1:n, ordered=FALSE)
  {
    ## We avoid calling factor(), for efficiency.

    ## Must set levels before class.
    ## That way, `levels<-` will pick up an invalid
    ## labels specification.

    f <- rep_len(rep.int(1:n, rep.int(k,n)), length)
    levels(f) <- as.character(labels)
    class(f) <- c(if (ordered) "ordered", "factor")
    f
  }
#  File src/library/base/R/grep.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

strsplit <-
function(x, split, fixed = FALSE, perl = FALSE, useBytes = FALSE)
    .Internal(strsplit(x, as.character(split), fixed, perl, useBytes))

grep <-
function(pattern, x, ignore.case = FALSE, perl = FALSE,
         value = FALSE, fixed = FALSE, useBytes = FALSE, invert = FALSE)
{
    ## when value = TRUE we return names
    if(!is.character(x)) x <- structure(as.character(x), names=names(x))
    .Internal(grep(as.character(pattern), x, ignore.case, value,
                   perl, fixed, useBytes, invert))
}

grepl <-
function(pattern, x, ignore.case = FALSE, perl = FALSE,
         fixed = FALSE, useBytes = FALSE)
{
    if(!is.character(x)) x <- as.character(x)
    .Internal(grepl(as.character(pattern), x, ignore.case, FALSE,
                    perl, fixed, useBytes, FALSE))
}

sub <-
function(pattern, replacement, x, ignore.case = FALSE,
         perl = FALSE, fixed = FALSE, useBytes = FALSE)
{
    if (!is.character(x)) x <- as.character(x)
     .Internal(sub(as.character(pattern), as.character(replacement), x,
                  ignore.case, perl, fixed, useBytes))
}

gsub <-
function(pattern, replacement, x, ignore.case = FALSE,
         perl = FALSE, fixed = FALSE, useBytes = FALSE)
{
    if (!is.character(x)) x <- as.character(x)
    .Internal(gsub(as.character(pattern), as.character(replacement), x,
                   ignore.case, perl, fixed, useBytes))
}

regexpr <-
function(pattern, text, ignore.case = FALSE, perl = FALSE,
         fixed = FALSE, useBytes = FALSE)
{
    if (!is.character(text)) text <- as.character(text)
    .Internal(regexpr(as.character(pattern), text,
                      ignore.case, perl, fixed, useBytes))
}

gregexpr <-
function(pattern, text, ignore.case = FALSE, perl = FALSE,
         fixed = FALSE, useBytes = FALSE)
{
    if (!is.character(text)) text <- as.character(text)
    .Internal(gregexpr(as.character(pattern), text,
                       ignore.case, perl, fixed, useBytes))
}

grepRaw <-
function(pattern, x, offset = 1L, ignore.case = FALSE, value = FALSE,
         fixed = FALSE, all = FALSE, invert = FALSE)
{
    if (!is.raw(pattern)) pattern <- charToRaw(as.character(pattern))
    if (!is.raw(x)) x <- charToRaw(as.character(x))
    .Internal(grepRaw(pattern, x, offset, ignore.case, fixed, value, all, invert))
}

regexec <-
function(pattern, text, ignore.case = FALSE,
         fixed = FALSE, useBytes = FALSE)
    .Internal(regexec(pattern, text, ignore.case, fixed, useBytes))

agrep <-
function(pattern, x, max.distance = 0.1, costs = NULL,
         ignore.case = FALSE, value = FALSE, fixed = TRUE,
         useBytes = FALSE)
{
    pattern <- as.character(pattern)
    if(!is.character(x))
        x <- as.character(x)

    ## TRE needs integer costs: coerce here for simplicity.
    costs <- as.integer(.amatch_costs(costs))
    bounds <- .amatch_bounds(max.distance)

    .Internal(agrep(pattern, x, ignore.case, value, costs, bounds,
                    useBytes, fixed))
}

.amatch_bounds <-
function(x = 0.1)
{
    ## Expand max match distance argument for agrep() et al into bounds
    ## for the TRE regaparams struct.

    ## Note that TRE allows for possibly different (integer) costs for
    ## insertions, deletions and substitions, and allows for specifying
    ## separate bounds for these numbers as well as the total number of
    ## "errors" (transformations) and the total cost.
    ##
    ## When using unit costs (and older versions of agrep() did not
    ## allow otherwise), the total number of errors is the same as the
    ## total cost, and bounds on the total number of errors imply the
    ## same bounds for the individual transformation counts.  This no
    ## longer holds when using possibly different costs.
    ##
    ## See ? agrep for details on handling the match distance argument.
    ##
    ## Older versions of agrep() expanded fractions (of the pattern
    ## length) in R code: but as the C code determines whether matching
    ## used bytes or characters, only the C code can determine the
    ## pattern length and hence expand fractions.
    ##
    ## Unspecified bounds are taken as NA_real_, and set to INT_MAX by
    ## the C code.

    if(!is.list(x)) {
        ## Sanity checks.
        if(!is.numeric(x) || (x < 0))
            stop("match distance components must be non-negative")
        bounds <- c(as.double(x), rep.int(NA_real_, 4L))
    } else {
        table <-
            c("cost", "insertions", "deletions", "substitutions", "all")
        ## Partial matching.
        pos <- pmatch(names(x), table)
        if(any(is.na(pos))) {
            warning("unknown match distance components ignored")
            x <- x[!is.na(pos)]
        }
        names(x) <- table[pos]
        ## Sanity checks.
        x <- unlist(x)
        if(!all(is.numeric(x)) || any(x < 0))
            stop("match distance components must be non-negative")
        ## Defaults.
        if(!is.na(x["cost"])) {
            bounds <- rep.int(NA_real_, 5L)
        } else {
            ## If 'cost' is missing: if 'all' is missing it is set to
            ## 0.1, and the other transformation number bounds default
            ## to 'all'.
            if(is.na(x["all"]))
                x["all"] <- 0.1
            bounds <- c(NA_real_, rep.int(x["all"], 4L))
        }
        names(bounds) <- table
        bounds[names(x)] <- x
    }

    bounds
}

.amatch_costs <-
function(x = NULL)
{
    costs <- c(insertions = 1, deletions = 1, substitutions = 1)
    if(!is.null(x)) {
        x <- as.list(x)
        ## Partial matching.
        pos <- pmatch(names(x), names(costs))
        if(any(is.na(pos))) {
            warning("unknown cost components ignored")
            x <- x[!is.na(pos)]
        }
        ## Sanity checks.
        x <- unlist(x)
        if(!all(is.numeric(x)) || any(x < 0))
            stop("cost components must be non-negative")
        costs[pos] <- x
    }

    costs
}

regmatches <-
function(x, m, invert = FALSE)
{
    if(length(x) != length(m))
        stop(gettextf("%s and %s must have the same length",
                      sQuote("x"), sQuote("m")),
             domain = NA)

    ili <- is.list(m)

    ## Handle useBytes/encoding issues.
    ## For regexpr() and gregexpr(), we get useBytes as TRUE if useBytes
    ## was given as TRUE, or all character string involved were ASCII.
    ## Hence, if useBytes is TRUE, we need to convert non-ASCII strings
    ## to a bytes encoding before computing match substrings.
    useBytes <- if(ili)
        any(unlist(lapply(m, attr, "useBytes")))
    else
        any(attr(m, "useBytes"))
    if(useBytes) {
        ## Cf. tools::showNonASCII():
        asc <- iconv(x, "latin1", "ASCII")
        ind <- is.na(asc) | (asc != x)
        ## Alternatively, could do as in tools:::.is_ASCII().
        if(any(ind))
            Encoding(x[ind]) <- "bytes"
    }

    ## What should we do about NA matches (from matching a non-NA
    ## pattern on an NA string)?  For now, let us always "drop" them so
    ## that extracting direct and inverse matches always gives nothing.

    if(!ili && !invert) {
        so <- m[ind <- (!is.na(m) & (m > -1L))]
        eo <- so + attr(m, "match.length")[ind] - 1L
        return(substring(x[ind], so, eo))
    }

    y <- if(invert) {
        Map(function(u, so, ml) {
            if((n <- length(so)) == 1L) {
                if(is.na(so))
                    return(character())
                else if(so == -1L)
                    return(u)
            }
            beg <- if(n > 1L) {
                ## regexec() could give overlapping matches.
                ## Matches are non-overlapping iff
                ##   eo[i] < so[i + 1], i = 1, ..., n - 1.
                eo <- so + ml - 1L
                if(any(eo[-n] >= so[-1L]))
                    stop(gettextf("need non-overlapping matches for %s",
                                  sQuote("invert = TRUE")),
                         domain = NA)
                c(1L, eo + 1L)
            } else {
                c(1L, so + ml)
            }
            end <- c(so - 1L, nchar(u))
            substring(u, beg, end)
        },
            x, m,
            if(ili)
            lapply(m, attr, "match.length")
            else
            attr(m, "match.length"),
            USE.NAMES = FALSE)
    } else {
        Map(function(u, so, ml) {
            if(length(so) == 1L) {
                if(is.na(so) || (so == -1L))
                    return(character())
            }
            substring(u, so, so + ml - 1L)
        },
            x, m,
            lapply(m, attr, "match.length"),
            USE.NAMES = FALSE)
    }

    names(y) <- names(x)
    y
}

`regmatches<-` <-
function(x, m, invert = FALSE, value)
{
    if(!length(x)) return(x)

    y <- regmatches(x, m, !invert)

    ili <- is.list(m)

    if(!ili && !invert) {
        ## For non-list m and invert = FALSE, we need a character vector
        ## of replacement values with length the number of matched
        ## elements.
        value <- as.character(value)
        if(any(is.na(value)))
            stop("missing replacement values are not allowed")
        ## Entries for matched elements have length 2.
        pos <- which(sapply(y, length) == 2L)
        np <- length(pos)
        nv <- length(value)
        if(np != nv) {
            if(!nv)
                stop("must have replacement values for matches")
            value <- rep_len(value, np)
        }
        y <- y[pos]
        x[pos] <- paste0(sapply(y, `[`, 1L), value, sapply(y, `[`, 2L))
        return(x)
    }

    ## We need a list of character vectors without missings, which has
    ## the same length as x.
    value <- lapply(value, as.character)
    if(any(is.na(unlist(value))))
        stop("missing replacement values are not allowed")
    if(!length(value))
        stop("value does not provide any replacement values")
    value <- rep_len(value, length(x))

    y <- if(invert) {
        ## Replace non-matches.
        ## An element of x with k matches has a corresponding y element
        ## of length k, and needs k + 1 replacement values.
        Map(function(u, v) {
            nu <- length(u)
            nv <- length(v)
            if(nv != (nu + 1L)) {
                if(!nv)
                    stop("must have replacements for non-matches")
                v <- rep_len(v, nu + 1L)
            }
            paste0(v, c(u, ""), collapse = "")
        },
            y, value, USE.NAMES = FALSE)
    } else {
        ## Replace matches.
        ## An element of x with k matches has a corresponding y element
        ## of length k + 1, and needs k replacement values.
        Map(function(u, v) {
            nu <- length(u)
            nv <- length(v)
            if(nv != (nu - 1L)) {
                if(!nv)
                    stop("must have replacements for matches")
                v <- rep_len(v, nu - 1L)
            }
            paste0(u, c(v, ""), collapse = "")
        },
            y, value, USE.NAMES = FALSE)
    }

    y <- unlist(y)
    names(y) <- names(x)

    y
}
#  File src/library/base/R/identical.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

identical <- function(x, y, num.eq = TRUE, single.NA = TRUE,
                      attrib.as.set = TRUE, ignore.bytecode = TRUE,
                      ignore.environment = FALSE)
    .Internal(identical(x,y, num.eq, single.NA, attrib.as.set,
                        ignore.bytecode, ignore.environment))

isTRUE <- function(x) identical(TRUE, x)
#  File src/library/base/R/ifelse.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

ifelse <- function (test, yes, no)
{
    if(is.atomic(test)) # do not lose attributes
	storage.mode(test) <- "logical"
    else ## typically a "class"; storage.mode<-() typically fails
	test <- if(isS4(test)) as(test, "logical") else as.logical(test)
    ans <- test
    ok <- !(nas <- is.na(test))
    if (any(test[ok]))
	ans[test & ok] <- rep(yes, length.out = length(ans))[test & ok]
    if (any(!test[ok]))
	ans[!test & ok] <- rep(no, length.out = length(ans))[!test & ok]
    ans[nas] <- NA
    ans
}
#  File src/library/base/R/interaction.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

### This is almost like the Primitive ":" for factors
### but with drop=TRUE, used in reshape
interaction <- function(..., drop = FALSE, sep = ".", lex.order = FALSE)
{
    args <- list(...)
    narg <- length(args)
    if (narg == 1L && is.list(args[[1L]])) {
	args <- args[[1L]]
	narg <- length(args)
    }
    for(i in narg:1L) {
        f <- as.factor(args[[i]])[, drop = drop]
        l <- levels(f)
        if1 <- as.integer(f) - 1L
        if(i == narg) {
            ans <- if1
            lvs <- l
        } else {
            if(lex.order) {
                ll <- length(lvs)
                ans <- ans + ll * if1
                lvs <- paste(rep(l, each = ll), rep(lvs, length(l)), sep=sep)
            } else {
                ans <- ans * length(l) + if1
                lvs <- paste(rep(l, length(lvs)),
                             rep(lvs, each = length(l)), sep=sep)
            }
            if(anyDuplicated(lvs)) { ## fix them up
                ulvs <- unique(lvs)
                while((i <- anyDuplicated(flv <- match(lvs, ulvs)))) {
                    lvs <- lvs[-i]
                    ans[ans+1L == i] <- match(flv[i], flv[1:(i-1)]) - 1L
                    ans[ans+1L > i] <- ans[ans+1L > i] - 1L
                }
                lvs <- ulvs
            }
            if(drop) {
                olvs <- lvs
                lvs <- lvs[sort(unique(ans+1L))]
                ans <- match(olvs[ans+1L], lvs) - 1L
            }
        }
    }
    structure(as.integer(ans+1L), levels=lvs, class = "factor")
}
#  File src/library/base/R/is.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

is.vector <- function(x, mode="any") .Internal(is.vector(x,mode))

`is.na<-` <- function(x, value) UseMethod("is.na<-")

`is.na<-.default` <- function(x, value)
{
    x[value] <- NA
    x
}

is.primitive <- function(x)
    switch(typeof(x), "special" = , "builtin" = TRUE, FALSE)
#  File src/library/base/R/jitter.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

### Unimplemented Idea {for amount = NULL ?}
### Really "optimal" (e.g. for rug()), use a non-constant amount,
### e.g. use "d" = diff(xx)  BEFORE  taking min()...

jitter <- function(x, factor = 1, amount=NULL)
{
    if(length(x) == 0L)
	return(x)
    if(!is.numeric(x))
        stop("'x' must be numeric")
    z <- diff(r <- range(x[is.finite(x)]))
    if(z == 0) z <- abs(r[1L])
    if(z == 0) z <- 1

    if(is.null(amount)) {		# default: Find 'necessary' amount
	d <- diff(xx <- unique(sort.int(round(x, 3 - floor(log10(z))))))
	d <- if(length(d)) min(d) else if(xx != 0) xx/10 else z/10
	amount <- factor/5 * abs(d)
    } else if(amount == 0)		# only then: S compatibility
	amount <- factor * (z/50)

    x + stats::runif(length(x),  - amount, amount)
}
#  File src/library/base/R/kappa.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1998 B. D. Ripley
#  Copyright (C) 1998-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

norm <- function(x, type = c("O", "I", "F", "M", "2")) {
    if(identical("2", type)) {
	svd(x, nu = 0L, nv = 0L)$d[1L]
	## *faster* at least on some platforms {but possibly less accurate}:
	##sqrt(eigen(crossprod(x), symmetric=TRUE, only.values=TRUE)$values[1L])
    } else
	.Internal(La_dlange(x, type))
} ## and define it as implicitGeneric, so S4 methods are consistent

kappa <- function(z, ...) UseMethod("kappa")

## Note that  all 4 Lapack version now work in the following
rcond <- function(x, norm = c("O","I","1"), triangular = FALSE, ...) {
    norm <- match.arg(norm)
    stopifnot(is.matrix(x))
    if({d <- dim(x); d[1L] != d[2L]})## non-square matrix -- use QR
        return(rcond(qr.R(qr(if(d[1L] < d[2L]) t(x) else x)), norm=norm, ...))

    ## x = square matrix :
    if(is.complex(x)) {
        if(triangular) .Internal(La_ztrcon(x, norm))
        else .Internal(La_zgecon(x, norm))
    } else {
        if(triangular) .Internal(La_dtrcon(x, norm))
        else .Internal(La_dgecon(x, norm))
    }
}

kappa.default <- function(z, exact = FALSE,
                          norm = NULL, method = c("qr", "direct"), ...)
{
    method <- match.arg(method)
    z <- as.matrix(z)
    norm <- if(!is.null(norm)) match.arg(norm, c("2", "1","O", "I")) else "2"
    if(exact && norm == "2") {
        s <- svd(z, nu = 0, nv = 0)$d
        max(s)/min(s[s > 0])
    }
    else { ## exact = FALSE or norm in "1", "O", "I"
	if(exact)
	    warning(gettextf("norm '%s' currently always uses exact = FALSE",
			     norm))
        d <- dim(z)
        if(method == "qr" || d[1L] != d[2L])
	    kappa.qr(qr(if(d[1L] < d[2L]) t(z) else z),
		     exact = FALSE, norm = norm, ...)
        else .kappa_tri(z, exact = FALSE, norm = norm, ...)
    }
}

kappa.lm <- function(z, ...) kappa.qr(z$qr, ...)

kappa.qr <- function(z, ...)
{
    qr <- z$qr
    R <- qr[1L:min(dim(qr)), , drop = FALSE]
    R[lower.tri(R)] <- 0
    .kappa_tri(R, ...)
}

.kappa_tri <- function(z, exact = FALSE, LINPACK = TRUE, norm=NULL, ...)
{
    if(exact) {
        stopifnot(is.null(norm) || identical("2", norm))
        kappa.default(z, exact = TRUE) ## using "2 - norm" !
    }
    else { ## norm is "1" ("O") or "I(nf)" :
        p <- as.integer(nrow(z))
        if(is.na(p)) stop("invalid nrow(x)")
	if(p != ncol(z)) stop("triangular matrix should be square")
	if(is.null(norm)) norm <- "1"
	if(is.complex(z)) 1/.Internal(La_ztrcon(z, norm))
	else if(LINPACK) {
	    if(norm == "I") # instead of "1" / "O"
		z <- t(z)
	    ##	dtrco  *differs* from Lapack's dtrcon() quite a bit
	    ## even though dtrco's doc also say to compute the
	    ## 1-norm reciprocal condition
            storage.mode(z) <- "double"
	    1 / .Fortran(.F_dtrco, z, p, p, k = double(1), double(p), 1L)$k
	}
	else 1/.Internal(La_dtrcon(z, norm))
    }
}
#  File src/library/base/R/kronecker.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

kronecker <- function (X, Y, FUN = "*", make.dimnames = FALSE, ...)
{
    ## This is principally to allow Matrix/SparseM to set S4 methods
    ## on %x%, which calls base::kronecker.
    if (.isMethodsDispatchOn() && (isS4(X) || isS4(Y))) {
        return(methods::kronecker(X, Y, FUN = FUN,
                                  make.dimnames = make.dimnames, ...))
    }
    .kronecker(X, Y, FUN = FUN, make.dimnames = make.dimnames, ...)
}

.kronecker <- function (X, Y, FUN = "*", make.dimnames = FALSE, ...)
{
    X <- as.array(X)
    Y <- as.array(Y)
    if (make.dimnames) {
	dnx <- dimnames(X)
	dny <- dimnames(Y)
    }
    dX <- dim(X)
    dY <- dim(Y)
    ld <- length(dX) - length(dY)
    if (ld < 0L)
	dX <- dim(X) <- c(dX, rep.int(1, -ld))
    else if (ld > 0L)
	dY <- dim(Y) <- c(dY, rep.int(1, ld))
    opobj <- outer(X, Y, FUN, ...)
    dp <- as.vector(t(matrix(1L:(2*length(dX)), ncol = 2)[, 2:1]))
    opobj <- aperm(opobj, dp)
    dim(opobj) <- dX * dY

    if (make.dimnames && !(is.null(dnx) && is.null(dny))) {
	if (is.null(dnx))
	    dnx <- vector("list", length(dX))
	else if (ld < 0L)
	    dnx <- c(dnx, vector("list", -ld))
	tmp <- which(sapply(dnx, is.null))
	dnx[tmp] <- lapply(tmp, function(i) rep.int("", dX[i]))

	if (is.null(dny))
	    dny <- vector("list", length(dY))
	else if (ld > 0)
	    dny <- c(dny, vector("list", ld))
	tmp <- which(sapply(dny, is.null))
	dny[tmp] <- lapply(tmp, function(i) rep.int("", dY[i]))

	k <- length(dim(opobj))
	dno <- vector("list", k)
	for (i in 1L:k) {
	    tmp <- outer(dnx[[i]], dny[[i]], FUN="paste", sep=":")
	    dno[[i]] <- as.vector(t(tmp))
	}
	dimnames(opobj) <- dno
    }
    opobj
}

## Binary operator, hence don't simply do "%x%" <- kronecker.
`%x%` <- function(X, Y) kronecker(X, Y)
#  File src/library/base/R/labels.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1998 B. D. Ripley
#  Copyright (C) 1998-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

labels <- function(object, ...) UseMethod("labels")

labels.default <- function(object, ...)
{
    if(length(d <- dim(object))) {	# array or data frame
	nt <- dimnames(object)
	if(is.null(nt)) nt <- vector("list", length(d))
	for(i in seq_along(d))
	    if(!length(nt[[i]])) nt[[i]] <- as.character(seq_len(d[i]))
    } else {
	nt <- names(object)
	if(!length(nt)) nt <- as.character(seq_along(object))
    }
    nt
}
#  File src/library/base/R/lapply.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

lapply <- function (X, FUN, ...)
{
    FUN <- match.fun(FUN)
    ## internal code handles all vector types, including expressions
    ## However, it would be OK to have attributes which is.vector
    ## disallows.
    if(!is.vector(X) || is.object(X)) X <- as.list(X)
    ## Note ... is not passed down.  Rather the internal code
    ## evaluates FUN(X[i], ...) in the frame of this function
    .Internal(lapply(X, FUN))
}

rapply <-
    function(object, f, classes = "ANY", deflt = NULL,
             how = c("unlist", "replace", "list"), ...)
{
    if(typeof(object) != "list")
        stop("'object' must be a list")
    how <- match.arg(how)
    res <- .Internal(rapply(object, f, classes, deflt, how))
    if(how == "unlist") unlist(res, recursive = TRUE) else res
}
#  File src/library/base/R/lazyload.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

## This code should be kept in step with code in ../baseloader.R
##
## This code has been factored in a somewhat peculiar way to allow the
## lazy load data base mechanism to be used for storing processed .Rd
## files. This isn't quite right as the .Rd use only uses the data
## base, not the lazy load part, but for now it will do. LT

lazyLoadDBexec <- function(filebase, fun, filter)
{
    ##
    ## bootstrapping definitions so we can load base
    ## - not that this version is actually used to load base
    ##
    glue <- function (..., sep = " ", collapse = NULL)
        .Internal(paste(list(...), sep, collapse))
    readRDS <- function (file) {
        halt <- function (message) .Internal(stop(TRUE, message))
        gzfile <- function (description, open)
            .Internal(gzfile(description, open, "", 6))
        close <- function (con) .Internal(close(con, "rw"))
        if (! is.character(file)) halt("bad file name")
        con <- gzfile(file, "rb")
        on.exit(close(con))
        .Internal(unserializeFromConn(con, baseenv()))
    }
    `parent.env<-` <-
        function (env, value) .Internal(`parent.env<-`(env, value))
    existsInFrame <- function (x, env) .Internal(exists(x, env, "any", FALSE))
    getFromFrame <- function (x, env) .Internal(get(x, env, "any", FALSE))
    set <- function (x, value, env) .Internal(assign(x, value, env, FALSE))
    environment <- function () .Internal(environment(NULL))
    mkenv <- function() .Internal(new.env(TRUE, baseenv(), 29L))

    ##
    ## main body
    ##
    mapfile <- glue(filebase, "rdx", sep = ".")
    datafile <- glue(filebase, "rdb", sep = ".")
    env <- mkenv()
    map <- readRDS(mapfile)
    vars <- names(map$variables)
    rvars <- names(map$references)
    compressed <- map$compressed
    for (i in seq_along(rvars))
        set(rvars[i], map$references[[i]], env)
    envenv <- mkenv()
    envhook <- function(n) {
        if (existsInFrame(n, envenv))
            getFromFrame(n, envenv)
        else {
            e <- mkenv()
            set(n, e, envenv)           # MUST do this immediately
            key <- getFromFrame(n, env)
            data <- lazyLoadDBfetch(key, datafile, compressed, envhook)
            ## comment from r41494
            ## modified the loading of old environments, so that those
            ## serialized with parent.env NULL are loaded with the
            ## parent.env=emptyenv(); and yes an alternative would have been
            ## baseenv(), but that was seldom the intention of folks that
            ## set the environment to NULL.
            if (is.null(data$enclos))
                parent.env(e) <- emptyenv()
            else
                parent.env(e) <- data$enclos
            vars <- names(data$bindings)
            for (i in seq_along(vars))
                set(vars[i], data$bindings[[i]], e)
            if (! is.null(data$attributes))
                attributes(e) <- data$attributes
            if (! is.null(data$isS4) && data$isS4)
                .Internal(setS4Object(e, TRUE, TRUE))
            if (! is.null(data$locked) && data$locked)
                .Internal(lockEnvironment(e, FALSE))
            e
        }
    }
    if (!missing(filter)) {
        use <- filter(vars)
        vars <- vars[use]
        vals <- map$variables[use]
        use <- NULL
    } else
        vals <-  map$variables

    res <- fun(environment())

    ## reduce memory use
    map <- NULL
    vars <- NULL
    vals <- NULL
    rvars <- NULL
    mapfile <- NULL
    readRDS <- NULL

    res
}

lazyLoad <- function(filebase, envir = parent.frame(), filter)
{
    fun <- function(db) {
        vals <- db$vals
        vars <- db$vars
        expr <- quote(lazyLoadDBfetch(key, datafile, compressed, envhook))
        .Internal(makeLazy(vars, vals, expr, db, envir))
    }
    lazyLoadDBexec(filebase, fun, filter)
}
#  File src/library/base/R/library.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

testPlatformEquivalence <-
function(built, run)
{
    ## args are "cpu-vendor-os", but os might be 'linux-gnu'!
    ## remove vendor field
    built <- gsub("([^-]*)-([^-]*)-(.*)", "\\1-\\3", built)
    run <- gsub("([^-]*)-([^-]*)-(.*)", "\\1-\\3", run)
    ## Mac OS X supports multiple CPUs by using 'universal' binaries
    if (length(grep("^universal-darwin", built)) && nzchar(.Platform$r_arch))
        built <- sub("^universal", R.version$arch, built)
    ## allow for small mismatches, e.g. OS version number and i686 vs i586.
    length(agrep(built, run)) > 0
}

library <-
function(package, help, pos = 2, lib.loc = NULL, character.only = FALSE,
         logical.return = FALSE, warn.conflicts = TRUE,
	 quietly = FALSE, verbose = getOption("verbose"))
{
    testRversion <- function(pkgInfo, pkgname, pkgpath)
    {
        if(is.null(built <- pkgInfo$Built))
            stop(gettextf("package %s has not been installed properly\n",
                          sQuote(pkgname)),
                 call. = FALSE, domain = NA)

        ## which version was this package built under?
        R_version_built_under <- as.numeric_version(built$R)
        if(R_version_built_under < "3.0.0")
            stop(gettextf("package %s was built before R 3.0.0: please re-install it",
                          sQuote(pkgname)), call. = FALSE, domain = NA)

        current <- getRversion()
        ## depends on R version?
        ## as it was installed >= 2.7.0 it will have Rdepends2
        if(length(Rdeps <- pkgInfo$Rdepends2)) {
            for(dep in Rdeps)
                if(length(dep) > 1L) {
                    target <- dep$version
                    res <- if(is.character(target)) {
                        do.call(dep$op, # these are both strings
                                list(as.numeric(R.version[["svn rev"]]),
                                     as.numeric(sub("^r", "", dep$version))))
                    } else {
                        do.call(dep$op,
                                list(current, as.numeric_version(target)))
##                        target <- as.numeric_version(dep$version)
##                        eval(parse(text=paste("current", dep$op, "target")))
                    }
                    if(!res)
                        stop(gettextf("This is R %s, package %s needs %s %s",
                                      current, sQuote(pkgname), dep$op, target),
                             call. = FALSE, domain = NA)
                }
        }
        ## warn if installed under a later version of R
        if(R_version_built_under > current)
            warning(gettextf("package %s was built under R version %s",
                             sQuote(pkgname), as.character(built$R)),
                    call. = FALSE, domain = NA)
        platform <- built$Platform
        r_arch <- .Platform$r_arch
        if(.Platform$OS.type == "unix") {
            ## allow mismatches if r_arch is in use, e.g.
            ## i386-gnu-linux vs x86-gnu-linux depending on
            ## build system.
            if(!nzchar(r_arch) && length(grep("\\w", platform)) &&
               !testPlatformEquivalence(platform, R.version$platform))
                stop(gettextf("package %s was built for %s",
                              sQuote(pkgname), platform),
                     call. = FALSE, domain = NA)
        } else {  # Windows
            ## a check for 'mingw' suffices, since i386 and x86_64
            ## have DLLs in different places.  This allows binary packages
            ## to be merged.
            if(nzchar(platform) && !grepl("mingw", platform))
                stop(gettextf("package %s was built for %s",
                              sQuote(pkgname), platform),
                     call. = FALSE, domain = NA)
        }
        ## if using r_arch subdirs, check for presence
        if(nzchar(r_arch)
           && file.exists(file.path(pkgpath, "libs"))
           && !file.exists(file.path(pkgpath, "libs", r_arch)))
            stop(gettextf("package %s is not installed for 'arch = %s'",
                          sQuote(pkgname), r_arch),
                 call. = FALSE, domain = NA)
    }

    checkLicense <- function(pkg, pkgInfo, pkgPath)
    {
        L <- tools:::analyze_license(pkgInfo$DESCRIPTION["License"])
        if(!L$is_empty && !L$is_verified) {
            site_file <- path.expand(file.path(R.home("etc"), "licensed.site"))
            if(file.exists(site_file) &&
               pkg %in% readLines(site_file)) return()
            personal_file <- path.expand("~/.R/licensed")
            if(file.exists(personal_file)) {
                agreed <- readLines(personal_file)
                if(pkg %in% agreed) return()
            } else agreed <- character()
            if(!interactive())
                stop(gettextf("package %s has a license that you need to accept in an interactive session", sQuote(pkg)), domain = NA)
            lfiles <- file.path(pkgpath, c("LICENSE", "LICENCE"))
            lfiles <- lfiles[file.exists(lfiles)]
            if(length(lfiles)) {
                message(gettextf("package %s has a license that you need to accept after viewing", sQuote(pkg)), domain = NA)
                readline("press RETURN to view license")
                encoding <- pkgInfo$DESCRIPTION["Encoding"]
                if(is.na(encoding)) encoding <- ""
                ## difR and EVER have a Windows' 'smart quote' LICEN[CS]E file
                if(encoding == "latin1") encoding <- "cp1252"
                file.show(lfiles[1L], encoding = encoding)
            } else {
                message(gettextf("package %s has a license that you need to accept:\naccording to the DESCRIPTION file it is", sQuote(pkg)), domain = NA)
                message(pkgInfo$DESCRIPTION["License"], domain = NA)
            }
            choice <- menu(c("accept", "decline"),
                           title = paste("License for", sQuote(pkg)))
            if(choice != 1)
                stop(gettextf("license for package %s not accepted",
                              sQuote(package)), domain = NA, call. = FALSE)
            dir.create(dirname(personal_file), showWarnings=FALSE)
            writeLines(c(agreed, pkg), personal_file)
        }
    }

    checkNoGenerics <- function(env, pkg)
    {
        nenv <- env
        ns <- .getNamespace(as.name(pkg))
        if(!is.null(ns)) nenv <- asNamespace(ns)
        if (exists(".noGenerics", envir = nenv, inherits = FALSE))
            TRUE
        else {
            ## A package will have created a generic
            ## only if it has created a formal method.
            length(objects(env, pattern="^\\.__[MT]", all.names=TRUE)) == 0L
        }
    }

    checkConflicts <- function(package, pkgname, pkgpath, nogenerics, env)
    {
        dont.mind <- c("last.dump", "last.warning", ".Last.value",
                       ".Random.seed", ".Last.lib", ".onDetach",
                       ".packageName", ".noGenerics", ".required",
                       ".no_S3_generics", ".Depends", ".requireCachedGenerics")
        sp <- search()
        lib.pos <- match(pkgname, sp)
        ## ignore generics not defined for the package
        ob <- objects(lib.pos, all.names = TRUE)
        if(!nogenerics) {
            ##  Exclude generics that are consistent with implicit generic
            ## from another package.  A better test would be to move this
            ## down into the loop and test against specific other package name
            ## but subtle conflicts like that are likely to be found elsewhere
	    these <- ob[substr(ob, 1L, 6L) == ".__T__"]
            gen <- gsub(".__T__(.*):([^:]+)", "\\1", these)
            from <- gsub(".__T__(.*):([^:]+)", "\\2", these)
            gen <- gen[from != package]
            ob <- ob[!(ob %in% gen)]
        }
        fst <- TRUE
	ipos <- seq_along(sp)[-c(lib.pos,
				 match(c("Autoloads", "CheckExEnv"), sp, 0L))]
        for (i in ipos) {
            obj.same <- match(objects(i, all.names = TRUE), ob, nomatch = 0L)
            if (any(obj.same > 0)) {
                same <- ob[obj.same]
                same <- same[!(same %in% dont.mind)]
                Classobjs <- grep("^\\.__", same)
                if(length(Classobjs)) same <- same[-Classobjs]
                ## report only objects which are both functions or
                ## both non-functions.
		same.isFn <- function(where)
		    vapply(same, exists, NA,
                           where = where, mode = "function", inherits = FALSE)
		same <- same[same.isFn(i) == same.isFn(lib.pos)]
		## if a package imports and re-exports, there's no problem
		not.Ident <- function(ch, TRAFO=identity, ...)
		    vapply(ch, function(.)
                           !identical(TRAFO(get(., i)),
                                      TRAFO(get(., lib.pos)), ...),
                           NA)
		if(length(same)) same <- same[not.Ident(same)]
		## if the package is 'base' it cannot be imported and re-exported,
		## allow a "copy":
		if(length(same) && identical(sp[i], "package:base"))
		    same <- same[not.Ident(same, ignore.environment = TRUE)]
                if(length(same)) {
                    if (fst) {
                        fst <- FALSE
                        packageStartupMessage(gettextf("\nAttaching package: %s\n",
                                                       sQuote(package)),
                                              domain = NA)
                    }

                    objs <- strwrap(paste(same, collapse=", "),
                                    indent = 4L, exdent = 4L)
                    msg <- sprintf(ngettext(length(same),
                                            "The following object is masked %s %s:\n\n%s\n",
                                            "The following objects are masked %s %s:\n\n%s\n"),
                                   if (i < lib.pos) "_by_" else "from",
                                   sQuote(sp[i]), paste(objs, collapse="\n"))
		    packageStartupMessage(msg)
                }
            }
        }
    }

    if(verbose && quietly)
	message("'verbose' and 'quietly' are both true; being verbose then ..")
    if(!missing(package)) {
        if (is.null(lib.loc)) lib.loc <- .libPaths()
        ## remove any non-existent directories
        lib.loc <- lib.loc[file.info(lib.loc)$isdir %in% TRUE]

	if(!character.only)
	    package <- as.character(substitute(package))
        if(length(package) != 1L)
            stop("'package' must be of length 1")
        if(is.na(package) || (package == ""))
            stop("invalid package name")

	pkgname <- paste("package", package, sep = ":")
	newpackage <- is.na(match(pkgname, search()))
	if(newpackage) {
            ## Check for the methods package before attaching this
            ## package.
            ## Only if it is _already_ here do we do cacheMetaData.
            ## The methods package caches all other pkgs when it is
            ## attached.

            pkgpath <- find.package(package, lib.loc, quiet = TRUE,
                                    verbose = verbose)
            if(length(pkgpath) == 0L) {
                txt <- if(length(lib.loc))
                    gettextf("there is no package called %s", sQuote(package))
                else
                    gettext("no library trees found in 'lib.loc'")
                if(logical.return) {
                    warning(txt, domain = NA)
		    return(FALSE)
		} else stop(txt, domain = NA)
            }
            which.lib.loc <- normalizePath(dirname(pkgpath), "/", TRUE)
            pfile <- system.file("Meta", "package.rds", package = package,
                                 lib.loc = which.lib.loc)
            if(!nzchar(pfile))
            	stop(gettextf("%s is not a valid installed package",
                              sQuote(package)), domain = NA)
            pkgInfo <- readRDS(pfile)
            testRversion(pkgInfo, package, pkgpath)
            ## avoid any bootstrapping issues by these exemptions
            if(!package %in% c("datasets", "grDevices", "graphics", "methods",
                               "splines", "stats", "stats4", "tcltk", "tools",
                               "utils") &&
               isTRUE(getOption("checkPackageLicense", FALSE)))
                checkLicense(package, pkgInfo, pkgpath)

            ## The check for inconsistent naming is now in find.package

            if(is.character(pos)) {
                npos <- match(pos, search())
                if(is.na(npos)) {
                    warning(gettextf("%s not found on search path, using pos = 2", sQuote(pos)), domain = NA)
                    pos <- 2
                } else pos <- npos
            }
            .getRequiredPackages2(pkgInfo, quietly = quietly)
            deps <- unique(names(pkgInfo$Depends))

            ## If the namespace mechanism is available and the package
            ## has a namespace, then the namespace loading mechanism
            ## takes over.
            if (packageHasNamespace(package, which.lib.loc)) {
                tt <- try({
                    ns <- loadNamespace(package, c(which.lib.loc, lib.loc))
                    env <- attachNamespace(ns, pos = pos, deps)
                })
                if (inherits(tt, "try-error"))
                    if (logical.return)
                        return(FALSE)
                    else stop(gettextf("package or namespace load failed for %s",
                                       sQuote(package)),
                              call. = FALSE, domain = NA)
                else {
                    on.exit(detach(pos = pos))
                    ## If there are S4 generics then the package should
                    ## depend on methods
                    nogenerics <-
                        !.isMethodsDispatchOn() || checkNoGenerics(env, package)
                    if(warn.conflicts && # never will with a namespace
                       !exists(".conflicts.OK", envir = env, inherits = FALSE))
                        checkConflicts(package, pkgname, pkgpath,
                                       nogenerics, ns)
                    on.exit()
                    if (logical.return)
                        return(TRUE)
                    else
                        return(invisible(.packages()))
                }
            } else
            stop(gettextf("package %s does not have a namespace and should be re-installed",
                          sQuote(package)), domain = NA)
	}
	if (verbose && !newpackage)
            warning(gettextf("package %s already present in search()",
                             sQuote(package)), domain = NA)

    }
    else if(!missing(help)) {
	if(!character.only)
	    help <- as.character(substitute(help))
        pkgName <- help[1L]            # only give help on one package
        pkgPath <- find.package(pkgName, lib.loc, verbose = verbose)
        docFiles <- c(file.path(pkgPath, "Meta", "package.rds"),
                      file.path(pkgPath, "INDEX"))
        if(file.exists(vignetteIndexRDS <-
                       file.path(pkgPath, "Meta", "vignette.rds")))
            docFiles <- c(docFiles, vignetteIndexRDS)
        pkgInfo <- vector("list", 3L)
        readDocFile <- function(f) {
            if(basename(f) %in% "package.rds") {
                txt <- readRDS(f)$DESCRIPTION
                if("Encoding" %in% names(txt)) {
                    to <- if(Sys.getlocale("LC_CTYPE") == "C") "ASCII//TRANSLIT"else ""
                    tmp <- try(iconv(txt, from=txt["Encoding"], to=to))
                    if(!inherits(tmp, "try-error"))
                        txt <- tmp
                    else
                        warning("'DESCRIPTION' has an 'Encoding' field and re-encoding is not possible", call.=FALSE)
                }
                nm <- paste0(names(txt), ":")
                formatDL(nm, txt, indent = max(nchar(nm, "w")) + 3)
            } else if(basename(f) %in% "vignette.rds") {
                txt <- readRDS(f)
                ## New-style vignette indices are data frames with more
                ## info than just the base name of the PDF file and the
                ## title.  For such an index, we give the names of the
                ## vignettes, their titles, and indicate whether PDFs
                ## are available.
                ## The index might have zero rows.
                if(is.data.frame(txt) && nrow(txt))
                    cbind(basename(gsub("\\.[[:alpha:]]+$", "",
                                        txt$File)),
                          paste(txt$Title,
                                paste0(rep.int("(source", NROW(txt)),
                                       ifelse(txt$PDF != "",
                                              ", pdf",
                                              ""),
                                       ")")))
                else NULL
            } else
            readLines(f)
        }
        for(i in which(file.exists(docFiles)))
            pkgInfo[[i]] <- readDocFile(docFiles[i])
        y <- list(name = pkgName, path = pkgPath, info = pkgInfo)
        class(y) <- "packageInfo"
        return(y)
    }
    else {
	## library():
        if(is.null(lib.loc))
            lib.loc <- .libPaths()
        db <- matrix(character(), nrow = 0L, ncol = 3L)
        nopkgs <- character()

        for(lib in lib.loc) {
            a <- .packages(all.available = TRUE, lib.loc = lib)
            for(i in sort(a)) {
                ## All packages installed under 2.0.0 should have
                ## 'package.rds' but we have not checked.
                file <- system.file("Meta", "package.rds", package = i,
                                    lib.loc = lib)
                title <- if(file != "") {
                    txt <- readRDS(file)
                    if(is.list(txt)) txt <- txt$DESCRIPTION
                    ## we may need to re-encode here.
                    if("Encoding" %in% names(txt)) {
                        to <- if(Sys.getlocale("LC_CTYPE") == "C") "ASCII//TRANSLIT" else ""
                        tmp <- try(iconv(txt, txt["Encoding"], to, "?"))
                        if(!inherits(tmp, "try-error"))
                            txt <- tmp
                        else
                            warning("'DESCRIPTION' has an 'Encoding' field and re-encoding is not possible", call.=FALSE)
                    }
                    txt["Title"]
                } else NA
                if(is.na(title))
                    title <- " ** No title available ** "
                db <- rbind(db, cbind(i, lib, title))
            }
            if(length(a) == 0L)
                nopkgs <- c(nopkgs, lib)
        }
        dimnames(db) <- list(NULL, c("Package", "LibPath", "Title"))
        if(length(nopkgs) && !missing(lib.loc)) {
            pkglist <- paste(sQuote(nopkgs), collapse = ", ")
            msg <- sprintf(ngettext(length(nopkgs),
                                    "library %s contains no packages",
                                    "libraries %s contain no packages"),
                           pkglist)
            warning(msg, domain=NA)
        }

        y <- list(header = NULL, results = db, footer = NULL)
        class(y) <- "libraryIQR"
        return(y)
    }

    if (logical.return)
	TRUE
    else invisible(.packages())
}

format.libraryIQR <-
function(x, ...)
{
    db <- x$results
    if(!nrow(db)) return(character())
    ## Split according to LibPath, preserving order of libraries.
    libs <- db[, "LibPath"]
    libs <- factor(libs, levels = unique(libs))
    out <- lapply(split(1 : nrow(db), libs),
                  function(ind) db[ind, c("Package", "Title"),
                                   drop = FALSE])
    c(unlist(Map(function(lib, sep) {
        c(gettextf("%sPackages in library %s:\n", sep, sQuote(lib)),
          formatDL(out[[lib]][, "Package"],
                   out[[lib]][, "Title"]))
    },
                 names(out),
                 c("", rep.int("\n", length(out) - 1L)))),
      x$footer)
}

print.libraryIQR <-
function(x, ...)
{
    s <- format(x)
    if(!length(s)) {
        message("no packages found")
    } else {
        outFile <- tempfile("RlibraryIQR")
        writeLines(s, outFile)
        file.show(outFile, delete.file = TRUE,
                  title = gettext("R packages available"))
    }
    invisible(x)
}

library.dynam <-
function(chname, package, lib.loc, verbose = getOption("verbose"),
         file.ext = .Platform$dynlib.ext, ...)
{
    dll_list <- .dynLibs()

    if(missing(chname) || !nzchar(chname)) return(dll_list)

    ## For better error messages, force these to be evaluated.
    package
    lib.loc

    r_arch <- .Platform$r_arch
    chname1 <- paste0(chname, file.ext)
    ## it is not clear we should allow this, rather require a single
    ## package and library.
    for(pkg in find.package(package, lib.loc, verbose = verbose)) {
        DLLpath <- if(nzchar(r_arch)) file.path(pkg, "libs", r_arch)
	else    file.path(pkg, "libs")
        file <- file.path(DLLpath, chname1)
        if(file.exists(file)) break else file <- ""
    }
    if(file == "")
        if(.Platform$OS.type == "windows")
            stop(gettextf("DLL %s not found: maybe not installed for this architecture?", sQuote(chname)), domain = NA)
        else
            stop(gettextf("shared object %s not found", sQuote(chname1)),
                 domain = NA)
    ## for consistency with library.dyn.unload:
    file <- file.path(normalizePath(DLLpath, "/", TRUE), chname1)
    ind <- vapply(dll_list, function(x) x[["path"]] == file, NA)
    if(length(ind) && any(ind)) {
        if(verbose)
            if(.Platform$OS.type == "windows")
                message(gettextf("DLL %s already loaded", sQuote(chname1)),
                        domain = NA)
            else
                message(gettextf("shared object '%s' already loaded",
                                 sQuote(chname1)), domain = NA)
        return(invisible(dll_list[[ seq_along(dll_list)[ind] ]]))
    }
    if(.Platform$OS.type == "windows") {
        ## Make it possible to find other DLLs in the same place as
        ## @code{file}, so that e.g. binary packages can conveniently
        ## provide possibly missing DLL dependencies in this place
        ## (without having to bypass the default package dynload
        ## mechanism).  Note that this only works under Windows, and a
        ## more general solution will have to be found eventually.
        ##
        ## 2.7.0: there's a more general mechanism in DLLpath=,
        ## so not clear if this is still needed.
        PATH <- Sys.getenv("PATH")
        Sys.setenv(PATH = paste(gsub("/", "\\\\", DLLpath), PATH, sep=";"))
        on.exit(Sys.setenv(PATH = PATH))
    }
    if(verbose)
        message(gettextf("now dyn.load(\"%s\") ...", file), domain = NA)
    dll <- if("DLLpath" %in% names(list(...))) dyn.load(file, ...)
    else dyn.load(file, DLLpath = DLLpath, ...)
    .dynLibs(c(dll_list, list(dll)))
    invisible(dll)
}

library.dynam.unload <-
function(chname, libpath, verbose = getOption("verbose"),
         file.ext = .Platform$dynlib.ext)
{
    dll_list <- .dynLibs()

    if(missing(chname) || (nc_chname <- nchar(chname, "c")) == 0L)
        if(.Platform$OS.type == "windows")
            stop("no DLL was specified")
        else
            stop("no shared object was specified")

    ## We need an absolute path here, and separators consistent with
    ## library.dynam
    libpath <- normalizePath(libpath, "/", TRUE)
    chname1 <- paste0(chname, file.ext)
    file <- if(nzchar(.Platform$r_arch))
             file.path(libpath, "libs", .Platform$r_arch, chname1)
     else    file.path(libpath, "libs", chname1)

    pos <- which(vapply(dll_list, function(x) x[["path"]] == file, NA))
    if(!length(pos))
        if(.Platform$OS.type == "windows")
            stop(gettextf("DLL %s was not loaded", sQuote(chname1)),
                 domain = NA)
        else
            stop(gettextf("shared object %s was not loaded", sQuote(chname1)),
                 domain = NA)

    if(!file.exists(file))
        if(.Platform$OS.type == "windows")
            stop(gettextf("DLL %s not found", sQuote(chname1)), domain = NA)
        else
            stop(gettextf("shared object '%s' not found", sQuote(chname1)),
                 domain = NA)
    if(verbose)
        message(gettextf("now dyn.unload(\"%s\") ...", file), domain = NA)
    dyn.unload(file)
    .dynLibs(dll_list[-pos])
    invisible(dll_list[[pos]])
}

require <-
function(package, lib.loc = NULL, quietly = FALSE, warn.conflicts = TRUE,
         character.only = FALSE)
{
    if(!character.only)
        package <- as.character(substitute(package)) # allowing "require(eda)"
    loaded <- paste("package", package, sep = ":") %in% search()

    if (!loaded) {
	if (!quietly)
            packageStartupMessage(gettextf("Loading required package: %s",
                                           package), domain = NA)
	value <- tryCatch(library(package, lib.loc = lib.loc,
                                  character.only = TRUE,
                                  logical.return = TRUE,
                                  warn.conflicts = warn.conflicts,
				  quietly = quietly),
                          error = function(e) e)
        if (inherits(value, "error")) {
            if (!quietly) {
                msg <- conditionMessage(value)
                cat("Failed with error:  ",
                    sQuote(msg), "\n", file = stderr(), sep = "")
                .Internal(printDeferredWarnings())
            }
            return(invisible(FALSE))
        }
        if (!value) return(invisible(FALSE))
    } else value <- TRUE
    invisible(value)
}

.packages <-
function(all.available = FALSE, lib.loc = NULL)
{
    if(is.null(lib.loc))
        lib.loc <- .libPaths()
    if(all.available) {
	ans <- character()
        for(lib in lib.loc[file.exists(lib.loc)]) {
            a <- list.files(lib, all.files = FALSE, full.names = FALSE)
            pfile <- file.path(lib, a, "Meta", "package.rds")
            ans <- c(ans, a[file.exists(pfile)])
        }
        return(unique(ans))
    } ## else
    s <- search()
    return(invisible(substring(s[substr(s, 1L, 8L) == "package:"], 9)))
}

path.package <-
function(package = NULL, quiet = FALSE)
{
    if(is.null(package)) package <- .packages()
    if(length(package) == 0L) return(character())
    s <- search()
    searchpaths <-
        lapply(seq_along(s), function(i) attr(as.environment(i), "path"))
    searchpaths[[length(s)]] <- system.file()
    pkgs <- paste("package", package, sep = ":")
    pos <- match(pkgs, s)
    if(any(m <- is.na(pos))) {
        if(!quiet) {
            if(all(m))
                stop("none of the packages are loaded")
            else
                warning(sprintf(ngettext(as.integer(sum(m)),
                                         "package %s is not loaded",
                                         "packages %s are not loaded"),
                                paste(package[m], collapse=", ")),
                        domain = NA)
        }
        pos <- pos[!m]
    }
    unlist(searchpaths[pos], use.names = FALSE)
}

## As from 2.9.0 ignore versioned installs
find.package <-
function(package = NULL, lib.loc = NULL, quiet = FALSE,
         verbose = getOption("verbose"))
{
    if(is.null(package) && is.null(lib.loc) && !verbose) {
        ## We only want the paths to the attached packages.
        return(path.package())
    }

    ## don't waste time looking for the standard packages:
    ## we know where they are and this can take a significant
    ## time with 1000+ packages installed.
    if(length(package) == 1L  &&
       package %in% c("base", "tools", "utils", "grDevices", "graphics",
                      "stats", "datasets", "methods", "grid", "parallel",
                      "splines", "stats4", "tcltk"))
        return(file.path(.Library, package))

    use_loaded <- FALSE
    if(is.null(package)) package <- .packages()
    if(is.null(lib.loc)) {
        use_loaded <- TRUE
        lib.loc <- .libPaths()
    }

    if(!length(package)) return(character())

    bad <- character()
    out <- character()

    for(pkg in package) {
        paths <- character()
        for(lib in lib.loc) {
            dirs <- list.files(lib,
                               pattern = paste0("^", pkg, "$"),
                               full.names = TRUE)
            ## Note that we cannot use tools::file_test() here, as
            ## cyclic namespace dependencies are not supported.  Argh.
            paths <- c(paths,
                       dirs[file.info(dirs)$isdir &
                            file.exists(file.path(dirs,
                                                  "DESCRIPTION"))])
        }
        if(use_loaded && pkg %in% loadedNamespaces()) {
            dir <- if (pkg == "base") system.file()
            else getNamespaceInfo(pkg, "path")
            paths <- c(dir, paths)
        }
        ## trapdoor for tools:::setRlibs
        if(length(paths) &&
           file.exists(file.path(paths[1], "dummy_for_check"))) {
            bad <- c(bad, pkg)
            next
        }
        if(length(paths)) {
            paths <- unique(paths)
            valid_package_version_regexp <-
                .standard_regexps()$valid_package_version
            db <- lapply(paths, function(p) {
                ## Note that this is sometimes used for source
                ## packages, e.g. by promptPackage from package.skeleton
                pfile <- file.path(p, "Meta", "package.rds")
                info <- if(file.exists(pfile))
                    ## this must have these fields to get installed
                    readRDS(pfile)$DESCRIPTION[c("Package", "Version")]
                else {
                    info <- tryCatch(read.dcf(file.path(p, "DESCRIPTION"),
                                              c("Package", "Version"))[1, ],
                                     error = identity)
                    if(inherits(info, "error")
                       || (length(info) != 2L)
                       || any(is.na(info)))
                        c(Package = NA, Version = NA) # need dimnames below
                    else
                        info
                }
            })
            db <- do.call("rbind", db)
            ok <- (apply(!is.na(db), 1L, all)
                   & (db[, "Package"] == pkg)
                   & (grepl(valid_package_version_regexp, db[, "Version"])))
            paths <- paths[ok]
        }

        if(length(paths) == 0L) {
            bad <- c(bad, pkg)
            next
        }
        if(length(paths) > 1L) {
            ## If a package was found more than once ...
            paths <- paths[1L]
            if(verbose)
                warning(gettextf("package %s found more than once,\nusing the one found in %s",
                                 sQuote(pkg), sQuote(paths)), domain = NA)
        }
        out <- c(out, paths)
    }

    if(!quiet && length(bad)) {
        if(length(out) == 0L) {
            if(length(bad) == 1L) {
                stop(gettextf("there is no package called %s", sQuote(pkg)),
                     domain = NA)
            } else {
                stop(ngettext(length(bad),
                              "there is no package called",
                              "there are no packages called"), " ",
                     paste(sQuote(bad), collapse = ", "), domain = NA)

            }
        }
        for(pkg in bad)
            warning(gettextf("there is no package called %s", sQuote(pkg)),
                    domain = NA)
    }

    out
}

format.packageInfo <-
function(x, ...)
{
    if(!inherits(x, "packageInfo")) stop("wrong class")
    vignetteMsg <-
        gettextf("Further information is available in the following vignettes in directory %s:",
                 sQuote(file.path(x$path, "doc")))
    headers <- sprintf("\n%s\n",
                       c(gettext("Description:"),
                         gettext("Index:"),
                         paste(strwrap(vignetteMsg), collapse = "\n")))
    formatDocEntry <- function(entry) {
        if(is.list(entry) || is.matrix(entry))
            formatDL(entry, style = "list")
        else
            entry
    }
    c(gettextf("\n\t\tInformation on package %s", sQuote(x$name)),
      unlist(lapply(which(!vapply(x$info, is.null, NA)),
                    function(i)
                        c(headers[i], formatDocEntry(x$info[[i]])))))

}

print.packageInfo <-
function(x, ...)
{
    outFile <- tempfile("RpackageInfo")
    writeLines(format(x), outFile)
    file.show(outFile, delete.file = TRUE,
              title =
              gettextf("Documentation for package %s", sQuote(x$name)))
    invisible(x)
}

.getRequiredPackages <-
function(file="DESCRIPTION", lib.loc = NULL, quietly = FALSE, useImports = FALSE)
{
    ## OK to call tools as only used during installation.
    pkgInfo <- tools:::.split_description(tools:::.read_description(file))
    .getRequiredPackages2(pkgInfo, quietly, lib.loc, useImports)
    invisible()
}

.getRequiredPackages2 <-
function(pkgInfo, quietly = FALSE, lib.loc = NULL, useImports = FALSE)
{
    pkgs <- unique(names(pkgInfo$Depends))
    if (length(pkgs)) {
        pkgname <- pkgInfo$DESCRIPTION["Package"]
        for(pkg in pkgs) {
            ## several packages 'Depends' on base!
            if (pkg == "base") next
            ## allow for multiple occurrences
            zs <- pkgInfo$Depends[names(pkgInfo$Depends) == pkg]
            have_vers <- any(vapply(zs, length, 1L) > 1L)
            if ( !paste("package", pkg, sep = ":") %in% search() ) {
                if (have_vers) {
                    pfile <- system.file("Meta", "package.rds",
                                         package = pkg, lib.loc = lib.loc)
                    if(!nzchar(pfile))
                        stop(gettextf("package %s required by %s could not be found",
                                      sQuote(pkg), sQuote(pkgname)),
                             call. = FALSE, domain = NA)
                    current <- readRDS(pfile)$DESCRIPTION["Version"]
                    for(z in zs)
                        if(length(z) > 1L) {
                            target <- as.numeric_version(z$version)
                            if (!do.call(z$op, list(as.numeric_version(current), target)))
##                            if (!eval(parse(text=paste("current", z$op, "target"))))
                                stop(gettextf("package %s %s was found, but %s %s is required by %s",
                                              sQuote(pkg), current, z$op,
                                              target, sQuote(pkgname)),
                                     call. = FALSE, domain = NA)
                        }
                }

                if (!quietly)
                    packageStartupMessage(gettextf("Loading required package: %s",
                                     pkg), domain = NA)
                library(pkg, character.only = TRUE, logical.return = TRUE,
                        lib.loc = lib.loc) ||
                stop(gettextf("package %s could not be loaded", sQuote(pkg)),
                     call. = FALSE, domain = NA)
            } else {
                ## check the required version number, if any
                if (have_vers) {
                    pfile <- system.file("Meta", "package.rds",
                                         package = pkg, lib.loc = lib.loc)
                    current <- readRDS(pfile)$DESCRIPTION["Version"]
                    for(z in zs)
                        if (length(z) > 1L) {
                            target <- as.numeric_version(z$version)
                            if (!do.call(z$op, list(as.numeric_version(current), target)))
##                            if (!eval(parse(text=paste("current", z$op, "target"))))
                                stop(gettextf("package %s %s is loaded, but %s %s is required by %s",
                                              sQuote(pkg), current, z$op,
                                              target, sQuote(pkgname)),
                                     call. = FALSE, domain = NA)
                        }
                }
            }
        }
    }
    if(useImports) {
        nss <- names(pkgInfo$Imports)
        for(ns in nss) loadNamespace(ns, lib.loc)
    }
}

.expand_R_libs_env_var <-
function(x)
{
    v <- paste(R.version[c("major", "minor")], collapse = ".")

    expand <- function(x, spec, expansion)
        gsub(paste0("(^|[^%])(%%)*%", spec),
             sprintf("\\1\\2%s", expansion), x)

    ## %V => version x.y.z
    x <- expand(x, "V", v)
    ## %v => version x.y
    x <- expand(x, "v", sub("\\.[^.]*$", "", v))
    ## %p => platform
    x <- expand(x, "p", R.version$platform)
    ## %a => arch
    x <- expand(x, "a", R.version$arch)
    ## %o => os
    x <- expand(x, "o", R.version$os)

    gsub("%%", "%", x)
}
#  File src/library/base/R/license.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

licence <- license <- function() {
    cat("\nThis software is distributed under the terms of the GNU General\n")
    cat("Public License, either Version 2, June 1991 or Version 3, June 2007.\n")
    cat("The terms of version 2 of the license are in a file called COPYING\nwhich you should have received with\n")
    cat("this software and which can be displayed by RShowDoc(\"COPYING\").\n")
    cat("Version 3 of the license can be displayed by RShowDoc(\"GPL-3\").\n")
    cat("\n")
    cat("Copies of both versions 2 and 3 of the license can be found\n")
    cat("at http://www.R-project.org/Licenses/.\n")
    cat("\n")
    cat("A small number of files (the API header files listed in\n")
    cat("R_DOC_DIR/COPYRIGHTS) are distributed under the\n")
    cat("LESSER GNU GENERAL PUBLIC LICENSE, version 2.1 or later.\n")
    cat("This can be displayed by RShowDoc(\"LGPL-2.1\"),\n")
    cat("or obtained at the URI given.\n")
    cat("Version 3 of the license can be displayed by RShowDoc(\"LGPL-3\").\n")
    cat("\n")
    cat("'Share and Enjoy.'\n\n")
}
#  File src/library/base/R/load.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

load <- function (file, envir = parent.frame(), verbose = FALSE)
{
    if (is.character(file)) {
        ## files are allowed to be of an earlier format
        ## gzfile can open gzip, bzip2, xz and uncompressed files.
        con <- gzfile(file)
        on.exit(close(con))
        ## Since the connection is not open this opens it in binary mode
        ## and closes it again.
        magic <- readChar(con, 5L, useBytes = TRUE)
	if (!length(magic)) stop("empty (zero-byte) input file")
	if (!grepl("RD[AX]2\n", magic)) {
            ## a check while we still know the call to load()
            if(grepl("RD[ABX][12]\r", magic))
                stop("input has been corrupted, with LF replaced by CR")
            ## Not a version 2 magic number, so try the pre-R-1.4.0 code
            warning(sprintf("file %s has magic number '%s'\n",
                            sQuote(basename(file)),
                            gsub("[\n\r]*", "", magic)),
                    "  ",
                    "Use of save versions prior to 2 is deprecated",
                    domain = NA, call. = FALSE)
            return(.Internal(load(file, envir)))
        }
    } else if (inherits(file, "connection")) {
        con <- if(inherits(file, "gzfile") || inherits(file, "gzcon")) file
               else gzcon(file)
    } else stop("bad 'file' argument")

    if (verbose)
    	cat("Loading objects:\n")
    	
    .Internal(loadFromConn2(con, envir, verbose))
}

save <- function(..., list = character(),
                 file = stop("'file' must be specified"),
                 ascii = FALSE, version = NULL, envir = parent.frame(),
                 compress = !ascii, compression_level,
                 eval.promises = TRUE, precheck = TRUE)
{
    opts <- getOption("save.defaults")
    if (missing(compress) && ! is.null(opts$compress))
        compress <- opts$compress
    if (missing(ascii) && ! is.null(opts$ascii))
        ascii <- opts$ascii
    if (missing(version)) version <- opts$version
    if (!is.null(version) && version < 2)
        warning("Use of save versions prior to 2 is deprecated", domain = NA)

    if(missing(list) && !length(list(...)))
	warning("nothing specified to be save()d")
    names <- as.character( substitute(list(...)))[-1L]
    list <- c(list, names)
    if (!is.null(version) && version == 1)
        .Internal(save(list, file, ascii, version, envir, eval.promises))
    else {
        if (precheck) {
            ## check for existence of objects before opening connection
            ## (and e.g. clobering file)
            ok <- unlist(lapply(list, exists, envir=envir))
            if(!all(ok)) {
                n <- sum(!ok)
                stop(sprintf(ngettext(n,
                                      "object %s not found",
                                      "objects %s not found"
                                      ),
                             paste(sQuote(list[!ok]), collapse = ", ")
                             ), domain = NA)
            }
        }
        if (is.character(file)) {
	    if(!nzchar(file)) stop("'file' must be non-empty string")
	    if(!is.character(compress)) {
		if(!is.logical(compress))
		    stop("'compress' must be logical or character")
		compress <- if(compress) "gzip" else "no compression"
	    }
	    con <- switch(compress,
			  "bzip2" = {
			      if (!missing(compression_level))
				  bzfile(file, "wb", compression = compression_level)
			      else bzfile(file, "wb")
			  }, "xz" = {
			      if (!missing(compression_level))
				  xzfile(file, "wb", compression = compression_level)
			      else xzfile(file, "wb", compression = 9)
			  }, "gzip" = {
			      if (!missing(compression_level))
				  gzfile(file, "wb", compression = compression_level)
			      else gzfile(file, "wb")
			  },
			  "no compression" = file(file, "wb"),

			  ## otherwise:
			  stop(gettextf("'compress = \"%s\"' is invalid", compress)))
	    on.exit(close(con))
	}
	else if (inherits(file, "connection"))
	    con <- file
	else stop("bad file argument")
	if(isOpen(con) && summary(con)$text != "binary")
	    stop("can only save to a binary connection")
	.Internal(saveToConn(list, con, ascii, version, envir, eval.promises))
    }
}

save.image <- function (file = ".RData", version = NULL, ascii = FALSE,
                        compress = !ascii, safe = TRUE)
{
    if (! is.character(file) || file == "")
        stop("'file' must be non-empty string")

    opts <- getOption("save.image.defaults")
    if(is.null(opts)) opts <- getOption("save.defaults")

    if (missing(safe) && ! is.null(opts$safe))
        safe <- opts$safe
    if (missing(ascii) && ! is.null(opts$ascii))
        ascii <- opts$ascii
    if (missing(compress) && ! is.null(opts$compress))
        compress <- opts$compress
    if (missing(version)) version <- opts$version

    if (safe) {
        ## find a temporary file name in the same directory so we can
        ## rename it to the final output file on success
        outfile <- paste0(file, "Tmp")
        i <- 0
        while (file.exists(outfile)) {
            i <- i + 1
            outfile <- paste0(file, "Tmp", i)
        }
    }
    else outfile <- file

    on.exit(file.remove(outfile))
    save(list = ls(envir = .GlobalEnv, all.names = TRUE), file = outfile,
         version = version, ascii = ascii, compress = compress,
         envir = .GlobalEnv, precheck = FALSE)
    if (safe)
        if (! file.rename(outfile, file)) {
            on.exit()
            stop(gettextf("image could not be renamed and is left in %s",
                          outfile), domain = NA)
        }
    on.exit()
}

sys.load.image <- function(name, quiet)
{
    if (file.exists(name)) {
        load(name, envir = .GlobalEnv)
        if (! quiet)
	    message("[Previously saved workspace restored]", "\n")
    }
}

sys.save.image <- function(name)
{
    ## Ensure that there is a reasonable chance that we can open a
    ## connection.
    closeAllConnections()
    save.image(name)
}

findPackageEnv <- function(info)
{
    if(info %in% search()) return(as.environment(info))
    message(gettextf("Attempting to load the environment %s", sQuote(info)),
            domain = NA)
    pkg <- substr(info, 9L, 1000L)
    if(require(pkg, character.only=TRUE, quietly = TRUE))
        return(as.environment(info))
    message("Specified environment not found: using '.GlobalEnv' instead")
    .GlobalEnv
}
#  File src/library/base/R/locales.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

Sys.getlocale <- function(category = "LC_ALL")
{
    category <- match(category, c("LC_ALL", "LC_COLLATE", "LC_CTYPE",
                                  "LC_MONETARY", "LC_NUMERIC", "LC_TIME",
                                  "LC_MESSAGES", "LC_PAPER", "LC_MEASUREMENT"))
    if(is.na(category)) stop("invalid 'category' argument")
    .Internal(Sys.getlocale(category))
}

Sys.setlocale <- function(category = "LC_ALL", locale = "")
{
    category <- match(category, c("LC_ALL", "LC_COLLATE", "LC_CTYPE",
                                  "LC_MONETARY", "LC_NUMERIC", "LC_TIME",
                                  "LC_MESSAGES", "LC_PAPER", "LC_MEASUREMENT"))
    if(is.na(category)) stop("invalid 'category' argument")
    .Internal(Sys.setlocale(category, locale))
}

Sys.localeconv <- function() .Internal(Sys.localeconv())
#  File src/library/base/R/lower.tri.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

lower.tri <- function(x, diag = FALSE)
{
    x <- as.matrix(x)
    if(diag) row(x) >= col(x)
    else row(x) > col(x)
}
#  File src/library/base/R/mapply.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

mapply <- function(FUN,..., MoreArgs = NULL, SIMPLIFY = TRUE, USE.NAMES = TRUE)
{
    FUN <- match.fun(FUN)
    dots <- list(...)

    answer <- .Internal(mapply(FUN, dots, MoreArgs))

    if (USE.NAMES && length(dots)) {
	if (is.null(names1 <- names(dots[[1L]])) && is.character(dots[[1L]]))
	    names(answer) <- dots[[1L]]
	else if (!is.null(names1))
	    names(answer) <- names1
    }
    if(!identical(SIMPLIFY, FALSE) && length(answer))
	simplify2array(answer, higher = (SIMPLIFY == "array"))
    else answer
}

.mapply <- function(FUN, dots, MoreArgs)
    .Internal(mapply(FUN, dots, MoreArgs))

Vectorize <- function(FUN, vectorize.args = arg.names, SIMPLIFY = TRUE,
                      USE.NAMES = TRUE)
{
    arg.names <- as.list(formals(FUN))
    arg.names[["..."]] <- NULL
    arg.names <- names(arg.names)

    vectorize.args <- as.character(vectorize.args)

    if (!length(vectorize.args)) return(FUN)

    if (!all(vectorize.args %in% arg.names))
    	stop("must specify names of formal arguments for 'vectorize'")

    FUNV <- function() { ## will set the formals below
        args <- lapply(as.list(match.call())[-1L], eval, parent.frame())
        names <- if(is.null(names(args))) character(length(args))
        else names(args)
        dovec <- names %in% vectorize.args
        do.call("mapply", c(FUN = FUN,
                            args[dovec],
                            MoreArgs = list(args[!dovec]),
                            SIMPLIFY = SIMPLIFY,
                            USE.NAMES = USE.NAMES))
    }
    formals(FUNV) <- formals(FUN)
    FUNV
}
#  File src/library/base/R/match.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

match <- function(x, table, nomatch = NA_integer_, incomparables = NULL)
    .Internal(match(x, table, nomatch, incomparables))

match.call <-
    function(definition=NULL, call=sys.call(sys.parent()), expand.dots=TRUE)
    .Internal(match.call(definition,call,expand.dots))

pmatch <- function(x, table, nomatch = NA_integer_, duplicates.ok = FALSE)
    .Internal(pmatch(as.character(x), as.character(table), nomatch,
                     duplicates.ok))

`%in%`  <- function(x, table) match(x, table, nomatch = 0L) > 0L

match.arg <- function (arg, choices, several.ok = FALSE)
{
    if (missing(choices)) {
	formal.args <- formals(sys.function(sys.parent()))
	choices <- eval(formal.args[[deparse(substitute(arg))]])
    }
    if (is.null(arg)) return(choices[1L])
    else if(!is.character(arg))
	stop("'arg' must be NULL or a character vector")
    if (!several.ok) { # most important (default) case:
        ## the arg can be the whole of choices as a default argument.
        if(identical(arg, choices)) return(arg[1L])
        if(length(arg) > 1L) stop("'arg' must be of length 1")
    } else if(length(arg) == 0L) stop("'arg' must be of length >= 1")

    ## handle each element of arg separately
    i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
    if (all(i == 0L))
	stop(gettextf("'arg' should be one of %s",
                      paste(dQuote(choices), collapse = ", ")),
             domain = NA)
    i <- i[i > 0L]
    if (!several.ok && length(i) > 1)
        stop("there is more than one match in 'match.arg'")
    choices[i]
}

charmatch <- function(x, table, nomatch = NA_integer_)
    .Internal(charmatch(as.character(x), as.character(table), nomatch))

char.expand <- function(input, target, nomatch = stop("no match"))
{
    if(length(input) != 1L)
	stop("'input' must have length 1")
    if(!(is.character(input) && is.character(target)))
	stop("'input' and 'target' must be character vectors")
    y <- .Internal(charmatch(input, target, NA_integer_))
    if(any(is.na(y))) eval(nomatch)
    target[y]
}
#  File src/library/base/R/match.fun.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

### clean up FUN arguments to *apply, outer, sweep, etc.
### note that this grabs two levels back and is not designed
### to be called at top level
match.fun <- function (FUN, descend = TRUE)
{
    if ( is.function(FUN) )
        return(FUN)
    if (!(is.character(FUN) && length(FUN) == 1L || is.symbol(FUN))) {
        ## Substitute in parent
        FUN <- eval.parent(substitute(substitute(FUN)))
        if (!is.symbol(FUN))
            stop(gettextf("'%s' is not a function, character or symbol",
                          deparse(FUN)), domain = NA)
    }
    envir <- parent.frame(2)
    if( descend )
        FUN <- get(as.character(FUN), mode = "function", envir = envir)
    else {
        FUN <- get(as.character(FUN), mode = "any", envir = envir)
        if( !is.function(FUN) )
           stop(gettextf("found non-function '%s'", FUN), domain = NA)
    }
    return(FUN)
}
#  File src/library/base/R/matrix.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

matrix <- function(data=NA, nrow=1, ncol=1, byrow=FALSE, dimnames=NULL)
{
    ## avoid copying to strip attributes in simple cases
    if (is.object(data) || !is.atomic(data)) data <- as.vector(data)
    ## NB: the defaults are not really nrow=1, ncol=1: missing values
    ## are treated differently, using length(data).
    .Internal(matrix(data, nrow, ncol, byrow, dimnames,
                     missing(nrow), missing(ncol)))
}

nrow <- function(x) dim(x)[1L]
ncol <- function(x) dim(x)[2L]

NROW <- function(x) if(length(d <- dim(x)))      d[1L] else length(x)
NCOL <- function(x) if(length(d <- dim(x)) > 1L) d[2L] else 1L

rownames <- function(x, do.NULL = TRUE, prefix = "row")
{
    dn <- dimnames(x)
    if(!is.null(dn[[1L]]))
	dn[[1L]]
    else {
        nr <- NROW(x)
	if(do.NULL) NULL
        else if(nr > 0L) paste0(prefix, seq_len(nr))
        else character()
    }
}

`rownames<-` <- function(x, value)
{
    if(is.data.frame(x)) {
        row.names(x) <- value
    } else {
        dn <- dimnames(x)
        if(is.null(dn)) {
            if(is.null(value)) return(x)
            if((nd <- length(dim(x))) < 1L)
                stop("attempt to set 'rownames' on an object with no dimensions")
            dn <- vector("list", nd)
        }
        if(length(dn) < 1L)
            stop("attempt to set 'rownames' on an object with no dimensions")
        if(is.null(value)) dn[1L] <- list(NULL) else dn[[1L]] <- value
        dimnames(x) <- dn
    }
    x
}

colnames <- function(x, do.NULL = TRUE, prefix = "col")
{
    if(is.data.frame(x) && do.NULL)
	return(names(x))
    dn <- dimnames(x)
    if(!is.null(dn[[2L]]))
	dn[[2L]]
    else {
        nc <- NCOL(x)
	if(do.NULL) NULL
        else if(nc > 0L) paste0(prefix, seq_len(nc))
        else character()
    }
}

`colnames<-` <- function(x, value)
{
    if(is.data.frame(x)) {
        names(x) <- value
    } else {
        dn <- dimnames(x)
        if(is.null(dn)) {
            if(is.null(value)) return(x)
            if((nd <- length(dim(x))) < 2L)
                stop("attempt to set 'colnames' on an object with less than two dimensions")
            dn <- vector("list", nd)
        }
        if(length(dn) < 2L)
            stop("attempt to set 'colnames' on an object with less than two dimensions")
        if(is.null(value)) dn[2L] <- list(NULL) else dn[[2L]] <- value
        dimnames(x) <- dn
    }
    x
}

row <- function(x, as.factor=FALSE)
{
    if(as.factor) {
        labs <- rownames(x, do.NULL=FALSE, prefix="")
        res <- factor(.Internal(row(dim(x))), labels=labs)
        dim(res) <- dim(x)
        res
    } else .Internal(row(dim(x)))
}

col <- function(x, as.factor=FALSE)
{
    if(as.factor) {
        labs <- colnames(x, do.NULL=FALSE, prefix="")
        res <- factor(.Internal(col(dim(x))), labels=labs)
        dim(res) <- dim(x)
        res
    } else .Internal(col(dim(x)))
}

crossprod <- function(x, y=NULL) .Internal(crossprod(x,y))
tcrossprod <- function(x, y=NULL) .Internal(tcrossprod(x,y))

t <- function(x) UseMethod("t")
## t.default is <primitive>
t.data.frame <- function(x)
{
    x <- as.matrix(x)
    NextMethod("t")
}
## as.matrix  is in "as"
#  File src/library/base/R/max.col.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

max.col <- function(m, ties.method = c("random", "first", "last"))
{
    ties.method <- match.arg(ties.method)
    tieM <- which(ties.method == eval(formals()[["ties.method"]]))
    .Internal(max.col(as.matrix(m), tieM))
}
#  File src/library/base/R/mean.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

mean <- function(x, ...) UseMethod("mean")

mean.default <- function(x, trim = 0, na.rm = FALSE, ...)
{
    if(!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
        warning("argument is not numeric or logical: returning NA")
        return(NA_real_)
    }
    if (na.rm)
	x <- x[!is.na(x)]
    if(!is.numeric(trim) || length(trim) != 1L)
        stop("'trim' must be numeric of length one")
    n <- length(x)
    if(trim > 0 && n) {
	if(is.complex(x))
	    stop("trimmed means are not defined for complex data")
        if(any(is.na(x))) return(NA_real_)
	if(trim >= 0.5) return(stats::median(x, na.rm=FALSE))
	lo <- floor(n*trim)+1
	hi <- n+1-lo
	x <- sort.int(x, partial = unique(c(lo, hi)))[lo:hi]
    }
    .Internal(mean(x))
}
#  File src/library/base/R/merge.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

merge <- function(x, y, ...) UseMethod("merge")

merge.default <- function(x, y, ...)
    merge(as.data.frame(x), as.data.frame(y), ...)

merge.data.frame <-
    function(x, y, by = intersect(names(x), names(y)), by.x = by, by.y = by,
             all = FALSE, all.x = all, all.y = all,
             sort = TRUE, suffixes = c(".x",".y"), incomparables = NULL,
             ...)
{
    fix.by <- function(by, df)
    {
        ## fix up 'by' to be a valid set of cols by number: 0 is row.names
        if(is.null(by)) by <- numeric()
        by <- as.vector(by)
        nc <- ncol(df)
        if(is.character(by)) {
            poss <- c("row.names", names(df))
            # names(df) are not necessarily unique, so check for multiple matches.
            if(any(bad <- !charmatch(by, poss, 0L)))
                stop(ngettext(sum(bad),
                              "'by' must specify a uniquely valid column",
                              "'by' must specify uniquely valid columns"),
                     domain = NA)
            by <- match(by, poss) - 1L
        } else if(is.numeric(by)) {
            if(any(by < 0L) || any(by > nc))
                stop("'by' must match numbers of columns")
        } else if(is.logical(by)) {
            if(length(by) != nc) stop("'by' must match number of columns")
            by <- seq_along(by)[by]
        } else stop("'by' must specify one or more columns as numbers, names or logical")
        if(any(bad <- is.na(by)))
            stop(ngettext(sum(bad),
                          "'by' must specify a uniquely valid column",
                          "'by' must specify uniquely valid columns"),
                 domain = NA)
         unique(by)
    }

    nx <- nrow(x <- as.data.frame(x)); ny <- nrow(y <- as.data.frame(y))
    by.x <- fix.by(by.x, x)
    by.y <- fix.by(by.y, y)
    if((l.b <- length(by.x)) != length(by.y))
        stop("'by.x' and 'by.y' specify different numbers of columns")
    if(l.b == 0L) {
        ## return the cartesian product of x and y, fixing up common names
        nm <- nm.x <- names(x)
        nm.y <- names(y)
        has.common.nms <- any(cnm <- nm.x %in% nm.y)
        if(has.common.nms) {
            names(x)[cnm] <- paste0(nm.x[cnm], suffixes[1L])
            cnm <- nm.y %in% nm
            names(y)[cnm] <- paste0(nm.y[cnm], suffixes[2L])
        }
        if (nx == 0L || ny == 0L) {
            res <- cbind(x[FALSE, ], y[FALSE, ])
        } else {
            ij <- expand.grid(seq_len(nx), seq_len(ny))
            res <- cbind(x[ij[, 1L], , drop = FALSE], y[ij[, 2L], , drop = FALSE])
        }
    }
    else {
        if(any(by.x == 0L)) {
            x <- cbind(Row.names = I(row.names(x)), x)
            by.x <- by.x + 1L
        }
        if(any(by.y == 0L)) {
            y <- cbind(Row.names = I(row.names(y)), y)
            by.y <- by.y + 1L
        }
        row.names(x) <- NULL
        row.names(y) <- NULL
        ## create keys from 'by' columns:
        if(l.b == 1L) {                  # (be faster)
            bx <- x[, by.x]; if(is.factor(bx)) bx <- as.character(bx)
            by <- y[, by.y]; if(is.factor(by)) by <- as.character(by)
        } else {
            ## Do these together for consistency in as.character.
            ## Use same set of names.
            bx <- x[, by.x, drop=FALSE]; by <- y[, by.y, drop=FALSE]
            names(bx) <- names(by) <- paste0("V", seq_len(ncol(bx)))
            bz <- do.call("paste", c(rbind(bx, by), sep = "\r"))
            bx <- bz[seq_len(nx)]
            by <- bz[nx + seq_len(ny)]
        }
        comm <- match(bx, by, 0L)
        bxy <- bx[comm > 0L]             # the keys which are in both
        xinds <- match(bx, bxy, 0L, incomparables)
        yinds <- match(by, bxy, 0L, incomparables)
        if(nx > 0L && ny > 0L)
            m <- .Internal(merge(xinds, yinds, all.x, all.y))
        else
            m <- list(xi = integer(), yi = integer(),
                      x.alone = seq_len(nx), y.alone = seq_len(ny))
        nm <- nm.x <- names(x)[-by.x]
        nm.by <- names(x)[by.x]
        nm.y <- names(y)[-by.y]
        ncx <- ncol(x)
        if(all.x) all.x <- (nxx <- length(m$x.alone)) > 0L
        if(all.y) all.y <- (nyy <- length(m$y.alone)) > 0L
        lxy <- length(m$xi)             # == length(m$yi)
        ## x = [ by | x ] :
        has.common.nms <- any(cnm <- nm.x %in% nm.y)
        if(has.common.nms && nzchar(suffixes[1L]))
            nm.x[cnm] <- paste0(nm.x[cnm], suffixes[1L])
        x <- x[c(m$xi, if(all.x) m$x.alone),
               c(by.x, seq_len(ncx)[-by.x]), drop=FALSE]
        names(x) <- c(nm.by, nm.x)
        if(all.y) { ## add the 'y.alone' rows to x[]
            ## need to have factor levels extended as well -> using [cr]bind
            ya <- y[m$y.alone, by.y, drop = FALSE]
            names(ya) <- nm.by
            ## this used to use a logical matrix, but that was not good
            ## enough as x could be zero-row.
            ya <- cbind(ya, x[rep.int(NA_integer_, nyy), nm.x, drop=FALSE ])
            x <- rbind(x, ya)
        }
        ## y (w/o 'by'):
        if(has.common.nms && nzchar(suffixes[2L])) {
            cnm <- nm.y %in% nm
            nm.y[cnm] <- paste0(nm.y[cnm], suffixes[2L])
        }
        y <- y[c(m$yi, if(all.x) rep.int(1L, nxx), if(all.y) m$y.alone),
               -by.y, drop = FALSE]
        if(all.x) {
            zap <- (lxy+1L):(lxy+nxx)
            for(i in seq_along(y)) {
                ## do it this way to invoke methods for e.g. factor
                if(is.matrix(y[[1]])) y[[1]][zap, ] <- NA
                else is.na(y[[i]]) <- zap
            }
        }

        if(has.common.nms) names(y) <- nm.y
        nm <- c(names(x), names(y))
        if(any(d <- duplicated(nm)))
            if(sum(d) > 1L)
                warning("column names ",
                        paste(sQuote(nm[d]), collapse = ", "),
                        " are duplicated in the result", domain = NA)
            else
                warning("column name ", sQuote(nm[d]),
                        " is duplicated in the result", domain = NA)
        res <- cbind(x, y)

        if (sort)
            res <- res[if(all.x || all.y) ## does NOT work
                       do.call("order", x[, seq_len(l.b), drop = FALSE])
            else sort.list(bx[m$xi]),, drop = FALSE]
    }
    attr(res, "row.names") <- .set_row_names(nrow(res))
    res
}
#  File src/library/base/R/message.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

simpleMessage <-
function(message, call = NULL)
    structure(list(message = message, call = call),
              class = c("simpleMessage", "message", "condition"))

suppressMessages <-
function(expr)
    withCallingHandlers(expr,
                        message = function(c)
                        invokeRestart("muffleMessage"))

message <-
function(..., domain = NULL, appendLF = TRUE)
{
    args <- list(...)
    cond <- if (length(args) == 1L && inherits(args[[1L]], "condition")) {
        if(nargs() > 1L)
            warning("additional arguments ignored in message()")
        args[[1L]]
    } else {
        msg <- .makeMessage(..., domain=domain, appendLF = appendLF)
        call <- sys.call()
        simpleMessage(msg, call)
    }
    defaultHandler <- function(c) {
        ## Maybe use special connection here?
        cat(conditionMessage(c), file=stderr(), sep = "")
    }
    withRestarts({
        signalCondition(cond)
        ## We don't get to the default handler if the signal
        ## is handled with a non-local exit, e.g. by
        ## invoking the muffleMessage restart.
        defaultHandler(cond)
    }, muffleMessage = function() NULL)
    invisible()
}

## also used by warning() and stop()
.makeMessage <- function(..., domain = NULL, appendLF = FALSE)
 {
    args <- list(...)
    msg <- if(length(args)) {
        args <- lapply(list(...), as.character)
        if(is.null(domain) || !is.na(domain))
            args <- .Internal(gettext(domain, unlist(args)))
        paste(args, collapse = "")
    } else ""
    if(appendLF) paste0(msg, "\n") else msg
}

.packageStartupMessage <- function (message, call = NULL)
    structure(list(message = message, call = call),
              class = c("packageStartupMessage", "condition", "message",
              "simpleMessage"))

suppressPackageStartupMessages <- function (expr)
    withCallingHandlers(expr, packageStartupMessage=function(c)
                        invokeRestart("muffleMessage"))

packageStartupMessage <- function(..., domain = NULL, appendLF = TRUE)
{
    call <- sys.call()
    msg <- .makeMessage(..., domain=domain, appendLF = appendLF)
    message(.packageStartupMessage(msg, call))
}
#  File src/library/base/R/methodsSupport.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

trace <- function(what, tracer, exit, at, print, signature, where = topenv(parent.frame()), edit = FALSE)
{
    needsAttach <- nargs() > 1L && !.isMethodsDispatchOn()
    if(needsAttach) {
        ns <- try(loadNamespace("methods"))
        if(isNamespace(ns))
            message("(loaded the methods namespace)", domain = NA)
        else
            stop("tracing functions requires the 'methods' package, but unable to load the 'methods' namespace")
    }
    else if(nargs() == 1L)
        return(.primTrace(what))
    tState <- tracingState(FALSE)
    on.exit(tracingState(tState))
    ## now call the version in the methods package, to ensure we get
    ## the correct namespace (e.g., correct version of class())
    call <- sys.call()
    call[[1L]] <- quote(methods::.TraceWithMethods)
    call$where <- where
    value <- eval.parent(call)
    on.exit() ## no error
    tracingState(tState)
    value
}

untrace <- function(what, signature = NULL, where = topenv(parent.frame())) {
    ## NOTE: following test is TRUE after loadNamespace("methods") (even if not in search())
    MethodsDispatchOn <- .isMethodsDispatchOn()
    if(MethodsDispatchOn) {
        tState <- tracingState(FALSE)
        on.exit(tracingState(tState))
    }
    if(!MethodsDispatchOn)
        return(.primUntrace(what)) ## can't have called trace except in primitive form
    ## at this point we can believe that the methods namespace was successfully loaded
    ## now call the version in the methods package, to ensure we get
    ## the correct namespace (e.g., correct version of class())
    call <- sys.call()
    call[[1L]] <- quote(methods::.TraceWithMethods)
    call$where <- where
    call$untrace <- TRUE
    value <- eval.parent(call)
    on.exit() ## no error
    tracingState(tState)
    invisible(value)
}


tracingState <- function(on = NULL) .Internal(traceOnOff(on))


asS4 <- function(object, flag = TRUE, complete = TRUE)
    .Internal(setS4Object(object, flag, complete))

asS3 <- function(object, flag = TRUE, complete = TRUE)
    .Internal(setS4Object(object, !as.logical(flag), complete))


.doTrace <- function(expr, msg) {
    on <- tracingState(FALSE)	   # turn it off QUICKLY (via a .Internal)
    if(on) {
	on.exit(tracingState(TRUE)) # restore on exit, keep off during trace
	if(!missing(msg)) {
	    call <- deparse(sys.call(sys.parent(1L)))
	    if(length(call) > 1L)
		call <- paste(call[[1L]], "....")
	    cat("Tracing", call, msg, "\n")
	}
	exprObj <- substitute(expr)
	eval.parent(exprObj)
    }
    NULL
}
#  File src/library/base/R/mode.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

mode <- function(x) {
    if(is.expression(x)) return("expression")
    if(is.call(x))
	return(switch(deparse(x[[1L]])[1L],
		      "(" = "(",
		      ## otherwise
		      "call"))
    if(is.name(x)) "name" else
    switch(tx <- typeof(x),
	   double =, integer = "numeric", # 'real=' dropped, 2000/Jan/14
	   closure =, builtin =, special = "function",
	   ## otherwise
	   tx)
}

`mode<-` <- function(x, value)
{
    if (storage.mode(x) == value) return(x)
    if(is.factor(x)) stop("invalid to change the storage mode of a factor")
    mde <- paste0("as.",value)
    atr <- attributes(x)
    isSingle <- !is.null(attr(x, "Csingle"))
    setSingle <- value == "single"
    x <- eval(call(mde,x), parent.frame())
    attributes(x) <- atr
    ## this avoids one copy
    if(setSingle != isSingle)
        attr(x, "Csingle") <- if(setSingle) TRUE # else NULL
    x
}

storage.mode <- function(x)
    switch(tx <- typeof(x),
	   closure = , builtin = , special = "function",
	   ## otherwise
	   tx)

### storage.mode<- is primitive as from R 2.6.0
#  File src/library/base/R/namespace.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2013 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

## give the base namespace a table for registered methods
`.__S3MethodsTable__.` <- new.env(hash = TRUE, parent = baseenv())

## NOTA BENE:
##  1) This code should work also when methods is not yet loaded
##  2) We use  ':::' instead of '::' inside the code below, for efficiency only

getNamespace <- function(name) {
    ns <- .Internal(getRegisteredNamespace(as.name(name)))
    if (! is.null(ns)) ns
    else tryCatch(loadNamespace(name), error = function(e) stop(e))
}

.getNamespace <- function(name) .Internal(getRegisteredNamespace(as.name(name)))

..getNamespace <- function(name, where) {
    ns <- .Internal(getRegisteredNamespace(as.name(name)))
    if (! is.null(ns)) ns
    else tryCatch(loadNamespace(name),
                  error = function(e) {
                      warning(gettextf("namespace %s is not available and has been replaced\nby .GlobalEnv when processing object %s",
                                       sQuote(name)[1L], sQuote(where)),
                              domain = NA, call. = FALSE, immediate. = TRUE)
                      .GlobalEnv
                  })
}

loadedNamespaces <- function()
    ls(.Internal(getNamespaceRegistry()), all.names = TRUE)

getNamespaceName <- function(ns) {
    ns <- asNamespace(ns)
    if (isBaseNamespace(ns)) "base"
    else getNamespaceInfo(ns, "spec")["name"]
}

getNamespaceVersion <- function(ns) {
    ns <- asNamespace(ns)
    if (isBaseNamespace(ns))
        c(version = paste(R.version$major, R.version$minor, sep = "."))
    else getNamespaceInfo(ns, "spec")["version"]
}

getNamespaceExports <- function(ns) {
    ns <- asNamespace(ns)
    if (isBaseNamespace(ns)) ls(.BaseNamespaceEnv, all.names = TRUE)
    else ls(getNamespaceInfo(ns, "exports"), all.names = TRUE)
}

getNamespaceImports <- function(ns) {
    ns <- asNamespace(ns)
    if (isBaseNamespace(ns)) NULL
    else getNamespaceInfo(ns, "imports")
}

getNamespaceUsers <- function(ns) {
    nsname <- getNamespaceName(asNamespace(ns))
    users <- character()
    for (n in loadedNamespaces()) {
        inames <- names(getNamespaceImports(n))
        if (match(nsname, inames, 0L))
            users <- c(n, users)
    }
    users
}

getExportedValue <- function(ns, name) {
    getInternalExportName <- function(name, ns) {
        exports <- getNamespaceInfo(ns, "exports")
        if (exists(name, envir = exports, inherits = FALSE))
            get(get(name, envir = exports, inherits = FALSE), envir = ns)
        else {
            ld <- getNamespaceInfo(ns, "lazydata")
            if (exists(name, envir = ld, inherits = FALSE))
                get(name, envir = ld, inherits = FALSE)
            else
                stop(gettextf("'%s' is not an exported object from 'namespace:%s'",
                              name, getNamespaceName(ns)),
                     call. = FALSE, domain = NA)
        }
    }
    ns <- asNamespace(ns)
    if (isBaseNamespace(ns)) get(name, envir = ns, inherits = FALSE)
    else getInternalExportName(name, ns)
}

`::` <- function(pkg, name) {
    pkg <- as.character(substitute(pkg))
    name <- as.character(substitute(name))
    getExportedValue(pkg, name)
}

`:::` <- function(pkg, name) {
    pkg <- as.character(substitute(pkg))
    name <- as.character(substitute(name))
    get(name, envir = asNamespace(pkg), inherits = FALSE)
}


attachNamespace <- function(ns, pos = 2L, depends = NULL)
{
    ## only used to run .onAttach
    runHook <- function(hookname, env, libname, pkgname) {
        if (exists(hookname, envir = env, inherits = FALSE)) {
            fun <- get(hookname, envir = env, inherits = FALSE)
            res <- tryCatch(fun(libname, pkgname), error = identity)
            if (inherits(res, "error")) {
                stop(gettextf("%s failed in %s() for '%s', details:\n  call: %s\n  error: %s",
                              hookname, "attachNamespace", nsname,
                              deparse(conditionCall(res))[1L],
                              conditionMessage(res)),
                     call. = FALSE, domain = NA)
            }
        } else if (exists(".First.lib", envir = env, inherits = FALSE) &&
                   nsname == Sys.getenv("R_INSTALL_PKG"))
            warning(sprintf("ignoring .First.lib() for package %s",
                            sQuote(nsname)), domain = NA, call. = FALSE)
    }
    runUserHook <- function(pkgname, pkgpath) {
        hook <- getHook(packageEvent(pkgname, "attach")) # might be list()
        for(fun in hook) try(fun(pkgname, pkgpath))
    }

    ns <- asNamespace(ns, base.OK = FALSE)
    nsname <- getNamespaceName(ns)
    nspath <- getNamespaceInfo(ns, "path")
    attname <- paste("package", nsname, sep = ":")
    if (attname %in% search())
        stop("namespace is already attached")
    env <- attach(NULL, pos = pos, name = attname)
    ## we do not want to run e.g. .onDetach here
    on.exit(.Internal(detach(pos)))
    attr(env, "path") <- nspath
    exports <- getNamespaceExports(ns)
    importIntoEnv(env, exports, ns, exports)
    ## always exists, might be empty
    dimpenv <- getNamespaceInfo(ns, "lazydata")
    dnames <- ls(dimpenv, all.names = TRUE)
    .Internal(importIntoEnv(env, dnames, dimpenv, dnames))
    if(length(depends)) assign(".Depends", depends, env)
    Sys.setenv("_R_NS_LOAD_" = nsname)
    on.exit(Sys.unsetenv("_R_NS_LOAD_"), add = TRUE)
    runHook(".onAttach", ns, dirname(nspath), nsname)
    lockEnvironment(env, TRUE)
    runUserHook(nsname, nspath)
    on.exit()
    Sys.unsetenv("_R_NS_LOAD_")
    invisible(env)
}

loadNamespace <- function (package, lib.loc = NULL,
                           keep.source = getOption("keep.source.pkgs"),
                           partial = FALSE, versionCheck = NULL)
{
    package <- as.character(package)[[1L]]

    ## check for cycles
    dynGet <- function(name,
                       notFound = stop(gettextf("%s not found", name),
                       domain = NA))
    {
        n <- sys.nframe()
        while (n > 1) {
            n <- n - 1
            env <- sys.frame(n)
            if (exists(name, envir = env, inherits = FALSE))
                return(get(name, envir = env, inherits = FALSE))
        }
        notFound
    }
    loading <- dynGet("__NameSpacesLoading__", NULL)
    if (match(package, loading, 0L))
        stop("cyclic namespace dependency detected when loading ",
             sQuote(package), ", already loading ",
             paste(sQuote(loading), collapse = ", "),
             domain = NA)
    "__NameSpacesLoading__" <- c(package, loading)

    ns <- .Internal(getRegisteredNamespace(as.name(package)))
    if (! is.null(ns)) {
        if(length(z <- versionCheck) == 3L) {
            current <- getNamespaceVersion(ns)
            if(!do.call(z$op, list(as.numeric_version(current), z$version)))
                stop(gettextf("namespace %s %s is already loaded, but %s %s is required",
                              sQuote(package), current, z$op, z$version),
                     domain = NA)
        }
        ns
    } else {
        ## only used here for .onLoad
        runHook <- function(hookname, env, libname, pkgname) {
            if (exists(hookname, envir = env, inherits = FALSE)) {
                fun <- get(hookname, envir = env, inherits = FALSE)
                res <- tryCatch(fun(libname, pkgname), error = identity)
                if (inherits(res, "error")) {
                    stop(gettextf("%s failed in %s() for '%s', details:\n  call: %s\n  error: %s",
                                  hookname, "loadNamespace", pkgname,
                                  deparse(conditionCall(res))[1L],
                                  conditionMessage(res)),
                         call. = FALSE, domain = NA)
                }
            }
        }
        runUserHook <- function(pkgname, pkgpath) {
            hooks <- getHook(packageEvent(pkgname, "onLoad")) # might be list()
            for(fun in hooks) try(fun(pkgname, pkgpath))
        }
        makeNamespace <- function(name, version = NULL, lib = NULL) {
            impenv <- new.env(parent = .BaseNamespaceEnv, hash = TRUE)
            attr(impenv, "name") <- paste("imports", name, sep = ":")
            env <- new.env(parent = impenv, hash = TRUE)
            name <- as.character(as.name(name))
            version <- as.character(version)
            info <- new.env(hash = TRUE, parent = baseenv())
            assign(".__NAMESPACE__.", info, envir = env)
            assign("spec", c(name = name, version = version), envir = info)
            setNamespaceInfo(env, "exports", new.env(hash = TRUE, parent = baseenv()))
            dimpenv <- new.env(parent = baseenv(), hash = TRUE)
            attr(dimpenv, "name") <- paste("lazydata", name, sep = ":")
            setNamespaceInfo(env, "lazydata", dimpenv)
            setNamespaceInfo(env, "imports", list("base" = TRUE))
            ## this should be an absolute path
            setNamespaceInfo(env, "path",
                             normalizePath(file.path(lib, name), "/", TRUE))
            setNamespaceInfo(env, "dynlibs", NULL)
            setNamespaceInfo(env, "S3methods", matrix(NA_character_, 0L, 3L))
            assign(".__S3MethodsTable__.",
                   new.env(hash = TRUE, parent = baseenv()),
                   envir = env)
            .Internal(registerNamespace(name, env))
            env
        }
        sealNamespace <- function(ns) {
            namespaceIsSealed <- function(ns)
               environmentIsLocked(ns)
            ns <- asNamespace(ns, base.OK = FALSE)
            if (namespaceIsSealed(ns))
                stop(gettextf("namespace %s is already sealed in 'loadNamespace'",
                              sQuote(getNamespaceName(ns))),
                     call. = FALSE, domain = NA)
            lockEnvironment(ns, TRUE)
            lockEnvironment(parent.env(ns), TRUE)
        }
        addNamespaceDynLibs <- function(ns, newlibs) {
            dynlibs <- getNamespaceInfo(ns, "dynlibs")
            setNamespaceInfo(ns, "dynlibs", c(dynlibs, newlibs))
        }

        bindTranslations <- function(pkgname, pkgpath)
        {
            ## standard packages are treated differently
            std <- c("compiler", "foreign", "grDevices", "graphics", "grid",
                     "methods", "parallel", "splines", "stats", "stats4",
                     "tcltk", "tools", "utils")
            popath <- if (pkgname %in% std) .popath else file.path(pkgpath, "po")
            if(!file.exists(popath)) return()
            bindtextdomain(pkgname, popath)
            bindtextdomain(paste("R", pkgname, sep = "-"), popath)
        }

        assignNativeRoutines <- function(dll, lib, env, nativeRoutines) {
            if(length(nativeRoutines) == 0L) return(NULL)

            if(nativeRoutines$useRegistration) {
                ## Use the registration information to register ALL the symbols
                fixes <- nativeRoutines$registrationFixes
                routines <- getDLLRegisteredRoutines.DLLInfo(dll, addNames = FALSE)
                lapply(routines,
                       function(type) {
                           lapply(type,
                                  function(sym) {
                                      varName <- paste0(fixes[1L], sym$name, fixes[2L])
                                      if(exists(varName, envir = env))
                                          warning(gettextf("failed to assign RegisteredNativeSymbol for %s to %s since %s is already defined in the %s namespace",
                                                           sym$name, varName, varName, sQuote(package)),
                                                  domain = NA)
                                      else
                                          assign(varName, sym, envir = env)
                                  })
                       })

            }

            symNames <- nativeRoutines$symbolNames
            if(length(symNames) == 0L) return(NULL)

            symbols <- getNativeSymbolInfo(symNames, dll, unlist = FALSE,
                                           withRegistrationInfo = TRUE)
            lapply(seq_along(symNames),
                   function(i) {
                       ## could vectorize this outside of the loop
                       ## and assign to different variable to
                       ## maintain the original names.
                       varName <- names(symNames)[i]
                       origVarName <- symNames[i]
                       if(exists(varName, envir = env))
                           if(origVarName != varName)
                               warning(gettextf("failed to assign NativeSymbolInfo for %s to %s since %s is already defined in the %s namespace",
                                                origVarName, varName, varName, sQuote(package)),
                                       domain = NA)
                           else
                               warning(gettextf("failed to assign NativeSymbolInfo for %s since %s is already defined in the %s namespace",
                                                origVarName, varName, sQuote(package)),
                                       domain = NA)
                       else
                           assign(varName, symbols[[origVarName]], envir = env)

                   })
            symbols
        }

        ## find package and check it has a namespace
        pkgpath <- find.package(package, lib.loc, quiet = TRUE)
        if (length(pkgpath) == 0L)
            stop(gettextf("there is no package called %s", sQuote(package)),
                 domain = NA)
        bindTranslations(package, pkgpath)
        package.lib <- dirname(pkgpath)
        package <- basename(pkgpath) # need the versioned name
        if (! packageHasNamespace(package, package.lib)) {
            hasNoNamespaceError <-
                function (package, package.lib, call = NULL) {
                class <- c("hasNoNamespaceError", "error", "condition")
                msg <- gettextf("package %s does not have a namespace",
                                sQuote(package))
                structure(list(message = msg, package = package,
                               package.lib = package.lib, call = call),
                          class = class)
            }
            stop(hasNoNamespaceError(package, package.lib))
        }

        ## create namespace; arrange to unregister on error
        ## Can we rely on the existence of R-ng 'nsInfo.rds' and
        ## 'package.rds'?
        ## No, not during builds of standard packages
        ## stats4 depends on methods, but exports do not matter
        ## whilst it is being built
        nsInfoFilePath <- file.path(pkgpath, "Meta", "nsInfo.rds")
        nsInfo <- if(file.exists(nsInfoFilePath)) readRDS(nsInfoFilePath)
        else parseNamespaceFile(package, package.lib, mustExist = FALSE)

        pkgInfoFP <- file.path(pkgpath, "Meta", "package.rds")
        if(file.exists(pkgInfoFP)) {
            pkgInfo <- readRDS(pkgInfoFP)
            version <- pkgInfo$DESCRIPTION["Version"]
            vI <- pkgInfo$Imports
            if(is.null(built <- pkgInfo$Built))
                stop(gettextf("package %s has not been installed properly\n",
                              sQuote(basename(pkgpath))),
                     call. = FALSE, domain = NA)
            R_version_built_under <- as.numeric_version(built$R)
            if(R_version_built_under < "3.0.0")
                stop(gettextf("package %s was built before R 3.0.0: please re-install it",
                             sQuote(basename(pkgpath))),
                     call. = FALSE, domain = NA)
            ## we need to ensure that S4 dispatch is on now if the package
            ## will require it, or the exports will be incomplete.
            dependsMethods <- "methods" %in% names(pkgInfo$Depends)
            if(dependsMethods) loadNamespace("methods")
            if(length(z <- versionCheck) == 3L &&
               !do.call(z$op, list(as.numeric_version(version), z$version)))
                stop(gettextf("namespace %s %s is being loaded, but %s %s is required",
                              sQuote(package), version, z$op, z$version),
                     domain = NA)
        }
        ns <- makeNamespace(package, version = version, lib = package.lib)
        on.exit(.Internal(unregisterNamespace(package)))

        ## process imports
        for (i in nsInfo$imports) {
            if (is.character(i))
                namespaceImport(ns,
                                loadNamespace(i, c(lib.loc, .libPaths()),
                                              versionCheck = vI[[i]]),
                                from = package)
            else
                namespaceImportFrom(ns,
                                    loadNamespace(j <- i[[1L]],
                                                  c(lib.loc, .libPaths()),
                                                  versionCheck = vI[[j]]),
                                    i[[2L]], from = package)
        }
        for(imp in nsInfo$importClasses)
            namespaceImportClasses(ns, loadNamespace(j <- imp[[1L]],
                                                     c(lib.loc, .libPaths()),
                                                     versionCheck = vI[[j]]),
                                   imp[[2L]], from = package)
        for(imp in nsInfo$importMethods)
            namespaceImportMethods(ns, loadNamespace(j <- imp[[1L]],
                                                     c(lib.loc, .libPaths()),
                                                     versionCheck = vI[[j]]),
                                   imp[[2L]], from = package)

        ## store info for loading namespace for loadingNamespaceInfo to read
        "__LoadingNamespaceInfo__" <- list(libname = package.lib,
                                           pkgname = package)

        env <- asNamespace(ns)
        ## save the package name in the environment
        assign(".packageName", package, envir = env)

        ## load the code
        codename <- strsplit(package, "_", fixed = TRUE)[[1L]][1L]
        codeFile <- file.path(pkgpath, "R", codename)
        if (file.exists(codeFile)) {
            res <- try(sys.source(codeFile, env, keep.source = keep.source))
            if(inherits(res, "try-error"))
                stop(gettextf("unable to load R code in package %s",
                              sQuote(package)), call. = FALSE, domain = NA)
        }
        # a package without R code currently is required to have a namespace
        # else warning(gettextf("package %s contains no R code",
        #                        sQuote(package)), call. = FALSE, domain = NA)

        ## partial loading stops at this point
        ## -- used in preparing for lazy-loading
        if (partial) return(ns)

        ## lazy-load any sysdata
        dbbase <- file.path(pkgpath, "R", "sysdata")
        if (file.exists(paste0(dbbase, ".rdb"))) lazyLoad(dbbase, env)

        ## load any lazydata into a separate environment
        dbbase <- file.path(pkgpath, "data", "Rdata")
        if(file.exists(paste0(dbbase, ".rdb")))
            lazyLoad(dbbase, getNamespaceInfo(ns, "lazydata"))

        ## register any S3 methods
        registerS3methods(nsInfo$S3methods, package, env)

        ## load any dynamic libraries
        dlls <- list()
        dynLibs <- nsInfo$dynlibs
        for (i in seq_along(dynLibs)) {
            lib <- dynLibs[i]
            dlls[[lib]]  <- library.dynam(lib, package, package.lib)
            assignNativeRoutines(dlls[[lib]], lib, env,
                                 nsInfo$nativeRoutines[[lib]])

            ## If the DLL has a name as in useDynLib(alias = foo),
            ## then assign DLL reference to alias.  Check if
            ## names() is NULL to handle case that the nsInfo.rds
            ## file was created before the names were added to the
            ## dynlibs vector.
            if(!is.null(names(nsInfo$dynlibs))
               && names(nsInfo$dynlibs)[i] != "")
                assign(names(nsInfo$dynlibs)[i], dlls[[lib]], envir = env)
            setNamespaceInfo(env, "DLLs", dlls)
        }
        addNamespaceDynLibs(env, nsInfo$dynlibs)


        ## used in e.g. utils::assignInNamespace
        Sys.setenv("_R_NS_LOAD_" = package)
        on.exit(Sys.unsetenv("_R_NS_LOAD_"), add = TRUE)
        ## run the load hook
        runHook(".onLoad", env, package.lib, package)

        ## process exports, seal, and clear on.exit action
        exports <- nsInfo$exports

        for (p in nsInfo$exportPatterns)
            exports <- c(ls(env, pattern = p, all.names = TRUE), exports)
        ##
        if(.isMethodsDispatchOn() && methods:::.hasS4MetaData(ns) &&
           !identical(package, "methods") ) {
            ## cache generics, classes in this namespace (but not methods itself,
            ## which pre-cached at install time
            methods:::cacheMetaData(ns, TRUE, ns)
            ## load actions may have added objects matching patterns
            for (p in nsInfo$exportPatterns) {
                expp <- ls(ns, pattern = p, all.names = TRUE)
                newEx <- !(expp %in% exports)
                if(any(newEx))
                    exports <- c(expp[newEx], exports)
            }
            ## process class definition objects
            expClasses <- nsInfo$exportClasses
            ##we take any pattern, but check to see if the matches are classes
            pClasses <- character()
            aClasses <- methods:::getClasses(ns)
            classPatterns <- nsInfo$exportClassPatterns
            ## defaults to exportPatterns
            if(!length(classPatterns))
                classPatterns <- nsInfo$exportPatterns
            for (p in classPatterns) {
                pClasses <- c(aClasses[grep(p, aClasses)], pClasses)
            }
            pClasses <- unique(pClasses)
            if( length(pClasses) ) {
                good <- vapply(pClasses, methods:::isClass, NA, where = ns)
                if( !any(good) && length(nsInfo$exportClassPatterns))
                    warning(gettextf("'exportClassPattern' specified in 'NAMESPACE' but no matching classes in package %s", sQuote(package)),
                            call. = FALSE, domain = NA)
                expClasses <- c(expClasses, pClasses[good])
            }
            if(length(expClasses)) {
                missingClasses <-
                    !vapply(expClasses, methods:::isClass, NA, where = ns)
                if(any(missingClasses))
                    stop(gettextf("in package %s classes %s were specified for export but not defined",
                                  sQuote(package),
                                  paste(expClasses[missingClasses],
                                        collapse = ", ")),
                         domain = NA)
                expClasses <- paste0(methods:::classMetaName(""), expClasses)
            }
            ## process methods metadata explicitly exported or
            ## implied by exporting the generic function.
            allGenerics <- unique(c(methods:::.getGenerics(ns),
                                   methods:::.getGenerics(parent.env(ns))))
            expMethods <- nsInfo$exportMethods
            ## check for generic functions corresponding to exported methods
            addGenerics <- expMethods[is.na(match(expMethods, exports))]
            if(length(addGenerics)) {
                nowhere <- sapply(addGenerics, function(what) !exists(what, mode = "function", envir = ns))
                if(any(nowhere)) {
                    warning(gettextf("no function found corresponding to methods exports from %s for: %s",
                                     sQuote(package),
                                     paste(sQuote(sort(unique(addGenerics[nowhere]))), collapse = ", ")),
                         domain = NA, call. = FALSE)
                    addGenerics <- addGenerics[!nowhere]
                }
                if(length(addGenerics)) {
                    ## skip primitives
                    addGenerics <- addGenerics[sapply(addGenerics, function(what) ! is.primitive(get(what, mode = "function", envir = ns)))]
                    ## the rest must be generic functions, implicit or local
                    ## or have been cached via a DEPENDS package
                    ok <- sapply(addGenerics, methods:::.findsGeneric, ns)
                    if(!all(ok)) {
                        bad <- sort(unique(addGenerics[!ok]))
                        msg <-
                            ngettext(length(bad),
                                     "Function found when exporting methods from the namespace %s which is not S4 generic: %s",
                                     "Functions found when exporting methods from the namespace %s which are not S4 generic: %s", domain = "R-base")
                        stop(sprintf(msg, sQuote(package),
                                     paste(sQuote(bad), collapse = ", ")),
                             domain = NA, call. = FALSE)
                    }
                    else if(any(ok > 1L))  #from the cache, don't add
                        addGenerics <- addGenerics[ok < 2L]
                }
### <note> Uncomment following to report any local generic functions
### that should have been exported explicitly.  But would be reported
### whenever the package is loaded, which is not when it is relevant.
### </note>
                ## local <- sapply(addGenerics, function(what) identical(as.character(get(what, envir = ns)@package), package))
                ## if(any(local))
                ##     message(gettextf("export(%s) from package %s generated by exportMethods()",
                ##        paste(addGenerics[local], collapse = ", ")),
                ##             domain = NA)
                exports <- c(exports, addGenerics)
            }
            expTables <- character()
            if(length(allGenerics)) {
                expMethods <-
                    unique(c(expMethods,
                             exports[!is.na(match(exports, allGenerics))]))
                missingMethods <- !(expMethods %in% allGenerics)
                if(any(missingMethods))
                    stop(gettextf("in %s methods for export not found: %s",
                                  sQuote(package),
                                  paste(expMethods[missingMethods],
                                        collapse = ", ")),
                         domain = NA)
                tPrefix <- methods:::.TableMetaPrefix()
                allMethodTables <-
                    unique(c(methods:::.getGenerics(ns, tPrefix),
                             methods:::.getGenerics(parent.env(ns), tPrefix)))
                needMethods <-
                    (exports %in% allGenerics) & !(exports %in% expMethods)
                if(any(needMethods))
                    expMethods <- c(expMethods, exports[needMethods])
                ## Primitives must have their methods exported as long
                ## as a global table is used in the C code to dispatch them:
                ## The following keeps the exported files consistent with
                ## the internal table.
                pm <- allGenerics[!(allGenerics %in% expMethods)]
                if(length(pm)) {
                    prim <- logical(length(pm))
                    for(i in seq_along(prim)) {
                        f <- methods:::getFunction(pm[[i]], FALSE, FALSE, ns)
                        prim[[i]] <- is.primitive(f)
                    }
                    expMethods <- c(expMethods, pm[prim])
                }
                for(i in seq_along(expMethods)) {
                    mi <- expMethods[[i]]
                    if(!(mi %in% exports) &&
                       exists(mi, envir = ns, mode = "function",
                              inherits = FALSE))
                        exports <- c(exports, mi)
                    pattern <- paste0(tPrefix, mi, ":")
                    ii <- grep(pattern, allMethodTables, fixed = TRUE)
                    if(length(ii)) {
			if(length(ii) > 1L) {
			    warning(gettextf("multiple methods tables found for %s",
				    sQuote(mi)), call. = FALSE, domain = NA)
			    ii <- ii[1L]
			}
                        expTables[[i]] <- allMethodTables[ii]
                     }
                    else { ## but not possible?
                      warning(gettextf("failed to find metadata object for %s",
                                       sQuote(mi)), call. = FALSE, domain = NA)
                    }
                }
            }
            else if(length(expMethods))
                stop(gettextf("in package %s methods %s were specified for export but not defined",
                              sQuote(package),
                              paste(expMethods, collapse = ", ")),
                     domain = NA)
            exports <- unique(c(exports, expClasses,  expTables))
        }
        ## certain things should never be exported.
        if (length(exports)) {
            stoplist <- c(".__NAMESPACE__.", ".__S3MethodsTable__.",
                          ".packageName", ".First.lib", ".onLoad",
                          ".onAttach", ".conflicts.OK", ".noGenerics")
            exports <- exports[! exports %in% stoplist]
        }
        namespaceExport(ns, exports)
        sealNamespace(ns)
        runUserHook(package, pkgpath)
        on.exit()
        Sys.unsetenv("_R_NS_LOAD_")
        ns
    }
}

## A version which returns TRUE/FALSE
requireNamespace <- function (package, ..., quietly = FALSE)
{
    package <- as.character(package)[[1L]] # like loadNamespace
    ns <- .Internal(getRegisteredNamespace(as.name(package)))
    res <- TRUE
    if (is.null(ns)) {
        packageStartupMessage(gettextf("Loading required namespace: %s",
                                       package), domain = NA)
        value <- tryCatch(loadNamespace(package, ...), error = function(e) e)
        if (inherits(value, "error")) {
            if (!quietly) {
                msg <- conditionMessage(value)
                cat("Failed with error:  ",
                    sQuote(msg), "\n", file = stderr(), sep = "")
                .Internal(printDeferredWarnings())
            }
            res <- FALSE
        }
    }
    invisible(res)
}

loadingNamespaceInfo <- function() {
    dynGet <- function(name,
                       notFound = stop(gettextf("%s not found", sQuote(name)),
                       domain = NA))
    {
        n <- sys.nframe()
        while (n > 1) {
            n <- n - 1
            env <- sys.frame(n)
            if (exists(name, envir = env, inherits = FALSE))
                return(get(name, envir = env, inherits = FALSE))
        }
        notFound
    }
    dynGet("__LoadingNamespaceInfo__", stop("not loading a namespace"))
}

topenv <- function(envir = parent.frame(),
                   matchThisEnv = getOption("topLevelEnvironment")) {
    while (! identical(envir, emptyenv())) {
        nm <- attributes(envir)[["names", exact = TRUE]]
        if ((is.character(nm) && length(grep("^package:" , nm))) ||
	    ## matchThisEnv is used in sys.source
            identical(envir, matchThisEnv) ||
            identical(envir, .GlobalEnv) ||
            identical(envir, baseenv()) ||
            .Internal(isNamespaceEnv(envir)) ||
	    ## packages except base and those with a separate namespace have .packageName
            exists(".packageName", envir = envir, inherits = FALSE))
            return(envir)
        else envir <- parent.env(envir)
    }
    return(.GlobalEnv)
}

unloadNamespace <- function(ns)
{
    ## only used to run .onUnload
    runHook <- function(hookname, env, ...) {
        if (exists(hookname, envir = env, inherits = FALSE)) {
            fun <- get(hookname, envir = env, inherits = FALSE)
            res <- tryCatch(fun(...), error=identity)
            if (inherits(res, "error")) {
                warning(gettextf("%s failed in %s() for '%s', details:\n  call: %s\n  error: %s",
                                 hookname, "unloadNamespace", nsname,
                                 deparse(conditionCall(res))[1L],
                                 conditionMessage(res)),
                        call. = FALSE, domain = NA)
            }
        }
    }
    ns <- asNamespace(ns, base.OK = FALSE)
    nsname <- getNamespaceName(ns)
    pos <- match(paste("package", nsname, sep = ":"), search())
    if (! is.na(pos)) detach(pos = pos)
    users <- getNamespaceUsers(ns)
    if (length(users))
        stop(gettextf("namespace %s is imported by %s so cannot be unloaded",
                      sQuote(getNamespaceName(ns)),
                      paste(sQuote(users), collapse = ", ")),
             domain = NA)
    nspath <- getNamespaceInfo(ns, "path")
    hook <- getHook(packageEvent(nsname, "onUnload")) # might be list()
    for(fun in rev(hook)) try(fun(nsname, nspath))
    runHook(".onUnload", ns, nspath)
    .Internal(unregisterNamespace(nsname))
    if(.isMethodsDispatchOn() && methods:::.hasS4MetaData(ns))
        methods:::cacheMetaData(ns, FALSE, ns)
    .Internal(lazyLoadDBflush(paste0(nspath, "/R/", nsname, ".rdb")))
    invisible()
}

isNamespace <- function(ns) .Internal(isNamespaceEnv(ns))

isBaseNamespace <- function(ns) identical(ns, .BaseNamespaceEnv)

getNamespaceInfo <- function(ns, which) {
    ns <- asNamespace(ns, base.OK = FALSE)
    info <- get(".__NAMESPACE__.", envir = ns, inherits = FALSE)
    get(which, envir = info, inherits = FALSE)
}

setNamespaceInfo <- function(ns, which, val) {
    ns <- asNamespace(ns, base.OK = FALSE)
    info <- get(".__NAMESPACE__.", envir = ns, inherits = FALSE)
    assign(which, val, envir = info)
}

asNamespace <- function(ns, base.OK = TRUE) {
    if (is.character(ns) || is.name(ns))
        ns <- getNamespace(ns)
    if (! isNamespace(ns))
        stop("not a namespace")
    else if (! base.OK && isBaseNamespace(ns))
        stop("operation not allowed on base namespace")
    else ns
}

namespaceImport <- function(self, ..., from = NULL)
    for (ns in list(...))
        namespaceImportFrom(self, asNamespace(ns), from = from)

namespaceImportFrom <- function(self, ns, vars, generics, packages, from = NULL)
{
    addImports <- function(ns, from, what) {
        imp <- structure(list(what), names = getNamespaceName(from))
        imports <- getNamespaceImports(ns)
        setNamespaceInfo(ns, "imports", c(imports, imp))
    }
    namespaceIsSealed <- function(ns)
       environmentIsLocked(ns)
    makeImportExportNames <- function(spec) {
        old <- as.character(spec)
        new <- names(spec)
        if (is.null(new)) new <- old
        else new[new == ""] <- old[new == ""]
        names(old) <- new
        old
    }
    whichMethodMetaNames <- function(impvars) {
        if(!.isMethodsDispatchOn())
            return(numeric())
        mm <- ".__T__"
        seq_along(impvars)[substr(impvars, 1L, nchar(mm, type = "c")) == mm]
    }
    if (is.character(self))
        self <- getNamespace(self)
    ns <- asNamespace(ns)
    nsname <- getNamespaceName(ns)
    impvars <- if (missing(vars)) {
        ## certain things should never be imported:
        ## but most of these are never exported (exception: .Last.lib)
        stoplist <- c(".__NAMESPACE__.", ".__S3MethodsTable__.",
                      ".packageName", ".First.lib", ".Last.lib",
                      ".onLoad", ".onAttach", ".onDetach",
                      ".conflicts.OK", ".noGenerics")
        vars <- getNamespaceExports(ns)
        vars <- vars[! vars %in% stoplist]
    } else vars
    impvars <- makeImportExportNames(impvars)
    impnames <- names(impvars)
    if (anyDuplicated(impnames)) {
        stop(gettextf("duplicate import names %s",
                      paste(sQuote(impnames[duplicated(impnames)]),
                            collapse = ", ")), domain = NA)
    }
    if (isNamespace(self) && isBaseNamespace(self)) {
        impenv <- self
        msg <- gettext("replacing local value with import %s when loading %s")
        register <- FALSE
    }
    else if (isNamespace(self)) {
        if (namespaceIsSealed(self))
            stop("cannot import into a sealed namespace")
        impenv <- parent.env(self)
        msg <- gettext("replacing previous import by %s when loading %s")
        register <- TRUE
    }
    else if (is.environment(self)) {
        impenv <- self
        msg <- gettext("replacing local value with import %s when loading %s")
        register <- FALSE
    }
    else stop("invalid import target")
    which <- whichMethodMetaNames(impvars)
    if(length(which)) {
	## If methods are already in impenv, merge and don't import
	delete <- integer()
	for(i in which) {
	    methodsTable <- .mergeImportMethods(impenv, ns, impvars[[i]])
	    if(is.null(methodsTable))
	    {} ## first encounter, just import it
	    else { ##
		delete <- c(delete, i)
		if(!missing(generics)) {
		    genName <- generics[[i]]
                    ## if(i > length(generics) || !nzchar(genName))
                    ##   {warning("got invalid index for importing ",mlname); next}
		    fdef <- methods:::getGeneric(genName,
                                                 where = impenv,
                                                 package = packages[[i]])
		    if(is.null(fdef))
			warning(gettextf("found methods to import for function %s but not the generic itself",
					 sQuote(genName)),
                                call. = FALSE, domain = NA)
		    else
			methods:::.updateMethodsInTable(fdef, ns, TRUE)
		}
	    }
	}
	if(length(delete)) {
	    impvars <- impvars[-delete]
	    impnames <- impnames[-delete]
	}
    }
    for (n in impnames)
	if (exists(n, envir = impenv, inherits = FALSE)) {
	    if (.isMethodsDispatchOn() && methods:::isGeneric(n, ns)) {
		## warn only if generic overwrites a function which
		## it was not derived from
		genNs <- get(n, envir = ns)@package
                genImp <- get(n, envir = impenv)
                if(methods::is(genImp, "genericFunction") &&
                   identical(genNs, genImp@package)) next # same generic
		genImpenv <- environmentName(environment(genImp))
                ## May call environment() on a non-function--an undocumented
                ## "feature" of environment() is that it returns a special
                ## attribute for non-functions, usually NULL
		if (!identical(genNs, genImpenv) ||
                    methods:::isGeneric(n, impenv)) {}
                else next
	    }
            ## this is always called from another function, so reporting call
            ## is unhelpful
            warning(sprintf(msg, sQuote(paste(nsname, n, sep = "::")),
                            sQuote(from)),
                    call. = FALSE, domain = NA)
	}
    importIntoEnv(impenv, impnames, ns, impvars)
    if (register)
        addImports(self, ns, if (missing(vars)) TRUE else impvars)
}

namespaceImportClasses <- function(self, ns, vars, from = NULL)
{
    for(i in seq_along(vars))
        vars[[i]] <- methods:::classMetaName(vars[[i]])
    namespaceImportFrom(self, asNamespace(ns), vars, from = from)
}

namespaceImportMethods <- function(self, ns, vars, from = NULL)
{
    allVars <- character()
    generics <- character()
    packages <- character()
    allFuns <- methods:::.getGenerics(ns) # all the methods tables in ns
    allPackages <- attr(allFuns, "package")
    pkg <- methods:::getPackageName(ns)
    if(!all(vars %in% allFuns)) {
        message(gettextf("No methods found in \"%s\" for requests: %s",
                         pkg, paste(vars[is.na(match(vars, allFuns))], collapse = ", ")),
                domain = NA)
        vars <- vars[vars %in% allFuns]
    }
    tPrefix <- methods:::.TableMetaPrefix()
    allMethodTables <- methods:::.getGenerics(ns, tPrefix)
    if(any(is.na(match(vars, allFuns))))
        stop(gettextf("requested methods not found in environment/package %s: %s",
                      sQuote(pkg),
                      paste(vars[is.na(match(vars, allFuns))],
                            collapse = ", ")), call. = FALSE, domain = NA)
    for(i in seq_along(allFuns)) {
        ## import methods tables if asked for
        ## or if the corresponding generic was imported
        g <- allFuns[[i]]
        p <- allPackages[[i]]
        if(exists(g, envir = self, inherits = FALSE) # already imported
           || g %in% vars) { # requested explicitly
            tbl <- methods:::.TableMetaName(g, p)
            if(is.null(.mergeImportMethods(self, ns, tbl))) { # a new methods table
               allVars <- c(allVars, tbl) # import it;else, was merged
               generics <- c(generics, g)
               packages <- c(packages, p)
            }
        }
        if(g %in% vars && !exists(g, envir = self, inherits = FALSE)) {
            if(exists(g, envir = ns) &&
               methods:::is(get(g, envir = ns), "genericFunction")) {
                allVars <- c(allVars, g)
                generics <- c(generics, g)
                packages <- c(packages, p)
            }
            else { # should be primitive
                fun <- methods::getFunction(g, mustFind = FALSE, where = self)
                if(is.primitive(fun) || methods::is(fun, "genericFunction")) {}
                else
                    warning(gettextf("No generic function found corresponding to requested imported methods for \"%s\" from package \"%s\" (malformed exports?)",
                                 g, pkg),
                        domain = NA)
            }
        }
    }
    namespaceImportFrom(self, asNamespace(ns), allVars, generics, packages,
                        from = from)
}

importIntoEnv <- function(impenv, impnames, expenv, expnames) {
    exports <- getNamespaceInfo(expenv, "exports")
    ex <- .Internal(ls(exports, TRUE))
    if(!all(expnames %in% ex)) {
        miss <- expnames[! expnames %in% ex]
        ## if called (indirectly) for namespaceImportClasses
        ## these are all classes
        if(all(grepl("^\\.__C__", miss))) {
            miss <- sub("^\\.__C__", "", miss)
            stop(sprintf(ngettext(length(miss),
                                  "class %s is not exported by 'namespace:%s'",
                                  "classes %s are not exported by 'namespace:%s'"),
                         paste(paste0('"', miss, '"'), collapse = ", "),
                         getNamespaceName(expenv)),
                 call. = FALSE, domain = NA)
        } else {
            stop(sprintf(ngettext(length(miss),
                                  "object %s is not exported by 'namespace:%s'",
                                  "objects %s are not exported by 'namespace:%s'"),
                         paste(sQuote(miss), collapse = ", "),
                         getNamespaceName(expenv)),
                 call. = FALSE, domain = NA)
        }
    }
    expnames <- unlist(lapply(expnames, get, envir = exports, inherits = FALSE))
    if (is.null(impnames)) impnames <- character()
    if (is.null(expnames)) expnames <- character()
    .Internal(importIntoEnv(impenv, impnames, expenv, expnames))
}

namespaceExport <- function(ns, vars) {
    namespaceIsSealed <- function(ns)
       environmentIsLocked(ns)
    if (namespaceIsSealed(ns))
        stop("cannot add to exports of a sealed namespace")
    ns <- asNamespace(ns, base.OK = FALSE)
    if (length(vars)) {
        addExports <- function(ns, new) {
            exports <- getNamespaceInfo(ns, "exports")
            expnames <- names(new)
            intnames <- new
            objs <- .Internal(ls(exports, TRUE))
            ex <- expnames %in% objs
            if(any(ex))
                warning(sprintf(ngettext(sum(ex),
                                         "previous export '%s' is being replaced",
                                         "previous exports '%s' are being replaced"),
                                paste(sQuote(expnames[ex]), collapse = ", ")),
                        call. = FALSE, domain = NA)
            for (i in seq_along(new))
                assign(expnames[i], intnames[i], envir = exports)
        }
        makeImportExportNames <- function(spec) {
            old <- as.character(spec)
            new <- names(spec)
            if (is.null(new)) new <- old
            else new[new == ""] <- old[new == ""]
            names(old) <- new
            old
        }
        new <- makeImportExportNames(unique(vars))
        ## calling exists each time is too slow, so do two phases
        undef <- new[! new %in% .Internal(ls(ns, TRUE))]
        undef <- undef[! vapply(undef, exists, NA, envir = ns)]
        if (length(undef)) {
            undef <- do.call("paste", as.list(c(undef, sep = ", ")))
            stop(gettextf("undefined exports: %s", undef), domain = NA)
        }
        if(.isMethodsDispatchOn()) .mergeExportMethods(new, ns)
        addExports(ns, new)
    }
}

.mergeExportMethods <- function(new, ns) {
##    if(!.isMethodsDispatchOn()) return(FALSE)
    mm <- methods:::methodsPackageMetaName("M","")
    newMethods <- new[substr(new, 1L, nchar(mm, type = "c")) == mm]
    nsimports <- parent.env(ns)
    for(what in newMethods) {
        if(exists(what, envir = nsimports, inherits = FALSE)) {
            m1 <- get(what, envir = nsimports)
            m2 <- get(what, envir = ns)
            assign(what, envir = ns, methods:::mergeMethods(m1, m2))
        }
    }
}


## NB this needs a decorated name, foo_ver, if appropriate
packageHasNamespace <- function(package, package.lib) {
    namespaceFilePath <- function(package, package.lib)
        file.path(package.lib, package, "NAMESPACE")
    file.exists(namespaceFilePath(package, package.lib))
}

parseNamespaceFile <- function(package, package.lib, mustExist = TRUE)
{
    namespaceFilePath <- function(package, package.lib)
        file.path(package.lib, package, "NAMESPACE")

    ## These two functions are essentially local to the parsing of
    ## the namespace file and don't need to be made available to
    ## users.  These manipulate the data from useDynLib() directives
    ## for the same DLL to determine how to map the symbols to R
    ## variables.

    nativeRoutineMap <-
        ## Creates a new NativeRoutineMap.
        function(useRegistration, symbolNames, fixes) {
            proto <- list(useRegistration = FALSE,
                          symbolNames = character())
            class(proto) <- "NativeRoutineMap"

            mergeNativeRoutineMaps(proto, useRegistration, symbolNames, fixes)
        }

    mergeNativeRoutineMaps <-
        ## Merges new settings into a NativeRoutineMap
        function(map, useRegistration, symbolNames, fixes) {
            if(!useRegistration)
                names(symbolNames) <-
                    paste0(fixes[1L],  names(symbolNames), fixes[2L])
            else
                map$registrationFixes <- fixes
            map$useRegistration <- map$useRegistration || useRegistration
            map$symbolNames <- c(map$symbolNames, symbolNames)
            map
        }

    nsFile <- namespaceFilePath(package, package.lib)
    descfile <- file.path(package.lib, package, "DESCRIPTION")
    enc <- if (file.exists(descfile)) {
        read.dcf(file = descfile, "Encoding")[1L]
    } else NA_character_
    if (file.exists(nsFile))
        directives <- if (!is.na(enc) &&
                          ! Sys.getlocale("LC_CTYPE") %in% c("C", "POSIX")) {
	    con <- file(nsFile, encoding=enc)
            on.exit(close(con))
	    parse(con, srcfile=NULL)
        } else parse(nsFile, srcfile=NULL)
    else if (mustExist)
        stop(gettextf("package %s has no 'NAMESPACE' file", sQuote(package)),
             domain = NA)
    else directives <- NULL
    exports <- character()
    exportPatterns <- character()
    exportClasses <- character()
    exportClassPatterns <- character()
    exportMethods <- character()
    imports <- list()
    importMethods <- list()
    importClasses <- list()
    dynlibs <- character()
    nS3methods <- 1000L
    S3methods <- matrix(NA_character_, nS3methods, 3L)
    nativeRoutines <- list()
    nS3 <- 0L
    parseDirective <- function(e) {
        ## trying to get more helpful error message:
	asChar <- function(cc) {
	    r <- as.character(cc)
	    if(any(r == ""))
		stop(gettextf("empty name in directive '%s' in 'NAMESPACE' file",
			      as.character(e[[1L]])),
		     domain = NA)
	    r
	}
        switch(as.character(e[[1L]]),
               "if" = if (eval(e[[2L]], .GlobalEnv))
               parseDirective(e[[3L]])
               else if (length(e) == 4L)
               parseDirective(e[[4L]]),
               "{" =  for (ee in as.list(e[-1L])) parseDirective(ee),
               "=" =,
               "<-" = {
                   parseDirective(e[[3L]])
                   if(as.character(e[[3L]][[1L]]) == "useDynLib")
                       names(dynlibs)[length(dynlibs)] <<- asChar(e[[2L]])
               },
               export = {
                   exp <- e[-1L]
                   exp <- structure(asChar(exp), names = names(exp))
                   exports <<- c(exports, exp)
               },
               exportPattern = {
                   pat <- asChar(e[-1L])
                   exportPatterns <<- c(pat, exportPatterns)
               },
               exportClassPattern = {
                   pat <- asChar(e[-1L])
                   exportClassPatterns <<- c(pat, exportClassPatterns)
               },
               exportClass = , exportClasses = {
                   exportClasses <<- c(asChar(e[-1L]), exportClasses)
               },
               exportMethods = {
                   exportMethods <<- c(asChar(e[-1L]), exportMethods)
               },
               import = imports <<- c(imports, as.list(asChar(e[-1L]))),
               importFrom = {
                   imp <- e[-1L]
                   ivars <- imp[-1L]
                   inames <- names(ivars)
                   imp <- list(asChar(imp[1L]),
                               structure(asChar(ivars), names = inames))
                   imports <<- c(imports, list(imp))
               },
               importClassFrom = , importClassesFrom = {
                   imp <- asChar(e[-1L])
                   pkg <- imp[[1L]]
                   impClasses <- imp[-1L]
                   imp <- list(asChar(pkg), asChar(impClasses))
                   importClasses <<- c(importClasses, list(imp))
               },
               importMethodsFrom = {
                   imp <- asChar(e[-1L])
                   pkg <- imp[[1L]]
                   impMethods <- imp[-1L]
                   imp <- list(asChar(pkg), asChar(impMethods))
                   importMethods <<- c(importMethods, list(imp))
               },
               useDynLib = {

                   ## This attempts to process as much of the
                   ## information as possible when NAMESPACE is parsed
                   ## rather than when it is loaded and creates
                   ## NativeRoutineMap objects to handle the mapping
                   ## of symbols to R variable names.

                   ## The name is the second element after useDynLib
                   dyl <- as.character(e[2L])
                   ## We ensure uniqueness at the end.
                   dynlibs <<-
                       structure(c(dynlibs, dyl),
                                 names = c(names(dynlibs),
                                 ifelse(!is.null(names(e)) &&
                                        names(e)[2L] != "", names(e)[2L], "" )))
                   if (length(e) > 2L) {
                       ## Author has specified some mappings for the symbols

                       symNames <- as.character(e[-c(1L, 2L)])
                       names(symNames) <- names(e[-c(1, 2)])

                       ## If there are no names, then use the names of
                       ## the symbols themselves.
                       if (length(names(symNames)) == 0L)
                           names(symNames) = symNames
                       else if (any(w <- names(symNames) == "")) {
                           names(symNames)[w] = symNames[w]
                       }

                       ## For each DLL, we build up a list the (R
                       ## variable name, symbol name) mappings. We do
                       ## this in a NativeRoutineMap object and we
                       ## merge potentially multiple useDynLib()
                       ## directives for the same DLL into a single
                       ## map.  Then we have separate NativeRoutineMap
                       ## for each different DLL.  E.g. if we have
                       ## useDynLib(foo, a, b, c) and useDynLib(bar,
                       ## a, x, y) we would maintain and resolve them
                       ## separately.

                       dup <- duplicated(names(symNames))
                       if (any(dup))
                           warning(gettextf("duplicate symbol names %s in useDynLib(\"%s\")",
                                            paste(sQuote(names(symNames)[dup]),
                                                  collapse = ", "), dyl),
                                   domain = NA)

                       symNames <- symNames[!dup]

                       ## Deal with any prefix/suffix pair.
                       fixes <- c("", "")
                       idx <- match(".fixes", names(symNames))
                       if(!is.na(idx)) {
                           ## Take .fixes and treat it as a call,
                           ## e.g. c("pre", "post") or a regular name
                           ## as the prefix.
                           if(symNames[idx] != "") {
                               e <- parse(text = symNames[idx], srcfile = NULL)[[1L]]
                               if(is.call(e))
                                   val <- eval(e)
                               else
                                   val <- as.character(e)
                               if(length(val))
                                   fixes[seq_along(val)] <- val
                           }
                           symNames <- symNames[-idx]
                       }

                       ## Deal with a .registration entry. It must be
                       ## .registration = value and value will be coerced
                       ## to a logical.
                       useRegistration <- FALSE
                       idx <- match(".registration", names(symNames))
                       if(!is.na(idx)) {
                           useRegistration <- as.logical(symNames[idx])
                           symNames <- symNames[-idx]
                       }

                       ## Now merge into the NativeRoutineMap.
                       nativeRoutines[[ dyl ]] <<-
                           if(dyl %in% names(nativeRoutines))
                               mergeNativeRoutineMaps(nativeRoutines[[ dyl ]],
                                                      useRegistration,
                                                      symNames, fixes)
                           else
                               nativeRoutineMap(useRegistration, symNames,
                                                fixes)
                   }
               },
               S3method = {
                   spec <- e[-1L]
                   if (length(spec) != 2L && length(spec) != 3L)
                       stop(gettextf("bad 'S3method' directive: %s",
                                     deparse(e)),
                            call. = FALSE, domain = NA)
                   nS3 <<- nS3 + 1L
                   if(nS3 > nS3methods) {
                       old <- S3methods
                       nold <- nS3methods
                       nS3methods <<- nS3methods * 2L
                       new <- matrix(NA_character_, nS3methods, 3L)
                       ind <- seq_len(nold)
                       for (i in 1:3) new[ind, i] <- old[ind, i]
                       S3methods <<- new
                       rm(old, new)
                   }
                   S3methods[nS3, seq_along(spec)] <<- asChar(spec)
               },
               stop(gettextf("unknown namespace directive: %s", deparse(e, nlines=1L)),
                    call. = FALSE, domain = NA)
               )
    }
    for (e in directives)
        parseDirective(e)

    ## need to preserve the names on dynlibs, so unique() is not appropriate.
    dynlibs <- dynlibs[!duplicated(dynlibs)]
    list(imports = imports, exports = exports,
         exportPatterns = unique(exportPatterns),
         importClasses = importClasses, importMethods = importMethods,
         exportClasses = unique(exportClasses),
         exportMethods = unique(exportMethods),
         exportClassPatterns = unique(exportClassPatterns),
         dynlibs = dynlibs, nativeRoutines = nativeRoutines,
         S3methods = unique(S3methods[seq_len(nS3), , drop = FALSE]) )
} ## end{parseNamespaceFile}

registerS3method <- function(genname, class, method, envir = parent.frame()) {
    addNamespaceS3method <- function(ns, generic, class, method) {
        regs <- getNamespaceInfo(ns, "S3methods")
        regs <- rbind(regs, c(generic, class, method))
        setNamespaceInfo(ns, "S3methods", regs)
    }
    groupGenerics <- c("Math", "Ops",  "Summary", "Complex")
    defenv <- if(genname %in% groupGenerics) .BaseNamespaceEnv
    else {
        genfun <- get(genname, envir = envir)
        if(.isMethodsDispatchOn() && methods:::is(genfun, "genericFunction"))
            genfun <- methods:::finalDefaultMethod(genfun@default)
        if (typeof(genfun) == "closure") environment(genfun)
	else .BaseNamespaceEnv
    }
    if (! exists(".__S3MethodsTable__.", envir = defenv, inherits = FALSE))
        assign(".__S3MethodsTable__.", new.env(hash = TRUE, parent = baseenv()),
               envir = defenv)
    table <- get(".__S3MethodsTable__.", envir = defenv, inherits = FALSE)
    if (is.character(method)) {
        assignWrapped <- function(x, method, home, envir) {
            method <- method            # force evaluation
            home <- home                # force evaluation
            delayedAssign(x, get(method, envir = home), assign.env = envir)
        }
        if(!exists(method, envir = envir)) {
            ## need to avoid conflict with message at l.1298
            warning(gettextf("S3 method %s was declared but not found",
                             sQuote(method)), call. = FALSE)
        } else {
	    assignWrapped(paste(genname, class, sep = "."), method, home = envir,
	    	    envir = table)
        }
    }
    else if (is.function(method))
        assign(paste(genname, class, sep = "."), method, envir = table)
    else stop("bad method")
    if (isNamespace(envir) && ! identical(envir, .BaseNamespaceEnv))
        addNamespaceS3method(envir, genname, class, method)
}


registerS3methods <- function(info, package, env)
{
    n <- NROW(info)
    if(n == 0L) return()

    assignWrapped <- function(x, method, home, envir) {
	method <- method            # force evaluation
	home <- home                # force evaluation
	delayedAssign(x, get(method, envir = home), assign.env = envir)
    }
    .registerS3method <- function(genname, class, method, nm, envir)
    {
        ## S3 generics should either be imported explicitly or be in
        ## the base namespace, so we start the search at the imports
        ## environment, parent.env(envir), which is followed by the
        ## base namespace.  (We have already looked in the namespace.)
        ## However, in case they have not been imported, we first
        ## look up where some commonly used generics are (including the
        ## group generics).
        defenv <- if(!is.na(w <- .knownS3Generics[genname])) asNamespace(w)
        else {
            if(!exists(genname, envir = parent.env(envir)))
                stop(gettextf("object '%s' not found whilst loading namespace '%s'",
                              genname, package), call. = FALSE, domain = NA)
            genfun <- get(genname, envir = parent.env(envir))
            if(.isMethodsDispatchOn() && methods:::is(genfun, "genericFunction"))
		genfun <- genfun@default  # nearly always, the S3 generic
            if (typeof(genfun) == "closure") environment(genfun)
            else .BaseNamespaceEnv
        }
        if (! exists(".__S3MethodsTable__.", envir = defenv, inherits = FALSE))
            assign(".__S3MethodsTable__.", new.env(hash = TRUE, parent = baseenv()),
                   envir = defenv)
        table <- get(".__S3MethodsTable__.", envir = defenv, inherits = FALSE)
	assignWrapped(nm, method, home = envir, envir = table)
    }

    methname <- paste(info[,1], info[,2], sep = ".")
    z <- is.na(info[,3])
    info[z,3] <- methname[z]
    Info <- cbind(info, methname)
    loc <- .Internal(ls(env, TRUE))
    notex <- !(info[,3] %in% loc)
    if(any(notex))
        warning(sprintf(ngettext(sum(notex),
                                 "S3 method %s was declared in NAMESPACE but not found",
                                 "S3 methods %s were declared in NAMESPACE but not found"),
                        paste(sQuote(info[notex, 3]), collapse = ", ")),
                call. = FALSE, domain = NA)
    Info <- Info[!notex, , drop = FALSE]

    ## Do local generics first (this could be load-ed if pre-computed).
    ## However, the local generic could be an S4 takeover of a non-local
    ## (or local) S3 generic.  We can't just pass S4 generics on to
    ## .registerS3method as that only looks non-locally (for speed).
    l2 <- localGeneric <- Info[,1] %in% loc
    if(.isMethodsDispatchOn())
        for(i in which(localGeneric)) {
            genfun <- get(Info[i, 1], envir = env)
            if(methods:::is(genfun, "genericFunction")) {
                localGeneric[i] <- FALSE
                registerS3method(Info[i, 1], Info[i, 2], Info[i, 3], env)
            }
        }
    if(any(localGeneric)) {
        lin <- Info[localGeneric, , drop = FALSE]
        S3MethodsTable <-
            get(".__S3MethodsTable__.", envir = env, inherits = FALSE)
        ## we needed to move this to C for speed.
        ## for(i in seq_len(nrow(lin)))
        ##    assign(lin[i,4], get(lin[i,3], envir = env),
        ##           envir = S3MethodsTable)
        .Internal(importIntoEnv(S3MethodsTable, lin[,4], env, lin[,3]))
    }

    ## now the rest
    fin <- Info[!l2, , drop = FALSE]
    for(i in seq_len(nrow(fin)))
        .registerS3method(fin[i, 1], fin[i, 2], fin[i, 3], fin[i, 4], env)

    setNamespaceInfo(env, "S3methods",
                     rbind(info, getNamespaceInfo(env, "S3methods")))
}

.mergeImportMethods <- function(impenv, expenv, metaname)
{
    expMethods <- get(metaname, envir = expenv)
    if(exists(metaname, envir = impenv, inherits = FALSE)) {
        impMethods <- get(metaname, envir = impenv)
        assign(metaname,
               methods:::.mergeMethodsTable2(impMethods,
                                             expMethods, expenv, metaname),
               envir = impenv)
        impMethods
    } else NULL
}
#  File src/library/base/R/notyet.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

.NotYetImplemented <- function ()
    stop(gettextf("'%s' is not implemented yet",
                  as.character(sys.call(sys.parent())[[1L]])), call. = FALSE)

.NotYetUsed <- function(arg, error = TRUE) {
    msg <- gettextf("argument '%s' is not used (yet)", arg)
    if(error) stop(msg, domain = NA, call. = FALSE)
    else warning(msg, domain = NA, call. = FALSE)
}
#  File src/library/base/R/octhex.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

format.octmode <- function(x, width = NULL, ...)
{
    isna <- is.na(x)
    y <- as.integer(x[!isna])
    fmt <- if(!is.null(width)) paste0("%0", width, "o") else "%o"
    ans <- rep.int(NA_character_, length(x))
    ans0 <- sprintf(fmt, y)
    if(is.null(width) && length(y) > 1L) {
        ## previous version padded with zeroes to a common field width
        nc <- max(nchar(ans0))
        ans0 <- sprintf(paste0("%0", nc, "o"), y)
    }
    ans[!isna] <- ans0
    dim(ans) <- dim(x)
    dimnames(ans) <- dimnames(x)
    names(ans) <- names(x)
    ans
}

as.character.octmode <- function(x, ...) format.octmode(x, ...)

print.octmode <- function(x, ...)
{
    print(format(x), ...)
    invisible(x)
}

`[.octmode` <- function (x, i)
{
    cl <- oldClass(x)
    y <- NextMethod("[")
    oldClass(y) <- cl
    y
}

as.octmode <- function(x)
{
    if(inherits(x, "octmode")) return(x)
    if(is.double(x) && x == as.integer(x)) x <- as.integer(x)
    if(is.integer(x)) return(structure(x, class="octmode"))
    if(is.character(x)) {
        z <- strtoi(x, 8L)
        if(!any(is.na(z) | z < 0)) return(structure(z, class="octmode"))
    }
    stop("'x' cannot be coerced to class \"octmode\"")
}

## BioC packages cellHTS2 and flowCore misuse this for doubles,
## hence the as.integer() call
format.hexmode <- function(x, width = NULL, upper.case = FALSE, ...)
{
    isna <- is.na(x)
    y <- as.integer(x[!isna])
    fmt0 <- if(upper.case) "X" else "x"
    fmt <- if(!is.null(width)) paste0("%0", width, fmt0) else paste0("%", fmt0)
    ans <- rep.int(NA_character_, length(x))
    ans0 <- sprintf(fmt, y)
    if(is.null(width) && length(y) > 1L) {
        ## previous version padded with zeroes to a common field width
        nc <- max(nchar(ans0))
        ans0 <- sprintf(paste0("%0", nc, fmt0), y)
    }
    ans[!isna] <- ans0
    dim(ans) <- dim(x)
    dimnames(ans) <- dimnames(x)
    names(ans) <- names(x)
    ans
}

as.character.hexmode <- function(x, ...) format.hexmode(x, ...)

print.hexmode <- function(x, ...)
{
    print(format(x), ...)
    invisible(x)
}

`[.hexmode` <- function (x, i)
{
    cl <- oldClass(x)
    y <- NextMethod("[")
    oldClass(y) <- cl
    y
}

as.hexmode <- function(x)
{
    if(inherits(x, "hexmode")) return(x)
    if(is.double(x) && (x == as.integer(x))) x <- as.integer(x)
    if(is.integer(x)) return(structure(x, class = "hexmode"))
    if(is.character(x)) {
        z <- strtoi(x, 16L)
        if(!any(is.na(z) | z < 0)) return(structure(z, class = "hexmode"))
    }
    stop("'x' cannot be coerced to class \"hexmode\"")
}


`!.octmode` <- function(a) as.octmode(bitwNot(as.octmode(a)))

`&.octmode` <- function(a, b) as.octmode(bitwAnd(as.octmode(a), as.octmode(b)))
`|.octmode` <- function(a, b) as.octmode(bitwOr(as.octmode(a), as.octmode(b)))
xor.octmode <- function(a, b) as.octmode(bitwXor(as.octmode(a), as.octmode(b)))

`!.hexmode` <- function(a) as.hexmode(bitwNot(as.hexmode(a)))

`&.hexmode` <- function(a, b) as.hexmode(bitwAnd(as.hexmode(a), as.hexmode(b)))
`|.hexmode` <- function(a, b) as.hexmode(bitwOr(as.hexmode(a), as.hexmode(b)))
xor.hexmode <- function(a, b) as.hexmode(bitwXor(as.hexmode(a), as.hexmode(b)))
#  File src/library/base/R/options.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

options <- function(...)
    .Internal(options(...))

getOption <- function(x, default = NULL)
{
    ## To avoid always performing the %in%,
    ## we use the original code if default is not specified.
    if(missing(default)) return(options(x)[[1L]])

    if(x %in% names(options())) options(x)[[1L]] else default
}
#  File src/library/base/R/outer.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2013 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

outer <- function (X, Y, FUN = "*", ...)
{
    if(is.array(X)) {
        dX <- dim(X)
        nx <- dimnames(X)
        no.nx <- is.null(nx)
    } else { # a vector
        dX <- length(X)  # cannot be long, as form a matrix below
        no.nx <- is.null(names(X))
        if(!no.nx) nx <- list(names(X))
    }
    if(is.array(Y)) {
        dY <- dim(Y)
        ny <- dimnames(Y)
        no.ny <- is.null(ny)
    } else { # a vector
        dY <- length(Y)
        no.ny <- is.null(names(Y))
        if(!no.ny) ny <- list(names(Y))
    }
    if (is.character(FUN) && FUN=="*") {
        if(!missing(...)) stop('using ... with FUN = "*" is an error')
        # this is for numeric vectors, so dropping attributes is OK
        robj <- as.vector(X) %*% t(as.vector(Y))
        dim(robj) <- c(dX, dY)
    } else {
        FUN <- match.fun(FUN)
        ## Y may have a class, so don't use rep.int
        Y <- rep(Y, rep.int(length(X), length(Y)))
        ##  length.out is not an argument of the generic rep()
        ##  X <- rep(X, length.out = length(Y))
        if(length(X))
            X <- rep(X, times = ceiling(length(Y)/length(X)))
        robj <- FUN(X, Y, ...)
        dim(robj) <- c(dX, dY) # careful not to lose class here
    }
    ## no dimnames if both don't have ..
    if(!(no.nx && no.ny)) {
	if(no.nx) nx <- vector("list", length(dX)) else
	if(no.ny) ny <- vector("list", length(dY))
	dimnames(robj) <- c(nx, ny)
    }
    robj
}

## Binary operator, hence don't simply do "%o%" <- outer.
`%o%` <- function(X, Y) outer(X, Y)
#  File src/library/base/R/pairlist.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

as.pairlist <- function(x) .Internal(as.vector(x, "pairlist"))
pairlist <- function(...) as.pairlist(list(...))
## This is now .Primitive:
##is.pairlist <- function(x) typeof(x) == "pairlist"
#  File src/library/base/R/parse.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

parse <- function(file = "", n = NULL, text = NULL, prompt = "?",
		  keep.source = getOption("keep.source"),
                  srcfile = NULL, encoding = "unknown")
{
    keep.source <- isTRUE(keep.source)
    if(!is.null(text)) {
    	if (length(text) == 0L) return(expression())
	if (missing(srcfile)) {
	    srcfile <- "<text>"
	    if (keep.source)
	       srcfile <- srcfilecopy(srcfile, text)
	}
	file <- stdin()
    } else {
	if(is.character(file)) {
            if(file == "") {
            	file <- stdin()
            	if (missing(srcfile))
            	    srcfile <- "<stdin>"
            } else {
		filename <- file
		file <- file(filename, "r")
            	if (missing(srcfile))
            	    srcfile <- filename
            	if (keep.source) {
		    text <- readLines(file, warn = FALSE)
		    if (!length(text)) text <- ""
            	    close(file)
            	    file <- stdin()
        	    srcfile <- srcfilecopy(filename, text, file.info(filename)[1,"mtime"],
        	                       isFile = TRUE)
                } else
		    on.exit(close(file))
	    }
	}
    }
    .Internal(parse(file, n, text, prompt, srcfile, encoding))
}
#  File src/library/base/R/paste.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

paste <- function (..., sep = " ", collapse = NULL)
    .Internal(paste(list(...), sep, collapse))
paste0 <- function(..., collapse = NULL)
    .Internal(paste0(list(...), collapse))

##=== Could we extend  paste(.) to (optionally) accept a
##    2-vector for collapse ?	 With the following functionality

##- paste.extra <- function(r, collapse=c(", "," and ")) {
##-	    n <- length(r)
##-	    if(n <= 1) paste(r)
##-	    else
##-	      paste(paste(r[-n],collapse=collapse[1L]),
##-		    r[n], sep=collapse[min(2,length(collapse))])
##- }
#  File src/library/base/R/pmax.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

pmin.int <- function(..., na.rm = FALSE) .Internal(pmin(na.rm, ...))
pmax.int <- function(..., na.rm = FALSE) .Internal(pmax(na.rm, ...))

pmax <- function (..., na.rm = FALSE)
{
    elts <- list(...)
    if(length(elts) == 0L) stop("no arguments")
    if(all(vapply(elts, function(x) is.atomic(x) && !is.object(x), NA))) {
        ## NB: NULL passes is.atomic
        mmm <- .Internal(pmax(na.rm, ...))
    } else {
        mmm <- elts[[1L]]
        attr(mmm, "dim") <- NULL  # dim<- would drop names
        has.na <- FALSE
        for (each in elts[-1L]) {
            attr(each, "dim") <- NULL
            l1 <- length(each); l2 <- length(mmm)
            if(l2 < l1) {
                if (l2 && l1 %% l2)
                    warning("an argument will be fractionally recycled")
                mmm <- rep(mmm, length.out = l1)
            } else if(l1 && l1 < l2) {
                if (l2 %% l1)
                    warning("an argument will be fractionally recycled")
                each <- rep(each, length.out = l2)
            }
            nas <- cbind(is.na(mmm), is.na(each))
            if(has.na || (has.na <- any(nas))) {
                mmm[nas[, 1L]] <- each[nas[, 1L]]
                each[nas[, 2L]] <- mmm[nas[, 2L]]
            }
            change <- mmm < each
            change <- change & !is.na(change)
            mmm[change] <- each[change]
            if (has.na && !na.rm) mmm[nas[, 1L] | nas[, 2L]] <- NA
        }
    }
    mostattributes(mmm) <- attributes(elts[[1L]])
    mmm
}

pmin <- function (..., na.rm = FALSE)
{
    elts <- list(...)
    if(length(elts) == 0L) stop("no arguments")
    if(all(vapply(elts, function(x) is.atomic(x) && !is.object(x), NA))) {
        mmm <- .Internal(pmin(na.rm, ...))
    } else {
        mmm <- elts[[1L]]
        attr(mmm, "dim") <- NULL  # dim<- would drop names
        has.na <- FALSE
        for (each in elts[-1L]) {
            attr(each, "dim") <- NULL
            l1 <- length(each); l2 <- length(mmm)
            if(l2 < l1) {
                if (l2 && l1 %% l2)
                    warning("an argument will be fractionally recycled")
                mmm <- rep(mmm, length.out = l1)
            } else if(l1 && l1 < l2) {
                if (l2 %% l1)
                    warning("an argument will be fractionally recycled")
                each <- rep(each, length.out = l2)
            }
            nas <- cbind(is.na(mmm), is.na(each))
            if(has.na || (has.na <- any(nas))) {
                mmm[nas[, 1L]] <- each[nas[, 1L]]
                each[nas[, 2L]] <- mmm[nas[, 2L]]
            }
            change <- mmm > each
            change <- change & !is.na(change)
            mmm[change] <- each[change]
            if (has.na && !na.rm) mmm[nas[, 1L] | nas[, 2L]] <- NA
        }
    }
    mostattributes(mmm) <- attributes(elts[[1L]])
    mmm
}
#  File src/library/base/R/pretty.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

pretty <- function(x, ...) UseMethod("pretty")


pretty.default <-
    function(x, n = 5, min.n = n %/% 3, shrink.sml = 0.75,
             high.u.bias = 1.5, u5.bias = .5 + 1.5*high.u.bias,
             eps.correct = 0, ...)
{
    x <- x[is.finite(x <- as.numeric(x))]
    if(!length(x)) return(x)
    z <- .Internal(pretty(min(x), max(x), n, min.n, shrink.sml,
                          c(high.u.bias, u5.bias), eps.correct))
    s <- seq.int(z$l, z$u, length.out = z$n + 1)
    if(!eps.correct && z$n) { # maybe zap smalls from seq() rounding errors
        ## better than zapsmall(s, digits = 14) :
        delta <- diff(range(z$l, z$u)) / z$n  # or abs(z$u - z$l)
        if(any(small <- abs(s) < 1e-14 * delta)) s[small] <- 0
    }
    s
}
#  File src/library/base/R/print.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

print <- function(x, ...) UseMethod("print")

##- Need '...' such that it can be called as  NextMethod("print", ...):
print.default <- function(x, digits = NULL, quote = TRUE, na.print = NULL,
                          print.gap = NULL, right = FALSE, max = NULL,
                          useSource = TRUE, ...)
{
    noOpt <- missing(digits) && missing(quote) && missing(na.print) &&
	missing(print.gap) && missing(right) && missing(max) &&
	missing(useSource) && missing(...)
    .Internal(print.default(x, digits, quote, na.print, print.gap, right, max,
			    useSource, noOpt))
}

prmatrix <-
    function (x, rowlab = dn[[1]], collab = dn[[2]],
              quote = TRUE, right = FALSE,
              na.print = NULL, ...)
{
    x <- as.matrix(x)
    dn <- dimnames(x)
    .Internal(prmatrix(x, rowlab, collab, quote, right, na.print))
}

noquote <- function(obj) {
    ## constructor for a useful "minor" class
    if(!inherits(obj,"noquote")) class(obj) <- c(attr(obj, "class"),"noquote")
    obj
}

as.matrix.noquote <- function(x, ...) noquote(NextMethod("as.matrix", x))
c.noquote <- function(..., recursive = FALSE)
    structure(NextMethod("c"), class = "noquote")

`[.noquote` <- function (x, ...) {
    attr <- attributes(x)
    r <- unclass(x)[...] ## shouldn't this be NextMethod?
    attributes(r) <- c(attributes(r),
		       attr[is.na(match(names(attr),
                                        c("dim","dimnames","names")))])
    r
}

print.noquote <- function(x, ...) {
    if(!is.null(cl <- attr(x, "class"))) {
	cl <- cl[cl != "noquote"]
        attr(x, "class") <-
          (if(length(cl)) cl else NULL)
      }
    print(x, quote = FALSE, ...)
}

## for alias:
print.listof <- function(x, ...)
{
    nn <- names(x)
    ll <- length(x)
    if(length(nn) != ll) nn <- paste("Component", seq.int(ll))
    for(i in seq_len(ll)) {
	cat(nn[i], ":\n"); print(x[[i]], ...); cat("\n")
    }
    invisible(x)
}

## formerly same as [.AsIs
`[.listof` <- function(x, i, ...) structure(NextMethod("["), class = class(x))

## used for version:
print.simple.list <- function(x, ...)
    print(noquote(cbind("_"=unlist(x))), ...)

`[.simple.list` <- `[.listof`

print.function <- function(x, useSource = TRUE, ...)
    .Internal(print.function(x, useSource, ...))
#  File src/library/base/R/qr.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

is.qr <- function(x) inherits(x, "qr")

qr <- function(x, ...) UseMethod("qr")

qr.default <- function(x, tol = 1e-07, LAPACK = FALSE, ...)
{
    x <- as.matrix(x)
    if(is.complex(x))
        return(structure(.Internal(La_qr_cmplx(x)), class = "qr"))
    ## otherwise :
    if(LAPACK)
        return(structure(.Internal(La_qr(x)), useLAPACK = TRUE, class = "qr"))

    p <- as.integer(ncol(x))
    if(is.na(p)) stop("invalid ncol(x)")
    n <- as.integer(nrow(x))
    if(is.na(n)) stop("invalid nrow(x)")
    if(1.0 * n * p > 2147483647) stop("too large a matrix for LINPACK")
    storage.mode(x) <- "double"
    res <- .Fortran(.F_dqrdc2,
	     qr = x,
	     n,
	     n,
	     p,
	     as.double(tol),
	     rank = integer(1L),
	     qraux = double(p),
	     pivot = as.integer(seq_len(p)),
	     double(2L*p))[c(1,6,7,8)]# c("qr", "rank", "qraux", "pivot")
    if(!is.null(cn <- colnames(x)))
        colnames(res$qr) <- cn[res$pivot]
    class(res) <- "qr"
    res
}

## + qr.lm  method defined in ../../stats/R/lm.R


qr.coef <- function(qr, y)
{
    if( !is.qr(qr) ) stop("first argument must be a QR decomposition")
    n <- as.integer(nrow(qr$qr)); if(is.na(n)) stop("invalid nrow(qr$qr)")
    p <- as.integer(ncol(qr$qr)); if(is.na(p)) stop("invalid ncol(qr$qr)")
    k <- as.integer(qr$rank); if(is.na(k)) stop("invalid ncol(qr$rank)")
    im <- is.matrix(y)
    if (!im) y <- as.matrix(y)
    ny <- as.integer(ncol(y))
    if(is.na(ny)) stop("invalid ncol(y)")
    if (p == 0L) return( if (im) matrix(0, p, ny) else numeric() )
    ix <- if ( p > n ) c(seq_len(n), rep(NA, p - n)) else seq_len(p)
    if(is.complex(qr$qr)) {
	coef <- matrix(NA_complex_, nrow = p, ncol = ny)
	coef[qr$pivot, ] <- .Internal(qr_coef_cmplx(qr, y))[ix, ]
	return(if(im) coef else c(coef))
    }
    ## else {not complex} :
    if(isTRUE(attr(qr, "useLAPACK"))) {
	coef <- matrix(NA_real_, nrow = p, ncol = ny)
	coef[qr$pivot, ] <- .Internal(qr_coef_real(qr, y))[ix, ]
	return(if(im) coef else c(coef))
    }
    if (k == 0L) return( if (im) matrix(NA, p, ny) else rep.int(NA, p))

    storage.mode(y) <- "double"
    if( nrow(y) != n )
	stop("'qr' and 'y' must have the same number of rows")
    z <- .Fortran(.F_dqrcf,
		  as.double(qr$qr),
		  n, k,
		  as.double(qr$qraux),
		  y,
		  ny,
		  coef = matrix(0, nrow = k,ncol = ny),
		  info = integer(1L),
		  NAOK = TRUE)[c("coef","info")]
    if(z$info) stop("exact singularity in 'qr.coef'")
    if(k < p) {
	coef <- matrix(NA_real_, nrow = p, ncol = ny)
	coef[qr$pivot[seq_len(k)], ] <- z$coef
    }
    else coef <- z$coef

    if(!is.null(nam <- colnames(qr$qr)))
        if(k < p) rownames(coef)[qr$pivot] <- nam
        else rownames(coef) <- nam

    if(im && !is.null(nam <- colnames(y)))
        colnames(coef) <- nam

    if(im) coef else drop(coef)
}

qr.qy <- function(qr, y)
{
    if(!is.qr(qr)) stop("argument is not a QR decomposition")
    if(is.complex(qr$qr))
        return(.Internal(qr_qy_cmplx(qr, as.matrix(y), FALSE)))
    if(isTRUE(attr(qr, "useLAPACK")))
        return(.Internal(qr_qy_real(qr, as.matrix(y), FALSE)))

    n <- as.integer(nrow(qr$qr))
    if(is.na(n)) stop("invalid nrow(qr$qr)")
    k <- as.integer(qr$rank)
    ny <- as.integer(NCOL(y))
    if(is.na(ny)) stop("invalid NCOL(y)")
    storage.mode(y) <- "double"
    if(NROW(y) != n)
	stop("'qr' and 'y' must have the same number of rows")
    .Fortran(.F_dqrqy,
	     as.double(qr$qr),
	     n, k,
	     as.double(qr$qraux),
	     y,
	     ny,
	     qy = y# incl. {dim}names
	     )$qy
}

qr.qty <- function(qr, y)
{
    if(!is.qr(qr)) stop("argument is not a QR decomposition")
    if(is.complex(qr$qr))
        return(.Internal(qr_qy_cmplx(qr, as.matrix(y), TRUE)))
    if(isTRUE(attr(qr, "useLAPACK")))
        return(.Internal(qr_qy_real(qr, as.matrix(y), TRUE)))

    n <- as.integer(nrow(qr$qr))
    if(is.na(n)) stop("invalid nrow(qr$qr)")
    k <- as.integer(qr$rank)
    ny <- as.integer(NCOL(y))
    if(is.na(ny)) stop("invalid NCOL(y)")
    if(NROW(y) != n)
	stop("'qr' and 'y' must have the same number of rows")
    storage.mode(y) <- "double"
    .Fortran(.F_dqrqty,
	     as.double(qr$qr),
	     n, k,
	     as.double(qr$qraux),
	     y,
	     ny,
	     qty = y# incl. {dim}names
             )$qty
}

qr.resid <- function(qr, y)
{
    if(!is.qr(qr)) stop("argument is not a QR decomposition")
    if(is.complex(qr$qr)) stop("not implemented for complex 'qr'")
    if(isTRUE(attr(qr, "useLAPACK"))) stop("not supported for LAPACK QR")
    k <- as.integer(qr$rank)
    if (k==0) return(y)
    n <- as.integer(nrow(qr$qr))
    if(is.na(n)) stop("invalid nrow(qr$qr)")
    ny <- as.integer(NCOL(y))
    if(is.na(ny)) stop("invalid NCOL(y)")
    if( NROW(y) != n )
	stop("'qr' and 'y' must have the same number of rows")
    storage.mode(y) <- "double"
    .Fortran(.F_dqrrsd,
	     as.double(qr$qr), n, k, as.double(qr$qraux), y, ny, rsd = y)$rsd
}

qr.fitted <- function(qr, y, k=qr$rank)
{
    if(!is.qr(qr)) stop("argument is not a QR decomposition")
    if(is.complex(qr$qr)) stop("not implemented for complex 'qr'")
    if(isTRUE(attr(qr, "useLAPACK"))) stop("not supported for LAPACK QR")
    n <- as.integer(nrow(qr$qr))
    if(is.na(n)) stop("invalid nrow(qr$qr)")
    k <- as.integer(k)
    if(k > qr$rank) stop("'k' is too large")
    ny <- as.integer(NCOL(y))
    if(is.na(ny)) stop("invalid NCOL(y)")
    if( NROW(y) != n )
	stop("'qr' and 'y' must have the same number of rows")
    storage.mode(y) <- "double"
    .Fortran(.F_dqrxb,
	     as.double(qr$qr), n, k, as.double(qr$qraux), y, ny, xb = y)$xb
}

## qr.solve is defined in  ./solve.R

##---- The next three are from Doug Bates ('st849'):
qr.Q <- function (qr, complete = FALSE, Dvec)
{
    if(!is.qr(qr)) stop("argument is not a QR decomposition")
    dqr <- dim(qr$qr)
    n <- dqr[1L]
    cmplx <- mode(qr$qr) == "complex"
    if(missing(Dvec))
	Dvec <- rep.int(if (cmplx) 1 + 0i else 1,
			if (complete) n else min(dqr))
    D <-
	if (complete) diag(Dvec, n)
	else {
	    ncols <- min(dqr)
	    diag(Dvec[seq_len(ncols)], nrow = n, ncol = ncols)
	}
    qr.qy(qr, D)
}

qr.R <- function (qr, complete = FALSE, ...)
{
    if(!is.qr(qr)) stop("argument is not a QR decomposition")
    if(!missing(...))
	warning(sprintf(ngettext(length(list(...)),
				 "extra argument %s will be disregarded",
				 "extra arguments %s will be disregarded"),
			 paste(sQuote(names(list(...))), collapse = ", ")),
		domain = NA)
    R <- qr$qr
    if (!complete)
	R <- R[seq.int(min(dim(R))), , drop = FALSE]
    R[row(R) > col(R)] <- 0
    R
}

qr.X <- function (qr, complete = FALSE,
		  ncol = if (complete) nrow(R) else min(dim(R)))
{
    if(!is.qr(qr)) stop("argument is not a QR decomposition")
    pivoted <- !identical(qr$pivot, ip <- seq_along(qr$pivot))
    R <- qr.R(qr, complete = TRUE)
    if(pivoted && ncol < length(qr$pivot))
        stop("need larger value of 'ncol' as pivoting occurred")
    cmplx <- mode(R) == "complex"
    p <- as.integer(dim(R)[2L])
    if(is.na(p)) stop("invalid NCOL(R)")
    if (ncol < p)
	R <- R[, seq_len(ncol), drop = FALSE]
    else if (ncol > p) {
	tmp <- diag(if (!cmplx) 1 else 1 + 0i, nrow(R), ncol)
	tmp[, seq_len(p)] <- R
	R <- tmp
    }
    res <- qr.qy(qr, R)
    cn <- colnames(res)
    if(pivoted) {# res may have more columns than length(qr$pivot)
	res[, qr$pivot] <- res[, ip]
        if(!is.null(cn)) colnames(res)[qr$pivot] <- cn[ip]
    }
    res
}
#  File src/library/base/R/quit.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

quit <- function(save = "default", status=0, runLast=TRUE)
    .Internal(quit(save, status, runLast))
q <- quit
#  File src/library/base/R/range.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

range.default <- function(..., na.rm = FALSE, finite = FALSE)
{
    x <- c(..., recursive = TRUE)
    if(is.numeric(x)) {
        if(finite) x <- x[is.finite(x)]
        else if(na.rm) x <- x[!is.na(x)]
	return (c(min(x), max(x)))
    }
    c(min(x, na.rm=na.rm), max(x, na.rm=na.rm))
}
#  File src/library/base/R/rank.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

rank <- function(x, na.last = TRUE,
		 ties.method = c("average", "first", "random", "max", "min"))
{
    nas <- is.na(x)
    nm <- names(x)
    ties.method <- match.arg(ties.method)
    ## To preserve past behaviour
    if(is.factor(x)) x <- as.integer(x)
    xx <- x[!nas]
    ## we pass length(xx) to allow
    y <- switch(ties.method,
		"average" = , "min" = , "max" =
		.Internal(rank(xx, length(xx), ties.method)),
		"first" = sort.list(sort.list(xx)),
		"random" = sort.list(order(xx, stats::runif(sum(!nas)))))
    if(!is.na(na.last) && any(nas)) {
	## the internal code has ranks in [1, length(y)]
	yy <- integer(length(x))
	storage.mode(yy) <- storage.mode(y) # integer or double
	yy <- NA
	NAkeep <- (na.last == "keep")
	if(NAkeep || na.last) {
	    yy[!nas] <- y
	    if(!NAkeep) yy[nas] <- (length(y) + 1L) : length(yy)
	} else {
	    len <- sum(nas)
	    yy[!nas] <- y + len
	    yy[nas] <- 1L : len
	}
	y <- yy
	names(y) <- nm
    } else names(y) <- nm[!nas]
    y
}
#  File src/library/base/R/raw.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

raw <- function(length = 0L) .Internal(vector("raw", length))

#as.raw <- function(x) .Internal(as.raw(x))

charToRaw <- function(x) .Internal(charToRaw(x))
rawToChar <- function(x, multiple=FALSE) .Internal(rawToChar(x, multiple))

rawShift <- function(x, n) .Internal(rawShift(x, n))

rawToBits <- function(x) .Internal(rawToBits(x))
intToBits <- function(x) .Internal(intToBits(x))

packBits <- function(x, type=c("raw", "integer"))
{
    type <- match.arg(type)
    .Internal(packBits(x, type))
}

utf8ToInt <- function(x) .Internal(utf8ToInt(x))
intToUtf8 <- function(x, multiple=FALSE) .Internal(intToUtf8(x, multiple))
#  File src/library/base/R/rep.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

rep.int <- function(x, times) .Internal(rep.int(x, times))

rep_len <- function(x, length.out) .Internal(rep_len(x, length.out))


rep.factor <- function(x, ...)
{
    y <- NextMethod()
    structure(y, class=class(x), levels=levels(x))
}
#  File src/library/base/R/replace.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

replace <-
    function (x, list, values)
{
    x[list] <- values
    x
}
#  File src/library/base/R/replicate.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

replicate <- function(n, expr, simplify = "array")
        sapply(integer(n),
           eval.parent(substitute(function(...)expr)), simplify = simplify)
#  File src/library/base/R/rev.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

rev <- function(x) UseMethod("rev")

rev.default <- function(x) if (length(x)) x[length(x):1L] else x
#  File src/library/base/R/rle.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

rle <- function(x)
{
    if (!is.vector(x) && !is.list(x))
        stop("'x' must be an atomic vector")
    n <- length(x)
    if (n == 0L)
	return(structure(list(lengths = integer(), values = x),
			 class = "rle"))
    y <- x[-1L] != x[-n]
    i <- c(which(y | is.na(y)), n)
    structure(list(lengths = diff(c(0L, i)), values = x[i]),
              class = "rle")
}

print.rle <- function(x, digits = getOption("digits"), prefix = "", ...)
{
    if(is.null(digits)) digits <- getOption("digits")
    cat("", "Run Length Encoding\n", "  lengths:", sep=prefix)
    utils::str(x$lengths)
    cat("", "  values :", sep=prefix)
    utils::str(x$values, digits.d = digits)
    invisible(x)
}

inverse.rle <- function(x, ...)
{
    if(is.null(le <- x$lengths) ||
       is.null(v  <- x$values) || length(le) != length(v))
        stop("invalid 'rle' structure")
    rep.int(v, le)
}
#  File src/library/base/R/rm.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

rm <-
    function (..., list = character(), pos = -1, envir = as.environment(pos),
              inherits = FALSE)
{
    dots <- match.call(expand.dots=FALSE)$...
    if(length(dots) &&
       !all(sapply(dots, function(x) is.symbol(x) || is.character(x))))
       stop("... must contain names or character strings")
    names <- sapply(dots, as.character)
    if (length(names) == 0L) names <- character()
    list <- .Primitive("c")(list, names)
    .Internal(remove(list, envir, inherits))
}

remove <- rm
#  File src/library/base/R/rowsum.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

rowsum <- function(x, group, reorder = TRUE, ...) UseMethod("rowsum")

rowsum.default <- function(x, group, reorder = TRUE, na.rm = FALSE, ...)
{
    if (!is.numeric(x)) stop("'x' must be numeric")
    if (length(group) != NROW(x)) stop("incorrect length for 'group'")
    if (any(is.na(group))) warning("missing values for 'group'")
    ugroup <- unique(group)
    if (reorder) ugroup <- sort(ugroup, na.last = TRUE, method = "quick")
    ## ugroup can be either a vector or a factor, so do as.character here
    .Internal(rowsum_matrix(x, group, ugroup, na.rm, as.character(ugroup)))
}

rowsum.data.frame <- function(x, group, reorder = TRUE, na.rm = FALSE, ...)
{
    if (!is.data.frame(x)) stop("not a data frame") ## make MM happy
    if (length(group) != NROW(x)) stop("incorrect length for 'group'")
    if (any(is.na(group))) warning("missing values for 'group'")
    ugroup <- unique(group)
    if (reorder) ugroup <- sort(ugroup, na.last = TRUE, method = "quick")
    .Internal(rowsum_df(x, group, ugroup, na.rm, as.character(ugroup)))
}
#  File src/library/base/R/sample.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

sample <- function(x, size, replace = FALSE, prob = NULL)
{
    if(length(x) == 1L && is.numeric(x) && x >= 1) {
	if(missing(size)) size <- x
	sample.int(x, size, replace, prob)
    } else {
	if(missing(size)) size <- length(x)
	x[sample.int(length(x), size, replace, prob)]
    }
}

sample.int  <- function(n, size = n, replace = FALSE, prob = NULL)
{
    if (!replace && is.null(prob) && n > 1e7 && size <= n/2)
        .Internal(sample2(n, size))
    else .Internal(sample(n, size, replace, prob))
}
#  File src/library/base/R/sapply.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

##' "Simplify" a list of commonly structured components into an array.
##'
##' @title simplify list() to an array if the list elements are structurally equal
##' @param x a list, typically resulting from lapply()
##' @param higher logical indicating if an array() of "higher rank"
##'  should be returned when appropriate, namely when all elements of
##' \code{x} have the same \code{\link{dim}()}ension.
##' @return x itself, or an array if the simplification "is sensible"
simplify2array <- function(x, higher = TRUE)
{
    if(length(common.len <- unique(unlist(lapply(x, length)))) > 1L)
        return(x)
    if(common.len == 1L)
        unlist(x, recursive = FALSE)
    else if(common.len > 1L) {
        n <- length(x)
        ## make sure that array(*) will not call rep() {e.g. for 'call's}:
        r <- as.vector(unlist(x, recursive = FALSE))
        if(higher && length(c.dim <- unique(lapply(x, dim))) == 1 &&
           is.numeric(c.dim <- c.dim[[1L]]) &&
           prod(d <- c(c.dim, n)) == length(r)) {

            iN1 <- is.null(n1 <- dimnames(x[[1L]]))
            n2 <- names(x)
            dnam <-
                if(!(iN1 && is.null(n2)))
                    c(if(iN1) rep.int(list(n1), length(c.dim)) else n1,
                      list(n2)) ## else NULL
            array(r, dim = d, dimnames = dnam)

        } else if(prod(d <- c(common.len, n)) == length(r))
            array(r, dim = d,
                  dimnames = if(!(is.null(n1 <- names(x[[1L]])) &
                  is.null(n2 <- names(x)))) list(n1,n2))
        else x
    }
    else x
}

sapply <- function(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE)
{
    FUN <- match.fun(FUN)
    answer <- lapply(X = X, FUN = FUN, ...)
    if(USE.NAMES && is.character(X) && is.null(names(answer)))
	names(answer) <- X
    if(!identical(simplify, FALSE) && length(answer))
	simplify2array(answer, higher = (simplify == "array"))
    else answer
}

vapply <- function(X, FUN, FUN.VALUE, ...,  USE.NAMES = TRUE)
{
    FUN <- match.fun(FUN)
    if(!is.vector(X) || is.object(X)) X <- as.list(X)
    .Internal(vapply(X, FUN, FUN.VALUE, USE.NAMES))
}


#  File src/library/base/R/scale.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

scale <- function(x, center = TRUE, scale = TRUE) UseMethod("scale")

scale.default <- function(x, center = TRUE, scale = TRUE)
{
    x <- as.matrix(x)
    nc <- ncol(x)
    if (is.logical(center)) {
	if (center) {
            center <- colMeans(x, na.rm=TRUE)
	    x <- sweep(x, 2L, center, check.margin=FALSE)
        }
    }
    else if (is.numeric(center) && (length(center) == nc))
	x <- sweep(x, 2L, center, check.margin=FALSE)
    else
	stop("length of 'center' must equal the number of columns of 'x'")
    if (is.logical(scale)) {
	if (scale) {
	    f <- function(v) {
		v <- v[!is.na(v)]
		sqrt(sum(v^2) / max(1, length(v) - 1L))
	    }
            scale <- apply(x, 2L, f)
	    x <- sweep(x, 2L, scale, "/", check.margin=FALSE)
	}
    }
    else if (is.numeric(scale) && length(scale) == nc)
	x <- sweep(x, 2L, scale, "/", check.margin=FALSE)
    else
	stop("length of 'scale' must equal the number of columns of 'x'")
    if(is.numeric(center)) attr(x, "scaled:center") <- center
    if(is.numeric(scale)) attr(x, "scaled:scale") <- scale
    x
}
#  File src/library/base/R/scan.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

scan <-
function(file = "", what = double(), nmax = -1L, n = -1L, sep = "",
         quote = if(identical(sep, "\n")) "" else "'\"",
         dec = ".", skip = 0L, nlines = 0L,
         na.strings = "NA", flush = FALSE, fill = FALSE,
         strip.white = FALSE, quiet = FALSE, blank.lines.skip = TRUE,
         multi.line = TRUE, comment.char = "", allowEscapes = FALSE,
         fileEncoding = "", encoding = "unknown", text)
{
    na.strings <- as.character(na.strings)# allow it to be NULL
    if(!missing(n)) {
        if(missing(nmax))
            nmax <- n / pmax(length(what), 1L)
        else
            stop("either specify 'nmax' or 'n', but not both.")
    }
    if (missing(file) && !missing(text)) {
	file <- textConnection(text)
	on.exit(close(file))
    }

    if(is.character(file))
        if(file == "") file <- stdin()
        else {
            file <- if(nzchar(fileEncoding))
                file(file, "r", encoding = fileEncoding) else file(file, "r")
            on.exit(close(file))
        }
    if(!inherits(file, "connection"))
        stop("'file' must be a character string or connection")
    .Internal(scan(file, what, nmax, sep, dec, quote, skip, nlines,
                   na.strings, flush, fill, strip.white, quiet,
                   blank.lines.skip, multi.line, comment.char,
                   allowEscapes, encoding))
}
#  File src/library/base/R/seq.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

seq <- function(...) UseMethod("seq")

seq.default <-
    function(from = 1, to = 1, by = ((to - from)/(length.out - 1)),
             length.out = NULL, along.with = NULL, ...)
{
    if((One <- nargs() == 1L) && !missing(from)) {
	lf <- length(from)
	return(if(mode(from) == "numeric" && lf == 1L) {
            if(!is.finite(from)) stop("'from' cannot be NA, NaN or infinite")
            1L:from
        } else if(lf) 1L:lf else integer())
    }
    if(!missing(along.with)) {
	length.out <- length(along.with)
	if(One) return(if(length.out) seq_len(length.out) else integer())
    }
    else if(!missing(length.out)) {
        len <- length(length.out)
        if(!len) stop("argument 'length.out' must be of length 1")
        if(len > 1L) {
            warning("first element used of 'length.out' argument")
            length.out <- length.out[1L]
        }
	length.out <- ceiling(length.out)
    }
    if(!missing(...))
        warning(sprintf(ngettext(length(list(...)),
                                 "extra argument %s will be disregarded",
                                 "extra arguments %s will be disregarded"),
			 paste(sQuote(names(list(...))), collapse = ", ")),
		domain = NA)
    if (!missing(from) && length(from) != 1L) stop("'from' must be of length 1")
    if (!missing(to) && length(to) != 1L) stop("'to' must be of length 1")
    if (!missing(from) && !is.finite(from))
        stop("'from' cannot be NA, NaN or infinite")
    if (!missing(to) && !is.finite(to))
        stop("'to' cannot be NA, NaN or infinite")
    if(is.null(length.out))
	if(missing(by))
	    from:to
	else { # dealing with 'by'
	    del <- to - from
	    if(del == 0 && to == 0) return(to)
	    n <- del/by
	    if(!(length(n) && is.finite(n))) {
		if(length(by) && by == 0 && length(del) && del == 0)
		    return(from)
		stop("invalid (to - from)/by in seq(.)")
	    }
	    if(n < 0L)
		stop("wrong sign in 'by' argument")
	    if(n > .Machine$integer.max)
		stop("'by' argument is much too small")

	    dd <- abs(del)/max(abs(to), abs(from))
	    if (dd < 100*.Machine$double.eps) return(from)
            if (is.integer(del) && is.integer(by)) {
                n <- as.integer(n) # truncates
                from + (0L:n) * by
            } else {
                n <- as.integer(n + 1e-10)
                x <- from + (0L:n) * by
                ## correct for possible overshot because of fuzz
                if(by > 0) pmin(x, to) else pmax(x, to)
            }
	}
    else if(!is.finite(length.out) || length.out < 0L)
	stop("length must be non-negative number")
    else if(length.out == 0L) integer()
    else if (One) seq_len(length.out)
    else if(missing(by)) {
	# if(from == to || length.out < 2) by <- 1
	if(missing(to))
	    to <- from + length.out - 1L
	if(missing(from))
	    from <- to - length.out + 1L
	if(length.out > 2L) # not clear why these have as.vector, and not others
	    if(from == to) rep.int(from, length.out)
	    else as.vector(c(from, from + seq_len(length.out - 2L) * by, to))
	else as.vector(c(from, to))[seq_len(length.out)]
    }
    else if(missing(to))
	from + (0L:(length.out - 1L)) * by
    else if(missing(from))
	to - ((length.out - 1L):0L) * by
    else stop("too many arguments")
}

## In reverence to the very first versions of R which already had sequence():
sequence <- function(nvec) unlist(lapply(nvec, seq_len))
#  File src/library/base/R/serialize.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

saveRDS <-
    function(object, file = "", ascii = FALSE, version = NULL,
             compress = TRUE, refhook = NULL)
{
    if(is.character(file)) {
        if(file == "") stop("'file' must be non-empty string")
        mode <- if(ascii) "w" else "wb"
        con <- if (identical(compress, "bzip2")) bzfile(file, mode)
            else if (identical(compress, "xz")) xzfile(file, mode)
            else if(compress) gzfile(file, mode) else file(file, mode)
        on.exit(close(con))
    }
    else if(inherits(file, "connection")) {
        if (!missing(compress))
            warning("'compress' is ignored unless 'file' is a file name")
        con <- file
    }
    else
        stop("bad 'file' argument")
    .Internal(serializeToConn(object, con, ascii, version, refhook))
}

readRDS <- function(file, refhook = NULL)
{
    if(is.character(file)) {
        con <- gzfile(file, "rb")
        on.exit(close(con))
    } else if(inherits(file, "connection"))
        con <- file
    else stop("bad 'file' argument")
    .Internal(unserializeFromConn(con, refhook))
}

serialize <-
    function(object, connection, ascii = FALSE, xdr = TRUE,
             version = NULL, refhook = NULL)
{
    if (!is.null(connection)) {
        if (!inherits(connection, "connection"))
            stop("'connection' must be a connection")
        if (missing(ascii)) ascii <- summary(connection)$text == "text"
    }
    if (!ascii && inherits(connection, "sockconn"))
        .Internal(serializeb(object, connection, xdr, version, refhook))
    else {
        if (!isTRUE(ascii) && !xdr) ascii <- NA
        .Internal(serialize(object, connection, ascii, version, refhook))
    }
}

unserialize <- function(connection, refhook = NULL)
{
    if (typeof(connection) != "raw" &&
        !is.character(connection) &&
        !inherits(connection, "connection"))
        stop("'connection' must be a connection")
    .Internal(unserialize(connection, refhook))
}
#  File src/library/base/R/sets.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

## See the help for why as.vector is used:
## it includes coercing factors.
union <- function(x, y) unique(c(as.vector(x), as.vector(y)))

intersect <- function(x, y)
{
    y <- as.vector(y)
    unique(y[match(as.vector(x), y, 0L)])
}

setdiff <- function(x, y)
{
    x <- as.vector(x)
    y <- as.vector(y)
    unique(if(length(x) || length(y)) x[match(x, y, 0L) == 0L] else x)
}
## Faster versions, see R-devel, Jan.4-6, 2000;  optimize later...
setequal <- function(x, y)
{
    x <- as.vector(x)
    y <- as.vector(y)
    all(c(match(x, y, 0L) > 0L, match(y, x, 0L) > 0L))
}

##  same as %in% ( ./match.R ) but different arg names:
is.element <- function(el, set) match(el, set, 0L) > 0L
#  File src/library/base/R/sink.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

sink <- function(file=NULL, append = FALSE, type = c("output", "message"),
                 split=FALSE)
{
    type <- match.arg(type)
    if(type == "message") {
        if(is.null(file)) file <- stderr()
        else if(!inherits(file, "connection") || !isOpen(file))
           stop("'file' must be NULL or an already open connection")
        if (split) stop("cannot split the message connection")
        .Internal(sink(file, FALSE, TRUE, FALSE))
    } else {
        closeOnExit <- FALSE
        if(is.null(file)) file <- -1L
        else if(is.character(file)) {
            file <- file(file, ifelse(append, "a", "w"))
            closeOnExit <- TRUE
        } else if(!inherits(file, "connection"))
            stop("'file' must be NULL, a connection or a character string")
        .Internal(sink(file, closeOnExit, FALSE,split))
    }
}

sink.number <- function(type = c("output", "message"))
{
    type <- match.arg(type)
    .Internal(sink.number(type != "message"))
}
#  File src/library/base/R/solve.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2013 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

solve.qr <- function(a, b, ...)
{
    if( !is.qr(a) )
	stop("this is the \"qr\" method for the generic function solve()")
    nc <- ncol(a$qr); nr <- nrow(a$qr)
    if( a$rank != min(nc, nr) )
    if( a$rank != nc )
	stop("singular matrix 'a' in 'solve'")
    if( missing(b) ) {
	if( nc != nr )
	    stop("only square matrices can be inverted")
	b <- diag(1, nc)
    }
    res <- qr.coef(a, b)
    res[is.na(res)] <- 0
    res
}

solve.default <-
    function(a, b, tol = .Machine$double.eps, LINPACK = FALSE, ...)
{
    if(is.complex(a) || (!missing(b) && is.complex(b))) {
	a <- as.matrix(a)
	if(missing(b)) {
	    b <- diag(1.0+0.0i, nrow(a))
	    colnames(b) <- rownames(a)
	}
        return(.Internal(La_solve_cmplx(a, b)))
    }

    if(is.qr(a)) {
	warning("solve.default called with a \"qr\" object: use 'qr.solve'")
	return(solve.qr(a, b, tol))
    }

    if(LINPACK)
        warning("LINPACK = TRUE is defunct and will be ignored", domain = NA)
    a <- as.matrix(a)
    if(missing(b)) {
        b <- diag(1.0, nrow(a))
        colnames(b) <- rownames(a)
    }
    .Internal(La_solve(a, b, tol))
}

solve <- function(a, b, ...) UseMethod("solve")

qr.solve <- function(a, b, tol = 1e-7)
{
    if( !is.qr(a) )
	a <- qr(a, tol = tol)
    nc <- ncol(a$qr); nr <- nrow(a$qr)
    if( a$rank != min(nc, nr) )
	stop("singular matrix 'a' in solve")
    if( missing(b) ) {
	if( nc != nr )
	    stop("only square matrices can be inverted")
	b <- diag(1, nc)
    }
    res <- qr.coef(a, b)
    res[is.na(res)] <- 0
    res
}

#  File src/library/base/R/sort.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

sort <- function(x, decreasing = FALSE, ...)
{
    if(!is.logical(decreasing) || length(decreasing) != 1L)
        stop("'decreasing' must be a length-1 logical vector.\nDid you intend to set 'partial'?")
    UseMethod("sort")
}

sort.default <- function(x, decreasing = FALSE, na.last = NA, ...)
{
    ## The first case includes factors.
    if(is.object(x)) x[order(x, na.last = na.last, decreasing = decreasing)]
    else sort.int(x, na.last = na.last, decreasing = decreasing, ...)
}

sort.int <-
    function(x, partial = NULL, na.last = NA, decreasing = FALSE,
             method = c("shell", "quick"), index.return = FALSE)
{
    if(isfact <- is.factor(x)) {
        if(index.return) stop("'index.return' only for non-factors")
	lev <- levels(x)
	nlev <- nlevels(x)
 	isord <- is.ordered(x)
        x <- c(x) # drop attributes
    } else if(!is.atomic(x))
        stop("'x' must be atomic")

    if(has.na <- any(ina <- is.na(x))) {
        nas <- x[ina]
        x <-  x[!ina]
    }
    if(index.return && !is.na(na.last))
        stop("'index.return' only for 'na.last = NA'")
    if(!is.null(partial)) {
        if(index.return || decreasing || isfact || !missing(method))
	    stop("unsupported options for partial sorting")
        if(!all(is.finite(partial))) stop("non-finite 'partial'")
        y <- if(length(partial) <= 10L) {
            partial <- .Internal(qsort(partial, FALSE))
            .Internal(psort(x, partial))
        } else if (is.double(x)) .Internal(qsort(x, FALSE))
        else .Internal(sort(x, FALSE))
    } else if(isfact && missing(method) && nlev < 100000) {
        o <- sort.list(x, decreasing = decreasing, method = "radix")
        y <- x[o]
    } else {
        nms <- names(x)
        method <- if(is.numeric(x)) match.arg(method) else "shell"
        switch(method,
               "quick" = {
                   if(!is.null(nms)) {
                       if(decreasing) x <- -x
                       y <- .Internal(qsort(x, TRUE))
                       if(decreasing) y$x <- -y$x
                       names(y$x) <- nms[y$ix]
                       if (!index.return) y <- y$x
                   } else {
                       if(decreasing) x <- -x
                       y <- .Internal(qsort(x, index.return))
                       if(decreasing)
                           if(index.return) y$x <- -y$x else y <- -y
                   }
               },
               "shell" = {
                   if(index.return || !is.null(nms)) {
                       o <- sort.list(x, decreasing = decreasing)
                       y <- if (index.return) list(x = x[o], ix = o) else x[o]
                   }
                   else
                       y <- .Internal(sort(x, decreasing))
               })
    }
    if(!is.na(na.last) && has.na)
	y <- if(!na.last) c(nas, y) else c(y, nas)
    if(isfact)
        y <- (if (isord) ordered else factor)(y, levels = seq_len(nlev),
                                              labels = lev)
    y
}

order <- function(..., na.last = TRUE, decreasing = FALSE)
{
    z <- list(...)
    if(any(unlist(lapply(z, is.object)))) {
        z <- lapply(z, function(x) if(is.object(x)) xtfrm(x) else x)
        if(!is.na(na.last))
            return(do.call("order", c(z, na.last = na.last,
                                      decreasing = decreasing)))
    } else if(!is.na(na.last)) {
        if (length(z) == 1L && is.factor(zz <- z[[1L]]) && nlevels(zz) < 100000)
            return(.Internal(radixsort(zz, na.last, decreasing)))
        else return(.Internal(order(na.last, decreasing, ...)))
    }

    ## na.last = NA case: remove nas
    if(any(diff(l.z <- vapply(z, length, 1L)) != 0L))
        stop("argument lengths differ")
    ans <- vapply(z, is.na, rep.int(NA, l.z[1L]))
    ok <- if(is.matrix(ans)) !apply(ans, 1, any) else !any(ans)
    if(all(!ok)) return(integer())
    z[[1L]][!ok] <- NA
    ans <- do.call("order", c(z, decreasing = decreasing))
    keep <- seq_along(ok)[ok]
    ans[ans %in% keep]
}

sort.list <- function(x, partial = NULL, na.last = TRUE, decreasing = FALSE,
                      method = c("shell", "quick", "radix"))
{
    if (missing(method) && is.factor(x) && nlevels(x) < 100000) method <-"radix"
    method <- match.arg(method)
    if(!is.atomic(x))
        stop("'x' must be atomic for 'sort.list'\nHave you called 'sort' on a list?")
    if(!is.null(partial))
        .NotYetUsed("partial != NULL")
    if(method == "quick") {
        if(is.factor(x)) x <- as.integer(x) # sort the internal codes
        if(is.numeric(x))
            return(sort(x, na.last = na.last, decreasing = decreasing,
                        method = "quick", index.return = TRUE)$ix)
        else stop("method = \"quick\" is only for numeric 'x'")
    }
    if(method == "radix") {
        if(!typeof(x) == "integer") # we do want to allow factors here
            stop("method = \"radix\" is only for integer 'x'")
        if(is.na(na.last))
            return(.Internal(radixsort(x[!is.na(x)], TRUE, decreasing)))
        else
            return(.Internal(radixsort(x, na.last, decreasing)))
    }
    ## method == "shell"
    if(is.na(na.last)) .Internal(order(TRUE, decreasing, x[!is.na(x)]))
    else .Internal(order(na.last, decreasing, x))
}


## xtfrm is now primitive
## xtfrm <- function(x) UseMethod("xtfrm")
xtfrm.default <- function(x)
    if(is.numeric(x)) unclass(x) else as.vector(rank(x, ties.method="min", na.last="keep"))
xtfrm.factor <- function(x) as.integer(x) # primitive, so needs a wrapper
xtfrm.Surv <- function(x)
    if(ncol(x) == 2L) order(x[,1L], x[,2L]) else order(x[,1L], x[,2L], x[,3L]) # needed by 'party'
xtfrm.AsIs <- function(x)
{
    if(length(cl <- class(x)) > 1) oldClass(x) <- cl[-1L]
    NextMethod("xtfrm")
}

## callback from C code for rank/order
.gt <- function(x, i, j)
{
    xi <- x[i]; xj <- x[j]
    if (xi == xj) 0L else if(xi > xj) 1L else -1L;
}

## callback for C code for is.unsorted, hence negation.
.gtn <- function(x, strictly)
{
    n <- length(x)
    if(strictly) !all(x[-1L] > x[-n]) else !all(x[-1L] >= x[-n])
}
#  File src/library/base/R/source.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

source <-
function(file, local = FALSE, echo = verbose, print.eval = echo,
	 verbose = getOption("verbose"),
	 prompt.echo = getOption("prompt"),
	 max.deparse.length = 150, chdir = FALSE,
         encoding = getOption("encoding"),
         continue.echo = getOption("continue"),
         skip.echo = 0, keep.source = getOption("keep.source"))
{
    envir <- if (isTRUE(local)) {
        parent.frame()
    } else if(identical(local, FALSE)) {
        .GlobalEnv
    } else if (is.environment(local)) {
        local
    } else stop("'local' must be TRUE, FALSE or an environment")
    have_encoding <- !missing(encoding) && encoding != "unknown"
    if (!missing(echo)) {
	if (!is.logical(echo))
	    stop("'echo' must be logical")
	if (!echo && verbose) {
	    warning("'verbose' is TRUE, 'echo' not; ... coercing 'echo <- TRUE'")
	    echo <- TRUE
	}
    }
    if (verbose) {
	cat("'envir' chosen:")
	print(envir)
    }
    ofile <- file # for use with chdir = TRUE
    from_file <- FALSE
    srcfile <- NULL
    if(is.character(file)) {
        if(identical(encoding, "unknown")) {
            enc <- utils::localeToCharset()
            encoding <- enc[length(enc)]
        } else enc <- encoding
        if(length(enc) > 1L) {
            encoding <- NA
            owarn <- options("warn"); options(warn = 2)
            for(e in enc) {
                if(is.na(e)) next
                zz <- file(file, encoding = e)
                res <- tryCatch(readLines(zz, warn = FALSE), error = identity)
                close(zz)
                if(!inherits(res, "error")) { encoding <- e; break }
            }
            options(owarn)
        }
        if(is.na(encoding))
            stop("unable to find a plausible encoding")
        if(verbose)
            cat(gettextf('encoding = "%s" chosen', encoding), "\n", sep = "")
        if(file == "") {
	    file <- stdin()
	    srcfile <- "<stdin>"
        } else {
            filename <- file
	    file <- file(filename, "r", encoding = encoding)
	    on.exit(close(file))
            if (isTRUE(keep.source)) {
	    	lines <- readLines(file, warn = FALSE)
	    	on.exit()
	    	close(file)
            	srcfile <- srcfilecopy(filename, lines, file.info(filename)[1,"mtime"],
            			       isFile = TRUE)
	    } else {
            	from_file <- TRUE
		srcfile <- filename
	    }

            ## We translated the file (possibly via a guess),
            ## so don't want to mark the strings.as from that encoding
            ## but we might know what we have encoded to, so
            loc <- utils::localeToCharset()[1L]
            encoding <- if(have_encoding)
                switch(loc,
                       "UTF-8" = "UTF-8",
                       "ISO8859-1" = "latin1",
                       "unknown")
            else "unknown"
	}
    } else {
    	lines <- readLines(file, warn = FALSE)
    	if (isTRUE(keep.source))
    	    srcfile <- srcfilecopy(deparse(substitute(file)), lines)
    	else
    	    srcfile <- deparse(substitute(file))
    }

    exprs <- if (!from_file) {
        if (length(lines))  # there is a C-level test for this
            .Internal(parse(stdin(), n = -1, lines, "?", srcfile, encoding))
        else expression()
    } else
    	.Internal(parse(file, n = -1, NULL, "?", srcfile, encoding))

    on.exit()
    if (from_file) close(file)

    Ne <- length(exprs)
    if (verbose)
	cat("--> parsed", Ne, "expressions; now eval(.)ing them:\n")

    if (chdir){
        if(is.character(ofile)) {
            isURL <- length(grep("^(ftp|http|file)://", ofile)) > 0L
            if(isURL)
                warning("'chdir = TRUE' makes no sense for a URL")
            if(!isURL && (path <- dirname(ofile)) != ".") {
                owd <- getwd()
                if(is.null(owd))
                    stop("cannot 'chdir' as current directory is unknown")
                on.exit(setwd(owd), add=TRUE)
                setwd(path)
            }
        } else {
            warning("'chdir = TRUE' makes no sense for a connection")
        }
    }

    if (echo) {
	## Reg.exps for string delimiter/ NO-string-del /
	## odd-number-of-str.del needed, when truncating below
	sd <- "\""
	nos <- "[^\"]*"
	oddsd <- paste0("^", nos, sd, "(", nos, sd, nos, sd, ")*", nos, "$")
        ## A helper function for echoing source.  This is simpler than the
        ## same-named one in Sweave
	trySrcLines <- function(srcfile, showfrom, showto) {
	    lines <- try(suppressWarnings(getSrcLines(srcfile, showfrom, showto)), silent=TRUE)
	    if (inherits(lines, "try-error")) lines <- character()
    	    lines
	}
    }
    yy <- NULL
    lastshown <- 0
    srcrefs <- attr(exprs, "srcref")
    for (i in seq_len(Ne+echo)) {
    	tail <- i > Ne
        if (!tail) {
	    if (verbose)
		cat("\n>>>> eval(expression_nr.", i, ")\n\t	 =================\n")
	    ei <- exprs[i]
	}
	if (echo) {
	    nd <- 0
	    srcref <- if(tail) attr(exprs, "wholeSrcref") else
		if(i <= length(srcrefs)) srcrefs[[i]] # else NULL
 	    if (!is.null(srcref)) {
	    	if (i == 1) lastshown <- min(skip.echo, srcref[3L]-1)
	    	if (lastshown < srcref[3L]) {
	    	    srcfile <- attr(srcref, "srcfile")
	    	    dep <- trySrcLines(srcfile, lastshown+1, srcref[3L])
	    	    if (length(dep)) {
			leading <- if(tail) length(dep) else srcref[1L]-lastshown
			lastshown <- srcref[3L]
			while (length(dep) && length(grep("^[[:blank:]]*$", dep[1L]))) {
			    dep <- dep[-1L]
			    leading <- leading - 1L
			}
			dep <- paste0(rep.int(c(prompt.echo, continue.echo),
					      c(leading, length(dep)-leading)),
				      dep, collapse="\n")
			nd <- nchar(dep, "c")
		    } else
		    	srcref <- NULL  # Give up and deparse
	    	}
	    }
	    if (is.null(srcref)) {
	    	if (!tail) {
		    # Deparse.  Must drop "expression(...)"
		    dep <- substr(paste(deparse(ei, control = "showAttributes"),
					collapse = "\n"), 12L, 1e+06L)
		    ## We really do want chars here as \n\t may be embedded.
		    dep <- paste0(prompt.echo,
				  gsub("\n", paste0("\n", continue.echo), dep))
		    nd <- nchar(dep, "c") - 1L
		}
	    }
	    if (nd) {
		do.trunc <- nd > max.deparse.length
		dep <- substr(dep, 1L, if (do.trunc) max.deparse.length else nd)
		cat("\n", dep, if (do.trunc)
		    paste(if (length(grep(sd, dep)) && length(grep(oddsd, dep)))
			  " ...\" ..." else " ....", "[TRUNCATED] "),
		    "\n", sep = "")
	    }
	}
	if (!tail) {
	    yy <- withVisible(eval(ei, envir))
	    i.symbol <- mode(ei[[1L]]) == "name"
	    if (!i.symbol) {
		## ei[[1L]] : the function "<-" or other
		curr.fun <- ei[[1L]][[1L]]
		if (verbose) {
		    cat("curr.fun:")
		    utils::str(curr.fun)
		}
	    }
	    if (verbose >= 2) {
		cat(".... mode(ei[[1L]])=", mode(ei[[1L]]), "; paste(curr.fun)=")
		utils::str(paste(curr.fun))
	    }
	    if (print.eval && yy$visible) {
		if(isS4(yy$value))
		    methods::show(yy$value)
		else
		    print(yy$value)
	    }
	    if (verbose)
		cat(" .. after ", sQuote(deparse(ei,
		    control = c("showAttributes","useSource"))), "\n", sep = "")
 	}
    }
    invisible(yy)
}

sys.source <-
function(file, envir = baseenv(), chdir = FALSE,
	 keep.source = getOption("keep.source.pkgs"))
{
    if(!(is.character(file) && file.exists(file)))
	stop(gettextf("'%s' is not an existing file", file))
    keep.source <- as.logical(keep.source)
    oop <- options(keep.source = keep.source,
		   topLevelEnvironment = as.environment(envir))
    on.exit(options(oop))
    if (keep.source) {
    	lines <- readLines(file, warn = FALSE)
    	srcfile <- srcfilecopy(file, lines, file.info(file)[1,"mtime"], isFile = TRUE)
    	exprs <- parse(text = lines, srcfile = srcfile, keep.source = TRUE)
    } else
    	exprs <- parse(n = -1, file = file, srcfile = NULL, keep.source = FALSE)
    if (length(exprs) == 0L)
	return(invisible())
    if (chdir && (path <- dirname(file)) != ".") {
	owd <- getwd()
        if(is.null(owd))
            stop("cannot 'chdir' as current directory is unknown")
	on.exit(setwd(owd), add = TRUE)
	setwd(path)
    }
    for (i in exprs) eval(i, envir)
    invisible()
}
#  File src/library/base/R/split.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2013 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

split <- function(x, f, drop = FALSE, ...) UseMethod("split")

split.default <- function(x, f, drop = FALSE, sep = ".", ...)
{
    if(!missing(...)) .NotYetUsed(deparse(...), error = FALSE)

    if (is.list(f)) f <- interaction(f, drop = drop, sep = sep)
    else if (!is.factor(f)) f <- as.factor(f) # docs say as.factor
    else if (drop) f <- factor(f) # drop extraneous levels
    storage.mode(f) <- "integer"  # some factors have had double in the past
    if (is.null(attr(x, "class")))
	return(.Internal(split(x, f)))
    ## else
    lf <- levels(f)
    y <- vector("list", length(lf))
    names(y) <- lf
    ind <- .Internal(split(seq_along(x), f))
    for(k in lf) y[[k]] <- x[ind[[k]]]
    y
}

## This is documented to work for matrices too
split.data.frame <- function(x, f, drop = FALSE, ...)
    lapply(split(x = seq_len(nrow(x)), f = f, drop = drop, ...),
           function(ind) x[ind, , drop = FALSE])

`split<-` <- function(x, f, drop = FALSE, ..., value) UseMethod("split<-")

`split<-.default` <- function(x, f, drop = FALSE, ..., value)
{
    ix <- split(seq_along(x), f, drop = drop, ...)
    n <- length(value)
    j <- 0
    for (i in ix) {
        j <- j %% n + 1
        x[i] <- value[[j]]
    }
    x
}

## This is documented to work for matrices too
`split<-.data.frame` <- function(x, f, drop = FALSE, ..., value)
{
    ix <- split(seq_len(nrow(x)), f, drop = drop, ...)
    n <- length(value)
    j <- 0
    for (i in ix) {
        j <- j %% n + 1
        x[i,] <- value[[j]]
    }
    x
}

unsplit <- function (value, f, drop = FALSE)
{
    len <- length(if (is.list(f)) f[[1L]] else f)
    if (is.data.frame(value[[1L]])) {
        x <- value[[1L]][rep(NA, len),, drop = FALSE]
        rownames(x) <- unsplit(lapply(value, rownames), f, drop = drop)
    } else
        x <- value[[1L]][rep(NA, len)]
    split(x, f, drop = drop) <- value
    x
}
#  File src/library/base/R/srcfile.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

# a srcfile is a file with a timestamp

srcfile <- function(filename, encoding = getOption("encoding"), Enc = "unknown")
{
    stopifnot(is.character(filename), length(filename) == 1L)

    ## This is small, no need to hash.
    e <- new.env(hash = FALSE, parent = emptyenv())

    e$wd <- getwd()
    e$filename <- filename

    # If filename is a URL, this will return NA
    e$timestamp <- file.info(filename)[1,"mtime"]

    if (identical(encoding, "unknown")) encoding <- "native.enc"
    e$encoding <- encoding
    e$Enc <- Enc

    class(e) <- "srcfile"
    return(e)
}

print.srcfile <- function(x, ...) {
    cat(x$filename, "\n")
    invisible(x)
}

summary.srcfile <- function(object, ...) {
    cat(utils:::.normalizePath(object$filename, object$wd), "\n")

    if (inherits(object$timestamp, "POSIXt"))
    	cat("Timestamp: ", format(object$timestamp, usetz=TRUE), "\n", sep="")

    cat('Encoding: "', object$encoding, '"', sep="")
    if (!is.null(object$Enc) && object$Enc != object$encoding && object$Enc != "unknown")
    	cat(', re-encoded to "', object$Enc, '"', sep="")
    cat("\n")

    invisible(object)
}

open.srcfile <- function(con, line, ...) {

    srcfile <- con

    oldline <- srcfile$line
    if (!is.null(oldline) && oldline > line) close(srcfile)

    conn <- srcfile$conn
    if (is.null(conn)) {
        if (!is.null(srcfile$wd)) {
	    olddir <- setwd(srcfile$wd)
	    on.exit(setwd(olddir))
	}
	timestamp <- file.info(srcfile$filename)[1,"mtime"]
	if (!is.null(srcfile$timestamp)
	    && !is.na(srcfile$timestamp)
	    && ( is.na(timestamp) || timestamp != srcfile$timestamp) )
            warning(gettextf("Timestamp of %s has changed", sQuote(srcfile$filename)),
                    call. = FALSE, domain = NA)
	if (is.null(srcfile$encoding)) encoding <- getOption("encoding")
	else encoding <- srcfile$encoding
	# Specifying encoding below means that reads will convert to the native encoding
	srcfile$conn <- conn <- file(srcfile$filename, open="rt", encoding=encoding)
	srcfile$line <- 1L
	oldline <- 1L
    } else if (!isOpen(conn)) {
	open(conn, open="rt")
	srcfile$line <- 1
	oldline <- 1L
    }
    if (oldline < line) {
	readLines(conn, line - oldline, warn = FALSE)
	srcfile$line <- line
    }
    invisible(conn)
}

close.srcfile <- function(con, ...) {
    srcfile <- con
    conn <- srcfile$conn
    if (is.null(conn)) return()
    else {
	close(conn)
	rm(list=c("conn", "line"), envir=srcfile)
    }
}

# srcfilecopy saves a copy of lines from a file

srcfilecopy <- function(filename, lines, timestamp = Sys.time(), isFile = FALSE) {
    stopifnot(is.character(filename), length(filename) == 1L)

    e <- new.env(parent=emptyenv())

    # Remove embedded newlines
    if (any(grepl("\n", lines, fixed=TRUE)))
	lines <- unlist(strsplit(sub("$", "\n", as.character(lines)), "\n"))

    e$filename <- filename
    e$wd <- getwd()
    e$isFile <- isFile
    e$lines <- as.character(lines)
    e$fixedNewlines <- TRUE  	# we have removed the newlines already
    e$timestamp <- timestamp
    e$Enc <- "unknown"

    class(e) <- c("srcfilecopy", "srcfile")
    return(e)
}

open.srcfilecopy <- function(con, line, ...) {

    srcfile <- con

    oldline <- srcfile$line
    if (!is.null(oldline) && oldline > line) close(srcfile)

    conn <- srcfile$conn
    if (is.null(conn)) {
	srcfile$conn <- conn <- textConnection(srcfile$lines, open="r")
	srcfile$line <- 1L
	oldline <- 1L
    } else if (!isOpen(conn)) {
	open(conn, open="r")
	srcfile$line <- 1L
	oldline <- 1L
    }
    if (oldline < line) {
	readLines(conn, line - oldline, warn = FALSE)
	srcfile$line <- line
    }
    invisible(conn)
}

srcfilealias <- function(filename, srcfile) {
    stopifnot(is.character(filename), length(filename) == 1L)

    e <- new.env(parent=emptyenv())

    e$filename <- filename
    e$original <- srcfile

    class(e) <- c("srcfilealias", "srcfile")
    return(e)
}

open.srcfilealias <- function(con, line, ...)
    open(con$original, line, ...)

close.srcfilealias <- function(con, ...)
    close(con$original, ...)

.isOpen <- function(srcfile) {
    conn <- srcfile$conn
    return( !is.null(conn) && isOpen(conn) )
}

getSrcLines <- function(srcfile, first, last) {
    if (first > last) return(character())
    if (inherits(srcfile, "srcfilealias"))
    	srcfile <- srcfile$original
    if (inherits(srcfile, "srcfilecopy")) {
	# Remove embedded newlines if we haven't done this already
	if (is.null(srcfile$fixedNewlines)) {
	    lines <- srcfile$lines
    	    if (any(grepl("\n", lines, fixed=TRUE)))
		srcfile$lines <- unlist(strsplit(sub("$", "\n", as.character(lines)), "\n"))
	    srcfile$fixedNewlines <- TRUE
	}
        last <- min(last, length(srcfile$lines))
        if (first > last) return(character())
        else return(srcfile$lines[first:last])
    }
    if (!.isOpen(srcfile)) on.exit(close(srcfile))
    conn <- open(srcfile, first)
    lines <- readLines(conn, n = last - first + 1L, warn = FALSE)
    # Re-encode from native encoding to specified one
    if (!is.null(Enc <- srcfile$Enc) && !(Enc %in% c("unknown", "native.enc")))
    	lines <- iconv(lines, "", Enc)
    srcfile$line <- first + length(lines)
    return(lines)
}

# a srcref gives start and stop positions of text
# lloc entries are first_line, first_byte, last_line, last_byte,
#  first_column, last_column, first_parse, last_parse
# all are inclusive

srcref <- function(srcfile, lloc) {
    stopifnot(inherits(srcfile, "srcfile"), length(lloc) %in% c(4L,6L,8L))
    if (length(lloc) == 4) lloc <- c(lloc, lloc[c(2,4)])
    if (length(lloc) == 6) lloc <- c(lloc, lloc[c(1,3)])
    structure(as.integer(lloc), srcfile=srcfile, class="srcref")
}

as.character.srcref <- function(x, useSource = TRUE, ...)
{
    srcfile <- attr(x, "srcfile")
    if (!is.null(srcfile) && !inherits(srcfile, "srcfile")) {
       cat("forcing class on") ## debug
	print(str(srcfile))
       class(srcfile) <- c("srcfilealias", "srcfile")
    }
    if (useSource) {
    	if (inherits(srcfile, "srcfilecopy") || inherits(srcfile, "srcfilealias"))
    	    lines <- try(getSrcLines(srcfile, x[7L], x[8L]), TRUE)
    	else
 	    lines <- try(getSrcLines(srcfile, x[1L], x[3L]), TRUE)
    }
    if (!useSource || inherits(lines, "try-error"))
    	lines <- paste("<srcref: file \"", srcfile$filename, "\" chars ",
                       x[1L],":",x[5L], " to ",x[3L],":",x[6L], ">", sep="")
    else if (length(lines)) {
    	enc <- Encoding(lines)
    	Encoding(lines) <- "latin1"  # so byte counting works
        if (length(lines) < x[3L] - x[1L] + 1L)
            x[4L] <- .Machine$integer.max
    	lines[length(lines)] <- substring(lines[length(lines)], 1L, x[4L])
    	lines[1L] <- substring(lines[1L], x[2L])
    	Encoding(lines) <- enc
    }
    lines
}

print.srcref <- function(x, useSource = TRUE, ...) {
    cat(as.character(x, useSource = useSource), sep="\n")
    invisible(x)
}

summary.srcref <- function(object, useSource = FALSE, ...) {
    cat(as.character(object, useSource = useSource), sep="\n")
    invisible(object)
}
#  File src/library/base/R/stop.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

stop <- function(..., call. = TRUE, domain = NULL)
{
    args <- list(...)
    if (length(args) == 1L && inherits(args[[1L]], "condition")) {
        cond <- args[[1L]]
        if(nargs() > 1L)
            warning("additional arguments ignored in stop()")
        message <- conditionMessage(cond)
        call <- conditionCall(cond)
        .Internal(.signalCondition(cond, message, call))
        .Internal(.dfltStop(message, call))
    } else
        .Internal(stop(as.logical(call.), .makeMessage(..., domain = domain)))
}

stopifnot <- function(...)
{
    n <- length(ll <- list(...))
    if(n == 0L)
	return(invisible())
    mc <- match.call()
    for(i in 1L:n)
	if(!(is.logical(r <- ll[[i]]) && !any(is.na(r)) && all(r))) {
	    ch <- deparse(mc[[i+1]], width.cutoff = 60L)
	    if(length(ch) > 1L) ch <- paste(ch[1L], "....")
            stop(sprintf(ngettext(length(r),
                                  "%s is not TRUE",
                                  "%s are not all TRUE"),
                         ch), call. = FALSE, domain = NA)
	}
    invisible()
}

warning <- function(..., call. = TRUE, immediate. = FALSE, domain = NULL)
{
    args <- list(...)
    if (length(args) == 1L && inherits(args[[1L]], "condition")) {
        cond <- args[[1L]]
        if(nargs() > 1L)
            cat(gettext("additional arguments ignored in warning()"),
                "\n", sep = "", file = stderr())
        message <- conditionMessage(cond)
        call <- conditionCall(cond)
        withRestarts({
                .Internal(.signalCondition(cond, message, call))
                .Internal(.dfltWarn(message, call))
            }, muffleWarning = function() NULL) #**** allow simpler form??
        invisible(message)
    } else
        .Internal(warning(as.logical(call.), as.logical(immediate.),
                          .makeMessage(..., domain = domain)))
}

gettext <- function(..., domain = NULL) {
    args <- lapply(list(...), as.character)
    .Internal(gettext(domain, unlist(args)))
}

bindtextdomain <- function(domain, dirname = NULL)
    .Internal(bindtextdomain(domain, dirname))

ngettext <- function(n, msg1, msg2, domain = NULL)
    .Internal(ngettext(n, msg1, msg2, domain))

gettextf <- function(fmt, ..., domain = NULL)
    sprintf(gettext(fmt, domain = domain), ...)
#  File src/library/base/R/structure.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

## This remaps special names as they are used by deparsing, but why are they?
##
structure <- function (.Data, ...)
{
    attrib <- list(...)
    if(length(attrib)) {
        specials <- c(".Dim", ".Dimnames", ".Names", ".Tsp", ".Label")
        replace <- c("dim", "dimnames", "names", "tsp", "levels")
	m <- match(names(attrib), specials)
	ok <- (!is.na(m) & m)
	names(attrib)[ok] <- replace[m[ok]]
        ## prior to 2.5.0 factors would deparse to double codes
	if("factor" %in% attrib[["class", exact = TRUE]]
           && typeof(.Data) == "double")
	   storage.mode(.Data) <- "integer"
	attributes(.Data) <- c(attributes(.Data), attrib)
    }
    return(.Data)
}
#  File src/library/base/R/strwrap.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

strtrim <- function(x, width)
{
    if(!is.character(x)) x <- as.character(x)
    .Internal(strtrim(x, width))
}

strwrap <-
function(x, width = 0.9 * getOption("width"), indent = 0, exdent = 0,
         prefix = "", simplify = TRUE, initial = prefix)
{
    if(!is.character(x)) x <- as.character(x)
    ## Useful variables.
    indentString <- paste(rep.int(" ", indent), collapse = "")
    exdentString <- paste(rep.int(" ", exdent), collapse = "")
    y <- list()                         # return value
    UB <- TRUE
    ## input need not be valid in this locale, e.g. from write.dcf
    ## but if x has UTF-8 encoding we want to preserve it, so
    if(all(Encoding(x) == "UTF-8")) UB <- FALSE
    else {
        ## Let's convert anything else with a marked encoding
        ## to the current locale
        enc <- Encoding(x) %in% c("latin1", "UTF-8")
        if(length(enc)) x[enc] <- enc2native(x[enc])
    }
    z <- lapply(strsplit(x, "\n[ \t\n]*\n", perl = TRUE, useBytes = UB),
                strsplit, "[ \t\n]", perl = TRUE, useBytes = UB)
    ## Now z[[i]][[j]] is a character vector of all "words" in
    ## paragraph j of x[i].

    for(i in seq_along(z)) {
        yi <- character()
        for(j in seq_along(z[[i]])) {
            ## Format paragraph j in x[i].
            words <- z[[i]][[j]]
            nc <- nchar(words, type="w")
	    if(any(is.na(nc))) {
		## use byte count as a reasonable substitute
		nc0 <- nchar(words, type="b")
		nc[is.na(nc)] <- nc0[is.na(nc)]
	    }

            ## Remove extra white space unless after a period which
            ## hopefully ends a sentence.
            ## Add ? ! as other possible ends, and there might be
            ## quoted and parenthesised sentences.
            ## NB, input could be invalid here.
            if(any(nc == 0L)) {
                zLenInd <- which(nc == 0L)
                zLenInd <- zLenInd[!(zLenInd %in%
                                     (grep("[.?!][)\"']{0,1}$", words,
                                           perl = TRUE, useBytes = TRUE) + 1L))]
                if(length(zLenInd)) {
                    words <- words[-zLenInd]
                    nc <- nc[-zLenInd]
                }
            }

            if(!length(words)) {
                yi <- c(yi, "", initial)
                next
            }

            currentIndex <- 0L
            lowerBlockIndex <- 1L
            upperBlockIndex <- integer()
            lens <- cumsum(nc + 1L)

            first <- TRUE
            maxLength <- width - nchar(initial, type="w") - indent

            ## Recursively build a sequence of lower and upper indices
            ## such that the words in line k are the ones in the k-th
            ## index block.
            while(length(lens)) {
                k <- max(sum(lens <= maxLength), 1L)
                if(first) {
                    first <- FALSE
                    maxLength <- width - nchar(prefix, type="w") - exdent
                }
                currentIndex <- currentIndex + k
                if(nc[currentIndex] == 0L)
                    ## Are we sitting on a space?
                    upperBlockIndex <- c(upperBlockIndex,
                                         currentIndex - 1L)
                else
                    upperBlockIndex <- c(upperBlockIndex,
                                         currentIndex)
                if(length(lens) > k) {
                    ## Are we looking at a space?
                    if(nc[currentIndex + 1L] == 0L) {
                        currentIndex <- currentIndex + 1L
                        k <- k + 1L
                    }
                    lowerBlockIndex <- c(lowerBlockIndex,
                                         currentIndex + 1L)
                }
                if(length(lens) > k)
                    lens <- lens[-seq_len(k)] - lens[k]
                else
                    lens <- NULL
            }

            nBlocks <- length(upperBlockIndex)
	    s <- paste0(c(initial, rep.int(prefix, nBlocks - 1L)),
			c(indentString, rep.int(exdentString, nBlocks - 1L)))
            initial <- prefix
            for(k in seq_len(nBlocks))
		s[k] <- paste0(s[k], paste(words[lowerBlockIndex[k] :
						 upperBlockIndex[k]],
					   collapse = " "))

            yi <- c(yi, s, prefix)
        }
        y <- if(length(yi))
            c(y, list(yi[-length(yi)]))
        else
            c(y, "")
    }

    if(simplify) y <- as.character(unlist(y))
    y
}

formatDL <-
function(x, y, style = c("table", "list"),
         width = 0.9 * getOption("width"), indent = NULL)
{
    if(is.list(x)) {
        if(length(x) == 2L && diff(vapply(x, length, 1L)) == 0L) {
            y <- x[[2L]]; x <- x[[1L]]
        }
        else
            stop("incorrect value for 'x'")
    }
    else if(is.matrix(x)) {
        if(NCOL(x) == 2L) {
            y <- x[, 2L]; x <- x[, 1L]
        }
        else
            stop("incorrect value for 'x'")
    }
    else if(missing(y) && !is.null(nms <- names(x))) {
        y <- x
        x <- nms
    }
    else if(length(x) != length(y))
        stop("'x' and 'y' must have the same length")
    x <- as.character(x)
    if(!length(x)) return(x)
    y <- as.character(y)

    style <- match.arg(style)

    if(is.null(indent))
        indent <- switch(style, table = width / 3, list = width / 9)
    if(indent > 0.5 * width)
        stop("incorrect values of 'indent' and 'width'")

    indentString <- paste(rep.int(" ", indent), collapse = "")

    if(style == "table") {
        i <- (nchar(x, type="w") > indent - 3L)
        if(any(i))
            x[i] <- paste0(x[i], "\n", indentString)
        i <- !i
        if(any(i))
            x[i] <- formatC(x[i], width = indent, flag = "-")
        y <- lapply(strwrap(y, width = width - indent,
                            simplify = FALSE),
                    paste,
                    collapse = paste0("\n", indentString))
        r <- paste0(x, unlist(y))
    }
    else if(style == "list") {
        y <- strwrap(paste0(x, ": ", y), exdent = indent,
                     width = width, simplify = FALSE)
        r <- unlist(lapply(y, paste, collapse = "\n"))
    }
    r
}
#  File src/library/base/R/summary.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

summary <- function (object, ...) UseMethod("summary")

summary.default <-
    function(object, ..., digits = max(3L, getOption("digits") - 3L))
{
    if(is.factor(object))
	return(summary.factor(object, ...))
    else if(is.matrix(object))
	return(summary.matrix(object, digits = digits, ...))

    value <- if(is.logical(object)) # scalar or array!
	c(Mode = "logical",
          {tb <- table(object, exclude = NULL) # incl. NA s
           if(!is.null(n <- dimnames(tb)[[1L]]) && any(iN <- is.na(n)))
               dimnames(tb)[[1L]][iN] <- "NA's"
           tb
           })
    else if(is.numeric(object)) {
	nas <- is.na(object)
	object <- object[!nas]
	qq <- stats::quantile(object)
	qq <- signif(c(qq[1L:3L], mean(object), qq[4L:5L]), digits)
	names(qq) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
	if(any(nas))
	    c(qq, "NA's" = sum(nas))
	else qq
    } else if(is.recursive(object) && !is.language(object) &&
	      (n <- length(object))) { # do not allow long dims
	sumry <- array("", c(n, 3L), list(names(object),
                                          c("Length", "Class", "Mode")))
	ll <- numeric(n)
	for(i in 1L:n) {
	    ii <- object[[i]]
	    ll[i] <- length(ii)
	    cls <- oldClass(ii)
	    sumry[i, 2L] <- if(length(cls)) cls[1L] else "-none-"
	    sumry[i, 3L] <- mode(ii)
	}
	sumry[, 1L] <- format(as.integer(ll))
	sumry
    }
    else c(Length = length(object), Class = class(object), Mode = mode(object))
    class(value) <- c("summaryDefault", "table")
    value
}

format.summaryDefault <- function(x, ...)
{
    xx <- if(is.numeric(x) || is.complex(x)) zapsmall(x) else x
    class(xx) <- class(x)[-1]
    m <- match("NA's", names(x), 0)
    if(inherits(x, "Date") || inherits(x, "POSIXct")) {
        if(length(a <- attr(x, "NAs")))
            c(format(xx, ...), "NA's" = as.character(a))
        else format(xx)
    } else if(m && !is.character(x))
        xx <- c(format(xx[-m], ...), "NA's" = as.character(xx[m]))
    else format(xx, ...)
}

print.summaryDefault <- function(x, ...)
{
    xx <- if(is.numeric(x) || is.complex(x)) zapsmall(x) else x
    class(xx) <- class(x)[-1] # for format
    m <- match("NA's", names(xx), 0)
    if(inherits(x, "Date") || inherits(x, "POSIXct")) {
        xx <- if(length(a <- attr(x, "NAs")))
            c(format(xx), "NA's" = as.character(a))
        else format(xx)
        print(xx, ...)
        return(invisible(x))
    } else if(m && !is.character(x))
        xx <- c(format(xx[-m]), "NA's" = as.character(xx[m]))
    print.table(xx, ...)
    invisible(x)
}

summary.factor <- function(object, maxsum = 100, ...)
{
    nas <- is.na(object)
    ll <- levels(object)
    if(any(nas)) maxsum <- maxsum - 1
    tbl <- table(object)
    tt <- c(tbl) # names dropped ...
    names(tt) <- dimnames(tbl)[[1L]]
    if(length(ll) > maxsum) {
	drop <- maxsum:length(ll)
	o <- sort.list(tt, decreasing = TRUE)
	tt <- c(tt[o[ - drop]], "(Other)" = sum(tt[o[drop]]))
    }
    if(any(nas)) c(tt, "NA's" = sum(nas)) else tt
}

summary.matrix <- function(object, ...) {
    ## we do want this changed into separate columns, so use matrix method
    summary.data.frame(as.data.frame.matrix(object), ...)
}

summary.data.frame <-
    function(object, maxsum = 7L, digits = max(3L, getOption("digits") - 3L), ...)
{
    ncw <- function(x) {
        z <- nchar(x, type="w")
        if (any(na <- is.na(z))) {
            # FIXME: can we do better
            z[na] <- nchar(encodeString(z[na]), "b")
        }
        z
    }
    # compute results to full precision.
    z <- lapply(X = as.list(object), FUN = summary,
                maxsum = maxsum, digits = 12L, ...)
    nv <- length(object)
    nm <- names(object)
    lw <- numeric(nv)
    nr <- max(unlist(lapply(z, NROW)))
    for(i in 1L:nv) {
        sms <- z[[i]]
        if(is.matrix(sms)) {
            ## need to produce a single column, so collapse matrix
            ## across rows
            cn <- paste(nm[i], gsub("^ +", "", colnames(sms), useBytes = TRUE),
                        sep=".")
            tmp <- format(sms)
            if(nrow(sms) < nr)
                tmp <- rbind(tmp, matrix("", nr - nrow(sms), ncol(sms)))
            sms <- apply(tmp, 1L, function(x) paste(x, collapse="  "))
            ## produce a suitable colname: undoing padding
            wid <- sapply(tmp[1L, ], nchar, type="w") # might be NA
            blanks <- paste(character(max(wid)), collapse = " ")
            wcn <- ncw(cn)
            pad0 <- floor((wid - wcn)/2)
            pad1 <- wid - wcn - pad0
            cn <- paste0(substring(blanks, 1L, pad0), cn,
                         substring(blanks, 1L, pad1))
            nm[i] <- paste(cn, collapse="  ")
            z[[i]] <- sms
        } else {
            sms <- format(sms, digits = digits) # may add NAs row
            lbs <- format(names(sms))
            sms <- paste0(lbs, ":", sms, "  ")
            lw[i] <- ncw(lbs[1L])
            length(sms) <- nr
            z[[i]] <- sms
        }
    }
    z <- unlist(z, use.names=TRUE)
    dim(z) <- c(nr, nv)
    if(any(is.na(lw)))
	warning("probably wrong encoding in names(.) of column ",
		paste(which(is.na(lw)), collapse = ", "))
    blanks <- paste(character(max(lw, na.rm=TRUE) + 2L), collapse = " ")
    pad <- floor(lw - ncw(nm)/2)
    nm <- paste0(substring(blanks, 1, pad), nm)
    dimnames(z) <- list(rep.int("", nr), nm)
    attr(z, "class") <- c("table") #, "matrix")
    z
}
#  File src/library/base/R/svd.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2013 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

svd <- function(x, nu = min(n,p), nv = min(n,p), LINPACK = FALSE)
{
    if(LINPACK)
        warning("LINPACK = TRUE is defunct and will be ignored", domain = NA)
    x <- as.matrix(x)
    if (any(!is.finite(x))) stop("infinite or missing values in 'x'")
    dx <- dim(x)
    n <- dx[1L]
    p <- dx[2L]
    if(!n || !p) stop("a dimension is zero")
    La.res <- La.svd(x, nu, nv)
    res <- list(d = La.res$d)
    if (nu) res$u <- La.res$u
    if (nv) {
	if (is.complex(x))
	    res$v <- Conj(t(La.res$vt))
	else
	    res$v <- t(La.res$vt)
    }
    res
}
#  File src/library/base/R/sweep.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

sweep <- function(x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ...)
{
    FUN <- match.fun(FUN)
    dims <- dim(x)
    if (check.margin) {
        dimmargin <- dims[MARGIN]
        dimstats <- dim(STATS)
        lstats <- length(STATS)
        if (lstats > prod(dimmargin)) {
            warning("STATS is longer than the extent of 'dim(x)[MARGIN]'")
        } else if (is.null(dimstats)) { # STATS is a vector
            cumDim <- c(1L, cumprod(dimmargin))
            upper <- min(cumDim[cumDim >= lstats])
            lower <- max(cumDim[cumDim <= lstats])
            if (lstats && (upper %% lstats != 0L || lstats %% lower != 0L))
                warning("STATS does not recycle exactly across MARGIN")
        } else {
            dimmargin <- dimmargin[dimmargin > 1L]
            dimstats <- dimstats[dimstats > 1L]
            if (length(dimstats) != length(dimmargin) ||
                any(dimstats != dimmargin))
                warning("length(STATS) or dim(STATS) do not match dim(x)[MARGIN]")
        }
    }
    perm <- c(MARGIN, seq_along(dims)[ - MARGIN])
    FUN(x, aperm(array(STATS, dims[perm]), order(perm)), ...)
}
#  File src/library/base/R/sys.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

sys.call <- function(which = 0L)
    .Internal(sys.call(which))

sys.calls <- function()
    .Internal(sys.calls())

sys.frame <- function(which = 0L)
    .Internal(sys.frame(which))

sys.function <- function(which = 0L)
    .Internal(sys.function(which))

sys.frames <- function()
    .Internal(sys.frames())

sys.nframe <- function()
    .Internal(sys.nframe())

sys.parent <- function(n = 1L)
    .Internal(sys.parent(n))

sys.parents <- function()
    .Internal(sys.parents())

sys.status <- function()
    list(sys.calls = sys.calls(), sys.parents = sys.parents(),
         sys.frames = sys.frames())

sys.on.exit <- function()
    .Internal(sys.on.exit())
#  File src/library/base/R/table.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

table <- function (..., exclude = if (useNA=="no") c(NA, NaN),
                   useNA = c("no", "ifany", "always"),
		   dnn = list.names(...), deparse.level = 1)
{
    list.names <- function(...) {
	l <- as.list(substitute(list(...)))[-1L]
	nm <- names(l)
	fixup <- if (is.null(nm)) seq_along(l) else nm == ""
	dep <- vapply(l[fixup], function(x)
		      switch(deparse.level + 1,
			     "", ## 0
			     if (is.symbol(x)) as.character(x) else "", ## 1
			     deparse(x, nlines=1)[1L] ## 2
			     ),
		      "")
	if (is.null(nm))
	    dep
	else {
	    nm[fixup] <- dep
	    nm
	}
    }
    if (!missing(exclude) && is.null(exclude))
        useNA <- "always"

    useNA <- match.arg(useNA)
    args <- list(...)
    if (!length(args))
	stop("nothing to tabulate")
    if (length(args) == 1L && is.list(args[[1L]])) {
	args <- args[[1L]]
	if (length(dnn) != length(args))
	    dnn <- if (!is.null(argn <- names(args)))
		 argn
	    else
		 paste(dnn[1L], seq_along(args), sep = ".")
    }
    # 0L, 1L, etc: keep 'bin' and 'pd' integer - as long as tabulate() requires it
    bin <- 0L
    lens <- NULL
    dims <- integer()
    pd <- 1L
    dn <- NULL
    for (a in args) {
	if (is.null(lens)) lens <- length(a)
	else if (length(a) != lens)
	    stop("all arguments must have the same length")
        cat <-
            if (is.factor(a)) {
                if (any(is.na(levels(a)))) # Don't touch this!
                    a
                else {
                    ## The logic here is tricky because it tries to do
                    ## something sensible if both 'exclude' and
                    ## 'useNA' are set.
                    ##
                    ## A non-null setting of 'exclude' sets the
                    ## excluded levels to missing, which is different
                    ## from the <NA> factor level. Excluded levels are
                    ## NOT tabulated, even if 'useNA' is set.
                    if (is.null(exclude) && useNA != "no")
                        addNA(a, ifany = (useNA == "ifany"))
                    else {
                        if (useNA != "no")
                            a <- addNA(a, ifany = (useNA == "ifany"))
                        ll <- levels(a)
                        a <- factor(a, levels = ll[!(ll %in% exclude)],
                               exclude = if (useNA == "no") NA)
                    }
                }
            }
            else { # NB: this excludes first, unlike the case above.
                a <- factor(a, exclude = exclude)
                if (useNA != "no")
                    addNA(a, ifany = (useNA == "ifany"))
                else
                    a
            }

	nl <- length(ll <- levels(cat))
	dims <- c(dims, nl)
        if (prod(dims) > .Machine$integer.max)
            stop("attempt to make a table with >= 2^31 elements")
	dn <- c(dn, list(ll))
	## requiring   all(unique(as.integer(cat)) == 1L:nlevels(cat))  :
	bin <- bin + pd * (as.integer(cat) - 1L)
	pd <- pd * nl
    }
    names(dn) <- dnn
    bin <- bin[!is.na(bin)]
    if (length(bin)) bin <- bin + 1L # otherwise, that makes bin NA
    y <- array(tabulate(bin, pd), dims, dimnames = dn)
    class(y) <- "table"
    y
}


## NB: NA in dimnames should be printed.
print.table <-
function (x, digits = getOption("digits"), quote = FALSE, na.print = "",
	  zero.print = "0",
	  justify = "none", ...)
{
    ## tables with empty extents have no contents and are hard to
    ## output in a readable way, so just say something descriptive and
    ## return.
    d <- dim(x)
    if (any(d == 0)) {
        cat ("< table of extent", paste(d, collapse=" x "), ">\n")
        return ( invisible(x) )
    }

    xx <- format(unclass(x), digits = digits, justify = justify)
    ## na.print handled here
    if(any(ina <- is.na(x)))
	xx[ina] <- na.print

    if(zero.print != "0" && any(i0 <- !ina & x == 0) && all(x == round(x)))
	## MM thinks this should be an option for many more print methods...
	xx[i0] <- sub("0", zero.print, xx[i0])

    ## Numbers get right-justified by format(), irrespective of 'justify'.
    ## We need to keep column headers aligned.
    if (is.numeric(x) || is.complex(x))
        print(xx, quote = quote, right = TRUE, ...)
    else
        print(xx, quote = quote, ...)
    invisible(x)
}

summary.table <- function(object, ...)
{
    if(!inherits(object, "table"))
	stop(gettextf("'object' must inherit from class %s",
                      dQuote("table")),
             domain = NA)
    n.cases <- sum(object)
    n.vars <- length(dim(object))
    y <- list(n.vars = n.vars,
	      n.cases = n.cases)
    if(n.vars > 1) {
	m <- vector("list", length = n.vars)
	relFreqs <- object / n.cases
	for(k in 1L:n.vars)
	    m[[k]] <- apply(relFreqs, k, sum)
	expected <- apply(do.call("expand.grid", m), 1L, prod) * n.cases
	statistic <- sum((c(object) - expected)^2 / expected)
	lm <- vapply(m, length, 1L)
	parameter <- prod(lm) - 1L - sum(lm - 1L)
	y <- c(y, list(statistic = statistic,
		       parameter = parameter,
		       approx.ok = all(expected >= 5),
		       p.value = stats::pchisq(statistic, parameter, lower.tail=FALSE),
		       call = attr(object, "call")))
    }
    class(y) <- "summary.table"
    y
}

print.summary.table <-
function(x, digits = max(1L, getOption("digits") - 3L), ...)
{
    if(!inherits(x, "summary.table"))
	stop(gettextf("'x' must inherit from class %s",
                      dQuote("summary.table")),
             domain = NA)
    if(!is.null(x$call)) {
	cat("Call: "); print(x$call)
    }
    cat("Number of cases in table:", x$n.cases, "\n")
    cat("Number of factors:", x$n.vars, "\n")
    if(x$n.vars > 1) {
	cat("Test for independence of all factors:\n")
	ch <- x$statistic
	cat("\tChisq = ",	format(round(ch, max(0, digits - log10(ch)))),
	    ", df = ",		x$parameter,
	    ", p-value = ",	format.pval(x$p.value, digits, eps = 0),
	    "\n", sep = "")
	if(!x$approx.ok)
	    cat("\tChi-squared approximation may be incorrect\n")
    }
    invisible(x)
}

as.data.frame.table <-
    function(x, row.names = NULL, ..., responseName = "Freq",
             stringsAsFactors = TRUE)
{
    ex <- quote(data.frame(do.call("expand.grid",
				   c(dimnames(provideDimnames(x)),
				     KEEP.OUT.ATTRS = FALSE,
                                     stringsAsFactors = stringsAsFactors)),
                           Freq = c(x),
                           row.names = row.names))
    names(ex)[3L] <- responseName
    eval(ex)
}

is.table <- function(x) inherits(x, "table")
as.table <- function(x, ...) UseMethod("as.table")
as.table.default <- function(x, ...)
{
    if(is.table(x)) return(x)
    else if(is.array(x) || is.numeric(x)) {
	x <- as.array(x)
	if(any(dim(x) == 0L)) stop("cannot coerce to a table")
	structure(class = c("table", oldClass(x)), provideDimnames(x))
    } else stop("cannot coerce to a table")
}

prop.table <- function(x, margin = NULL)
{
    if(length(margin))
	sweep(x, margin, margin.table(x, margin), "/", check.margin=FALSE)
    else
	x / sum(x)
}

margin.table <- function(x, margin = NULL)
{
    if(!is.array(x)) stop("'x' is not an array")
    if (length(margin)) {
	z <- apply(x, margin, sum)
	dim(z) <- dim(x)[margin]
	dimnames(z) <- dimnames(x)[margin]
    }
    else return(sum(x))
    class(z) <- oldClass(x) # avoid adding "matrix"
    z
}
#  File src/library/base/R/tabulate.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

tabulate <- function(bin, nbins = max(1L, bin, na.rm = TRUE))
{
    if(!is.numeric(bin) && !is.factor(bin))
	stop("'bin' must be numeric or a factor")
    ## avoid a copy for factors, since as.integer strips attributes
    if (typeof(bin) != "integer") bin <- as.integer(bin)
    if (nbins > .Machine$integer.max)
        stop("attempt to make a table with >= 2^31 elements")
    nbins <- as.integer(nbins)
    if (is.na(nbins)) stop("invalid value of 'nbins'")
    .Internal(tabulate(bin, nbins))
}
#  File src/library/base/R/tapply.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

tapply <- function (X, INDEX, FUN = NULL, ..., simplify = TRUE)
{
    FUN <- if (!is.null(FUN)) match.fun(FUN)
    if (!is.list(INDEX)) INDEX <- list(INDEX)
    nI <- length(INDEX)
    if (!nI) stop("'INDEX' is of length zero")
    namelist <- vector("list", nI)
    names(namelist) <- names(INDEX)
    extent <- integer(nI)
    nx <- length(X)
    one <- 1L
    group <- rep.int(one, nx) #- to contain the splitting vector
    ngroup <- one
    for (i in seq_along(INDEX)) {
	index <- as.factor(INDEX[[i]])
	if (length(index) != nx)
	    stop("arguments must have same length")
	namelist[[i]] <- levels(index)#- all of them, yes !
	extent[i] <- nlevels(index)
	group <- group + ngroup * (as.integer(index) - one)
	ngroup <- ngroup * nlevels(index)
    }
    if (is.null(FUN)) return(group)
    ans <- lapply(X = split(X, group), FUN = FUN, ...)
    index <- as.integer(names(ans))
    if (simplify && all(unlist(lapply(ans, length)) == 1L)) {
	ansmat <- array(dim = extent, dimnames = namelist)
	ans <- unlist(ans, recursive = FALSE)
    } else {
	ansmat <- array(vector("list", prod(extent)),
			dim = extent, dimnames = namelist)
    }
    if(length(index)) {
        names(ans) <- NULL
        ansmat[index] <- ans
    }
    ansmat
}





#  File src/library/base/R/taskCallback.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

addTaskCallback <- function(f, data = NULL, name = character())
{
    if(!is.function(f))
        stop("handler must be a function")
    val <- .Call(.C_R_addTaskCallback, f, data, !missing(data),
                 as.character(name))
    val + 1L
}

removeTaskCallback <- function(id)
{
    if(!is.character(id))
        id <- as.integer(id)

    .Call(.C_R_removeTaskCallback, id)
}

getTaskCallbackNames <- function() .Call(.C_R_getTaskCallbackNames)


taskCallbackManager <-
  #
  #
  #
function(handlers = list(), registered = FALSE, verbose = FALSE)
{
    suspended <- FALSE
    .verbose <- verbose

    add <-
    #
    # this is used to register a callback.
    # It has the same call sequence and semantics
    # as addTaskCallback but provides an optional
    # name by which to identify the element.
    # This can be used to remove the value in the future.
    # The default name is the next available position in the
    # list.
    # The result is stored in the `handlers' list using the
    # name.
    #
    # The element in the list contains the function
    # in the `f' slot,  and optionally a data field
    # to store the `data' argument.
    #
    # This could arrange to register itself using
    # addTaskCallback() if the size of the handlers list
    # becomes 1.
        function(f, data = NULL, name = NULL, register = TRUE)
        {

      # generate default name if none supplied
            if(is.null(name))
                name <- as.character(length(handlers) + 1L)

      # Add to handlers, replacing any element with that name
      # if needed.
            handlers[[name]] <<- list(f = f)

      # If data was specified, add this to the new element
      # so that it will be included in the call for this function
            if(!missing(data))
                handlers[[name]][["data"]] <<- data

      # We could arrange to register the evaluate function
      # so that the handlers list would be active. However,
      # we would have to unregister it in the remove()
      # function when there were no handlers.
            if(!registered && register) {
                register()
            }

            name
        }

    remove <- function(which)
    {
        if(is.character(which)) {
            tmp <- seq_along(handlers)[!is.na(match(which, names(handlers)))]
            if(length(tmp))
                stop(gettextf("no such element '%s'", which), domain = NA)
            which <- tmp
        } else
        which <- as.integer(which)

        handlers <<- handlers[-which]

        return(TRUE)
    }


    evaluate <-
    #
    # This is the actual callback that is registered with the C-level
    # mechanism. It is invoked by R when a top-level task is completed.
    # It then calls each of the functions in the handlers list
    # passing these functions the arguments it received and any
    # user-level data for those functions registered in the call to
    # add() via the `data' argument.
    #
    # At the end of the evaluation, any function that returned FALSE
    # is discarded.
        function(expr, value, ok, visible)
        {
            if(suspended)
                return(TRUE)
            discard <- character()
            for(i in names(handlers)) {
                h <- handlers[[i]]
                if(length(h) > 1L) {
                    val <- h[["f"]](expr, value, ok, visible, i[["data"]])
                } else {
                    val <- h[["f"]](expr, value, ok, visible)
                }
                if(!val) {
                    discard <- c(discard, i)
                }
            }
            if(length(discard)) {
                if(.verbose)
                    cat(gettextf("Removing %s", paste(discard, collapse=", ")), "\n")
                idx <- is.na(match(names(handlers), discard))
                if(length(idx))
                    handlers <<- handlers[idx]
                else
                    handlers <<- list()
            }
            return(TRUE)
        }

    suspend <-
        function(status = TRUE) {
            suspended <<- status
        }

    register <-
        function(name = "R-taskCallbackManager", verbose = .verbose)
        {
            if(verbose)
                cat(gettext("Registering 'evaluate' as low-level callback\n"))
            id <- addTaskCallback(evaluate, name = name)
            registered <<- TRUE
            id
        }

    list(add = add,
         evaluate = evaluate,
         remove = remove,
         register = register,
         suspend = suspend,
         callbacks = function()
         handlers
         )
}
#  File src/library/base/R/temp.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

tempfile <- function(pattern = "file", tmpdir = tempdir(), fileext = "")
    .Internal(tempfile(pattern, tmpdir, fileext))

tempdir <- function() .Internal(tempdir())
#  File src/library/base/R/time.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

system.time <- function(expr, gcFirst = TRUE)
{
    ppt <- function(y) {
        if(!is.na(y[4L])) y[1L] <- y[1L] + y[4L]
        if(!is.na(y[5L])) y[2L] <- y[2L] + y[5L]
        y[1L:3L]
    }
    if(!exists("proc.time")) return(rep(NA_real_, 5L))
    if(gcFirst)  gc(FALSE)
    time <- proc.time()
    ## need on.exit after 'time' has been set:
    ## on some systems proc.time throws an error.
    on.exit(cat("Timing stopped at:", ppt(proc.time() - time), "\n"))
    expr # evaluated here because of lazy evaluation
    new.time <- proc.time()
    on.exit()
    structure(new.time - time, class="proc_time")
}
unix.time <- system.time

date <- function() .Internal(date())

summary.proc_time <-
function(object, ...)
{
    if(!is.na(object[4L])) 
        object[1L] <- object[1L] + object[4L]
    if(!is.na(object[5L])) 
        object[2L] <- object[2L] + object[5L]
    object <- object[1L : 3L]
    names(object) <-
        c(gettext("user"), gettext("system"), gettext("elapsed"))
    object
}

print.proc_time <-
function(x, ...)
{
    print(summary(x, ...))
    invisible(x)
}
#  File src/library/base/R/toString.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

#functions to convert their first argument to strings
toString <- function(x, ...) UseMethod("toString")

toString.default <- function(x, width = NULL, ...)
{
    string <- paste(x, collapse=", ")
    if( missing(width) || is.null(width) || width == 0) return(string)
    if( width < 0 ) stop("'width' must be positive")
    if(nchar(string, type = "w") > width) {
        width <- max(6, width) ## Leave something!
        string <- paste0(strtrim(string, width - 4), "....")
    }
    string
}
#  File src/library/base/R/traceback.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

traceback <-
function(x = NULL, max.lines = getOption("deparse.max.lines"))
{
    if(is.null(x) && (exists(".Traceback", envir = baseenv())))
	x <- get(".Traceback", envir = baseenv())
    else if (is.numeric(x))
    	x <- .Internal(traceback(x))
    n <- length(x)
    if(n == 0L)
        cat(gettext("No traceback available"), "\n")
    else {
        for(i in 1L:n) {
            label <- paste0(n-i+1L, ": ")
            m <- length(x[[i]])
            if (!is.null(srcref <- attr(x[[i]], "srcref"))) {
            	srcfile <- attr(srcref, "srcfile")
            	x[[i]][m] <- paste0(x[[i]][m], " at ",
				    basename(srcfile$filename), "#", srcref[1L])
            }
            if(m > 1)
                label <- c(label, rep(substr("          ", 1L,
                                             nchar(label, type="w")),
                                      m - 1L))
            if(is.numeric(max.lines) && max.lines > 0L && max.lines < m) {
                cat(paste0(label[1L:max.lines], x[[i]][1L:max.lines]),
                    sep = "\n")
                cat(label[max.lines+1L], " ...\n")
            } else
            	cat(paste0(label, x[[i]]), sep="\n")
        }
    }
    invisible(x)
}
#  File src/library/base/R/unix/system.unix.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

system <- function(command, intern = FALSE,
                   ignore.stdout = FALSE, ignore.stderr = FALSE,
                   wait = TRUE, input = NULL,
                   show.output.on.console = TRUE, minimized = FALSE,
                   invisible = TRUE)
{
    if(!missing(show.output.on.console) || !missing(minimized)
       || !missing(invisible))
        message("arguments 'show.output.on.console', 'minimized' and 'invisible' are for Windows only")

    if(!is.logical(intern) || is.na(intern))
        stop("'intern' must be TRUE or FALSE")
    if(!is.logical(ignore.stdout) || is.na(ignore.stdout))
        stop("'ignore.stdout' must be TRUE or FALSE")
    if(!is.logical(ignore.stderr) || is.na(ignore.stderr))
        stop("'ignore.stderr' must be TRUE or FALSE")
    if(!is.logical(wait) || is.na(wait))
        stop("'wait' must be TRUE or FALSE")

    if(ignore.stdout) command <- paste(command, ">/dev/null")
    if(ignore.stderr) command <- paste(command, "2>/dev/null")
    if(!is.null(input)) {
        if(!is.character(input))
            stop("'input' must be a character vector or 'NULL'")
        f <- tempfile()
        on.exit(unlink(f))
        writeLines(input, f)
        # cat(input, file=f, sep="\n")
        command <- paste(command, "<", shQuote(f))
    }
    if(!wait && !intern) command <- paste(command, "&")
    .Internal(system(command, intern))
}

system2 <- function(command, args = character(),
                    stdout = "", stderr = "", stdin = "", input = NULL,
                    env = character(),
                    wait = TRUE, minimized = FALSE, invisible = TRUE)
{
    if(!missing(minimized) || !missing(invisible))
        message("arguments 'minimized' and 'invisible' are for Windows only")
    if(!is.logical(wait) || is.na(wait))
        stop("'wait' must be TRUE or FALSE")

    intern <- FALSE
    command <- paste(c(env, shQuote(command), args), collapse = " ")

    if(is.null(stdout)) stdout <- FALSE
    if(is.null(stderr)) stderr <- FALSE

    if (isTRUE(stderr)) {
        if (!isTRUE(stdout)) warning("setting stdout = TRUE")
        stdout <- TRUE
    }
    if (identical(stdout, FALSE))
        command <- paste(command, ">/dev/null")
    else if(isTRUE(stdout))
        intern <- TRUE
    else if(is.character(stdout)) {
        if(length(stdout) != 1L) stop("'stdout' must be of length 1")
        if(nzchar(stdout)) {
            command <- if (identical(stdout, stderr))
                paste(command, ">", shQuote(stdout), "2>&1")
            else command <- paste(command, ">", shQuote(stdout))
        }
    }
    if (identical(stderr, FALSE))
        command <- paste(command, "2>/dev/null")
    else if(isTRUE(stderr)) { # stdout == TRUE
        command <- paste(command, "2>&1")
    } else if(is.character(stderr)) {
        if(length(stderr) != 1L) stop("'stderr' must be of length 1")
        if(nzchar(stderr) && !identical(stdout, stderr))
            command <- paste(command, "2>", shQuote(stderr))
    }
    if(!is.null(input)) {
        if(!is.character(input))
            stop("'input' must be a character vector or 'NULL'")
        f <- tempfile()
        on.exit(unlink(f))
        writeLines(input, f)
        # cat(input, file=f, sep="\n")
        command <- paste(command, "<", shQuote(f))
    } else if (nzchar(stdin)) command <- paste(command, "<", stdin)
    if(!wait && !intern) command <- paste(command, "&")
    .Internal(system(command, intern))
}

## Some people try to use this with NA inputs (PR#15147)
Sys.which <- function(names)
{
    res <- character(length(names)); names(res) <- names
    ## hopefully configure found [/usr]/bin/which
    which <- "/usr/bin/which"
    if (!nzchar(which)) {
        warning("'which' was not found on this platform")
        return(res)
    }
    for(i in seq_along(names)) {
        if(is.na(names[i])) {res[i] <- NA; next}
        ## Quoting was added in 3.0.0
        ans <- suppressWarnings(system(paste(which, shQuote(names[i])),
                                       intern = TRUE, ignore.stderr = TRUE))
        ## Solaris' which gives 'no foo in ...' message on stdout,
        ## GNU which does it on stderr
        if(grepl("solaris", R.version$os)) {
            tmp <- strsplit(ans[1], " ", fixed = TRUE)[[1]]
            if(identical(tmp[1:3], c("no", i, "in"))) ans <- ""
        }
        res[i] <- if(length(ans)) ans[1] else ""
        ## final check that this is a real path and not an error message
        if(!file.exists(res[i])) res[i] <- ""
    }
    res
}
#  File src/library/base/R/unlist.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

unlist <- function(x, recursive=TRUE, use.names=TRUE)
{
    if(.Internal(islistfactor(x, recursive))) {
        lv <- unique(.Internal(unlist(lapply(x, levels), recursive, FALSE)))
        nm <- if(use.names) names(.Internal(unlist(x, recursive, use.names)))
        res <- .Internal(unlist(lapply(x, as.character), recursive, FALSE))
        res <- match(res, lv)
        ## we cannot make this ordered as level set may have been changed
        structure(res, levels=lv, names=nm, class="factor")
    } else .Internal(unlist(x, recursive, use.names))
}
#  File src/library/base/R/unname.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

unname <- function (obj, force = FALSE)
{
    if (!is.null(names(obj)))
        names(obj) <- NULL
    if (!is.null(dimnames(obj)) && (force || !is.data.frame(obj)))
        dimnames(obj) <- NULL
    obj
}
#  File src/library/base/R/upper.tri.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

upper.tri <- function(x, diag = FALSE)
{
    x <- as.matrix(x)
    if(diag) row(x) <= col(x)
    else row(x) < col(x)
}
#  File src/library/base/R/userhooks.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2013 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

## presumed small
.userHooksEnv <- new.env(hash = FALSE, parent = baseenv())

packageEvent <-
    function(pkgname, event=c("onLoad", "attach", "detach", "onUnload"))
{
    event <- match.arg(event)
    pkgname <- strsplit(pkgname, "_", fixed=TRUE)[[1L]][1L]
    paste("UserHook", pkgname, event, sep = "::")
}

getHook <- function(hookName)
{
    if (exists(hookName, envir = .userHooksEnv, inherits = FALSE))
        get(hookName, envir = .userHooksEnv, inherits = FALSE)
    else list()
}

setHook <- function(hookName, value,
                    action = c("append", "prepend", "replace"))
{
    action <- match.arg(action)
    old <- getHook(hookName)
    new <- switch(action,
                  "append" = c(old, value),
                  "prepend" = c(value, old),
                  "replace" = if (is.null(value) || is.list(value)) value else list(value))
    if (length(new))
        assign(hookName, new, envir = .userHooksEnv, inherits = FALSE)
    else if(exists(hookName, envir = .userHooksEnv, inherits = FALSE))
        remove(list=hookName, envir = .userHooksEnv, inherits = FALSE)
    invisible()
}
#  File src/library/base/R/utilities.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

mat.or.vec <- function(nr,nc)
    if(nc == 1L) numeric(nr) else matrix(0, nr, nc)

## Use  'version' since that exists in all S dialects :
is.R <-
    function() exists("version") && !is.null(vl <- version$language) && vl == "R"

#  File src/library/base/R/utils.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

shQuote <- function(string, type = c("sh", "csh", "cmd"))
{
    cshquote <- function(x) {
        xx <- strsplit(x, "'", fixed = TRUE)[[1L]]
        paste(paste0("'", xx, "'"), collapse="\"'\"")
    }
    if(missing(type) && .Platform$OS.type == "windows") type <- "cmd"
    type <- match.arg(type)
    if(type == "cmd") {
        paste0('"', gsub('"', '\\\\"', string), '"')
    } else {
        if(!length(string)) return("")
        has_single_quote <- grep("'", string)
        if(!length(has_single_quote))
            return(paste0("'", string, "'"))
        if(type == "sh")
            paste0('"', gsub('(["$`\\])', "\\\\\\1", string), '"')
        else {
            if(!length(grep("([$`])", string))) {
                paste0('"', gsub('(["!\\])', "\\\\\\1", string), '"')
            } else vapply(string, cshquote, "")
        }
    }
}

.standard_regexps <-
function()
{
    list(valid_package_name = "[[:alpha:]][[:alnum:].]*[[:alnum:]]",
         valid_package_version = "([[:digit:]]+[.-]){1,}[[:digit:]]+",
         valid_R_system_version =
         "[[:digit:]]+\\.[[:digit:]]+\\.[[:digit:]]+",
         valid_numeric_version = "([[:digit:]]+[.-])*[[:digit:]]+")
}
#  File src/library/base/R/vector.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2013 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

vector <- function(mode = "logical", length = 0L) .Internal(vector(mode, length))
logical <- function(length = 0L) .Internal(vector("logical", length))
character <- function(length = 0L) .Internal(vector("character", length))
integer <- function(length = 0L) .Internal(vector("integer", length))
numeric <- double <-
    function(length = 0L) .Internal(vector("double", length))

complex <- function(length.out = 0L,
		    real = numeric(), imaginary = numeric(),
		    modulus = 1, argument = 0) {
    if(missing(modulus) && missing(argument)) {
	## assume 'real' and 'imaginary'
	.Internal(complex(length.out, real, imaginary))
    } else {
	n <- max(length.out, length(argument), length(modulus))
	rep_len(modulus, n) * exp(1i * rep_len(argument, n))
    }
}

single <- function(length = 0L)
    structure(vector("double", length), Csingle=TRUE)
#  File src/library/base/R/version.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

## A simple S3 class for numeric versions (including package versions),
## and associated methods.

## We represent "vectors" of numeric versions as lists of sequences of
## integers, as obtained by splitting the version strings on the
## separators.  By default, only valid version specs (sequences of
## integers of suitable length), separated by '.' or '-', are allowed.
## If strictness is turned off, invalid specs result in integer()
## (rather than NA) to keep things simple.  (Note: using NULL would make
## subscripting more cumbersome ...)

## (In fact, the underlying mechanism could easily be extended to more
## general alphanumberic version specs.  E.g., one could allow "letters"
## in version numbers by replacing the non-sep characters in the version
## string by their ASCII codes.  However, this is not straightforward:
## alternatively, one could use an extended scheme with special markup
## for alpha, beta, release candidate, release, and patch versions, as
## used by many open source programs.  See e.g. the version::AlphaBeta
## module on CPAN.)

.make_numeric_version <-
function(x, strict = TRUE, regexp, classes = NULL)
{
    ## Internal creator for numeric version objects.

    nms <- names(x)
    x <- as.character(x)
    y <- rep.int(list(integer()), length(x))
    valid_numeric_version_regexp <- sprintf("^%s$", regexp)
    if(length(x)) {
        ok <- grepl(valid_numeric_version_regexp, x)
        if(!all(ok) && strict)
            stop(gettextf("invalid version specification %s",
                          paste(sQuote(unique(x[!ok])), collapse = ", ")),
                 call. = FALSE, domain = NA)
        y[ok] <- lapply(strsplit(x[ok], "[.-]"), as.integer)
    }
    names(y) <- nms
    class(y) <- unique(c(classes, "numeric_version"))
    y
}

## Basic numeric versions.

numeric_version <-
function(x, strict = TRUE)
    .make_numeric_version(x, strict,
                          .standard_regexps()$valid_numeric_version)

is.numeric_version <-
function(x)
    inherits(x, "numeric_version")

as.numeric_version <-
function(x)
{
    if(is.numeric_version(x)) x
    else if(is.package_version(x)) {
        ## Pre 2.6.0 is.package_version() compatibility code ...
        ## Simplify eventually ...
        structure(x, class = c(class(x), "numeric_version"))
    }
    else numeric_version(x)
}

## Package versions must have at least two integers, corresponding to
## major and minor.

package_version <-
function(x, strict = TRUE)
{
    ## Special-case R version lists.
    ## Currently, do this here for backward compatibility.
    ## Should this be changed eventually?
    if(is.list(x) && all(c("major", "minor") %in% names(x)))
        return(R_system_version(paste(x[c("major", "minor")],
                                      collapse = ".")))
    .make_numeric_version(x, strict,
                          .standard_regexps()$valid_package_version,
                          "package_version")
}

is.package_version <- function(x) inherits(x, "package_version")

as.package_version <- function(x)
    if(is.package_version(x)) x else package_version(x)

## R system versions must have exactly three integers.
## (Not sure if reduced strictness makes a lot of sense here.)

R_system_version <-
function(x, strict = TRUE)
    .make_numeric_version(x, strict,
                          .standard_regexps()$valid_R_system_version,
                          c("R_system_version", "package_version"))

getRversion <-
function()
    package_version(R.version)

## Workhorses.

## <NOTE>
## Could use this for or as as.double.numeric_version() ...
## </NOTE>

.encode_numeric_version <-
function(x, base = NULL)
{
    if(!is.numeric_version(x)) stop("wrong class")
    if(is.null(base)) base <- max(unlist(x), 0, na.rm = TRUE) + 1
    classes <- class(x)
    nms <- names(x)
    x <- unclass(x)
    lens <- vapply(x, length, 1L)
    ## We store the lengths so that we know when to stop when decoding.
    ## Alternatively, we need to be smart about trailing zeroes.  One
    ## approach is to increment all numbers in the version specs and
    ## base by 1, and when decoding only retain the non-zero entries and
    ## decrement by 1 one again.
    x <- vapply(x, function(t)
		sum(t / base^seq.int(0, length.out = length(t))), 1.)
    structure(ifelse(lens > 0L, x, NA_real_),
              base = base, lens = lens, .classes = classes, names = nms)
}

## <NOTE>
## Currently unused.
## Is there any point in having a 'base' argument?
## </NOTE>
.decode_numeric_version <-
function(x, base = NULL)
{
    if(is.null(base)) base <- attr(x, "base")
    if(!is.numeric(base)) stop("wrong argument")
    lens <- attr(x, "lens")
    y <- vector("list", length = length(x))
    for(i in seq_along(x)) {
        n <- lens[i]
        encoded <- x[i]
        decoded <- integer(n)
        for(k in seq_len(n)) {
            decoded[k] <- encoded %/% 1
            encoded <- base * (encoded %% 1)
        }
        y[[i]] <- as.integer(decoded)
    }
    class(y) <- unique(c(attr(x, ".classes"), "numeric_version"))
    y
}

## Methods.

`[.numeric_version` <-
function(x, i, j)
{
    y <- if(missing(j))
        unclass(x)[i]
    else
        lapply(unclass(x)[i], "[", j)
    ## Change sequences which are NULL or contains NAs to integer().
    bad <- vapply(y, function(t) is.null(t) || any(is.na(t)), NA)
    if(any(bad))
        y[bad] <- rep.int(list(integer()), length(bad))
    class(y) <- class(x)
    y
}

`[[.numeric_version` <-
function(x, ..., exact = NA)
{
   if(length(list(...)) < 2L)
      structure(list(unclass(x)[[..., exact=exact]]), class = oldClass(x))
   else
      unclass(x)[[..1, exact=exact]][..2]
}

## allowed forms
## x[[i]] <- "1.2.3"; x[[i]] <- 1L:3L; x[[c(i,j)]] <- <single integer>
## x[[i,j]] <- <single integer>
`[[<-.numeric_version` <-
function(x, ..., value)
{
   z <- unclass(x)
   if(nargs() < 4L) {
       if(length(..1) < 2L) {
           if(is.character(value) && length(value) == 1L)
               value <- unclass(as.numeric_version(value))[[1L]]
           else if(!is.integer(value)) stop("invalid 'value'")
       } else {
           value <- as.integer(value)
           if(length(value) != 1L) stop("invalid 'value'")
       }
       z[[..1]] <- value
   } else {
       value <- as.integer(value)
       if(length(value) != 1L) stop("invalid 'value'")
       z[[..1]][..2] <- value
   }
   structure(z, class = oldClass(x))
}


Ops.numeric_version <-
function(e1, e2)
{
    if(nargs() == 1L)
        stop(gettextf("unary '%s' not defined for \"numeric_version\" objects",
                      .Generic), domain = NA)
    boolean <- switch(.Generic, "<" = , ">" = , "==" = , "!=" = ,
        "<=" = , ">=" = TRUE, FALSE)
    if(!boolean)
        stop(gettextf("'%s' not defined for \"numeric_version\" objects",
                      .Generic), domain = NA)
    if(!is.numeric_version(e1)) e1 <- as.numeric_version(e1)
    if(!is.numeric_version(e2)) e2 <- as.numeric_version(e2)
    base <- max(unlist(e1), unlist(e2), 0) + 1
    e1 <- .encode_numeric_version(e1, base = base)
    e2 <- .encode_numeric_version(e2, base = base)
    NextMethod(.Generic)
}

Summary.numeric_version <-
function(..., na.rm)
{
    ok <- switch(.Generic, max = , min = , range = TRUE, FALSE)
    if(!ok)
        stop(gettextf("%s not defined for \"numeric_version\" objects",
                      .Generic), domain = NA)
    x <- do.call("c", lapply(list(...), as.numeric_version))
    v <- .encode_numeric_version(x)
    if(!na.rm && length(pos <- which(is.na(v)))) {
        y <- x[pos[1L]]
        if(as.character(.Generic) == "range")
            c(y, y)
        else
            y
    }
    else
        switch(.Generic,
               max = x[which.max(v)],
               min = x[which.min(v)],
               range = x[c(which.min(v), which.max(v))])
}

as.character.numeric_version <-
function(x, ...)
    as.character(format(x))

as.data.frame.numeric_version <- as.data.frame.vector

as.list.numeric_version <-
function(x, ...)
{
    nms <- names(x)
    names(x) <- NULL
    y <- lapply(seq_along(x), function(i) x[i])
    names(y) <- nms
    y
}

c.numeric_version <-
function(..., recursive = FALSE)
{
    x <- lapply(list(...), as.numeric_version)
    ## Try to preserve common extension classes.
    ## Note that this does not attempt to turn character strings into
    ## *package* versions if possible.
    classes <- if(length(unique(lapply(x, class))) == 1L)
        class(x[[1L]])
    else
        "numeric_version"
    structure(unlist(x, recursive = FALSE), class = classes)
}

duplicated.numeric_version <-
function(x, incomparables = FALSE, ...)
{
    x <- .encode_numeric_version(x)
    NextMethod("duplicated")
}

format.numeric_version <-
function(x, ...)
{
    x <- unclass(x)
    y <- rep.int(NA_character_, length(x))
    names(y) <- names(x)
    ind <- vapply(x, length, 1L) > 0L
    y[ind] <- unlist(lapply(x[ind], paste, collapse = "."))
    y
}

is.na.numeric_version <-
function(x)
    is.na(.encode_numeric_version(x))

print.numeric_version <-
function(x, ...)
{
    y <- as.character(x)
    if(!length(y))
        writeLines(gettext("<0 elements>"))
    else
        print(noquote(ifelse(is.na(y), NA_character_, sQuote(y))), ...)
    invisible(x)
}

rep.numeric_version <-
function(x, ...)
    structure(NextMethod("rep"), class = oldClass(x))

unique.numeric_version <-
function(x, incomparables = FALSE, ...)
    x[!duplicated(x, incomparables, ...)]

xtfrm.numeric_version <-
function(x)
    .encode_numeric_version(x)

## <NOTE>
## Versions of R prior to 2.6.0 had only a package_version class.
## We now have package_version extend numeric_version.
## We only provide named subscripting for package versions.
## </NOTE>

`$.package_version` <-
function(x, name)
{
    name <- pmatch(name, c("major", "minor", "patchlevel"))
    x <- unclass(x)
    switch(name,
	   major = vapply(x, "[", 0L, 1L),
	   minor = vapply(x, "[", 0L, 2L),
	   patchlevel = vapply(x, "[", 0L, 3L))
}
#  File src/library/base/R/warnings.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

warnings <- function(...)
{
    if(!exists("last.warning", envir=baseenv())) return()
    last.warning <- get("last.warning", envir=baseenv())
    if(!(n <- length(last.warning))) return()
    structure(last.warning, dots=list(...), class="warnings")
}

`[.warnings` <- function(x, ...)
    structure(NextMethod("["), class="warnings")

print.warnings <- function(x, ...)
{
    if(n <- length(x)) {
        cat(ngettext(n, "Warning message:\n", "Warning messages:\n"))
        msgs <- names(x)
        for(i in seq_len(n)) {
            ind <- if(n == 1L) "" else paste0(i, ": ")
            out <- if(length(x[[i]])) {
                ## deparse can overshoot cutoff
                temp <- deparse(x[[i]], width.cutoff = 50L, nlines = 2L)
                ## Put on one line if narrow enough.
                sm <- strsplit(msgs[i], "\n")[[1L]]
                nl <- if(nchar(ind, "w") + nchar(temp[1L], "w") +
                         nchar(sm[1L], "w") <= 75L)
                    " " else "\n  "
                paste(ind, "In ",
                      temp[1L], if(length(temp) > 1L) " ...",
                      " :", nl, msgs[i], sep = "")
            } else paste0(ind, msgs[i])
            do.call("cat", c(list(out), attr(x, "dots"), fill=TRUE))
        }
    }
    invisible(x)
}
#  File src/library/base/R/which.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

which <- function(x, arr.ind = FALSE, useNames = TRUE)
{
    wh <- .Internal(which(x))
    if (arr.ind && !is.null(d <- dim(x)))
	arrayInd(wh, d, dimnames(x), useNames=useNames) else wh
}

arrayInd <- function(ind, .dim, .dimnames = NULL, useNames = FALSE) {
    ##-- return a matrix  length(ind) x rank == length(ind) x length(.dim)
    m <- length(ind)
    rank <- length(.dim)
    wh1 <- ind - 1L
    ind <- 1L + wh1 %% .dim[1L]
    ind <- matrix(ind, nrow = m, ncol = rank,
		  dimnames = if(useNames)
		  list(.dimnames[[1L]][ind],
		       if(rank == 2L) c("row", "col") # for matrices
		       else paste0("dim", seq_len(rank))))
    if(rank >= 2L) {
	denom <- 1L
	for (i in 2L:rank) {
	    denom <- denom * .dim[i-1L]
	    nextd1 <- wh1 %/% denom	# (next dim of elements) - 1
	    ind[,i] <- 1L + nextd1 %% .dim[i]
	}
    }
    storage.mode(ind) <- "integer"
    ind
}

which.min <- function(x) .Internal(which.min(x))
which.max <- function(x) .Internal(which.max(x))
#  File src/library/utils/R/withVisible.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

withVisible <- function(x) .Internal(withVisible(x))
#  File src/library/base/R/write.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

write <- function(x, file = "data", ncolumns = if(is.character(x)) 1 else 5,
                  append = FALSE, sep = " ")
    cat(x, file = file, sep = c(rep.int(sep, ncolumns-1), "\n"),
        append = append)
#  File src/library/base/R/xor.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

xor <- function(x, y) { (x | y) & !(x & y) }
#  File src/library/base/R/zapsmall.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

zapsmall <- function(x, digits = getOption("digits"))
{
    if (length(digits) == 0L)
        stop("invalid 'digits'")
    if (all(ina <- is.na(x)))
        return(x)
    mx <- max(abs(x[!ina]))
    round(x, digits = if(mx > 0) max(0L, digits - log10(mx)) else digits)
}
#  File src/library/base/R/zdatetime.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

## needs to run after paste()
.leap.seconds <- local({
    .leap.seconds <-
        c("1972-6-30", "1972-12-31", "1973-12-31", "1974-12-31",
          "1975-12-31", "1976-12-31", "1977-12-31", "1978-12-31",
          "1979-12-31", "1981-6-30", "1982-6-30", "1983-6-30",
          "1985-6-30", "1987-12-31", "1989-12-31", "1990-12-31",
          "1992-6-30", "1993-6-30", "1994-6-30","1995-12-31",
          "1997-6-30", "1998-12-31", "2005-12-31", "2008-12-31", "2012-6-30")
    .leap.seconds <- strptime(paste(.leap.seconds , "23:59:60"),
                              "%Y-%m-%d %H:%M:%S")
    c(as.POSIXct(.leap.seconds, "GMT")) # lose the timezone
})
#  File src/library/base/R/zdynvars.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

## Need to ensure this comes late enough ...
## Perhaps even merge it into the common profile?

.dynLibs <- local({
    ## <NOTE>
    ## Versions of R prior to 1.4.0 had .Dyn.libs in .AutoloadEnv
    ## (and did not always ensure getting it from there).
    ## Until 1.6.0, we consistently used the base environment.
    ## Now we have a dynamic variable instead.
    ## </NOTE>
    .Dyn.libs <- structure(list(), class = "DLLInfoList")
    function(new) {
        if(!missing(new)) {
            class(new) <- "DLLInfoList"
            .Dyn.libs <<- new
        }
        else
            .Dyn.libs
    }
})

.libPaths <- local({
    .lib.loc <- character()            # Profiles need to set this.
    function(new) {
        if(!missing(new)) {
            ## paths don't really need to be unique, but searching
            ## large library trees repeatedly would be inefficient.
            ## Use normalizePath for display: but also does path.expand
            new <- Sys.glob(path.expand(new))
            paths <- unique(normalizePath(c(new, .Library.site, .Library), '/'))
            .lib.loc <<- paths[file.info(paths)$isdir %in% TRUE]
        }
        else
            .lib.loc
    }
})
#  File src/library/base/R/zzz.R
#  Part of the R package, http://www.R-project.org
#
#  Copyright (C) 1995-2013 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

## top-level assignments that need to be copied to baseloader.R
as.numeric <- as.double
is.name <- is.symbol


## extracted from existing NAMESPACE files in Dec 2003
.knownS3Generics <- local({

    ## include the S3 group generics here
    baseGenerics <- c("Math", "Ops", "Summary", "Complex",
        "as.character", "as.data.frame", "as.environment", "as.matrix", "as.vector",
        "cbind", "labels", "print", "rbind", "rep", "seq", "seq.int",
        "solve", "summary", "t")

    utilsGenerics <- c("edit", "str")

    graphicsGenerics <- c("contour", "hist", "identify", "image",
        "lines", "pairs", "plot", "points", "text")

    statsGenerics <- c( "add1", "AIC", "anova", "biplot", "coef",
        "confint", "deviance", "df.residual", "drop1", "extractAIC",
        "fitted", "formula", "logLik", "model.frame", "model.matrix",
        "predict", "profile", "qqnorm", "residuals", "se.contrast",
        "terms", "update", "vcov")

    tmp <- rep.int(c("base", "utils", "graphics", "stats"),
                   c(length(baseGenerics), length(utilsGenerics),
                     length(graphicsGenerics), length(statsGenerics)))
    names(tmp) <-
        c(baseGenerics, utilsGenerics, graphicsGenerics, statsGenerics)
    tmp
})

.ArgsEnv <- new.env(hash = TRUE, parent = emptyenv())

assign("%*%", function(x, y) NULL, envir = .ArgsEnv)
assign(".C", function(.NAME, ..., NAOK = FALSE, DUP = TRUE, PACKAGE,
                      ENCODING) NULL,
       envir = .ArgsEnv)
assign(".Fortran",
       function(.NAME, ..., NAOK = FALSE, DUP = TRUE, PACKAGE, ENCODING) NULL,
       envir = .ArgsEnv)
assign(".Call", function(.NAME, ..., PACKAGE) NULL, envir = .ArgsEnv)
assign(".Call.graphics", function(.NAME, ..., PACKAGE) NULL, envir = .ArgsEnv)
assign(".External", function(.NAME, ..., PACKAGE) NULL, envir = .ArgsEnv)
assign(".External2", function(.NAME, ..., PACKAGE) NULL, envir = .ArgsEnv)
assign(".External.graphics", function(.NAME, ..., PACKAGE) NULL,
       envir = .ArgsEnv)
assign(".Internal", function(call) NULL, envir = .ArgsEnv)
assign(".Primitive", function(name) NULL, envir = .ArgsEnv)
assign(".isMethodsDispatchOn", function(x, onOff = NULL) NULL, envir = .ArgsEnv)
assign(".primTrace", function(obj) NULL, envir = .ArgsEnv)
assign(".primUntrace", function(obj) NULL, envir = .ArgsEnv)
assign(".subset", function(x, ...) NULL, envir = .ArgsEnv)
assign(".subset2", function(x, ...) NULL, envir = .ArgsEnv)
assign("UseMethod", function(generic, object) NULL, envir = .ArgsEnv)
assign("as.call", function(x) NULL, envir = .ArgsEnv)
assign("attr", function(x, which, exact = FALSE) NULL, envir = .ArgsEnv)
assign("attr<-", function(x, which, value) NULL, envir = .ArgsEnv)
assign("attributes", function(obj) NULL, envir = .ArgsEnv)
assign("attributes<-", function(obj, value) NULL, envir = .ArgsEnv)
assign("baseenv", function() NULL, envir = .ArgsEnv)
assign("browser",
       function(text="", condition=NULL, expr = TRUE, skipCalls = 0L) NULL,
       envir = .ArgsEnv)
assign("call", function(name, ...) NULL, envir = .ArgsEnv)
assign("class", function(x) NULL, envir = .ArgsEnv)
assign("class<-", function(x, value) NULL, envir = .ArgsEnv)
assign(".cache_class", function(class, extends) NULL, envir = .ArgsEnv)
assign("emptyenv", function() NULL, envir = .ArgsEnv)
assign("enc2native", function(x) NULL, envir = .ArgsEnv)
assign("enc2utf8", function(x) NULL, envir = .ArgsEnv)
assign("environment<-", function(fun, value) NULL, envir = .ArgsEnv)
assign("expression", function(...) NULL, envir = .ArgsEnv)
assign("gc.time", function(on = TRUE) NULL, envir = .ArgsEnv)
assign("globalenv", function() NULL, envir = .ArgsEnv)
assign("interactive", function() NULL, envir = .ArgsEnv)
assign("invisible", function(x) NULL, envir = .ArgsEnv)
assign("is.atomic", function(x) NULL, envir = .ArgsEnv)
assign("is.call", function(x) NULL, envir = .ArgsEnv)
assign("is.character", function(x) NULL, envir = .ArgsEnv)
assign("is.complex", function(x) NULL, envir = .ArgsEnv)
assign("is.double", function(x) NULL, envir = .ArgsEnv)
assign("is.environment", function(x) NULL, envir = .ArgsEnv)
assign("is.expression", function(x) NULL, envir = .ArgsEnv)
assign("is.function", function(x) NULL, envir = .ArgsEnv)
assign("is.integer", function(x) NULL, envir = .ArgsEnv)
assign("is.language", function(x) NULL, envir = .ArgsEnv)
assign("is.list", function(x) NULL, envir = .ArgsEnv)
assign("is.logical", function(x) NULL, envir = .ArgsEnv)
assign("is.name", function(x) NULL, envir = .ArgsEnv)
assign("is.null", function(x) NULL, envir = .ArgsEnv)
assign("is.object", function(x) NULL, envir = .ArgsEnv)
assign("is.pairlist", function(x) NULL, envir = .ArgsEnv)
assign("is.raw", function(x) NULL, envir = .ArgsEnv)
assign("is.recursive", function(x) NULL, envir = .ArgsEnv)
assign("is.single", function(x) NULL, envir = .ArgsEnv)
assign("is.symbol", function(x) NULL, envir = .ArgsEnv)
assign("isS4", function(object) NULL, envir = .ArgsEnv)
assign("list", function(...) NULL, envir = .ArgsEnv)
assign("lazyLoadDBfetch", function(key, file, compressed, hook) NULL,
       envir = .ArgsEnv)
assign("missing", function(x) NULL, envir = .ArgsEnv)
assign("nargs", function() NULL, envir = .ArgsEnv)
assign("nzchar", function(x) NULL, envir = .ArgsEnv)
assign("oldClass", function(x) NULL, envir = .ArgsEnv)
assign("oldClass<-", function(x, value) NULL, envir = .ArgsEnv)
assign("on.exit", function(expr = NULL, add = FALSE) NULL, envir = .ArgsEnv)
assign("pos.to.env", function(x) NULL, envir = .ArgsEnv)
assign("proc.time", function() NULL, envir = .ArgsEnv)
assign("quote", function(expr) NULL, envir = .ArgsEnv)
assign("retracemem", function(x, previous = NULL) NULL, envir = .ArgsEnv)
assign("seq_along", function(along.with) NULL, envir = .ArgsEnv)
assign("seq_len", function(length.out) NULL, envir = .ArgsEnv)
assign("standardGeneric", function(f, fdef) NULL, envir = .ArgsEnv)
assign("storage.mode<-", function(x, value) NULL, envir = .ArgsEnv)
assign("substitute", function(expr, env) NULL, envir = .ArgsEnv)
assign("switch", function(EXPR, ...) NULL, envir = .ArgsEnv)
assign("tracemem", function(x) NULL, envir = .ArgsEnv)
assign("unclass", function(x) NULL, envir = .ArgsEnv)
assign("untracemem", function(x) NULL, envir = .ArgsEnv)


.S3PrimitiveGenerics <-
  c("as.character", "as.complex", "as.double", "as.environment",
    "as.integer", "as.logical", "as.numeric", "as.raw",
    "c", "dim", "dim<-", "dimnames", "dimnames<-",
    "is.array", "is.finite",
    "is.infinite", "is.matrix", "is.na", "is.nan", "is.numeric",
    "length", "length<-", "levels<-", "names", "names<-", "rep",
    "seq.int", "xtfrm")

.GenericArgsEnv <- local({
    env <- new.env(hash = TRUE, parent = emptyenv())
    ## those with different arglists are overridden below.
    for(f in .S3PrimitiveGenerics) {
        fx <- function(x) {}
        body(fx) <- substitute(UseMethod(ff), list(ff=f))
        environment(fx) <- .BaseNamespaceEnv
        assign(f, fx, envir = env)
    }
    ## now add the group generics
    ## round, signif, log, trunc are handled below
    fx <- function(x) {}
    for(f in c("abs", "sign", "sqrt", "floor", "ceiling",
               "exp", "expm1", "log1p", "log10", "log2",
               "cos", "sin", "tan", "acos", "asin", "atan", "cosh", "sinh",
               "tanh", "acosh", "asinh", "atanh",
               "gamma", "lgamma", "digamma", "trigamma",
               "cumsum", "cumprod", "cummax", "cummin")) {
        body(fx) <- substitute(UseMethod(ff), list(ff=f))
        environment(fx) <- .BaseNamespaceEnv
        assign(f, fx, envir = env)
    }

    ## ! is unary and handled below
    fx <- function(e1, e2) {}
    for(f in c("+", "-", "*", "/", "^", "%%", "%/%", "&", "|",
               "==", "!=", "<", "<=", ">=", ">")) {
        body(fx) <- substitute(UseMethod(ff), list(ff=f))
        environment(fx) <- .BaseNamespaceEnv
        assign(f, fx, envir = env)
    }

    for(f in c("all", "any", "sum", "prod", "max", "min", "range")) {
        fx <- function(..., na.rm = FALSE) {}
        body(fx) <- substitute(UseMethod(ff), list(ff=f))
        environment(fx) <- .BaseNamespaceEnv
        assign(f, fx, envir = env)
    }

    for(f in c("Arg", "Conj", "Im", "Mod", "Re")) {
        fx <- function(z) {}
        body(fx) <- substitute(UseMethod(ff), list(ff=f))
        environment(fx) <- .BaseNamespaceEnv
        assign(f, fx, envir = env)
    }
    env
})
### do these outside to get the base namespace as the environment.
assign("!", function(x) UseMethod("!"), envir = .GenericArgsEnv)
assign("as.character", function(x, ...) UseMethod("as.character"),
       envir = .GenericArgsEnv)
assign("as.complex", function(x, ...) UseMethod("as.complex"),
       envir = .GenericArgsEnv)
assign("as.double", function(x, ...) UseMethod("as.double"),
       envir = .GenericArgsEnv)
assign("as.integer", function(x, ...) UseMethod("as.integer"),
       envir = .GenericArgsEnv)
assign("as.logical", function(x, ...) UseMethod("as.logical"),
       envir = .GenericArgsEnv)
#assign("as.raw", function(x) UseMethod("as.raw"), envir = .GenericArgsEnv)
assign("c", function(..., recursive = FALSE) UseMethod("c"),
       envir = .GenericArgsEnv)
#assign("dimnames", function(x) UseMethod("dimnames"), envir = .GenericArgsEnv)
assign("dim<-", function(x, value) UseMethod("dim<-"), envir = .GenericArgsEnv)
assign("dimnames<-", function(x, value) UseMethod("dimnames<-"),
       envir = .GenericArgsEnv)
assign("length<-", function(x, value) UseMethod("length<-"),
       envir = .GenericArgsEnv)
assign("levels<-", function(x, value) UseMethod("levels<-"),
       envir = .GenericArgsEnv)
assign("log", function(x, base=exp(1)) UseMethod("log"),
       envir = .GenericArgsEnv)
assign("names<-", function(x, value) UseMethod("names<-"),
       envir = .GenericArgsEnv)
assign("rep", function(x, ...) UseMethod("rep"), envir = .GenericArgsEnv)
assign("round", function(x, digits=0) UseMethod("round"),
       envir = .GenericArgsEnv)
assign("seq.int", function(from, to, by, length.out, along.with, ...)
       UseMethod("seq.int"), envir = .GenericArgsEnv)
assign("signif", function(x, digits=6) UseMethod("signif"),
       envir = .GenericArgsEnv)
assign("trunc", function(x, ...) UseMethod("trunc"), envir = .GenericArgsEnv)
#assign("xtfrm", function(x) UseMethod("xtfrm"), envir = .GenericArgsEnv)

## make this the same object as as.double
assign("as.numeric", get("as.double", envir = .GenericArgsEnv),
       envir = .GenericArgsEnv)
