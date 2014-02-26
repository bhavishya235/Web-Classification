#  File src/library/parallel/R/RngStream.R
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

nextRNGStream <- function(seed)
{
    if(!is.integer(seed) || seed[1L] %% 100L != 7L)
        stop("invalid value of 'seed'")
    .Call(C_nextStream, seed)
}

nextRNGSubStream <- function(seed)
{
    if(!is.integer(seed) || seed[1L] %% 100L != 7L)
        stop("invalid value of 'seed'")
    .Call(C_nextSubStream, seed)
}

## Different from snow's RNG code
clusterSetRNGStream <- function(cl = NULL, iseed = NULL)
{
    cl <- defaultCluster(cl)
    oldseed <-
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        else NULL
    RNGkind("L'Ecuyer-CMRG")
    if(!is.null(iseed)) set.seed(iseed)
    nc <- length(cl)
    seeds <- vector("list", nc)
    seeds[[1L]] <- .Random.seed
    for(i in seq_len(nc-1L)) seeds[[i+1L]] <- nextRNGStream(seeds[[i]])
    ## Reset the random seed in the master.
    if(!is.null(oldseed))
        assign(".Random.seed", oldseed, envir = .GlobalEnv)
    else rm(.Random.seed, envir = .GlobalEnv)
    for (i in seq_along(cl)) {
        expr <- substitute(assign(".Random.seed", seed, envir = .GlobalEnv),
                           list(seed = seeds[[i]]))
        sendCall(cl[[i]], eval, list(expr))
    }
    checkForRemoteErrors(lapply(cl, recvResult))
    invisible()
}

RNGenv <- new.env()

mc.reset.stream <- function() {
    if (RNGkind()[1L] == "L'Ecuyer-CMRG") {
        if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            runif(1)
        assign("LEcuyer.seed",
               get(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
               envir = RNGenv)
    }
}

## For use in the master before forking
mc.advance.stream <- function(reset = FALSE)
{
    if (RNGkind()[1L] == "L'Ecuyer-CMRG") {
        if (reset ||
            !exists("LEcuyer.seed", envir = RNGenv, inherits = FALSE)) {
            if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
                runif(1)
            assign("LEcuyer.seed",
                   get(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
                   envir = RNGenv)
        } else {
            assign("LEcuyer.seed",
                   nextRNGStream(get("LEcuyer.seed", envir = RNGenv)),
                   envir = RNGenv)
        }
    }
}

## For use in the child
mc.set.stream <- function()
{
    if (RNGkind()[1L] == "L'Ecuyer-CMRG") {
            assign(".Random.seed", get("LEcuyer.seed", envir = RNGenv),
                   envir = .GlobalEnv)
    } else {
        ## It is random to simply unset the seed
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
}
#  File src/library/parallel/R/clusterApply.R
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

## Derived from snow 0.3-6 by Luke Tierney

staticClusterApply <- function(cl = NULL, fun, n, argfun) {
    cl <- defaultCluster(cl)
    p <- length(cl)
    if (n > 0L && p) {
        val <- vector("list", n)
        start <- 1L
        while (start <= n) {
            end <- min(n, start + p - 1L)
	    jobs <- end - start + 1L
            for (i in 1:jobs)
                sendCall(cl[[i]], fun, argfun(start + i - 1L))
            val[start:end] <- lapply(cl[1:jobs], recvResult)
            start <- start + jobs
        }
        checkForRemoteErrors(val)
    }
}

dynamicClusterApply <- function(cl = NULL, fun, n, argfun) {
    cl <- defaultCluster(cl)
    p <- length(cl)
    if (n > 0L && p) {
        submit <- function(node, job)
            sendCall(cl[[node]], fun, argfun(job), tag = job)
        for (i in 1:min(n, p)) submit(i, i)
        val <- vector("list", n)
        for (i in 1:n) {
            d <- recvOneResult(cl)
            j <- i + min(n, p)
            if (j <= n) submit(d$node, j)
            val[d$tag] <- list(d$value)
        }
        checkForRemoteErrors(val)
    }
}

## exported and documented from here down unless otherwise stated.

clusterCall  <- function(cl = NULL, fun, ...)
{
    cl <- defaultCluster(cl)
    for (i in seq_along(cl)) sendCall(cl[[i]], fun, list(...))
    checkForRemoteErrors(lapply(cl, recvResult))
}


clusterEvalQ <- function(cl = NULL, expr)
    clusterCall(cl, eval, substitute(expr), env=.GlobalEnv)

clusterExport <- local({
    gets <- function(n, v) { assign(n, v, envir = .GlobalEnv); NULL }
    function(cl = NULL, varlist, envir = .GlobalEnv) {
        ## do this with only one clusterCall--loop on workers?
        for (name in varlist) {
            clusterCall(cl, gets, name, get(name, envir = envir))
        }
    }
})

clusterApply <- function(cl = NULL, x, fun, ...)
{
    ## **** this closure is sending all of x to all nodes
    argfun <- function(i) c(list(x[[i]]), list(...))
    staticClusterApply(cl, fun, length(x), argfun)
}

clusterApplyLB <- function(cl = NULL, x, fun, ...)
{
    ## **** this closure is sending all of x to all nodes
    argfun <- function(i) c(list(x[[i]]), list(...))
    dynamicClusterApply(cl, fun, length(x), argfun)
}

clusterMap <- function (cl = NULL, fun, ..., MoreArgs = NULL, RECYCLE = TRUE,
                        SIMPLIFY = FALSE, USE.NAMES = TRUE,
                        .scheduling = c("static", "dynamic"))
{
    cl <- defaultCluster(cl)
    args <- list(...)
    if (length(args) == 0) stop("need at least one argument")
    .scheduling <- match.arg(.scheduling)
    n <- sapply(args, length)
    if (RECYCLE) {
        vlen <- max(n)
        if(vlen && min(n) == 0L)
            stop("zero-length inputs cannot be mixed with those of non-zero length")
        if (!all(n == vlen))
            for (i in seq_along(args)) # why not lapply?
                args[[i]] <- rep(args[[i]], length.out = vlen)
    }
    else vlen <- min(n)
    ## **** this closure is sending all of ... to all nodes
    argfun <- function(i) c(lapply(args, function(x) x[[i]]), MoreArgs)
    answer <-
        if(.scheduling == "dynamic") dynamicClusterApply(cl, fun, vlen, argfun)
    else staticClusterApply(cl, fun, vlen, argfun)
    ## rest matches mapply(): with a different default for SIMPLIFY
    if (USE.NAMES && length(args)) {
        if (is.null(names1 <- names(args[[1L]])) && is.character(args[[1L]]))
            names(answer) <- args[[1L]]
        else if (!is.null(names1))
            names(answer) <- names1
    }
    if (!identical(SIMPLIFY, FALSE) && length(answer))
        simplify2array(answer, higher = (SIMPLIFY == "array"))
    else answer
}

## splitIndices <- function(nx, ncl)
## {
##     i <- seq_len(nx)
##     if (ncl == 1L) i
##     else structure(split(i, cut(i, ncl)), names = NULL)
## }

# The fuzz used by cut() is too small when nx and ncl are both large
# and causes some groups to be empty. The definition below avoids that
# while minimizing changes from the results produced by the definition
# above.
splitIndices <- function(nx, ncl) {
    i <- 1L:nx
    if (ncl == 1L || nx == 1L) i
    else {
        fuzz <- min((nx - 1L) / 1000, 0.4 * nx / ncl)
        breaks <- seq(1 - fuzz, nx + fuzz, length = ncl + 1L)
        structure(split(i, cut(i, breaks)), names = NULL)
    }
}

clusterSplit <- function(cl = NULL, seq) {
    cl <- defaultCluster(cl)
    lapply(splitIndices(length(seq), length(cl)), function(i) seq[i])
}

#internal
splitList <- function(x, ncl)
    lapply(splitIndices(length(x), ncl), function(i) x[i])

#internal
splitRows <- function(x, ncl)
    lapply(splitIndices(nrow(x), ncl), function(i) x[i, , drop=FALSE])

#internal
splitCols <- function(x, ncl)
    lapply(splitIndices(ncol(x), ncl), function(i) x[, i, drop=FALSE])

parLapply <- function(cl = NULL, X, fun, ...)
{
    cl <- defaultCluster(cl)
    do.call(c,
            clusterApply(cl, x = splitList(X, length(cl)),
                         fun = lapply, fun, ...),
            quote = TRUE)
}

parLapplyLB <- function(cl = NULL, X, fun, ...)
{
    cl <- defaultCluster(cl)
    do.call(c,
            clusterApplyLB(cl, x = splitList(X, length(cl)),
                           fun = lapply, fun, ...),
            quote = TRUE)
}

parRapply <- function(cl = NULL, x, FUN, ...)
{
    cl <- defaultCluster(cl)
    do.call(c,
            clusterApply(cl = cl, x = splitRows(x, length(cl)),
                         fun = apply, MARGIN = 1L, FUN = FUN, ...),
            quote = TRUE)
}

parCapply <- function(cl = NULL, x, FUN, ...) {
    cl <- defaultCluster(cl)
    do.call(c,
            clusterApply(cl = cl, x = splitCols(x, length(cl)),
                         fun = apply, MARGIN = 2L, FUN = FUN, ...),
            quote = TRUE)
}


parSapply <-
    function (cl = NULL, X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE)
{
    FUN <- match.fun(FUN) # should this be done on worker?
    answer <- parLapply(cl, X = as.list(X), fun = FUN, ...)
    if(USE.NAMES && is.character(X) && is.null(names(answer)))
	names(answer) <- X
    if(!identical(simplify, FALSE) && length(answer))
	simplify2array(answer, higher = (simplify == "array"))
    else answer
}

parSapplyLB <-
    function (cl = NULL, X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE)
{
    FUN <- match.fun(FUN) # should this be done on worker?
    answer <- parLapplyLB(cl, X = as.list(X), fun = FUN, ...)
    if(USE.NAMES && is.character(X) && is.null(names(answer)))
	names(answer) <- X
    if(!identical(simplify, FALSE) && length(answer))
	simplify2array(answer, higher = (simplify == "array"))
    else answer
}


parApply <- function(cl = NULL, X, MARGIN, FUN, ...)
{
    cl <- defaultCluster(cl) # initial sanity check
    FUN <- match.fun(FUN) # should this be done on worker?

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
    d.ans  <- d[MARGIN]
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
    arglist <- if(length(d.call) < 2L) {# vector
        if (length(dn.call)) dimnames(newX) <- c(dn.call, list(NULL))
        lapply(seq_len(d2), function(i) newX[,i])
    } else
        lapply(seq_len(d2), function(i) array(newX[,i], d.call, dn.call))
    ans <- parLapply(cl = cl, X = arglist, fun = FUN, ...)

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

#  File src/library/parallel/R/detectCores.R
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

## In part based on code in package multicore 0.1-6 by Simon Urbanek

detectCores <-
    if(.Platform$OS.type == "windows") {
        function(all.tests = FALSE, logical = TRUE) {
            ## result is # cores, logical processors.
            res <- .Call(C_ncpus, FALSE)
            ifelse(logical, res[2L], res[1L]);
        }
    } else {
        function(all.tests = FALSE, logical = FALSE) {
            systems <-
                list(darwin = "/usr/sbin/sysctl -n hw.ncpu 2>/dev/null",
                     freebsd = "/sbin/sysctl -n hw.ncpu 2>/dev/null",
                     linux = "grep processor /proc/cpuinfo 2>/dev/null | wc -l",
                     irix  = c("hinv | grep Processors | sed 's: .*::'",
                     "hinv | grep '^Processor '| wc -l"),
                     solaris = if(logical) "/usr/sbin/psrinfo -v | grep 'Status of.*processor' | wc -l" else "/bin/kstat -p -m cpu_info | grep :core_id | cut -f2 | uniq | wc -l")
            for (i in seq(systems))
                if(all.tests ||
		   length(grep(paste0("^", names(systems)[i]), R.version$os)))
                    for (cmd in systems[i]) {
                        a <- gsub("^ +","", system(cmd, TRUE)[1])
                        if (length(grep("^[1-9]", a))) return(as.integer(a))
                    }
            NA_integer_
        }
    }
#  File src/library/parallel/R/snow.R
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

## Derived from snow 0.3-6 by Luke Tierney

.reg <-  new.env()
assign("default", NULL, envir = .reg)

defaultCluster <- function(cl = NULL)
{
    if(is.null(cl)) cl <- get("default", envir = .reg)
    if(is.null(cl)) stop("no cluster 'cl' supplied and none is registered")
    checkCluster(cl)
    cl
}

setDefaultCluster <- function(cl = NULL)
{
    if(!is.null(cl)) checkCluster(cl)
    assign("default", cl, envir = .reg)
}

#
# Checking and subsetting
#

checkCluster <- function(cl)
    if (!inherits(cl, "cluster")) stop("not a valid cluster");

`[.cluster` <- function(cl, ...) {
    v <- NextMethod()
    class(v) <- class(cl)
    v
}


#
# Higher-Level Node Functions
#

closeNode <- function(node) UseMethod("closeNode")
closeNode.default <- function(node) {}

## These have SOCK methods
sendData <- function(node, data) UseMethod("sendData")
recvData <- function(node) UseMethod("recvData")
recvOneData <- function(cl) UseMethod("recvOneData")

postNode <- function(con, type, value = NULL, tag = NULL)
    sendData(con, list(type = type, data = value, tag = tag))

stopNode <- function(n) {
    postNode(n, "DONE")
    closeNode(n)
}



#
#  Cluster Creation and Destruction
#

defaultClusterOptions <- NULL

#**** check valid cluster option

initDefaultClusterOptions <- function(libname)
{
    rscript <- file.path(R.home("bin"), "Rscript")
    port <- Sys.getenv("R_PARALLEL_PORT")
    port <- if (identical(port, "random")) NA else as.integer(port)
    if (is.na(port))
        port <- 11000 + 1000 * ((stats::runif(1L) + unclass(Sys.time())/300) %% 1)
    options <- list(port = as.integer(port),
                    timeout = 60 * 60 * 24 * 30, # 30 days
                    master =  Sys.info()["nodename"],
                    homogeneous = TRUE,
                    type = "PSOCK",
                    outfile = "/dev/null",
                    rscript = rscript,
                    user = Sys.info()["user"],
                    rshcmd = "ssh",
                    manual = FALSE,
                    methods = TRUE,
                    renice = NA_integer_,
                    ## rest are unused in parallel
                    rhome = R.home(),
                    rlibs = Sys.getenv("R_LIBS"),
                    scriptdir = file.path(libname, "parallel"),
                    rprog = file.path(R.home("bin"), "R"),
                    snowlib = .libPaths()[1],
                    useRscript = TRUE, useXDR = TRUE)
    defaultClusterOptions <<- addClusterOptions(emptyenv(), options)
}

addClusterOptions <- function(options, new) {
    if (!is.null(new)) {
        options <- new.env(parent = options)
        names <- names(new)
        for (i in seq_along(new))
            assign(names[i], new[[i]], envir = options)
    }
    options
}

getClusterOption <- function(name, options = defaultClusterOptions)
    get(name, envir = options)

setDefaultClusterOptions <- function(...) {
    list <- list(...)
    names <- names(list)
    for (i in seq_along(list))
        assign(names[i], list[[i]], envir = defaultClusterOptions)
}


makeCluster <-
    function (spec, type = getClusterOption("type"), ...)
{
    switch(type,
           PSOCK = makePSOCKcluster(spec, ...),
           FORK = makeForkCluster(spec, ...),
           SOCK = snow::makeSOCKcluster(spec, ...),
           MPI = snow::makeMPIcluster(spec, ...),
           NWS = snow::makeNWScluster(spec, ...),
           stop("unknown cluster type"))
}


stopCluster <- function(cl = NULL)
{
    cl <- defaultCluster(cl)
    if(identical(cl, get("default", envir = .reg)))
        assign("default", NULL, envir = .reg)
    UseMethod("stopCluster")
}

stopCluster.default <- function(cl) for (n in cl) stopNode(n)


#
# Cluster Functions
#

sendCall <- function (con, fun, args, return = TRUE, tag = NULL)
{
    timing <-  .snowTimingData$running()
    postNode(con, "EXEC",
             list(fun = fun, args = args, return = return, tag = tag))
    if (timing)
        .snowTimingData$enterSend(con$rank, start, proc.time()[3L])
    NULL
}

recvResult <- function(con)
{
    if (.snowTimingData$running()) {
        start <- proc.time()[3L]
        r <- recvData(con)
        end <- proc.time()[3L]
        .snowTimingData$enterRecv(con$rank, start, end, r$time[3L])
    }
    else r <- recvData(con)
    r$value
}

checkForRemoteErrors <- function(val)
{
    count <- 0
    firstmsg <- NULL
    for (v in val) {
        if (inherits(v, "try-error")) {
            count <- count + 1
            if (count == 1) firstmsg <- v
        }
    }
    ## These will not translate
    if (count == 1)
        stop("one node produced an error: ", firstmsg, domain = NA)
    else if (count > 1)
        stop(count, " nodes produced errors; first error: ", firstmsg, domain = NA)
    val
}

recvOneResult <- function (cl) {
    if (.snowTimingData$running()) {
        start <- proc.time()[3]
        v <- recvOneData(cl)
        end <- proc.time()[3]
        .snowTimingData$enterRecv(v$node, start, end, v$value$time[3])
    }
    else v <- recvOneData(cl)
    list(value = v$value$value, node = v$node, tag = v$value$tag)
}

findRecvOneTag <- function(cl, anytag) {
    rtag <- NULL
    for (node in cl) {
        if (is.null(rtag))
            rtag <- node$RECVTAG
        else if (rtag != node$RECVTAG) {
            rtag <- anytag
            break;
        }
    }
    rtag
}

### ========== snow support ===========

## place holder for now.
.snowTimingData <-
    list(running = function() FALSE,
         enterSend = function(...) {},
         enterRecv = function(...) {})


closeNode.NWSnode <- function(node) snow::closeNode.NWSnode(node)

recvData.MPInode <- function(node) snow::recvData.MPInode(node)
recvData.NWSnode <- function(node) snow::recvData.NWSnode(node)

recvOneData.MPIcluster <- function(cl) snow::recvOneData.MPIcluster(cl)
recvOneData.NWScluster <- function(cl) snow::recvOneData.NWScluster(cl)

sendData.MPInode <- function(node, data) snow::sendData.MPInode(node, data)
sendData.NWSnode <- function(node, data) snow::sendData.NWSnode(node, data)

## these use NextMethod() so need copies.
stopCluster.MPIcluster <- function(cl) {
    NextMethod()
    snow::setMPIcluster(NULL)
}

stopCluster.spawnedMPIcluster <- function(cl) {
    comm <- 1
    NextMethod()
    Rmpi::mpi.comm.disconnect(comm)
}

stopCluster.NWScluster <- function(cl) {
    NextMethod()
    nws::nwsDeleteWs(cl[[1]]$wsServer, nws::nwsWsName(cl[[1]]$ws))
    close(cl[[1]]$wsServer)
}

#  File src/library/parallel/R/snowSOCK.R
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

## Derived from snow 0.3-6 by Luke Tierney
## Uses solely Rscript, and a function in the package rather than scripts.

newPSOCKnode <- function(machine = "localhost", ...,
                         options = defaultClusterOptions, rank)
{
    options <- addClusterOptions(options, list(...))
    if (is.list(machine)) {
        options <- addClusterOptions(options, machine)
        machine <- machine$host
    }
    outfile <- getClusterOption("outfile", options)
    master <- if (machine == "localhost") "localhost"
    else getClusterOption("master", options)
    port <- getClusterOption("port", options)
    manual <- getClusterOption("manual", options)
    timeout <- getClusterOption("timeout", options)
    methods <- getClusterOption("methods", options)
    useXDR <- getClusterOption("useXDR", options)

    ## build the local command for starting the worker
    env <- paste0("MASTER=", master,
                 " PORT=", port,
                 " OUT=", outfile,
                 " TIMEOUT=", timeout,
                 " METHODS=", methods,
                 " XDR=", useXDR)
    arg <- "parallel:::.slaveRSOCK()"
    rscript <- if (getClusterOption("homogeneous", options)) {
        shQuote(getClusterOption("rscript", options))
    } else "Rscript"

    cmd <- paste(rscript, "-e", shQuote(arg), env)

    ## We do redirection of connections at R level once the process is
    ## running.  We could instead do it at C level here, at least on
    ## a Unix-alike.
    renice <- getClusterOption("renice", options)
    if(!is.na(renice) && renice) ## ignore 0
        cmd <- sprintf("nice +%d %s", as.integer(renice), cmd)

    if (manual) {
        cat("Manually start worker on", machine, "with\n    ", cmd, "\n")
        flush.console()
    } else {
        ## add the remote shell command if needed
        if (machine != "localhost") {
            ## This assumes an ssh-like command
            rshcmd <- getClusterOption("rshcmd", options)
            user <- getClusterOption("user", options)
            ## this assume that rshcmd will use a shell, and that is
            ## the same shell as on the master.
            cmd <- shQuote(cmd)
            cmd <- paste(rshcmd, "-l", user, machine, cmd)
        }

        if (.Platform$OS.type == "windows") {
            ## snow said:
            ## On Windows using input = something seems needed to
            ## disconnect standard input of an ssh process when run
            ## from Rterm (at least using putty's plink).  In
            ## principle this could also be used for supplying a
            ## password, but that is probably a bad idea. So, for now
            ## at least, on Windows password-less authentication is
            ## necessary.
            ##
            ## (Not clear if that is the current behaviour: works for me)
            system(cmd, wait = FALSE, input = "")
        }
        else system(cmd, wait = FALSE)
    }

    con <- socketConnection("localhost", port = port, server = TRUE,
                            blocking = TRUE, open = "a+b", timeout = timeout)
    structure(list(con = con, host = machine, rank = rank),
              class = if(useXDR) "SOCKnode" else "SOCK0node")
}

closeNode.SOCKnode <- closeNode.SOCK0node <- function(node) close(node$con)

sendData.SOCKnode <- function(node, data) serialize(data, node$con)
sendData.SOCK0node <- function(node, data) serialize(data, node$con, xdr = FALSE)

recvData.SOCKnode <- recvData.SOCK0node <- function(node) unserialize(node$con)

recvOneData.SOCKcluster <- function(cl)
{
    socklist <- lapply(cl, function(x) x$con)
    repeat {
        ready <- socketSelect(socklist)
        if (length(ready) > 0) break;
    }
    n <- which.max(ready) # may need rotation or some such for fairness
    list(node = n, value = unserialize(socklist[[n]]))
}

makePSOCKcluster <- function(names, ...)
{
    if (is.numeric(names)) names <- rep('localhost', names[1])
    options <- addClusterOptions(defaultClusterOptions, list(...))
    cl <- vector("list", length(names))
    for (i in seq_along(cl))
        cl[[i]] <- newPSOCKnode(names[[i]], options = options, rank = i)
    class(cl) <- c("SOCKcluster", "cluster")
    cl
}

print.SOCKcluster <- function(x, ...)
{
    nc <- length(x)
    hosts <- unique(sapply(x, "[[", "host"))
    msg <- sprintf(ngettext(length(hosts),
                            "socket cluster with %d nodes on host %s",
                            "socket cluster with %d nodes on hosts %s"),
                   nc, paste(sQuote(hosts), collapse = ", "))
    cat(msg, "\n", sep = "")
    invisible(x)
}

print.SOCKnode <- print.SOCK0node <- function(x, ...)
{
    sendCall(x, eval, list(quote(Sys.getpid())))
    pid <- recvResult(x)

    msg <- gettextf("node of a socket cluster on host %s with pid %d",
                    sQuote(x[["host"]]), pid)
    cat(msg, "\n", sep = "")
    invisible(x)
}

.slaveRSOCK <- function()
{
    makeSOCKmaster <- function(master, port, timeout, useXDR)
    {
        port <- as.integer(port)
        ## maybe use `try' and sleep/retry if first time fails?
        con <- socketConnection(master, port = port, blocking = TRUE,
                                open = "a+b", timeout = timeout)
        structure(list(con = con),
                  class = if(useXDR) "SOCKnode" else "SOCK0node")
    }

    ## set defaults in case run manually without args.
    master <- "localhost"
    port <- NA_integer_ # no point in getting option on worker.
    outfile <- Sys.getenv("R_SNOW_OUTFILE") # defaults to ""
    methods <- TRUE
    useXDR <- TRUE

    for (a in commandArgs(TRUE)) {
        ## Or use strsplit?
        pos <- regexpr("=", a)
        name <- substr(a, 1L, pos - 1L)
        value <- substr(a, pos + 1L, nchar(a))
        switch(name,
               MASTER = {master <- value},
               PORT = {port <- value},
               OUT = {outfile <- value},
               TIMEOUT = {timeout <- value},
               METHODS = {methods <- value},
               XDR = {useXDR <- as.logical(value)})
    }
    if (is.na(port)) stop("PORT must be specified")

    if(as.logical(methods)) library("methods") ## because Rscript does not load methods by default
    ## We should not need to attach parallel, as we are running in the namespace.

    sinkWorkerOutput(outfile)
    msg <- sprintf("starting worker pid=%d on %s at %s\n",
                   Sys.getpid(), paste(master, port, sep = ":"),
                   format(Sys.time(), "%H:%M:%OS3"))
    cat(msg)
    slaveLoop(makeSOCKmaster(master, port, timeout, useXDR))
}
#  File src/library/parallel/R/unix/forkCluster.R
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

makeForkCluster <- function(nnodes = getOption("mc.cores", 2L), ...)
{
    if(nnodes < 1L) stop("'nnodes' must be >= 1")
    cl <- vector("list", nnodes)
    for (i in seq_along(cl)) cl[[i]] <- newForkNode(..., rank = i)
    class(cl) <- c("SOCKcluster", "cluster")
    cl
}


newForkNode <- function(..., options = defaultClusterOptions, rank)
{
    options <- addClusterOptions(options, list(...))
    outfile <- getClusterOption("outfile", options)
    port <- getClusterOption("port", options)
    timeout <- getClusterOption("timeout", options)
    renice <- getClusterOption("renice", options)

    f <- mcfork()
    if (inherits(f, "masterProcess")) { # the slave
        on.exit(mcexit(1L, structure("fatal error in wrapper code",
                                  class = "try-error")))
        # closeStdout()
        master <- "localhost"
        makeSOCKmaster <- function(master, port, timeout)
        {
            port <- as.integer(port)
            ## maybe use `try' and sleep/retry if first time fails?
            con <- socketConnection(master, port = port, blocking = TRUE,
                                    open = "a+b", timeout = timeout)
            structure(list(con = con), class = "SOCK0node")
        }
        sinkWorkerOutput(outfile)
        msg <- sprintf("starting worker pid=%d on %s at %s\n",
                       Sys.getpid(), paste(master, port, sep = ":"),
                       format(Sys.time(), "%H:%M:%OS3"))
        cat(msg)
        ## allow this to quit when the loop is done.
        tools::pskill(Sys.getpid(), tools::SIGUSR1)
        if(!is.na(renice) && renice) ## ignore 0
            tools::psnice(Sys.getpid(), renice)
        slaveLoop(makeSOCKmaster(master, port, timeout))
        mcexit(0L)
    }

    con <- socketConnection("localhost", port = port, server = TRUE,
                            blocking = TRUE, open = "a+b", timeout = timeout)
    structure(list(con = con, host = "localhost", rank = rank),
              class = c("forknode", "SOCK0node"))
}
#  File src/library/parallel/R/unix/mcfork.R
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

### Derived from multicore version 0.1-6 by Simon Urbanek

### --- multicore --- low-level functions ---

## all not exported in parallel.

mc_pids <- new.env()
assign("pids", integer(), envir = mc_pids)
clean_pids <- function(e)
    if(length(pids <- get("pids", envir = e))) tools::pskill(pids, tools::SIGKILL)

mcfork <- function() {
    r <- .Call(C_mc_fork)
    assign("pids", c(get("pids",envir = mc_pids), r[1L]), envir = mc_pids)
    structure(list(pid = r[1L], fd = r[2:3]),
              class = c(if(r[1L]) "childProcess"
                        else "masterProcess", "process"))
}

## not used
readChildren <- function(timeout = 0)
    .Call(C_mc_read_children, as.double(timeout))

## used by mccollect, mclapply
readChild <- function(child)
{
    if (inherits(child, "process")) child <- processID(child)
    if (!is.numeric(child)) stop("invalid 'child' argument")
    .Call(C_mc_read_child, as.integer(child))
}

## used by mccollect, mclapply
selectChildren <- function(children = NULL, timeout = 0)
{
    if (!length(children)) children <- integer()
    if (inherits(children, "process")) children <- processID(children)
    if (is.list(children))
        children <- unlist(lapply(children, function(x) if (inherits(x, "process")) x$pid
        else stop("'children' must be a list of processes or a single process")))
    if (!is.numeric(children))
        stop("'children' must be a list of processes or a single process")
    .Call(C_mc_select_children, as.double(timeout), as.integer(children))
}

rmChild <- function(child)
{
    if (inherits(child, "process")) child <- processID(child)
    if (!is.numeric(child)) stop("invalid 'child' argument")
    .Call(C_mc_rm_child, as.integer(child))
}

## used in pvec, mclapply
mckill <- function(process, signal = 2L)
{
    process <- processID(process)
    ## or simply tools::pskill(process, signal)
    unlist(lapply(process, function(p)
                  .Call(C_mc_kill, as.integer(p), as.integer(signal))))
}

## used by mcparallel, mclapply
sendMaster <- function(what)
{
    # This is talking to the same machine, so no point in using xdr.
    if (!is.raw(what)) what <- serialize(what, NULL, xdr = FALSE)
    .Call(C_mc_send_master, what)
}

processID <- function(process) {
    if (inherits(process, "process")) process$pid
    else if (is.list(process)) unlist(lapply(process, processID))
    else stop(gettextf("'process' must be of class %s", dQuote("process")),
              domain = NA)
}

# unused
sendChildStdin <- function(child, what)
{
    if (inherits(child, "process") || is.list(child)) child <- processID(child)
    if (!is.numeric(child) || !length(child))
        stop("'child' must be a valid child process")
    child <- as.integer(child)
    if (is.character(what)) what <- charToRaw(paste(what, collapse='\n'))
    if (!is.raw(what)) stop("'what' must be a character or raw vector")
    invisible(unlist(lapply(child, function(p)
                            .Call(C_mc_send_child_stdin, p, what))))
}

## used by mcparallel, mclapply
mcexit <- function(exit.code = 0L, send = NULL)
{
    if (!is.null(send)) try(sendMaster(send), silent = TRUE)
    .Call(C_mc_exit, as.integer(exit.code))
}

## used by mccollect, mclapply
children <- function(select)
{
    p <- .Call(C_mc_children)
    if (!missing(select)) p <- p[p %in% processID(select)]
    ## FIXME: this is not the meaning of this class as returned by mcfork
    lapply(p, function(x)
           structure(list(pid = x), class = c("childProcess", "process")))
}

childrenDescriptors <- function(index = 0L)
    .Call(C_mc_fds, as.integer(index))

masterDescriptor <- function() .Call(C_mc_master_fd)

isChild <- function() .Call(C_mc_is_child)

closeStdout <- function(to.null=FALSE) .Call(C_mc_close_stdout, to.null)
closeStderr <- function(to.null=FALSE) .Call(C_mc_close_stderr, to.null)
closeFD <- function(fds) .Call(C_mc_close_fds, as.integer(fds))

closeAll <- function(includeStd = FALSE)
{
    if (!isChild()) {
        warning("closeAll() is a no-op in the master process", domain = NA)
        return(invisible(FALSE))
    }
    fds <- masterDescriptor()
    if (identical(fds, -1L)) fds <- integer(0)
    if (includeStd) fds <- c(1L, 2L, fds)
    mf <- max(fds) + 16L # take a few more ...
    ## close all but those that we actually use
    closeFD((1:mf)[-fds])
}
#  File src/library/parallel/R/unix/mclapply.R
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

### Derived from multicore version 0.1-6 by Simon Urbanek

mclapply <- function(X, FUN, ..., mc.preschedule = TRUE, mc.set.seed = TRUE,
                     mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L),
                     mc.cleanup = TRUE, mc.allow.recursive = TRUE)
{
    env <- parent.frame()
    cores <- as.integer(mc.cores)
    if(cores < 1L) stop("'mc.cores' must be >= 1")

    if (isChild() && !isTRUE(mc.allow.recursive))
        return(lapply(X = X, FUN = FUN, ...))

    if(mc.set.seed) mc.reset.stream()

    jobs <- list()
    cleanup <- function() {
        ## kill children if cleanup is requested
        if (length(jobs) && mc.cleanup) {
            ## first take care of uncollected children
            mccollect(children(jobs), FALSE)
            mckill(children(jobs),
                   if (is.integer(mc.cleanup)) mc.cleanup else tools::SIGTERM)
            mccollect(children(jobs))
        }
        if (length(jobs)) {
            ## just in case there are zombies
            mccollect(children(jobs), FALSE)
        }
    }
    on.exit(cleanup())
    ## Follow lapply
    if(!is.vector(X) || is.object(X)) X <- as.list(X)

    if (!mc.preschedule) {              # sequential (non-scheduled)
        FUN <- match.fun(FUN)
        if (length(X) <= cores) { # we can use one-shot parallel
            jobs <- lapply(seq_along(X),
                           function(i) mcparallel(FUN(X[[i]], ...),
                                                  name = names(X)[i],
                                                  mc.set.seed = mc.set.seed,
                                                  silent = mc.silent))
            res <- mccollect(jobs)
            if (length(res) == length(X)) names(res) <- names(X)
            has.errors <- sum(sapply(res, inherits, "try-error"))
        } else { # more complicated, we have to wait for jobs selectively
            sx <- seq_along(X)
            res <- vector("list", length(sx))
            names(res) <- names(X)
            ent <- rep(FALSE, length(X)) # values entered (scheduled)
            fin <- rep(FALSE, length(X)) # values finished
            jobid <- seq_len(cores)
            jobs <- lapply(jobid,
                           function(i) mcparallel(FUN(X[[i]], ...),
                                                  mc.set.seed = mc.set.seed,
                                                  silent = mc.silent))
            jobsp <- processID(jobs)
            ent[jobid] <- TRUE
            has.errors <- 0L
            while (!all(fin)) {
                s <- selectChildren(jobs, 0.5)
                if (is.null(s)) break   # no children -> no hope
                if (is.integer(s))
                    for (ch in s) {
                        ji <- which(jobsp == ch)[1]
                        ci <- jobid[ji]
                        r <- readChild(ch)
                        if (is.raw(r)) {
                            child.res <- unserialize(r)
                            if (inherits(child.res, "try-error"))
                                has.errors <- has.errors + 1L
                            ## we can't just assign it since a NULL
                            ## assignment would remove it from the list
                            if (!is.null(child.res)) res[[ci]] <- child.res
                        } else {
                            fin[ci] <- TRUE
                            if (!all(ent)) { # still something to do,
                                             # spawn a new job
                                nexti <- which(!ent)[1]
                                jobid[ji] <- nexti
                                jobs[[ji]] <- mcparallel(FUN(X[[nexti]], ...),
                                                         mc.set.seed = mc.set.seed,
                                                         silent = mc.silent)
                                jobsp[ji] <- processID(jobs[[ji]])
                                ent[nexti] <- TRUE
                            }
                        }
                    }
            }
        }
        if (has.errors)
            warning(gettextf("%d function calls resulted in an error",
                             has.errors), domain = NA)
        return(res)
    }

    ## mc.preschedule = TRUE from here on.
    if (length(X) < cores) cores <- length(X)
    if (cores < 2L) return(lapply(X = X, FUN = FUN, ...))
    sindex <- lapply(seq_len(cores),
                     function(i) seq(i, length(X), by = cores))
    schedule <- lapply(seq_len(cores),
                       function(i) X[seq(i, length(X), by = cores)])
    ch <- list()
    res <- vector("list", length(X))
    names(res) <- names(X)
    cp <- rep(0L, cores)
    fin <- rep(FALSE, cores)
    dr <- rep(FALSE, cores)
    inner.do <- function(core) {
        S <- schedule[[core]]
        f <- mcfork()
        if (isTRUE(mc.set.seed)) mc.advance.stream()
        if (inherits(f, "masterProcess")) { # this is the child process
            on.exit(mcexit(1L, structure("fatal error in wrapper code", class="try-error")))
            if (isTRUE(mc.set.seed)) mc.set.stream()
            if (isTRUE(mc.silent)) closeStdout(TRUE)
            sendMaster(try(lapply(X = S, FUN = FUN, ...), silent = TRUE))
            mcexit(0L)
        }
        jobs[[core]] <<- ch[[core]] <<- f
        cp[core] <<- f$pid
        NULL
    }
    job.res <- lapply(seq_len(cores), inner.do)
    ac <- cp[cp > 0]
    has.errors <- integer(0)
    while (!all(fin)) {
        s <- selectChildren(ac, 1)
        if (is.null(s)) break # no children -> no hope we get anything
        if (is.integer(s))
            for (ch in s) {
                a <- readChild(ch)
                if (is.integer(a)) {
                    core <- which(cp == a)
                    fin[core] <- TRUE
                } else if (is.raw(a)) {
                    core <- which(cp == attr(a, "pid"))
                    job.res[[core]] <- ijr <- unserialize(a)
                    if (inherits(ijr, "try-error"))
                        has.errors <- c(has.errors, core)
                    dr[core] <- TRUE
                }
            }
    }
    for (i in seq_len(cores)) {
        this <- job.res[[i]]
        if (inherits(this, "try-error")) { ## length-1 result
            for (j in sindex[[i]]) res[[j]] <- this
        } else res[sindex[[i]]] <- this
    }
    if (length(has.errors)) {
        if (length(has.errors) == cores)
            warning("all scheduled cores encountered errors in user code")
        else
            warning(sprintf(ngettext(has.errors,
                                     "scheduled core %s encountered error in user code, all values of the job will be affected",
                                     "scheduled cores %s encountered errors in user code, all values of the jobs will be affected"),
                            paste(has.errors, collapse = ", ")),
                    domain = NA)
    }
    res
}
#  File src/library/parallel/R/unix/mcmapply.R
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

mcmapply <-
    function(FUN, ..., MoreArgs = NULL, SIMPLIFY = TRUE, USE.NAMES = TRUE,
             mc.preschedule = TRUE, mc.set.seed = TRUE,
             mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L),
             mc.cleanup = TRUE)
{
    FUN <- match.fun(FUN)
    dots <- list(...)
    if(!length(dots)) return(list())
    lens <- sapply(dots, length)
    n <- max(lens)
    if(n && min(lens) == 0L)
        stop("Zero-length inputs cannot be mixed with those of non-zero length")
    answer <- if(n <= mc.cores) .mapply(FUN, dots, MoreArgs)
    else {
        ## recycle shorter vectors
        X <- if (!all(lens == n))
            lapply(dots, function(x) rep(x, length.out = n))
        else dots
        do_one <- function(indices, ...) {
            dots <- lapply(X, function(x) x[indices])
            .mapply(FUN, dots, MoreArgs)
        }
        answer <- mclapply(seq_len(n), do_one, mc.preschedule = mc.preschedule,
                           mc.set.seed = mc.set.seed, mc.silent = mc.silent,
                           mc.cores = mc.cores, mc.cleanup = mc.cleanup)
        do.call(c, answer)
    }
    if (USE.NAMES && length(dots)) {
        if (is.null(names1 <- names(dots[[1L]])) && is.character(dots[[1L]]))
            names(answer) <- dots[[1L]]
        else if (!is.null(names1))
            names(answer) <- names1
    }
    if (!identical(SIMPLIFY, FALSE) && length(answer))
        simplify2array(answer, higher = (SIMPLIFY == "array"))
    else answer
}

mcMap <- function (f, ...)
{
    f <- match.fun(f)
    mcmapply(f, ..., SIMPLIFY = FALSE, mc.silent = TRUE)
}
#  File src/library/parallel/R/unix/mcparallel.R
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

### Derived from multicore version 0.1-6 by Simon Urbanek

mcaffinity <- function(affinity = NULL) .Call(C_mc_affinity, affinity)

mcparallel <- function(expr, name, mc.set.seed = TRUE, silent = FALSE, mc.affinity = NULL, mc.interactive = FALSE)
{
    f <- mcfork()
    env <- parent.frame()
    if (isTRUE(mc.set.seed)) mc.advance.stream()
    if (inherits(f, "masterProcess")) {
        on.exit(mcexit(1L, structure("fatal error in wrapper code",
                                  class = "try-error")))
        if (isTRUE(mc.set.seed)) mc.set.stream()
        mc.interactive <- as.logical(mc.interactive)
        if (isTRUE(mc.interactive)) .Call(C_mc_interactive, TRUE)
        if (isTRUE(!mc.interactive)) .Call(C_mc_interactive, FALSE)
        if (!is.null(mc.affinity)) .Call(C_mc_affinity, mc.affinity)
        if (isTRUE(silent)) closeStdout(TRUE)
        sendMaster(try(eval(expr, env), silent = TRUE))
        mcexit(0L)
    }
    if (!missing(name) && !is.null(name)) f$name <- as.character(name)[1L]
    class(f) <- c("parallelJob", class(f))
    f
}

mccollect <- function(jobs, wait = TRUE, timeout = 0, intermediate = FALSE)
{
    if (missing(jobs)) jobs <- children()
    if (!length(jobs)) return (NULL)
    if (isTRUE(intermediate)) intermediate <- str
    if (!wait) {
        s <- selectChildren(jobs, timeout)
        if (is.logical(s) || !length(s)) return(NULL)
        lapply(s, function(x) {
            r <- readChild(x)
            if (is.raw(r)) unserialize(r) else NULL
        })
    } else {
        pids <- if (inherits(jobs, "process") || is.list(jobs))
            processID(jobs) else jobs
        if (!length(pids)) return(NULL)
        if (!is.numeric(pids)) stop("invalid 'jobs' argument")
        pids <- as.integer(pids)
        pnames <- as.character(pids)
        if (!inherits(jobs, "process") && is.list(jobs))
            for(i in seq(jobs))
                if (!is.null(jobs[[i]]$name))
                    pnames[i] <- as.character(jobs[[i]]$name)
        res <- lapply(pids, function(x) NULL)
        names(res) <- pnames
        fin <- rep(FALSE, length(jobs))
        while (!all(fin)) {
            s <- selectChildren(pids, 0.5)
            if (is.integer(s)) {
                for (pid in s) {
                    r <- readChild(pid)
                    if (is.integer(r) || is.null(r)) fin[pid == pids] <- TRUE
                    if (is.raw(r)) # unserialize(r) might be null
                        res[which(pid == pids)] <- list(unserialize(r))
                }
                if (is.function(intermediate)) intermediate(res)
            } else if (all(is.na(match(pids, processID(children()))))) break
        }
        res
    }
}
#  File src/library/parallel/R/unix/pvec.R
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

### Derived from multicore version 0.1-6 by Simon Urbanek

pvec <- function(v, FUN, ..., mc.set.seed = TRUE, mc.silent = FALSE,
                 mc.cores = getOption("mc.cores", 2L), mc.cleanup = TRUE)
{
    if (!is.vector(v)) stop("'v' must be a vector")

    env <- parent.frame()
    cores <- as.integer(mc.cores)
    if(cores < 1L) stop("'mc.cores' must be >= 1")
    if(cores == 1L) return(FUN(v, ...))

    if(mc.set.seed) mc.reset.stream()

    n <- length(v)
    l <- if (n <= cores) as.list(v) else {
        ## compute the scheduling, making it as fair as possible
        il <- as.integer(n / cores)
        xc <- n - il * cores
        sl <- rep(il, cores)
        if (xc) sl[1:xc] <- il + 1L
        si <- cumsum(c(1L, sl))
        se <- si + c(sl, 0L) - 1L
        lapply(seq_len(cores), function(ix) v[si[ix]:se[ix]])
    }
    jobs <- NULL
    cleanup <- function() {
        ## kill children if cleanup is requested
        if (length(jobs) && mc.cleanup) {
            ## first take care of uncollected children
            mccollect(children(jobs), FALSE)
            mckill(children(jobs),
                   if (is.integer(mc.cleanup)) mc.cleanup else 15L)
            mccollect(children(jobs))
        }
        if (length(jobs)) {
            ## just in case there are zombies
            mccollect(children(jobs), FALSE)
        }
    }
    on.exit(cleanup())
    FUN <- match.fun(FUN)
    ## may have more cores than tasks ....
    jobs <- lapply(seq_len(min(n, cores)),
                   function(i) mcparallel(FUN(l[[i]], ...), name = i,
                                          mc.set.seed = mc.set.seed,
                                          silent = mc.silent))
    res <- mccollect(jobs)
    names(res) <- NULL
    res <- do.call(c, res)
    if (length(res) != n)
        warning("some results may be missing, folded or caused an error")
    res
}
#  File src/library/parallel/R/worker.R
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

## Derived from snow 0.3-6 by Luke Tierney

slaveLoop <- function(master)
{
    repeat
        tryCatch({
            msg <- recvData(master)
            # cat(paste("Type:", msg$type, "\n"))

            if (msg$type == "DONE") {
                closeNode(master)
                break;
            } else if (msg$type == "EXEC") {
                success <- TRUE
                ## This uses the message rather than the exception since
                ## the exception class/methods may not be available on the
                ## master.
                handler <- function(e) {
                    success <<- FALSE
                    structure(conditionMessage(e),
                              class = c("snow-try-error","try-error"))
                }
                t1 <- proc.time()
                value <- tryCatch(do.call(msg$data$fun, msg$data$args, quote = TRUE),
                                  error = handler)
                t2 <- proc.time()
                value <- list(type = "VALUE", value = value, success = success,
                              time = t2 - t1, tag = msg$data$tag)
                sendData(master, value)
            }
        }, interrupt = function(e) NULL)
}

## NB: this only sinks the connections, not C-level stdout/err.
sinkWorkerOutput <- function(outfile)
{
    if (nzchar(outfile)) {
        if (.Platform$OS.type == "windows" && outfile == "/dev/null")
            outfile <- "nul:"
        ## all the workers log to the same file.
        outcon <- file(outfile, open = "a")
        sink(outcon)
        sink(outcon, type = "message")
    }
}

#  File src/library/parallel/R/zzz.R
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

.noGenerics <- TRUE

if (.Platform$OS.type == "windows")
    utils::globalVariables(c("mc_pids", "clean_pids"), add = TRUE)

.onLoad <- function(libname, pkgname)
{
    initDefaultClusterOptions(libname)
    cores <- getOption("mc.cores", NULL)
    if(is.null(cores) && !is.na(nc <- as.integer(Sys.getenv("MC_CORES"))))
        options("mc.cores" = nc)
    if(.Platform$OS.type == "unix") reg.finalizer(mc_pids, clean_pids, TRUE)
}

.onUnload <-
function(libpath)
    library.dynam.unload("parallel", libpath)
