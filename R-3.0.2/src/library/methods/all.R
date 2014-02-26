.packageName <- "methods"
#  File src/library/methods/R/BasicClasses.R
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

## a few class name definitions needed elsewhere
.anyClassName <- structure("ANY", package = "methods")
.signatureClassName <- structure("signature", package = "methods")



.InitBasicClasses <- function(envir)
{
    ## setClass won't allow redefining basic classes,
    ## so make the list of these empty for now.
    assign(".BasicClasses", character(), envir)
    ## hide some functions that would break because the basic
    ## classes are not yet defined
    real.reconcileP <- reconcilePropertiesAndPrototype
    assign("reconcilePropertiesAndPrototype",
           function(name, properties, prototype, extends, where) {
               list(properties=properties, prototype = prototype, extends = extends)
           }, envir)
    clList = character()
    setClass("VIRTUAL", where = envir); clList <- c(clList, "VIRTUAL")
    setClass("ANY", where = envir); clList <- c(clList, "ANY")
    setClass("vector", where = envir); clList <- c(clList, "vector")
    setClass("missing", where = envir); clList <- c(clList, "missing")
    ## "numeric" is the class returned by class() for double vectors
    vClasses <- c("logical", "numeric", "character",
                  "complex", "integer", "raw",
                  "expression", "list")
    ## now some pseudo-classes in base, marked specially for new()
    for(.class in vClasses) {
        .setBaseClass(.class, prototype = newBasic(.class), where = envir)
    }
    .setBaseClass("expression", prototype = expression(), where = envir)
    clList <- c(clList, vClasses)
    nullF <- function()NULL; environment(nullF) <- .GlobalEnv
    attr(nullF, "source") <- NULL
    .setBaseClass("function", prototype = nullF, where = envir); clList <- c(clList, "function")

    setClass("language", where = envir); clList <- c(clList, "language")
    .setBaseClass("environment", prototype = new.env(), where = envir); clList <- c(clList, "environment")

    .setBaseClass("externalptr", prototype = .newExternalptr(), where = envir); clList <- c(clList, "externalptr")

    .setBaseClass("builtin", prototype = `<-`, where = envir); clList <- c(clList, "builtin")

    .setBaseClass("special", prototype = `if`, where = envir); clList <- c(clList, "special")

    ## S4, S3 are basic classes that are used to define methods related to being S4, S3 object
    for(cl in c("S4", "S3")) {
        tmp <- newClassRepresentation(className=cl, prototype = defaultPrototype(), virtual=TRUE, package = "methods")
        assignClassDef(cl, tmp, where = envir); clList <- c(clList, cl)
    }

    ## NULL is weird in that it has NULL as a prototype, but is not virtual
    tmp <- newClassRepresentation(className="NULL", prototype = NULL, virtual=FALSE, package = "methods")
    assignClassDef("NULL", tmp, where = envir); clList <- c(clList, "NULL")
    ## the pseudo-NULL used to store NULL as a slot
    ## must match the C code in attrib.c (would be better to use that
    ## code to create .pseudoNULL)
    assign(".pseudoNULL", as.name("\001NULL\001"), envir = envir)


    setClass("structure", where = envir); clList <- c(clList, "structure")
    setClass("nonStructure",  where = envir); #NOT a basic class
    stClasses <- c("matrix", "array") # classes that have attributes, but no class attr.
    for(.class in stClasses) {
        .setBaseClass(.class, prototype = newBasic(.class), where = envir)
    }
    ## "ts" will be defined below as an S3 class, but it is still
    ## included in .BasicClasses, to allow its coerce() method to use
    ## as.ts().  This decision may be revisited.
    clList <- c(clList, stClasses, "ts")
    assign(".BasicClasses", clList, envir)

    ## Now we can define the SClassExtension class and use it to instantiate some
    ## is() relations.
    .InitExtensions(envir)

    for(.class in vClasses)
        setIs(.class, "vector", where = envir)

    setIs("integer", "numeric", where = envir)

    setIs("structure", "vector", coerce = .gblEnv(function(object) as.vector(object)),
          replace = .gblEnv(function(from, to, value) {
              attributes(value) <- attributes(from)
              value
          }),
          where = envir)

    setIs("array", "structure", where = envir)
    setIs("matrix", "array", where = envir)
### Rather want a simple  setAs("array", "matrix", ..) method..
    ## setIs("array", "matrix", test = .gblEnv(function(object) length(dim(object)) == 2),
    ##       replace = .gblEnv(function(from, to, value) {
    ##           if(is(value, "matrix"))
    ##               value
    ##           else
    ##               stop("replacement value is not a matrix")
    ##       }),
    ##       where = envir)

    ## Some class definitions extending "language", delayed to here so
    ## setIs will work.
    .setBaseClass("name", "language", prototype = as.name("<UNDEFINED>"), where = envir); clList <- c(clList, "name")
    .setBaseClass("call", "language", prototype = quote("<undef>"()), where = envir); clList <- c(clList, "call")
    .setBaseClass("{", "language", prototype = quote({}), where = envir); clList <- c(clList, "{")
    .setBaseClass("if", "language", prototype = quote(if(NA) TRUE else FALSE), where = envir); clList <- c(clList, "if")
    .setBaseClass("<-", "language", prototype = quote("<undef>"<-NULL), where = envir); clList <- c(clList, "<-")
    .setBaseClass("for", "language", prototype = quote(for(NAME in logical()) NULL), where = envir); clList <- c(clList, "for")
    .setBaseClass("while", "language", prototype = quote(while(FALSE) NULL), where = envir); clList <- c(clList, "while")
    .setBaseClass("repeat", "language", prototype = quote(repeat{break}), where = envir); clList <- c(clList, "repeat")
    .setBaseClass("(", "language", prototype = quote((NULL)), where = envir); clList <- c(clList, "(")

    ## a virtual class used to allow NULL as an indicator that a possible function
    ## is not supplied (used, e.g., for the validity slot in classRepresentation
    setClass("OptionalFunction", where = envir)
    setIs("function", "OptionalFunction", where = envir)
    setIs("NULL", "OptionalFunction")
    assign(".BasicClasses", clList, envir)
    assign(".SealedClasses", clList, envir)
    ## restore the true definition of the hidden functions
    assign("reconcilePropertiesAndPrototype", real.reconcileP, envir)
}

.InitS3Classes <- function(envir) {
    ## create a virtual class from which all S3 classes will inherit the .S3Class slot
    setClass("oldClass", representation(.S3Class = "character"),
             contains = "VIRTUAL", prototype = prototype(.S3Class = character()),
             where = envir)
    ## call setOldClass on some known old-style classes.  Ideally this would be done
    ## in the code that uses the classes, but that code doesn't know about the methods
    ## package.
    ## Two steps; first, those classes with a known prototype.  These
    ## can be non-Virtual
    clList <- get(".SealedClasses", envir = envir)
    for(i in seq_along(.OldClassesPrototypes)) {
        el <- .OldClassesPrototypes[[i]]
        if(is.list(el) && length(el) > 1)
            setOldClass(el[[1L]], prototype = el[[2L]],  where = envir)
        else
            warning(gettextf("OOPS: something wrong with line %d in '.OldClassesPrototypes'", i), domain = NA)
    }
    setGeneric("slotsFromS3", where = envir)
    ## the method for "oldClass" is really a constant, just hard to express that way
    setMethod("slotsFromS3", "oldClass", function(object) getClass("oldClass")@slots,
              where = envir)

    setClass("ts", contains = "structure", representation(tsp = "numeric"),
             prototype = prototype(NA, tsp = rep(1,3)), where = envir)

    setOldClass("ts", S4Class = "ts", where = envir)

    setClass("mts", contains=c("matrix", "ts"), prototype =
             prototype(matrix(NA,1,1), tsp = rep(1,3), .S3Class = c("mts", "ts")))
    .init_ts <-	 function(.Object,  ...) {
	if(nargs() < 2) # guaranteed to be called with .Object from new
	    return(.Object)
	args <- list(...)
	argnames <- names(args)
	slotnames <- if(is.null(argnames)) FALSE else {
            nzchar(argnames) & is.na(match(argnames, .tsArgNames)) }
	if(any(slotnames)) {
	    value <- do.call(stats::ts, args[!slotnames])
	    .mergeAttrs(value, .Object, args[slotnames])
	}
	else
	    .mergeAttrs(stats::ts(...), .Object)
    }
    setMethod("initialize", "ts", .init_ts, where = envir)
    setMethod("initialize", "mts", .init_ts, where = envir) #else, it's ambiguous
    ## the following mimics settings for other basic classes ("ts" was
    ## not defined at the time these are done).
    setMethod("coerce", c("ANY", "ts"), function (from, to, strict = TRUE)
              {
                  value <- as.ts(from)
                  if(strict) {
                      attrs <- attributes(value)
                      if(length(attrs) > 2)
                        attributes(value) <- attrs[c("class", "tsp")]
                      value <- .asS4(value)
                  }
                  value
              },
              where = envir)
    setClass("factor", contains = "integer", representation(levels = "character"),
	     validity = function(object) {
		 levs <- levels(object)
		 if (!is.character(levs))
		     return("factor levels must be \"character\"")
		 if (d <- anyDuplicated(levs))
		     return(sprintf("duplicated level [%d] in factor", d))
		 ## 'else'	ok :
		 TRUE
	     },
	     where = envir)
    setOldClass("factor", S4Class = "factor", where = envir)
    if(!isGeneric("show", envir))
        setGeneric("show", where = envir, simpleInheritanceOnly = TRUE)
    setMethod("show", "oldClass", function(object) {
        if(!isS4(object))  {
            print(object)
            return(invisible())
        }
        cl <- as.character(class(object))
        S3Class <- object@.S3Class
        if(length(S3Class)) S3Class <- S3Class[[1L]]
        else S3Class <- "oldClass"      # or error?
        cat("Object of class \"", cl, "\"\n", sep = "")
        print(S3Part(object, strictS3 = TRUE))
        otherSlots <- slotNames(cl)
        S3slots <- slotNames(S3Class)
        otherSlots <- otherSlots[is.na(match(otherSlots, S3slots))]
        for(what in otherSlots) {
            cat('Slot "', what, '":\n', sep = "")
            show(slot(object, what))
            cat("\n")
        }
        NULL
    }, where = envir)
   .initS3 <- function(.Object, ...) {
         if(nargs() < 2)
           return(.Object)
         Class <- class(.Object)
         ClassDef <- getClass(Class)
         S3Class <- attr(ClassDef@prototype, ".S3Class")
         if(is.null(S3Class)) # not a class set up by setOldClass()
             return(callNextMethod())
        S3ClassP <- S3Class[[1L]]
         args <- list(...)
        ## separate the slots, superclass objects
        snames <- allNames(args)
        which <- nzchar(snames)
        elements <- args[which]
        supers <- args[!which]
        thisExtends <- names(ClassDef@contains)
        slotDefs <- ClassDef@slots
        dataPart <- elNamed(slotDefs, ".Data")
        if(is.null(dataPart))
          dataPart <- "missing" # nothing will extend this => no data part args allowed
        if(length(supers) > 0) {
            for(i in rev(seq_along(supers))) {
                obj <- el(supers, i)
                Classi <- class(obj)
                defi <- getClassDef(Classi)
                if(is.null(defi))
                  stop(gettextf("unnamed argument to initialize() for S3 class must have a class definition; %s does not",
                                dQuote(Classi)),
                       domain = NA)
                if(is(obj, S3ClassP)) {
                    ## eligible to be the S3 part; merge other slots from prototype;
                    ## obj then becomes the object, with its original class as the S3Class
                    if(is.null(attr(obj, ".S3Class"))) # must be an S3 object; use its own class
                       attr(obj, ".S3Class") <- Classi
                    .Object <- .asS4(.mergeAttrs(obj, .Object))
                }
                else if(is(obj, dataPart)) {
                    ## the S3Class stays from the prototype
                    .Object <- .mergeAttrs(obj, .Object)
                }
                else stop(gettextf("unnamed argument must extend either the S3 class or the class of the data part; not true of class %s", dQuote(Classi)), domain = NA)

            }
        }
        ## named slots are done as in the default method, which will also call validObject()
        if(length(elements)>0) {
            elements <- c(list(.Object), elements)
            .Object <- do.call(`callNextMethod`, elements)
        }
         else
           validObject(.Object)
         .Object
    }
    setMethod("initialize", "oldClass", .initS3, where = envir)
    ## Next, miscellaneous S3 classes.
    for(cl in .OldClassesList)
        setOldClass(cl, where = envir)
    ## special mess for "maov"; see comment in .OldClassesList
    setIs("maov", "aov")
    setClassUnion("data.frameRowLabels", c("character", "integer"), where = envir)
    setClass("data.frame",
             representation(names = "character", row.names = "data.frameRowLabels"),
             contains = "list", prototype = unclass(data.frame()), where = envir) # the S4 version
    setOldClass("data.frame", S4Class = "data.frame", where = envir)
    ## the S3 method for $<- does some stupid things to class()
    ## This buffers the effect from S4 classes
    setMethod("$<-", "data.frame", where = envir,
              function(x, name, value) {
                  x@.Data <- as.list(`$<-.data.frame`(structure(x@.Data, names = x@names,
                         row.names = x@row.names, class = "data.frame"),
                     name, value))
                  ## Assert:  the only slot/attribute that can change
                  ## in $<-.data.frame is "names", and the assignment
                  ## of the .Data "slot" copies in the new names
                  x
              })
    ## methods to go from S4 to S3; first, using registered class; second, general S4 object
    setMethod("coerce", c("oldClass", "S3"), function (from, to, strict = TRUE)
              {
                  from <- .notS4(from) # not needed? ensures that class() can return >1 string
                  cl <- class(from)
                  cl1 <- .class1(from)
                  classDef <- getClassDef(cl1)
                  S3Class <- attr(classDef@prototype, ".S3Class")
                  if(length(S3Class) > length(cl))  #add S3 inheritance
                      attr(from, "class") <- S3Class
                  from
              },
              where = envir)
    setMethod("coerce", c("ANY", "S3"), function (from, to, strict = TRUE)
              {
                  switch(typeof(from),
                         S4 =
                         stop(gettextf("class %s does not have an S3 data part, and so is of type \"S4\"; no S3 equivalent",
                                       dQuote(class(from))),
                              domain = NA),
                         .notS4(from) )
              },
              where = envir)
    setMethod("coerce", c("ANY", "S4"), function (from, to, strict = TRUE)
              {
                  if(isS4(from)) {
                      value <- from
                  }
                  else {
                      cl <- .class1(from)
                      classDef <- getClass(cl)
                      if(identical(classDef@virtual, TRUE))
                        stop(gettextf("class %s is VIRTUAL; not meaningful to create an S4 object from this class",
                                      dQuote(cl)),
                             domain = NA)
                      pr <- classDef@prototype
                      value <- new(cl)
                      slots <- classDef@slots
                      if(match(".Data", names(slots), 0L) > 0L) {
                          data <- unclass(from)
                          if(!is(data, slots[[".Data"]]))
                            stop(gettextf("object must be a valid data part for class %s; not true of type %s", dQuote(cl), dQuote(class(data))),
                                 domain = NA)
                          value@.Data <- unclass(from)
                      }
                      ## copy attributes:  Note that this copies non-slots as well
                      ## but checks the slots for validity
                      anames <- names(attributes(from))
                      isSlot <- anames %in% names(slots)
                      for(i in seq_along(anames)) {
                          what <- anames[[i]]
                          if(isSlot[[i]])
                            slot(value, what) <- attr(from, what)
                          else
                            attr(value, what) <- attr(from, what)
                      }
                  }
                  if(strict)
                    ## validate.  If we created S4 object, slots were tested; else, not
                    ## so complete= is set accordingly.
                      validObject(value, complete = isS4(from))
                  value
              })
    assign(".SealedClasses", c(clList,unique(unlist(.OldClassesList))),  envir)
}

### create a class definition for one of the pseudo-classes in base
### The class name does _not_ have a package attribute, which signals
### the C coded for new() to return an object w/o explicit class
### attribute, to be consistent with older R code
.setBaseClass <- function(cl, ..., where) {
    setClass(cl, ..., where = where)
    def <- getClassDef(cl, where)
    def@className <- as.character(def@className)
    def@prototype <- .notS4(def@prototype)
    assignClassDef(cl, def, where = where)
}


.tsArgNames <- names(formals(stats::ts))

### The following methods are now activated
### via the last line of the function .InitMethodDefinitions in ./MethodsListClass.R
###
### Tradeoff between intuition of users that
### new("matrix", ...) should be like matrix(...) vs consistency of new().
### Relevant when new class has basic class as its data part.
.InitBasicClassMethods <- function(where) {
    ## methods to initialize "informal" classes by using the
    ## functions of the same name.

    ## These methods are designed to be inherited or extended
    initMatrix <- function(.Object, data = NA, nrow = 1, ncol = 1,
                           byrow = FALSE, dimnames = NULL, ...) {
        na <- nargs()
        if(length(dots <- list(...)) && ".Data" %in% names(dots)) {
            if(na == 2)
              .Object <- .mergeAttrs(dots$.Data, .Object)
            else {
                dat <- dots$.Data
                dots <- dots[names(dots) != ".Data"]
                if(na == 2 + length(dots)) {
                    .Object <- .mergeAttrs(as.matrix(dat), .Object, dots)
                }
                else
                  stop("cannot specify matrix() arguments when specifying '.Data'")
            }
        }
        else if(is.matrix(data) && na == 2 + length(dots))
          .Object <- .mergeAttrs(data, .Object, dots)
        else {
            if (missing(nrow))
              nrow <- ceiling(length(data)/ncol)
            else if (missing(ncol))
              ncol <- ceiling(length(data)/nrow)
            value <- matrix(data, nrow, ncol, byrow, dimnames)
            .Object <- .mergeAttrs(value, .Object, dots)
        }
        validObject(.Object)
        .Object
    }
    .matrixExtends <- unique(c("matrix", names(getClass("matrix")@contains)))
    setMethod("initialize", "matrix", where = where,
              function(.Object, ...) {
		  if(nargs() < 2) # guaranteed to be called with .Object from new
                    return(.Object)
		  else {
                      if(isMixin(getClass(class(.Object)))) # other superclasses
                          callNextMethod()
                      else
                          initMatrix(.Object, ...)
                  }
              }
	      )
    initArray <- function(.Object, data = NA, dim = length(data),
                          dimnames = NULL, ...) {
        na <- nargs()
        if(length(dots <- list(...)) && ".Data" %in% names(dots)) {
            if(na == 2)
              .Object <- .mergeAttrs(dots$.Data, .Object)
            else {
                dat <- dots$.Data
                dots <- dots[names(dots) != ".Data"]
                if(na == 2 + length(dots)) {
                    .Object <- .mergeAttrs(as.array(dat), .Object, dots)
                }
                else
                  stop("cannot specify array() arguments when specifying '.Data'")
            }
        }
        else if(is.array(data) && na == 2 + length(dots))
          .Object <- .mergeAttrs(data, .Object, dots)
        else {
            value <- array(data, dim, dimnames)
            .Object <- .mergeAttrs(value, .Object, dots)
        }
        validObject(.Object)
        .Object
    }
    .arrayExtends <- unique(c("array", names(getClass("array")@contains)))
    setMethod("initialize", "array", where = where,
              function(.Object, ...) {
		  if(nargs() < 2) # guaranteed to be called with .Object from new
                    .Object
		  else {
                      if(isMixin(getClass(class(.Object)))) # other superclasses
                          callNextMethod()
                      else
                          initArray(.Object, ...)
                  }
              }
	      )
    ## following should not be needed if data_class2 returns "array",...
##     setMethod("[", # a method to avoid invalid objects from an S4 class
##               signature(x = "array"), where = where,
##               function (x, i, j, ..., drop = TRUE)
##               {
##                 value <- callNextMethod()
##                 if(is(value, class(x)))
##                   value@.Data
##                 else
##                   value
##               })

}

## .OldClassList is a purely heuristic list of known old-style classes, with emphasis
## on old-style class inheritiance.  Used in .InitBasicClasses to call setOldClass for
## each known class pattern.
## .OldClassesPrototypes is a list of S3 classes for which prototype
## objects are known & reasonable.  The classes will reappear in
## .OldClassesList, but will have been initialized first in
## .InitBasicClasses.  NB:  the methods package will NOT set up
## prototypes for S3 classes except those in package base and for "ts"
## (and would rather not do those either).  The package that owns the
## S3 class should have code to call setOldClass in its
## initialization.
.OldClassesPrototypes <-
  list(
       list("data.frame",  data.frame(), "data.frame"),
       list("factor",  factor()),
       list(c("ordered", "factor"), ordered(character())),
       list("table",  table(factor())),
       list("summary.table",  summary.table(table(factor())))
       , list("ts", stats::ts())
       )
.OldClassesList <-
    list(
         c("anova", "data.frame"),
         c("mlm", "lm"),
         c("aov", "lm"),
         ## note:  definition of "maov" below differs from the
         ## current S3 attribute, which has an inconsistent combination
         ## of "aov" and "mlm" (version 2.12 devel, rev. 51984)
         c("maov", "mlm", "lm"),
         c("POSIXct", "POSIXt"),
         c("POSIXlt", "POSIXt"),
         "Date",
         "dump.frames",
         c("ordered", "factor"),
         c("glm.null", "glm", "lm"),
         c("anova.glm.null", "anova.glm"),
         "hsearch",
         "integrate",
         "packageInfo",
         "libraryIQR",
         "packageIQR",
         "mtable",
         "table",
         c("summaryDefault","table"),
         "summary.table",
         "recordedplot",
         "socket",
         "packageIQR",
         "density",
         "formula",
         "logLik",
         "rle"
)

.InitSpecialTypesAndClasses <- function(where) {
    if(!exists(".S3MethodsClasses", envir = where, inherits = FALSE)) {
      S3table <- new.env()
      assign(".S3MethodsClasses", S3table, envir = where)
    }
    else S3table <- get(".S3MethodsClasses", envir = where)
    specialClasses <- .indirectAbnormalClasses
    specialTypes <- .AbnormalTypes # only part matching classes used
    for(i in seq_along(specialClasses)) {
        cl <- specialTypes[[i]]
      ncl <- specialClasses[[i]]
      setClass(ncl, representation(.xData = cl), where = where)
      setIs(ncl, cl, coerce = function(from) from@.xData,
        replace = function(from, value){ from@.xData <- value; from},
        where = where)
      ## these classes need explicit coercion for S3 methods
      assign(cl, getClass(cl, where), envir = S3table)
    }
    ## a few other special classes
    setClass("namedList", representation(names = "character"),
             contains = "list", where = where)
    if(!isGeneric("show", where))
        setGeneric("show", where = where, simpleInheritanceOnly = TRUE)
    setMethod("show", "namedList", function(object) {
        cat("An object of class ", dQuote(class(object)), "\n")
        print(structure(object@.Data, names=object@names))
        showExtraSlots(object, getClass("namedList"))
    })
    setClass("listOfMethods", representation(  arguments = "character",
                                                 signatures = "list", generic = "genericFunction"), contains = "namedList",
             where = where)
    specialClasses <- c(specialClasses, "namedList", "listOfMethods")
    assign(".SealedClasses", c(get(".SealedClasses", where), specialClasses), where)
    setMethod("initialize", ".environment", # for simple subclasses of "environment"
              function(.Object, ...) {
                  args <- list(...)
                  objs <- names(args)
                  hasEnvArg <- length(args) && !all(nzchar(objs))
                  if(hasEnvArg) {
                      ii <- seq_along(args)[!nzchar(objs)]
                      i <- integer()
                      for(iii in ii) {
                          if(is(args[[iii]], "environment"))
                              i <- c(i, iii)
                      }
                      if(length(i)>1)
                          stop("cannot have more than one unnamed argument as environment")
                      if(length(i) == 1) {
                          selfEnv <- args[[i]]
                          args <- args[-i]
                          objs <- objs[-i]
                          if(!is(selfEnv, "environment"))
                              stop("unnamed argument to new() must be an environment for the new object")
                          selfEnv <- as.environment(selfEnv)
                      }
                      ## else, no environment superclasses
                      else
                          selfEnv <- new.env()
                  }
                  else
                      selfEnv <- new.env()
                  if(length(objs)) {
                      ## don't assign locally named slots of subclasses
                      ClassDef <- getClass(class(.Object))
                      slots <- slotNames(ClassDef)
                      localObjs <- is.na(match(objs, slots))
                      if(any(localObjs)) {
                          for(what in objs[localObjs])
                              assign(what, elNamed(args, what), envir = selfEnv)
                          objs <- objs[!localObjs]
                          args <- args[!localObjs]
                      }
                  }
                  .Object@.xData <- selfEnv
                  if(length(objs)) # call next method with remaining args
                      .Object <- do.call(callNextMethod, c(.Object, args))
                  .Object
              }, where = where)
}
#  File src/library/methods/R/BasicFunsList.R
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

## Lists of functions and expressions used in dispatch of functions
## defined internally (as .Primitive's) for which formal argument lists
## are not available, or for which a generic, if created,
## needs to have a special form (e.g., belonging to one of the
## predefined groups of functions).

## The list is expanded in .makeBasicFuns by adding the S4 group generics
## and the remaining primitives.

.BasicFunsList <-
list(
### subset/subassignment ops are regarded as language elements
"$" = structure(function(x, name)
{
    name <- as.character(substitute(name))
    standardGeneric("$")
}, signature = c("x"))
, "$<-" = structure(function(x, name, value)
{
    name <- as.character(substitute(name))
    standardGeneric("$<-")
}, signature = c("x", "value"))
, "[" = function(x, i, j, ..., drop = TRUE) standardGeneric("[")
, "[<-" = function(x, i, j, ..., value) standardGeneric("[<-")
, "[[" = function(x, i, j, ...) standardGeneric("[[")
, "[[<-" = function(x, i, j, ..., value) standardGeneric("[[<-")
### S4 generic via R_possible_dispatch in do_matprod
, "%*%" = function(x, y) standardGeneric("%*%")
, "xtfrm" = function(x) standardGeneric("xtfrm")
### these have a different arglist from the primitives
, "c" = function(x, ..., recursive = FALSE) standardGeneric("c")
, "all" = function(x, ..., na.rm = FALSE) standardGeneric("all")
, "any" = function(x, ..., na.rm = FALSE) standardGeneric("any")
, "sum" = function(x, ..., na.rm = FALSE) standardGeneric("sum")
, "prod" = function(x, ..., na.rm = FALSE) standardGeneric("prod")
, "max" = function(x, ..., na.rm = FALSE) standardGeneric("max")
, "min" = function(x, ..., na.rm = FALSE) standardGeneric("min")

, "range" = function(x, ..., na.rm = FALSE) standardGeneric("range")
## , "!" = function(e1) standardGeneric("!")
)

## the names of the basic funs with the style of "["
## R implements these in an inconsistent call mechanism, in which missing arguments
## are allowed, and significant, but argument names are not used.  See callNextMethod

.BasicSubsetFunctions <- c("[", "[[", "[<-", "[[<-")

## create generic functions corresponding to the basic (primitive) functions
## but don't leave them as generics in the package.  Instead store them in
## a named list to be used by setMethod, w/o forcing method dispatch on these
## functions.

.addBasicGeneric <-
    function(funslist, f, fdef, group = list())
{
    signature <- attr(fdef, "signature") #typically NULL, but see the case for "$"
    deflt <- get(f, "package:base")
    ## use the arguments of the base package function
    ##FIXME:  should also deal with the functions having ... as the first
    ## argument, but needs to create a generic with different args from the deflt
    ## => constructing a call to the base function from the default
    if(is.primitive(deflt)) {
        body(fdef, envir = globalenv()) <-
            substitute(standardGeneric(FNAME, DEFLT), list(FNAME=f, DEFLT=deflt))
    }
    else {
        fdef <- deflt
        body(fdef, envir = globalenv()) <-
            substitute(standardGeneric(FNAME), list(FNAME=f))
    }
    deflt <- .derivedDefaultMethod(deflt)
    elNamed(funslist, f) <- makeGeneric(f, fdef, deflt, group = group, package = "base",
                                        signature = signature)
    funslist
}

.ShortPrimitiveSkeletons <-
    list( quote(f(x,i)), quote(fgets(x,i,value=value)))

.EmptyPrimitiveSkeletons <-
    list( quote(f(x)), quote(fgets(x,value=value)))

## utilities to get and set the primitive generics.
## Version below uses the environment, not the list
## in order to work with namespace for methods package
# genericForPrimitive <- function(f, where = topenv(parent.frame())) {
#     what <- methodsPackageMetaName("G", f)
#     if(exists(what, where))
#         get(what, where)
#     else
#         NULL
# }

# setGenericForPrimitive <-function(f, value, where = topenv(parent.frame()))
#     assign(methodsPackageMetaName("G", f), value, where)

## temporary versions while primitives are still handled by a global table

genericForPrimitive <- function(f, where = topenv(parent.frame()), mustFind = TRUE) {
#    if(.matchBasic(f, .ExcludePrimitiveGenerics, FALSE))
#        stop(gettextf("methods may not be defined for primitive function %s in this version of R", sQuote(f)), domain = NA)
    env <- .findBasicFuns(where)
    funs <- get(".BasicFunsList", envir = env)
    ans <- elNamed(funs, f)
    ## this element may not exist (yet, during loading), dom't test null
    if(mustFind && identical(ans, FALSE))
        stop(gettextf("methods may not be defined for primitive function %s in this version of R",
                      sQuote(f)),
             domain = NA)
    ans
}

.findBasicFuns <- function(where) {
    allWhere <- .findAll(".BasicFunsList", where = where)
    if(length(allWhere) == 0)
        .methodsNamespace
    else
        as.environment(allWhere[[1L]])
}
#  File src/library/methods/R/ClassExtensions.R
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

.InitExtensions <- function(where) {
    ## to be called from the initialization
    setClass("SClassExtension",
	     representation(subClass = "character", superClass = "character",
			    package = "character", coerce = "function",
			    test = "function", replace = "function",
			    simple = "logical", by = "character",
			    dataPart = "logical", distance = "numeric"),
	     where = where)
    ## a class for conditional extensions, so they will not break the hierarchical
    ## structure.
    setClass("conditionalExtension", contains = "SClassExtension")
    assign(".SealedClasses", c(get(".SealedClasses", where), "SClassExtension",
                               "conditionalExtension"),
	   where)
}

.simpleExtCoerce <- function(from, strict = TRUE)from
.simpleIsCoerce <- function(from)from
.simpleExtTest <- function(object)TRUE
## TO DO:  the simple replace below only handles the case of classes with slots.
## There are some other simple relations (e.g., with virtual classes).  Replacing in
## these cases is less likely but needs to be tested (below) and a suitable
## replace function inserted.
.simpleExtReplace <- function(from, to, value){
    for(what in .InhSlotNames(to))
        slot(from, what) <- slot(value, what)
    from
}
## slot names for inheritance (to be used in replace methods).  Extends slots to implicit
## .Data for basic classes.
.InhSlotNames <- function(Class) {
   ClassDef <- getClass(Class)
    value <- names(ClassDef@slots)
    if(length(value)==0 && (Class %in% .BasicClasses || extends(ClassDef, "vector")))
        ## No slots, but extends "vector" => usually a basic class; treat as data part
        value <- ".Data"
   value
}
.dataPartReplace <- list(f1 = function(from, to, value){
    from@.Data <- value
    from
},

f2 = function(from, to, value){
    from@.Data <- as(value, THISCLASS, strict = FALSE)
    from
},

## and a version of dataPartReplace w/o the unused `to' argument
f2args = function(from, value) {
    from@.Data <- value
    from
})

S3Part <- function(object, strictS3 = FALSE, S3Class) {
    if(!isS4(object))
      return(object)
    classDef <- getClass(class(object))
    oldClassCase <- extends(classDef, "oldClass")
    defltS3Class <- missing(S3Class)
    if(oldClassCase) {
        if(defltS3Class)
            S3Class <- .S3Class(object)
        keepSlots <- slotNames(S3Class[[1L]])
     }
    else {
        if(all(is.na(match(extends(classDef), .BasicClasses))))
          stop(gettextf("S3Part() is only defined for classes set up by setOldCLass(), basic classes or subclasses of these:  not true of class %s", dQuote(class(object))), domain = NA)
        if(missing(S3Class)) {
            S3Class <- classDef@slots$.Data
            if(is.null(S3Class)) # is this an error?
              S3Class <- typeof(object)
            keepSlots <- character()
        }
        else
          keepSlots <- slotNames(S3Class[[1L]])
    }
    if(!(defltS3Class || extends(classDef, S3Class)))
      stop(gettextf("the 'S3Class' argument must be a superclass of %s:  not true of class %s", dQuote(class(object)), dQuote(S3Class)), domain = NA)
    if(strictS3)
      keepSlots <- keepSlots[is.na(match(keepSlots, ".S3Class"))]
    deleteSlots = slotNames(classDef)
    deleteSlots <- deleteSlots[is.na(match(deleteSlots,keepSlots))]
    for(slot in deleteSlots)
      attr(object, slot) <- NULL
    if(strictS3) {
        object <- .notS4(object)
        class(object) <- S3Class
    }
    else
      class(object) <- S3Class[[1L]]
    object
}

"S3Part<-" <- function(object, strictS3 = FALSE, needClass = .S3Class(object) , value) {
    S3Class <- .S3Class(value)
    def <- getClassDef(S3Class[[1L]])
    if(is.null(def) || !extends(def, needClass[[1L]]))
      stop(gettextf("replacement value must extend class %s, got %s", dQuote(needClass), dQuote(S3Class[[1L]])), domain = NA)
    slots <- slotNames(class(object))
    if(!strictS3) {
        fromValue <- names(attributes(value))
        slots <- slots[is.na(match(slots, fromValue))]
    }
    slots <- c("class", slots)  # always preserve class(object)
    for(slot in slots)
      attr(value, slot) <- attr(object, slot)
    if(extends(def, "oldClass"))
      attr(value, ".S3Class") <- S3Class
    if(isS4(object))
      value <- .asS4(value)
    value
}

## templates for replacement methods for S3 classes in classes that extend oldClass
.S3replace <-
    list(e1 =
         quote( {
             S3Part(from, needClass = NEED) <- value
             from
         }),
         e2 = quote( {
             if(is(value, CLASS)) {
                 S3Part(from,  needClass = NEED) <- value
                 from
             }
             else
                 stop(gettextf("replacement value must be of class %s, got one of class %s",
                               dQuote(CLASS),
                               dQuote(class(value)[[1L]])))

         })
         )

.S3coerce <- function(from, to) {
    S3Part(from)
}

.ErrorReplace <- function(from, to, value)
    stop(gettextf("no 'replace' method was defined for 'as(x, \"%s\") <- value' for class %s",
                  to, dQuote(class(from))), domain = NA)

.objectSlotNames <- function(object) {
    ## a quick version that makes no attempt to check the class definition
    value <- names(attributes(object))
    if(is.null(value)) ## not possible with methods package?
        character()
    else
        value[-match("class", value, 0L)]
}

makeExtends <- function(Class, to,
                        coerce = NULL, test = NULL, replace = NULL,
                        by = character(), package,
                        slots = getSlots(classDef1),
                        classDef1 = getClass(Class), classDef2) {
    ## test for datapart class:  must be the data part class, except
    ## that extensions within the basic classes are allowed (numeric, integer)
    dataEquiv <- function(cl1, cl2) {
        .identC(cl1, cl2) ||
          (extends(cl1, cl2) && !any(is.na(match(c(cl1, cl2), .BasicClasses))))
    }
    packageEnv <- .requirePackage(package)
    class1Defined <- missing(slots) # only at this time can we construct methods
    simple <- is.null(coerce) && is.null(test) && is.null(replace) && (length(by)==0)
    distance <- 1
    ##FIX ME:  when by is supplied, should use the existing extension information
    ## to compute distance
    dataPartClass <- elNamed(slots, ".Data")
    dataPart <- FALSE
    if(simple && !is.null(dataPartClass)) {
        if(!(is.null(getClassDef(dataPartClass)) || is.null(getClassDef(to)))) {
            ## note that dataPart, to are looked up in the methods package & parents,
            ## because the default in getClassDef is the topenv of the caller (this fun.):
            ## Assertion is that only these classes are allowed as data slots
            dataPart <- dataEquiv(dataPartClass, to)
        }
    }
    if(is.null(coerce)) {
        coerce <- .simpleExtCoerce
        if(isXS3Class(classDef2)) {
            allNames <- names(slots)
            body(coerce, envir = packageEnv) <-
                substitute({
                    if(strict) S3Part(from, S3Class = S3CLASS)
                    else from
                }, list(S3CLASS =  to))
        }
        else if(!isVirtualClass(classDef2))
            body(coerce, envir = packageEnv) <-
                 .simpleCoerceExpr(Class, to, names(slots), classDef2)
    }
    else if(is(coerce, "function")) {
        ## we allow definitions with and without the `strict' argument
        ## but create a  function that can be called with the argument
        if(length(formals(coerce)) == 1) {
            coerce <- .ChangeFormals(coerce, .simpleIsCoerce, "'coerce' argument to setIs ")
            tmp <- .simpleExtCoerce
            body(tmp, envir = environment(coerce)) <- body(coerce)
            coerce <- tmp
        }
        else
            coerce <- .ChangeFormals(coerce, .simpleExtCoerce, "'coerce' argument to setIs ")

    }
    else stop(gettextf("the 'coerce' argument to 'setIs' should be a function of one argument, got an object of class %s",
                       dQuote(class(coerce))), domain = NA)
    if(is.null(test)) {
        test <- .simpleExtTest
        extClass <- "SClassExtension"
    }
    else {
        test <- .ChangeFormals(test, .simpleExtTest, "'test' argument to setIs ")
        extClass <- "conditionalExtension"
    }
    if(is.null(replace)) {
        if(dataPart) {
            extn <- elNamed(classDef2@contains, dataPartClass)
            if(is(extn, "SClassExtension"))
                easy <- extn@simple
            else
                easy <- FALSE
            if(easy)
                replace <- .dataPartReplace$f1
            else {
                replace <- .dataPartReplace$f2
                bdy <- body(replace)
                body(replace, envir = environment(replace)) <-
                    substituteDirect(bdy, list(THISCLASS = dataPartClass))
            }
        }
        else if(simple) {
            replace <- .simpleExtReplace
            if(isXS3Class(classDef2)) {  # replace the S3 part & slots in class to
                S3Class <- attr(classDef2@prototype, ".S3Class")
                if(is.null(S3Class)) # the setOldClass case ?
                  S3Class <- to
                body(replace, envir = packageEnv) <-
                  quote({
                      S3Part(from) <- value
                      from
                  })
            }
            else if(isVirtualClass(classDef2)) {  # a simple is to a virtual class => a union
                body(replace, envir = packageEnv) <-
                    substitute({
                        if(!is(value, TO))
                            stop(gettextf("the computation: 'as(object,\"%s\") <- value' is valid when object has class %s only if 'is(value, \"%s\")' is TRUE ('class(value)' was %s)\n",
                                 TO, dQuote(FROM), TO, dQuote(class(value))), domain = NA)
                        value
                    }, list(FROM = Class, TO = to))
            }
            else if(class1Defined && length(slots) == 0) {
                ## check for the classes having the same representation
                ## (including the case of no slots)
                ext <- getAllSuperClasses(classDef1, TRUE)
                toSlots <- classDef2@slots
                sameSlots <- TRUE
                for(eclass in ext) {
                    ## does any superclass other than "to" have slots?
                    if(.identC(eclass, to))
                        next
                    edef <- getClassDef(eclass, where = packageEnv)
                    if(!is.null(edef) && length(edef@slots) > 0) {
                        sameSlots <- FALSE
                        break
                    }
                }
                if(sameSlots)
                    body(replace, envir = packageEnv) <-
                        substitute({class(value) <- FROM; value}, list(FROM = Class))
                else if(length(toSlots) == 0) # seems replacement not defined in this case?
                    replace <- .ErrorReplace
            }
            else
                body(replace, envir = packageEnv) <-
                    .simpleReplaceExpr(classDef2)
        }
        else
            replace <- .ErrorReplace
        if(identical(replace, .ErrorReplace))
            warning(gettextf("there is no automatic definition for 'as(object, \"%s\") <- value' when object has class %s and no 'replace' argument was supplied; replacement will be an error",
                             to, dQuote(Class)), domain = NA)
    }
    else if(is(replace, "function")) {
        ## turn function of two or three arguments into correct 3-arg form
        if(length(formals(replace)) == 2) {
            replace <- .ChangeFormals(replace, .dataPartReplace$f2args, "'replace' argument to setIs ")
            tmp  <- .ErrorReplace
            body(tmp, envir = environment(replace)) <- body(replace)
            replace <- tmp
        }
        else
            replace <- .ChangeFormals(replace, .ErrorReplace, "'replace' argument to setIs ")
    }
    else
        stop(gettextf("the 'replace' argument to setIs() should be a function of 2 or 3 arguments, got an object of class %s",
                      dQuote(class(replace))), domain = NA)

    new(extClass, subClass = Class, superClass = to, package = package,
	coerce = coerce, test = test, replace = replace, simple = simple,
	by = by, dataPart = dataPart, distance = distance)
}

.findAll <- function(what, where = topenv(parent.frame())) {
    ## search in envir. & parents thereof
    ## For namespaces, this follows R's soft namespace policy
    ## by not stopping when it reaches the basenamespace
    ## The code used to do so and then had a kludge for looking
    ## in the methods namespace.  But that failed anyway on
    ## non-namespace (package) environments and was inconsistent
    ## with the normal R lookup with namespace environments.
    value <- list()
    if(is.environment(where)) {
        if(isNamespace(where)) repeat {
            if(exists(what, where, inherits = FALSE))
                value <- c(value, list(where))
            if(identical(where, emptyenv()))
                break
            where <- parent.env(where)
        }
        else {  # typically, a package environment: look here, then in the search list
            if(exists(what, where, inherits = FALSE))
                value <- c(value, list(where))
            for(i in seq_along(search())) {
                if(exists(what, i, inherits = FALSE)) {
                    evi <- as.environment(i)
                    addMe<- TRUE
                    for(other in value)
                        if(identical(other, evi)) {
                            addMe <- FALSE
                            break
                        }
                    if(addMe)
                        value <- c(value, list(evi))
                }
            }
        }
    }
    else
        for(i in where) {
            if(exists(what, i, inherits = FALSE))
                value <- c(value, list(i))
        }
    value
}

.S4inherits <- function(x, what, which) {
    superClasses <- extends(getClass(class(x)))
    if(which)
       match(what, superClasses, 0L)
    else
      what %in% superClasses
}

## find the S3 classes or their extensions in the indirect superclasses
## and give them the correct coerce and replacement methods
.S3Extends <- function(ClassDef, exts, where) {
    superClasses <- names(exts)
    S3Class <- attr(ClassDef@prototype, ".S3Class")
    need <- S3Class[[1L]]
    for(i in seq_along(exts)) {
        exti <- exts[[i]]
        if(exti@distance == 1)
            next # asserted that this was done by makeExtends
        what <- superClasses[[i]]
        whatDef <- getClassDef(what, where)
        if(is.null(whatDef) # but shouldn't happen,
           || !isXS3Class(whatDef))
            next
        coerce <- exti@coerce
        body(coerce, environment(coerce))<- body(.S3coerce)
        exti@coerce <- coerce
        replace <- exti@replace
        pos <- match(what, S3Class, 0L)
        if(pos > 1) # not the complete S3 class, probably an error
          body(replace, environment(replace)) <-
            substituteDirect(.S3replace$e2, list(CLASS = what, NEED = need))
        else
          body(replace, environment(replace))  <-
            substituteDirect(.S3replace$e1, list(NEED = need))
        exti@replace <- replace
        exts[[i]] <- exti
    }
    exts
}

#  File src/library/methods/R/ClassUnion.R
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

.InitClassUnion <- function(where) {
    setClass("ClassUnionRepresentation",  "classRepresentation",
             validity =function(object) {
                 if(identical(object@virtual, TRUE) && length(object@slots)==0 &&
                    is.null(object@prototype))
                     TRUE
                 else
                     "Class must be an empty virtual class with NULL prototype"
             }, where = where)
    ## some classes in methods package are unions--now they can be officially
    setClassUnion("OptionalFunction", c("function", "NULL"), where)
    setClassUnion("PossibleMethod", c("function", "MethodDefinition"), where)
    clList <- c("ClassUnionRepresentation", "OptionalFunction",
                "PossibleMethod")
    assign(".SealedClasses", c(get(".SealedClasses", where), clList), where)
}

setClassUnion <- function(name, members = character(), where = topenv(parent.frame())) {
    if(length(members)>0) {
        membersDefined <- sapply(members, isClass, where = as.environment(where))
        if(!all(membersDefined))
            stop(gettextf("the member classes must be defined: not true of %s",
                          paste(.dQ(as(members[!membersDefined], "character")), collapse=", ")), domain = NA)
    }
    def <- new("ClassUnionRepresentation",
               makeClassRepresentation(name, package = getPackageName(where), where = where))
    prev <- getClassDef(name, where = where)
    value <- setClass(name, def, where = where)
    failed <- character()
    ## the prototype of the union will be from the first non-virtual
    ## subclass, except that we prefer NULL if "NULL" is a subclass
    hasNull <- match("NULL", members, 0)
    if(hasNull)
        members <- c("NULL", members[-hasNull])
    for(what in members) {
        if(is(try(setIs(what, name, where = where)), "try-error")) {
            if(!is.character(what))
                what <- getClass(what, TRUE, where)@className
            failed <- c(failed, what)
        }
    }
    if(length(failed)>0) {
        if(is.null(prev))
            try(removeClass(name, where = where))
        else
            try(setClass(name, prev, where = where))
        stop(gettextf("unable to create union class:  could not set members %s",
                      paste(.dQ(failed), collapse=", ")), domain = NA)
    }
    invisible(value)
}

isClassUnion <- function(Class) {
    ## test the class DEFINITION for representing a union
    if(is.character(Class))
        Class <- getClass(Class, TRUE) # the real def. or a dummy
    extends(class(Class), "ClassUnionRepresentation")
}
#  File src/library/methods/R/Defunct.R
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

## old access functions that do nothing but get a slot
## likely to be made defunct & disappear.

## The above comment is there since May 2003, R 1.8.0
## similarly in the help file ../man/oldGet.Rd
## Finally (2008-02) made deprecated formally for 2.7.0, defunct in 2.8.0

getAccess <- function (ClassDef) .Defunct()

getClassName <- function (ClassDef) .Defunct()

getClassPackage <- function (ClassDef) .Defunct()

getExtends <- function (ClassDef) .Defunct()

getProperties <- function (ClassDef)  .Defunct()

getPrototype <- function (ClassDef) .Defunct()

getSubclasses <- function (ClassDef) .Defunct()

getVirtual <- function (ClassDef) .Defunct()

## Deprecated in 2.7.0, defunct in 2.8.0
getAllMethods <- function(f, fdef, where = topenv(parent.frame())) .Defunct()

mlistMetaName <- function(name = "", package = "") .Defunct()

removeMethodsObject <- function(f, where = topenv(parent.frame())) .Defunct()

seemsS4Object <- function(object) .Defunct("isS4")

## Deprecated in 2.8.0, defunct in 2.9.0

allGenerics <- function(...)
    ## this is used nowhere, and we already have too many functions
    .Defunct("getGenerics")

#  File src/library/methods/R/Methods.R
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


setGeneric <-
  ## Define `name' to be a generic  function, for which methods will be defined.
  ##
  ## If there is already a non-generic function of this name, it will be used
  ## to define the generic unless `def' is supplied, and the current
  ## function will become the default method for the generic.
  ##
  ## If `def' is supplied, this defines the generic function.  The
  ## default method for a new generic will usually be an existing
  ## non-generic.  See the .Rd page
  ##
    function(name, def = NULL, group = list(), valueClass = character(),
             where = topenv(parent.frame()),
             package = NULL, signature = NULL,
             useAsDefault = NULL, genericFunction = NULL,
             simpleInheritanceOnly = NULL)
{
    if(is.character(.isSingleName(name)))
        stop(gettextf("invalid argument 'name': %s",
                      .isSingleName(name)), domain = NA)
    if(exists(name, "package:base") &&
       is.primitive(get(name, "package:base"))) { # primitives

        name <- switch(name, "as.double" = "as.numeric", name)
        fdef <- getGeneric(name) # will fail if this can't have methods
        if(nargs() <= 1) {
            ## generics for primitives are global, so can & must always be cached
            .cacheGeneric(name, fdef)
            return(name)
        }
        ## you can only conflict with a primitive if you supply
        ## useAsDefault to signal you really mean a different function
        if(!is.function(useAsDefault) && !identical(useAsDefault, FALSE)) {
            msg <- gettextf("%s is a primitive function;  methods can be defined, but the generic function is implicit, and cannot be changed.", sQuote(name))
            stop(msg, domain = NA)
        }
    }
    simpleCall <- { nargs() < 2 ||
		    all(missing(def), missing(group), missing(valueClass),
			missing(package), missing(signature), missing(useAsDefault),
			missing(genericFunction), missing(simpleInheritanceOnly)) }
    stdGenericBody <- substitute(standardGeneric(NAME), list(NAME = name))
    ## get the current function which may already be a generic
    fdef <-
	if(is.null(package))
	    getFunction(name, mustFind = FALSE, where = where)
	else {
	    ev <- .NamespaceOrPackage(package)
	    if(simpleCall)
		implicitGeneric(name, ev) # generic or NULL
	    else
		getFunction(name, mustFind = FALSE, where = ev)
	}
    if(simpleCall) {
        if(is(fdef, "genericFunction"))
          return(.GenericAssign(name, fdef, where))
    }
    if(is.null(fdef)) {
        if(isNamespace(where))
            fdef <- .getFromStandardPackages(name)
        else
            fdef <- getFunction(name, mustFind = FALSE)
    }
    if(is.null(fdef) && is.function(useAsDefault))
        fdef <- useAsDefault
    ## Use the previous function definition to get the default
    ## and to set the package if not supplied.
    doUncache <- FALSE
    if(is.object(fdef) && is(fdef, "genericFunction")) {
        doUncache <- TRUE
        oldDef <- fdef
        prevDefault <- finalDefaultMethod(fdef@default)
        if(is.null(package))
            package <- fdef@package
    }
    else if(is.function(fdef)) {
        prevDefault <- fdef
        if(is.primitive(fdef)) package <- "base"
        if(is.null(package))
            package <- getPackageName(environment(fdef))
    }
    else
        prevDefault <- NULL

    if(is.primitive(fdef)) ## get the pre-defined version
        fdef <- getGeneric(name, where = where)
    else if(is.function(fdef))
        body(fdef, envir = as.environment(where)) <- stdGenericBody
    if(!is.null(def)) {
        if(is.primitive(def) || !is.function(def))
            stop(gettextf("if the 'def' argument is supplied, it must be a function that calls standardGeneric(\"%s\") or is the default",
                          name), domain = NA)
        nonstandardCase <- .NonstandardGenericTest(body(def), name, stdGenericBody)
        if(is.na(nonstandardCase)) {
            if(is.null(useAsDefault)) {# take this as the default
                useAsDefault <- def
            }
            body(def, envir = as.environment(where)) <- stdGenericBody
            nonstandardCase <- FALSE
        }
        fdef <- def
        if(is.null(genericFunction) && nonstandardCase)
            genericFunction <- new("nonstandardGenericFunction") # force this class for fdef
    }
    thisPackage <- getPackageName(where)
    if(is.null(package) || !nzchar(package))
        ## either no previous def'n or failed to find its package name
        package <- thisPackage
    if(is.null(fdef))
        stop(gettextf("must supply a function skeleton for %s, explicitly or via an existing function", sQuote(name)), domain = NA)
    ensureGeneric.fdef <- function(sig = signature) {
        if(!(is.object(fdef) && is(fdef, "genericFunction"))) {
            fdeflt <-
                if(is.function(useAsDefault)) useAsDefault
                else if(identical(useAsDefault, FALSE)) NULL
                else if(is.function(prevDefault) &&
                        !identical(formalArgs(prevDefault), formalArgs(fdef)) &&
                        !is.primitive(prevDefault))
                    NULL
                else prevDefault
            if(is.function(fdeflt))
                fdeflt <- .derivedDefaultMethod(fdeflt)
            fdef <<-
                makeGeneric(name, fdef, fdeflt, group=group, valueClass=valueClass,
                            package = package, signature = sig,
                            genericFunction = genericFunction,
                            simpleInheritanceOnly = simpleInheritanceOnly)
        }
    }
    if(identical(package, thisPackage)) {
        ensureGeneric.fdef()
    } else {
        ## setting a generic for a function in another package.
        ## In this case, the generic definition must agree with the implicit
        ## generic for the given function and package
        implicit <- implicitGeneric(name, .NamespaceOrPackage(package))
        if(is.null(implicit)) { # New function, go ahead
            ensureGeneric.fdef()
        }
        else {
	    ## possibly take the signature from the *implicit* generic:
	    ensureGeneric.fdef(if(is.null(signature) && is.null(def))
			       implicit@signature else signature)
	    cmp <- .identicalGeneric(fdef, implicit,
				     allow.extra.dots =
				     !nzchar(Sys.getenv("R_SETGENERIC_PICKY_DOTS")))
            if(identical(cmp, TRUE)) {
                fdef <- implicit
            }  # go ahead silently
            else if(is.function(implicit)) {
                thisPName <- if(identical(thisPackage, ".GlobalEnv"))
                    "the global environment" else paste("package", sQuote(thisPackage))
                ## choose the implicit unless an explicit def was given
                if(is.null(def) && is.null(signature)) {
                    message(gettextf(
                       "Creating a generic function for %s from %s in %s\n    (from the saved implicit definition)",
                                     sQuote(name), sQuote(package),
                                     thisPName), domain = NA)
                    fdef <- implicit
                }
                else {
                    message(gettextf(
                         "Creating a new generic function for %s in %s",
                                     sQuote(name), thisPName),
                        domain = NA)
                    fdef@package <- attr(fdef@generic, "package") <- thisPackage
                }
            }
            else { # generic prohibited
                warning(gettextf(
			"no generic version of %s on package %s is allowed;\n   a new generic will be assigned for %s",
                                 sQuote(name), sQuote(package),
                                 thisPName),
                        domain = NA)
                fdef@package <- attr(fdef@generic, "package") <- thisPackage
            }
        }
    }
    if(identical(fdef@signature, "..."))
	fdef <- .dotsGeneric(fdef)
    if(doUncache)
	.uncacheGeneric(name, oldDef)
    groups <- fdef@group
    for(group in groups) { # add as member of group generic(s) if not there
        gdef <- getGeneric(group)
        if(is(gdef, "groupGenericFunction") &&
           is.na(match(fdef@generic, as.character(gdef@groupMembers)))) {
            gwhere <- .genEnv(group, where)
            gdef@groupMembers <- c(gdef@groupMembers, list(fdef@generic))
            assign(group, gdef, gwhere)
        }
    }
    .GenericAssign(name, fdef, where)
}

.GenericAssign <- function(name, fdef, where) {
    assign(name, fdef, where)
    .cacheGeneric(name, fdef)
    methods <- fdef@default # empty or containing the default
    assignMethodsMetaData(name, methods, fdef, where, finalDefaultMethod(fdef@default))
    .assignMethodsTableMetaData(name, fdef, where)
    name
}

## Mimic the search for a function in the standard search() list for packages
## with namespace, to be consistent with the evaluator's search for objects
.standardPackageNamespaces <- new.env()
.standardPackages <- c("stats", "graphics", "grDevices", "utils", "datasets", "methods")
.getFromStandardPackages <- function(name) {
    where <- objects(.standardPackageNamespaces)
    if(!length(where)) { # initialize the table of namespaces
        for(pkg in .standardPackages) {
            ## the tryCatch nonsense is needed because this code gets called
            ## while installing methods (why?) and throws an error on that namespace
            ns <- tryCatch(loadNamespace(pkg), error = function(e) new.env())
            assign(pkg, ns, envir = .standardPackageNamespaces)
        }
        where <- .standardPackages
    }
    ## search
    for(pkg in where) {
        ns <- get(pkg, envir = .standardPackageNamespaces)
        if(exists(name, envir = ns, inherits = FALSE)) {
            obj <- get(name, envir = ns)
            if(is.function(obj))
                return(obj)
        }
    }
    return(NULL)
}

##
## make a generic function object corresponding to the given function name.
##

isGeneric <-
  ## Is there a function named `f', and if so, is it a generic?
  ##
  ## If the `fdef' argument is supplied, take this as the definition of the
  ## generic, and test whether it is really a generic, with `f' as the name of
  ## the generic.  (This argument is not available in S-Plus.)
  function(f, where = topenv(parent.frame()), fdef = NULL, getName = FALSE)
{
    if(is.null(fdef) && missing(where)) {
         fdef <- .getGenericFromCache(f, where)
        ## a successful search will usually end here w/o other tests
         if(!is.null(fdef))
           return(if(getName) fdef@generic else TRUE)
     }
    if(is.null(fdef))
        fdef <- getFunction(f, where=where, mustFind = FALSE)
    if(is.null(fdef))
      return(FALSE)
    ## check primitives. These are never found as explicit generic functions.
    if(is.primitive(fdef)) {
        if(is.character(f) && f %in% "as.double") f <- "as.numeric"
        ## the definition of isGeneric() for a primitive is that methods are defined
        ## (other than the default primitive)
        gen <- genericForPrimitive(f, mustFind = FALSE)
        return(is.function(gen) && length(objects(.getMethodsTable(gen), all.names=TRUE)) > 1L)
    }
    if(!is(fdef, "genericFunction"))
        return(FALSE)
    gen <- fdef@generic # the name with package attribute
    if(missing(f) || .identC(gen, f)) {
	if(getName)
	    gen
	else
	    TRUE
    }
    else {
        warning(gettextf("function %s appears to be a generic function, but with generic name %s",
                         sQuote(f), sQuote(gen)),
                domain = NA)
        FALSE
    }
}

removeGeneric <-
  ## Remove the generic function of this name, specifically the first version
  ## encountered from environment where
  ##
    function(f, where = topenv(parent.frame()))
{
    ev <- fdef <- NULL
    allEv <- findFunction(f, where = where)
    for(maybeEv in allEv) {
        fdef <- get(f, maybeEv)
        if(is(fdef, "genericFunction")) {
            ev <- maybeEv
            break
        }
    }
    found <- is(fdef, "genericFunction")
    if(found) {
         .removeMethodsMetaTable(fdef, where)
         oldMetaName <- methodsPackageMetaName("M",fdef@generic, fdef@package)
         if(exists(oldMetaName, where, inherits = FALSE))
           rm(list = oldMetaName, pos = where)
        .uncacheGeneric(f, fdef)
        rm(list = fdef@generic, pos = where)
    }
    else {
        if(!is.character(f))
            f <- deparse(f)
        warning(gettextf("generic function %s not found for removal",
                         sQuote(f)),
                domain = NA)
    }
    return(found)
}

getMethods <-
    ## The list of methods for the specified generic.  If the function is not
    ## a generic function, returns NULL.
    ## The `f' argument can be either the character string name of the generic
    ## or the object itself.
    ##
    ## The `where' argument optionally says where to look for the function, if
    ## `f' is given as the name.
    ## This function returns a MethodsList object, no longer used for method dispatch
    ## A better structure for most purposes is the linear methods list returned by findMethods()
    ## There are no plans currently to make getMethods defunct, but it will be less
    ## efficient than findMethods()  both for creating the object and using it.

  ##  The function getMethods continues to
  ## return a methods list object, but now this is the metadata from where,
  ## or is converted from the internal table if where is missing
  ## or Mlists are dummies.

    function(f, where = topenv(parent.frame()), table = FALSE)
{
    if(!table)
      .MlistDeprecated("getMethods", "findMethods")
    nowhere <- missing(where)
    fdef <- getGeneric(f, where = where)
    f <- fdef@generic
    if(!is.null(fdef)) {
        if(table)
          return(getMethodsForDispatch(fdef, TRUE))
        value <-
            (if(nowhere) {
                if(is(fdef, "genericFunction"))
                    .makeMlistFromTable(fdef) # else NULL
            }
            else if(.noMlists())
		.makeMlistFromTable(fdef, where)
            else getMethodsMetaData(f, where = where)
             )

        if(is.null(value)) ## return empty methods list
            new("MethodsList", argument = fdef@default@argument) # but deprecated from 2.11.0
        else
            value
    }
    else
      NULL
}

getMethodsForDispatch <- function(fdef, inherited = FALSE)
{
    .getMethodsTable(fdef, environment(fdef), inherited = inherited)
}

## Some functions used in MethodsListSelect, that must be safe against recursive
## method selection.

.setIfBase <- function(f, fdef, mlist) {
    if(is.null(f))
        FALSE
    else {
        found <- base::exists(f, "package:base")
	if(found) {
	    ## force (default) computation of mlist in MethodsListSelect
	    base::assign(".Methods", envir = base::environment(fdef),
			 base::get(f, "package:base"))
	}
        found
    }
}

##NB used internally in MethodsListSelect.  Must NOT use the standard version
## to prevent recursion
.getMethodsForDispatch <- function(fdef) {
    ev <- base::environment(fdef)
    if(base::exists(".Methods", envir = ev))
        base::get(".Methods", envir = ev)
    ## else NULL
}

.setMethodsForDispatch <- function(f, fdef, mlist) {
    ev <- environment(fdef)
    if(!is(fdef, "genericFunction") ||
       !exists(".Methods", envir = ev, inherits = FALSE))
        stop(sprintf("internal error: did not get a valid generic function object for function %s",
                      sQuote(f)),
             domain = NA)
    assign(".Methods", envir = ev, mlist)
}

cacheMethod <-
  ## cache the given definition in the method metadata for f
  ## Support function:  DON'T USE DIRECTLY (does no checking)
  function(f, sig, def, args = names(sig), fdef, inherited = FALSE) {
      ev <- environment(fdef)

      .cacheMethodInTable(fdef, sig, def,
			  .getMethodsTable(fdef, ev, inherited = inherited))
      ## if this is not an inherited method, update the inherited table as well
      ## TODO:	in this case, should uncache inherited methods, though the callin
      ##  function will normally have done this.
      if(!inherited)
	  .cacheMethodInTable(fdef, sig, def,
			      .getMethodsTable(fdef, ev, inherited = TRUE))
  }

.removeCachedMethod <- function(f, sig, fdef = getGeneric(f))
    cacheMethod(f, sig, NULL, names(sig), fdef)


setMethod <-
    ## Define a method for the specified combination of generic function and signature.
    ## The method is stored in the methods meta-data of the specified database.
    ##
    ## Note that assigning methods anywhere but the global environment (`where==1') will
    ## not have a permanent effect beyond the current R session.
    function(f, signature = character(), definition,
	     where = topenv(parent.frame()), valueClass = NULL,
	     sealed = FALSE)
{
    ## Methods are stored in metadata in database where.  A generic function will be
    ## assigned if there is no current generic, and the function is NOT a primitive.
    ## Primitives are dispatched from the main C code, and an explicit generic NEVER
    ## is assigned for them.
    if(is.function(f) && is(f, "genericFunction")) {
        ## (two-part test to deal with bootstrapping of methods package)
        fdef <- f
        f <- fdef@generic
        gwhere <- .genEnv(f)
    }
    else if(is.function(f)) {
        if(is.primitive(f)) {
            f <- .primname(f)
            fdef <- genericForPrimitive(f)
            gwhere <- .genEnv(f)
        }
        else
            stop("a function for argument 'f' must be a generic function")
    }
    ## slight subtlety:  calling getGeneric vs calling isGeneric
    ## For primitive functions, getGeneric returns the (hidden) generic function,
    ## even if no methods have been defined.  An explicit generic MUST NOT be
    ## for these functions, dispatch is done inside the evaluator.
    else {
        where <- as.environment(where)
        gwhere <- .genEnv(f, where)
        f <- switch(f, "as.double" = "as.numeric", f)
        fdef <- getGeneric(f, where = if(identical(gwhere, baseenv())) where else gwhere)
    }
    if(.lockedForMethods(fdef, where))
        stop(gettextf("the environment %s is locked; cannot assign methods for function %s",
                      sQuote(getPackageName(where)),
                      sQuote(f)),
             domain = NA)
    hasMethods <- !is.null(fdef)
    deflt <- getFunction(f, generic = FALSE, mustFind = FALSE, where = where)
    ## where to insert the methods in generic
    if(identical(gwhere, baseenv())) {
        allWhere <- findFunction(f, where = where)
        generics <- logical(length(allWhere))
        if(length(allWhere)) { # put methods into existing generic
            for(i in seq_along(allWhere)) {
                fi <- get(f, allWhere[[i]])
                geni <- is(fi, "genericFunction")
                generics[[i]] <- geni
                if(!geni && is.null(deflt))
                    deflt <- fi
            }
        }
        if(any(generics)) {
            ## try to add method to the existing generic, but if the corresponding
            ## environment is sealed, must create a new generic in where
            gwhere <- as.environment(allWhere[generics][[1L]])
            if(.lockedForMethods(fdef, gwhere)) {
                if(identical(as.environment(where), gwhere))
                    stop(gettextf("the 'where' environment (%s) is a locked namespace; cannot assign methods there",
                                  getPackageName(where)), domain = NA)
                msg <-
                    gettextf("Copying the generic function %s to environment %s, because the previous version was in a sealed namespace (%s)",
                             sQuote(f),
                             sQuote(getPackageName(where)),
                             sQuote(getPackageName(gwhere)))
                message(strwrap(msg), domain = NA)
                assign(f, fdef, where)
                gwhere <- where
            }
        }
    }
    if(!hasMethods)
        fdef <- deflt
    if(is.null(fdef))
        stop(gettextf("no existing definition for function %s",
                      sQuote(f)),
             domain = NA)
    if(!hasMethods) {
        ## create using the visible non-generic as a pattern and default method
        setGeneric(f, where = where)
        fdef <- getGeneric(f, where = where)
        thisPackage <- getPackageName(where)
        thisPName <- if(identical(thisPackage, ".GlobalEnv"))
            "the global environment" else paste("package", sQuote(thisPackage))
        if(identical(as.character(fdef@package), thisPackage))
          message(gettextf("Creating a generic function from function %s in %s",
                           sQuote(f), thisPName), domain = NA)
        else
          message(gettextf("Creating a generic function for %s from package %s in %s",
                           sQuote(f), sQuote(fdef@package), thisPName),
                  domain = NA)
    }
    else if(identical(gwhere, NA)) {
        ## better be a primitive since getGeneric returned a generic, but none was found
        if(is.null(elNamed(.BasicFunsList, f)))
            stop(sprintf("apparent internal error: a generic function was found for \"%s\", but no corresponding object was found searching from \"%s\"",
                          f, getPackageName(where)), domain = NA)
        if(!isGeneric(f))
            setGeneric(f) # turn on this generic and cache it.
    }
    if(isSealedMethod(f, signature, fdef, where=where))
        stop(gettextf("the method for function %s and signature %s is sealed and cannot be re-defined",
                      sQuote(f),
                      .signatureString(fdef, signature)),
             domain = NA)
    signature <- matchSignature(signature, fdef, where)
    createMethod <- FALSE # TRUE for "closure" only
    switch(typeof(definition),
	   "closure" = {
	       fnames <- formalArgs(fdef)
	       mnames <- formalArgs(definition)
	       if(!identical(mnames, fnames)) {
		   ## fix up arg name for single-argument generics
		   ## useful for e.g. '!'
		   if(length(fnames) == length(mnames) && length(mnames) == 1L) {
		       warning(gettextf("For function %s, signature %s: argument in method definition changed from (%s) to (%s)",
					sQuote(f),
                                        sQuote(signature),
                                        mnames,
                                        fnames),
                               domain = NA, call. = FALSE)
		       formals(definition) <- formals(fdef)
		       ll <- list(as.name(formalArgs(fdef))); names(ll) <- mnames
		       body(definition) <- substituteDirect(body(definition), ll)
		       mnames <- fnames
		   }
		   else {
		       ## omitted arguments (classes) in method => "missing"
		       fullSig <- conformMethod(signature, mnames, fnames, f, fdef, definition)
		       if(!identical(fullSig, signature)) {
			   formals(definition, envir = environment(definition)) <- formals(fdef)
			   signature <- fullSig
		       }
		       ## extra arguments (classes) in method => use "..." to rematch
		       definition <- rematchDefinition(definition, fdef, mnames, fnames, signature)
		   }
	       }
	       definition <- matchDefaults(definition, fdef) # use generic's defaults if none in method
               createMethod <- TRUE
	   },
	   "builtin" = , "special" = {
	       ## the only primitive methods allowed are those equivalent
	       ## to the default, for generics that were primitives before
	       ## and will be dispatched by C code.
	       if(!identical(definition, deflt))
		   stop("primitive functions cannot be methods; they must be enclosed in a regular function")
	   },
	   "NULL" = {

	   },
           stop(gettextf("invalid method definition: expected a function, got an object of class %s",
			 dQuote(class(definition))), domain = NA)
	   )
    fenv <- environment(fdef)
    ## check length against active sig. length, reset if necessary in .addToMetaTable
    nSig <- .getGenericSigLength(fdef, fenv, TRUE)
    signature <- .matchSigLength(signature, fdef, fenv, TRUE)
    margs <- (fdef@signature)[seq_along(signature)]
    if(createMethod) {
        definition <- asMethodDefinition(definition, signature, sealed, fdef)
        definition@generic <- fdef@generic
    }
    is.not.base <- !identical(where, baseenv())
##    do.mlist <- is.not.base && (!.noMlists() || all(signature == "ANY"))
    do.mlist <- is.not.base && !.noMlists()
    if(do.mlist)
        whereMethods <- insertMethod(.getOrMakeMethodsList(f, where, fdef),
                                     signature, margs, definition)
    else
        whereMethods <- NULL
    mtable <- getMethodsForDispatch(fdef)
    if(cacheOnAssign(where)) { # will be FALSE for sourceEnvironment's
        ## cache in both direct and inherited tables
        .cacheMethodInTable(fdef, signature, definition, mtable) #direct
        .cacheMethodInTable(fdef, signature, definition) # inherited, by default
        if(is.not.base)
            .addToMetaTable(fdef, signature, definition, where, nSig)
        resetGeneric(f, fdef, mtable, gwhere, deflt) # Note: gwhere not used by resetGeneric
    }
    ## assigns the methodslist object
    ## and deals with flags for primitives & for updating group members
    assignMethodsMetaData(f, whereMethods, fdef, where, deflt)
    f
}

removeMethod <- function(f, signature = character(), where = topenv(parent.frame())) {
    if(is.function(f)) {
      if(is(f, "genericFunction"))
         { fdef <- f; f <- f@generic}
      else if(is.primitive(f))
        { f <- .primname(f); fdef <- genericForPrimitive(f)}
      else
        stop("function supplied as argument 'f' must be a generic")
    }
    else
      fdef <- getGeneric(f, where = where)
    if(is.null(fdef)) {
        warning(gettextf("no generic function %s found", sQuote(f)),
                domain = NA)
        return(FALSE)
    }
    if(is.null(getMethod(fdef, signature, optional=TRUE))) {
        warning(gettextf("no method found for function %s and signature %s",
                         sQuote(fdef@generic),
                         paste(.dQ(signature), collapse =", ")),
                domain = NA)
        return(FALSE)
    }
    setMethod(f, signature, NULL, where = where)
    TRUE
}

## an extension to removeMethod that resets inherited methods as well
.undefineMethod <- function(f, signature = character(), where = topenv(parent.frame())) {
    fdef <- getGeneric(f, where = where)
    if(is.null(fdef)) {
        warning(gettextf("no generic function %s found",
                         sQuote(f)),
                domain = NA)
        return(FALSE)
    }
    if(!is.null(getMethod(fdef, signature, optional=TRUE)))
      setMethod(f, signature, NULL, where = where)
  }

findMethod <- function(f, signature, where = topenv(parent.frame())) {
    if(is(f, "genericFunction")) {
        fdef <- f
        f <- fdef@generic
    }
    else
      fdef <- getGeneric(f, where = where)
    if(is.null(fdef)) {
        warning(gettextf("no generic function %s found",
                         sQuote(f)),
                domain = NA)
        return(character())
    }
    fM <- .TableMetaName(fdef@generic, fdef@package)
    where <- .findAll(fM, where)
    found <- logical(length(where))
    for(i in seq_along(where)) {
        wherei <- where[[i]]
        table <- get(fM, wherei, inherits=FALSE)
        mi <- .findMethodInTable(signature, table, fdef)
        found[i] <- !is.null(mi)
    }
    value <- where[found]
    ## to conform to the API, try to return a numeric or character vector
    ## if possible
    what <- sapply(value, class)
    if(identical(what, "numeric") || identical(what, "character"))
        unlist(value)
    else
        value
}

getMethod <-
  ## Return the function that is defined as the method for this generic function and signature
  ## (classes to be matched to the arguments of the generic function).
  function(f, signature = character(), where = topenv(parent.frame()), optional = FALSE,
           mlist, fdef )
{
    if(!missing(where)) {
        env <- .NamespaceOrEnvironment(where)
        if(is.null(env))
          stop(gettextf("no environment or package corresponding to argument where=%s",
               deparse(where)), domain = NA)
        where <- env
    }
    if(missing(fdef)) {
        if(missing(where))
          fdef <-  getGeneric(f, FALSE)
        else {
            fdef <-  getGeneric(f, FALSE, where = where)
            if(is.null(fdef))
              fdef <- getGeneric(f, FALSE)
        }
    }
    if(!is(fdef, "genericFunction")) {
	if(optional)
	    return(NULL)
	## else
	if(!is.character(f)) f <- deparse(substitute(f))
	stop(gettextf("no generic function found for '%s'", f), domain = NA)
    }
    if(missing(mlist)) {
	if(missing(where))
	    mlist <- getMethodsForDispatch(fdef)
	else
	    mlist <- .getMethodsTableMetaData(fdef, where, optional)
    }
    if(is.environment(mlist)) {
	signature <- matchSignature(signature, fdef)
	value <- .findMethodInTable(signature, mlist, fdef)
	if(is.null(value) && !optional) {
	    if(!is.character(f)) f <- deparse(substitute(f))
	    stop(gettextf("no method found for function '%s' and signature %s",
			  f, paste(signature, collapse = ", ")))
	}
        return(value)
    }
    else if(is.null(mlist)) return(mlist)

    ## the rest of the code will be executed only if a methods list object is supplied
    ## as an argument.  Should be deleted from 2.8.0
    message("Warning: using defunct methods list search", domain = NA)
    i <- 1
    argNames <- fdef@signature
    signature <- matchSignature(signature, fdef)
    Classes <- signature # a copy just for possible error message
    while(length(signature) && is(mlist, "MethodsList")) {
        if(!identical(argNames[[i]], as.character(mlist@argument)))
            stop(sprintf("apparent inconsistency in the methods for function %s; argument %s in the signature corresponds to %s in the methods list object",
                          sQuote(.genericName(f)),
                          sQuote(argNames[[i]]),
                          sQuote(as.character(mlist@argument))),
                 domain = NA)
        Class <- signature[[1L]]
        signature <- signature[-1L]
        methods <- slot(mlist, "methods")
        mlist <- elNamed(methods, Class)# may be function, MethodsList or NULL
        i <- i + 1
    }
    if(length(signature) == 0L) {
        ## process the implicit remaining "ANY" elements
        if(is(mlist, "MethodsList"))
            mlist <- finalDefaultMethod(mlist)
        if(is(mlist, "function"))
            return(mlist) # the only successful outcome
    }
    if(optional)
        mlist                           ## may be NULL or a MethodsList object
    else {
        ## for friendliness, look for (but don't return!) an S3 method
        if(length(Classes) == 1L && exists(paste(.genericName(f), Classes, sep="."), where))
            stop(sprintf("no S4 method for function %s and signature %s; consider getS3method() if you wanted the S3 method",
                         sQuote(.genericName(f)), Classes),
                 domain = NA)
        if(length(Classes)) {
            length(argNames) <- length(Classes)
            Classes <- paste(argNames," = \"", unlist(Classes),
                             "\"", sep = "", collapse = ", ")
        }
        else
            Classes <- "\"ANY\""
        stop(sprintf("no method defined in methods list object for function %s and signature %s",
                     sQuote(.genericName(f)), Classes),
             domain = NA)
    }
}

dumpMethod <-
  ## Dump the method for this generic function and signature.
  ## The resulting source file will recreate the method.
  function(f, signature=character(), file = defaultDumpName(f, signature),
           where = topenv(parent.frame()),
           def = getMethod(f, signature, where=where, optional = TRUE))
{
    if(!is.function(def))
        def <- getMethod(f, character(), where=where, optional = TRUE)

    ## sink() handling as general as possible -- unbelievably unpretty coding:
    closeit <- TRUE ; isSTDOUT <- FALSE
    if (is.character(file)) {
        if(!(isSTDOUT <- file == "")) ## stdout() -- no sink() needed
            file <- file(file, "w")
    }
    else if (inherits(file, "connection")) {
	if (!isOpen(file)) open(file, "w") else closeit <- FALSE
    } else stop("'file' must be a character string or a connection")
    if(!isSTDOUT){ sink(file); on.exit({sink(); if(closeit) close(file)}) }

    cat("setMethod(\"", f, "\", ", deparse(signature), ",\n", sep="")
    dput(def@.Data)
    cat(")\n", sep="")
    if(!isSTDOUT) { on.exit(); sink(); if(closeit) close(file) }
    invisible(file)
}

dumpMethods <- function(f, file = "", signature = NULL, methods= findMethods(f, where = where),
                        where = topenv(parent.frame()) )
{
    ## The signature argument was used in recursive calls to dumpMethods()
    ## using the old MethodsList objects.  It is not meaningful with
    ## the current listOfMethods class
    if(length(signature) > 0)
        warning("argument 'signature' is not meaningful with the current implementation and is ignored \n(extract a subset of the methods list instead)")

    ## sink() handling as general as possible -- unbelievably unpretty coding:
    closeit <- TRUE ; isSTDOUT <- FALSE
    if (is.character(file)) {
        if(!(isSTDOUT <- file == "")) ## stdout() -- no sink() needed
            file <- file(file, "w")
    }
    else if (inherits(file, "connection")) {
	if (!isOpen(file)) open(file, "w") else closeit <- FALSE
    } else stop("'file' must be a character string or a connection")
    if(!isSTDOUT){ sink(file); on.exit({sink(); if(closeit) close(file)}) }
    sigs <- methods@signatures
    for(i in seq_along(methods))
        dumpMethod(f, sigs[[i]], file = "", def = methods[[i]])
}


selectMethod <-
    ## Returns the method (a function) that R would use to evaluate a call to
    ## generic 'f' with arguments corresponding to the specified signature.
    function(f, signature, optional = FALSE, useInherited = TRUE,
	     mlist = if(!is.null(fdef)) getMethodsForDispatch(fdef),
	     fdef = getGeneric(f, !optional), verbose = FALSE, doCache = FALSE)
{
    if(is.environment(mlist))  {# using methods tables
        fenv <- environment(fdef)
        nsig <- .getGenericSigLength(fdef, fenv, FALSE)
        if(verbose)
            cat("* mlist environment with", length(mlist),"potential methods\n")
        if(length(signature) < nsig)
            signature[(length(signature)+1):nsig] <- "ANY"
        if(identical(fdef@signature, "...")) {
            method <- .selectDotsMethod(signature, mlist,
                 if(useInherited) getMethodsForDispatch(fdef, inherited = TRUE))
            if(is.null(method) && !optional)
              stop(gettextf("no method for %s matches class %s",
                            sQuote("..."), dQuote(signature)),
                   domain = NA)
            return(method)
        }
        method <- .findMethodInTable(signature, mlist, fdef)
	if(is.null(method)) {
	    if(missing(useInherited))
		useInherited <- (is.na(match(signature, "ANY")) & # -> vector
				 if(identical(fdef, coerce))# careful !
				 c(TRUE,FALSE) else TRUE)
	    if(verbose) cat("  no direct match found to signature (",
			    paste(signature, collapse=", "),")\n", sep="")
	    methods <-
		if(any(useInherited)) {
		    allmethods <- .getMethodsTable(fdef, fenv, check=FALSE,
                                                   inherited=TRUE)
		    ## look in the supplied (usually standard) table
		    .findInheritedMethods(signature, fdef,
					  mtable = allmethods, table = mlist,
					  useInherited = useInherited,
                                          verbose = verbose, doCache = doCache)
		    ##MM: TODO? allow 'excluded' to be passed
		}
		## else list() : just look in the direct table

	    if(length(methods))
		return(methods[[1L]])
	    else if(optional)
		return(NULL)
	    else stop(gettextf("no method found for signature %s",
			       paste(signature, collapse=", ")))
	}
	else
	  return(method)
    }
    else if(is.null(mlist)) {
	if(optional)
	    return(mlist)
	else
	    stop(gettextf("%s has no methods defined",
                          sQuote(f)),
                 domain = NA)
    }
    else ## mlist not an environment nor NULL :
	stop("selectMethod(): mlist is not an environment or NULL :\n",
	     "** should no longer happen!", domain = NA)
}

hasMethod <-
  ## returns `TRUE' if `f' is the name of a generic function with an (explicit or inherited) method for
  ## this signature.
  function(f, signature = character(), where = .genEnv(f, topenv(parent.frame())))
{
    fdef <- getGeneric(f, where = where)
    if(is.null(fdef))
        FALSE
    else
        !is.null(selectMethod(f, signature, optional = TRUE, fdef = fdef))
}

existsMethod <-
  ## returns `TRUE' if `f' is the name of a generic function with an (explicit) method for
  ## this signature.
  function(f, signature = character(), where = topenv(parent.frame()))
{
        if(missing(where))
          method <- getMethod(f, signature,  optional = TRUE)
        else
          method <- getMethod(f, signature, where = where, optional = TRUE)
        !is.null(method)
}

signature <-
  ## A named list of classes to be matched to arguments of a generic function.
  ## It is recommended to supply signatures to `setMethod' via a call to `signature',
  ## to make clear which arguments are being used to select this method.
  ## It works, however, just to give a vector of character strings, which will
  ## be associated with the formal arguments of the function, in order.  The advantage
  ## of using `signature' is to provide a check on which arguments you meant, as well
  ## as clearer documentation in your method specification.  In addition, `signature'
  ## checks that each of the elements is a single character string.
  function(...)
{
    value <- list(...)
    names <- names(value)
    for(i in seq_along(value)) {
        sigi <- el(value, i)
        if(!is.character(sigi) || length(sigi) != 1L)
            stop(gettextf("bad class specified for element %d (should be a single character string)", i), domain = NA)
    }
      value <- as.character(value)
      names(value) <- names
      value
}

showMethods <-
    ## Show all the methods for the specified function.
    ##
    ## If `where' is supplied, the definition from that database will
    ## be used; otherwise, the current definition is used (which will
    ## include inherited methods that have arisen so far in the
    ## session).
    ##
    ## The output style is different from S-Plus in that it does not
    ## show the database from which the definition comes, but can
    ## optionally include the method definitions, if `includeDefs == TRUE'.
    ##
    function(f = character(), where = topenv(parent.frame()), classes = NULL,
             includeDefs = FALSE, inherited = !includeDefs,
             showEmpty, printTo = stdout(), fdef = getGeneric(f, where = where))
{
    if(missing(showEmpty))
	showEmpty <- !missing(f)
    if(identical(printTo, FALSE))
        con <- textConnection(NULL, "w")
    else
        con <- printTo
    ## must resolve showEmpty in line; using an equivalent default
    ## fails because R resets the "missing()" result for f later on (grumble)
    if(is(f, "function")) {
        fdef <- f ## note that this causes missing(fdef) to be FALSE below
        if(missing(where))
            where <- environment(f)
        f <- deparse(substitute(f))
        if(length(f) > 1L) f <- paste(f, collapse = "; ")
    }
    if(!is(f, "character"))
        stop(gettextf("first argument should be the names of one of more generic functions (got object of class %s)",
                      dQuote(class(f))), domain = NA)
    if(length(f) ==  0L) { ## usually, the default character()
        f <- if(missing(where)) getGenerics() else getGenerics(where)
    }
    if(length(f) == 0L)
	cat(file = con, "no applicable functions\n")
    else if(length(f) > 1L) {
	for(ff in f) { ## recall for each
            ffdef <- getGeneric(ff, where = where)
            if(missing(where)) {
                if(isGeneric(ff))
		    Recall(ff, classes=classes,
			   includeDefs=includeDefs, inherited=inherited,
			   showEmpty=showEmpty, printTo=con, fdef = ffdef)
            }
            else if(isGeneric(ff, where)) {
                Recall(ff, where=where, classes=classes,
                       includeDefs=includeDefs, inherited=inherited,
                       showEmpty=showEmpty, printTo=con, fdef = ffdef)
            }
	}
    }
    else { ## f of length 1 --- the "workhorse" :
        out <- paste0("\nFunction \"", f, "\":\n")
        if(!is(fdef, "genericFunction"))
            cat(file = con, out, "<not an S4 generic function>\n")
        else
            ## maybe no output for showEmpty=FALSE
            .showMethodsTable(fdef, includeDefs, inherited,
                              classes = classes, showEmpty = showEmpty,
                              printTo = con)
    }
    if(identical(printTo, FALSE)) {
        txtOut <- textConnectionValue(con)
        close(con)
        txtOut
    }
    else
        invisible(printTo)
}



removeMethods <-
  ## removes all the methods defined for this generic function.  Returns `TRUE' if
  ## `f' was a generic function, `FALSE' (silently) otherwise.
  ##
  ## If there is a default method, the function will be re-assigned as
  ## a simple function with this definition; otherwise, it will be removed.  The
  ## assignment or removal can be controlled by optional argument `where', which
  ## defaults to the first element of the search list having a function called `f'.
  function(f, where = topenv(parent.frame()), all = missing(where))
{
    ## NOTE:  The following is more delicate than one would like, all because of
    ## methods for primitive functions.  For those, no actual generic function exists,
    ## but isGeneric(f) is TRUE if there are methods.  We have to get the default from
    ## the methods object BEFORE calling removeMethodsObject, in case there are no more
    ## methods left afterwards. AND we can't necessarily use the same default "where"
    ## location for methods object and generic, for the case of primitive functions.
    ## And missing(where) only works in R BEFORE the default is calculated.  Hence
    ## the peculiar order of computations and the explicit use of missing(where).
    fdef <- getGeneric(f, where = where)
    if(!is(fdef, "genericFunction")) {
        warning(gettextf("%s is not an S4 generic function in %s; methods not removed",
                         sQuote(f),
                         sQuote(getPackageName(where))),
                domain = NA)
        return(FALSE)
    }

    methods <- getMethodsForDispatch(fdef)
    default <- getMethod(fdef, "ANY", optional = TRUE)
    fMetaName <- .TableMetaName(fdef@generic, fdef@package)
    oldMetaName <- methodsPackageMetaName("M",fdef@generic, fdef@package)
    allWhere <- .findAll(fMetaName, where)
    if(!all)
        allWhere <- allWhere[1L]
    value <- rep(TRUE, length(allWhere))
    ## cacheGenericsMetaData is called to clear primitive methods if there
    ## are none for this generic on other databases.
    cacheGenericsMetaData(f, fdef, FALSE, where)
    .uncacheGeneric(f, fdef) # in case it gets removed or re-assigned
    doGeneric <- TRUE # modify the function
    for(i in seq_along(allWhere)) {
        db <- as.environment(allWhere[[i]])
        if(environmentIsLocked(db)) {
                warning(gettextf("cannot remove methods for %s in locked environment/package %s",
                                 sQuote(f), sQuote(getPackageName(db))),
                        domain = NA)
                value[[i]] <- FALSE
                next
            }
            if(exists(fMetaName, db, inherits = FALSE)) {
                ## delete these methods from the generic
                theseMethods <- get(fMetaName, db)
                .mergeMethodsTable(fdef, methods, theseMethods, FALSE)
                rm(list = fMetaName, pos = db)
                if(exists(oldMetaName, db, inherits = FALSE))
                  rm(list = oldMetaName, pos = db)
            }
    }
    all <- all && base::all(value) # leave methods on any locked packages
    # now find and reset the generic function
    for(i in seq_along(allWhere)) {
        db <- as.environment(allWhere[[i]])
        if(doGeneric && isGeneric(f, db)) {
            ## restore the original function if one was used as default
            if(all && is(default, "derivedDefaultMethod")) {
                default <- as(default, "function") # strict, removes slots
                rm(list=f, pos = db)
                if(!existsFunction(f, FALSE, db)) {
                    message(gettextf("Restoring default function definition of %s",
                                     sQuote(f)),
                            domain = NA)
                    assign(f, default, db)
                }
                ## else the generic is removed, nongeneric will be found elsewhere
            }
            ## else, leave the generic in place, with methods removed
            ## and inherited methods reset
            else {
                resetGeneric(f, fdef, where = db, deflt = default)
            }
            doGeneric <- FALSE
        }
    }
    any(value)
}


resetGeneric <- function(f, fdef = getGeneric(f, where = where),
			 mlist = getMethodsForDispatch(fdef),
			 where = topenv(parent.frame()),
			 deflt = finalDefaultMethod(mlist))
{
    if(!is(fdef, "genericFunction")) {
            stop(gettextf("error in updating S4 generic function %s; the function definition is not an S4 generic function (class %s)", sQuote(f), dQuote(class(fdef))),
                 domain = NA)
        }
    ## reset inherited methods
    .updateMethodsInTable(fdef, attach = "reset")
    f
}

setReplaceMethod <-
  function(f, ..., where = topenv(parent.frame()))
  setMethod(paste0(f, "<-"), ..., where = where)

setGroupGeneric <-
    ## create a group generic function for this name.
    function(name, def = NULL, group = list(), valueClass = character(),
             knownMembers = list(), package = getPackageName(where), where = topenv(parent.frame()))
{
    if(is.null(def)) {
        def <- getFunction(name, where = where)
        if(isGroup(name, fdef = def)) {
            if(nargs() == 1) {
                message(gettextf("Function %s is already a group generic; no change",
                                 sQuote(name)),
                        domain = NA)
                return(name)
            }
        }
    }
    ## By definition, the body must generate an error.
    body(def, envir = environment(def)) <- substitute(
              stop(MSG, domain = NA),
              list(MSG =
                   gettextf("Function %s is a group generic; do not call it directly",
                            sQuote(name))))
    if(is.character(knownMembers))
        knownMembers <- as.list(knownMembers) # ? or try to find them?
    setGeneric(name, def, group = group, valueClass = valueClass,
               package = package, useAsDefault = FALSE,
               genericFunction =
                 new("groupGenericFunction", def, groupMembers = knownMembers),
               where = where)
    .MakeImplicitGroupMembers(name, knownMembers, where)
    name
}

isGroup <-
  function(f, where = topenv(parent.frame()), fdef = getGeneric(f, where = where))
  {
    is(fdef, "groupGenericFunction")
  }

callGeneric <- function(...)
{
    frame <- sys.parent()
    envir <- parent.frame()
    call <- sys.call(frame)

    ## localArgs == is the evaluation in a method that adds special arguments
    ## to the generic.  If so, look back for the call to generic.  Also expand  "..."
    localArgs <- FALSE
    ## the  lines below this comment do what the previous version
    ## did in the expression fdef <- sys.function(frame)
    if(exists(".Generic", envir = envir, inherits = FALSE))
	fname <- get(".Generic", envir = envir)
    else { # in a local method (special arguments), or	an error
        localArgs <- identical(as.character(call[[1L]]), ".local")
	if(localArgs)
	    call <- sys.call(sys.parent(2))
	fname <- as.character(call[[1L]])
    }
    fdef <- get(fname, envir = envir)

    if(is.primitive(fdef)) {
        if(nargs() == 0)
            stop("'callGeneric' with a primitive needs explicit arguments (no formal args defined)")
        else {
            fname <- as.name(fname)
            call <- substitute(fname(...))
        }
    }
    else {
        env <- environment(fdef)
        if(!exists(".Generic", env, inherits = FALSE))
            stop("'callGeneric' must be called from a generic function or method")
        f <- get(".Generic", env, inherits = FALSE)
        fname <- as.name(f)
        if(nargs() == 0) {
            call[[1L]] <- as.name(fname) # in case called from .local
            ## if ... appears as an arg name, must be a nested callGeneric()
            ##  or callNextMethod?  If so, leave alone so "..." will be evaluated
            if("..." %in% names(call)) {  }
            else {
                ## expand the ... if this is  a locally modified argument list.
                ## This is a somewhat ambiguous case and may not do what the
                ## user expects.  Not clear there is a single solution.  Should we warn?
                call <- match.call(fdef, call, expand.dots = localArgs)
                anames <- names(call)
                matched <- !is.na(match(anames, names(formals(fdef))))
                for(i in seq_along(anames))
                  if(matched[[i]])
                    call[[i]] <- as.name(anames[[i]])
            }
        }
        else {
            call <- substitute(fname(...))
        }
    }
    eval(call, sys.frame(sys.parent()))
}

## This uses 'where' to record the methods namespace: default may not be that
initMethodDispatch <- function(where = topenv(parent.frame()))
    .Call(C_R_initMethodDispatch, as.environment(where))# C-level initialization

### dummy version for booting
isSealedMethod <- function(f, signature, fdef = getGeneric(f, FALSE, where = where),
			   where = topenv(parent.frame())) FALSE

### real version
.isSealedMethod <- function(f, signature, fdef = getGeneric(f, FALSE, where = where),
			   where = topenv(parent.frame()))
{
    ## look for the generic to see if it is a primitive
    fGen <- getFunction(f, TRUE, FALSE, where = where)
    if(!is.primitive(fGen)) {
        mdef <- getMethod(f, signature, optional = TRUE, where = where, fdef = fGen)
        return(is(mdef, "SealedMethodDefinition"))
    }
    ## else, a primitive
    if(is(fdef, "genericFunction"))
        signature <- matchSignature(signature, fdef)
    if(length(signature) == 0L)
        TRUE # default method for primitive
    else if(f %in% .subsetFuns)
        ## primitive dispatch requires some argument to be an S4 object.
        ## This does not quite guarantee an S4 object; e.g., a class union might have only basic types in it.
        !any(is.na(match(signature, .BasicClasses)))
    else {
        sealed <- !is.na(match(signature[[1L]], .BasicClasses))
        if(sealed &&
           (!is.na(match("Ops", c(f, getGroup(f, TRUE))))
            || !is.na(match(f, c("%*%", "crossprod")))))
            ## Ops methods are only sealed if both args are basic classes
            sealed <- sealed && (length(signature) > 1L) &&
                      !is.na(match(signature[[2L]], .BasicClasses))
        sealed
    }
}

.subsetFuns <- c("[", "[[","[<-","[[<-")

.lockedForMethods <- function(fdef, env) {
    ## the env argument is NULL if setMethod is only going to assign into the
    ## table of the generic function, and not to assign methods list object
    if(is.null(env) || !environmentIsLocked(env))
        return(FALSE) #? can binding be locked and envir. not?
    if(!is(fdef, "genericFunction"))
      return(TRUE)
    name <- fdef@generic
    package <- fdef@package
    objs <- c(name, .TableMetaName(name, package))
    for(obj in objs) {
        hasIt <- exists(obj, env, inherits = FALSE)
        ## the method object may be bound, or a new one may be needed
        ## in which case the env. better not be locked
        if((!hasIt || bindingIsLocked(obj, env)))
            return(TRUE)
    }
    FALSE
}

implicitGeneric <- function(...) NULL

## real version, installed after methods package initialized

.implicitGeneric <- function(name, where = topenv(parent.frame()),
                             generic = getGeneric(name, where = where))
### Add the named function to the table of implicit generics in environment where.
###
### If there is a generic function of this name, it is saved to the
### table.  This is the reccomended approach and is required if you
### want the saved generic to include any non-default methods.
###
  {
      if(!nzchar(name))
        stop(gettextf('expected a non-empty character string for argument name'), domain = NA)
      if(!missing(generic) && is(generic, "genericFunction") && !.identC(name, generic@generic))
        stop(gettextf('generic function supplied was not created for %s',
                      sQuote(name)),
             domain = NA)
      createGeneric <- (missing(generic) || !is(generic, "genericFunction")) && !isGeneric(name, where)
      if(createGeneric) {
          fdefault <- getFunction(name, where = where, mustFind = FALSE)
          if(is.null(fdefault))
            return(NULL)  # no implicit generic
          env <- environment(fdefault) # the environment for an implicit generic table
          fdefault <- .derivedDefaultMethod(fdefault)
          if(is.primitive(fdefault)) {
              value <- genericForPrimitive(name)
              if(!missing(generic) && !identical(value, generic))
                  stop(gettextf("%s is a primitive function; its generic form cannot be redefined",
                                sQuote(name)),
                       domain = NA)
              generic <- value
              package <- "base"
          }
          else
              package <- getPackageName(env)
          ## look for a group
          group <-
              .getImplicitGroup(name,
                                if(identical(package,"base"))
                                .methodsNamespace else environment(fdefault))
          if(missing(generic)) {
            generic <- .getImplicitGeneric(name, env, package)
            if(is.null(generic))  { # make a new one
                generic <- makeGeneric(name, fdefault = fdefault, package = package,
                                       group = group)
                .cacheImplicitGeneric(name, generic)
            }
          }
          else {
            generic <- makeGeneric(name, generic, fdefault, package = package,
                                   group = group)
            .cacheImplicitGeneric(name, generic)
        }
      }
      generic
  }

setGenericImplicit <- function(name, where = topenv(parent.frame()), restore = TRUE) {
    if(!isGeneric(name, where)) {
        warning(gettextf("%s is not currently a generic:  define it first to create a non-default implicit form",
                         sQuote(name)),
                domain = NA)
        return(FALSE)
    }
    generic <- getGeneric(name, where = where)
    if(restore)
        removeMethods(name, where, TRUE)
    else
        removeGeneric(name, where)
    .saveToImplicitGenerics(name, generic, where)
}

prohibitGeneric <- function(name, where = topenv(parent.frame()))
### store a definition in the implicit generic table that explicitly prohibits
### a function from being made generic
  {
      .saveToImplicitGenerics(name, FALSE, where)
  }

registerImplicitGenerics <- function(what = .ImplicitGenericsTable(where),
                                     where = topenv(parent.frame()))
{
    if(!is.environment(what))
        stop(gettextf("must provide an environment table; got class %s",
                      dQuote(class(what))), domain = NA)
    objs <- objects(what, all.names = TRUE)
    for(f in objs)
        .cacheImplicitGeneric(f, get(f, envir = what))
    NULL
}


### the metadata name for the implicit generic table
.ImplicitGenericsMetaName <- ".__IG__table" # methodsPackageMetaName("IG", "table")

.ImplicitGenericsTable <- function(where)
  {
### internal utility to add a function to the implicit generic table
      if(!exists(.ImplicitGenericsMetaName, where, inherits = FALSE))
        assign(.ImplicitGenericsMetaName, new.env(TRUE), where)
      get(.ImplicitGenericsMetaName, where)
  }

.saveToImplicitGenerics <- function(name, def, where)
  .cacheGenericTable(name, def, .ImplicitGenericsTable(where))

.getImplicitGeneric <- function(name, where, pkg = "")
{
    value <- .getImplicitGenericFromCache(name, where, pkg)
    if(is.null(value) && exists(.ImplicitGenericsMetaName, where, inherits = FALSE)) {
        tbl <-  get(.ImplicitGenericsMetaName, where)
        value <- .getGenericFromCacheTable(name, where, pkg, tbl)
    }
    value
}

## only called from setGeneric, f1 = supplied, f2 = implicit
.identicalGeneric <- function(f1, f2, allow.extra.dots = FALSE)
{
    gpString <- function(gp) {
	if(length(gp))
	    paste(as.character(gp), collapse = ", ")
	else
	    "<none>"
    }
    if(identical(f2, FALSE))
	return(gettext("original function is prohibited as a generic function"))
    if(!(is.function(f2) && is.function(f1)))
	return(gettext("not both functions!"))
    ## environments will be different
    if(!identical(class(f1), class(f2)))
	return(sprintf("classes: %s, %s",
                       .dQ(class(f1)), .dQ(class(f2))))
    if(!isS4(f1)) return(gettextf("argument %s is not S4",
                                  deparse(substitute(f1))))
    if(!isS4(f2)) return(gettextf("argument %s is not S4",
                                  deparse(substitute(f2))))
    f1d <- f1@.Data
    f2d <- f2@.Data
    ## xtra... <- FALSE
    if(!identical(formals(f1d), formals(f2d))) {
	a1 <- names(formals(f1d)); a2 <- names(formals(f2d))
	if(identical(a1, a2))
	    return(gettext("formal arguments differ (in default values?)"))
	else if(identical(c(a1, "..."), a2) && allow.extra.dots)
            ## silently accept an extra "..."
            { } ## xtra... <- TRUE
	    ## and continue
	else
	    return(gettextf("formal arguments differ: (%s), (%s)",
			    paste(a1, collapse = ", "),
			    paste(a2, collapse = ", ")))
    }
    if(!identical(f1@valueClass, f2@valueClass))
	return(gettextf("value classes differ: %s, %s",
                        .dQ(gpString(f1@valueClass)),
                        .dQ(gpString(f2@valueClass))))
    if(!identical(body(f1d), body(f2d)))
	return("function body differs")
    if(!identical(f1@signature, f2@signature))
	return(gettextf("signatures differ:  (%s), (%s)",
                        paste(f1@signature, collapse = ", "),
                        paste(f2@signature, collapse = ", ")))
    if(!identical(f1@package, f2@package))
	return(gettextf("package slots  differ: %s, %s",
                        .dQ(gpString(f1@package)),
                        .dQ(gpString(f2@package))))
    if(!identical(f1@group, f2@group)) {
	return(gettextf("groups differ: %s, %s",
                        .dQ(gpString(f1@group)),
                        .dQ(gpString(f2@group))))
    }
    if(!identical(as.character(f1@generic), as.character(f2@generic)))
	return(gettextf("generic names differ: %s, %s",
                        .dQ(f1@generic), .dQ(f2@generic)))
    TRUE
}

.ImplicitGroupMetaName <- ".__IGM__table"
.MakeImplicitGroupMembers <- function(group, members, where) {
    if(!exists(.ImplicitGroupMetaName, where, inherits = FALSE))
        assign(.ImplicitGroupMetaName, new.env(TRUE), where)
    tbl <- get(.ImplicitGroupMetaName, where)
    for(what in members)
        assign(what, as.list(group), envir = tbl)
    NULL
}

.getImplicitGroup <- function(name, where) {
    if(exists(.ImplicitGroupMetaName, where, inherits = FALSE)) {
        tbl <- get(.ImplicitGroupMetaName, where)
        if(exists(name, envir = tbl, inherits = FALSE))
            return(get(name, envir = tbl))
    }
    list()
}

findMethods <- function(f, where, classes = character(), inherited = FALSE, package = "") {
    if(is(f, "genericFunction")) {
        fdef <- f
        f <- fdef@generic
    }
    else if(.isSingleString(f)) {
        if(missing(where))
            fdef <- getGeneric(f, package = package)
        else { # the generic may not be in the where= environment
            ##  but we prefer that version if it is
            fdef <- getGeneric(f, where = where, package = package)
            if(is.null(fdef))
                fdef <- getGeneric(f, package = package)
        }
    }
    else if(!is(f, "function"))
        stop(gettextf("argument %s must be a generic function or a single character string; got an object of class %s",
                      sQuote("f"), dQuote(class(f))),
             domain = NA)
    else {
        fdef <- f
        f <- deparse(substitute(f))
    }
    if(!is(fdef, "genericFunction")) {
        warning(gettextf("non-generic function '%s' given to findMethods()", f),
                domain = NA)
        return(list())
    }
    object <- new("listOfMethods", arguments = fdef@signature,
                  generic = fdef) # empty list of methods
    if(missing(where))
      table <- get(if(inherited) ".AllMTable" else ".MTable", envir = environment(fdef))
    else {
        if(!identical(inherited, FALSE))
          stop("only FALSE is meaningful for 'inherited', when 'where' is supplied (got ", inherited, "\"")
        where <- as.environment(where)
        what <- .TableMetaName(f, fdef@package)
        if(exists(what, envir = where, inherits = FALSE))
          table <- get(what, envir = where)
        else
          return(object)
    }
    objNames <- objects(table, all.names = TRUE)
    if(length(classes)) {
        classesPattern <- paste0("#", classes, "#", collapse = "|")
        which <- grep(classesPattern, paste0("#",objNames,"#"))
        objNames <- objNames[which]
    }
    object@.Data <- lapply(objNames, function(x)get(x, envir = table))
    object@names <- objNames
    object@signatures <- strsplit(objNames, "#", fixed = TRUE)
    object
}

findMethodSignatures <- function(..., target = TRUE, methods = findMethods(...))
{
    what <- methods@arguments
    if(target)
      sigs <- methods@signatures
    else {
        anySig <- rep("ANY", length(what))
        ## something of a kludge for the case of some primitive
        ## default methods to get a vector of "ANY" of right length
        for(m in methods)
          if(!is.primitive(m)) {
              length(anySig) <- length(m@defined)
              break
          }
        sigs <- lapply(methods, function(x)
                       if(is.primitive(x)) anySig else as.character(x@defined))
    }
    lens <- unique(sapply(sigs, length))
    if(length(lens) == 0)
        return(matrix(character(), 0, length(methods@arguments)))
    if(length(lens) > 1L) {
        lens <- max(lens)
        anys <- rep("ANY", lens)
        sigs <- lapply(sigs, function(x) {
            if(length(x) < lens) {
              anys[seq_along(x)] <- x
              anys
          } else x
        })
    }
    length(what) <- lens # if not all possible arguments used
    t(matrix(unlist(sigs), nrow = lens, dimnames = list(what, NULL)))
}

hasMethods <- function(f, where, package = "")
{
    fdef <- NULL
    nowhere <- missing(where) # because R resets this if where is assigned
    if(is(f, "genericFunction")) {
        fdef <- f
        f <- fdef@generic
        if(missing(package))
            package <- fdef@package
    }
    else if(!.isSingleString(f))
        stop(gettextf("argument 'f' must be a generic function or %s",
                      .notSingleString(f)), domain = NA)
    else if(missing(package)) {
        package <- packageSlot(f) # maybe a string with package slot
	if(is.null(package)) {
            if(missing(where))
                fdef <- getGeneric(f)
            else { # the generic may not be in this package, but prefer it if so
                fdef <- getGeneric(f, where = where)
                if(is.null(fdef))
                    fdef <- getGeneric(f)
            }
            if(is(fdef, "genericFunction"))
                package <- fdef@package
	    else
		stop(gettextf("'%s' is not a known generic function {and 'package' not specified}",
			      f),
		     domain = NA)
	}
    }
    what <- .TableMetaName(f, package)
    testEv <- function(ev)
      exists(what, envir = ev, inherits = FALSE) &&
    length(objects(get(what, envir = ev), all.names = TRUE))
    if(nowhere) {
        for(i in seq_along(search())) {
            if(testEv(as.environment(i)))
              return(TRUE)
        }
        return(FALSE)
    }
    else
      testEv(as.environment(where))
}
## returns TRUE if the argument is a non-empty character vector of length 1
## otherwise, returns a diagnostic character string reporting the non-conformance
.isSingleName <- function(x) {
    if(!is.character(x))
      return(paste0('required to be a character vector, got an object of class "', class(x)[[1L]], '"'))
    if(length(x) != 1)
      return(paste0("required to be a character vector of length 1, got length ",length(x)))
    if(is.na(x) || !nzchar(x))
      return(paste0('required a non-empty string, got "',x, '"'))
    TRUE
}
#  File src/library/methods/R/MethodsList.R
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

MethodsList <-
  ## Create a MethodsList object out of the arguments.
  ##
  ## Conceptually, this object is a named collection of methods to be
  ## dispatched when the (first) argument in a function call matches the
  ## class corresponding to one of the names.  A final, unnamed element
  ## (i.e., with name `""') corresponds to the default method.
  ##
  ## The elements can be either a function, or another MethodsList.  In
  ## the second case, this list implies dispatching on the second
  ## argument to the function using that list, given a selection of this
  ## element on the first argument.  Thus, method dispatching on an
  ## arbitrary number of arguments is defined.
  ##
  ## MethodsList objects are used primarily to dispatch OOP-style
  ## methods and, in R, to emulate S4-style methods.
  function(.ArgName, ...)
{
    value <- makeMethodsList(list(...))
    if(is.name(.ArgName)){}
    else if(is.character(.ArgName) && length(.ArgName) == 1)
        .ArgName <- as.name(.ArgName)
    else stop("invalid first argument: should be the name of the first argument in the dispatch")
    slot(value, "argument") <- .ArgName
    value
}

makeMethodsList <- function(object, level=1)
{
    mnames <- allNames(object)
    if(.noMlists()) {
        keep <- mnames %in% c("", "ANY")
        mnames <- mnames[keep]
        object <- object[keep]
    }
    value <- new("MethodsList")
    i <- match("", mnames)
    if(!is.na(i)) {
        ## convert to ANY
        el(mnames, i) <- "ANY"
        names(object) <- mnames
    }
    if(anyDuplicated(mnames))
        stop(gettextf("duplicate element names in 'MethodsList' at level %d: %s",
             level, paste("\"", unique(mnames[duplicated(mnames)]), "\"",
                          collapse=", ")), domain = NA)
    for(i in seq_along(object)) {
        eli <- el(object, i)
        if(is(eli, "function")
           || is(eli, "MethodsList")) {}
        else if(is(eli, "list") ||
                is(eli, "named"))
            el(object, i) <- Recall(eli, NULL, level+1)
        else
            stop(gettextf("element %d at level %d (class %s) cannot be interpreted as a function or named list",
                          i, level, dQuote(class(eli))),
                 domain = NA)
    }
    slot(value, "methods") <- object
    value
}

SignatureMethod <-
  ## construct a MethodsList object containing (only) this method, corresponding
  ## to the signature; i.e., such that signature[[1L]] is the match for the first
  ## argument, signature[[2L]] for the second argument, and so on.  The string
  ## "missing" means a match for a missing argument, and "ANY" means use this as the
  ## default setting at this level.
  ##
  ## The first argument is the argument names to be used for dispatch corresponding to
  ## the signatures.
  function(names, signature, definition)
{
    n <- length(signature)
    if(n > length(names))
        stop("arguments 'names' and 'signature' must have the same length")
    if(n == 0)
        return(definition)
    Class <- el(signature,n)
    name <- el(names, n)
    m <- MethodsList(name)
    elNamed(slot(m, "methods"), Class) <- definition
    slot(m, "argument") <- as.name(name)
    SignatureMethod(names[-n], signature[-n], m)
}


insertMethod <-
  ## insert the definition `def' into the MethodsList object, `mlist', corresponding to
  ## the signature, and return the modified MethodsList.
  function(mlist, signature, args, def, cacheOnly = FALSE)
{
    if(.noMlists() && !identical(unique(signature), "ANY"))
      return(mlist)
    ## Checks for assertions about valid calls.
    ## See rev. 1.17 for the code before the assertions added.
    if(identical(args[1L], "...") && !identical(names(signature), "...")) {
        if(identical(signature[[1L]], "ANY"))
           stop(gettextf("inserting method with invalid signature matching argument '...' to class %s",
                         dQuote(signature[[1L]])),
                domain = NA)
        args <- args[-1L]
        signature <- signature[-1L]
        if(length(signature) == 0L)
            return(mlist)
    }
    if(length(signature) == 0L)
        stop("inserting method corresponding to empty signature")
    if(!is(mlist, "MethodsList"))
        stop(gettextf("inserting method into non-methods-list object (class %s)",
                      dQuote(.class1(mlist))),
             domain = NA)
    if(length(args) > 1 && !cacheOnly)
        mlist <- balanceMethodsList(mlist, args)
    Class <- el(signature, 1)
    methods <- if(cacheOnly) mlist@allMethods else mlist@methods
    current <- elNamed(methods, Class)
    if(is(current, "MethodsList")) {
        nextArg <- as.character(current@argument)
        sigArgs <- args
        n <- length(signature)
        length(sigArgs) <- n
        if(is.na(match(nextArg, sigArgs))) {
            n <- match(nextArg, args) - n
            if(is.na(n)) { ## not in args eitiher
                n <- 1
                args <- c(args, nextArg)
            }
            ## make explicit the trailing ANY's needed
            signature <- c(signature, rep("ANY", n))
        }
    }
    if(length(signature) == 1) {
        if(is.null(current)) {
            if(!is.null(def))
                elNamed(methods, Class) <- def
            ## else, no change
        }
        else {
            which <- match(Class, names(methods))
            if(is.null(def))
                ## delete the method
                methods <- methods[-which]
            else
                el(methods, which) <- def
        }
    }
    else { ## recursively merge, initializing current if necessary
        if(is.null(current))
            current <- new("MethodsList", argument = as.name(args[2L]))
        else if(is.function(current))
            current <- new("MethodsList", argument = as.name(args[2L]),
			   methods = list(ANY = current))
        elNamed(methods, Class) <-
            Recall(current, signature[-1L], args[-1L], def, cacheOnly)
    }
    mlist@allMethods <- methods
    if(!cacheOnly)
        mlist@methods <-  methods
    mlist
}


MethodsListSelect <-
  ## select the element of a MethodsList object corresponding to the
  ## actual arguments (as defined in the suppled environment),
  ## and return the object, extended to include that method if necessary.
  ##
  ## Works recursively.  At each level finds an argument name from the current `mlist'
  ## object, and evaluates this argument (if it is not missing), then uses the
  ## class of the result to select an element of `mlist'.  If such an element
  ## exists and is another `MethodsList' object, `MethodsListSelect'  calls itself recursively
  ## to resolve using further arguments.  Matching includes using a default selection or
  ## a method specifically linked to class `"missing"'.  Once a function is found, it
  ## is returned as the value.  If matching fails,  NULL is returned.
    function(f, env,
             mlist = NULL,
             fEnv = if(is(fdef, "genericFunction")) environment(fdef) else baseenv(),
             finalDefault = finalDefaultMethod(mlist),
             evalArgs = TRUE,
             useInherited = TRUE,  ## supplied when evalArgs is FALSE
             fdef = getGeneric(f, where = env), # MUST BE SAFE FROM RECUSIVE METHOD SELECTION
             resetAllowed = TRUE # FALSE when called from selectMethod, .findNextMethod
 )
{
    if(!resetAllowed) # ensure we restore the real methods for this function
	resetMlist <- .getMethodsForDispatch(fdef)
    ## look for call from C dispatch code during another call to MethodsListSelect
    if(is.null(f)) {} # Recall, not from C
    else {
	fMethods <- .getMethodsForDispatch(fdef)
        if(is.null(mlist) || (evalArgs && is.function(fMethods)))
            mlist <- fMethods
    }
    resetNeeded <- .setIfBase(f, fdef, mlist) # quickly protect against recursion -- see Methods.R
    if(resetNeeded) {
        on.exit(.setMethodsForDispatch(f, fdef, mlist))
    }
    if(!is(mlist, "MethodsList")) {
        if(is.function(mlist)) # call to f, inside MethodsListSelect
            {on.exit(); return(mlist)}
        if(is.null(f)) # recursive recall of MethodsListSelect
            stop("invalid method sublist")
        else if(!is.null(mlist)) # NULL => 1st call to genericFunction
            stop(gettextf("%f is not a valid generic function: methods list was an object of class %s",
                          sQuote(f), dQuote(class(mlist))),
                 domain = NA)
    }
    if(!is.logical(useInherited))
        stop(gettextf("%s must be TRUE, FALSE, or a named logical vector of those values; got an object of class %s",
                      sQuote("useInherited"),
                      dQuote(class(useInherited))),
             domain = NA)
    if(identical(mlist, .getMethodsForDispatch(fdef))) {
        resetNeeded <- TRUE
        ## On the initial call:
        ## turn off any further method dispatch on this function, to avoid recursive
        ## loops if f is a function used in MethodsListSelect.
        ## TODO: Using namespaces in the methods package would eliminate the need for this
        .setMethodsForDispatch(f, fdef, finalDefault)
        if(is(mlist, "MethodsList")) {
            on.exit(.setMethodsForDispatch(f, fdef, mlist))
        }
    }
    argName <- slot(mlist, "argument")
    arg <- NULL ## => don't use instance-specific inheritance
    if(evalArgs) {
        ## check for missing argument. NB: S sense, not that of R base missing()
        if(missingArg(argName, env, TRUE))
            thisClass <- "missing"
        else {
            arg <- eval(as.name(argName), env) ## DO use instance-specific inheritance
	    if(missing(arg)) ## S3 weird R code? Bail out!
		return(finalDefault)
            thisClass <- .class1(arg)
        }
    }
    else
        thisClass <- get(as.character(argName), envir = env, inherits = FALSE)
    if(identical(useInherited, TRUE) || identical(useInherited, FALSE))
        thisInherit <- nextUseInherited <- useInherited
    else {
        which <- match(as.character(argName), names(useInherited))
        if(is.na(which)) {
            nextUseInherited <- useInherited
            thisInherit <- TRUE
        }
        else {
            thisInherit <- useInherited[[which]]
            nextUseInherited <- useInherited[-which]
        }
    }
    fromClass <- thisClass ## will mark the class actually providing the method
    allMethods <- mlist@allMethods
    which <- match(thisClass, names(allMethods))
    inherited <- is.na(which)
    selection <- if(inherited) NULL else allMethods[[which]]
    if(!inherited) {
        if(is(selection, "function")) {
            if(is.null(f)) {
              ## An inherited method at the next level up.
              ## only the inherited method should be added
              mlist <- .trimMlist(mlist, fromClass)
            }
            value <- mlist ## no change
          }
        else {
            ## recursive call with NULL function name, to allow search to fail &
            ## to suppress any reset actions.
            method <- Recall(NULL, env, selection, finalDefault = finalDefault,
                   evalArgs = evalArgs, useInherited = nextUseInherited, fdef = fdef,
                             )
            if(is(method, "EmptyMethodsList"))
                value <- method
            else {
                mlist@allMethods[[which]] <- method
                value <- mlist
            }
        }
    }
    if(inherited || is(value, "EmptyMethodsList"))  {
        ## direct selection failed at this level or below
        method <- NULL
        if(thisInherit)  {
            allSelections <- inheritedSubMethodLists(arg, fromClass, mlist, env)
            allClasses <- names(allSelections)
            for(i in seq_along(allSelections)) {
                selection <- allSelections[[i]]
                fromClass <- allClasses[[i]]
                if(is(selection, "function"))
                    method <- selection
                else if(is(selection, "MethodsList")) {
                    ## go on to try matching further arguments
                    method <- Recall(NULL, env, selection, finalDefault = finalDefault,
                                     evalArgs = evalArgs,
                                     useInherited = nextUseInherited, fdef = fdef)
                    if(is(method, "EmptyMethodsList"))
                        selection <- method   ## recursive selection failed
                }
                if(!is(selection, "EmptyMethodsList"))
                    break
            }
        }
        if((is.null(selection) || is(selection, "EmptyMethodsList"))
           && !is.null(f) && !is.null(finalDefault)) {
            ## only use the final default method after exhausting all
            ## other possibilities, at all levels.
            method <- finalDefault
            fromClass <- "ANY"
        }
        if(is.null(method) || is(method, "EmptyMethodsList"))
            value <- emptyMethodsList(mlist, thisClass) ## nothing found
        else {
            method <- MethodAddCoerce(method, argName, thisClass, fromClass)
            value <- .insertCachedMethods(mlist, as.character(argName), thisClass, fromClass,
                                         method)
        }
    }
    if(!is.null(f)) {
        ## top level
        if(is(value, "EmptyMethodsList")) ## selection failed
            value <- NULL
        if(resetNeeded) {
            on.exit() # cancel the restore of the original mlist
            if(resetAllowed) {
                if(is.null(value)) resetMlist <- mlist else resetMlist <- value
            }
            .setMethodsForDispatch(f, fdef, resetMlist)
            if(is.primitive(finalDefault))
                setPrimitiveMethods(f, finalDefault, "set", fdef, resetMlist)
        }

    }
    value
}

emptyMethodsList <- function(mlist, thisClass = "ANY", sublist = list()) {
    sublist[thisClass] <- list(NULL)
    new("EmptyMethodsList", argument = mlist@argument, sublist = sublist)
}

insertMethodInEmptyList <- function(mlist, def) {
    value <- new("MethodsList", argument = mlist@argument)
    sublist <- mlist@sublist
    submethods <- sublist[[1L]]
    if(is.null(submethods))
        sublist[[1L]] <- def
    else
        sublist[[1L]] <- Recall(submethods, def)
    value@allMethods <- sublist
    value
}




finalDefaultMethod <-
  ## Return the default method from the generic (it may be NULL, a method object or a primitive.
  ## this previously searched in a MethodsList object.  Once those are gone, the loop should
  ## be irrelevant except as an error check.
  function(method)
{
    repeat {
        if(is.function(method) #somewhat liberal, but catches both methods and primitives
           || is.null(method))
          break
        value <- NULL
        if(is(method, "MethodsList"))
            method <-  elNamed(slot(method, "methods"), "ANY")
        else
          stop(gettextf("default method must be a method definition, a primitive or NULL: got an object of class %s", dQuote(class(method))),
               domain = NA)
    }
    method
}


inheritedSubMethodLists <-
  ## Utility function to match the object to the elements of a methods list.
  ##
  ## The function looks only for an inherited match, and only among
  ## the methods that are not themselves inherited.  (Inherited methods when found are
  ## stored in the session copy of the methods list, but they themselves should not be
  ## used for finding inherited matches, because an erroneous match could be found depending
  ## on which methods were previously used.  See the detailed discussion of methods.)
  function(object, thisClass, mlist, ev)
{
  methods <- slot(mlist, "methods")## only direct methods
  defaultMethod <- elNamed(methods, "ANY")## maybe NULL
  classes <- names(methods)
  value <- list()
  if(.identC(thisClass, "missing")) {
        ## no superclasses for "missing"
  }
  else {
      ## search in the superclasses, but don't use inherited methods
      ## There are two cases:  if thisClass is formally defined & unsealed, use its
      ## superclasses.  Otherwise, look in the subclasses of those classes for
      ## which methods exist.
      classDef <- getClassDef(thisClass, ev)
      useSuperClasses <- !is.null(classDef) && !classDef@sealed
      if(useSuperClasses) {
          ## for consistency, order the available methods by
          ## the ordering of the superclasses of thisClass
          superClasses <- names(classDef@contains)
          classes <- superClasses[!is.na(match(superClasses, classes))]
          for(which in seq_along(classes)) {
              tryClass <- el(classes, which)
              ## TODO:  There is potential bug here:  If the is relation is conditional,
              ## we should not cache this selection.  Needs another trick in the environment
              ## to FORCE no caching regardless of what happens elsewhere; e.g., storing a
              ## special object in .Class
              if(is.null(object) || is(object, tryClass)) {
                  elNamed(value, tryClass) <- elNamed(methods, tryClass)
              }
          }
      }
      else {
          for(which in seq_along(classes)) {
              tryClass <- el(classes, which)
              tryClassDef <- getClassDef(tryClass, ev)
              if(!is.null(tryClassDef) &&
                 !is.na(match(thisClass, names(tryClassDef@subclasses))))
                  elNamed(value, tryClass) <- el(methods, which)
          }
      }
  }
  if(!is.null(defaultMethod))
      elNamed(value, "ANY") <- defaultMethod
  value
}


matchSignature <-
  ## Match the signature object (a partially or completely named subset of the
  ## arguments of `fun', and return a vector of all the classes in the order specified
  ## by the signature slot of the generic.  The classes not specified by `signature
  ##' will be `"ANY"' in the value.
  function(signature, fun, where = baseenv())
{
    if(!is(fun, "genericFunction"))
        stop(gettextf("trying to match a method signature to an object (of class %s) that is not a generic function",
                      dQuote(class(fun))),
             domain = NA)
    anames <- fun@signature
    if(length(signature) == 0L)
        return(character())
    if(is(signature,"character")) {
        pkgs <- packageSlot(signature) # includes case of  "ObjectsWithPackage"
        if(is.null(pkgs))
            pkgs <- character(length(signature))
        else if(length(pkgs) != length(signature))
            stop("invalid 'package' slot or attribute, wrong length")
        sigClasses <- as.character(signature)
    }
    else if(is(signature, "list")) {
        sigClasses <- pkgs <- character(length(signature))
        for(i in seq_along(signature)) {
            cli <- signature[[i]]
            if(is(cli, "classRepresentation")) {
                sigClasses[[i]] <- cli@className
                pkgs[[i]] <- cli@package
            }
            else if(is(cli, "character") && length(cli) == 1) {
                sigClasses[[i]] <- cli
                pkgi <- packageSlot(cli)
                if(is.character(pkgi))
                    pkgs[[i]] <- pkgi
            }
            else
                stop(gettextf("invalid element in a list for \"signature\" argument; element %d is neither a class definition nor a class name",
                     i), domain = NA)
        }
    }
    else
        stop(gettextf("trying to match a method signature of class %s; expects a list or a character vector",
                      dQuote(class(signature))),
             domain = NA)
    if(!identical(where, baseenv())) {
        ## fill in package information, warn about undefined classes
        unknown <- !nzchar(pkgs)
        for(i in seq_along(sigClasses)[unknown]) {
            cli <- getClassDef(sigClasses[[i]], where)
            if(!is.null(cli)) {
                pkgs[[i]] <- cli@package
                unknown[[i]] <- FALSE
            }
        }
        if(any(unknown)) {
            unknown <- unique(sigClasses[unknown])
            ## coerce(), i.e., setAs() may use *one* unknown class
	    MSG <- if(identical(as.vector(coerce@generic), "coerce") &&
		      length(unknown) == 1) message
	    else function(...) warning(..., call. = FALSE)
	    MSG(.renderSignature(fun@generic, signature),
		sprintf(ngettext(length(unknown),
				 "no definition for class %s",
				 "no definition for classes %s"),
			paste(dQuote(unknown), collapse = ", ")),
		domain = NA)
        }
    }
    signature <- as.list(signature)
    if(length(sigClasses) != length(signature))
        stop(gettextf("object to use as a method signature for function %s does not look like a legitimate signature (a vector of single class names): there were %d class names, but %d elements in the signature object",
                      sQuote(fun@generic),
                      length(sigClasses),
                      length(signature)),
             domain = NA)
    if(is.null(names(signature))) {
        which <- seq_along(signature)
        if(length(which) > length(anames))
          stop(gettextf("more elements in the method signature (%d) than in the generic signature (%d) for function %s",
	       length(which), length(anames), sQuote(fun@generic)), domain = NA)
    }
    else {
    ## construct a function call with the same naming pattern  &
      ## values as signature
    sigList <- signature
    for(i in seq_along(sigList))
        sigList[[i]] <- c(sigClasses[[i]], pkgs[[i]])
    fcall <- do.call("call", c("fun", sigList))
    ## match the call to the formal signature (usually the formal args)
    if(identical(anames, formalArgs(fun)))
        smatch <- match.call(fun, fcall)
    else {
        fmatch <- fun
        ff <- as.list(anames); names(ff) <- anames
        formals(fmatch, envir = environment(fun)) <- ff
        smatch <- match.call(fmatch, fcall)
    }
    snames <- names(smatch)[-1L]
    which <- match(snames, anames)
    ## Assertion:  match.call has permuted the args into the order of formal args,
    ## and carried along the values.  Get the supplied classes in that
    ## order, from the matched args in the call object.
    if(any(is.na(which)))
        stop(sprintf(ngettext(sum(is.na(which)),
                              "in the method signature for function %s invalid argument name in the signature: %s",
                              "in the method signature for function %s invalid argument names in the signature: %s"),
                     sQuote(fun@generic),
                     paste(snames[is.na(which)], collapse = ", ")),
             domain = NA)
    smatch <- smatch[-1]
    for(i in seq_along(smatch)) {
        eli <- smatch[[i]]
        sigClasses[[i]] <- eli[[1]]
        pkgs[[i]] <- eli[[2]]
    }
}
    n <- length(anames)
    value <- rep("ANY", n)
    valueP <- rep("methods", n)
    names(value) <- anames
    value[which] <- sigClasses
    valueP[which] <- pkgs
    unspec <- value == "ANY"
    ## remove the trailing unspecified classes
    while(n > 1 && unspec[[n]])
        n <- n-1
    length(value) <- length(valueP) <- n
    attr(value, "package") <- valueP
    ## <FIXME> Is there a reason (bootstrapping?) why this
    ## is not an actual object from class "signature"?
    ## See .MakeSignature() </FIXME>
    value
}

showMlist <-
  ## Prints the contents of the MethodsList.  If `includeDefs' the signatures and the
  ## corresponding definitions will be printed; otherwise, only the signatures.
  ##
  ## If `includeDefs' is `TRUE', the currently known inherited methods are included;
  ## otherwise, only the directly defined methods.
function(mlist, includeDefs = TRUE, inherited = TRUE, classes = NULL, useArgNames = TRUE,
         printTo = stdout())
{
    if(identical(printTo, FALSE)) {
        tmp <- tempfile()
        con <- file(tmp, "w")
    }
    else
        con <- printTo
  object <- linearizeMlist(mlist, inherited)
  methods <- object@methods
  signatures <- object@classes
  args <- object@arguments
  if(!is.null(classes) && length(signatures)>0) {
    keep <- !sapply(signatures, function(x, y)all(is.na(match(x, y))), classes)
    methods <- methods[keep]
    signatures <- signatures[keep]
    args <- args[keep]
  }
  if(length(methods) == 0)
    cat(file=con, "<Empty Methods List>\n")
  else {
   n <- length(methods)
    labels <- character(n)
    if(useArgNames) {
      for(i in 1L:n) {
        sigi <- signatures[[i]]
        labels[[i]] <- paste(args[[i]], " = \"", sigi, "\"",
                             sep = "", collapse = ", ")
      }
    }
    else {
      for(i in 1L:n)
        labels[[i]] <- paste(signatures[[i]], collapse = ", ")
    }
    for(i in seq_along(methods)) {
      cat(file=con, (if(includeDefs) "## Signature:" else ""), labels[[i]])
      method <- methods[[i]]
      if(includeDefs) {
        cat(file=con, ":\n")
        if(is(method, "MethodDefinition")) ## really an assertion
          cat(file=con, deparse(method@.Data), sep="\n")
        else
          cat(file=con, deparse(method), sep="\n")
      }
      if(is(method, "MethodDefinition") &&
         !identical(method@target, method@defined)) {
          defFrom <- method@defined
          cat(file = con, if(includeDefs) "##:" else "\n",
              "    (inherited from ",
              paste(names(defFrom), " = \"", as.character(defFrom),
                    "\"", sep = "", collapse = ", "),
               ")", if(includeDefs) "\n", sep="")
      }
      cat(file=con, "\n")
    }
  }
    if(identical(printTo, FALSE)) {
        close(con)
        value <- readLines(tmp)
        unlink(tmp)
        value
    }
}

promptMethods <- function(f, filename = NULL, methods)
{
    ## Generate information in the style of 'prompt' for the methods of
    ## the generic named 'f'.
    ##
    ## 'filename' can be a logical or NA or the name of a file to print
    ## to.  If it 'FALSE', the methods skeleton is returned, to be
    ## included in other printing (typically, the output from 'prompt').

    escape <- function(txt) gsub("%", "\\\\%", txt)
    packageString <- ""

    fdef <- getGeneric(f)
    if(!isGeneric(f, fdef=fdef))
	stop(gettextf("no generic function found corresponding to %s",
                      sQuote(f)),
	     domain = NA)
    if(missing(methods)) {
	methods <- findMethods(fdef)
	## try making  packageString
	where <- .genEnv(fdef, topenv(parent.frame()))
	if(!identical(where, .GlobalEnv))
	    packageString <-
                sprintf("in Package \\pkg{%s}", getPackageName(where))
    }
    fullName <- utils:::topicName("methods", f)
    n <- length(methods)
    labels <- character(n)
    aliases <- character(n)
    signatures <- findMethodSignatures(methods = methods, target=TRUE)
    args <- colnames(signatures) # the *same* for all
    for(i in seq_len(n)) {
        sigi <- signatures[i, ]
	labels[[i]] <-
            sprintf("\\code{signature(%s)}",
                    paste(sprintf("%s = \"%s\"", args, escape(sigi)),
                          collapse = ", "))
	aliases[[i]] <-
	    paste0("\\alias{",
		   utils:::topicName("method", c(f, signatures[i,])),
		   "}")
    }
    text <- paste0("\n\\item{", labels,
                   "}{\n%%  ~~describe this method here~~\n}")
    text <- c("\\section{Methods}{\n\\describe{", text, "}}")
    aliasText <- c(paste0("\\alias{", escape(fullName), "}"), escape(aliases))
    if(identical(filename, FALSE))
        return(c(aliasText, text))

    if(is.null(filename) || identical(filename, TRUE))
        filename <- paste0(fullName, ".Rd")

    Rdtxt <-
        list(name = paste0("\\name{", fullName, "}"),
             type = "\\docType{methods}",
             aliases = aliasText,
             ## <FIXME>
             ## Title and description are ok as auto-generated: should
             ## they be flagged as such (via '~~' which are quite often
             ## left in by authors)?
             title =
             sprintf("\\title{ ~~ Methods for Function \\code{%s} %s ~~}",
                     f, packageString),
             description =
             paste0("\\description{\n ~~ Methods for function",
                    " \\code{", f, "} ",
                    sub("^in Package", "in package", packageString),
                    " ~~\n}"),
             ## </FIXME>
             "section{Methods}" = text,
             keywords = c("\\keyword{methods}",
             "\\keyword{ ~~ other possible keyword(s) ~~ }"))

    if(is.na(filename)) return(Rdtxt)

    cat(unlist(Rdtxt), file = filename, sep = "\n")
    .message("A shell of methods documentation has been written",
             .fileDesc(filename), ".\n")
    invisible(filename)
}

linearizeMlist <-
    ## Undo the recursive nature of the methods list, making a list of
    ## function definitions, with the names of the list being the
    ## corresponding signatures (designed for printing; for looping over
    ## the methods, use `listFromMlist' instead).
    ##
    ## The function calls itself recursively.  `prev' is the previously
    ## selected class names.
    ##
    ## If argument `classes' is provided, only signatures containing one
    ## of these classes will be included.
    function(mlist, inherited = TRUE) {
        methods <- mlist@methods
        allMethods <- mlist@allMethods
        if(inherited && length(allMethods) >= length(methods)) {
            anames <- names(allMethods)
            inh <- is.na(match(anames, names(methods)))
            methods <- allMethods
        }
        preC <- function(y, x)c(x,y) # used with lapply below
        cnames <- names(methods)
        value <- list()
        classes <- list()
        arguments <- list()
        argname <- as.character(mlist@argument)
        for(i in seq_along(cnames)) {
            mi <- methods[[i]]
            if(is.function(mi)) {
                value <- c(value, list(mi))
                classes <- c(classes, list(cnames[[i]]))
                arguments <- c(arguments, list(argname))
            }
            else if(is(mi, "MethodsList")) {
                mi <- Recall(mi, inherited)
                value <- c(value, mi@methods)
                classes <- c(classes, lapply(mi@classes, preC, cnames[[i]]))
                arguments <- c(arguments, lapply(mi@arguments, preC, argname))
            }
            else
                warning(gettextf("skipping methods list element %s of unexpected class %s\n\n",
                                 paste(cnames[i], collapse = ", "),
                                 dQuote(.class1(mi))),
                        domain = NA)
        }
        new("LinearMethodsList", methods = value, classes = classes, arguments = arguments)
    }

print.MethodsList <- function(x, ...)
    showMlist(x)


listFromMlist <-
  ## linearizes the MethodsList object into list(sigs, methods); `prefix' is the partial
  ## signature (a named list of classes) to be prepended to the signatures in this object.
  ##
  ## A utility function used to iterate over all the individual methods in the object.
  function(mlist, prefix = list(), sigs. = TRUE, methods. = TRUE)
{
    methodSlot <- slot(mlist, "methods")
    mnames <- names(methodSlot)
    argName <- as.character(slot(mlist, "argument"))
    sigs <- list()
    methods <- list()
    for(i in seq_along(methodSlot)) {
        thisMethod <- el(methodSlot, i)
        thisClass <- el(mnames, i)
        elNamed(prefix, argName) <- thisClass
        if(is.function(thisMethod)) {
            if(sigs.) sigs <- c(sigs, list(prefix))
            if(methods.) methods <- c(methods, list(thisMethod))
        }
        else {
            more <- Recall(thisMethod, prefix)
            if(sigs.) sigs <- c(sigs, el(more, 1))
            if(methods.) methods <- c(methods, el(more, 2))
        }
    }
    list(sigs, methods)
}

.insertCachedMethods <- function(mlist, argName, Class, fromClass, def) {
    if(is(def, "MethodsList")) {
        ## insert all the cached methods in def
        newArg <- c(argName, as.character(def@argument))
        newDefs <- def@allMethods
        newSigs <- as.list(names(newDefs))
        for(j in seq_along(newDefs))
            mlist <- Recall(mlist, newArg, c(Class, newSigs[[j]]), fromClass,
                            newDefs[[j]])
    }
    else {
        def <- .addMethodFrom(def, argName[1L], Class[1L], fromClass)
        mlist <- insertMethod(mlist, Class, argName, def, TRUE)
    }
    mlist
}

.addMethodFrom <- function(def, arg, Class, fromClass) {
    if(is(def, "MethodDefinition")) {
        ## eventually, we may enforce method definition objects
        ## If not, just leave raw functions alone (NextMethod won't work)
        def@target[[arg]] <- Class
        def@defined[[arg]] <- fromClass
    }
    def
}

## Define a trivial version of asMethodDefinition for bootstrapping.
## The real version requires several class definitions as well as
## methods for as<-
asMethodDefinition <- function(def, signature = list(.anyClassName), sealed = FALSE, fdef = def) {
  if(is.primitive(def))
    def
  else {
    value = new("MethodDefinition")
    value@.Data <- def
    classes <- .MakeSignature(new("signature"),  def, signature, fdef)
        value@target <- classes
        value@defined <- classes
    value
  }
  }

.trimMlist <- function(mlist, fromClass) {
  mlist@methods <- mlist@methods[fromClass]
  mlist@allMethods <- mlist@allMethods[fromClass]
  mlist
}

.noMlistsFlag <- TRUE
.noMlists <- function() {
   ## if this were to be dynamically variable, but
  ## it can't, IMO
  ## identical(getOption("noMlists"), TRUE)
  ## so instead
  .noMlistsFlag
}

.MlistDepTable <- new.env()
.MlistDeprecated <- function(this = "<default>", instead) {
    if(is.character(this)) {
        if(exists(this, envir = .MlistDepTable, inherits = FALSE))
            return()
        else
            assign(this, TRUE, envir = .MlistDepTable)
    }
    if(missing(this))
        msg <-"Use of the \"MethodsList\" meta data objects is deprecated."
    else if(is.character(this))
        msg <- gettextf("%s, along with other use of the \"MethodsList\" metadata objects, is deprecated.", dQuote(this))
    else
        msg <- gettextf("in %s: use of \"MethodsList\" metadata objects is deprecated.", deparse(this))
    if(!missing(instead))
      msg <- paste(msg, gettextf("use %s instead.", dQuote(instead)))
    msg <- paste(msg, "see ?MethodsList. (This warning is shown once per session.)")
    base::.Deprecated(msg = msg)
}

#  File src/library/methods/R/MethodsListClass.R
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

.InitMethodsListClass <- function(envir)
{
    if(exists(classMetaName("MethodsList"), envir))
        return(FALSE)
    clList <- character()
    setClass("MethodsList",
             representation(methods = "list", argument = "name", allMethods = "list"),
             where = envir); clList <- c(clList, "MethodsList")
    setClass("EmptyMethodsList", representation(argument = "name", sublist = "list"),
             where = envir); clList <- c(clList, "EmptyMethodsList")

    ## the classes for method definitions
    setClass("PossibleMethod", where = envir); clList <- c(clList, "PossibleMethod")
    ## functions (esp. primitives) are methods
    setIs("function", "PossibleMethod", where = envir)

    ## the default slot of a generic function can be a method, primitive or NULL
    setClass("optionalMethod", where = envir); clList <- c(clList, "optionalMethod")
    setIs("PossibleMethod", "optionalMethod", where = envir)
    setIs("NULL", "optionalMethod", where = envir)
    ## prior to 2.11.0, the default slot in generic function objects was a MethodsList or NULL
    setIs("MethodsList", "optionalMethod", where = envir) #only until MethodsList class is defunct

    ## signatures -- multiple class names w. package slot in ||
    setClass("signature", representation("character", names = "character", package = "character"), where = envir); clList <- c(clList, "signature")

    ## className -- a single class name with package
    setClass("className", contains = "character",
             representation(package = "character"))

    ## formal method definition for all but primitives
    setClass("MethodDefinition", contains = "function",
             representation(target = "signature", defined = "signature", generic = "character"),
             where = envir); clList <- c(clList, "MethodDefinition")
    ## class for default methods made from ordinary functions
    setClass("derivedDefaultMethod", "MethodDefinition")
    ## class for methods with precomputed information for callNextMethod
    setClass("MethodWithNext",
             representation("MethodDefinition", nextMethod = "PossibleMethod", excluded = "list"), where = envir); clList <- c(clList, "MethodWithNext")
    setClass("SealedMethodDefinition", contains = "MethodDefinition"); clList <- c(clList, "SealedMethodDefinition")
    setClass("genericFunction", contains = "function",
             representation( generic = "character", package = "character",
                            group = "list", valueClass = "character",
                            signature = "character", default = "optionalMethod",
                            skeleton = "call"), where = envir); clList <- c(clList, "genericFunction")
    ## standard generic function -- allows immediate dispatch
    setClass("standardGeneric",  contains = "genericFunction")
    setClass("nonstandardGeneric", # virtual class to mark special generic/group generic
             where = envir); clList <- c(clList, "nonstandardGeneric")
    setClass("nonstandardGenericFunction",
             representation("genericFunction", "nonstandardGeneric"),
             where = envir); clList <- c(clList, "nonstandardGenericFunction")
    setClass("groupGenericFunction",
             representation("genericFunction", groupMembers = "list"),
             where = envir); clList <- c(clList, "groupGenericFunction")
    setClass("nonstandardGroupGenericFunction",
             representation("groupGenericFunction", "nonstandardGeneric"),
             where = envir); clList <- c(clList, "nonstandardGroupGenericFunction")
    setClass("LinearMethodsList", representation(methods = "list", arguments = "list",
                                                 classes = "list", generic = "genericFunction"),
             where = envir); clList <- c(clList, "LinearMethodsList")
    setClass("ObjectsWithPackage", representation("character", package = "character"),
             where = envir); clList <- c(clList, "ObjectsWithPackage")
    assign(".SealedClasses", c(get(".SealedClasses", envir), clList), envir)
    TRUE
}

## some initializations that need to be done late
.InitMethodDefinitions <- function(envir) {
    assign("asMethodDefinition",
           function(def, signature = list(.anyClassName), sealed = FALSE, fdef = def) {
        ## primitives can't take slots, but they are only legal as default methods
        ## and the code will just have to accomodate them in that role, w/o the
        ## MethodDefinition information.
        ## NULL is a valid def, used to remove methods.
        switch(typeof(def),
               "builtin" = , "special" = , "NULL" = return(def),
               "closure" = {},
               stop(gettextf("invalid object for formal method definition: type %s",
                             dQuote(typeof(def))),
                    domain = NA)
               )
        if(is(def, "MethodDefinition")) {
            value <- def
            if(missing(signature))
                signature <- value@defined
        }
        else
            value <- new("MethodDefinition", def)

        if(sealed)
            value <- new("SealedMethodDefinition", value)
        if(is(signature, "signature"))
            classes <- signature
        else
            classes <- .MakeSignature(new("signature"),  def, signature, fdef)
        value@target <- classes
        value@defined <- classes
        value
    }, envir = envir)
    setGeneric("loadMethod", where = envir)
    setMethod("loadMethod", "MethodDefinition",
              function(method, fname, envir) {
                  assign(".target", method@target, envir = envir)
                  assign(".defined", method@defined, envir = envir)
                  assign(".Method", method, envir = envir)
                  method
              }, where = envir)
    setMethod("loadMethod", "MethodWithNext",
              function(method, fname, envir) {
                  callNextMethod()
                  assign(".nextMethod", method@nextMethod, envir = envir)
                  method
              }, where = envir)
    setGeneric("addNextMethod", function(method, f = "<unknown>",
                                         mlist, optional = FALSE, envir)
               standardGeneric("addNextMethod"), where = envir)
    setMethod("addNextMethod", "MethodDefinition",
	      function(method, f, mlist, optional, envir) {
		  .findNextFromTable(method, f, optional, envir)
	      }, where = envir)
    setMethod("addNextMethod", "MethodWithNext",
	      function(method, f, mlist, optional, envir) {
		  .findNextFromTable(method, f, optional, envir, method@excluded)
	      }, where = envir)

    .initGeneric <- function(.Object, ...) {
            value <- standardGeneric("initialize")
            if(!identical(class(value), class(.Object))) {
                cv <- class(value)
                co <- class(.Object)
                if(.identC(cv[[1L]], co)) {
                  ## ignore S3 with multiple classes  or basic classes
                    if(is.na(match(cv, .BasicClasses)) &&
                       length(cv) == 1L) {
                        warning(gettextf("missing package slot (%s) in object of class %s (package info added)",
                                         packageSlot(co),
                                         dQuote(class(.Object))),
                                domain = NA)
                        class(value) <- class(.Object)
                    }
                    else
                        return(value)
                }
                else
                    stop(gettextf("'initialize' method returned an object of class %s instead of the required class %s",
                                  paste(dQuote(class(value)), collapse=", "),
                                  dQuote(class(.Object))),
                         domain = NA)
            }
            value
        }
    if(!isGeneric("initialize", envir)) {
        ## save the default method
        assign(".initialize", initialize, envir)
        setGeneric("initialize",  .initGeneric, where = envir, useAsDefault = TRUE, simpleInheritanceOnly = TRUE)
    }
    setMethod("initialize", "signature",
              function(.Object, functionDef, ...) {
                  if(nargs() < 2)
                      .Object
                  else if(missing(functionDef))
                      .MakeSignature(.Object, , list(...))
                  else if(!is(functionDef, "function"))
                      .MakeSignature(.Object, , list(functionDef, ...))
                  else
                      .MakeSignature(.Object, functionDef, list(...))
              }, where = envir)
    setMethod("initialize", "environment", # only for new("environment",...); see .InitSpecialTypesAndClasses for subclasses
              function(.Object, ...) {
                  value <- new.env()
                  args <- list(...)
                  objs <- names(args)
                  for(what in objs)
                      assign(what, elNamed(args, what), envir = value)
                  value
              }, where = envir)
    ## from 2.11.0, the MethodsList classs is deprecated
    setMethod("initialize", "MethodsList", function(.Object, ...) {
        .MlistDeprecated()
        callNextMethod()
    }, where = envir)

    ## make sure body(m) <- .... leaves a method as a method
    setGeneric("body<-", where = envir)
    setMethod("body<-", "MethodDefinition", function (fun, envir, value) {
        ff <- as(fun, "function")
        body(ff, envir = envir) <- value
        fun@.Data <- ff
        fun
    }, where = envir)
    ## a show method for lists of generic functions, etc; see metaNameUndo
    if(!isGeneric("show", envir))
        setGeneric("show", where = envir, simpleInheritanceOnly = TRUE)
    setMethod("show", "ObjectsWithPackage",
              function(object) {
                  pkg <- object@package
                  data <- as(object, "character")
                  cat("An object of class \"", class(object), "\":\n", sep="")
                  if(length(unique(pkg))==1) {
                      show(data)
                      cat("(All from \"", unique(pkg), "\")\n", sep="")
                  }
                  else {
                      mat <- rbind(data, pkg)
		      dimnames(mat) <- list(c("Object:", "Package:"),
					    rep("", length(data)))
                      show(mat)
                  }
              }, where = envir)
    ## show method for reports of method selection ambiguities; see MethodsTable.R
    setMethod("show", "MethodSelectionReport", where = envir,
              function(object) {
                  nreport <- length(object@target)
                  cat(sprintf(ngettext(nreport,
                                       "Reported %d ambiguous selection out of %d for function %s\n",
                                       "Reported %d ambiguous selections out of %d for function %s\n"),
                              nreport, length(object@allSelections), object@generic))
                  target <- object@target; selected = object@selected
                  candidates <- object@candidates; note <- object@note
                  for(i in seq_len(nreport)) {
                      these <- candidates[[i]]; notei <- note[[i]]
                      these <- these[is.na(match(these, selected[[i]]))]
                      cat(gettextf(
                                   '%d: target "%s": chose "%s" (others: %s)',
                                   i,target[[i]], selected[[i]], paste0('"', these, '"', collapse =", ")))
                      if(nzchar(notei))
                          cat(gettextf("\n    Notes: %s.\n", notei))
                      else
                          cat(".\n")
                  }
                  NULL
              })
    setMethod("show", "classGeneratorFunction", where = envir,
              function(object) {
                  cat(gettextf("class generator function for class %s from package %s\n",
                               dQuote(object@className),
                               sQuote(object@package)))
                  show(as(object, "function"))
              })

    setGeneric("cbind2", function(x, y, ...) standardGeneric("cbind2"),
	       where = envir)
    ## and its default methods:
    setMethod("cbind2", signature(x = "ANY", y = "ANY"),
	      function(x,y) .__H__.cbind(deparse.level = 0, x, y) )
    setMethod("cbind2", signature(x = "ANY", y = "missing"),
	      function(x,y) .__H__.cbind(deparse.level = 0, x) )

    setGeneric("rbind2", function(x, y, ...) standardGeneric("rbind2"),
	       where = envir)
    ## and its default methods:
    setMethod("rbind2", signature(x = "ANY", y = "ANY"),
	      function(x,y) .__H__.rbind(deparse.level = 0, x, y) )
    setMethod("rbind2", signature(x = "ANY", y = "missing"),
	      function(x,y) .__H__.rbind(deparse.level = 0, x) )

    setGeneric("kronecker", where = envir)

    setMethod("kronecker", signature(X = "ANY", Y = "ANY"),
	      function(X, Y, FUN = "*", make.dimnames = FALSE, ...)
              .kronecker(X, Y, FUN = FUN, make.dimnames = make.dimnames, ...))

    .InitStructureMethods(envir)
### Uncomment next line if we want special initialize methods for basic classes
    .InitBasicClassMethods(envir)
}

.InitStructureMethods <- function(where) {
    ## these methods need to be cached (for the sake of the primitive
    ## functions in the group) if a class is loaded that extends
    ## one of the classes in `needed` (other classes than "structure" now
    ## also require generics for some primitives).
    if(!exists(".NeedPrimitiveMethods", where))
      needed <- list()
    else
      needed <- get(".NeedPrimitiveMethods", where)
    needed <- c(needed, list(structure = "Ops", vector = "Ops",
          array = "Ops", nonStructure = "Ops"),
          array = "[", structure = "[", nonStructure = "[",
          structure = "Math", nonStructure = "Math",
          refClass = "$", refClass = "$<-", data.frame = "$<-"
                )
    assign(".NeedPrimitiveMethods", needed, where)
    setMethod("Ops", c("structure", "vector"), where = where,
              function(e1, e2) {
                  value <- callGeneric(e1@.Data, e2)
                  if(length(value) == length(e1)) {
                      e1@.Data <- value
                      e1
                  }
                  else
                    value
              })
    setMethod("Ops", c("vector", "structure"), where = where,
              function(e1, e2) {
                  value <- callGeneric(e1, e2@.Data)
                  if(length(value) == length(e2)) {
                      e2@.Data <- value
                      e2
                  }
                  else
                    value
              })
    setMethod("Ops", c("structure", "structure"), where = where,
              function(e1, e2)
                 callGeneric(e1@.Data, e2@.Data)
              )
    ## We need some special cases for matrix and array.
    ## Although they extend "structure", their .Data "slot" is the matrix/array
    ## So op'ing them with a structure gives the matrix/array:  Not good?
    ## Following makes them obey the structure rule.
    setMethod("Ops", c("structure", "array"), where = where,
              function(e1, e2)
                 callGeneric(e1@.Data, as.vector(e2))
              )
    setMethod("Ops", c("array", "structure"), where = where,
              function(e1, e2)
                 callGeneric(as.vector(e1), e2@.Data)
              )
    ## but for two array-based strucures, we let the underlying
    ## code for matrix/array stand.
    setMethod("Ops", c("array", "array"), where = where,
              function(e1, e2)
                 callGeneric(e1@.Data, e2@.Data)
              )


    setMethod("Math", "structure", where = where,
              function(x) {
                  x@.Data <- callGeneric(x@.Data)
                  x
              })
    setMethod("Math2", "structure", where = where,
              function(x, digits) {
                  value <- x
                  x <- x@.Data
                  value@Data  <- callGeneric()
                  value
              })
    ## some methods for nonStructure, ensuring that the class and slots
    ## will be discarded
    setMethod("Ops", c("nonStructure", "vector"), where = where,
              function(e1, e2) {
                  callGeneric(e1@.Data, e2)
              })
    setMethod("Ops", c("vector", "nonStructure"), where = where,
              function(e1, e2) {
                  callGeneric(e1, e2@.Data)
              })
    setMethod("Ops", c("nonStructure", "nonStructure"), where = where,
              function(e1, e2)
                 callGeneric(e1@.Data, e2@.Data)
              )
    setMethod("Math", "nonStructure", where = where,
              function(x) {
                  callGeneric(x@.Data)
              })
    setMethod("Math2", "nonStructure", where = where,
              function(x, digits) {
                  x <- x@.Data
                  callGeneric()
              })
    setMethod("[", "nonStructure", where = where,
                        function (x, i, j, ..., drop = TRUE)
                        {
                          value <- callNextMethod()
                          value@.Data
                        })

}


.MakeSignature <- function(object, def = NULL, signature, fdef = def) {
    ## fill in the signature information in object
    ## In effect, object must come from class "signature" or a subclass
    ## but the only explicit requirement is that it has compatible
    ## .Data and "package" slots
    signature <- unlist(signature)
    if(length(signature)>0) {
        classes <- as.character(signature)
        sigArgs <- names(signature)
        pkgs <- attr(signature, "package")
        if(is.null(pkgs))
            pkgs <- character(length(signature))
        if(is(fdef, "genericFunction"))
            formalNames <- fdef@signature
        else if(is.function(def)) {
            if(!is(fdef, "function")) fdef <- def
            formalNames <- formalArgs(fdef)
            dots <- match("...", formalNames)
            if(!is.na(dots))
                formalNames <- formalNames[-dots]
        }
        else formalNames <- character()
        if(length(formalNames) > 0) {
            if(is.null(sigArgs))
              names(signature) <- formalNames[seq_along(classes)]
            else if(length(sigArgs) && any(is.na(match(sigArgs, formalNames))))
              stop(gettextf("the names in signature for method (%s) do not match %s's arguments (%s)",
                            paste(sigArgs, collapse = ", "),
                            if(is(fdef, "genericFunction")) fdef@generic else "function",
                            paste(formalNames, collapse = ", ")),
                   domain = NA)
        }
        object@.Data <- signature
        object@package <- pkgs
    }
    object
}
#  File src/library/methods/R/NextMethod.R
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

callNextMethod <- function(...) {
    method <- nextMethod <-  NULL
    dotNextMethod <- as.name(".nextMethod")
    ## 2 environments are used here:  callEnv, from which the .nextMethod call
    ## takes place; and methodEnv, the method environment used to find the next method
    ## Because of the .local mechanism used to allow variable argument lists
    ## in methods (see rematchDefinition) these may be different.
    parent <- sys.parent(1)
    maybeMethod <- sys.function(parent)
    if(is(maybeMethod, "MethodDefinition")) {
        callEnv <- methodEnv <- parent.frame(1)
        mcall <- sys.call(parent)
        i <- 1
    }
    else {
        callEnv <- parent.frame(1)
        methodEnv <- parent.frame(2)
        mcall <- sys.call(sys.parent(2))
        i <- 2
    }
    ## set up the nextMethod object, load it
    ## into the calling environment, and cache it
    if(exists(".Method", envir = methodEnv, inherits = FALSE)) {
        ## call to standardGeneric(f)
        method <- get(".Method", envir = methodEnv, inherits = FALSE)
        if(exists(".nextMethod", envir = callEnv, inherits = FALSE))
            nextMethod <- get(".nextMethod", envir = callEnv)
        f <- get(".Generic", envir = methodEnv)
    }
    else if(identical(mcall[[1L]], dotNextMethod)) {
        ## a call from another callNextMethod()
        nextMethodEnv <- parent.frame(i+1)
        nextMethod <- get(".nextMethod", nextMethodEnv)
        f <- get(".Generic", envir = nextMethodEnv)
    }
    else {
        ## may be a method call for a primitive; not available as .Method
        f <- as.character(mcall[[1L]])
        fdef <- genericForPrimitive(f)
        ## check that this could be a basic function with methods
        if(is.null(fdef))
            stop(gettextf("a call to callNextMethod() appears in a call to %s, but the call does not seem to come from either a generic function or another 'callNextMethod'",
                          sQuote(f)),
                 domain = NA)
        f <- fdef@generic
        method <- maybeMethod
    }
    if(is(method, "MethodDefinition")) {
        if(is.null(nextMethod)) {
            if(!is(method, "MethodWithNext")) {
                method <- addNextMethod(method, f, envir=methodEnv)
                ## cache the method with the nextMethod included,
                ## so later calls will load this information.
                cacheMethod(f, method@target, method, fdef = getGeneric(f), inherited = TRUE)
            }
            nextMethod <- method@nextMethod
            assign(".nextMethod", nextMethod, envir = callEnv)
            assign(".Generic", f, envir = callEnv)
        }
    }
    else if(is.null(method)) {
        if(is.null(nextMethod))
            stop("call to 'callNextMethod' does not appear to be in a 'method' or 'callNextMethod' context")
        ## else, callNextMethod() from another callNextMethod
        method <- nextMethod
        if(!is(method, "MethodWithNext")) {
            method <- addNextMethod(method, f, envir=methodEnv)
        }
        nextMethod <- method@nextMethod
        ## store the nextmethod in the previous nextmethod's
        assign(".nextMethod", nextMethod, envir = callEnv)
        assign(".Generic", f, envir = callEnv)
        assign(".nextMethod", method, envir = nextMethodEnv)
        assign(".Generic", f, envir = nextMethodEnv)
    }
    else
        stop(gettextf("bad object found as method (class %s)",
                      dQuote(class(method))), domain = NA)
    subsetCase <- !is.na(match(f, .BasicSubsetFunctions))
    if(nargs()>0) {
      call <- sys.call()
      call[[1L]] <- as.name(".nextMethod")
      eval(call, callEnv)
      }
    else {
        if(subsetCase) {
            ## don't use match.call, because missing args will screw up for "[", etc.
            call <- as.list(mcall)
            ## don't test with identical(), there may  be a package attr.
            if((f ==  "[") && length(names(call)>0))
                call <- .doSubNextCall(call, method) # [ with a drop= arg.
            else {
               fnames <- c("", formalArgs(method))
               i <- match("...",fnames)
               if(is.na(i) || i > length(call))
                   length(fnames) <- length(call)
               else {
                   i <- i-1
                   length(fnames) <- i
                   fnames <- c(fnames, rep("", length(call) - i))
               }
               names(call) <- fnames
               call <- as.call(call)
           }
        }
        else
            call <- match.call(method, mcall, expand.dots = FALSE)
        .Call(C_R_nextMethodCall, call, callEnv)
    }
}

loadMethod <- function(method, fname, envir)
    method

.doSubNextCall <- function(call, method) {
    idrop <- match("drop", names(call))
    hasDrop <- !is.na(idrop)
    if(hasDrop) {
        drop <- call$drop
        call <- call[-idrop]
    }
    fnames <- c("", formalArgs(method))
    i <- match("...",fnames)
    if(is.na(i) || i > length(call))
        length(fnames) <- length(call)
    else {
        i <- i-1
        length(fnames) <- i
        fnames <- c(fnames, rep("", length(call) - i))
    }
    names(call) <- fnames
    if(hasDrop)
        call$drop <- drop
    as.call(call)
}
#  File src/library/methods/R/RClassUtils.R
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

testVirtual <-
  ## Test for a Virtual Class.
  ## Figures out, as well as possible, whether the class with these properties,
  ## extension, and prototype is a virtual class.
  ## Can be forced to be virtual by extending "VIRTUAL".  Otherwise, a class is
  ## virtual only if it has no slots, extends no non-virtual classes, and has a
  ## NULL Prototype
  function(properties, extends, prototype, where)
{
    if(length(extends)) {
        en <- names(extends)
        if(!is.na(match("VIRTUAL", en)))
            return(TRUE)
        ## does the class extend a known non-virtual class?
        for(what in en) {
            enDef <- getClassDef(what, where)
            if(!is.null(enDef) && identical(enDef@virtual, FALSE))
                return(FALSE)
        }
    }
    (length(properties) == 0L && is.null(prototype))
}

makePrototypeFromClassDef <-
  ## completes the prototype implied by
  ## the class definition.
  ##
  ##  The following three rules are applied in this order.
  ##
  ## If the class has slots, then the prototype for each
  ## slot is used by default, but a corresponding element in the explicitly supplied
  ## prototype, if there is one, is used instead (but it must be coercible to the
  ## class of the slot).
  ##
  ## If there are no slots but a non-null prototype was specified, this is returned.
  ##
  ## If there is a single non-virtual superclass (a class in the extends list),
  ## then its prototype is used.
  ##
  ## If all three of the above fail, the prototype is `NULL'.
  function(slots, ClassDef, extends, where)
{
    className <- ClassDef@className
    snames <- names(slots)
    ## try for a single superclass that is not virtual
    supers <- names(extends)
    virtual <- NA
    dataPartClass <- elNamed(slots, ".Data")
    prototype <- ClassDef@prototype
    dataPartDone <- is.null(dataPartClass)  || is(prototype, dataPartClass)# don't look for data part in supreclasses
    ## check for a formal prototype object (TODO:  sometime ensure that this happens
    ## at setClass() time, so prototype slot in classRepresentation can have that class
    if(!.identC(class(prototype), className) && .isPrototype(prototype)) {
        pnames <- prototype@slots
        prototype <- prototype@object
    }
    else
        pnames <- names(attributes(prototype))
    if(length(slots) == 0L && !is.null(prototype))
            return(prototype)
    for(i in seq_along(extends)) {
        what <- el(supers, i)
        exti <- extends[[i]]
        if(identical(exti@simple, FALSE))
            next ## only simple contains rel'ns give slots
        if(identical(what, "VIRTUAL"))
            ## the class is virtual, and the prototype usually NULL
            virtual <- TRUE
        else if(isClass(what, where = where)) {
            cli <- getClass(what, where = where)
            slotsi <- names(cli@slots)
            pri <- cli@prototype
            ## once in a while
            if(is.null(prototype)) {
                prototype <- pri
                pnames <- names(attributes(prototype))
                fromClass <- what
            }
            else if(length(slots)) {
                for(slotName in slotsi) {
                    if(identical(slotName, ".Data")) {
                        if(!dataPartDone) {
                            prototype <- setDataPart(prototype, getDataPart(pri), FALSE)
                            dataPartDone <- TRUE
                        }
                    }
                    else if(is.na(match(slotName, pnames))) {
                        ## possible that the prototype already had this slot specified
                        ## If not, add it now.
                        attr(prototype, slotName) <- attr(pri, slotName)
                        pnames <- c(pnames, slotName)
                    }
                }
            }
            else if(!dataPartDone && extends(cli, dataPartClass)) {
                 prototype <- setDataPart(prototype, pri, FALSE)
                 dataPartDone <- TRUE
            }
        }
    }
    if(length(slots) == 0L)
        return(prototype)
    if(is.null(prototype))
        prototype <- defaultPrototype()
    pnames <- names(attributes(prototype))
    ## watch out for a prototype of this class.  Not supposed to happen, but will
    ## at least for the basic class "ts", and can lead to inf. recursion
    pslots <-
        if(.identC(class(prototype), className))
            names(attributes(unclass(prototype)))
        else if(isClass(class(prototype)))
            names(getSlots(getClass(class(prototype))))
        ## else NULL

    ## now check that all the directly specified slots have corresponding elements
    ## in the prototype--the inherited slots were done in the loop over extends
    if(!is.na(match(".Data", snames))) {
        dataPartClass <- elNamed(slots, ".Data")

        ## check the data part
        if(!(isVirtualClass(dataPartClass))) {
            if(isClass(class(prototype), where = where)) {
                prototypeClass <- getClass(class(prototype), where = where)
                OK <- extends(prototypeClass, dataPartClass)
            }
            else
                OK <- FALSE
            if(identical(OK, FALSE))
                stop(gettextf("in constructing the prototype for class %s: prototype has class %s, but the data part specifies class %s",
                              dQuote(className),
                              dQuote(.class1(prototype)),
                              dQuote(dataPartClass)),
                     domain = NA)
        }
        iData <- -match(".Data", snames)
        snames <- snames[iData]
        slots <- slots[iData]
    }
    for(j in seq_along(snames)) {
        name <- el(snames, j)
        i <- match(name, pnames)
        if(is.na(i)) {
            ## if the class of the j-th element of slots is defined and non-virtual,
            ## generate an object from it; else insert NULL
            slot(prototype, name, check = FALSE) <- tryNew(el(slots, j), where)
        }
    }
    extra <- pnames[is.na(match(pnames, snames)) & !is.na(match(pnames, pslots))]
    if(length(extra) && is.na(match("oldClass", supers)))
        warning(gettextf("in constructing the prototype for class %s, slots in prototype and not in class: %s",
                         dQuote(className),
                         paste(extra, collapse=", ")),
                domain = NA)
    ## now check the elements of the prototype against the class definition
    slotDefs <- getSlots(ClassDef); slotNames <- names(slotDefs)
    pnames <- names(attributes(prototype))
    pnames <- pnames[!is.na(match(pnames, slotNames))]
    check <- rep.int(FALSE, length(pnames))
    for(what in pnames) {
        pwhat <- slot(prototype, what)
        slotClass <- getClassDef(slotDefs[[what]], where)
        if(is.null(slotClass) || !extends(class(pwhat), slotClass)) {
            if(is.null(pwhat)) { # does this still apply??
            }
            else if(is(slotClass, "classRepresentation") &&
                    slotClass@virtual) {} # no nonvirtual prototype;e.g. S3 class
            else
                check[match(what, pnames)] <- TRUE
        }
    }
    if(any(check))
        stop(gettextf("in making the prototype for class %s elements of the prototype failed to match the corresponding slot class: %s",
                      dQuote(className),
                      paste(pnames[check],
                            "(class",
                            .dQ(slotDefs[match(pnames[check], slotNames)]),
                            ")",
                            collapse = ", ")),
             domain = NA)
    prototype
}

newEmptyObject <-
  ## Utility function to create an empty object into which slots can be
  ## set.  Currently just creates an empty list with class "NULL"
  ##
  ## Later version should create a special object reference that marks an
  ## object currently with no slots and no data.
  function()
{
    value <- list()
    value
}


completeClassDefinition <-
  ## Completes the definition of Class, relative to the current environment
  ##
  ## The completed definition is stored in the session's class metadata,
  ## to be retrieved the next time that getClass is called on this class,
  ## and is returned as the value of the call.
  function(Class, ClassDef = getClassDef(Class), where, doExtends = TRUE)
{
    ClassDef <- .completeClassSlots(ClassDef, where)
    immediate <- ClassDef@contains
    properties <- ClassDef@slots
    prototype <- makePrototypeFromClassDef(properties, ClassDef, immediate, where)
    virtual <- ClassDef@virtual
    validity <- ClassDef@validity
    access <- ClassDef@access
    package <- ClassDef@package
    extends    <- if(doExtends) completeExtends   (ClassDef, where = where) else ClassDef@contains
    subclasses <- if(doExtends) completeSubclasses(ClassDef, where = where) else ClassDef@subclasses
    if(is.na(virtual))
        ## compute it from the immediate extensions, but all the properties
        virtual <- testVirtual(properties, immediate, prototype, where)
    ## modify the initial class definition object, rather than creating
    ## a new one, to allow extensions of "classRepresentation"
    ## Done by a separate function to allow a bootstrap version.
    ClassDef <- .mergeClassDefSlots(ClassDef,
                                    slots = properties,
                                    contains = extends,
                                    prototype = prototype,
                                    virtual = virtual,
                                    subclasses = subclasses)
    if(any(!is.na(match(names(ClassDef@subclasses), names(ClassDef@contains))))
       && getOption("warn") > 0 ## NEEDED:  a better way to turn on strict testing
       ) {
        bad <- names(ClassDef@subclasses)[!is.na(match(names(ClassDef@subclasses), names(ClassDef@contains)))]
        warning(gettextf("potential cycle in class inheritance: %s has duplicates in superclasses and subclasses (%s)",
                         dQuote(Class),
                         paste(bad, collapse = ", ")),
                domain = NA)
    }
    ClassDef
}

.completeClassSlots <- function(ClassDef, where) {
        properties <- ClassDef@slots
        simpleContains <- ClassDef@contains
        Class <- ClassDef@className
        package <- ClassDef@package
        ext <- getAllSuperClasses(ClassDef, TRUE)
        ## ext has the names of all the direct and indirect superClasses but NOT those that do
        ## an explicit coerce (we can't conclude anything about slots, etc. from them)
        if(length(ext)) {
            superProps <- vector("list", length(ext)+1L)
            superProps[[1L]] <- properties
            for(i in seq_along(ext)) {
                eClass <- ext[[i]]
                if(isClass(eClass, where = where))
                    superProps[[i+1]] <- getClassDef(eClass, where = where)@slots
            }
            properties <- unlist(superProps, recursive = FALSE)
            ## check for conflicting slot names
            if(anyDuplicated(allNames(properties))) {
                duped <- duplicated(names(properties))
#TEMPORARY -- until classes are completed in place & we have way to match non-inherited slots
                properties <- properties[!duped]
#                 dupNames <- unique(names(properties)[duped])
#                 if(!is.na(match(".Data", dupNames))) {
#                     dataParts <- seq_along(properties)[names(properties) == ".Data"]
#                     dupNames <- dupNames[dupNames != ".Data"]
#                     ## inherited data part classes are OK but should be consistent
#                     dataPartClasses <- unique(as.character(properties[dataParts]))
#                     if(length(dataPartClasses)>1)
#                         warning("Inconsistent data part classes inherited (",
#                                 paste(dataPartClasses, collapse = ", "),
#                                 "): coercion to some may fail")
#                     ## remove all but the first .Data
#                     properties <- properties[-dataParts[-1L]]
#                 }
#                 if(length(dupNames)>0) {
#                     dupClasses <- logical(length(superProps))
#                     for(i in seq_along(superProps)) {
#                         dupClasses[i] <- !all(is.na(match(dupNames, names(superProps[[i]]))))
#                     }
#                     stop(paste("Duplicate slot names: slots ",
#                                paste(dupNames, collapse =", "), "; see classes ",
#                                paste0(c(Class, ext)[dupClasses], collapse = ", ")))
#                }
            }
        }
        ## ensure that each element of the slots is a valid class reference
        undefClasses <- rep.int(FALSE, length(properties))
        for(i in seq_along(properties)) {
            cli <- properties[[i]]
            if(is.null(packageSlot(cli))) {
                cliDef <- getClassDef(cli, where)
                if(is.null(cliDef))
                    undefClasses[[i]] <- TRUE
                else
                    packageSlot(properties[[i]]) <- cliDef@package
            }
            else {
                cliDef <- getClassDef(cli)
                if(is.null(cliDef))
                    undefClasses[[i]] <- TRUE
            }
        }
        if(any(undefClasses))
            warning(gettextf("undefined slot classes in definition of %s: %s",
                             .dQ(ClassDef@className),
                             paste(names(properties)[undefClasses], "(class ",
                                   .dQ(unlist(properties, recursive = FALSE)[undefClasses]),
                                   ")", collapse = ", ", sep = "")),
                    call. = FALSE, domain = NA)
        ClassDef@slots <- properties
        ClassDef
}

.uncompleteClassDefinition <- function(ClassDef, slotName) {
    if(missing(slotName)) {
        ClassDef <- Recall(ClassDef, "contains")
        Recall(ClassDef, "subclasses")
    }
    else {
        prev <- slot(ClassDef, slotName)
        if(length(prev)) {
            indir <- sapply(prev, .isIndirectExtension)
            slot(ClassDef, slotName) <- slot(ClassDef, slotName)[!indir]
        }
        ClassDef
    }
}

.isIndirectExtension <- function(object) {
    is(object, "SClassExtension") && length(object@by) > 0
}

.mergeSlots <- function(classDef1, classDef2) {

}

.directSubClasses <- function(ClassDef) {
    ## no checks for input here:
    if(length(sc <- ClassDef@subclasses)) {
        names(sc)[sapply(sc, function(cc) cc@distance == 1L)]
    } ## else NULL
}

getAllSuperClasses <-
  ## Get the names of all the classes that this class definition extends.
  ##
  ## A utility function used to complete a class definition.  It
  ## returns all the superclasses reachable from this class, in
  ## depth-first order (which is the order used for matching methods);
  ## that is, the first direct superclass followed by all its
  ## superclasses, then the next, etc.  (The order is relevant only in
  ## the case that some of the superclasses have multiple inheritance.)
  ##
  ## The list of superclasses is stored in the extends property of the
  ## session metadata.  User code should not need to call
  ## getAllSuperClasses directly; instead, use getClass()@contains
  ## (which will complete the definition if necessary).
  function(ClassDef, simpleOnly = TRUE) {
    temp <- superClassDepth(ClassDef, simpleOnly = simpleOnly)
    unique(temp$label[sort.list(temp$depth)])
  }

superClassDepth <-
    ## all the superclasses of ClassDef, along with the depth of the relation
    ## Includes the extension definitions, but these are not currently used by
    ## getAllSuperClasses
  function(ClassDef, soFar = ClassDef@className, simpleOnly = TRUE)
{
    ext <- ClassDef@contains
    ## remove indirect and maybe non-simple superclasses (latter for inferring slots)
    ok <- rep.int(TRUE, length(ext))
    for(i in seq_along(ext)) {
        exti <- ext[[i]]
        if(.isIndirectExtension(exti) ||
           (simpleOnly && ! exti @simple))
            ok[i] <- FALSE
    }
    ext <- ext[ok]
    immediate <- names(ext)
    notSoFar <- is.na(match(immediate, soFar))
    immediate <- immediate[notSoFar]
    super <- list(label = immediate, depth = rep.int(1, length(immediate)),
                  ext = ext)
    for(i in seq_along(immediate)) {
        what <- immediate[[i]]
        if(!is.na(match(what, soFar)))
           ## watch out for loops (e.g., matrix/array have mutual is relationship)
           next
        exti <- ext[[i]]
        soFar <- c(soFar, what)
        if(!is(exti, "SClassExtension"))
            stop(gettextf("in definition of class %s, information for superclass %s is of class %s (expected \"SClassExtension\")",
                          dQuote(ClassDef@className),
                          dQuote(what),
                          dQuote(class(exti))),
                 domain = NA)
        superClass <-  getClassDef(exti@superClass, package = exti@package)
            if(is.null(superClass)) {
                warning(gettextf("class %s extends an undefined class, %s",
                                 dQuote(ClassDef@className),
                                 dQuote(what)),
                        domain = NA)
                next
            }
            more <- Recall(superClass, soFar)
            whatMore <- more$label
            if(!all(is.na(match(whatMore, soFar)))) {
                ## elminate classes reachable by more than one path
                ## (This is allowed in the model, however)
                ok <- is.na(match(whatMore, soFar))
                more$depth <- more$depth[ok]
                more$label <- more$label[ok]
                more$ext <- more$ext[ok]
                whatMore <- whatMore[ok]
            }
            if(length(whatMore)) {
                soFar <- c(soFar, whatMore)
                super$depth <- c(super$depth, 1+more$depth)
                super$label <- c(super$label, more$label)
                super$ext <- c(super$ext, more$ext)
            }
    }
    super
}

selectSuperClasses <-
    function(Class, dropVirtual = FALSE, namesOnly = TRUE,
             directOnly = TRUE, simpleOnly = directOnly,
             where = topenv(parent.frame()))
{
    ext <- if(isClassDef(Class))
        Class@contains
    else if(isClass(Class, where = where))
        getClass(Class, where = where)@contains
    else stop("'Class' must be a valid class definition or class")

    .selectSuperClasses(ext, dropVirtual = dropVirtual, namesOnly = namesOnly,
                        directOnly = directOnly, simpleOnly = simpleOnly)
}

.selectSuperClasses <- function(ext, dropVirtual = FALSE, namesOnly = TRUE,
                                directOnly = TRUE, simpleOnly = directOnly)
{
    ## No argument checking here
    addCond <- function(xpr, prev)
        if(length(prev)) substitute(P && N, list(P = prev, N = xpr)) else xpr
    C <- if(dropVirtual) {
        ## NB the default 'where' in getClass() may depend on specific superClass:
        isVirtualExt <- function(x) getClass(x@superClass)@virtual
        quote(!isVirtualExt(exti))
    } else expression()
    if(directOnly) C <- addCond(quote(length(exti@by) == 0), C)
    if(simpleOnly) C <- addCond(quote(exti@simple), C)
    if(length(C)) {
      F <- function(exti){}; body(F) <- C
      ext <- ext[unlist(lapply(ext, F), use.names=FALSE)]
    }
    if(namesOnly) names(ext) else ext
}

inheritedSlotNames <- function(Class, where = topenv(parent.frame()))
{
    ext <- if(isClassDef(Class))
        Class@contains
    else if(isClass(Class, where = where))
        getClass(Class, where = where)@contains
    supcl <- .selectSuperClasses(ext) ## maybe  simpleOnly = FALSE or use as argument?
    unique(unlist(lapply(lapply(supcl, getClassDef), slotNames), use.names=FALSE))
    ## or just the non-simplified part (*with* names):
    ##     lapply(sapply(supcl, getClassDef, simplify=FALSE), slotNames)
}


isVirtualClass <-
  ## Is the named class a virtual class?  A class is virtual if explicitly declared to
  ## be, and also if the class is not formally defined.
  function(Class, where = topenv(parent.frame())) {
      if(isClassDef(Class))
          Class@virtual
      else if(isClass(Class, where = where))
          getClass(Class, where = where)@virtual
      else
          TRUE
  }


assignClassDef <-
  ## assign the definition of the class to the specially named object
  function(Class, def, where = .GlobalEnv, force = FALSE) {
      if(!is(def,"classRepresentation"))
          stop(gettextf("trying to assign an object of class %s as the definition of class %s: must supply a \"classRepresentation\" object",
                        dQuote(class(def)),
                        dQuote(Class)),
               domain = NA)
      clName <- def@className; attributes(clName) <- NULL
      if(!.identC(Class, clName))
          stop(gettextf("assigning as %s a class representation with internal name %s",
                        dQuote(Class),
                        dQuote(def@className)),
               domain = NA)
      where <- as.environment(where)
      mname <- classMetaName(Class)
      if(exists(mname, envir = where, inherits = FALSE) && bindingIsLocked(mname, where)) {
          if(force)
            .assignOverBinding(mname, def, where, FALSE)
          else
            stop(gettextf("class %s has a locked definition in package %s",
                          dQuote(Class), sQuote(getPackageName(where))))
      }
      else
          assign(mname, def, where)
      if(cacheOnAssign(where)) # will be FALSE for sourceEnvironment's
          .cacheClass(clName, def, is(def, "ClassUnionRepresentation"), where)
  }


.InitClassDefinition <- function(where) {
    defSlots <- list(slots = "list", contains = "list", virtual = "logical",
                     prototype = "ANY", validity = "OptionalFunction", access = "list",
                     ## the above are to conform to the API; now some extensions
                     className = "character", package = "character",
                     subclasses = "list", versionKey = "externalptr", ## or "integer"??
                     sealed = "logical")
    ## the prototype of a new class def'n:  virtual class with NULL prototype
    protoSlots <- list(slots=list(), contains=list(), virtual=NA,
                  prototype = NULL, validity = NULL,
                  access = list(), className = character(), package = character(),
                  subclasses = list(), versionKey = .newExternalptr(),
                  sealed = FALSE)
    proto <- defaultPrototype()
    pnames <- names(protoSlots)
    for(i in seq_along(protoSlots))
        slot(proto, pnames[[i]], FALSE) <- protoSlots[[i]]
    classRepClass <- .classNameFromMethods("classRepresentation")
    class(proto) <- classRepClass
    object <- defaultPrototype()
    class(object) <- classRepClass
    slot(object, "slots", FALSE) <- defSlots
    slot(object, "className", FALSE) <- classRepClass
    slot(object, "virtual", FALSE) <- FALSE
    slot(object, "prototype", FALSE) <- proto
    for(what in c("contains", "validity", "access", "hasValidity", "subclasses",
                  "versionKey"))
        slot(object, what, FALSE) <- elNamed(protoSlots, what)
    slot(object, "sealed", FALSE) <- TRUE
    slot(object, "package", FALSE) <- getPackageName(where)
##    assignClassDef("classRepresentation", object, where)
    assign(classMetaName("classRepresentation"), object, where)
    ## the list of needed generics, initially empty (see .InitStructureMethods)
    assign(".NeedPrimitiveMethods", list(), where)
}

.classNameFromMethods <- function(what) {
    packageSlot(what) <- "methods"
    what
  }

.initClassSupport <- function(where) {
    setClass("classPrototypeDef", representation(object = "ANY", slots = "character", dataPart = "logical"),
             sealed = TRUE, where = where)
    setClass(".Other", representation(label = "character"),
             sealed = TRUE, where = where)  # nonvirtual, nobody's subclass, see testInheritedMethods
    ## a class and a method for reporting method selection ambiguities
    setClass("MethodSelectionReport",
         representation(generic = "character", allSelections = "character", target = "character", selected = "character", candidates = "list", note = "character"),
             sealed = TRUE, where = where)
    setClass("classGeneratorFunction",
             representation(className = "character", package = "character"),
             contains = "function")
}


newBasic <-
  ## the implementation of the function `new' for basic classes.
  ##
  ## See `new' for the interpretation of the arguments.
  function(Class, ...) {
      msg <- NULL
      value <- switch(Class,
               "NULL" = return(NULL), ## can't set attr's of NULL in R
               "logical" =,
               "numeric" =,
               "character" =,
               "complex" =,
               "integer" =,
               "raw" =,
               "list" =  as.vector(c(...), Class),
               "expression" = eval(substitute(expression(...))),
               "externalptr" = {
                   if(nargs() > 1)
                       stop("'externalptr' objects cannot be initialized from new()")
                   .newExternalptr()
               },
               "single" = as.single(c(...)),
                  ## note on array, matrix:  not possible to be compatible with
                  ## S-Plus on array, unless R allows 0-length .Dim attribute
               "array" = if(!missing(...)) array(...) else structure(numeric(), .Dim =0L),
               "matrix" = if (!missing(...)) matrix(...) else matrix(0, 0L, 0L),
#               "ts" = ts(...),
# break dependence on package stats
	       "ts" = if(!missing(...)) stats::ts(...) else
		      structure(NA, .Tsp = c(1, 1, 1), class = "ts"),

                ## otherwise:
                  {
                      args <- list(...)
                      if(length(args) == 1L && is(args[[1L]], Class)) {
                          value <- as(args[[1L]], Class)
                      }
                      else if(is.na(match(Class, .BasicClasses)))
                          msg <- paste("Calling new() on an undefined and non-basic class (\"",
                               Class, "\")", sep="")
                      else
                          msg <-
                              gettextf("initializing objects from class %s with these arguments is not supported",
                                       dQuote(Class))
                  }
                  )
  if(is.null(msg))
      value
  else
      stop(msg, domain = NA)
}


## this non-exported function turns on or off
## the use of the S4 type as class prototype
.useS4Prototype <- function(on = TRUE, where  = .methodsNamespace) {
    if(on)
     pp <- .Call(C_Rf_allocS4Object)
    else
     pp <-  list()
    .assignOverBinding(".defaultPrototype", where=where, pp, FALSE)
}

defaultPrototype <-
    ## the starting prototype for a non-virtual class
    ## Should someday be a non-vector sexp type
    function()
    .defaultPrototype

reconcilePropertiesAndPrototype <-
  ## makes a list or a structure look like a prototype for the given class.
  ##
  ## Specifically, returns a structure with attributes corresponding to the slot
  ## names in properties and values taken from prototype if they exist there, from
  ## `new(classi)' for the class, `classi' of the slot if that succeeds, and `NULL'
  ## otherwise.
  ##
  ## The prototype may imply slots not in the properties list.  It is not required that
    ## the extends classes be define at this time.  Should it be?
  function(name, properties, prototype, superClasses, where) {
      ## the StandardPrototype should really be a type that doesn't behave like
      ## a vector.  But none of the existing SEXP types work.  Someday ...
      StandardPrototype <- defaultPrototype()
      slots <-  validSlotNames(allNames(properties))
      dataPartClass <- elNamed(properties, ".Data")
      dataPartValue <- FALSE
      if(!is.null(dataPartClass) && is.null(.validDataPartClass(dataPartClass, where)))
          stop(gettextf("in defining class %s, the supplied data part class, %s is not valid (must be a basic class or a virtual class combining basic classes)",
                        dQuote(name), dQuote(dataPartClass)),
               domain = NA)
      prototypeClass <- getClass(class(prototype), where = where)
      if((!is.null(dataPartClass) || length(superClasses))
         && is.na(match("VIRTUAL", superClasses))) {
          ## Look for a data part in the superclasses, either an inherited
          ## .Data slot, or a basic class.  Uses the first possibility, warns of conflicts
          for(cl in superClasses) {
              clDef <- getClassDef(cl, where = where)
              if(is.null(clDef))
                stop(gettextf("no definition was found for superclass %s in the specification of class %s",
                              dQuote(cl), dQuote(name)),
                     domain = NA)
              thisDataPart <-  .validDataPartClass(clDef, where, dataPartClass)
              if(!is.null(thisDataPart)) {
                    dataPartClass <- thisDataPart
                    if(!is.null(clDef@prototype)) {
                      newObject <- clDef@prototype
                      dataPartValue <- TRUE
                    }
                  }
          }
          if(length(dataPartClass)) {
              if(is.na(match(".Data", slots))) {
                  properties <- c(list(".Data"= dataPartClass), properties)
                  slots <- names(properties)
              }
              else if(!extends(elNamed(properties, ".Data"), dataPartClass))
                  stop(gettextf("conflicting definition of data part: .Data = %s, superclass implies %s",
                                dQuote(elNamed(properties, ".Data")),
                                dQuote(dataPartClass)),
                       domain = NA)
              pslots <- NULL
              if(is.null(prototype)) {
                  if(dataPartValue)
                      prototype <- newObject
                  else if(isVirtualClass(dataPartClass, where = where))
                      ## the equivalent of new("vector")
                      prototype <- newBasic("logical")
                  else
                      prototype <- new(dataPartClass)
                  prototypeClass <- getClass(class(prototype), where = where)
              }
              else {
                  if(extends(prototypeClass, "classPrototypeDef")) {
                      hasDataPart <- identical(prototype@dataPart, TRUE)
                      if(!hasDataPart) {
                          if(!dataPartValue) # didn't get a .Data object
                            newObject <- new(dataPartClass)
                          pobject <- prototype@object
                          ## small amount of head-standing to preserve
                          ## any attributes in newObject & not in pobject
                          anames <- names(attributes(pobject))
                          attributes(newObject)[anames] <- attributes(pobject)
                          prototype@object <- newObject
                      }
                      else if(!extends(getClass(class(prototype@object), where = where)
                                       , dataPartClass))
                          stop(gettextf("a prototype object was supplied with object slot of class %s, but the class definition requires an object that is class %s",
                                        dQuote(class(prototype@object)),
                                        dQuote(dataPartClass)),
                               domain = NA)
                  }
                  else if(!extends(prototypeClass, dataPartClass))
                      stop(gettextf("a prototype was supplied of class %s, but the class definition requires an object that is class %s",
                                    dQuote(class(prototype)),
                                    dQuote(dataPartClass)),
                           domain = NA)
              }
          }
          if(is.null(prototype)) { ## non-vector (may extend NULL)
              prototype <- StandardPrototype
          }
      }
      ## check for conflicts in the slots
      allProps <- properties
      for(i in seq_along(superClasses)) {
          cl <- superClasses[[i]]
          clDef <- getClassDef(cl, where)
          if(is(clDef, "classRepresentation")) {
              theseProperties <- getSlots(clDef)
              theseSlots <- names(theseProperties)
              theseSlots <- theseSlots[theseSlots == ".Data"] # handled already
              dups <- !is.na(match(theseSlots, allProps))
              for(dup in theseSlots[dups])
                  if(!extends(elNamed(allProps, dup), elNamed(theseProperties, dup)))
                      stop(gettextf("slot %s in class %s currently defined (or inherited) as \"%s\", conflicts with an inherited definition in class %s",
                                    sQuote(dup),
                                    dQuote(name),
                                    elNamed(allProps, dup),
                                    dQuote(cl)),
                           domain = NA)
              theseSlots <- theseSlots[!dups]
              if(length(theseSlots))
                  allProps[theseSlots] <- theseProperties[theseSlots]
          }
          else
              stop(gettextf("class %s extends an undefined class (%s)",
                            dQuote(name), dQuote(cl)),
                   domain = NA)
      }
      if(is.null(dataPartClass)) {
          if(extends(prototypeClass, "classPrototypeDef"))
          {}
          else {
              if(is.list(prototype))
               prototype <- do.call("prototype", prototype)
              if(is.null(prototype))
                  prototype <- StandardPrototype
          }
      }
      else {
          dataPartDef <- getClass(dataPartClass)
          checkDataPart <- !isXS3Class(dataPartDef)
          if(checkDataPart)
            checkDataPart  <-
              ((is.na(match(dataPartClass, .BasicClasses)) &&
                !isVirtualClass(dataPartDef)) || length(dataPartDef@slots))
          if(checkDataPart)
              stop(gettextf("%s is not eligible to be the data part of another class (must be a basic class or a virtual class with no slots)",
                            dQuote(dataPartClass)),
                   domain = NA)
          if(extends(prototypeClass, "classPrototypeDef"))
          {}
          else if(extends(prototypeClass, dataPartClass)) {
              if(extends(prototypeClass, "list") && length(names(prototype)))
                  warning("prototype is a list with named elements (could be ambiguous):  better to use function prototype() to avoid trouble.")
          }
          else if(is.list(prototype))
              prototype <- do.call("prototype", prototype)
      }
      ## pnames will be the names explicitly defined in the prototype
      if(extends(prototypeClass, "classPrototypeDef")) {
          pnames <- prototype@slots
          prototype <- prototype@object
          if(length(superClasses) == 0L && any(is.na(match(pnames, slots))))
              stop(sprintf(ngettext(sum(is.na(match(pnames, slots))),
                                    "named elements of prototype do not correspond to slot name: %s",
                                    "named elements of prototype do not correspond to slot names: %s"),
                           paste(.dQ(pnames[is.na(match(pnames, slots))]),
                                 collapse =", ")),
                   domain = NA)
      }
      else
          pnames <- allNames(attributes(prototype))
       ## now set the slots not yet in the prototype object.
      ## An important detail is that these are
      ## set using slot<- with check=FALSE (because the slot will not be there already)
      ## what <- is.na(match(slots, pnames))
      what <- seq_along(properties)
      props <- properties[what]
      what <- slots[what]
      nm <- names(attributes(prototype))
      for(i in seq_along(what)) {
          propName <- el(what, i)
          if(!identical(propName, ".Data") && !propName %in% nm)
#             is.null(attr(prototype, propName)))
              slot(prototype, propName, FALSE) <- tryNew(el(props, i), where)
      }
      list(properties = properties, prototype = prototype)
  }

tryNew <-
    ## Tries to generate a new element from this class, but if
    ## the class is undefined just returns NULL.
    ##
    ## For virtual classes, returns the class prototype
    ## so that the object is valid member of class.
    ## Otherwise tries to generate a new() object, but in rare
    ## cases, this might fail if the install() method required
    ## an argument, so this case is trapped as well.
  function(Class, where)
{
    ClassDef <- getClassDef(Class, where)
    if(is.null(ClassDef))
        return(NULL)
    else if(identical(ClassDef@virtual, TRUE))
        ClassDef@prototype
    else tryCatch(new(ClassDef),
                  error = function(e) {
                      value <- ClassDef@prototype
                      class(value) <- ClassDef@className
                      value
                  })
}

empty.dump <- function() list()

isClassDef <- function(object) is(object, "classRepresentation")

showClass <-
    ## print the information about a class definition.
    ## If complete==TRUE, include the indirect information about extensions.
    function(Class, complete = TRUE, propertiesAreCalled = "Slots")
{
    if(isClassDef(Class)) {
        ClassDef <- Class
        Class <- ClassDef@className
    }
    else if(complete)
        ClassDef <- getClass(Class)
    else
        ClassDef <- getClassDef(Class)
    cat(if(identical(ClassDef@virtual, TRUE)) "Virtual ",
	"Class ", .dQ(Class),
	## Show the package if that is non-trivial:
	if(nzchar(pkg <- ClassDef@package))
	c(" [", if(pkg != ".GlobalEnv") "package" else "in", " \"", pkg,"\"]"),
	"\n", sep="")
    x <- ClassDef@slots
    if(length(x)) {
        printPropertiesList(x, propertiesAreCalled)
    }
    else
        cat("\nNo ", propertiesAreCalled, ", prototype of class \"",
            .class1(ClassDef@prototype), "\"\n", sep="")
    ext <- ClassDef@contains
    if(length(ext)) {
        cat("\nExtends: ")
        showExtends(ext)
    }
    ext <- ClassDef@subclasses
    if(length(ext)) {
        cat("\nKnown Subclasses: ")
        showExtends(ext)
    }
}

printPropertiesList <- function(x, propertiesAreCalled) {
    if(length(x)) {
        n <- length(x)
        cat("\n",propertiesAreCalled, ":\n", sep="")
        text <- format(c(names(x), as.character(x)), justify="right")
        text <- matrix(text, nrow = 2L, ncol = n, byrow = TRUE)
        dimnames(text) <- list(c("Name:", "Class:"), rep.int("", n))
        print(text, quote = FALSE)
    }
}

showExtends <-
    ## print the elements of the list of extensions.  Also used to print
    ## extensions recorded in the opposite direction, via a subclass list
    function(ext, printTo = stdout())
{
    what <- names(ext)
    how <- character(length(ext))
    for(i in seq_along(ext)) {
        eli <- el(ext, i)
        if(is(eli, "SClassExtension")) {
            how[i] <-
                if(length(eli@by))
		    paste("by class", paste0("\"", eli@by, "\", distance ",
					     eli@distance, collapse = ", "))
                else if(identical(eli@dataPart, TRUE))
                    "from data part"
                else "directly"
            if(!eli@simple) {
                if(is.function(eli@test) && !identical(body(eli@test), TRUE)) {
                    how[i] <-
                        paste(how[i], if(is.function(eli@coerce))
                              ", with explicit test and coerce" else
                              ", with explicit test", sep="")
                }
                else if(is.function(eli@coerce))
                    how[i] <- paste0(how[i], ", with explicit coerce")
            }
        }
    }
    if(identical(printTo, FALSE))
        list(what = what, how = how)
    else if(all(!nzchar(how)) ||  all(how == "directly")) {
        what <- paste0('"', what, '"')
        if(length(what) > 1L)
            what <- c(paste0(what[-length(what)], ","), what[[length(what)]])
        cat(file = printTo, what, fill=TRUE)
    }
    else cat(file = printTo, "\n",
	     paste0("Class \"", what, "\", ", how, "\n"), sep = "")
}



printClassRepresentation <-
  function(x, ...)
  showClass(x, propertiesAreCalled="Slots")

## bootstrap definition to be used before getClass() works
possibleExtends <- function(class1, class2, ClassDef1, ClassDef2)
    .identC(class1, class2) || .identC(class2, "ANY")

## "Real" definition (assigned in ./zzz.R )
.possibleExtends <-
    ## Find the information that says whether class1 extends class2,
    ## directly or indirectly.  This can be either a logical value or
    ## an object containing various functions to test and/or coerce the relationship.
    ## TODO:  convert into a generic function w. methods WHEN dispatch is really fast!
    function(class1, class2, ClassDef1 = getClassDef(class1),
             ClassDef2 = getClassDef(class2, where = .classEnv(ClassDef1)))
{
    if(.identC(class1[[1L]], class2) || .identC(class2, "ANY"))
        return(TRUE)
    ext <- TRUE # may become a list of extends definitions
    if(is.null(ClassDef1)) # class1 not defined
        return(FALSE)
    ## else
    ext <- ClassDef1@contains
    nm1 <- names(ext)
    i <- match(class2, nm1)
    if(is.na(i)) {
        ## look for class1 in the known subclasses of class2
        if(!is.null(ClassDef2)) {
            ext <- ClassDef2@subclasses
            ## check for a classUnion definition, not a plain "classRepresentation"
            if(!.identC(class(ClassDef2), "classRepresentation") &&
               isClassUnion(ClassDef2))
                ## a simple TRUE iff class1 or one of its superclasses belongs to the union
		i <- as.logical(anyDuplicated(c(class1, unique(nm1),
						names(ext))))
            else {
                ## class1 could be multiple classes here.
                ## I think we want to know if any extend
                i <- match(class1, names(ext))
                ii <- i[!is.na(i)]
                i <- if(length(ii))  ii[1L] else i[1L]
            }
        }
    }
    if(is.na(i))
        FALSE
    else if(is.logical(i))
        i
    else
        el(ext, i)
}

  ## complete the extends information in the class definition, by following
  ## transitive chains.
  ##
  ## Elements in the immediate extends list may be added and current elements may be
  ## replaced, either by replacing a conditional relation with an unconditional
  ## one, or by adding indirect relations.
  ##
completeExtends <-    function(ClassDef, class2, extensionDef, where) {
    ## check for indirect extensions => already completed
    ext <- ClassDef@contains
    for(i in seq_along(ext)) {
        if(.isIndirectExtension(ext[[i]])) {
            ClassDef <- .uncompleteClassDefinition(ClassDef, "contains")
            break
        }
    }
    exts <- .walkClassGraph(ClassDef, "contains", where, attr(ext, "conflicts"))
    if(length(exts)) {
##         ## sort the extends information by depth (required for method dispatch)
##         superClassNames <- getAllSuperClasses(ClassDef, FALSE)
##         ## FIXME:  getAllSuperClassses sometimes misses.  Why?
##         if(length(superClassNames) == length(exts))
##             exts <- exts[superClassNames]
        if("oldClass" %in% names(exts) &&
           length(ClassDef@slots) > 1L) # an extension of an S3 class
          exts <- .S3Extends(ClassDef, exts, where)
    }
    if(!missing(class2) && length(ClassDef@subclasses)) {
        strictBy <- TRUE # FIXME:  would like to make this conditional but a safe condition is unknown
        subclasses <-
            .transitiveSubclasses(ClassDef@className, class2, extensionDef, ClassDef@subclasses, strictBy)
        ## insert the new is relationship, but without any recursive completion
        ## (asserted not to be needed if the subclass slot is complete)
        for(i in seq_along(subclasses)) {
            obji <- subclasses[[i]]
            ## don't override existing relations
            ## TODO:  have a metric that picks the "closest" relationship
            if(!extends(obji@subClass, class2))
                setIs(obji@subClass, class2, extensionObject = obji, doComplete = FALSE,
                      where = where)
        }
    }
## TODO:  move these checks to a tool used by check & conditional on no .S3Class slot
##     S3Class <- attr(ClassDef@prototype, ".S3Class")
##     if(!is.null(S3Class)) {
##       others <- c(ClassDef@className, names(exts))
##       others <- others[is.na(match(others, S3Class))]
##       if(length(others)>0)
##         .checkS3forClass(ClassDef@className, where, others)
##     }
    exts
}

completeSubclasses <-
    function(classDef, class2, extensionDef, where, classDef2 = getClassDef(class2, where)) {
    ## check for indirect extensions => already completed
    ext <- classDef@subclasses
    for(i in seq_along(ext)) {
        if(.isIndirectExtension(ext[[i]])) {
            classDef <- .uncompleteClassDefinition(classDef, "subclasses")
            break
        }
    }
    subclasses <- .walkClassGraph(classDef, "subclasses", where)
    if(!missing(class2) && length(classDef@contains)) {
        strictBy <- TRUE
        contains <-
            .transitiveExtends(class2, classDef@className, extensionDef, classDef@contains, strictBy)
        ## insert the new is relationship, but without any recursive completion
        ## (asserted not to be needed if the subclass slot is complete)
        for(i in seq_along(contains)) {
            obji <- contains[[i]]
            cli <- contains[[i]]@superClass
            cliDef <- getClassDef(cli, where)
            ## don't override existing relations
            ## TODO:  have a metric that picks the "closest" relationship
            if(!extends(classDef2, cliDef))
                setIs(class2, cli, extensionObject = obji,
                      doComplete = FALSE, where = where)
        }
    }
    subclasses
}


## utility function to walk the graph of super- or sub-class relationships
## in order to incorporate indirect relationships
.walkClassGraph <-  function(ClassDef, slotName, where,  conflicts = character())
{
    ext <- slot(ClassDef, slotName)
    if(length(ext) == 0)
        return(ext)
    className <- ClassDef@className
    ## the super- vs sub-class is identified by the slotName
    superClassCase <- identical(slotName, "contains")
    what <- names(ext)
    for(i in seq_along(ext)) { # note that this loops only over the original ext
        by <- what[[i]]
        if(isClass(by, where = where)) {
            byDef <- getClass(by, where = where)
            exti <-  slot(byDef, slotName)
            coni <- attr(exti, "conflicts") # .resolveSuperclasses makes this
            if(superClassCase && length(coni) > 0) {
                conflicts <- unique(c(conflicts, coni))
              }
            ## add in those classes not already known to be super/subclasses
            exti <- exti[is.na(match(names(exti), what))]
            if(length(exti)) {
                if(superClassCase) {
                    strictBy <- TRUE  # FIXME:  need to find some safe test allowing non-strict
                      exti <- .transitiveExtends(className, by, ext[[i]], exti, strictBy)
                }
                else {
                    strictBy <- TRUE
                    exti <- .transitiveSubclasses(by, className, ext[[i]], exti, strictBy)
                }
                ext <- c(ext, exti)
            }
        }
        else
            stop(gettextf("the '%s' list for class %s, includes an undefined class %s",
                          if(superClassCase) "superClass" else "subClass",
                          dQuote(className),
                          dQuote(.className(by))),
                 domain = NA)
    }
    what <- names(ext)  ## the direct and indirect extensions
    if(!all(is.na(match(what, className)))) {
        ok <- is.na(match(what, className))
        ## A class may not contain itself, directly or indirectly
        ## but a non-simple cyclic relation, involving setIs, is allowed
        for(i in seq_along(what)[!ok]) {
            exti <- ext[[i]]
            if(!is(exti, "conditionalExtension")) {
                if(superClassCase) {
                    whatError <-  "contain itself"
                }
                else {
                    whatError <- "have itself as a subclass"
                }
                ## this is not translatable
                stop(sprintf("class %s may not %s: it contains class %s, with a circular relation back to %s",
                             dQuote(className), whatError,
                             dQuote(exti@by),
                             dQuote(className)),
                     domain = NA)
            }
        }
        ext <- ext[ok]
    }
    ## require superclasses to be sorted by distance
    distOrder <- sort.list(sapply(ext, function(x)x@distance))
    ext <- ext[distOrder]
    if(superClassCase && (anyDuplicated(what) || length(conflicts) > 0))
        ext <- .resolveSuperclasses(ClassDef, ext, where, conflicts)
    ext
}

.reportSuperclassConflicts <- function(className, ext, where) {
    what <- names(ext)
    conflicts <- character()
    for(i in seq_along(ext)) {
        by <- what[[i]]
        ## report only the direct superclass from which inconsistencies are inherited
        if(identical(ext[[i]]@distance, 1) && isClass(by, where = where)) {
            byDef <- getClass(by, where = where)
            exti <-  byDef@contains
            coni <- attr(exti, "conflicts") # .resolveSuperclasses makes this
            if( length(coni) > 0) {
                warning(gettextf("class %s is inheriting an inconsistent superclass structure from class %s, inconsistent with %s",
                                 .dQ(className), .dQ(by),
                                 paste(.dQ(coni), collapse = ", ")),
                        call. = FALSE, domain = NA)
                conflicts <- unique(c(conflicts, coni))
              }
          }
      }
          newconflicts <- attr(ext, "conflicts")
        if(length(newconflicts) > length(conflicts))
          warning(gettextf("unable to find a consistent ordering of superclasses for class %s: order chosen is inconsistent with the superclasses of %s",
                           .dQ(className),
                           paste(.dQ(setdiff(newconflicts, conflicts)),
                                 collapse = ", ")),
                  call. = FALSE, domain = NA)
        }


.resolveSuperclasses <- function(classDef, ext, where, conflicts = attr(ext, "conflicts")) {
  ## find conditional extensions, ignored in superclass ordering
  .condExts <- function(contains)
      sapply(contains, function(x) is(x, "conditionalExtension" ))
  .noncondExtsClass <- function(cl) {
    if(isClass(cl, where = where) ) {
      contains <- getClass(cl, where = where)@contains
      names(contains)[!.condExts(contains)]
    }
    else cl
  }
  what <- names(ext)
  dups <- unique(what[duplicated(what)])
  if(length(dups) > 0) {
    ## First, eliminate all conditional relations, which never override non-conditional
    affected <- match(what, dups, 0) > 0
    conditionals <- .condExts(ext)
    if(any(conditionals)) {
      affected[conditionals] <- FALSE
      what2 <- what[affected]
      dups <- unique(what2[duplicated(what2)])
      if(length(dups) == 0) {
        ##  eliminating conditonal relations removed duplicates
        if(length(conflicts) > 0)
          attr(ext, "conflicts") <- unique(c(conflicts, attr(ext, "conflicts")))
        return(ext)
      }
      ## else, go on with conditionals eliminated
    }
    directSupers <- sapply(classDef@contains, function(x) identical(x@distance, 1))
    directSupers <- unique(names(classDef@contains[directSupers]))
    ## form a list of the superclass orderings of the direct superclasses
    ## to check consistency with each way to eliminate duplicates
    ## Once again, conditional relations are eliminated
    superExts <- lapply(directSupers, .noncondExtsClass)
    names(superExts) <- directSupers
    retain = .choosePos(classDef@className, what, superExts, affected)
    if(is.list(retain)) {
      these <- retain[[2]]
      conflicts <- unique(c(conflicts, these)) # append the new conflicts
      retain <- retain[[1]]
    }
    ## eliminate the affected & not retained
    affected[retain] <- FALSE
    ext <- ext[!affected]
  }
  ## even if no dups here, may have inherited some conflicts,
  ## which will be copied to the contains list.
  ## FUTURE NOTE (7/09):  For now, we are using an attribute for conflicts,
  ## rather than promoting the ext list to a new class, which may be desirable
  ## if other code comes to depend on the conflicts information.
  attr(ext, "conflicts") <- conflicts
  ext
}

classMetaName <-
  ## a name for the object storing this class's definition
  function(name)
  methodsPackageMetaName("C", name)

# regexp for matching class metanames; semi-general but assumes the
# meta pattern starts with "." and has no other special characters
.ClassMetaPattern <- function()
    paste0("^[.]",substring(methodsPackageMetaName("C",""),2))

##FIXME:  C code should take multiple strings in name so paste() calls could  be avoided.
methodsPackageMetaName <-
  ## a name mangling device to simulate the meta-data in S4
  function(prefix, name, package = "")
  ## paste(".", prefix, name, sep="__") # too slow
    .Call(C_R_methodsPackageMetaName, prefix, name, package)

## a  non-exported regexp that matches  methods metanames
## This is quite general and matches all patterns that could be generated
## by calling methodsPackageMetaName() with a sequence of capital Latin letters
## Used by package.skeleton in utils
.methodsPackageMetaNamePattern <- "^[.]__[A-Z]+__"

requireMethods <-
  ## Require a subclass to implement methods for the generic functions, for this signature.
  ##
  ## For each generic, `setMethod' will be called to define a method that throws an error,
  ## with the supplied message.
  ##
  ## The `requireMethods' function allows virtual classes to require actual classes that
  ## extend them to implement methods for certain functions, in effect creating an API
  ## for the virtual class.  Otherwise, default methods for the corresponding function would
  ## be called, resulting in less helpful error messages or (worse still) silently incorrect
  ## results.
  function(functions, signature,
           message = "", where = topenv(parent.frame()))
{
    for(f in functions) {
        method <- getMethod(f, optional = TRUE)
        if(!is.function(method))
            method <- getGeneric(f, where = where)
        body(method) <- substitute(stop(methods:::.missingMethod(FF, MESSAGE, if(exists(".Method")).Method else NULL), domain=NA), list(FF=f, MESSAGE=message))
        environment(method) <- .GlobalEnv
        setMethod(f, signature, method, where = where)
    }
    NULL
}

## Construct an error message for an unsatisfied required method.
.missingMethod <- function(f, message = "", method) {
    if(nzchar(message))
        message <- paste0("(", message, ")")
    message <- paste("for function", f, message)
    if(is(method, "MethodDefinition")) {
        target <-  paste(.dQ(method@target), collapse=", ")
        defined <- paste(.dQ(method@defined), collapse=", ")
        message <- paste("Required method", message, "not defined for signature",
                         target)
        if(!identical(target, defined))
            message <- paste(message, ", required for signature", defined)
    }
    else message <- paste("Required method not defined", message)
    message
}

getSlots <- function(x) {
    classDef <- if(isClassDef(x)) x else getClass(x)
    props <- classDef@slots
    value <- as.character(props)
    names(value) <- names(props)
    value
}


## check for reserved slot names.  Currently only "class" is reserved
validSlotNames <- function(names) {
    if(is.na(match("class", names)))
        names
    else
        stop("\"class\" is a reserved slot name and cannot be redefined")
}

### utility function called from primitive code for "@"
getDataPart <- function(object) {
    if(identical(typeof(object),"S4")) {
        ## explicit .Data or .xData slot
        ## Some day, we may merge both of these as .Data
        value <- attr(object, ".Data")
        if(is.null(value)) {
            value <- attr(object, ".xData")
            if(is.null(value))
              stop("Data part is undefined for general S4 object")
          }
        if(identical(value, .pseudoNULL))
          return(NULL)
        else
          return(value)
    }
    temp <- getClass(class(object))@slots
    if(length(temp) == 0L)
        return(object)
    if(is.na(match(".Data", names(temp))))
       stop(gettextf("no '.Data' slot defined for class %s",
                     dQuote(class(object))),
            domain = NA)
    dataPart <- temp[[".Data"]]
    switch(dataPart,
           ## the common cases, for efficiency
           numeric = , vector = , integer = , character = , logical = ,
           complex = , list =
              attributes(object) <- NULL,
           matrix = , array = {
               value <- object
               attributes(value) <- NULL
               attr(value, "dim") <- attr(object, "dim")
               attr(value, "dimnames") <- attr(object, "dimnames")
               object <- value
           },
           ts = {
               value <- object
               attributes(value) <- NULL
               attr(value, "ts") <- attr(object, "ts")
               object <- value
           },
           ## default:
           if(is.na(match(dataPart, .BasicClasses))) {
               ## keep attributes not corresponding to slots
               attrVals <- attributes(object)
               attrs <- names(attrVals)
               attrs <- attrs[is.na(match(attrs, c("class", names(temp))))]
               attributes(object) <- attrVals[attrs]
           }
           else
           ## other basic classes have no attributes
               attributes(object) <- NULL
           )
    object
}

setDataPart <- function(object, value, check = TRUE) {
    if(check || identical(typeof(object), "S4")) {
        classDef <- getClass(class(object))
        slots <- getSlots(classDef)
        dataSlot <- .dataSlot(names(slots))
        if(length(dataSlot) == 1)
          dataClass <- elNamed(slots, dataSlot)
        else if(check)
          stop(gettextf("class %s does not have a data part (a .Data slot) defined",
                        dQuote(class(object))),
               domain = NA)
        else # this case occurs in making the methods package. why?
          return(.mergeAttrs(value, object))
        value <- as(value, dataClass)  # note that this is strict as()
        if(identical(typeof(object), "S4")) {
            if(is.null(value))
              value <- .pseudoNULL
            attr(object, dataSlot) <- value
            return(object)
        }
    }
    .mergeAttrs(value, object)
}

.validDataPartClass <- function(cl, where, prevDataPartClass = NULL) {
    if(is(cl, "classRepresentation")) {
        ClassDef <- cl
        cl <- ClassDef@className
    }
    else
        ClassDef <- getClass(cl, TRUE)

    switch(cl, matrix = , array = value <- cl,
           value <- elNamed(ClassDef@slots, ".Data"))
    if(is.null(value)) {
        if(.identC(cl, "structure"))
            value <- "vector"
        else if((extends(cl, "vector") || !is.na(match(cl, .BasicClasses))))
            value <- cl
        else if(extends(cl, "oldClass") && isVirtualClass(cl)) {
        }
        else if(identical(ClassDef@virtual, TRUE) &&
               length(ClassDef@slots) == 0L &&
               length(ClassDef@subclasses) ) {
                ## look for a union of basic classes
                subclasses <- ClassDef@subclasses
                what <- names(subclasses)
                value <- cl
                for(i in seq_along(what)) {
                    ext <- subclasses[[i]]
                    ##TODO:  the following heuristic test for an "original"
                    ## subclass should be replaced by a suitable class (extending SClassExtension)
                    if(length(ext@by) == 0L && ext@simple && !ext@dataPart &&
                       is.na(match(what[i], .BasicClasses))) {
                        value <- NULL
                        break
                    }
                }
            }
    }
    if(!(is.null(value) || is.null(prevDataPartClass) || extends(prevDataPartClass, value) ||
         isVirtualClass(value, where = where))) {
      warning(gettextf("more than one possible class for the data part: using %s rather than %s",
                  .dQ(prevDataPartClass), .dQ(value)), domain = NA)
      value <- NULL
    }
    value
}

.dataSlot <- function(slotNames) {
    dataSlot <- c(".Data", ".xData")
    dataSlot <- dataSlot[match(dataSlot, slotNames, 0)>0]
    if(length(dataSlot) > 1)
      stop("class cannot have both an ordinary and hidden data type")
    dataSlot
  }


.mergeAttrs <- function(value, object, explicit = NULL) {
    supplied <- attributes(object)
    if(length(explicit))
        supplied[names(explicit)] <- explicit
    valueAttrs <- attributes(value)
    ## names are special.
    if(length(supplied$names) && length(valueAttrs$names) == 0L) {
        if(length(value) != length(object))
            length(supplied$names) <- length(value)
    }
    if(length(valueAttrs)) {	 ## don't overwrite existing attrs
	valueAttrs$class <- NULL ## copy in class if it's supplied
	supplied[names(valueAttrs)] <- valueAttrs
    } ## else --  nothing to protect
    attributes(value) <- supplied
    if(isS4(object))
        .asS4(value)
    else
        value
}

.newExternalptr <- function()
    .Call(C_R_externalptr_prototype_object)

## modify the list moreExts, currently from class `by', to represent
## extensions instead from an originating class; byExt is the extension
## from that class to `by'
.transitiveExtends <- function(from, by, byExt, moreExts, strictBy) {
    what <- names(moreExts)
###    if(!strictBy) message("Extends: ",from, ": ", paste(what, collapse = ", "))
    for(i in seq_along(moreExts)) {
        toExt <- moreExts[[i]]
        to <- what[[i]]
        toExt <- .combineExtends(byExt, toExt, by, to, strictBy)
        moreExts[[i]] <- toExt
    }
    moreExts
###    if(!strictBy) message("Done")
}

.transitiveSubclasses <- function(by, to, toExt, moreExts, strictBy) {
    what <- names(moreExts)
###    if(!strictBy) message("Subclasses: ",by, ": ", paste(what, collapse = ", "))
    for(i in seq_along(moreExts)) {
        byExt <- moreExts[[i]]
        byExt <- .combineExtends(byExt, toExt, by, to, strictBy)
        moreExts[[i]] <- byExt
    }
    moreExts
###    if(!strictBy) message("Done")
}

.combineExtends <- function(byExt, toExt, by, to, strictBy) {
        ## construct the composite coerce method, taking into account the strict=
        ## argument.
        f <- toExt@coerce
        fR <- toExt@replace
            toExpr <- body(f)
            fBy <- byExt@coerce
            byExpr <- body(fBy)
        ## if both are simple extensions, so is the composition
        if(byExt@simple && toExt@simple) {
            expr <- (if(byExt@dataPart)
                     substitute({if(strict) from <- from@.Data; EXPR},
                                list(EXPR = toExpr))
                   else if(toExt@dataPart)
                     substitute({from <- EXPR;  if(strict) from@.Data},
                                list(EXPR = byExpr))
                   else  (if(identical(byExpr, quote(from)) && identical(toExpr, quote(from)))
                           quote(from)
                         else
                           substitute({from <- E1; E2}, list(E1 = byExpr, E2 = toExpr))
                         )
                     )
            body(f, envir = environment(f)) <- expr
        }
        else {
            toExt@simple <- FALSE
            if(!identical(byExpr, quote(from)))
                body(f, envir = environment(f)) <-
                    substitute( {from <- as(from, BY, strict = strict); TO},
                               list(BY = by, TO = toExpr))
        }
        toExt@coerce <- f
        f <- toExt@test
        toExpr <- body(f)
        byExpr <- body(byExt@test)
        ## process the test code
        if(!identical(byExpr, TRUE)) {
            if(!identical(toExpr, TRUE))
                body(f, envir = environment(f)) <- substitute((BY) && (TO),
                              list(BY = byExpr, TO = toExpr))
            else
                body(f, envir = environment(f)) <- byExpr
        }
        toExt@test <- f
        f <- byExt@replace
        byExpr <- body(f)
        if(!strictBy) {
            toDef <- getClassDef(to)
            byDef <- getClassDef(by)
            strictBy <- is.null(toDef) || is.null(byDef) || toDef@virtual || byDef@virtual
        }
        ## Is there a danger of infinite loop below?
        expr <- substitute({.value <- as(from, BY, STRICT); as(.value, TO) <- value; value <- .value; BYEXPR},
                           list(BY=by, TO = to, BYEXPR = byExpr, STRICT = strictBy))
        body(f, envir = environment(f)) <- expr
        toExt@replace <- f
        toExt@by <- toExt@subClass
        toExt@subClass <- byExt@subClass
        toExt@distance <- toExt@distance + byExt@distance
        ## the combined extension is conditional if either to or by is conditional
        if(is(byExt, "conditionalExtension") && !is(toExt, "conditionalExtension"))
          class(toExt) <- class(byExt)
        toExt
}

## construct the expression that implements the computations for coercing
## an object to one of its superclasses
## The fromSlots argument is provided for calls from makeClassRepresentation
## and completeClassDefinition,
## when the fromClass is in the process of being defined, so slotNames() would fail
.simpleCoerceExpr <- function(fromClass, toClass, fromSlots, toDef) {
    toSlots <- names(toDef@slots)
    sameSlots <- (length(fromSlots) == length(toSlots) &&
		  !any(is.na(match(fromSlots, toSlots))))
    if(!isVirtualClass(toDef))
        toClass <- class(new(toDef)) # get it with the package slot correct
    if(sameSlots)
	substitute({class(from) <- CLASS; from}, list(CLASS = toClass))
    else if(length(toSlots) == 0L) {
	## either a basic class or something with the same representation
	if(is.na(match(toClass, .BasicClasses)))
	    substitute({ attributes(from) <- NULL; class(from) <- CLASS; from},
		       list(CLASS = toClass))
	else if(isVirtualClass(toDef))
	    quote(from)
	else {
	    ## a basic class; a vector type, matrix, array, or ts
	    switch(toClass,
		   matrix = , array = {
		       quote({.dm <- dim(from); .dn <- dimnames(from)
			      attributes(from) <- NULL; dim(from) <- .dm
			      dimnames(from) <- .dn; from})
		   },
		   ts = {
		       quote({.tsp <- tsp(from); attributes(from) <- NULL
			      tsp(from) <- .tsp; class(from) <- "ts"; from})
		   },
		   quote({attributes(from) <- NULL; from})
		   )
	}
    }
    else {
	substitute({ value <- new(CLASS)
		     for(what in TOSLOTS)
			 slot(value, what) <- slot(from, what)
		     value },
		   list(CLASS = toClass, TOSLOTS = toSlots))
    }
}

.simpleReplaceExpr <- function(toDef) {
    toSlots <- names(toDef@slots)
    substitute({
        for(what in TOSLOTS)
            slot(from, what) <- slot(value, what)
        from
    }, list(TOSLOTS = toSlots))
}

## the boot version of newClassRepresentation (does no checking on slots to avoid
## requiring method selection on coerce).

newClassRepresentation <- function(...) {
    value <- new("classRepresentation")
    slots <- list(...)
    slotNames <- names(slots)
    for(i in seq_along(slotNames))
        slot(value, slotNames[[i]], FALSE) <- slots[[i]]
    value
}

## create a temporary definition of a class, but one that is distinguishable
## (by its class) from the real thing.  See comleteClassDefinition
.tempClassDef <- function(...) {
    value <- new("classRepresentation")
    slots <- list(...)
    slotNames <- names(slots)
    for(i in seq_along(slotNames))
        slot(value, slotNames[[i]], FALSE) <- slots[[i]]
    value
}

## the real version of newClassRepresentation, assigned in ..First.lib
.newClassRepresentation <- function(...)
    new("classRepresentation", ...)

.insertExpr <- function(expr, el) {
    if(!is(expr, "{"))
        expr <- substitute({EXPR}, list(EXPR = expr))
    expr[3L:(length(expr)+1)] <- expr[2L:length(expr)]
    expr[[2L]] <- el
    expr
}

## utility guaranteed to return only the first string of the class.
## Would not be needed if we dis-allowed S3 classes with multiple strings (or
## if the methods package version of class dropped the extra strings).
.class1 <- function(x) {
    cl <- class(x)
    if(length(cl) > 1L)
        cl[[1L]]
    else
        cl
}

substituteFunctionArgs <-
    function(def, newArgs, args = formalArgs(def), silent = FALSE,
             functionName = "a function")
{
    if(!identical(args, newArgs)) {
        if( !missing(functionName) ) # this style does not allow translation
            functionName <- paste("for", functionName)

        n <- length(args)
        if(n != length(newArgs))
            stop(sprintf("trying to change the argument list of %s with %d arguments to have arguments (%s)",
                         functionName, n, paste(newArgs, collapse = ", ")),
                 domain = NA)
        bdy <- body(def)
        ## check for other uses of newArgs
        checkFor <- newArgs[is.na(match(newArgs, args))]
        locals <- all.vars(bdy)
        if(length(checkFor) && any(!is.na(match(checkFor, locals))))
            stop(sprintf("get rid of variables in definition %s (%s); they conflict with the needed change to argument names (%s)",
                         functionName,
                         paste(checkFor[!is.na(match(checkFor, locals))], collapse = ", "),
                         paste(newArgs, collapse = ", ")), domain = NA)
        ll <- vector("list", 2L*n)
        for(i in seq_len(n)) {
            ll[[i]] <- as.name(args[[i]])
            ll[[n+i]] <- as.name(newArgs[[i]])
        }
        names(ll) <- c(args, newArgs)
        body(def, envir = environment(def)) <- substituteDirect(bdy, ll)
        if(!silent) {
            msg <-
                sprintf("NOTE: arguments in definition %s changed from (%s) to (%s)",
                        functionName,
                        paste(args, collapse = ", "),
                        paste(newArgs, collapse = ", "))
            message(msg, domain = NA)
        }
    }
    def
}

.makeValidityMethod <- function(Class, validity) {
    if(!is.null(validity)) {
        if(!is(validity, "function"))
            stop(gettextf("a validity method must be a function of one argument, got an object of class %s",
                          dQuote(class(validity))),
                 domain = NA)
        validity <- substituteFunctionArgs(validity, "object", functionName = sprintf("validity method for class '%s'", Class))
    }
    validity
}

# the bootstrap version of setting slots in completeClassDefinition
.mergeClassDefSlots <- function(ClassDef, ...) {
    slots <- list(...); slotNames <- names(slots)
    for(i in seq_along(slots))
        slot(ClassDef, slotNames[[i]], FALSE) <- slots[[i]]
    ClassDef
}

## the real version:  differs only in checking the slot values
..mergeClassDefSlots <- function(ClassDef, ...) {
    slots <- list(...); slotNames <- names(slots)
    for(i in seq_along(slots))
        slot(ClassDef, slotNames[[i]]) <- slots[[i]]
    ClassDef
}

### fix the annoying habit of R giving function definitions the local environment by default
.gblEnv <- function(f) {
    environment(f) <- .GlobalEnv
    f
}

## a utility for makePrototypeFromClassDef that causes inf. recursion if used too early
..isPrototype <- function(p)is(p, "classPrototypeDef")
## the simple version
.isPrototype <- function(p) .identC(class(p), "classPrototypeDef")

.className <- function(cl) if(is(cl, "classRepresentation")) cl@className else as(cl, "character")

## bootstrap version:  all classes and methods must be in the version of the methods
## package being built in the toplevel environment: MUST avoid require("methods") !
.requirePackage <- function(package, mustFind = TRUE)
    topenv(parent.frame())

.PackageEnvironments <- new.env(hash=TRUE) # caching for required packages

## real version of .requirePackage
..requirePackage <- function(package, mustFind = TRUE) {
    value <- package
    if(nzchar(package)) {
        if(package %in% loadedNamespaces())
            value <- getNamespace(package)
        else {
            if(identical(package, ".GlobalEnv"))
                return(.GlobalEnv)
            if(identical(package, "methods"))
                return(topenv(parent.frame())) # booting methods
            if(exists(package, envir = .PackageEnvironments, inherits = FALSE))
                return(get(package, envir = .PackageEnvironments)) #cached, but only if no namespace
        }
    }
    if(is.environment(value))
        return(value)
    topEnv <- options()$topLevelEnvironment
    if(is.null(topEnv))
        topEnv <- .GlobalEnv
    if(exists(".packageName", topEnv, inherits=TRUE) &&
       .identC(package, get(".packageName", topEnv)))
        return(topEnv) # kludge for source'ing package code
    if(nzchar(package) && require(package, character.only = TRUE)) {}
    else {
        if(mustFind)
          stop(gettextf("unable to find required package %s",
                        sQuote(package)),
               domain = NA)
        else
          return(NULL)
    }
    value <- .asEnvironmentPackage(package)
    assign(package, value, envir = .PackageEnvironments)
    value
}

.classDefEnv <- function(classDef) {
    .requirePackage(classDef@package)
}


.asEnvironmentPackage <- function(package) {
    if(identical(package, ".GlobalEnv"))
        .GlobalEnv
    else {
        ##FIXME:  the paste should not be needed
        pkg <- paste("package", package, sep=":")
        ## need to allow for versioned installs: prefer exact match.
        m <- charmatch(pkg, search())
        if(is.na(m)) # not attached, better be an available namespace
            getNamespace(package)
        else
          as.environment(search()[m])
    }
}

## bootstrap version, mustn't fail
.classEnv <- function(Class, default = .requirePackage("methods"), mustFind = TRUE) {
         package <- packageSlot(Class)
        if(is.null(package)) {
            ## unconditionally use the methods package
            default
        }
        else
            .requirePackage(package)
     }


## to be .classEnv()  --- currently used in 'Matrix'  (via wrapper)
..classEnv <- function(Class, default = .requirePackage("methods"), mustFind = TRUE) {
    package <- { if(is.character(Class)) packageSlot(Class) else
		 ## must then be a class definition
		 Class@package }
    if(is.null(package)) {
	## use the default, but check that the class is there, and if not
	## try a couple of other heuristics
	value <- default
	def <- getClassDef(Class, value, NULL)
	if(is.null(def)) {
	    value <- .GlobalEnv
	    def <- getClassDef(Class, value, NULL)
	    if(is.null(def)) {
		value <- .requirePackage("methods")
		if(!identical(default, value)) # user supplied default
		    def <- getClassDef(Class, value, NULL)
	    }
	}
	if(is.null(def) && mustFind)
	    stop(gettextf("unable to find an environment containing class %s",
			  dQuote(Class)),
                 domain = NA)
	value
    }
    else
	.requirePackage(package)
}

## find a generic function reference, using the package slot if present
## FIXME:  this and .classEnv should be combined and implemented in C for speed
## They differ in that  .classEnv uses the class metaname when it searches; i.e.,
## they use getClassDef and .getGeneric resp.  Also, .getEnv returns baseenv() rather
## than generating an error if no generic found (so getGeneric can return gen'c for prim'ves)

.genEnv <-  function(f, default = .requirePackage("methods"), package = "")
{
    if(!nzchar(package))
        package <- packageSlot(f)
    if(is.null(package)) {
        ## use the default, but check that the object is there, and if not
        ## try a couple of other heuristics
        value <- default
        def <- .getGeneric(f, value)
        if(is.null(def)) {
            value <- .GlobalEnv
            def <- .getGeneric(f, value)
            if(is.null(def)) {
                value <- .requirePackage("methods")
                if(!identical(default, value)) # user supplied default
                    def <- .getGeneric(f, value)
            }
        }
        if(is.null(def))
            baseenv()
        else
            value
    }
    else
        .requirePackage(package)
}

## cache and retrieve class definitions  If there is a conflict with
## packages a list of  classes will be cached
## See .cacheGeneric, etc. for analogous computations for generics
.classTable <- new.env(TRUE, baseenv())
assign("#HAS_DUPLICATE_CLASS_NAMES", FALSE, envir = .classTable)
.duplicateClassesExist <- function(on) {
    value <- get("#HAS_DUPLICATE_CLASS_NAMES", envir = .classTable)
    if(nargs())
        assign("#HAS_DUPLICATE_CLASS_NAMES", on, envir = .classTable)
    value
}

.cacheClass <- function(name, def, doSubclasses = FALSE, env) {
    if(!identical(doSubclasses, FALSE))
      .recacheSubclasses(def@className, def, doSubclasses, env)
    if(exists(name, envir = .classTable, inherits = FALSE)) {
        newpkg <- def@package
        prev <- get(name, envir = .classTable)
        if(is(prev, "classRepresentation")) {
            if(identical(prev, def))
               return()
            pkg <- prev@package # start a per-package list
            if(identical(pkg, newpkg)) { # redefinition
                ## cache for S3, to override possible previous cache
                base:::.cache_class(name, .extendsForS3(def))
##                base:::.cache_class(name, extends(def))
                return(assign(name, def, envir = .classTable))
            }
            else if(.simpleDuplicateClass(def, prev))
                return()
            prev <- list(prev)
            names(prev) <- pkg
        }
        i <- match(newpkg, names(prev))
        if(is.na(i))
           prev[[newpkg]] <- def
        else if(identical(def, prev[[i]]))
          return()
        else
            prev[[i]] <- def
        def <- prev
        .duplicateClassesExist(TRUE)
    }
    assign(name, def, envir = .classTable)
}

## test for identical def, prev class definitions
## An exhaustive test would be very complicated, having to test
## superclasses in detail, prototypes for the slots, etc.
.simpleDuplicateClass <- function(def, prev) {
    supers <- names(def@contains)
    prevSupers <- names(prev@contains)
    if(length(supers) != length(prevSupers) ||
       any(is.na(match(supers, prevSupers))))
        return(FALSE)
    warnLevel <- getOption("warn")
    S3 <- "oldClass" %in% supers
    if(S3) {
        ## it is possible one  of these is inconsistent, but unlikely
        ## and we will get here often from multiple setOldClass(...)'s
        if(warnLevel)
            message(gettextf("Note: the specification for S3 class %s in package %s seems equivalent to one from package %s: not turning on duplicate class definitions for this class.",
                             dQuote(def@className),
                             sQuote(def@package),
                             sQuote(prev@package)),
                    domain = NA)
        return(TRUE)
    }
    ## if there are already duplicate classes, we check duplicates
    ## for the superclasses
    dupsExist <- .duplicateClassesExist()
    if(dupsExist) {
        dups <- match(supers, multipleClasses(), 0) > 0
        if(any(dups)) {
            if(warnLevel)
                message(gettextf("Note: some superclasses of class %s in package %s have duplicate definitions.  This definition is not being treated as equivalent to that from package %s",
                                 dQuote(def@className),
                                 sQuote(def@package),
                                 sQuote(prev@package)),
                    domain = NA)
            return(FALSE)
        }
    }
    ## now check the slots
    slots <- names(def@slots)
    prevSlots <- names(prev@slots)
    if(length(slots) != length(prevSlots) ||
       any(is.na(match(slots, prevSlots))))
        return(FALSE)
    for(what in slots) {
        slotClasses <- def@slots
        prevClasses <- prev@slots
        clWhat <- slotClasses[[what]]
        prevWhat <- prevClasses[[what]]
        if(!identical(as.character(clWhat), as.character(prevWhat)) ||
           (dupsExist && !identical(as.character(packageSlot(clWhat)),
              as.character(packageSlot(prevWhat)))))
            return(FALSE)
    }
    if(warnLevel)
        message(gettextf("Note: the specification for class %s in package %s seems equivalent to one from package %s: not turning on duplicate class definitions for this class.",
                         dQuote(def@className),
                         sQuote(def@package),
                         sQuote(prev@package)),
                    domain = NA)
    TRUE
}

.uncacheClass <- function(name, def) {
    if(exists(name, envir = .classTable, inherits = FALSE)) {
        if(is(def, "classRepresentation")) # paranoia: should only be called this way
            newpkg <- def@package
        else
            newpkg <- ""
        prev <- get(name, envir = .classTable)
        if(is(prev, "classRepresentation") &&
           identical(prev@package, newpkg) )
            return(remove(list = name, envir = .classTable))
         i <- match(newpkg, names(prev))
        if(!is.na(i))
           prev[[i]] <- NULL
        else # we might warn about unchaching more than once
          return()
        if(length(prev) == 0L)
          return(remove(list = name, envir = .classTable))
        else if(length(prev) == 1L)
          prev <- prev[[1L]]
        assign(name, prev, envir  = .classTable)
    }
}

## the workhorse of class access
## The underlying C code will return name if it is not a character vector
## in the assumption this is a classRepresentation or subclass of that.
## In principle, this could replace the checks on class(name) in getClassDef
## and new(), which don't work for subclasses of classRepresentation anyway.
.getClassFromCache <- function(name, where) {
	value <- .Call(C_R_getClassFromCache, name, .classTable)
	if(is.list(value)) { ## multiple classes with this name
	    pkg <- packageSlot(name)
	    if(is.null(pkg))
		pkg <- if(is.character(where)) where else getPackageName(where, FALSE) # may be ""
	    pkgs <- names(value)
	    i <- match(pkg, pkgs, 0L)
	    if(i == 0L) ## try 'methods':
		i <- match("methods", pkgs, 0L)
	    if(i > 0L) value[[i]]
            else NULL
	}
	else #either a class definition or NULL
	    value
}

### insert superclass information into all the subclasses of this
### class.  Used to incorporate inheritance information from
### ClassUnions
.recacheSubclasses <- function(class, def, doSubclasses, env) {
    subs <- def@subclasses
    subNames <- names(subs)
    for(i in seq_along(subs)) {
        what <- subNames[[i]]
        subDef <- getClassDef(what, env)
        if(is.null(subDef))
            warning(gettextf("undefined subclass %s of class %s; definition not updated",
                             .dQ(what), .dQ(def@className)))
        else if(is.na(match(what, names(subDef@contains)))) {
            ## insert the new superclass to maintain order by distance
            cntns <- subDef@contains
            cntns[[class]] <- subs[[i]]
            cntns <- cntns[sort.list(sapply(cntns, function(x)x@distance))]
            subDef@contains <- cntns
            .cacheClass(what, subDef, FALSE, env)
        }
    }
    NULL
}

## alternative to .recacheSubclasses, only needed for non-unions
## Inferior in that nonlocal subclasses will not be updated, hence the
## warning when the subclass is not in where
.checkSubclasses <- function(class, def, class2, def2, where, where2) {
    where <- as.environment(where)
    where2 <- as.environment(where2)
   subs <- def@subclasses
    subNames <- names(subs)
    extDefs <- def2@subclasses
    for(i in seq_along(subs)) {
        what <- subNames[[i]]
        if(.identC(what, class2))
          next # catch recursive relations
        cname <- classMetaName(what)
        if(exists(cname, envir = where, inherits = FALSE)) {
            subDef <- get(cname, envir = where)
            cwhere <- where
        }
        else if(exists(cname, envir = where2, inherits = FALSE)) {
            subDef <- get(cname, envir = where2)
            cwhere <- where2
        }
        else {
          warning(gettextf("subclass %s of class %s is not local and cannot be updated for new inheritance information; consider setClassUnion()",
                           .dQ(what), .dQ(class)),
                  call. = FALSE, domain = NA)
          next
        }
        extension <- extDefs[[what]]
        if(is.null(extension)) # not possible if the setIs behaved?
          warning(gettextf("no definition of inheritance from %s to %s, though the relation was implied by the setIs() from %s",
                           .dQ(what), .dQ(def2@className), .dQ(class)),
                  call. = FALSE, domain = NA)
        else if(is.na(match(class2, names(subDef@contains)))) {
            subDef@contains[[class2]] <- extension
            assignClassDef(what, subDef, cwhere, TRUE)
        }
    }
    NULL
}

.removeSuperclassBackRefs <- function(Class, classDef, classWhere)
{
    if(length(classDef@contains)) {
        superclasses <- names(classDef@contains)
        for(what in superclasses) {
            superWhere <- findClass(what, classWhere)
            if(length(superWhere)) {
                superWhere <- superWhere[[1L]]
                .removeSubClass(what, Class, superWhere)
            } else if(! what %in% c(.BasicClasses, "oldClass"))
                warning(gettextf("could not find superclass %s to clean up when removing subclass references to class %s",
                                 .dQ(what), .dQ(Class)))
        }
    }
    NULL
}


## remove subclass from the known subclasses of class
## both in the package environment and in the cache
.removeSubClass <- function(class, subclass, where) {
    mname <- classMetaName(class)
    where <- as.environment(where)
    if(exists(mname, envir = where, inherits = FALSE)) {
        cdef <- get(mname, envir = where)
        newdef <- .deleteSubClass(cdef, subclass)
        if(!is.null(newdef))
          assignClassDef(class, newdef,  where, TRUE)
        else { # check the cache
            cdef <- .getClassFromCache(cdef@className, where)
            if(is.null(cdef)) {}
            else {
                newdef <- .deleteSubClass(cdef, subclass)
                if(!is.null(newdef))
                  .cacheClass(class, newdef, FALSE, where)
            }
        }
        sig <- signature(from=subclass, to=class)
        if(existsMethod("coerce", sig))
          .removeCachedMethod("coerce", sig)
        if(existsMethod("coerce<-", sig))
          .removeCachedMethod("coerce<-", sig)
        if(is(cdef, "classRepresentation"))
            .uncacheClass(class, cdef)
    }
    else
      warning(gettextf("no class %s found as expected in removing subclass %s",
                       .dQ(class), .dQ(subclass)))
}

.deleteSubClass <- function(cdef, subclass) {
        subclasses <- cdef@subclasses
        ii <- match(subclass, names(subclasses), 0)
        ## the subclass may not be there, e.g., if that class has been
        ## unloaded.
        if(ii > 0) {
            cdef@subclasses <- subclasses[-ii]
            cdef
        }
        else
          NULL
    }

## remove superclass from  definition of class in the cache & in environments
## on search list
.removeSuperClass <- function(class, superclass) {
    cdef <- .getClassFromCache(class, where)
    if(is.null(cdef)) {}
    else {
        newdef <- .deleteSuperClass(cdef, superclass)
        if(!is.null(newdef))
          .cacheClass(class, newdef, FALSE, where)
    }
    sig <- signature(from=class, to=superclass)
    if(existsMethod("coerce", sig))
      .removeCachedMethod("coerce", sig)
    if(existsMethod("coerce<-", sig))
      .removeCachedMethod("coerce<-", sig)
    evv <- findClass(class, .GlobalEnv) # what about hidden classes?  how to find them?
    mname <- classMetaName(class)
    for(where in evv) {
        if(exists(mname, envir = where, inherits = FALSE)) {
            cdef <- get(mname, envir = where)
            newdef <- .deleteSuperClass(cdef, superclass)
            if(!is.null(newdef)) {
              assignClassDef(class, newdef,  where, TRUE)
              ## message("deleted ",superclass, " from ",class, "in environment")
          }
        }
    }
    NULL
}

.deleteSuperClass <- function(cdef, superclass) {
        superclasses <- cdef@contains
        ii <- match(superclass, names(superclasses), 0)
        if(ii > 0) {
            cdef@contains <- superclasses[-ii]
            for(subclass in names(cdef@subclasses))
              .removeSuperClass(subclass, superclass)
            cdef
        }
        else
          NULL
    }

classesToAM <- function(classes, includeSubclasses = FALSE,
                        abbreviate = 2) {
  .mergeMatrices <- function(m1, m2) {
    if(nrow(m1) == 0)
      return(m2)
    dn1 <- dimnames(m1)
    dn2 <- dimnames(m2)
    rows <- unique(c(dn1[[1]], dn2[[1]]))
    columns <- unique(c(dn1[[2]], dn2[[2]]))
    value <- matrix(0, length(rows), length(columns), dimnames = list(rows, columns))
    value[dn1[[1]], dn1[[2]] ] <- m1
    value[dn2[[1]], dn2[[2]] ] <- m2
    value
  }
  if(length(includeSubclasses) == 1)
    includeSubclasses <- rep.int(includeSubclasses, length(classes))
  if(!is(includeSubclasses, "logical") || length(includeSubclasses) != length(classes))
    stop("argument 'includeSubclasses' must be a logical, either one value or a vector of the same length as argument 'classes'")
  value <- matrix(0,0,0)
  for(i in seq_along(classes)) {
    class <- classes[[i]] # to allow for package attribute
    classDef <- getClass(class) # throws an error if undefined.  Make a warning?
    value <- .mergeMatrices(value, .oneClassToAM(classDef, includeSubclasses[[i]]))
  }
  abbr <- match(as.integer(abbreviate), 0:3)-1
  if(length(abbr) != 1 || is.na(abbr))
    stop("argument 'abbreviate' must be 0, 1, 2, or 3")
  if(abbr %% 2)
    dimnames(value)[[1]] <- base::abbreviate(dimnames(value)[[1]])
  if(abbr %/% 2)
    dimnames(value)[[2]] <- base::abbreviate(dimnames(value)[[2]])
  value
}

.oneClassToAM <- function(classDef, includeSubclasses = FALSE, short = FALSE) {
    findEdges <- function(extensions) {
        superclasses <- names(extensions)
        edges <- numeric()
        for(what in superclasses) {
            whatDef <- getClassDef(what)
            ifrom <- match(what, nodes)
            if(is.null(whatDef) || is.na(ifrom))
              next
            exts <- whatDef@contains
            whatedges <- names(exts)
            ito <- match(whatedges, nodes, 0)
            for(i in seq_along(exts))
              if(ito[[i]] >0 && exts[[i]]@distance == 1)
                edges <- c(edges, ifrom, ito[[i]])
        }
        edges
    }
    nodes <- c(classDef@className, names(classDef@contains))
    if(includeSubclasses)
      nodes <- c(nodes, names(classDef@subclasses))
    nodes <- unique(nodes)
    labels <-
        if(isTRUE(short)) abbreviate(nodes)
        else if(is.character(short)) {
            if(length(short) != length(nodes))
                stop(gettextf("needed the supplied labels vector of length %d, got %d",
                              length(nodes), length(short)), domain = NA)
            else short
        } else nodes
    size <- length(nodes)
    value <- matrix(0, size, size, dimnames = list(labels, labels))
    ifrom <- match(classDef@className, nodes) # well, 1, but just for consistency
    ## the following could use the current fact that direct superclasses come
    ## first, but the efficiency gain is minor, so we use the findEdges logic
    extensions <- classDef@contains
    superclasses <- names(extensions)
    ito <- match(superclasses, nodes)
    edges <- numeric()
    for(i in seq_along(extensions)) {
        exti <- extensions[[i]]
        if(exti@distance == 1)
            edges <- c(edges, ifrom, ito[[i]])
    }
    edges <- c(edges, findEdges(classDef@contains))
    if(includeSubclasses) {
        edges <- c(edges, findEdges(classDef@subclasses))
    }
    edges <- t(matrix(edges, nrow=2))
    value[edges] <- 1
    value
}

.choosePos <- function (thisClass, superclasses, subNames, affected)
  ## find if possible a set of superclass relations that gives a consistent
  ## ordering and eliminates any duplicates in the affected relations
  ## Note that the returned indices are against the index of superclasses
  ## If no successful selection is possible, return (one of) the best
  ## attempt, and the superclass(es) inconsistently embedded
{
    candidates <- list()
    allNames <- c(thisClass, superclasses)
    dups <- unique(superclasses[affected])
    whichCase <- names(subNames)
    for(what in dups) {
        where <- seq_along(allNames)[match( allNames, what,0)>0]
        ## make a list of all the subsets to remove duplicates
        whatRemove <- lapply(-seq_along(where), function(x,y) y[x], y=where)
        if(length(candidates) == 0)
          candidates <- whatRemove
        else # all the pairwise combinations with the previous
          candidates <- outer(candidates, whatRemove,
                              function(x,y)mapply(c,x,y, SIMPLIFY=FALSE))
    }
    ## check each way to make the list unique against each superclass extension
    problems <- function(x,y) any(diff(match(y, x))<0)
    possibles <- lapply(candidates, function(x, names)names[-x], names=allNames)
    ## the next could be vectorized, but here we choose instead to exit early.
    scores <- vector("list", length(possibles))
    for(i in seq_along(possibles)) {
        score <- sapply(subNames, problems, x=possibles[[i]])
        scores[[i]] <- whichCase[score]
        if(!any(score))
          return(-candidates[[i]]+1)
    }
    # the first min. scoring possibility and its score
    i <- which.min(sapply(scores, length))
    list(-candidates[[i]]+1, scores[[i]])
}

.checkGeneric <- function(what, where) {
  .checkFun <-  function(x) {
      maybe <- (if(exists(x, where)) {
        f <- get(x, where)
        is.function(f)
      }
      else
        FALSE)
      if(maybe)
        maybe <- is(f, "genericFunction") ||
              (length(grep("UseMethod", deparse(f))) > 0) ||
              is.primitive(f)
      maybe
    }
  sapply(what, .checkFun)
}


S3forS4Methods <- function(where, checkClasses = character()) {
  allClasses <- getClasses(where)
  if(length(checkClasses) > 0)
    allClasses <- allClasses[match(allClasses, checkClasses, 0) > 0]
  if(length(allClasses) == 0)
    return(allClasses)
  pattern <- paste0("([.]",allClasses, "$)", collapse="|")
  allObjects <- objects(where, all.names = TRUE)
  allObjects <- allObjects[-grep("^[.][_][_]", allObjects)] # remove meta data
  allObjects <- grep(pattern, allObjects, value = TRUE)
  if(length(allObjects) > 0) {
    badMethods <- allObjects
    funs <- sub(pattern, "", badMethods)
    uniqueFuns <- unique(funs)
    uniqueFuns <- uniqueFuns[nzchar(uniqueFuns)]
    possible <- .checkGeneric(uniqueFuns, where)
    if(!any(possible))
      return(character())
    uniqueFuns <- uniqueFuns[possible]
    badMethods <- badMethods[match(funs, uniqueFuns, 0) > 0]
    allObjects <- badMethods
    attr(allObjects, "functions") <- uniqueFuns
  }
  allObjects
}

## ## this function warns of S3 methods for S4 classes, but only once per package
## ## per session.
## .checkS3forS4 <- function(method) {
##   envir <- environment(method)
##   pkg <- getPackageName(envir)
##   if(!nzchar(pkg)) pkg <- getPackageName(parent.env(pkg)) #? if generic function
##   if(!nzchar(pkg)) pkg <- format(envir)
##   if(!exists(".WarnedS3forS4", .GlobalEnv, inherits = FALSE))
##     assign(".WarnedS3forS4", character(), envir = .GlobalEnv)
##   if(is.na(match(pkg, .WarnedS3forS4))) {
##       methods <-   S3forS4Methods(envir)
##       .WarnedS3forS4 <<- c(.WarnedS3forS4, pkg)
##       if(length(methods) > 0) {
##         warning("S3 methods written for S4 classes will fail inheritance!\nPackage ", pkg, " apparently has ",
##             length(methods), " such methods  for the functions ", paste(attr(methods, "functions"), collapse = ", "), "\n\n",
##         "Possible dangerous methods: ", paste(methods, collapse =", "),
##                 "\n\n(Warnings generated once per package per session)")
##       }
##   }
## }

## a warning when a class is defined that extends classes with S3 methods.
## .checkS3forClass <- function(className, where, what = className) {
##   badMethods <- S3forS4Methods(where, what)
##   if(length(badMethods) > 0) {
##     msg <- paste0("The apparent methods are ", paste('"',badMethods, '"', collapse = ", "))
##     warning("Some of the superclasses in the definition of class \"",
##             className, "\" have apparent S3 methods.\n\nThese will be hidden by the S3 class that this class contains. (See ?Methods)\n\n", msg)
##   }
## }

## a utility to detect mixin classes:  meant to be fast for use in
## initialize methods (cf the "matrix" method in BasicClasses.R)
isMixin <- function(classDef) {
    val <- 0
    cc <- classDef@contains
    ## relies on the superclasses in contains slot being ordered by distance
    for(cl in cc) {
        if(cl@distance > 1 || val > 1)
          break
        val <- val + 1
    }
    val > 1
}

.classDefIsLocked <- function(classDef) {
    what <- classMetaName(classDef@className)
    env <- .NamespaceOrEnvironment(classDef@package)
    is.environment(env) && exists(what, envir = env, inherits = FALSE) &&
       bindingIsLocked(what, env)
}

#  File src/library/methods/R/RMethodUtils.R
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

## The real version of makeGeneric, to be installed after there are some
## generic functions to boot the definition (in particular, coerce and coerce<-)

.makeGeneric <-
## Makes a generic function object corresponding to the given function name.
## and definition.
  function(f, fdef,
           fdefault = fdef,
           group = list(),
           valueClass = character(),
           package = getPackageName(environment(fdef)),
           signature = NULL,
           genericFunction = NULL,
           simpleInheritanceOnly = NULL)
{
    checkTrace <- function(fun, what, f) {
        if(is(fun, "traceable")) {
            warning(gettextf("the function being used as %s in making a generic function for %s is currently traced; the function used will have tracing removed",
                             what,
                             sQuote(f)),
                    domain = NA)
            .untracedFunction(fun)
        }
        else
            fun
    }
    if(missing(fdef)) {
        if(missing(fdefault))
            stop(gettextf("must supply either a generic function or a function as default for %s",
                          sQuote(f)),
                 domain = NA)
        else if(is.primitive(fdefault)) {
            return(genericForPrimitive(f))
        }
        fdef <- fdefault
        body(fdef) <- substitute(standardGeneric(NAME), list(NAME = f))
        environment(fdef) <- .methodsNamespace
    }
    ## give the function a new environment, to cache methods later
    ev <- new.env()
    parent.env(ev) <- environment(fdef)
    environment(fdef) <- ev
    packageSlot(f) <- package
    assign(".Generic", f, envir = ev)
    fdef <- checkTrace(fdef)
    if(length(valueClass))
        fdef <- .ValidateValueClass(fdef, f, valueClass)
    group <- .asGroupArgument(group)
    if(is.null(genericFunction))
        value <- new("standardGeneric")
    else if(is(genericFunction, "genericFunction"))
        value <- genericFunction
    else
        stop(gettextf("the %s argument must be NULL or a generic function object; got an object of class %s",
                      sQuote("genericFunction"),
                      dQuote(class(genericFunction))),
             domain = NA)
    value@.Data <- fdef
    value@generic <- f
    value@group <- group
    value@valueClass <- valueClass
    value@package <- package
    args <- formalArgs(fdef)
    if(is.null(signature))
        signature <- args
    else if(any(is.na(match(signature, args))))
        stop(sprintf(ngettext(sum(is.na(match(signature, args))),
                              "non-argument found in the signature: %s",
                              "non-arguments found in the signature: %s"),
                     paste(signature[is.na(match(signature, args))], collapse = ", ")),
             domain = NA)
    dots <- match("...", signature)
    if(!is.na(dots)) { # remove "..." unless it is the only element of the signature
        if(length(signature) > 1L)
            signature <- signature[-dots]
    }
    if(length(signature) == 0L)
        stop("no suitable arguments to dispatch methods in this function")
    attr(signature, "simpleOnly") <- simpleInheritanceOnly # usually NULL
    value@signature <- signature
    name <- signature[[1L]]
    if(is.null(fdefault))
        {} # pre 2.11.0: methods <- MethodsList(name)
    else {
        fdefault <- checkTrace(fdefault)
        if(!identical(formalArgs(fdefault), formalArgs(fdef)) &&
           !is.primitive(fdefault))
            stop(sprintf(ngettext(length(fdef),
                                  "the formal argument of the generic function for %s (%s) differs from that of the non-generic to be used as the default (%s)",
                                  "the formal arguments of the generic function for %s (%s) differ from those of the non-generic to be used as the default (%s)"),
                          paste(formalArgs(fdef), collapse = ", "),
                          paste(formalArgs(fdefault), collapse = ", ")),
                 domain = NA)
        fdefault <- asMethodDefinition(fdefault, fdef = value)
        if(is(fdefault, "MethodDefinition"))
            fdefault@generic <- value@generic
        ## pre 2.11.0 methods <- MethodsList(name, fdefault)
    }
    value@default <- fdefault # pre 2.11.0 methods
    assign(".Methods", fdefault, envir = ev) ## ? why
    .setupMethodsTables(value, TRUE)
    value@skeleton <- generic.skeleton(f, fdef, fdefault)
    value
}

## stripped down version of asS4 in base (asS4 can't be used until the methods
## namespace is available -- no longer true)
.asS4 <- function (object)
    asS4(object, TRUE, 0L)

.notS4 <- function (object)
    asS4(object, FALSE, 0L)


## the bootstrap version: "#----" brackets lines that replace parts of the real version
makeGeneric <-
      function(f, fdef,
           fdefault = getFunction(f, generic = FALSE, mustFind = FALSE),
           group = list(), valueClass = character(), package, signature = NULL,
           genericFunction = NULL, simpleInheritanceOnly = NULL)
{
    ## give the function a new environment, to cache methods later
    ev <- new.env()
    parent.env(ev) <- environment(fdef)
    environment(fdef) <- ev
    packageSlot(f) <- package
    assign(".Generic", f, envir = ev)
    if(length(valueClass))
        fdef <- .ValidateValueClass(fdef, f, valueClass)
    group <- .asGroupArgument(group)
###--------
    value <- .asS4(fdef)
    if(is.null(genericFunction))
        class(value) <- .classNameFromMethods("standardGeneric")
    else
        class(value) <- class(genericFunction)
    slot(value, "generic", FALSE) <- f
    slot(value, "group", FALSE) <- group
    slot(value, "valueClass", FALSE) <- valueClass
    slot(value, "package", FALSE) <- package
###--------
    args <- formalArgs(fdef)
    if(is.null(signature))
        signature <- args
    else if(any(is.na(match(signature, args))))
        stop(sprintf(ngettext(sum(is.na(match(signature, args))),
                              "non-argument found in the signature: %s",
                              "non-arguments found in the signature: %s"),
                     paste(signature[is.na(match(signature, args))], collapse = ", ")),
             domain = NA)
    attr(signature, "simpleOnly") <- simpleInheritanceOnly # usually NULL
    dots <- match("...", signature)
    if(!is.na(dots)) ## ... is not currently supported in method signatures
        signature <- signature[-dots]
    if(length(signature) == 0L)
        stop("no suitable arguments to dispatch methods in this function")
###--------
    slot(value, "signature", FALSE) <- signature
###--------
    name <- signature[[1L]]
    if(is.null(fdefault))
      {}
    else
        fdefault <- asMethodDefinition(fdefault, fdef = value)
        if(is(fdefault, "MethodDefinition"))
            fdefault@generic <- value@generic
        ## pre 2.11.0 methods <- MethodsList(name, fdefault)
###--------
    assign(".Methods", fdefault, envir = ev)
    slot(value, "default", FALSE) <- fdefault
    slot(value, "skeleton", FALSE) <- generic.skeleton(f, fdef, fdefault)
###--------
    value
}


makeStandardGeneric <-
  ## a utility function that makes a valid function calling standardGeneric for name f
  ## Works (more or less) even if the actual definition, fdef, is not a proper function,
  ## that is, it is a primitive or internal
  function(f, fdef)
{
    fgen <- fdef
    body(fgen) <- substitute(standardGeneric(FNAME), list(FNAME=f))
    ## detect R specials and builtins:  these don't provide an argument list
    if(typeof(fdef) != "closure") {
        ## Look in a list of pre-defined functions (and also of functions for which
        ## methods are prohibited)
        fgen <- genericForPrimitive(f)
        if(identical(fgen, FALSE))
            stop(gettextf("special function %s is not permitted to have methods",
                          sQuote(f)),
                 domain = NA)
        if(is.null(fgen)) {
            warning(gettextf("special function %s has no known argument list; will assume '(x, ...)'",
                             sQuote(f)),
                    domain = NA)
            ## unknown
            fgen <- function(x, ...) {}
        }
        else {
            message(gettextf("making a generic for special function %s",
                             sQuote(f)),
                    domain = NA)
            setPrimitiveMethods(f, fdef, "reset", fgen, NULL)
        }
        ## Note that the body of the function comes from the list.  In a few cases ("$"),
        ## this body is not just a call to standardGeneric
    }
    fgen
}

generic.skeleton <- function(name, fdef, fdefault)
{
    anames <- formalArgs(fdef)
    skeleton <- lapply(as.list(c(name, anames)), as.name)
    ## any arguments after "..." have to be named
    dots <- match("...", anames)
    if(!is.na(dots) && dots < length(anames)) {
        anames[1L:dots] <- ""
        names(skeleton) <- c("", anames)
    }
    if(is.null(fdefault)) {
        fdefault <- fdef
        body(fdefault) <- substitute(stop(MESSAGE, domain = NA), list(MESSAGE=
                                                   gettextf("invalid call in method dispatch to '%s' (no default method)", name)))
        environment(fdefault) <- baseenv()
    }
    skeleton[[1L]] <- fdefault
    as.call(skeleton)
}


defaultDumpName <-
  ## the default name to be used for dumping a method.
  function(generic, signature)
{
    if(missing(signature))
        paste(generic, "R", sep=".", collapse =".")
    else
        paste(generic, paste(signature, collapse ="."), "R", sep=".")
}


mergeMethods <-
    ## merge the methods in the second MethodsList object into the first,
    ## and return the merged result.
    function(m1, m2, genericLabel = character())
{
    if(length(genericLabel) && is(m2, "MethodsList"))
        m2 <- .GenericInPrimitiveMethods(m2, genericLabel)
    if(is.null(m1) || is(m1, "EmptyMethodsList"))
        return(m2)
    tmp <- listFromMlist(m2)
    sigs <- el(tmp, 1)
    methods <- el(tmp, 2)
    for(i in seq_along(sigs)) {
        sigi <- el(sigs, i)
        if(.noMlists() && !identical(unique(sigi), "ANY"))
          next
        args <- names(sigi)
        m1 <- insertMethod(m1, as.character(sigi), args, el(methods, i), FALSE)
    }
    m1
}

doPrimitiveMethod <-
  ## do a primitive call to builtin function 'name' the definition and call
  ## provided, and carried out in the environment 'ev'.
  ##
  ## A call to 'doPrimitiveMethod' is used when the actual method is a .Primitive.
  ##  (because primitives don't behave correctly as ordinary functions,
  ## not having either formal arguments nor a function body).
  function(name, def, call = sys.call(sys.parent()), ev = sys.frame(sys.parent(2)))
{
    cat("called doPrimitiveMethod\n\n")
    ## Store a local version of function 'name' back where the current version was
    ## called.  Restore the previous state there on exit, either removing or re-assigning.
    if(exists(name, envir=ev, inherits=FALSE)) {
        prev <- get(name, envir=ev)
        on.exit(assign(name, prev, envir = ev))
    }
    else
        on.exit(rm(list=name, envir=ev))
    assign(name, def, envir = ev)
    eval(call, ev)
}

.renderSignature <- function(f, signature)
{
    nm <- names(signature)
    nm[nzchar(nm)] <- paste0(nm[nzchar(nm)], "=")
    msig <- paste0(nm, '"', as.vector(signature), '"')
    msig <- paste(msig, collapse = ",")
    gettextf("in method for %s with signature %s: ", sQuote(f), sQuote(msig))
}

conformMethod <- function(signature, mnames, fnames,
			  f = "<unspecified>", fdef, method)
{
    sig0 <- signature
    fsig <- fdef@signature
    if(is.na(match("...", mnames)) && !is.na(match("...", fnames)))
        fnames <- fnames[-match("...", fnames)]
    imf <- match(fnames, mnames)
    omitted <- is.na(imf)
    if(is.unsorted(imf[!omitted]))
	stop(.renderSignature(f, signature),
             "formal arguments in method and generic do not appear in the same order",
             call. = FALSE)
    if(!any(omitted)) ## i.e. mnames contains all fnames
        return(signature)
    sigNames <- names(signature)
    omittedSig <- sigNames %in% fnames[omitted] #  names in signature & generic but not in method defn
### FIXME:  the test below is too broad, with all.names().  Would be nice to have a test
### for something like assigning to one of the omitted arguments.
    ##     missingFnames <- fnames[omitted]
    ##     foundNames <- missingFnames %in% all.names(body(method), unique = TRUE)
    ##     if(any(foundNames))
    ##         warning(gettextf("%s function arguments omitted from method arguments, (%s), were found in method definition",
    ##                       label, paste(missingFnames[foundNames], collapse = ", ")),
    ##              domain = NA)
    if(!any(omittedSig))
      return(signature)
    if(any(is.na(match(signature[omittedSig], c("ANY", "missing"))))) {
        bad <- omittedSig & is.na(match(signature[omittedSig], c("ANY", "missing")))
        bad2 <- paste0(fnames[bad], " = \"", signature[bad], "\"", collapse = ", ")
        stop(.renderSignature(f, sig0),
             gettextf("formal arguments (%s) omitted in the method definition cannot be in the signature", bad2),
             call. = TRUE, domain = NA)
    }
    else if(!all(signature[omittedSig] == "missing")) {
        omittedSig <- omittedSig && (signature[omittedSig] != "missing")
        .message("Note: ", .renderSignature(f, sig0),
                 gettextf("expanding the signature to include omitted arguments in definition: %s",
                          paste(sigNames[omittedSig], "= \"missing\"",collapse = ", ")))
        omittedSig <- seq_along(omittedSig)[omittedSig] # logical index will extend signature!
        signature[omittedSig] <- "missing"
    }
    ## remove trailing "ANY"'s
    n <- length(signature)
    while(.identC(signature[[n]], "ANY"))
        n <- n - 1L
    length(signature) <- n
    length(fsig) <- n
    names(signature) <- fsig
    signature
}

rematchDefinition <- function(definition, generic, mnames, fnames, signature)
{
    added <- any(is.na(match(mnames, fnames)))
    keepsDots <- !is.na(match("...", mnames))
    if(!added && keepsDots) {
        ## the formal args of the method must be identical to generic
        formals(definition, envir = environment(definition)) <- formals(generic)
        return(definition)
    }
    dotsPos <- match("...", fnames)
    if(added && is.na(dotsPos))
        stop(gettextf("methods can add arguments to the generic %s only if '...' is an argument to the generic", sQuote(generic@generic)),
             call. = TRUE)
    ## pass down all the names in common between method & generic,
    ## plus "..."  even if the method doesn't have it.  But NOT any
    ## arguments having class "missing" implicitly (see conformMethod),
    ## i.e., are not among 'mnames':
    useNames <- !is.na(imf <- match(fnames, mnames)) | fnames == "..."
    newCall <- lapply(c(".local", fnames[useNames]), as.name)

    ## Should not be needed, if conformMethod() has already been called:
    if(is.unsorted(imf[!is.na(imf)]))
	stop(.renderSignature(generic@generic, signature),
             "formal arguments in method and generic do not appear in the same order",
             call. = FALSE)

    ## leave newCall as a list while checking the trailing args
    if(keepsDots && dotsPos < length(fnames)) {
	## Trailing arguments are required to match.  This is a little
	## stronger than necessary, but this is a dicey case, because
	## the argument-matching may not be consistent otherwise (in
	## the generic, such arguments have to be supplied by name).
	## The important special case is replacement methods, where
	## value is the last argument.

	ntrail <- length(fnames) - dotsPos
	trailingArgs <- fnames[seq.int(to = length(fnames), length.out = ntrail)]
	if(!identical(	mnames[seq.int(to = length(mnames), length.out = ntrail)],
		      trailingArgs))
	    stop(gettextf("arguments (%s) after '...' in the generic must appear in the method, in the same place at the end of the argument list",
			  paste(trailingArgs, collapse=", ")),
                 call. = TRUE, domain = NA)
	newCallNames <- character(length(newCall))
	newCallNames[seq.int(to = length(newCallNames), length.out = ntrail)] <-
	    trailingArgs
	names(newCall) <- newCallNames
    }
    newCall <- as.call(newCall)
    newBody <- substitute({.local <- DEF; NEWCALL},
			  list(DEF = definition, NEWCALL = newCall))
    generic <- .copyMethodDefaults(generic, definition)
    body(generic, envir = environment(definition)) <- newBody
    generic
}

unRematchDefinition <- function(definition)
{
    ## undo the effects of rematchDefiniition, if it was used.
    ## Has the obvious disadvantage of depending on the implementation.
    ## If we considered the rematching part of the API, a cleaner solution
    ## would be to include the "as given to setMethod" definition as a slot
    bdy <- body(definition)
    if(.identC(class(bdy),"{") && length(bdy) > 1L) {
        bdy <- bdy[[2L]]
        if(.identC(class(bdy), "<-") &&
           identical(bdy[[2L]], as.name(".local")))
            definition <- bdy[[3L]]
    }
    definition
}

getGeneric <-
  ## return the definition of the function named f as a generic.
  ##
  ## If there is no definition, throws an error or returns
  ## NULL according to the value of mustFind.
  function(f, mustFind = FALSE, where, package = "")
{
    if(is.function(f)) {
        if(is(f, "genericFunction"))
            return(f)
        else if(is.primitive(f))
            return(genericForPrimitive(.primname(f)))
        else
            stop("argument 'f' must be a string, generic function, or primitive: got an ordinary function")
    }
    value <- if(missing(where))
	.getGeneric(f, , package) else
    .getGeneric(f, where, package)
    if(is.null(value) && exists(f, envir = baseenv(), inherits = FALSE)) {
        ## check for primitives
        baseDef <- get(f, envir = baseenv())
        if(is.primitive(baseDef)) {
            value <- genericForPrimitive(f)
            if(!is.function(value) && mustFind)
                stop(gettextf("methods cannot be defined for the primitive function %s",
                              sQuote(f)), domain = NA)
            if(is(value, "genericFunction"))
                value <- .cacheGeneric(f, value)
        }
    }
    if(is.function(value))
        value
    else {
        if(nzchar(package) && is.na(match(package, c("methods", "base")))) {
            ## try to load package, or attach it if necessary
            ev <- tryCatch(loadNamespace(package), error = function(e)e)
            if(is(ev, "error") &&
               require(package, character.only =TRUE))
                ev <- as.environment(paste("package",package,sep=":"))
            if(is.environment(ev))
                value <- .getGeneric(f, ev, package)
        }
        if(is.function(value))
            value
        else if(mustFind)
            ## the C code will have thrown an error if f is not a single string
            stop(gettextf("no generic function found for %s",
                          sQuote(f)),
                 domain = NA)
        else
            NULL
    }
}

## low-level version
.getGeneric <- function(f, where = .GlobalEnv, # default only for C search
                        package = "")
{
    ## do not search the cache if getGeneric() was called with explicit where=
    if(missing(where))
        value <- .getGenericFromCache(f, where,  package)
    else
        value <- NULL
    if(is.null(value)) {
        if(is.character(f) && f %in% "as.double") f <- "as.numeric"
        if(is.character(f) && !nzchar(f)) {
            message("Empty function name in .getGeneric")
            dput(sys.calls())
        }
        value <- .Call(C_R_getGeneric, f, FALSE, as.environment(where), package)
        ## cache public generics (usually these will have been cached already
        ## and we get to this code for non-exported generics)
        if(!is.null(value) && exists(f, .GlobalEnv) &&
           identical(get(f, .GlobalEnv), value))
            .cacheGeneric(f, value)
    }
    ##     if(is.null(value) && nzchar(package) && !identical(package, "base")) {
    ##         env <- .requirePackage(package, FALSE)
    ##         if(is.environment(env))
    ##           value <- .Call("R_getGeneric", f, FALSE, env, package,
    ##                      PACKAGE = "methods")
    ##     }
    value
}

## cache and retrieve generic functions.  If the same generic name
## appears for multiple packages, a named list of the generics is cached.
.genericTable <- new.env(TRUE, baseenv())

.implicitTable <- new.env(TRUE, baseenv())

.cacheGeneric <- function(name, def)
  .cacheGenericTable(name, def, .genericTable)

.cacheImplicitGeneric <- function(name, def)
   .cacheGenericTable(name, def, .implicitTable)

.cacheGenericTable <- function(name, def, table)
{
    fdef <- def
    if(exists(name, envir = table, inherits = FALSE)) {
        newpkg <- def@package
        prev <- get(name, envir = table)
        if(is.function(prev)) {
            if(identical(prev, def))
                return(fdef)
            ## the following makes the cached version != package
            ##  fdef <- def <- .makeGenericForCache(def)
            pkg <- prev@package
            if(identical(pkg, newpkg)) { # redefinition
                assign(name, def, envir = table)
                return(fdef)
            }
            prev <- list(prev)          # start a per-package list
            names(prev) <- pkg
        }
        i <- match(newpkg, names(prev))
        if(is.na(i))
            prev[[newpkg]] <- def # or, .makeGenericForCache(def) as above
        else if(identical(def, prev[[i]]))
            return(fdef)
        else
            prev[[i]] <- def  # or, .makeGenericForCache(def) as above
        def <- prev
    }

    .getMethodsTable(fdef)              # force initialization
    assign(name, def, envir = table)
    fdef
}

.uncacheGeneric <- function(name, def)
  .uncacheGenericTable(name, def, .genericTable)

.uncacheImplicitGeneric <- function(name, def)
  .uncacheGenericTable(name, def, .implicitTable)

.uncacheGenericTable <- function(name, def, table)
{
    if(exists(name, envir = table, inherits = FALSE)) {
        newpkg <- def@package
        prev <- get(name, envir = table)
        if(is.function(prev))  # we might worry if  prev not identical
            return(remove(list = name, envir = table))
        i <- match(newpkg, names(prev))
        if(!is.na(i))
            prev[[i]] <- NULL
        else           # we might warn about unchaching more than once
            return()
        if(length(prev) == 0L)
            return(remove(list = name, envir = table))
        else if(length(prev) == 1L)
            prev <- prev[[1L]]
        assign(name, prev, envir  = table)
    }
}

.getGenericFromCache <- function(name, where,  pkg = "")
   .getGenericFromCacheTable(name,where, pkg, .genericTable)

.getImplicitGenericFromCache <- function(name, where,  pkg = "")
   .getGenericFromCacheTable(name,where, pkg, .implicitTable)

.getGenericFromCacheTable <- function(name, where, pkg = "", table)
{
    if(exists(name, envir = table, inherits = FALSE)) {
        value <- get(name, envir = table)
        if(is.list(value)) {        # multiple generics with this name
            ## force a check of package name, even if argument is ""
            if(!nzchar(pkg)) {
                if(is.character(where))
                    pkg <- where
                else {
                    pkg <- attr(name, "package")
                    if(is.null(pkg))
                        pkg <- getPackageName(where, FALSE)
                    if(identical(pkg, ".GlobalEnv"))
                        pkg <- ""
                }
            }
            pkgs <- names(value)
            i <- match(pkg, pkgs, 0L)
            if(i > 0L)
                return(value[[i]])
            i <- match("methods", pkgs, 0L)
            if(i > 0L)
                return(value[[i]])
            i <- match("base", pkgs, 0L)
            if(i > 0L)
                return(value[[i]])
            else
                return(NULL)
        }
        else if(nzchar(pkg) && !identical(pkg, value@package))
            NULL
        else
            value
    }
    else
        NULL
}

.genericOrImplicit <- function(name, pkg, env)
{
    fdef <- .getGenericFromCache(name, env, pkg)
    if(is.null(fdef)) {
	penv <- tryCatch(getNamespace(pkg), error = function(e)e)
	if(!isNamespace(penv))	{      # no namespace--should be rare!
	    pname <- paste0("package:", pkg)
	    penv <- if(pname %in% search()) as.environment(pname) else env
	}
        fdef <- getFunction(name, TRUE, FALSE, penv)
        if(!is(fdef, "genericFunction")) {
            if(is.primitive(fdef))
                fdef <- genericForPrimitive(name, penv)
            else
                fdef <- implicitGeneric(name, penv)
        }
    }
    fdef
}


## copy the environments in the generic function so later merging into
## the cached generic will not modify the generic in the package.
## NOT CURRENTLY USED: see comments in .getGeneric()
.makeGenericForCache <- function(fdef)
{
    value <- fdef
    ev <- environment(fdef)
    environment(value) <- newEv <- new.env(TRUE, parent.env(ev))
    for(what in objects(ev, all.names=TRUE)) {
        obj <- get(what, envir = ev)
        if(is.environment(obj))
            obj <- .copyEnv(obj)
        assign(what, obj, envir = newEv)
    }
    value
}

.copyEnv <- function(env)
{
    value <- new.env(TRUE, parent.env(env))
    for(what in objects(env, all.names = TRUE))
        assign(what, get(what, envir = env), envir = value)
    value
}

getGroup <-
  ## return the groups to which this generic belongs.  If 'recursive=TRUE', also all the
  ## group(s) of these groups.
  function(fdef, recursive = FALSE, where = topenv(parent.frame()))
{
    if(is.character(fdef))
        fdef <- getGeneric(fdef, where = where)
    if(is(fdef, "genericFunction"))
        group <- fdef@group
    else
        group <- list()
    if(recursive && length(group)) {
        allGroups <- group
        for(gp in group) {
            fgp <- getGeneric(gp, where = where)
            if(is(fgp, "groupGenericFunction"))
                allGroups <- c(allGroups, Recall(fgp, TRUE, where))
        }
        if(length(allGroups) > 1L) {
            ids <- sapply(allGroups, function(x) {
                pkg <- packageSlot(x)
                if(is.null(pkg)) x
                else paste(x, pkg, sep=":")
            })
            allGroups <- allGroups[!duplicated(ids)]
        }
        allGroups
    }
    else
        group
}

getMethodsMetaData <- function(f, where = topenv(parent.frame()))
{
    fdef <- getGeneric(f, where = where)
    if(is.null(fdef))
        return(NULL)
    if(.noMlists()) {
        warning(sprintf("Methods list objects are not maintained in this version of R:  request for function %s may return incorrect information",
                        sQuote(fdef@generic)),
                domain = NA)
    }
    mname <- methodsPackageMetaName("M",fdef@generic, fdef@package)
    if (exists(mname, where = where, inherits = missing(where)))
        get(mname, where)
    else if(missing(where))
        .makeMlistFromTable(fdef)
    else
        .makeMlistFromTable(fdef, where)
}

assignMethodsMetaData <-
  ## assign value to be the methods metadata for generic f on database where.
  ## as of R 2.7.0 the mlist metadata is deprecated.
  ## If value is not a MethodsList,  only turns on primitives & groups
  function(f, value, fdef, where, deflt)
{
    where <- as.environment(where)
    if(is(value, "MethodsList")) {
        mname <- methodsPackageMetaName("M",fdef@generic, fdef@package)
        if(exists(mname, envir = where, inherits = FALSE) &&
           bindingIsLocked(mname, where))
          {}        # may be called from trace() with locked binding; ignore
        else
          assign(mname, value, where)
    }
    if(is.primitive(deflt))
        setPrimitiveMethods(f, deflt, "reset", fdef, NULL)
    if(is(fdef, "groupGenericFunction")) # reset or turn on members of group
        cacheGenericsMetaData(f, fdef, where = where, package = fdef@package)
}


## utility for getGenerics to return package(s)
.packageForGeneric <- function(object)
{
    if(is.list(object))                 # a list of objects
        lapply(object, .packageForGeneric)
    else if(is(object, "genericFunction"))
        object@package
    else ## ?? possibly a primitive
        "base"
}

getGenerics <- function(where, searchForm = FALSE)
{
    if(missing(where)) {
        ## all the packages cached ==? all packages with methods
        ## globally visible.  Assertion based on cacheMetaData + setMethod
        fnames <- as.list(objects(.genericTable, all.names=TRUE))
        packages <- vector("list", length(fnames))
        for(i in seq_along(fnames)) {
            obj <- get(fnames[[i]], envir = .genericTable)
            if(is.list(obj))
                fnames[[i]] <-  names(obj)
            packages[[i]] <- .packageForGeneric(obj)
        }
        new("ObjectsWithPackage", unlist(fnames), package=unlist(packages))
    }
    else {
        if(is.environment(where)) where <- list(where)
        these <- character()
        for(i in where)
            these <- c(these, objects(i, all.names=TRUE))
        metaNameUndo(unique(these), prefix = "T", searchForm = searchForm)
    }
}

## Find the pattern for methods lists or tables
## Currently driven by mlists, but eventually these will go away
## in favor of tables.

## always returns a compatible list, with an option of  prefix
.getGenerics <- function(where, trim = TRUE)
{
    if(missing(where)) where <- .envSearch(topenv(parent.frame()))
    else if(is.environment(where)) where <- list(where)
    these <- character()
    for(i in where) these <- c(these, objects(i, all.names=TRUE))
    these <- allThese <- unique(these)
    these <- these[substr(these, 1L, 6L) == ".__T__"]
    if(length(these) == 0L)
        return(character())
    funNames <- gsub(".__T__(.*):([^:]+)", "\\1", these)
    if(length(funNames) == 0L &&
       length(these[substr(these, 1L, 6L) == ".__M__"]))
        warning(sprintf("package %s seems to have out-of-date methods; need to reinstall from source",
                         sQuote(getPackageName(where[[1L]]))))
    packageNames <- gsub(".__T__(.*):([^:]+(.*))", "\\2", these)
    attr(funNames, "package") <- packageNames
    ## Would prefer following, but may be trouble bootstrapping methods
    ## funNames <- new("ObjectsWithPackage", funNames, package = packageNames)
    if(identical(trim, TRUE))
        funNames
    else {
        if(identical(trim, FALSE))
            these
        else
            gsub(".__T__", as.character(trim), these)
    }
}

cacheMetaData <-
    function(where, attach = TRUE, searchWhere = as.environment(where),
             doCheck = TRUE)
{
    ## a collection of actions performed on attach or detach
    ## to update class and method information.
    pkg <- getPackageName(where)
    classes <- getClasses(where)
    for(cl in classes) {
        cldef <- (if(attach) get(classMetaName(cl), where) # NOT getClassDef, it will use cache
                  else  getClassDef(cl, searchWhere))
        if(is(cldef, "classRepresentation")) {
            if(attach) {
                .cacheClass(cl, cldef, is(cldef, "ClassUnionRepresentation"), where)
            }
            else if(identical(cldef@package, pkg)) {
                .uncacheClass(cl, cldef)
                .removeSuperclassBackRefs(cl, cldef, searchWhere)
            }
        }
    }
    generics <- .getGenerics(where)
    packages <- attr(generics, "package")
    if(length(packages) <  length(generics))
        packages <- rep(packages, length.out = length(generics))
    if(attach && exists(".requireCachedGenerics", where, inherits = FALSE)) {
        others <- get(".requireCachedGenerics", where)
        generics <- c(generics, others)
        packages <- c(packages, attr(others, "package"))
    }
    ## check for duplicates
    dups <- duplicated(generics) & duplicated(packages)
    generics <- generics[!dups]
    for(i in seq_along(generics)) {
        f <- generics[[i]]
        fpkg <- packages[[i]]
        if(!identical(fpkg, pkg) && doCheck) {
            if(attach) {
                env <- as.environment(where)
                ## All instances of this generic in different attached packages must
                ## agree with the cached version of the generic for consistent
                ## method selection.
                if(exists(f, envir = env, inherits = FALSE)) {
                    def <- get(f, envir = env)
                    fdef <- .genericOrImplicit(f, fpkg, env)
                    if(is.function(def)) {
                        ## exclude a non-function of the same name as a primitive with methods (!)
                        if(identical(environment(def), environment(fdef)))
                            next        # the methods are identical
                        else if( is(fdef, "genericFunction")) {
                            .assignOverBinding(f, fdef,  env, FALSE)
                        }
                    }     # else, go ahead to update primitive methods
                }
                else          # either imported generic or a primitive
                    fdef <- getGeneric(f, FALSE, searchWhere, fpkg)
            }
            else
                fdef <- getGeneric(f, FALSE, searchWhere, fpkg)
        }
        else
            fdef <- getGeneric(f, FALSE, searchWhere, fpkg)
        if(!is(fdef, "genericFunction"))
            next ## silently ignores all generics not visible from searchWhere
        if(attach)
            .cacheGeneric(f, fdef)
        else
            .uncacheGeneric(f, fdef)
        methods <- .updateMethodsInTable(fdef, where, attach)
        cacheGenericsMetaData(f, fdef, attach, where, fdef@package, methods)
    }
    .doLoadActions(where, attach)
    invisible(NULL) ## as some people call this at the end of functions
}


cacheGenericsMetaData <- function(f, fdef, attach = TRUE,
                                  where = topenv(parent.frame()),
                                  package, methods)
{
    if(!is(fdef, "genericFunction")) {
	warning(gettextf("no methods found for %s; cacheGenericsMetaData() will have no effect",
			 sQuote(f)),
		domain = NA)
	return(FALSE)
    }
    if(missing(package))
        package <- fdef@package
### Assertion: methods argument unused except for primitives
### and then only for the old non-table case.
    deflt <- finalDefaultMethod(fdef@default) #only to detect primitives
    if(is.primitive(deflt)) {
	if(missing(methods)) ## "reset"
	    setPrimitiveMethods(f, deflt, "reset", fdef, NULL)
	else ## "set"
	    setPrimitiveMethods(f, deflt, "set", fdef, methods)
    }
    else if(isGroup(f, fdef = fdef)) {
	members <- fdef@groupMembers
	## do the computations for the members as well; important if the
	## members are primitive functions.
	for(ff in members) {
	    ffdef <- getGeneric(ff, where = where)
	    if(is(ffdef, "genericFunction"))
		Recall(ff, ffdef, attach, where,
                       methods = .getMethodsTable(ffdef))
	}
    }
    TRUE
}

setPrimitiveMethods <-
  function(f, fdef, code, generic, mlist = get(".Methods", envir = environment(generic)))
    .Call(C_R_M_setPrimitiveMethods, f, fdef, code, generic, mlist)

### utility to turn ALL primitive methods on or off (to avoid possible inf. recursion)
.allowPrimitiveMethods <-
  function(onOff) {
      if(onOff) code <- "SET"
      else code <- "CLEAR"
      .Call(C_R_M_setPrimitiveMethods, "", NULL, code, NULL, NULL)
  }


findUnique <- function(what, message, where = topenv(parent.frame()))
{
    where <- .findAll(what, where = where)
    if(length(where) > 1L) {
        if(missing(message))
            message <- sQuote(what)
        if(is.list(where))
            where <- unlist(where)
        if(is.numeric(where))
            where <- search()[where]
        warning(message,
                sprintf(" found on: %s; using the first one",
                        paste(sQuote(where), collapse = ", ")),
                domain = NA)
        where <- where[1L]
    }
    where
}

MethodAddCoerce <- function(method, argName, thisClass, methodClass)
{
    if(.identC(thisClass, methodClass))
        return(method)
    ext <- possibleExtends(thisClass, methodClass)
    ## if a non-simple coerce is required to get to the target class for
    ## dispatch, insert it in the method.
    if(is.logical(ext) || ext@simple)
        return(method)
    methodInsert <- function(method, addExpr) {
        if(is.function(method)) {
            newBody <- substitute({firstExpr; secondExpr},
                                  list(firstExpr = addExpr,
                                       secondExpr = body(method)))
            body(method, envir = environment(method)) <- newBody
        }
        else if(is(method, "MethodsList")) {
            methods <- method@allMethods
            for(i in seq_along(methods))
                methods[[i]] <- Recall(methods[[i]], addExpr)
            method@allMethods <- methods
        }
        method
    }
    addExpr <- substitute(XXX <- as(XXX, CLASS),
                          list(XXX = argName, CLASS = methodClass))
    methodInsert(method, addExpr)
}

missingArg <- function(symbol, envir = parent.frame(), eval = FALSE)
    .Call(C_R_missingArg, if(eval) symbol else substitute(symbol), envir)

balanceMethodsList <- function(mlist, args, check = TRUE)
{
    moreArgs <- args[-1L]
    if(length(moreArgs) == 0L)
        return(mlist)
    methods <- mlist@methods
    if(check && length(methods)) {
        ## check whether the current depth is enough (i.e.,
        ## whether a method with this no. of args or more was set before
        depth <- 0
        el <- methods[[1L]]
        while(is(el, "MethodsList")) {
            mm <- el@methods
            if(length(mm) == 0L)
                break
            depth <- depth+1L
            el <- mm[[1L]]
        }
        if(depth >= length(args))
            ## already balanced to this length: An assertion
            ## relying on balance having been used consistently,
            ## which in turn relies on setMethod being called to
            ## add methods.  If you roll your own, tough luck!
            return(mlist)
    }
    for(i in seq_along(methods)) {
        el <- methods[[i]]
        if(is(el, "MethodsList"))
            el <- Recall(el, moreArgs, FALSE)
        else {
            if(is(el, "MethodDefinition")) {
                el@target[moreArgs] <- "ANY"
                el@defined[moreArgs] <- "ANY"
            }
            for(what in rev(moreArgs))
                el <- new("MethodsList", argument = as.name(what),
                          methods = list(ANY = el))
        }
        methods[[i]] <- el
    }
    mlist@methods <- methods
    mlist
}


sigToEnv <- function(signature, generic)
{
    genericSig <- generic@signature
    package <- packageSlot(signature)
    if(is.null(package))
        parent <- environment(generic)
    else
        parent <- .requirePackage(package)
    value <- new.env(parent = parent)
    classes <- as.character(signature)
    args <- names(signature)
    for(i in seq_along(args))
        assign(args[[i]], classes[[i]], envir = value)
    ## missing args in signature have class "ANY"
    if(length(args) < length(genericSig))
        for(other in genericSig[is.na(match(genericSig, args))])
            assign(other, "ANY", envir = value)
    value
}

methodSignatureMatrix <- function(object, sigSlots = c("target", "defined"))
{
    if(length(sigSlots)) {
        allSlots <- lapply(sigSlots, slot, object = object)
        mm <- unlist(allSlots)
        mm <- matrix(mm, nrow = length(allSlots), byrow = TRUE)
        dimnames(mm) <- list(sigSlots, names(allSlots[[1L]]))
        mm
    }
    else matrix(character(), 0L, 0L)
}

.valueClassTest <- function(object, classes, fname)
{
    if(length(classes)) {
        for(Cl in classes)
            if(is(object, Cl)) return(object)
        stop(gettextf("invalid value from generic function %s, class %s, expected %s",
                      sQuote(fname),
                      dQuote(class(object)),
                      paste(dQuote(classes), collapse = " or ")),
             domain = NA)
    }
    ## empty test is allowed
    object
}


.getOrMakeMethodsList <- function(f, where, genericFun)
{
    allMethods <- getMethodsMetaData(f, where = where)
    if(is.null(allMethods)) {
        argName <- genericFun@signature[[1L]]
        allMethods <- new("MethodsList", argument = as.name(argName))
#         other <- getMethodsMetaData(f)
#         if(is.null(other))
#             ## this utility is called AFTER ensuring the existence of a generic for f
#             ## Therefore, the case below can only happen for a primitive for which
#             ## no methods currently are attached.  Make the primitive the default
#             deflt <- getFunction(f, generic = FALSE, mustFind = FALSE)
#         else
#             ## inherit the default method, if any
#             deflt <- finalDefaultMethod(other)
#         if(!is.null(deflt))
#             allMethods <- insertMethod(allMethods, "ANY", argName, deflt)
        }
    allMethods
}

.makeCallString <- function(def, name = substitute(def), args = formalArgs(def))
{
    if(is.character(def)) {
        if(missing(name))
            name <- def
        def <- getFunction(def)
    }
    if(is(def, "function"))
        paste0(name, "(", paste(args, collapse=", "), ")")
    else
        ""
}

.ValidateValueClass <- function(fdef, name, valueClass)
{
    ## include tests for value
    fbody <- body(fdef)
    body(fdef, envir = environment(fdef)) <-
        substitute(.valueClassTest(EXPR, VALUECLASS, FNAME),
                   list(EXPR = fbody, VALUECLASS = valueClass, FNAME = name))
    fdef
}

## interpret the group= argument to makeGeneric, allowing for char. argument
## and "" for compatibility.
## TO DO:  make it possible for this argument to be a group generic function
## (it may in fact work now).
.asGroupArgument <- function(group)
{
    if(is.character(group)) {
	if(identical(group, ""))
	    list()
	else
	    as.list(group) ## should we allow c(group, package) ?
    }
    else
	group
}

metaNameUndo <- function(strings, prefix, searchForm = FALSE)
{
    pattern <- methodsPackageMetaName(prefix, "")
    n <- nchar(pattern, "c")
    matched <- substr(strings, 1L, n) == pattern
    value <- substring(strings[matched], n+1L)
    pkg <- sub("^[^:]*", "", value)   # will be "" if no : in the name
    if(searchForm) {
        global <- grep(".GlobalEnv", value)
        if(length(global)) {
            pkg[-global] <- paste0("package", pkg[-global])
            pkg[global] <- substring(pkg[global], 2L)
        }
    }
    else
        pkg <- substring(pkg, 2L)
    value <- sub(":.*","", value)
    new("ObjectsWithPackage", value, package = pkg)
}

.recursiveCallTest <- function(x, fname)
{
    if(is(x, "call")) {
        if(identical(x[[1L]], quote(standardGeneric))) {
            if(!identical(x[[2L]], fname))
                warning(gettextf("the body of the generic function for %s calls 'standardGeneric' to dispatch on a different name (\"%s\")!",
                                 sQuote(fname),
                                 paste(as.character(x[[2L]]), collapse = "\n")),
                        domain = NA)
            TRUE
        }
        else {
            for(i in seq.int(from=2L, length.out = length(x)-1L)) {
                if(Recall(x[[i]], fname))
                    return(TRUE)
            }
            FALSE
        }
    }
    else if(is(x, "language")) {
        for(i in seq.int(from=2L, length.out = length(x)-1L)) {
            if(Recall(x[[i]], fname))
                return(TRUE)
        }
        FALSE
    }
    else
        FALSE
}

.NonstandardGenericTest <- function(body, fname, stdBody)
{
    if(identical(body, stdBody))
        FALSE
    else if(.recursiveCallTest(body, fname))
        TRUE
    else
        NA
}

.GenericInPrimitiveMethods <- function(mlist, f)
{
    methods <- mlist@methods
    for(i in seq_along(methods)) {
        mi <- methods[[i]]
        if(is(mi, "function")) {
            body(mi, envir = environment(mi)) <-
                substitute({.Generic <- FF; BODY},
                           list(FF = f,BODY = body(mi)))
        }
        else if(is(mi, "MethodsList"))
            mi <- Recall(mi, f)
        else
            stop(sprintf("internal error: Bad methods list object in fixing methods for primitive function %s",
                          sQuote(f)),
                 domain = NA)
        methods[[i]] <- mi
    }
    mlist@methods <- methods
    mlist
}

.signatureString <- function(fdef, signature)
{
    snames <- names(signature)
    if(is.null(snames)) {
        if(is(fdef, "genericFunction")) {
            snames <- fdef@signature
            signature <- matchSignature(signature, fdef)
            if(length(snames) > length(signature))
                length(snames) <- length(signature)
        }
        else                            # shouldn't happen,...
            return(paste(signature, collapse=", "))
    }
    else
        signature <- as.character(signature)
    paste(paste0(snames, "=\"", signature, "\""), collapse = ", ")
}

.ChangeFormals <- function(def, defForArgs, msg = "<unidentified context>")
{
    if(!is(def, "function"))
        stop(gettextf("trying to change the formal arguments in %s in an object of class %s; expected a function definition",
                      msg, dQuote(class(def))),
             domain = NA)
    if(!is(defForArgs, "function"))
        stop(gettextf("trying to change the formal arguments in %s, but getting the new formals from an object of class %s; expected a function definition",
                      msg, dQuote(class(def))),
             domain = NA)
    old <- formalArgs(def)
    new <- formalArgs(defForArgs)
    if(length(old) < length(new))
        stop(gettextf("trying to change the formal arguments in %s, but the number of existing arguments is less than the number of new arguments: (%s) vs (%s)",
                      msg, paste0("\"", old, "\"", collapse=", "),
                      paste0("\"", new, "\"", collapse=", ")),
             domain = NA)
    if(length(old) > length(new))
        warning(gettextf("trying to change the formal arguments in %s, but the number of existing arguments is greater than the number of new arguments (the extra arguments won't be used): (%s) vs (%s)",
                         msg, paste0("\"", old, "\"", collapse=", "),
                         paste0("\"", new, "\"", collapse=", ")),
                domain = NA)
    if(identical(old, new))           # including the case of 0 length
        return(def)
    dlist <- as.list(def)
    slist <- lapply(c(old, new), as.name)
    names(slist) <- c(new, old)
    vlist <- dlist
    for(i in seq_along(vlist))
        vlist[[i]] <- do.call("substitute", list(vlist[[i]], slist))
    dnames <- names(dlist)
    whereNames <- match(old, dnames)
    if(any(is.na(whereNames)))
	stop(gettextf("in changing formal arguments in %s, some of the old names are not in fact arguments: %s",
		      msg, paste0("\"", old[is.na(match(old, names(dlist)))], "\"", collapse=", ")),
	     domain = NA)
    dnames[whereNames] <- new
    names(vlist) <- dnames
    as.function(vlist, envir = environment(def))
}

## The search list, or a namespace's static search list, or an environment
.envSearch <- function(env = topenv(parent.frame()))
{
    if(identical(env, .GlobalEnv))
        seq_along(search())
    else if(isNamespace(env) && !isBaseNamespace(env)) {
        ## the static environments for this namespace, ending with the base namespace
        value <- list(env)
        repeat {
            if(identical(env, emptyenv()))
                stop("botched namespace: failed to find 'base' namespace in its parents", domain = NA)
            env <- parent.env(env)
            value <- c(value, list(env))
            if(isBaseNamespace(env))
                break
        }
        value
    }
    else
        list(env)
}

.genericName <- function(f)
{
    if(is(f, "genericFunction"))
        f@generic
    else
        as.character(f)
}

## the environment in which to start searching for methods, etc. related
## to this generic function.  Will normally be the namespace of the generic's
## home package, or else the global environment
.genericEnv <- function(fdef)
    parent.env(environment(fdef))

## the default environment in which to start searching for methods, etc. relative to this
## call to a methods package utility.  In the absence of other information, the current
## strategy is to look at the function _calling_ the methods package utility.
##TODO:  this utility can't really work right until the methods package itself has a
## namespace, so that calls from within the package can be detected.  The
## heuristic is that all callers are skipped as long as their enviornment is  identical
## to .methodsNamespace.  But that is currently initialized to .GlobalEnv.
##
## The logic will fail if a function in a package with a namespace calls a (non-methods)
## function in a package with no namespace, and that function then calls a methods package
## function.  The right answer then is .GlobalEnv, but we will instead get the package
## namespace.
.externalCallerEnv <- function(n = 2, nmax = sys.nframe() - n + 1)
{
    ## start n generations back; by default the caller of the caller to this function
    ## go back nmax at most (e.g., a function in the methods package that knows it's never
    ## called more than nmax levels in could supply this argument
    if(nmax < 1) stop("got a negative maximum number of frames to look at")
    ev <- topenv(parent.frame()) # .GlobalEnv or the environment in which methods is being built.
    for(back in seq.int(from = -n, length.out = nmax)) {
        fun <- sys.function(back)
        if(is(fun, "function")) {
            ## Note that "fun" may actually be a method definition, and still will be counted.
            ## This appears to be the correct semantics, in
            ## the sense that, if the call came from a method, it's the method's environment
            ## where one would expect to start the search (for a class definition, e.g.)
            ev <- environment(fun)
            if(!identical(ev, .methodsNamespace))
                break
        }
    }
    ev
}

## a list of environments, starting from ev, going back to the base package,
## or else terminated by finding a namespace
.parentEnvList <- function(ev)
{
    ev <- as.environment(ev)
    value <- list(ev)
    while(!isNamespace(ev)) {
        if(identical(ev, baseenv())) {
            value[[length(value)]] <- .BaseNamespaceEnv
            break
        } else if(identical(ev, emptyenv())) {
            break
        }
        ev <- parent.env(ev)
        value <- c(value, list(ev))
    }
    value
}

.genericAssign <- function(f, fdef, methods, where, deflt)
{
    ev <- environment(fdef)
    assign(".Methods", methods, ev)
}

## Mark the method as derived from a non-generic.
.derivedDefaultMethod <- function(fdef)
{
    if(is.function(fdef) && !is.primitive(fdef)) {
        value <- new("derivedDefaultMethod")
        value@.Data <- fdef
        value@target <- value@defined <- .newSignature(.anyClassName, formalArgs(fdef))
        value
    }
    else
        fdef
}

.identC <- function(c1 = NULL, c2 = NULL)
{
    ## are the two objects identical class or genric function string names?
    .Call(C_R_identC, c1, c2)
}

## match default exprs in the method to those in the generic
## if the method does not itself specify a default, and the
## generic does
matchDefaults <- function(method, generic)
{
    changes <- FALSE
    margs <- formals(method)
    gargs <- formals(generic)
    for(arg in names(margs)) {
        ##!! weird use of missing() here is required by R's definition
        ## of a missing arg as a name object with empty ("") name
        ## This is dangerously kludgy code but seems the only way
        ## to avoid spurious errors ("xxx missing with no default")
        marg <- margs[[arg]]
        garg <- gargs[[arg]]
        if(missing(marg) && !missing(garg)) {
            changes <- TRUE
            margs[arg] <- gargs[arg] # NOT  [[]], which woud fail for NULL element
        }
    }
    if(changes)
        formals(method, envir = environment(method)) <- margs
    method
}

getGroupMembers <- function(group, recursive = FALSE, character = TRUE)
{
    .recMembers <- function(members, where) {
        all = vector("list", length(members))
        for(i in seq_along(members)) {
            what <- members[[i]]
            f <- getGeneric(what, FALSE, where)
            if(!is.null(f))
                all[[i]] <- what
            if(is(f, "groupGenericFunction")) {
                newMem <- f@groupMembers
                all <- c(all, Recall(newMem, where))
            }
        }
        all
    }
    f <- getGeneric(group)
    if(is.null(f)) {
        warning(gettextf("%s is not a generic function (or not visible here)",
                         sQuote(f)),
                domain = NA)
        return(character())
    }
    else if(!is(f, "groupGenericFunction"))
        character()
    else {
        members <- f@groupMembers
        if(recursive) {
            where <- f@package
            if(identical(where, "base")) {
                where <- "methods"      # no generics actually on base
                members <- .recMembers(members, .methodsNamespace)
            }
            else
                members <- .recMembers(members, .asEnvironmentPackage(where))
        }
        if(character)
            sapply(members, function(x){
                if(is(x, "character"))
                    x
                else if(is(x, "genericFunction"))
                    x@generic
                else
		    stop(gettextf("invalid element in the \"groupMembers\" slot (class %s)",
				  dQuote(class(x))),
                         domain = NA)
            })
        else
            members
    }
}

.primname <- function(object)
{
    ## the primitive name is 'as.double', but S4 methods are
    ## traditionally set on 'as.numeric'
    f <- .Call(C_R_get_primname, object)
    if(f == "as.double") "as.numeric" else f
}

.copyMethodDefaults <- function(generic, method)
{
    emptyDefault <- function(value) missing(value) ||
    (is.name(value) && nzchar(as.character(value)) )
    fg <- formals(generic)
    mg <- formals(method)
    mgn <- names(mg)
    changed <- FALSE
    for(what in names(fg)) {
        i <- match(what, mgn, 0L)
        if(i > 0L) {
            deflt <- mg[[i]]
            if(!(emptyDefault(deflt) || identical(deflt, fg[[what]]))) {
                fg[[what]] <- deflt
                changed <- TRUE
            }
        }
    }
    if(changed)
        formals(generic) <- fg
    generic
}

.NamespaceOrPackage <- function(what)
{
    name <- as.name(what)
    ns <-  .getNamespace(name)
    if(!is.null(ns))
        asNamespace(ns)
    else {
        i <- match(paste("package", what, sep=":"), search())
        if(is.na(i))
            .GlobalEnv
        else
            as.environment(i)
    }
}

.NamespaceOrEnvironment <- function(where)
{
    value <- NULL
    if(is.environment(where))
        value <- where
    else if(is.character(where) && nzchar(where)) {
        ns <- .getNamespace(where)
        if(isNamespace(ns))
            value <- ns
        else if(where %in% search())
            value <- as.environment(where)
        else {
            where <- paste0("package:", where)
            if(where %in% search())
                value <- as.environment(where)
        }
    }
    else if(is.numeric(where) && where %in% seq_along(search()))
        value <- as.environment(where)
    value
}


.hasS4MetaData <- function(env)
  (length(objects(env, all.names = TRUE,
                          pattern = "^[.]__[CTA]_")))

## turn ordinary generic into one that dispatches on "..."
## currently only called in one place from setGeneric()
.dotsGeneric <- function(f)
{
    if(!is(f, "genericFunction"))
        f <- getGeneric(f)
    if(!is(f, "genericFunction") || !identical(f@signature, "..."))
        stop("argument f must be a generic function with signature \"...\"")
    def <- .standardGenericDots
    fenv <- environment(f)
    environment(def) <- fenv
    assign("standardGeneric", def, envir = fenv)
    assign(".dotsCall", .makeDotsCall(formalArgs(f)), envir = fenv)
    f
}

utils::globalVariables(c(".MTable", ".AllMTable", ".dotsCall"))

.standardGenericDots <- function(name)
{
    env <- sys.frame(sys.parent())
    dots <- eval(quote(list(...)), env)
    classes <- unique(unlist(lapply(dots, methods:::.class1)))
    method <- methods:::.selectDotsMethod(classes, .MTable, .AllMTable)
    if(is.null(method))
        stop(gettextf("no method or default matching the \"...\" arguments in %s",
                      deparse(sys.call(sys.parent()), nlines = 1)), domain = NA)
    assign(".Method", method, envir = env)
    eval(.dotsCall, env)
}


.quoteCall <- quote(.Method(...))
.makeDotsCall <- function(formals)
{
    call <- methods:::.quoteCall
    if(length(formals)  > 1L) {
        idots <- match("...", formals)
        for(what in formals[-idots]) {
            ## the following nonsense is required to get the names in the call
            ## expression to be empty for ... and there for other args
            eval(substitute(call$NAME <- as.name(WHAT),
                            list(NAME = as.name(what), WHAT = what)))
        }
    }
    call
}

.selectDotsMethod <- function(classes, mtable, allmtable)
{
    .pasteC <- function(names) paste0('"', names, '"', collapse = ", ")
    found <- character()
    distances <- numeric()
    methods <- objects(mtable, all.names = TRUE)
    direct <- match(classes, methods, 0L) > 0L
    if(all(direct)) {
        if(length(classes) > 1L) {
            warning(gettextf("multiple direct matches: %s; using the first of these", .pasteC(classes)), domain = NA)
            classes <- classes[1L]
        }
        else if(length(classes) == 0L)
            return( if(is.na(match("ANY", methods))) NULL else get("ANY", envir = mtable))
        return(get(classes,envir = mtable))
    }
    if(is.null(allmtable))
        return(NULL)

    ## Else, look for an acceptable inherited method, which must match or be a superclass
    ## of the class of each of the arguments.
    classes <- sort(classes) # make slection depend only on the set of classes
    label <- .sigLabel(classes)
    if(exists(label, envir = allmtable, inherits = FALSE))
        ## pre-cached, but possibly NULL to indicate no match
        return(get(label, envir = allmtable))
    for(i in seq_along(classes)) {
        classi <- classes[[i]]
        defi <- getClassDef(classi)
        if(is.null(defi)) next
        extendsi <- defi@contains
        namesi <- c(classi, names(extendsi))
        if(i == 1)
            namesi <- namesi[match(namesi, methods, 0L) > 0L]
        else { # only the superclass methods matching all arguments are kept
            namesi <- namesi[match(namesi, found, 0L) > 0L]
            found <- namesi
            if(length(found) == 0L) break # no possible non-default match
        }
        for(namei in namesi) {
            disti <- if(identical(namei, classi)) 0 else extendsi[[namei]]@distance
            prev <- match(namei, found)
            if(is.na(prev)) {           # must be the 1st element
                found <- c(found, namei)
                distances <- c(distances, disti)
            }
            else if(disti < distances[[prev]])
                distances[[prev]] <- disti
        }
    }
    if(length(found) == 0L)
        method <-  if(is.na(match("ANY", methods))) NULL else get("ANY", envir = mtable)
    else {
        classes <- found[which.min(distances)]
        if(length(classes) > 1L) {
            warning(gettextf("multiple equivalent inherited matches: %s; using the first of these",
                             .pasteC(classes)), domain = NA)
            classes <- classes[1L]
        }
        method <- get(classes,envir = mtable)
    }
    if(!is.null(method))
        method@target <- new("signature", ... = label) # ?? not a legal class name if > 1 classes
    assign(label, method, allmtable)
    method
}

.isSingleString <- function(what)
  is.character(what) && identical(nzchar(what), TRUE)

.notSingleString <- function(what)
{
    if(identical(what, ""))
        "non-empty string; got \"\""
    else if(is.character(what))
        paste("single string; got a character vector of length", length(what))
    else
        gettextf("single string; got an object of class %s",
                 dQuote(class(what)[[1L]]))
}

.dotsClass <- function(...) {
    if(missing(..1))
      "missing"
    else
      class(..1)
}

## a utility to exclude various annoying glitches during
## loading of the methods package
.methodsIsLoaded <- function()
    identical(.saveImage, TRUE)

if(FALSE) {
## Defined but not currently used:
## utilitity to test well-defined classes in signature,
## for setMethod(), setAs() [etc.?], the result to be
## assigned in package where=
## Returns a list of signature, messages and level of error

## Has undefined ns an package
 .validSignature <- function(signature, generic, where) {
    thisPkg <- getPackageName(where, FALSE)
    checkDups <- .duplicateClassesExist()
    if(is(signature, "character")) { # including class "signature"
        classes <- as.character(signature)
        names <- allNames(signature)
        pkgs <- attr(signature, "package")
    }
    else if(is(signature, "list")) {
        classes <- sapply(signature, as.character)
        names <- names(signature)
        pkgs <- character(length(signature))
        for(i in seq_along(pkgs)) {
            pkgi <- attr(signature[[i]], "package")
            pkgs[[i]] <- if(is.null(pkgi)) "" else pkgi
        }
    }
    msgs <- character(); level <- integer()
    for(i in seq_along(classes)) {
        ## classes must be defined
        ## if duplicates exist check for them
        ## An ambiguous duplicate is a warning if it can match thisPkg
        ## else, an error
        classi <- classes[[i]]
        pkgi <- pkgs[[i]]
        classDefi <- getClass(classi, where = where)
        if(checkDups &&
           classi %in% mulipleClasses()) { # hardly ever, we hope
            clDefsi <- get(classi, envir = .classTable)
            if(nzchar(pkgi) && pkgi %in% names(clDefsi))
                ## use the chosen class, no message
                classDefi <- clDefsi[[pkgi]]
            else if(nzchar(pkgi)){
                ## this is only a warning because it just might
                ## be the result of identical class defs (e.g., from setOldClass()
                msgs <- c(msgs,
                          gettextf("multiple definitions exist for class %s, but the supplied package (%s) is not one of them (%s)",
                                   dQuote(classi), sQuote(pkgi),
                                   paste(dQuote(get(classi, envir = .classTable)), collapse = ", ")))
                level <- c(level, 2) #warn
            }
            else {
                msgs <- c(msgs,
                          gettextf("multiple definitions exist for class %s; should specify one of them (%s), e.g. by className()",
                                   dQuote(classi),
                                   paste(dQuote(get(classi, envir = .classTable)), collapse = ", ")))
            }
        }
        else {
            ## just possibly the first reference to an available
            ## package not yet loaded.  It's an error to specify
            ## a non-loadable package
            if(nzchar(pkgi)) {
                loadNamespace(pkgi)
                classDefi <- getClass(classi, where = ns)
            }
            if(is.null(classDefi)) {
                classDefi <- getClassDef
                msgi <- gettextf("no definition found for class %s",
                                 dQuote(classi))
                ## ensure only one error message
                if(length(level) && any(level == 3))
                    msgs[level == 3] <- paste(msgs[level == 3], msgi, sep = "; ")
                else
                    msgs <- c(msgs, msgi)
                level <- c(level, 3)
            }
            ## note that we do not flag a pkgi different from
            ## the package of the def., mainly because of setOldClass()
            ## which currently generates potentially multiple versions
            ## of the same S3 class.
        }
        ## except for the obscure multiple identical class case
        ## we should not get here w/o a valid class def.
        if(is.null(classDefi)) {}
        else
            pkgs[[i]] <- classDefi@package
    }
    signature <- .MakeSignature(new("signature"), generic,
                                structure(classes, names = names, package = package))
    if(length(msgs) > 1) {
        ## sort by severity, to get all messages before errror
        ii <- sort.list(level)
        msgs <- msgs[ii]; level <- level[ii]
    }
    list(signature = signature, message = msgs, level = level)
}
}

.ActionMetaPattern <- function()
    paste0("^[.]",substring(methodsPackageMetaName("A",""),2))

.actionMetaName <- function(name)
    methodsPackageMetaName("A", name)


.doLoadActions <- function(where, attach) {
    ## at the moment, no unload actions
    if(!attach)return()
    actionListName <- .actionMetaName("")
    if(!exists(actionListName, envir = where, inherits = FALSE))
        return(list())
    actions <- get(actionListName, envir = where)
    ## check sanity:  methods must be loaded
    if(! "package:methods" %in% search()) {
        warning("trying to execute load actions without 'methods' package")
        library(methods)
    }
    for(what in actions) {
        aname <- .actionMetaName(what)
        if(!exists(aname, envir = where, inherits = FALSE)) {
            warning(gettextf("missing function for load action: %s", what))
            next
        }
        f <- get(aname, envir = where)
        value <- eval(substitute(tryCatch(FUN(WHERE), error = function(e)e),
                            list(FUN = f, WHERE = where)), where)
        if(is(value, "error")) {
            callString <- deparse(value$call)[[1]]
            stop(gettextf("error in load action %s for package %s: %s: %s",
                          aname, getPackageName(where), callString, value$message))
        }
    }
}

setLoadAction <- function(action,
              aname = "",
              where = topenv(parent.frame())) {
    currentAnames <- .assignActionListNames(where)
    if(!nzchar(aname))
        aname <- paste0(".", length(currentAnames)+1)
    .assignActions(list(action), aname, where)
    if(is.na(match(aname, currentAnames))) {
        actionListName <- .actionMetaName("")
        assign(actionListName, c(currentAnames, aname), envir = where)
    }
}

.assignActions <- function(actions, anames, where) {
    ## check all the actions before assigning any
    for(i in seq_along(actions)) {
        f <- actions[[i]]
        fname <- anames[[i]]
        if(!is(f, "function"))
            stop(gettextf("non-function action: %s",
                          sQuote(fname)),
                 domain = NA)
        if(length(formals(f)) == 0)
            stop(gettextf("action function %s has no arguments, should have at least 1",
                          sQuote(fname)),
                 domain = NA)
    }
    for(i in seq_along(actions))
        assign(.actionMetaName(anames[[i]]), actions[[i]], envir = where)
}

.assignActionListNames <- function(where) {
    actionListName <- .actionMetaName("")
    if(exists(actionListName, envir = where, inherits = FALSE))
        get(actionListName, envir = where)
    else
        character()
}

setLoadActions <- function(..., .where = topenv(parent.frame())) {
    actionListName <- .actionMetaName("")
    currentAnames <- .assignActionListNames(.where)
    actions <- list(...)
    anames <- allNames(actions)
    ## first, replacements
    previous <- anames %in% currentAnames
    if(any(previous)) {
        .assignActions(actions[previous], anames[previous], .where)
        if(all(previous))
            return(list())
        anames <- anames[!previous]
        actions <- actions[!previous]
    }
    anon <- !nzchar(anames)
    if(any(anon)) {
        n <- length(currentAnames)
        deflts <- paste0(".",seq(from = n+1, length.out = length(actions)))
        anames[anon] <- deflts[anon]
    }
    .assignActions(actions, anames, .where)
    assign(actionListName, c(currentAnames, anames), envir = .where)
}

hasLoadAction <- function(aname, where = topenv(parent.frame()))
    exists(.actionMetaName(aname), envir = where, inherits = FALSE)

getLoadActions <- function(where = topenv(parent.frame())) {
    actionListName <- .actionMetaName("")
    if(!exists(actionListName, envir = where, inherits = FALSE))
        return(list())
    actions <- get(actionListName, envir = where)
    if(length(actions)) {
        allExists <- sapply(actions, function(what) exists(.actionMetaName(what), envir = where, inherits = FALSE))
        if(!all(allExists)) {
            warning(gettextf("some actions are missing: %s",
                             paste(actions[!allExists], collapse =", ")),
                    domain = NA)
            actions <- actions[allExists]
        }
        allFuns <- lapply(actions, function(what) get(.actionMetaName(what), envir = where))
        names(allFuns) <- actions
        allFuns
    }
    else
        list()
}

evalOnLoad <- function(expr, where = topenv(parent.frame()), aname = "") {
    f <- function(env)NULL
    body(f, where) <- substitute(eval(EXPR,ENV), list(EXPR = expr, ENV = where))
    setLoadAction(f, aname, where)
}

evalqOnLoad <- function(expr, where = topenv(parent.frame()), aname = "")
    evalOnLoad(substitute(expr), where, aname)

## a utility function used to flag non-generics at the loadNamespace phase
## The calculation there used to ignore the generic cache, which is wrong logic
## if the package being loaded had a DEPENDS on a package containing the generic
## version of the function.
.findsGeneric <- function(what, ns) {
    if(is(get(what, mode = "function", envir = ns), "genericFunction"))
        1L
    else if(!is.null(.getGenericFromCache(what, ns)))
        2L
    else
        0L
}

#  File src/library/methods/R/SClasses.R
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

setClass <-
    ## Define Class to be an S4 class.
    function(Class, representation = list(), prototype = NULL,
             contains = character(), validity = NULL, access = list(),
             where = topenv(parent.frame()), version = .newExternalptr(),
             sealed = FALSE, package = getPackageName(where),
             S3methods = FALSE, slots)
{
    oldDef <- getClassDef(Class, where)
    if(is(oldDef, "classRepresentation") && oldDef@sealed)
        stop(gettextf("%s has a sealed class definition and cannot be redefined",
                      dQuote(Class)),
             domain = NA)
    if(!missing(slots)) {
        ## The modern version consistent with reference classes
        ## Arguments slots= and contains= are used, representation must not be
        if(!missing(representation))
            stop("Argument \"representation\" cannot be used if argument \"slots\" is supplied")
        properties <- inferProperties(slots, "slot")
        classDef <- makeClassRepresentation(Class, properties,contains, prototype, package,
                                             validity, access, version, sealed, where = where)
        superClasses <- names(classDef@contains)
    }
    else if(is(representation, "classRepresentation")) {
        ## supplied a class definition object
        classDef <- representation
        if(!(missing(prototype) && missing(contains) && missing(validity) && missing(access)
             && missing(version) && missing(package)))
            stop("only arguments 'Class' and 'where' can be supplied when argument 'representation' is a 'classRepresentation' object")
        if(length(classDef@package) == 0L)
            classDef@package <- package # the default
        superClasses <- allNames(classDef@contains)
    }
    else {
        ## catch the special case of a single class name as the representation
        if(is.character(representation) && length(representation) == 1L &&
           is.null(names(representation)))
            representation <- list(representation)
        slots <- nzchar(allNames(representation))
        superClasses <- c(as.character(representation[!slots]), contains)
        properties <- representation[slots]
        classDef <- makeClassRepresentation(Class, properties,superClasses, prototype, package,
                                             validity, access, version, sealed, where = where)
        superClasses <- names(classDef@contains)
    }
    classDef <- completeClassDefinition(Class, classDef, where, doExtends = FALSE)
    ## uncache an old definition for this package, if one is cached
    .uncacheClass(Class, classDef)
    if(length(superClasses) > 0L) {
        sealed <- classDef@sealed
        classDef@sealed <- FALSE # to allow setIs to work anyway; will be reset later
        assignClassDef(Class, classDef, where)
        badContains <- character()
        for(class2 in superClasses) {
            if(is(try(setIs(Class, class2, classDef = classDef, where = where)), "try-error"))
                badContains <- c(badContains, class2)
            else { # update class definition
                classDef <- getClassDef(Class, where = where)
                if(is.null(classDef))
                  stop(sprintf("internal error: definition of class %s not properly assigned",
                                dQuote(Class)),
                       domain = NA)
            }
          }
        if(length(badContains)) {
            msg <- paste(.dQ(badContains), collapse = ", ")
            if(is(try(removeClass(Class, where)), "try-error"))
                stop(gettextf("error in contained classes (%s) for class %s and unable to remove definition from %s",
                              msg, dQuote(Class),
                              sQuote(getPackageName(where))),
                     domain = NA)
            if(is.null(oldDef))
                stop(gettextf("error in contained classes (%s) for class %s; class definition removed from %s",
                              msg, dQuote(Class),
                              sQuote(getPackageName(where))),
                     domain = NA)
            else if(is(try(setClass(Class, oldDef, where=where)), "try-error"))
                stop(gettextf("error in contained classes (%s) for class %s and unable to restore previous definition from %s",
                              msg, dQuote(Class),
                              sQuote(getPackageName(where))),
                     domain = NA)
            else
                stop(gettextf("error in contained classes (%s) for class %s; previous definition restored to %s",
                              msg, dQuote(Class),
                              sQuote(getPackageName(where))),
                     domain = NA)
        }
        if(length(attr(classDef@contains, "conflicts")) > 0)
          .reportSuperclassConflicts(Class, classDef@contains, where)
        .checkRequiredGenerics(Class, classDef, where)
        if(sealed) {
            classDef@sealed <- TRUE
        }
    }
    if(S3methods)
      classDef <- .setS3MethodsOn(classDef)
    assignClassDef(Class, classDef, where)
    invisible(classGeneratorFunction(classDef, where))
}

representation <-
  ## Representation of a class; that is,
  ## a list of named slots and unnamed classes to be included in a class
  ## definition.
  function(...)
{
    value <- list(...)
    ## unlike the S-Plus function, this does not form the class representation,
    ## since set SClass works separately with the slots and extends arguments.
    anames <- allNames(value)
    for(i in seq_along(value)) {
        ei <- el(value, i)
        if(!is.character(ei) || length(ei) != 1L)
            stop(gettextf("element %d of the representation was not a single character string", i), domain = NA)
    }
    includes <- as.character(value[!nzchar(anames)])
    if(anyDuplicated(includes))
        stop(gettextf("duplicate class names among superclasses: %s",
                      paste(.dQ(includes[duplicated(includes)]),
                            collapse = ", ")),
             domain = NA)
    slots <- anames[nzchar(anames)]
    if(anyDuplicated(slots)) {
        dslots <- slots[duplicated(slots)]
        stop(sprintf(ngettext(length(dslots),
                              "duplicated slot name: %s",
                              "duplicated slot names: %s"),
                     paste(sQuote(dslots), collapse="")),
             domain = NA)
    }
    value
}

### the version called prototype is the external interface.  But functions with argument
### named prototype in R cannot call the prototype function (until there is a methods namespace
### to allow methods::prototype(...)
prototype <- function(...)
    .prototype(...)

.prototype <- function(...) {
    props <- list(...)
    names <- allNames(props)
    data <- !nzchar(names)
    dataPart <- any(data)
    if(dataPart) {
        if(sum(data) > 1)
            stop("only one data object (unnamed argument to prototype) allowed")
        obj <- unclass(props[[seq_along(data)[data] ]])
        props <- props[!data]
        names <- names[!data]
    }
    else
        obj <- defaultPrototype()
    for(i in seq_along(names))
        slot(obj, names[[i]], FALSE) <- props[[i]]
    new("classPrototypeDef", object = obj, slots = names, dataPart = dataPart)
}

makeClassRepresentation <-
  ## Set the Class Definition.
  ## The formal definition of the class is set according to the arguments.
  ##
  ## Users should call setClass instead of this function.
  function(name, slots = list(), superClasses = character(), prototype = NULL, package, validity = NULL, access = list(), version = .newExternalptr(), sealed = FALSE, virtual = NA, where)
{
    if(any(superClasses %in% .AbnormalTypes))
        superClasses <- .addAbnormalDataType(superClasses)
    if(!is.null(prototype) || length(slots) || length(superClasses)) {
        ## collect information about slots, create prototype if needed
        pp <- reconcilePropertiesAndPrototype(name, slots, prototype, superClasses, where)
        slots <- pp$properties
        prototype <- pp$prototype
    }
    contains <- list()
    if(nzchar(package))
        packageSlot(name) <- package
    for(what in superClasses) {
        if(is(what, "classRepresentation"))
            whatClassDef <- what
        else if(is.null(packageSlot(what)))
            whatClassDef <- getClass(what, where = where)
        else
            whatClassDef <- getClass(what)
        what <- whatClassDef@className # includes package name as attribute
        ## Create the SClassExtension objects (will be simple, possibly dataPart).
        ## The slots are supplied explicitly, since `name' is currently an undefined class
        elNamed(contains, what) <- makeExtends(name, what, slots = slots,
                                              classDef2 = whatClassDef, package = package)
    }
    validity <- .makeValidityMethod(name, validity)
    if(is.na(virtual)) {
        virtual <- testVirtual(slots, contains, prototype, where)
        if(virtual && !is.na(match("VIRTUAL", superClasses)))
            elNamed(contains, "VIRTUAL") <- NULL
    }
    # new() must return an S4 object, except perhaps for basic classes
    if(!is.null(prototype) && is.na(match(name, .BasicClasses)))
      prototype <- .asS4(prototype)
    if(".S3Class" %in% names(slots))
      prototype <- .addS3Class(name, prototype, contains, where)
    newClassRepresentation(className = name, slots = slots,
                           contains = contains,
                           prototype = prototype,
                           virtual = virtual,
                           validity = validity,
                           access = access,
                           package = package,
                           versionKey = version,
                           sealed = sealed)
}

getClassDef <-
  ## Get the definition of the class supplied as a string.
  function(Class, where = topenv(parent.frame()), package = packageSlot(Class),
           inherits = TRUE)
{
    if(inherits) #includes both the lookup and Class being alread a definition
      value <- .getClassFromCache(Class, where)
    else # want to force a search for the metadata in this case (Why?)
      value <- NULL
    if(is.null(value)) {
	cname <-
	    classMetaName(if(length(Class) > 1L)
			  ## S3 class; almost certainly has no packageSlot,
			  ## but we'll continue anyway
			  Class[[1L]] else Class)
	## a string with a package slot strongly implies the class definition
	## should be in that package.
	if(identical(nzchar(package), TRUE)) {
	    whereP <- .requirePackage(package)
	    if(exists(cname, whereP, inherits = inherits))
		value <- get(cname, whereP)
	}
	if(is.null(value) && exists(cname, where, inherits = inherits))
	    value <- get(cname, where)
    }
    value
}

getClass <-
  ## Get the complete definition of the class supplied as a string,
  ## including all slots, etc. in classes that this class extends.
  function(Class, .Force = FALSE,
	   where = .classEnv(Class, topenv(parent.frame()), FALSE))
{
    value <- .getClassFromCache(Class, where) # the quick way
    if(is.null(value)) {
        value <- getClassDef(Class, where) # searches
        if(is.null(value)) {
            if(!.Force)
                stop(gettextf("%s is not a defined class",
                              dQuote(Class)),
                     domain = NA)
            else
                value <- makeClassRepresentation(Class, package = "base",
                                                 virtual = TRUE, where = where)
        }
    }
    value
}

slot <-
  ## Get the value of the named slot.  This function does exact, not partial, matching of names,
  ## and the name must be one of the slot names specified in the class's definition.
  ##
  ## Because slots are stored as attributes, the validity check is not 100% guaranteed,
  ## but should be OK if nobody has "cheated" (e.g., by setting other attributes directly).
  function(object, name)
    .Call(C_R_get_slot, object, name)

"slot<-" <-
  ## Set the value of the named slot.  Must be one of the slots in the class's definition.
  function(object, name, check = TRUE, value) {
      if(check)
          value <- checkSlotAssignment(object, name, value)
      .Call(C_R_set_slot, object, name, value)
      ## currently --> R_do_slot_assign() in ../../../main/attrib.c
  }

## ". - hidden" since one should typically rather use is(), extends() etc:
.hasSlot <- function(object, name)
    .Call(C_R_hasSlot, object, name)

checkSlotAssignment <- function(obj, name, value)
{
    cl <- class(obj)
    ClassDef <- getClass(cl) # fails if cl not a defined class (!)
    slotClass <- elNamed(ClassDef@slots, name)
    if(is.null(slotClass))
        stop(gettextf("%s is not a slot in class %s",
                      sQuote(name), dQuote(cl)),
             domain = NA)
    valueClass <- class(value)
    if(.identC(slotClass, valueClass))
       return(value)
    ## check the value, but be careful to use the definition of the slot's class from
    ## the class environment of obj (change validObject too if a better way is found)
    ok <- possibleExtends(valueClass, slotClass,
                          ClassDef2 = getClassDef(slotClass, where = .classEnv(ClassDef)))
    if(identical(ok, FALSE))
       stop(gettextf("assignment of an object of class %s is not valid for slot %s in an object of class %s; is(value, \"%s\") is not TRUE",
		     dQuote(valueClass), sQuote(name), dQuote(cl), slotClass),
            domain = NA)
    else if(identical(ok, TRUE))
        value
    else
       as(value, slotClass, strict=FALSE, ext = ok)
}

## slightly simpler verison to be called from do_attrgets()
checkAtAssignment <- function(cl, name, valueClass)
{
    ClassDef <- getClass(cl) # fails if cl not a defined class (!)
    slotClass <- elNamed(ClassDef@slots, name)
    if(is.null(slotClass))
        stop(gettextf("%s is not a slot in class %s",
                      sQuote(name), dQuote(cl)),
             domain = NA)
    if(.identC(slotClass, valueClass))
       return(TRUE)
    ## check the value, but be careful to use the definition of the slot's class from
    ## the class environment of obj (change validObject too if a better way is found)
    ok <- possibleExtends(valueClass, slotClass,
                          ClassDef2 = getClassDef(slotClass, where = .classEnv(ClassDef)))
    if(identical(ok, FALSE))
       stop(gettextf("assignment of an object of class %s is not valid for @%s in an object of class %s; is(value, \"%s\") is not TRUE",
		     dQuote(valueClass), sQuote(name), dQuote(cl), slotClass),
            domain = NA)
    TRUE
}

## Now a primitive in base
## "@<-" <-
##    function(object, name, value) {
##      arg <- substitute(name)
##      if(is.name(arg))
##        name <- as.character(arg)
##      "slot<-"(object, name, TRUE, value)
##    }

##  The names of the class's slots.  The argument is either the name
##  of a class, or an object from the relevant class.

## NOTA BENE:  .slotNames() shouldn't be needed,
##             rather slotNames() should be changed (to work like .slotNames())!
slotNames <- function(x)
    if(is(x, "classRepresentation")) names(x@slots) else .slotNames(x)

.slotNames <- function(x)
{
    classDef <-
	getClassDef(if(is.character(x) && length(x) == 1L) x else class(x))
    if(is.null(classDef))
	character()
    else
	names(classDef@slots)
}


removeClass <-  function(Class, where = topenv(parent.frame())) {
    if(missing(where)) {
       classEnv <- .classEnv(Class, where, FALSE)
        classWhere <- findClass(Class, where = classEnv)
        if(length(classWhere) == 0L) {
            warning(gettextf("class definition for %s not found (no action taken)",
                             dQuote(Class)),
                    domain = NA)
            return(FALSE)
        }
        if(length(classWhere) > 1L)
            warning(gettextf("class %s has multiple definitions visible; only the first removed",
                             dQuote(Class)),
                    domain = NA)
        classWhere <- classWhere[[1L]]
    }
    else classWhere <- where
    classDef <- getClassDef(Class, where=classWhere)
    if(length(classDef@subclasses)) {
      subclasses <- names(classDef@subclasses)
      found <- sapply(subclasses, isClass, where = where)
      for(what in subclasses[found])
          .removeSuperClass(what, Class)
    }
    .removeSuperclassBackRefs(Class, classDef, classWhere)
    .uncacheClass(Class, classDef)
    .undefineMethod("initialize", Class, classWhere)
    what <- classMetaName(Class)
    rm(list=what, pos=classWhere)
    TRUE
}


isClass <-
  ## Is this a formally defined class?
  function(Class, formal=TRUE, where = topenv(parent.frame()))
    ## argument formal is for Splus compatibility & is ignored.  (All classes that
    ## are defined must have a class definition object.)
    !is.null(getClassDef(Class, where))

### TODO   s/Class/._class/  -- in order to allow 'Class' as regular slot name
new <-
  ## Generate an object from the specified class.
  ##
  ## Note that the basic vector classes, `"numeric"', etc. are implicitly defined,
  ## so one can use `new' for these classes.
  ##
  function(Class, ...)
{
    ClassDef <- getClass(Class, where = topenv(parent.frame()))
    value <- .Call(C_new_object, ClassDef)
    initialize(value, ...)
}

getClasses <-
  ## The names of all the classes formally defined on `where'.
  ## If called with no argument, all the classes currently known in the session
  ## (which does not include classes that may be defined on one of the attached
  ## libraries, but have not yet been used in the session).
  function(where = .externalCallerEnv(), inherits = missing(where))
{
    pat <- paste0("^",classMetaName(""))
    if(inherits) {
        evList <- .parentEnvList(where)
        clNames <- character()
        for(ev in evList)
            clNames <- c(clNames, objects(ev, pattern = pat, all.names = TRUE))
        clNames <- unique(clNames)
    }
    else
        clNames <- objects(where, pattern = pat, all.names = TRUE)
    ## strip off the leading pattern (this implicitly assumes the characters
    ## in classMetaName("") are either "." or not metacharacters
    substring(clNames, nchar(pat, "c"))
}


validObject <- function(object, test = FALSE, complete = FALSE)
{
    Class <- class(object)
    classDef <- getClassDef(Class)
    where <- .classEnv(classDef)
    anyStrings <- function(x) if(identical(x, TRUE)) character() else x
    ## perform, from bottom up, the default and any explicit validity tests
    ## First, validate the slots.
    errors <- character()
    slotTypes <- classDef@slots
    slotNames <- names(slotTypes)
    attrNames <- c(".Data", ".S3Class", names(attributes(object)))
    if(any(is.na(match(slotNames, attrNames)))) {
        badSlots <- is.na(match(slotNames, attrNames))
	errors <-
	    c(errors,
	      paste("slots in class definition but not in object:",
		    paste0('"', slotNames[badSlots], '"', collapse = ", ")))
        slotTypes <- slotTypes[!badSlots]
        slotNames <- slotNames[!badSlots]
    }
    for(i in seq_along(slotTypes)) {
	classi <- slotTypes[[i]]
	classDefi <- getClassDef(classi, where = where)
	if(is.null(classDefi)) {
	    errors <- c(errors,
			paste0("undefined class for slot \"", slotNames[[i]],
			       "\" (\"", classi, "\")"))
	    next
	}
        namei <- slotNames[[i]]
        sloti <- try(switch(namei,
                            ## .S3Class for S3 objects (e.g., "factor")
                            .S3Class = S3Class(object),
                            slot(object, namei)
                            ), silent = TRUE)
        if(inherits(sloti, "try-error")) {
           errors <- c(errors, sloti)
           next
        }
	## note that the use of possibleExtends is shared with checkSlotAssignment(), in case a
	## future revision improves on it!
	ok <- possibleExtends(class(sloti), classi, ClassDef2 = classDefi)
	if(identical(ok, FALSE)) {
	    errors <- c(errors,
			paste0("invalid object for slot \"", slotNames[[i]],
			       "\" in class \"", Class,
			       "\": got class \"", class(sloti),
			       "\", should be or extend class \"", classi, "\""))
	    next
	}
	if(!complete)
          next
        errori <- anyStrings(Recall(sloti, TRUE, TRUE))
        if(length(errori)) {
	    errori <- paste0("In slot \"", slotNames[[i]],
			     "\" of class \"", class(sloti), "\": ", errori)
            errors <- c(errors, errori)
        }
    }
    extends <- rev(classDef@contains)
    for(i in seq_along(extends)) {
	exti <- extends[[i]]
	superClass <- exti@superClass
	if(!exti@simple && !is(object, superClass))
	    next ## skip conditional relations that don't hold for this object
	superDef <- getClassDef(superClass, where = where)
	if(is.null(superDef)) {
	    errors <- c(errors,
			paste0("superclass \"", superClass,
			       "\" not defined in the environment of the object's class"))
	    break
	}
	validityMethod <- superDef@validity
	if(is(validityMethod, "function")) {
	    errors <- c(errors, anyStrings(validityMethod(as(object, superClass))))
	    if(length(errors))
		break
	}
    }
    validityMethod <- classDef@validity
    if(length(errors) == 0L && is(validityMethod, "function")) {
	errors <- c(errors, anyStrings(validityMethod(object)))
    }
    if(length(errors)) {
	if(test)
	    errors
	else {
	    msg <- gettextf("invalid class %s object", dQuote(Class))
	    if(length(errors) > 1L)
		stop(paste(paste0(msg, ":"),
                           paste(seq_along(errors), errors, sep=": "),
			   collapse = "\n"), domain = NA)
	    else stop(msg, ": ", errors, domain = NA)
	}
    }
    else
	TRUE
}

setValidity <- function(Class, method, where = topenv(parent.frame())) {
    if(isClassDef(Class)) {
	ClassDef <- Class
	Class <- ClassDef@className
    }
    else {
	ClassDef <- getClassDef(Class, where)
    }
    method <- .makeValidityMethod(Class, method)
    if(is.null(method) ||
       (is(method, "function") && length(formalArgs(method)) == 1L))
	ClassDef@validity <- method
    else
	stop("validity method must be NULL or a function of one argument")
    ## TO DO:  check the where argument against the package of the class def.
    assignClassDef(Class, ClassDef, where = where)
    resetClass(Class, ClassDef, where = where)
}

getValidity <- function (ClassDef) {
    ## "needed" according to ../man/validObject.Rd
    ClassDef@validity
}


resetClass <- function(Class, classDef, where) {
        if(is(Class, "classRepresentation")) {
            classDef <- Class
            Class <- Class@className
            if(missing(where))
                where <- .classDefEnv(classDef)
        }
        else {
            if(missing(where)) {
                if(missing(classDef))
                    where <- findClass(Class, unique = "resetting the definition")[[1L]]
                else
                    where <- .classDefEnv(classDef)
            }
            if(missing(classDef)) {
                classDef <- getClassDef(Class, where)
                if(is.null(classDef)) {
                    warning(gettextf("class %s not found on %s; 'resetClass' will have no effect",
                                     dQuote(Class),
                                     sQuote(getPackageName(where))),
                            domain = NA)
                    return(classDef)
                }
            }
            else if(!is(classDef, "classRepresentation"))
                stop(gettextf("argument 'classDef' must be a string or a class representation; got an object of class %s",
                              dQuote(class(classDef))),
                     domain = NA)
            package <- getPackageName(where)
        }
        if(classDef@sealed)
            warning(gettextf("class %s is sealed; 'resetClass' will have no effect",
                             dQuote(Class)),
                    domain = NA)
        else {
            classDef <-  .uncompleteClassDefinition(classDef)
            classDef <- completeClassDefinition(Class, classDef, where)
            assignClassDef(Class, classDef, where)
        }
        classDef
    }

## the (default) initialization:  becomes the default method when the function
## is made a generic by .InitMethodDefinitions

initialize <- function(.Object, ...) {
    args <- list(...)
    if(length(args)) {
        Class <- class(.Object)
        ## the basic classes have fixed definitions
        if(!is.na(match(Class, .BasicClasses)))
            return(newBasic(Class, ...))
        ClassDef <- getClass(Class)
        ## separate the slots, superclass objects
        snames <- allNames(args)
        which <- nzchar(snames)
        elements <- args[which]
        supers <- args[!which]
        thisExtends <- names(ClassDef@contains)
        slotDefs <- ClassDef@slots
        dataPart <- elNamed(slotDefs, ".Data")
        if(is.null(dataPart)) dataPart <- "missing"
        if(length(supers)) {
            for(i in rev(seq_along(supers))) {
                obj <- el(supers, i)
                Classi <- class(obj)
                if(length(Classi) > 1L)
                    Classi <- Classi[[1L]] #possible S3 inheritance
                ## test some cases that let information be copied into the
                ## object, ordered from more to less:  all the slots in the
                ## first two cases, some in the 3rd, just the data part in 4th
                if(.identC(Classi, Class))
                    .Object <- obj
                else if(extends(Classi, Class))
                    .Object <- as(obj, Class, strict=FALSE)
                else if(extends(Class, Classi))
                    as(.Object, Classi) <- obj
                else if(extends(Classi, dataPart))
                    .Object@.Data <- obj
                else {
                    ## is there a class to which we can coerce obj
                    ## that is then among the superclasses of Class?
                    extendsi <- extends(Classi)[-1L]
                    ## look for the common extensions, choose the first
                    ## one in the extensions of Class
                    which <- match(thisExtends, extendsi)
                    which <- seq_along(which)[!is.na(which)]
                    if(length(which)) {
                        Classi <- thisExtends[which[1L]]
###                    was:    as(.Object, Classi) <- as(obj, Classi, strict = FALSE)
                        ## but   as<- does an as(....) to its value argument
                        as(.Object, Classi) <- obj
                    }
                    else
                        stop(gettextf("cannot use object of class %s in new():  class %s does not extend that class",
                                      dQuote(Classi),
                                      dQuote(Class)),
                             domain = NA)
                }
            }
        }
        if(length(elements)) {
            snames <- names(elements)
	    if(anyDuplicated(snames))
                stop(gettextf("duplicated slot names: %s",
                              paste(sQuote(snames[duplicated(snames)]),
                                    collapse = ", ")), domain = NA)
            which  <- match(snames, names(slotDefs))
            if(any(is.na(which)))
                stop(sprintf(ngettext(sum(is.na(which)),
                                      "invalid name for slot of class %s: %s",
                                      "invalid names for slots of class %s: %s"),
                              dQuote(Class),
                              paste(snames[is.na(which)], collapse=", ")),
                     domain = NA)
            firstTime <- TRUE
            for(i in seq_along(snames)) {
                slotName <- el(snames, i)
                slotClass <- elNamed(slotDefs, slotName)
                slotClassDef <- getClassDef(slotClass, package=ClassDef@package)
                slotVal <- el(elements, i)
                ## perform non-strict coercion, but leave the error messages for
                ## values not conforming to the slot definitions to validObject(),
                ## hence the check = FALSE argument in the slot assignment
                if(!.identC(class(slotVal), slotClass)
                   && !is.null(slotClassDef) ) {
                    valClass <- class(slotVal)
                    valClassDef <- getClassDef(valClass, package = ClassDef@package)
                    if(!identical(possibleExtends(valClass, slotClass,
                                         valClassDef, slotClassDef), FALSE))
                        slotVal <- as(slotVal, slotClass, strict = FALSE)
                }
                if (firstTime) {
                    ## force a copy of .Object
                    slot(.Object, slotName, check = FALSE) <- slotVal
                    firstTime <- FALSE
                } else {
                    ## XXX: do the assignment in-place
                    "slot<-"(.Object, slotName, check = FALSE, slotVal)
                }
            }
        }
        validObject(.Object)
     }
    .Object
}

findClass <- function(Class, where = topenv(parent.frame()), unique = "") {
    if(is(Class, "classRepresentation")) {
        pkg <- Class@package
        classDef <- Class
        Class <- Class@className
    }
    else {
        pkg <- packageSlot(Class)
        if(is.null(pkg))
          pkg <- ""
        classDef <- getClassDef(Class, where, pkg)
    }
    if(missing(where) && nzchar(pkg))
            where <- .requirePackage(pkg)
    else
        where <- as.environment(where)
    what <- classMetaName(Class)
    where <- .findAll(what, where)
    if(length(where) > 1L && nzchar(pkg)) {
        pkgs <- sapply(where, function(db)get(what, db)@package)
        where <- where[match(pkg, pkgs, 0L)]
    }
    else
      pkgs <- pkg
    if(length(where) == 0L) {
        if(is.null(classDef))
            classDef <- getClassDef(Class) # but won't likely succeed over previous
        if(nzchar(unique)) {
            if(is(classDef, "classRepresentation"))
                stop(gettextf("class %s is defined, with package %s, but no corresponding metadata object was found (not exported?)",
                              dQuote(Class),
                              sQuote(classDef@package)),
                     domain = NA)
            else
                stop(gettextf("no definition of %s to use for %s",
                              dQuote(Class),
                              unique),
                     domain = NA)
        }
    }
    else if(length(where) > 1L) {
        pkgs <- sapply(where, getPackageName)
        where <- where[!duplicated(pkgs)]
        if(length(where) > 1L)
            if(nzchar(unique)) {
                pkgs <- base::unique(pkgs)
                where <- where[1L]
                ## problem: 'unique'x is text passed in, so do not translate
                warning(sprintf(ngettext(length(pkgs),
                                         "multiple definition of class %s visible (%s); using the definition\n   in package %s for %s",
                                         "multiple definitions of class %s visible (%s); using the definition\n   in package %s for %s"),
                                dQuote(Class),
                                paste(sQuote(pkgs), collapse = ", "),
                                sQuote(pkgs[[1L]]),
                                unique),
                        domain = NA)
            }
            ## else returns a list of >1 places, for the caller to sort out (e.g., .findOrCopyClass)
    }
    where
}

isSealedClass <- function(Class, where = topenv(parent.frame())) {
    if(is.character(Class))
            Class <- getClass(Class, TRUE, where)
    if(!is(Class, "classRepresentation"))
        FALSE
    else
        Class@sealed
}

sealClass <- function(Class, where = topenv(parent.frame())) {
    if(missing(where))
        where <- findClass(Class, unique = "sealing the class", where = where)
    classDef <- getClassDef(Class, where)
    if(!classDef@sealed) {
        classDef@sealed <- TRUE
        assignClassDef(Class, classDef, where)
    }
    invisible(classDef)
}

## see $RHOME/src/main/duplicate.c for the corresponding datatypes
## not copied by duplicate1
.AbnormalTypes <- c("environment", "name", "externalptr",  "NULL")


.indirectAbnormalClasses <- paste0(".", .AbnormalTypes)
names(.indirectAbnormalClasses) <- .AbnormalTypes

## the types not supported by indirect classes (yet)
.AbnormalTypes <- c(.AbnormalTypes,
                    "special","builtin", "weakref", "bytecode")

.addAbnormalDataType <- function(classes) {
  types <- match(classes, .AbnormalTypes, 0) > 0
  type = classes[types]
  if(length(type) == 0)
    return(classes)
  if(length(type) > 1)
    stop(gettextf("class definition cannot extend more than one of these data types: %s",
		  paste0('"',type, '"', collapse = ", ")),
         domain = NA)
  class <- .indirectAbnormalClasses[type]
  if(is.na(class))
    stop(gettextf("abnormal type %s is not supported as a superclass of a class definition",
                  dQuote(type)),
         domain = NA)
  ## this message USED TO BE PRINTED: reminds programmers that
  ## they will see an unexpected superclass
  ## message(gettextf('Defining type "%s" as a superclass via class "%s"',
  ##                 type, class), domain = NA)
  c(class, classes[!types])
}

.checkRequiredGenerics <- function(Class, classDef, where) {}

..checkRequiredGenerics <- function(Class, classDef, where) {
  ## If any of the superclasses are in the .NeedPrimitiveMethods
  ## list, cache the corresponding generics now and also save their names in
  ## .requireCachedGenerics to be used when the environment
  ## where= is loaded.
  supers <- names(classDef@contains)
  allNeeded <- get(".NeedPrimitiveMethods", envir = .methodsNamespace)
  specials <- names(allNeeded)
  needed <- match(specials, supers, 0) > 0
  if(any(needed)) {
    generics <- unique(allNeeded[needed])
    packages <- character()
    for(g in generics) {
      def <- getGeneric(g)
      packages <- c(packages, def@package) # must be "methods" ?
      cacheGenericsMetaData(g, def, TRUE, where, def@package)
    }
    if(exists(".requireCachedGenerics",  where, inherits = FALSE))
      previous <- get(".requireCachedGenerics",  where)
    else
      previous <- character()
    packages <- c(attr(previous, "package"), packages)
    gg <- c(previous, generics)
    attr(gg, "package") <- packages
    assign(".requireCachedGenerics", gg, where)
  }
}

.setS3MethodsOn <- function(classDef) {
    ext <- extends(classDef)
    slots <- classDef@slots
    if(is.na(match(".S3Class", names(slots)))) {
        ## add the slot if it's not there
        slots$.S3Class <- getClass("oldClass")@slots$.S3Class
        classDef@slots <- slots
    }
    ## in any case give the prototype the full extends as .S3Class
    proto <- classDef@prototype
    if(is.null(proto)) # simple virtual class--unlikely but valid
        proto <- defaultPrototype()
    attr(proto, ".S3Class") <- ext
    classDef@prototype <- proto
    classDef
  }

multipleClasses <- function(details = FALSE) {
    ctable <- .classTable
    cnames <- objects(ctable, all.names = TRUE)
    dups <- sapply(cnames, function(x) is.list(get(x, envir = ctable)))
    if(details) {
        value <- lapply(cnames[dups], function(x) get(x, envir = ctable))
        names(value) <- cnames[dups]
        value
    }
    else
        cnames[dups]
}

className <- function(class, package) {
    if(is(class, "character")) {
        className <- as.character(class)
        if(missing(package))
            package <- packageSlot(class)
        if(is.null(package)) {
            if(exists(className, envir = .classTable, inherits = FALSE))
                classDef <- get(className, envir = .classTable)
            else {
                classDef <- findClass(className, topenv(parent.frame()))
                if(length(classDef) == 1)
                    classDef <- classDef[[1]]
            }
            ## at this point, classDef is the definition if
            ## unique, otherwise a list of 0 or >1 definitions
            if(is(classDef, "classRepresentation"))
                package <- classDef@package
            else if(length(classDef) > 1L) {
                pkgs <- sapply(classDef, function(cl)cl@package)
                warning(gettextf("multiple class definitions for %s from packages: %s; picking the first",
                                 dQuote(className),
                                 paste(sQuote(pkgs), collapse = ", ")),
                        domain = NA)
                package <- pkgs[[1L]]
            }
            else
                stop(gettextf("no package name supplied and no class definition found for %s",
                              dQuote(className)),
                     domain = NA)
        }
    }
    else if(is(class, classDef)) {
        className <- class@className
        if(missing(package))
            package <- class@package
    }
    new("className", .Data = className, package = package)
}

## bootstrap version before the class is defined
classGeneratorFunction <- function(classDef, env = topenv(parent.frame())) {
    fun <- function(...)NULL
    ## put the class name with package attribute into new()
    body(fun) <- substitute(new(CLASS, ...),
                            list(CLASS = classDef@className))
    environment(fun) <- env
    fun
}

.classGeneratorFunction <- function(classDef, env = topenv(parent.frame())) {
    if(is(classDef, "classRepresentation")) {}
    else if(is(classDef, "character")) {
        if(is.null(packageSlot(classDef)))
            classDef <- getClass(classDef, where = env)
        else
            classDef <- getClass(classDef)
    }
    else
        stop("argument 'classDef' must be a class definition or the name of a class")
    fun <- function(...)NULL
    ## put the class name with package attribute into new()
    body(fun) <- substitute(new(CLASS, ...),
                            list(CLASS = classDef@className))
    environment(fun) <- env
    fun <- as(fun, "classGeneratorFunction")
    fun@className <- classDef@className
    fun@package <- classDef@package
    fun
}

## grammar: 'what' is an adjective, so not plural ....
inferProperties <- function(props, what) {
    .validPropNames <- function(propNames) {
        n <- length(props)
        if(!n)
            return(character())
        else if(is.null(propNames))
            stop(gettextf("No %s names supplied", what),
                 domain = NA, call. = FALSE)
        else if(!all(nzchar(propNames)))
            stop(gettextf("All %s names must be nonempty in:\n(%s)", what,
                          paste(sQuote(propNames), collapse = ", ")),
                 domain = NA, call. = FALSE)
        else if(any(duplicated(propNames))) # NB: not translatable because of plurals
            stop(gettextf("All %s names must be distinct in:\n(%s)", what,
                          paste(sQuote(propNames), collapse = ", ")),
                 domain = NA, call. = FALSE)
        propNames
    }
    if(is.character(props)) {
        propNames <- names(props)
        if(is.null(propNames)) {
            propNames <- .validPropNames(props) # the text is the names
            ## treat as "ANY"
            props <- as.list(rep("ANY", length(props)))
            names(props) <- propNames
        }
        else {
            .validPropNames(propNames)
            props <- as.list(props)
        }
    }
    else if(is.list(props)) {
        if(length(props) > 0) # just validate them
            .validPropNames(names(props))
    }
    else
        stop(gettextf("argument %s must be a list or a character vector; got an object of class %s",
                      dQuote(what), dQuote(class(fields))),
             domain = NA)
    props
}


#  File src/library/methods/R/addedFunctions.R
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


functionBody <- base::body #was get("body", mode = "function")

`functionBody<-` <- base::`body<-`
## was
## .ff <- function(fun, envir = environment(fun), value) fun
## body(.ff, envir = .GlobalEnv) <- body(get("body<-"))
## "functionBody<-" <- .ff
## rm(.ff)

allNames <-
  ## the character vector of names (unlike names(), never returns NULL)
  function(x)
{
    value <- names(x)
    if(is.null(value))
        character(length(x))
    else
        value
}

getFunction <- function(name, generic = TRUE, mustFind = TRUE,
                        where = topenv(parent.frame()))
      ## find the object as a function.
{
    if(!nzchar(name))
        stop(gettextf('expected a non-empty character string for argument name'), domain = NA)
    found <- FALSE
    where <- as.environment(where)
    f <- NULL
    ## parent.env sequence of a namespace ends in the base package namespace,
    ## of a non-namespace ends in NULL (equiv. to base environment) [sigh]
    lastEnv <- if(isNamespace(where)) function(where) isBaseNamespace(where) else
    function(where) identical(where, baseenv())
    repeat {
        if(exists(name, envir = where, mode = "function", inherits = FALSE)) {
            f <- get(name, envir = where)
            found <- generic || !is(f, "genericFunction")
        }
        if(found || lastEnv(where))
            break
        where <- parent.env(where)
    }
    if(!found && mustFind)
	stop(if(generic)
	     gettextf("no function %s found", sQuote(name)) else
	     gettextf("no non-generic function %s found", sQuote(name)),
	     domain = NA)
    f
}

el <-
  function(object, where)
  ## element of a vector; numeric index only.
  ##
  ## the definition allows indexing beyond current length of vector
  ## (consistent with [[]] in S but not in R).
  object[where][[1L]]

"el<-" <-
  ## set the element of a vector; numeric index only.
  .Primitive("[[<-")

elNamed <-
  ## get the element of the vector corresponding to name.  No partial matching.
  function(x, name, mustFind=FALSE)
{
    i <- match(name, names(x))
    if(is.na(i)) {
        if(mustFind)
            stop(gettextf("%s is not one of the element names",
                          sQuote(name)),
                 domain = NA)
        else NULL
    }
    else
        el(x,i)
}

"elNamed<-" <-
    ## set the element of the vector corresponding to name.
    function(x, name, value)
{
    x[[name]] <- value
    x
}

formalArgs <-
    ## Returns the names of the formal arguments of this function.
    function(def) names(formals(def))


findFunction <-
    ## return a list of all the places where a function
    ## definition for `name' exists.  If `generic' is FALSE, ignore generic
    ## functions.
    function(f, generic = TRUE, where = topenv(parent.frame()))
{
    allWhere <- .findAll(f, where) # .findAll() in ./ClassExtensions.R
    ok <- logical(length(allWhere))
    for(i in seq_along(allWhere)) {
	wherei <- allWhere[[i]]
	if(exists(f, wherei, inherits = FALSE)) {
	    fdef <- get(f, wherei)
	    ok[i] <- is.function(fdef) && (generic || is.primitive(fdef) || !isGeneric(f, wherei, fdef))
	}## else ok[i] <- FALSE
    }
    allWhere[ok]
}

existsFunction <- function(f, generic=TRUE, where = topenv(parent.frame()))
    length(findFunction(f, generic, where)) > 0L

Quote <- base::quote #was get("quote" , mode = "function")

.message <- function(..., domain = NULL, appendLF = TRUE) {
    ## Output all the arguments, pasted together with no intervening spaces,
    ## wrapping long lines
    text <- paste0(..., collapse="")
    lines <- strwrap(text, width = max(20, 7 * getOption("width") %/% 8))
    message(paste(lines, collapse="\n"), domain = domain, appendLF = appendLF)
}

hasArg <- function(name) {
    aname <- as.character(substitute(name))
    fnames <- names(formals(sys.function(sys.parent())))
    if(is.na(match(aname, fnames))) {
        if(is.na(match("...", fnames)))
            FALSE
        else {
            dotsCall <- eval(quote(substitute(list(...))), sys.parent())
            !is.na(match(aname, names(dotsCall)))
        }
    }
    else
        eval(substitute(!missing(name)), sys.frame(sys.parent()))
}
#  File src/library/methods/R/as.R
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

as <-
  ## Returns the version of this object coerced to be the given `Class'.
  ## If the corresponding `is' relation is true, it will be used.  In particular,
  ## if the relation has a coerce method, the method will be invoked on `object'.
  ##
  ## If the `is' relation is FALSE, and `coerceFlag' is `TRUE',
  ## the coerce function will be called (which will throw an error if there is
  ## no valid way to coerce the two objects).  Otherwise, `NULL' is returned.
  function(object, Class, strict = TRUE, ext = possibleExtends(thisClass, Class))
{
    ## prior to 2.7.0 there was a pseudo-class "double"
    if(.identC(Class, "double")) Class <- "numeric"
    thisClass <- .class1(object)
    if(.identC(thisClass, Class) || .identC(Class, "ANY"))
        return(object)
    where <- .classEnv(thisClass, mustFind = FALSE)
    coerceFun <- getGeneric("coerce", where = where)
    ## get the methods table, use inherited table
    coerceMethods <- .getMethodsTable(coerceFun,environment(coerceFun),inherited= TRUE)
    asMethod <- .quickCoerceSelect(thisClass, Class, coerceFun, coerceMethods, where)
    if(is.null(asMethod)) {
        sig <-  c(from=thisClass, to = Class)
        ## packageSlot(sig) <- where
        ## try first for an explicit (not inherited) method
        ## ?? Can this ever succeed if .quickCoerceSelect failed?
        asMethod <- selectMethod("coerce", sig, optional = TRUE,
                                 useInherited = FALSE, #optional, no inheritance
                                 fdef = coerceFun, mlist = getMethodsForDispatch(coerceFun))
        if(is.null(asMethod)) {
            canCache <- TRUE
            inherited <- FALSE
            if(is(object, Class)) {
                ClassDef <- getClassDef(Class, where)
                ## use the ext information, computed or supplied
                if(identical(ext, FALSE))
                    stop(sprintf("internal problem in as(): %s is(object, \"%s\") is TRUE, but the metadata asserts that the 'is' relation is FALSE",
                                 dQuote(thisClass), Class),
                         domain = NA)
                else if(identical(ext, TRUE))
                    asMethod <- .makeAsMethod(quote(from), TRUE, Class, ClassDef, where)
                else {
                  test <- ext@test
                  asMethod <- .makeAsMethod(ext@coerce, ext@simple, Class, ClassDef, where)
                  canCache <- (!is(test, "function")) || identical(body(test), TRUE)
                 }
            }
            if(is.null(asMethod) && extends(Class, thisClass)) {
                ClassDef <- getClassDef(Class, where)
                asMethod <- .asFromReplace(thisClass, Class, ClassDef, where)
            }
            ## if none of these applies, look for an inherited method
            ## but only on the from argument
            if(is.null(asMethod)) {
                asMethod <- selectMethod("coerce", sig, optional = TRUE,
                                         c(from = TRUE, to = FALSE),
                                         fdef = coerceFun, mlist = coerceMethods)
                inherited <- TRUE
            }
            else if(canCache)  # make into method definition
                asMethod <- .asCoerceMethod(asMethod, thisClass, ClassDef, FALSE, where)
	    if(is.null(asMethod))
		stop(gettextf("no method or default for coercing %s to %s",
			      dQuote(thisClass),
                              dQuote(Class)),
                     domain = NA)
	    else if(canCache) {
		## cache in the coerce function's environment
		cacheMethod("coerce", sig, asMethod, fdef = coerceFun,
			    inherited = inherited)
	    }
        }
    }
    if(strict)
        asMethod(object)
    else
        asMethod(object, strict = FALSE)
}

.quickCoerceSelect <- function(from, to, fdef, methods, where) {
    if(is.null(methods))
        return(NULL)
    else if(is.environment(methods)) {
        method <- .findMethodInTable(c(from, to), methods)
        if(is.environment(method))
            NULL # FIXME:  should resolve by checking package
        else
            method
    }
    else {
        allMethods <- methods@allMethods
        i <- match(from, names(allMethods))
        if(is.na(i))
          NULL
        else {
            methodsi <- allMethods[[i]]
            j <- match(to, names(methodsi))
            if(is.na(j))
              NULL
            else
              methodsi[[j]]
        }
    }
}

.asFromReplace <- function(fromClass, toClass, ClassDef, where) {
    ## toClass extends fromClass, so an asMethod will
    ## be the equivalent of new("toClass", fromObject)
    ## But must check that replacement is defined, in the case
    ## of nonstandard superclass relations
    replaceMethod <- elNamed(ClassDef@contains, fromClass)
    if(is(replaceMethod, "SClassExtension") &&
       !identical(as(replaceMethod@replace, "function"), .ErrorReplace)) {
        f <- function(from, to) NULL
        body(f, envir = where) <-
            substitute({obj <- new(TOCLASS); as(obj, FROMCLASS) <- from; obj},
                       list(FROMCLASS = fromClass, TOCLASS = toClass))
        f
    }
    else
        NULL

}



"as<-" <-
  ## Set the portion of the object defined by the right-hand side.
  ##
  ## Typically, the object being modified extends the class of the right-hand side object,
  ## and contains the slots of that object. These slots (only) will then be replaced.
  function(object, Class, value) {
    thisClass <- .class1(object)
    if(!.identC(.class1(value), Class))
      value <- as(value, Class, strict = FALSE)
    where <- .classEnv(class(object))
    coerceFun <- getGeneric("coerce<-", where = where)
    coerceMethods <- getMethodsForDispatch(coerceFun)
    asMethod <- .quickCoerceSelect(thisClass, Class, coerceFun, coerceMethods, where)
    if(is.null(asMethod)) {
        sig <-  c(from=thisClass, to = Class)
        canCache <- TRUE
        inherited <- FALSE
        asMethod <- selectMethod("coerce<-", sig, TRUE, FALSE, #optional, no inheritance
                                 fdef = coerceFun, mlist = coerceMethods)
        if(is.null(asMethod)) {
            if(is(object, Class)) {
                asMethod <- possibleExtends(thisClass, Class)
                if(identical(asMethod, TRUE)) {# trivial, probably identical classes
                    class(value) <- class(object)
                    return(value)
                }
                else {
                    test <- asMethod@test
                    asMethod <- asMethod@replace
                    canCache <- (!is(test, "function")) || identical(body(test), TRUE)
                    if(canCache) { ##the replace code is a bare function
                        ClassDef <- getClassDef(Class, where)
                        asMethod <- .asCoerceMethod(asMethod, thisClass, ClassDef, TRUE, where)
                    }
                }
            }
            else { # search for inherited method
              asMethod <- selectMethod("coerce<-", sig, TRUE, c(from = TRUE, to = FALSE), doCache = TRUE)
              inherited <- TRUE
            }
        }
        ## cache for next call
        if(canCache && !is.null(asMethod))
                 cacheMethod("coerce<-", sig, asMethod, fdef = coerceFun,
                             inherited = inherited)
     }
    if(is.null(asMethod))
        stop(gettextf("no method or default for as() replacement of %s with Class=\"%s\"",
                      dQuote(thisClass),
                      Class),
             domain = NA)
    asMethod(object, Class, value)
}



setAs <-
  function(from, to, def, replace = NULL, where = topenv(parent.frame()))
{
    ## where there is an "is" relation, modify it
    fromDef <- getClassDef(from, where)
    extds <- possibleExtends(from, to, fromDef)
    if(is(extds, "SClassExtension")) {
        test <- extds@test
        if(is.null(replace))
            replace <- extds@replace
        test <- NULL
        setIs(from, to, test = test, coerce = def, replace = replace, where = where)
    }
    else if(identical(extds, TRUE)) {
        if(.identC(from, to))
            stop(gettextf("trying to set an 'as' relation from %s to itself",
                          dQuote(.class1(from))),
                 domain = NA)
        ## usually to will be a class union, where setAs() is not
        ## allowed by the definition of a union
        toDef <- getClassDef(to, where=where)
        if(is.null(toDef))
            stop(gettextf("class %s is not defined in this environment",
                          dQuote(to)),
                 domain = NA)
        if(isClassUnion(toDef))
            stop(gettextf("class %s is a class union: 'coerce' relations to a class union are not meaningful",
                          dQuote(to)),
                 domain = NA)
        ## else go ahead (but are there any cases here where extds is TRUE?)
        setIs(from, to, coerce = def, replace = replace, where = where)
    }
    ## else extds is FALSE -- no is() action
        args <- formalArgs(def)
        if(!is.na(match("strict", args))) args <- args[-match("strict", args)]
        if(length(args) == 1)
            def <- substituteFunctionArgs(def, "from", functionName = "coerce")
        else  if(length(args) != 2 || !identical(args, c("from", "to")))
               stop(gettextf("'as' method should have one argument, or match the arguments of coerce(): got  (%s)",
                           paste(formalArgs(def), collapse = ", ")),
                  domain = NA)
    ## coerce@.Data is the "prototype" from which we construct the method
        method <- as.list(coerce@.Data) # the function def'n, just to get arguments correct
        method$to <- to
        method <- as.function(method)
        body(method, envir = environment(def)) <- body(def)
        setMethod("coerce", c(from, to), method, where = where)
        if(!is.null(replace)) {
            args <- formalArgs(replace)
            if(identical(args, c("from", "to", "value")))
                method <- replace
            else {
                ## if not from an extends object, process the arguments
                if(length(args) != 2)
                    stop(gettextf("a 'replace' method definition in 'setAs' must be a function of two arguments, got %d", length(args)), domain = NA)
                replace <- body(replace)
                if(!identical(args, c("from", "value"))) {
                    ll <- list(quote(from), quote(value))
                    names(ll) <- args
                    replace <- substituteDirect(replace, ll)
                    warning(gettextf("argument names in 'replace' changed to agree with 'coerce<-' generic:\n%s", paste(deparse(replace), sep="\n    ")),
                            domain = NA)
                }
                method <- eval(function(from, to, value)NULL)
                functionBody(method, envir = .GlobalEnv) <- replace
            }
            setMethod("coerce<-", c(from, to), method, where = where)
        }
}

.setCoerceGeneric <- function(where) {
  ## create the initial version of the coerce function, with methods that convert
  ## arbitrary objects to the basic classes by calling the corresponding as.<Class>
  ## functions.
  setGeneric("coerce", function(from, to, strict = TRUE) {
      if(TRUE) {
          warning("direct use of coerce() is deprecated:  use as(from, class(to)) instead", domain = NA)
          return(as(from, class(to), strict = strict))
      }
      standardGeneric("coerce")
      },
             where = where)
  setGeneric("coerce<-", function(from, to, value) {
      if(TRUE) {
          warning("direct use of coerce() is deprecated:  use as(from, class(to)) <- value instead", domain = NA)
          return(`as<-`(from, class(to), value))
      }
      standardGeneric("coerce<-")
      }, where = where)
  basics <- c(
 "POSIXct",  "POSIXlt", "Date",  "array",  "call",  "character",  "complex",  "data.frame",
 "environment",  "expression",  "factor",  "formula",  "function",  "integer",
 "list",  "logical",  "matrix",  "name",  "numeric",  "ordered",
  "single",  "table",   "vector")
  basics <- basics[!is.na(match(basics,.BasicClasses))]
  for(what in basics) {
      ## if the class is a basic class and there exists an as.<class> function,
      ## use it as the coerce method.
      method  <- .basicCoerceMethod
      switch(what,
	     array =, matrix = body(method, envir = environment(method)) <-
	     substitute({
		 value <- AS(from)
		 if(strict) {
		     dm <- dim(value)
		     dn <- dimnames(value)
		     attributes(value) <- NULL
		     dim(value) <- dm
		     dimnames(value) <- dn
		 }
		 value
	     }, list(AS = as.name(paste0("as.", what)))),
	     ##
	     ts = body(method, envir = environment(method)) <- quote({
		 value <- as.ts(from)
		 if(strict) {
		     attributes(value) <- NULL
		     class(value) <- class(new("ts"))
		     tsp(value) <- tsp(from)
		 }
		 value
	     }),
	     ## default: no attributes
	     body(method, envir = environment(method)) <- substitute({
		 value <- AS(from)
		 if(strict)
		     attributes(value) <- NULL
		 value
	     }, list(AS = as.name(paste0("as.", what))))
	     )
      setMethod("coerce", c("ANY", what), method, where = where)
  }
  ## and some hand-coded ones
  body(method) <- quote(as.null(from))
  setMethod("coerce", c("ANY", "NULL"), method, where = where)
  body(method) <- quote({
            if(length(from) != 1)
                warning("ambiguous object (length != 1) to coerce to \"name\"")
            as.name(from)
        })
  setMethod("coerce", c("ANY","name"), method, where = where)
  ## not accounted for and maybe not needed:  real, pairlist, double
}

.basicCoerceMethod <- function(from, to, strict = TRUE)
    stop("undefined 'coerce' method")

.makeAsMethod <- function(expr, simple, Class, ClassDef, where) {
    if(is(expr, "function")) {
        where <- environment(expr)
        args <- formalArgs(expr)
        if(!identical(args, "from"))
            expr <- .ChangeFormals(expr,
                    if(length(args) > 1) .simpleExtCoerce else .simpleIsCoerce)
        expr <- body(expr)
    }
    ## commented code below is needed if we don't assume asMethod sets the class correctly
#     if(isVirtualClass(ClassDef))
#         value <- expr
#     else if(identical(expr, quote(from)))
#         value <- substitute({class(from) <- CLASS; from},
#                            list(CLASS = Class))
#     else value <- substitute({from <- EXPR; class(from) <- CLASS; from},
#                            list(EXPR = expr, CLASS = Class) )
    ## else
    value <- expr
    if(simple && !identical(expr, quote(from)))
        value <- substitute(if(strict) EXPR else from,
                           list(EXPR = expr))
    f <- .simpleExtCoerce
    body(f, envir = where) <- value
    f
}

## check for and remove a previous coerce method.  Called from setIs
## We warn if the previous method seems to be from a
## setAs(), indicating a conflicting setIs() A previous
## version of setIs is OK, but we remove the methods anyway to be safe.
## Definitions are only removed from the current function's environment,
## not from a permanent copy.
.removePreviousCoerce <- function(from, to, where, prevIs) {
    sig <- c(from, to)
    cdef <- getGeneric("coerce", where = where)
    if(is.null(cdef))
        return(FALSE) # only for booting the methods package?
    prevCoerce <- !is.null(selectMethod("coerce", sig, TRUE, FALSE,
                                        fdef = cdef))
    rdef <- getGeneric("coerce<-", where = where)
    if(is.null(rdef))
        return(FALSE) # only for booting the methods package?
    prevRepl <- !is.null(selectMethod("coerce<-", sig, TRUE, FALSE,
                                      fdef = rdef))
    if(prevCoerce || prevRepl) {
        if(!prevIs)
            warning(gettextf("methods currently exist for coercing from %s to %s; they will be replaced.",
                             dQuote(from),
                             dQuote(to)),
                    domain = NA)
        if(prevCoerce)
            setMethod(cdef, sig, NULL, where = baseenv())
        if(prevRepl)
            setMethod(rdef, sig, NULL, where = baseenv())
        TRUE
    }
    else
        FALSE

}

canCoerce <- function(object, Class) {
    is(object, Class) ||
    !is.null(selectMethod("coerce", c(class(object), Class),
			  optional = TRUE,
			  useInherited = c(from=TRUE, to=FALSE)))
}

## turn raw function into method for coerce() or coerce<-()
## Cheats a little to get past booting the methods package
## (mainly in knowing the slots of the "signature" class).
.asCoerceMethod <- function(def, thisClass, ClassDef, replace, where) {
    fdef <-
	if(replace) quote(function(from, to = TO, value) NULL)
	else	    quote(function(from, to = TO, strict = TRUE) NULL)
    fdef[[2L]]$to <- ClassDef@className
    fdef <- eval(fdef)
    body(fdef, environment(def)) <- body(def)
    attr(fdef, "source") <- deparse(fdef) # because it's wrong from the quote()
    sig <- new("signature")
    sig@.Data <- c(thisClass, ClassDef@className)
    sig@names <- c("from", "to")
    thisPackage <- packageSlot(thisClass)
    sig@package <- if(is.null(thisPackage))
        c(getPackageName(where, FALSE), ClassDef@package) else
        c(thisPackage, ClassDef@package)
    value <- new("MethodDefinition")
    value@.Data <- fdef
    value@target <- sig
    value@defined <- sig
    value@generic <- structure( #FIXME: there should be a genericName()
                if(replace) "coerce<-" else "coerce", package = "methods")
    value
}
#  File src/library/methods/R/cbind.R
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

#### cbind() and rbind()  which build on  cbind2() / rbind2()
#### --------------------------------------------------------
### NOTE: We rely on
### o	dim(.) working reliably for all arguments of [cr]bind2()
### o	All [cr]bind2() methods are assumed to
###	       1) correctly set (row/col)names
###	       2) correctly set (col/row)names for *matrix*(like) arguments

### Note that this
### 1) is namespace-hidden usually,
### 2) cbind / rbind are almost never called in 'methods' itself
### hence, the following has almost no effect unless ``activated'' (see below)

## rbind() is in ./rbind.R  {so it's easier to keep them 100% - synchronized !}

cbind <- function(..., deparse.level = 1)
{
    na <- nargs() - !missing(deparse.level)
    deparse.level <- as.integer(deparse.level)
    stopifnot(0 <= deparse.level, deparse.level <= 2)

    argl <- list(...)
    ## remove trailing 'NULL's:
    while(na > 0 && is.null(argl[[na]])) { argl <- argl[-na]; na <- na - 1 }
    if(na == 0) return(NULL)
    if(na == 1) {
	if(isS4(..1)) return(cbind2(..1))
	else return(.__H__.cbind(..., deparse.level = deparse.level))
    }

    ## else :  na >= 2

    if(deparse.level) {
	symarg <- as.list(sys.call()[-1L])[1L:na] # the unevaluated arguments
	## For default 'deparse.level = 1', cbind(a, b) has to give *names*!
	Nms <- function(i) { # possibly 'deparsed' names of argument  i
	    if(is.null(r <- names(symarg[i])) || r == "") {
		if(is.symbol(r <- symarg[[i]]) || deparse.level == 2)
		    deparse(r)		# else NULL
	    } else r
	}
    }
    if(na == 2) {
	r <- ..2
	fix.na <- FALSE
    }
    else { ## na >= 3 arguments: -- RECURSION -- with care
	## determine nrow(<result>)  for e.g.,	cbind(diag(2), 1, 2)
	## only when the last two argument have *no* dim attribute:
	nrs <- unname(lapply(argl, nrow)) # of length na
	iV <- sapply(nrs, is.null)# is 'vector'
	fix.na <- identical(nrs[(na-1):na], list(NULL,NULL))
	if(fix.na) {
	    ## "fix" last argument, using 1-column `matrix' of proper nrow():
	    nr <- max(if(all(iV)) sapply(argl, length) else unlist(nrs[!iV]))
	    argl[[na]] <- cbind(rep(argl[[na]], length.out = nr),
				deparse.level = 0)
	    ## and since it's a 'matrix' now, cbind() below may not name it
	}
	## need to pass argl, the evaluated arg list to do.call();
	## OTOH, these may have lost their original 'symbols'
	if(deparse.level) {
	    if(fix.na)
		fix.na <- !is.null(Nna <- Nms(na))
	    if(!is.null(nmi <- names(argl))) iV <- iV & (nmi == "")
	    ## attach `symbols' to argl[-1L] for 'vectors'[iV]
	    ii <- if(fix.na) # need to fix later ([na] is 'matrix')
		2:(na-1) else 2:na
	    if(any(iV[ii])) {
		for(i in ii[iV[ii]])
		    if (!is.null(nmi <- Nms(i))) names(argl)[i] <- nmi
	    }
	}
	r <- do.call(cbind, c(argl[-1L], list(deparse.level=deparse.level)))
    }

    d2 <- dim(r)
    r <- cbind2(..1, r)
    if(deparse.level == 0)
	return(r)
    ism1 <- !is.null(d1 <- dim(..1)) && length(d1) == 2L
    ism2 <- !is.null(d2)	     && length(d2) == 2L && !fix.na
    if(ism1 && ism2) ## two matrices
	return(r)

    ## else -- Setting colnames correctly
    ##	       when one was not a matrix [needs some diligence!]
    Ncol <- function(x) {
	d <- dim(x); if(length(d) == 2L) d[2L] else as.integer(length(x) > 0L) }
    nn1 <- !is.null(N1 <- if((l1 <- Ncol(..1)) && !ism1) Nms(1)) # else NULL
    nn2 <- !is.null(N2 <- if(na == 2 && Ncol(..2) && !ism2) Nms(2))
    if(nn1 || nn2 || fix.na) {
	if(is.null(colnames(r)))
	    colnames(r) <- rep.int("", ncol(r))
	setN <- function(i, nams)
	    colnames(r)[i] <<- if(is.null(nams)) "" else nams
	if(nn1) setN(1,	 N1)
	if(nn2) setN(1+l1, N2)
	if(fix.na) setN(ncol(r), Nna)
    }
    r
}

## To be active, the above cbind() must "replace" cbind() in "base".
## This may be called on loading methods, see ./zzz.R
bind_activation <- function(on = TRUE)
{
    inBase <- function(x, value, ns)
    {
        unlockBinding(x, ns)
        assign(x, value, envir = ns, inherits = FALSE)
        w <- options("warn")
        on.exit(options(w))
        options(warn = -1)
        lockBinding(x, ns)
        invisible(NULL)
    }

    ## 'bind' : cbind && rbind
    ## as from 2.4.0 this saving is done in base, so could simplify code
    base.ns <- getNamespace("base")
    saved <- exists(".__H__.cbind", envir = base.ns, inherits = FALSE)
    if(was.on <- saved)
        was.on <- !identical(base::cbind, base::.__H__.cbind)
    if(on) {
        if(!saved) {
            inBase(".__H__.cbind", base::cbind, base.ns)
            inBase(".__H__.rbind", base::rbind, base.ns)
        }
	inBase("cbind", cbind, base.ns)
	inBase("rbind", rbind, base.ns)
    }
    else if(!on && was.on) { ## turn it off
        inBase("cbind", base::.__H__.cbind, ns = base.ns)
        inBase("rbind", base::.__H__.rbind, ns = base.ns)
    }
    was.on
}

### cbind2 () :	 Generic and methods need to be "method-bootstrapped"
### --------   --> ./MethodsListClass.R
#  File src/library/methods/R/fixPrevious.R
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

## fixPre1.8(names)
##   The objects in names should have been loaded from a version of R
##   previous to 1.8.0
##   The classes of these objects must be defined in the current session.
##   Objects are modified to have the correct version of its class,
##   and re-assigned.
##
##   The function checks for existence and consistency of class definitions.
fixPre1.8 <- function(names, where = topenv(parent.frame())) {
    done <- character()
    for(what in names) {
        objWhere <- methods:::.findAll(what, where)
        if(length(objWhere) == 0) {
            warning(gettextf("object %s not found",
                             sQuote(what)),
                    domain = NA)
            next
        }
        objWhere <- objWhere[[1L]]
        obj <- get(what, objWhere)
        ## don't fix up basic datatypes with no explicit class
        if(is.null(attr(obj, "class")))
            next
        Class <- class(obj)
        if(is.null(attr(Class, "package"))) {
            if(isClass(Class, where = where)) {
                ClassDef <- getClass(Class, where = where)
                ok <- !(isVirtualClass(ClassDef) ||
			!isTRUE(validObject(obj, test=TRUE)))
                if(ok) {
                    class(obj) <- ClassDef@className
                    assign(what, obj, objWhere)
                    done <- c(done, what)
                }
                else
                    warning(gettextf("object %s not changed (it is not consistent with the current definition of class %s from %s)",
                                     sQuote(what),
                                     dQuote(Class),
                                     sQuote(ClassDef@package)),
                            domain = NA)
            }
            else
                warning(gettextf("no definition for the class of %s (class %s) found",
                                 sQuote(what),
                                 dQuote(class)),
                        domain = NA)
        }
        else
            warning(gettextf("object %s not changed (it does not appear to be from a version of R earlier than 1.8.0)",
                             sQuote(what)),
                    domain = NA)
    }
    done
}


#  File src/library/methods/R/is.R
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

is <-
  # With two arguments, tests whether `object' can be treated as from `class2'.
  #
  # With one argument, returns all the super-classes of this object's class.
function(object, class2)
{
    cl <- class(object)
    S3Case <- length(cl) > 1L
    if(S3Case)
      cl <- cl[[1L]]
    if(missing(class2))
        return(extends(cl))
    class1Def <- getClassDef(cl)
    if(is.null(class1Def)) # an unregistered S3 class
      return(inherits(object, class2))
    if(is.character(class2))
      class2Def <- getClassDef(class2, .classDefEnv(class1Def))
    else {
        class2Def <- class2
        class2 <- class2Def@ className
    }
    ## S3 inheritance is applied if the object is not S4 and class2 is either a basic
    ## class or an S3 class (registered or not)
    S3Case <- S3Case || (is.object(object) && !isS4(object)) # first requirement
    S3Case <- S3Case && (is.null(class2Def) || class2 %in% .BasicClasses ||
                         extends(class2Def, "oldClass"))
    if(S3Case)
        return(inherits(object, class2))
    if(.identC(cl, class2) || .identC(class2, "ANY"))
        return(TRUE)
    ext <- possibleExtends(cl, class2, class1Def, class2Def)
    if(is.logical(ext))
        ext
    else if(ext@simple)
        TRUE
    else
       ext@test(object)
}

extends <-
  ## Does the first class extend the second class?
  ## Returns `maybe' if the extension includes a non-trivial test.
  function(class1, class2, maybe = TRUE, fullInfo = FALSE)
{
    if(is.character(class1)) {
        if(length(class1) > 1L)
            class1 <- class1[[1L]]
	classDef1 <- getClassDef(class1)
    } else if(is(class1, "classRepresentation")) {
	classDef1 <- class1
	class1 <- classDef1@className
    }
    else
	stop("'class1' must be the name of a class or a class definition")
    if(missing(class2)) {
        if(is.null(classDef1))
            return(class1)
        ext <- classDef1@contains
        if(!identical(maybe, TRUE) && length(ext) > 0)
        {
            noTest <- sapply(ext, function(obj)identical(body(obj@test), TRUE))
            ext <- ext[noTest]
        }
        if(fullInfo) {
            elNamed(ext, class1) <- TRUE
            return(ext)
        }
        else
            return(c(class1,names(ext)))
    }
    value <- NULL
    if(is.character(class2) && length(class2) == 1L) { ## fast first checks
	## the [[1L]] below handles old-style classes & throws away package attributes
	if(.identC(class1[[1L]], class2) || .identC(class2, "ANY"))
          return(TRUE)
        if(!is.null(classDef1) && class2 %in% names(classDef1@contains))
	    value <- classDef1@contains[[class2]]
        else
          classDef2 <- getClassDef(class2)
    }
    else if(is(class2, "classRepresentation")) {
	classDef2 <- class2
	class2 <- class2@className
    }
    else
	stop("'class2' must be the name of a class or a class definition")
    if(is.null(value))
      value <- possibleExtends(class1, class2, classDef1, classDef2)
    if(fullInfo)
        value
    else if(is.logical(value))
        value
    else if(value@simple || identical(body(value@test), TRUE))
        TRUE
    else
        maybe
}

.specialVirtual <- c("oldClass")

setIs <-
  ## Defines class1 to be an extension of class2.
  ## The relationship can be conditional, if a function is supplied as the `test'
  ## argument.  If a function is supplied as the `coerce' argument, this function will
  ## be applied to any `class1' object in order to turn it into a `class2' object.
  ##
  ## Extension may imply that a `class1' object contains a `class2' object.  The default
  ## sense of containing is that all the slots of the simpler class are found in the
  ## more elaborate one.  If the `replace' argument is supplied as an S replacement
  ## function, this function will be used to implement `as(obj, class2) <- value'.
  function(class1, class2, test = NULL, coerce = NULL,
           replace = NULL, by = character(), where = topenv(parent.frame()),
           classDef = getClass(class1, TRUE, where = where), extensionObject = NULL, doComplete = TRUE)
{
    ## class2 should exist
    where <- as.environment(where)
    classDef2 <- getClassDef(class2, where)
    if(is.null(classDef2))
        stop(gettextf("class %s has no visible definition from package or environment %s",
                      dQuote(class2),
                      sQuote(getPackageName(where))),
             domain = NA)
    ## check some requirements:
    ## One of the classes must be on the target environment (so that the relation can
    ## be retained by saving the corresponding image)
    m1 <- classMetaName(class1)
    local1 <- exists(m1, where, inherits = FALSE) &&
    !(classDef@sealed || bindingIsLocked(m1, where))
    m2 <- classMetaName(class2)
    local2 <- exists(m2, where, inherits = FALSE) &&
    !(classDef2@sealed || bindingIsLocked(m2, where))
    if(!(local1 || local2) )
        stop(gettextf("cannot create a 'setIs' relation when neither of the classes (%s and %s) is local and modifiable in this package",
                      dQuote(class1),
                      dQuote(class2)),
             domain = NA)
    if(classDef@sealed && !isClassUnion(classDef2))
        stop(gettextf("class %s is sealed; new superclasses can not be defined, except by 'setClassUnion'",
                      dQuote(class1)),
             domain = NA)
    prevIs <- !identical(possibleExtends(class1, class2,classDef, classDef2),
                         FALSE) # used in checking for previous coerce
    if(is.null(extensionObject))
        obj <- makeExtends(class1, class2, coerce, test, replace, by,
                           classDef1 = classDef, classDef2 = classDef2,
                           package = getPackageName(where))
    else
        obj <- extensionObject
    ## revise the superclass/subclass info in the stored class definition
    ok <- .validExtends(class1, class2, classDef,  classDef2, obj@simple)
    if(!identical(ok, TRUE))
      stop(ok)
    where2 <- .findOrCopyClass(class2, classDef2, where, "subclass")
    elNamed(classDef2@subclasses, class1) <- obj
    if(doComplete)
        classDef2@subclasses <- completeSubclasses(classDef2, class1, obj, where)
    ## try to provide a valid prototype for virtual classes
    if(classDef2@virtual && is.na(match(class2, .specialVirtual))) {
        ## For simplicity, we prefer NULL prototype if "NULL"
        ## is a subclass of a virtual class; otherwise the
        ## prototype is an element of class1 or its prototype if VIRTUAL
        if(extends(classDef, "NULL"))
            classDef2@prototype <- NULL
        else if(is.null(classDef2@prototype)
                && is.na(match("NULL", names(classDef2@subclasses)))) {
            if(classDef@virtual)
                classDef2@prototype <- classDef@prototype
            else # new(), but without intialize(), which may require an arg.
                classDef2@prototype <- .Call(C_new_object, classDef)
        }
    }
    assignClassDef(class2, classDef2, where2, TRUE)
    .removePreviousCoerce(class1, class2, where, prevIs)
    where1 <- .findOrCopyClass(class1, classDef, where, "superClass")
    ## insert the direct contains information in a valid spot
    .newDirectSuperclass(classDef@contains, class2, names(classDef2@contains)) <- obj
    if(doComplete) {
      classDef@contains <- completeExtends(classDef, class2, obj, where = where)
      if(!is(classDef, "ClassUnionRepresentation")) #unions are handled in assignClassDef
        .checkSubclasses(class1, classDef, class2, classDef2, where1, where2)
    }
    assignClassDef(class1, classDef, where1, TRUE)
    invisible(classDef)
 }

.findOrCopyClass <- function(class, classDef, where, purpose) {
    whereIs <- findClass(classDef, where)
    if(length(whereIs))
      whereIs[[1L]]
    else {
        warning(gettextf("class %s is defined (with package slot %s) but no metadata object found to revise %s information---not exported?  Making a copy in package %s",
                         .dQ(class), sQuote(classDef@package), purpose,
                         sQuote(getPackageName(where, FALSE))),
                call. = FALSE, domain = NA)
        where
    }
}


.validExtends <- function(class1, class2, classDef1,  classDef2, slotTests) {
    .msg <- function(class1, class2)
        gettextf("class %s cannot extend class %s",
                 dQuote(class1),
                 dQuote(class2))
    if((is.null(classDef1) || is.null(classDef2)) &&
       !(isVirtualClass(class1) && isVirtualClass(class2)))
        return(c(.msg(class1, class2), ": ",
             gettext("both classes must be defined")))
    if(slotTests) {
        slots2 <- classDef2@slots
        if(length(slots2)) {
            n2 <- names(slots2)
            slots1 <- classDef1@slots
            n1 <- names(slots1)
            if(any(is.na(match(n2, n1))))
                return(c(.msg(class1, class2), ": ",
                         sprintf(ngettext(sum(is.na(match(n2, n1))),
                                          "class %s is missing slot from class %s (%s), and no coerce method was supplied",
                                          "class %s is missing slots from class %s (%s), and no coerce method was supplied"),
                                 dQuote(class1),
                                 dQuote(class2),
                                 paste(n2[is.na(match(n2, n1))], collapse = ", "))))
            bad <- character()
            for(what in n2)
                if(!extends(elNamed(slots1, what), elNamed(slots2, what)))
                    bad <- c(bad, what)
            if(length(bad))
                return(c(.msg(class1, class2), ": ",
                         sprintf(ngettext(length(bad),
                                          "slot in class %s must extend corresponding slot in class %s: fails for %s",
                                          "slots in class %s must extend corresponding slots in class %s: fails for %s"),
                                 dQuote(class1),
                                 dQuote(class2),
                                 paste(bad, collapse = ", "))))
        }
    }
    TRUE
}

".newDirectSuperclass<-" <- function(contains, class2, superclasses2, value) {
    superclasses <- names(contains)
    if(length(superclasses2) == 0 || length(superclasses) == 0 ||
       all(is.na(match(superclasses2, superclasses))))
      elNamed(contains, class2) <- value
    else {
        sq <- seq_along(superclasses)
        before <- (sq[match(superclasses, superclasses2,0)>0])[[1]]
        contains <- c(contains[sq < before], value, contains[sq >= before])
        superclasses <- c(superclasses[sq < before], class2, superclasses[sq >= before])
        names(contains) <- superclasses
    }
    contains
}

#  File src/library/methods/R/languageEl.R
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

languageEl <-
  ## extract an element of a language object, consistently
  ## for different kinds of objects.
  ##
  ## The 1st., etc. elements of a function are the corresponding
  ## formal arguments, with the default expression if any as value.
  ##
  ## The first element of a call is the name or the function object being
  ## called.  The 2nd, 3rd, etc. elements are the 1st, 2nd, etc. arguments expressions.
  ## Note that the form of the extracted name is different for R and S-Plus.
  ## When the name (the first element) of a call is replaced, the languageEl replacement
  ## function coerces a character string to the internal form for each system.
  ##
  ## The 1st, 2nd, 3rd elements of an `if' expression are the test, first, and second branch.
  ##
  ## The 1st element of a `for' object is the name (symbol) being used in the loop,
  ## the second is the expression for the range of the loop, the third is the body of the loop.
  ##
  ## The first element of a `while' object is the loop test, and the second the body of
  ## the loop.
  function(object, which)
{
    data <- as.list(object)
    if(is.character(which))
        data[[which]]
    else if(typeof(object) == "language") {
        if(isGrammarSymbol(data[[1L]]))
            data[[which + 1]]
        else
            data[[which]]               ## other calls
    }
    else data[[which]]
}

"languageEl<-" <-
  ## replace an element of a language object, see "languageEl" for meaning.
  function(object, which, value)
{
    data <- as.list(object)
    n <- length(data)
    type <- typeof(object)
    if(type == "closure") {
        ev <- environment(object)
        if(is.character(which)) {
            if(is.na(match(which, names(data)))) {
                body <- data[[n]]
                data <- data[-n]
                data[[which]] <- value
                data[[n+1]] <- body
            }
            else
                data[[which]] <- value
        }
        else {
            if(which < 1 || which > n)
                stop("invalid index for function argument")
            ## we don't warn if this is used to replace the body (which == n)
            ## but maybe we should.
            data[[which]] <- value
        }
        object <- as.function(data)
        environment(object) <- ev
        object
    }
    else if(type == "language") {
        if(is.character(which))
            data[[which]] <- value
        else if(isGrammarSymbol(data[[1L]]))
            data[[which+1]] <- value
        else {
            if(identical(which, 1) && is.character(value))
                value <- as.symbol(value)
            data[[which]] <- value
        }
        as.call(data)
    }
    else {
        object[[which]] <- value
        object
    }
}


isGrammarSymbol <-
  function(symbol)
{
    if(typeof(symbol) != "symbol")
        FALSE
    else  switch(as.character(symbol),
                 ## the grammatical constructions
                 "{" =, "if" = , "for"= ,
                 "while" = , "repeat" = ,
                 "return" = , "next" = ,
                 "break" = , "<-" = , "<<-" = TRUE,
                 FALSE)
}
#  File src/library/methods/R/makeBasicFunsList.R
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

## the executable code to complete the generics corresponding to primitives,
## and to define the group generics for these functions.

## uses the primitive list and the function .addBasicGeneric
## defined (earlier) in BasicFunsList.R

utils::globalVariables(".addBasicGeneric")

.makeBasicFuns <- function(where)
{
    funs <- get(".BasicFunsList", envir=where)

    ## First, set up the existing functions in the list as valid generics.
    ## This will override anything except the S4 group generics.
    curNames <- names(funs)
    for(i in seq_along(funs)) {
	val <- funs[[i]]
        if (is.function(val))
            funs <- .addBasicGeneric(funs, curNames[[i]], val, "")
    }

    ## Next, add the remaining primitive generics
    prims <- ls(.GenericArgsEnv, all.names=TRUE)
    new_prims <- prims[!prims %in% names(funs)]
    for(nm in new_prims) {
        f <- get(nm, envir = .GenericArgsEnv)
        body(f) <- substitute(standardGeneric(ff), list(ff=val))
        funs <- .addBasicGeneric(funs, nm, f, "")
    }

    ## Then add all the primitives that are not already there.
    ff <- ls("package:base", all.names=TRUE)
    prims <- ff[sapply(ff, function(x) is.primitive(get(x, "package:base")))]
    new_prims <- prims[!prims %in% names(funs)]
    add <- rep(list(FALSE), length(new_prims))
    names(add) <- new_prims
    funs <- c(funs, add)

    ## the Math group.
    members <- c("abs", "sign", "sqrt",
		 "ceiling", "floor", "trunc",
		 "cummax", "cummin", "cumprod", "cumsum",
		 "exp", "expm1",
		 "log", "log10", "log2", "log1p",
		 "cos", "cosh", "sin", "sinh", "tan", "tanh",
		 "acos", "acosh", "asin", "asinh", "atan", "atanh",
		 "gamma", "lgamma", "digamma", "trigamma"
		 )
    for(f in members) {
	funs <-
	    .addBasicGeneric(funs, f,
			     if(f %in% c("log", "trunc")) {
				 function(x, ...) standardGeneric("")
			     } else   function(x) standardGeneric(""),
			     "Math")
    }

    setGroupGeneric(where=where, "Math", function(x)NULL,
		    knownMembers = members, package = "base")

    ## The Math2 group.
    funs <- .addBasicGeneric(funs, "round",
			     function(x, digits = 0) standardGeneric(""),
			     "Math2")
    funs <- .addBasicGeneric(funs, "signif",
			     function(x, digits = 6) standardGeneric(""),
			     "Math2")

    setGroupGeneric(where = where, "Math2", function(x, digits) NULL,
		    knownMembers = c("round", "signif"), package = "methods")

    ## The Arith group
    members <- c("+", "-", "*", "^", "%%", "%/%", "/")
    for(f in members)
	funs <- .addBasicGeneric(funs, f, function(e1, e2) standardGeneric(""),
				 "Arith")

    setGroupGeneric(where = where, "Arith", function(e1, e2)NULL,
		    group = "Ops", knownMembers = members, package = "base")

    ## the Compare group
    members <- c("==", ">", "<", "!=", "<=", ">=")
    for(f in members)
	funs <- .addBasicGeneric(funs, f, function(e1, e2) standardGeneric(""),
				 "Compare")

    setGroupGeneric(where = where, "Compare", function(e1, e2)NULL,
		    group = "Ops", knownMembers = members, package = "methods")

    ## The Logic group
    members <- c("&", "|") ## *not*  "!" since that has only one argument
    for(f in members)
	funs <- .addBasicGeneric(funs, f, function(e1, e2) standardGeneric(""),
				 "Logic")
    setGroupGeneric(where = where, "Logic", function(e1, e2) NULL,
		    group = "Ops", knownMembers = members, package = "base")

    ## the Ops group generic

    setGroupGeneric(where = where,"Ops", function(e1, e2) NULL,
		    knownMembers = c("Arith", "Compare", "Logic"),
                    package = "base")


    ## The Summary group

    ## These are a bit problematic, since they essentially have "..."
    ## as their only data-related formal argument.  The treatment
    ## since S3 has been to define the generic with a special first
    ## argument, to allow method dispatch.  But the method had better
    ## invoke the generic recursively or perform some other special
    ## computations, in order to avoid unintended anomalies, such as
    ## !identical(max(x,y), max(y,x))

    members <- c("max", "min", "range", "prod", "sum", "any", "all")
    for(f in members)
	funs <- .addBasicGeneric(funs, f, function (x, ..., na.rm = FALSE)
				 standardGeneric(""),
				 "Summary")

    setGroupGeneric(where = where, "Summary",
		    function(x, ..., na.rm = FALSE) NULL,
		    knownMembers = members, package = "base")

    ## The Complex group

    ## R adds this group to the previous S language function groups,
    ## for all the operations defined on complex numbers.  For
    ## applications wanting to define a new class that extends the
    ## concept of complex numbers, a function group is likely to be
    ## useful since all these functions may operate in a similar
    ## manner (by analogy, e.g., with the Math group).

    members <- c("Arg", "Conj", "Im", "Mod", "Re")
    for(f in members)
	funs <- .addBasicGeneric(funs, f, function(z) standardGeneric(""),
				 "Complex")

    setGroupGeneric(where=where,"Complex", function(z)NULL,
		    knownMembers = members, package = "base")

    assign(".BasicFunsList", funs, envir=where)
    rm(.addBasicGeneric, envir=where)
}


.initImplicitGenerics <- function(where)
{
    ## create implicit generics & possibly  methods for the functions in .BasicFunsList.

    setGeneric("with", signature = "data", where = where)
    setGenericImplicit("with", where, FALSE)

    ## when setMethod()ing on chol2inv, one should *not* have to deal with
    ## arguments  'size' and 'LINPACK' :
    setGeneric("chol2inv", function(x, ...) standardGeneric("chol2inv"),
	       useAsDefault = function(x, ...) base::chol2inv(x, ...),
	       signature = "x", where = where)
    setGenericImplicit("chol2inv", where, FALSE)

    setGeneric("rcond", function(x, norm, ...) standardGeneric("rcond"),
	       useAsDefault = function(x, norm, ...) base::rcond(x, norm, ...),
	       signature = c("x", "norm"), where = where)
    setGenericImplicit("rcond", where, FALSE)

    setGeneric("norm", function(x, type, ...) standardGeneric("norm"),
	       useAsDefault = function(x, type, ...) base::norm(x, type, ...),
	       signature = c("x", "type"), where = where)
    setGenericImplicit("norm", where, FALSE)

    setGeneric("backsolve", function(r, x, k = ncol(r), upper.tri = TRUE, transpose = FALSE, ...)
	       standardGeneric("backsolve"),
	       useAsDefault =
	       function(r, x, k = ncol(r), upper.tri = TRUE, transpose = FALSE, ...)
	       base::backsolve(r, x, k = k,
			       upper.tri = upper.tri, transpose = transpose, ...),
	       signature = c("r", "x"), where = where)
    setGenericImplicit("backsolve", where, FALSE)

    setGeneric("colMeans", function(x, na.rm = FALSE, dims = 1, ...)
			standardGeneric("colMeans"),
	       useAsDefault = function(x, na.rm = FALSE, dims = 1, ...)
			base::colMeans(x, na.rm=na.rm, dims=dims, ...),
	       signature = c("x", "na.rm", "dims"), where = where)
    setGeneric("colSums", function(x, na.rm = FALSE, dims = 1, ...)
			standardGeneric("colSums"),
	       useAsDefault = function(x, na.rm = FALSE, dims = 1, ...)
			base::colSums(x, na.rm=na.rm, dims=dims, ...),
	       signature = c("x", "na.rm", "dims"), where = where)
    setGeneric("rowMeans", function(x, na.rm = FALSE, dims = 1, ...)
			standardGeneric("rowMeans"),
	       useAsDefault = function(x, na.rm = FALSE, dims = 1, ...)
			base::rowMeans(x, na.rm=na.rm, dims=dims, ...),
	       signature = c("x", "na.rm", "dims"), where = where)
    setGeneric("rowSums", function(x, na.rm = FALSE, dims = 1, ...)
			standardGeneric("rowSums"),
	       useAsDefault = function(x, na.rm = FALSE, dims = 1, ...)
			base::rowSums(x, na.rm=na.rm, dims=dims, ...),
	       signature = c("x", "na.rm", "dims"), where = where)
    setGenericImplicit("colMeans", where, FALSE)
    setGenericImplicit("colSums",  where, FALSE)
    setGenericImplicit("rowMeans", where, FALSE)
    setGenericImplicit("rowSums",  where, FALSE)

    setGeneric("sample", function(x, size, replace = FALSE, prob = NULL, ...)
			standardGeneric("sample"),
	       useAsDefault = function(x, size, replace = FALSE, prob = NULL, ...)
			base::sample(x, size, replace=replace, prob=prob, ...),
	       signature = c("x", "size"), where = where)
    setGenericImplicit("sample", where, FALSE)

    ## not implicitGeneric() which is not yet available "here"
    registerImplicitGenerics(where = where)
}
#  File src/library/methods/R/method.skeleton.R
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

method.skeleton <- function (generic, signature, file, external = FALSE,
			     where = topenv(parent.frame()))
{
    fdef <- getGeneric(generic, where = where)
    if (is.null(fdef)) {
        fdef <- implicitGeneric(generic, where = where)
        if(is.null(fdef))
            stop(gettextf("no function definition found for %s",
                          sQuote(generic)),
                 domain = NA)
    }
    else {
        generic <- fdef@generic
    }
    signature <- matchSignature(signature, fdef)
    if (length(signature) == 0)
        signature <- "ANY"
    sigNames <- fdef@signature
    length(sigNames) <- length (signature)
    method <- function() {
    }
    formals(method) <- formals(fdef)
    body(method) <- quote({
        stop("need a definition for the method here")
    })
    methodName <- paste(c(generic, signature), collapse = "_")
    if (missing(file))
        file <- paste0(methodName, ".R")
    output <- c(paste0("setMethod(\"", generic, '",'),
		paste0("    signature(", paste0(sigNames, ' = "', signature, '"',
						collapse = ", "), "),"))
    method <- deparse(method)
    if (identical(external, FALSE))
        output <- c(output, paste0("    ", method), ")")
    else {
        if(is(external, "character") )
            methodName <- toString(external)
        method[[1L]] <- paste0("`", methodName, "` <- ", method[[1L]])
        output <- c(method, "", output, paste0("  `", methodName, "`)"))
    }
    writeLines(output, file)
    message(gettextf("Skeleton of method written to %s",
                     if (is.character(file)) file else "connection"),
            domain = NA)
    invisible(file)
}
#  File src/library/methods/R/methods-deprecated.R
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


## <entry>
## Deprecated in 2.10.0
## Defunct in 2.11.0
## Removed in 3.0.0
## trySilent <- function(expr) .Defunct("try(silent = TRUE)")
## </entry>

## <entry>
## Defunct in 3.0.0
traceOn <- function(what, tracer = browseAll, exit = NULL) {
    browseAll <- function() .Defunct()
    .Defunct("trace")
}
traceOff <- function(whatL) .Defunct("untrace")
## </entry>
#  File src/library/methods/R/methods-deprecated.R
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


## <entry>
## Deprecated in 2.10.0
## trySilent <- function(expr) {
##     .Deprecated("try(*, silent=TRUE)  or {more efficiently}\n tryCatch(*, error=function(e) e)")
##     try(expr, silent = TRUE)
## }
## </entry>
#  File src/library/methods/R/methodsTable.R
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

### merge version called from namespace imports code.  Hope to avoid using generic
.mergeMethodsTable2 <- function(table, newtable, envir, metaname) {
    old <- objects(envir=table, all.names=TRUE)
    mm <- 1
    for( what in old) {
      mm <- get(what, envir =table)
      if(is(mm, "MethodDefinition")) {
          mm <- length(mm@defined)
          break
      }
    }
    new <- objects(envir=newtable, all.names=TRUE)
    ## check that signature length doesn't change
    canStore <- TRUE
    for(what in new) {
        obj <- get(what, envir = newtable)
        if(is(obj, "MethodDefinition") &&
           length(obj@defined) != mm) {
            canStore <- FALSE
            break
        }
    }
    if(canStore) {
        for(what in new)
          assign(what, get(what, envir = newtable), envir = table)
        table
    }
    else { # rats! have to get the generic function
        f <- gsub(".__T__(.*):([^:]+)", "\\1", metaname)
        package <- gsub(".__T__(.*):([^:]+(.*))", "\\2", metaname)
        generic <- getGeneric(f, TRUE, envir, package)
        .mergeMethodsTable(generic, table, newtable, TRUE)
        table
    }
}

## action on attach, detach to merge methods tables
.mergeMethodsTable <- function(generic, table, newtable, add = TRUE) {
  fenv <- environment(generic)
  signature <- generic@signature
  if(!exists(".SigLength", envir = fenv, inherits = FALSE))
     .setupMethodsTables(generic)
  if(add)
      allTable <- NULL # .AllMTable but only if required
  else
      allTable <- get(".AllMTable", envir = fenv)
  n <- get(".SigLength", envir = fenv)
  anySig <- rep("ANY", n) # assert doesn't need to be a real signature
  anyLabel <- .sigLabel(anySig)
  newMethods <- objects(envir=newtable, all.names=TRUE)
  for(what in newMethods) {
    obj <- get(what, envir = newtable)
    if(is.primitive(obj))
      sig <- anySig
    else if(is(obj, "MethodDefinition"))
      sig <- obj@defined
    else if(is.environment(obj)) {
       objsWhat <- objects(obj, all.names=TRUE)
       if(length(objsWhat) == 0)
           next # empty environment, ignore
       sig <- NULL
       for(ww in objsWhat) {
           objw <- get(ww, envir = obj)
           if(is(objw, "MethodDefinition"))
               sig <- objw@defined
       }
       if(is.null(sig))
           sig <- anySig
    }
    else
      stop(gettextf("invalid object in meta table of methods for %s, label %s, had class %s",
                    sQuote(generic@generic),
                    sQuote(what),
                    dQuote(class(obj))),
           domain = NA)
    ns <- length(sig)
    if(ns == n) {}
    else {
      if(ns < n) {
        nadd <- n - ns
        sigPackage <- packageSlot(sig)
        if(length(sigPackage)< ns) # probably out of date?
            sigPackage <- c(sigPackage, rep("", ns - length(sigPackage)))
        sig <-  .simpleSignature(c(sig, rep("ANY", nadd)),
                    names = generic@signature[1:n],
                    packages = c(sigPackage, rep("methods", nadd)))
        obj <- .xpdSignature(obj, sig, n-ns)
        what <- .sigLabel(sig)
        ns <- n
      }
      else if(add) { # and ns > n
        signames <- generic@signature
        length(signames) <- ns
        .resetTable(table, ns, signames)
        assign(".SigLength", ns, envir = fenv)
        n <- ns
      }
    }
    if(add) {
        if(exists(what, envir = table, inherits = FALSE)) {
            obj <- .newOrMultipleMethod(obj, what, table)
            ## must replace in .AllMTable also
            if(is.null(allTable))
                allTable <- get(".AllMTable", envir = fenv)
            assign(what, obj, envir = allTable)
        }
        assign(what, obj, envir = table)
    }
    else if(exists(what, envir = table, inherits = FALSE) &&
            !all(obj@defined == "ANY") ) {
        ## remove methods, but not the default
        remove(list = what, envir = allTable)
        remove(list = what, envir = table)
    }
    ## else warning?
  }
  NULL
}

.xpdSignature <- function(obj, sig, nadd) {
    if(is(obj, "MethodDefinition")) {
        obj@defined <- sig
        obj@target <- sig
    }
    else if(is.environment(obj)) {
        xtrPkg <- rep("methods", nadd)
        for(what in objects(obj)) {
            objw <- get(what, envir = obj)
            if(is(objw, "MethodDefinition")) {
                sigw <- objw@defined
                pkgw <- packageSlot(sigw)
                if(length(pkgw) < length(sigw))
                    pkgw <- c(pkgw, rep("", length(sigw) - length(pkgw)))
                sigw <- .simpleSignature( c(sigw, rep("ANY", nadd)),
                        names = names(sig),
                        packages = c(pkgw, rep("methods", nadd)))
                objw@defined <- objw@target <- sigw
                remove(list = what, envir = obj)
                var <- .pkgMethodLabel(objw)
                if(nzchar(var)) assign(var, objw, envir = obj)
            }
        }
    }
    obj
}

## a simpler version of setting up a signature object
## For better or worse, the initialize() method expects
## a function definition and calls .MakeSignature()
.simpleSignature <- function(classes, names, packages) {
    object <- new("signature")
    object@.Data <- classes
    object@names <- names
    object@package <- packages
    object
}

.newOrMultipleMethod <- function(obj, what, table) {
    if(!.duplicateClassesExist())
        return(obj)
    current <- get(what, envir = table)
    if(is.environment(current)) {
        if(is.environment(obj))
            for(whatObj in objects(obj, all.names = TRUE))
                assign(whatObj, get(whatObj, envir = obj),
                       envir = current)
        else if(is(obj, "MethdodDefinition")) {
            var <- .pkgMethodLabel(obj)
            if(nzchar(var)) assign(var, obj, envir = current)
        }
        current
    }
    else if(is(current, "MethodDefinition")) {
        curPkg <- packageSlot(current@defined)
        if(is(obj, "MethodDefinition")) {
            objPkg <- packageSlot(obj@defined)
            if(is.null(curPkg) || is.null(objPkg) ||
               identical(curPkg, objPkg))
                return(obj)
            else {
                merge <- new.env()
                var <- .pkgMethodLabel(obj)
                if(nzchar(var)) assign(var, obj, envir = merge)
                var <- .pkgMethodLabel(current)
                if(nzchar(var)) assign(var, current, envir = merge)
                return(merge)
            }
        }
        else if(is.environment(obj)) {
            merge <- new.env()
            assign(.pkgMethodLabel(current), current, envir = merge)
            for(whatObj in objects(obj, all.names = TRUE))
                assign(whatObj, get(whatObj, envir = obj),
                       envir = merge)
            return(merge)
        }
        ## else adding a primitive, should do nothing
        else
            current
    }
}


.mlistAddToTable <- function(generic, mlist, table = new.env(TRUE, fenv), add = TRUE) {
  fenv <- environment(generic)
  signature <- generic@signature
    if(!exists(".SigLength", envir = fenv, inherits = FALSE))
     .setupMethodsTables(generic)
  n <- get(".SigLength", envir = fenv, inherits = FALSE)
  .storeMlist(table, rep("ANY", n), mlist, 1,  add, fenv)
  ## check for more args in the mlist than in the table
  nNow <- get(".SigLength", envir = fenv, inherits = FALSE)
  if(nNow > n) {
    length(signature) <- nNow
    .resetTable(table, nNow, signature)
  }
  table
}

.storeMlist <- function(table, sig, mlist, i, add, fenv) {
    ## once generic functions are installed from 2.11.0 or later, this should
    ## only be called with mlist a method or NULL.
    if(is.null(mlist)) return(table)
    m <- (if(is(mlist, "MethodsList")) mlist@methods
        else list(ANY=mlist)
        )
  ## once MethodsList is defunct, this should be rewritten (and renamed!)

  ## the methods slot is a list named by class, with elements either
  ## method definitions or mlists
  classes <- names(m)
  for(j in seq_along(m)) {
    el <- m[[j]]
    sig[[i]] <- classes[[j]]
    if(is(el, "MethodDefinition") || is.primitive(el)) {
      if(add)
        assign(.sigLabel(sig), el, envir = table)
      else
        remove(list = .sigLabel(sig), envir = table)
    }
    else if(is(el,"MethodsList")) {
      i1 <- i+1
      if(i1 >= length(sig)) {
        ## a reset of the labels will be needed
        assign(".SigLength", i1, envir = fenv)
        sig <- c(sig, rep("ANY", i1-length(sig)))
      }
      Recall(table, sig, el, i1, add, fenv)
    }
    else
      stop(gettextf("invalid mlist element for signature %s at level %d (should be methods list or method, had class %s)",
                    sQuote(classes[[j]]),
                    i,
                    dQuote(class(el))),
           domain = NA)
  }
  table
}

.cacheMethodInTable <- function(fdef, sig, def,
                                table = get(".AllMTable", envir = fenv)) {
    ## store method in cache table.
    ## called from setMethod()
    ## also Called from cacheMethod (from as(),
  ## as<-())
  fenv <- environment(fdef)
  if(missing(table) && !exists(".AllMTable", envir = fenv, inherits = FALSE))
    .setupMethodsTables(fdef)
  sig <- .matchSigLength(sig, fdef, fenv, TRUE)
  label <- .sigLabel(sig)
  isCurrent <- exists(label, envir = table, inherits = FALSE)
  if(is.null(def)) { # remove the method (convention for setMethod)
      if(isCurrent)
          remove(list = label, envir = table)
  }
  else {
      dupl <- .duplicateClassesExist()
      ## ensure that a valid object is assigned:  if duplicate classes
      ## exist, may need a table by package label; else, make sure
      ## the target and defined slots are complete
      ## IF we believed all methods up to date, the call could be conditional
      ##    if(dupl || isCurrent)
      def <- .methodPackageSlots(def, label, table, dupl, isCurrent)
      assign(label, def, envir = table)
  }
}

## check for duplicate classes and embed method in an environment if so
.methodPackageSlots <- function(def, ...) def

## the real version
..methodPackageSlots <- function(def, label, table, duplicatesExist, isCurrent) {
    sig <- def@target
    dups <- FALSE
    if(duplicatesExist) {
        def <- .fixPackageSlot(def, sig)
        for(cl in sig) {
            if(exists(cl, envir = .classTable, inherits = FALSE) && is.list(get(cl, envir = .classTable))) {
                dups <- TRUE
                break
            }
        }
        if(isCurrent) { # check that this is overwriting identical signature
            current <- get(label, envir = table)
            dups <- dups || !identical(current@target, sig)
        }
        if(dups) {
            if(isCurrent) {
                if(is(current, "MethodDefinition")) {
                    pkg <- attr(current@target, "package")
                    if(length(pkg) == 0)
                        current <- .fixPackageSlot(current, current@target)
                    env <- new.env()
                    ## zero-length seen 2011-07-29
                    var <- .pkgMethodLabel(current)
                    if(nzchar(var)) assign(var, current, envir = env)
                }
                else if(is.environment(current))
                    env <- current
                else
                    stop(
                         gettextf("bad method object stored in method table, class %s",
                                  dQuote(class(current))),
                         domain = NA)
            }
            else
                env <- new.env()
            assign(.pkgMethodLabel(def), def, envir = env)
            env
        }
        else # no change
            def
    }
    else # no duplicate classes
        def
}

.fixPackageSlot <- function(def, sig) {
    ## check the pkg slot
    pkgs <- attr(sig, "package")
    if(is.null(pkgs))
        pkgs <- character(length(sig))
    fixme <- !nzchar(pkgs)
    if(any(fixme)) {
        for(i in seq_along(pkgs)[fixme])
            pkgs[[i]] <- getClass(sig[[i]], .Force = TRUE)@package
        attr(sig, "package") <- pkgs
        def@target <- sig
        ## check the defined signature as well
        sig <- def@defined
        pkgs <- attr(sig, "package")
        if(is.null(pkgs))
            pkgs <- character(length(sig))
        fixme <- !nzchar(pkgs)
        if(any(fixme)) {
            for(i in seq_along(pkgs)[fixme])
                pkgs[[i]] <- getClass(sig[[i]], .Force = TRUE)@package
            attr(sig, "package") <- pkgs
            def@defined <- sig
        }
    }
    def
}

.okMethodLabel <- function(method) {
    if(is(method, "MethodDefinition")) {
        pkgs <- packageSlot(method@target)
        length(pkgs) > 0 && all(nzchar(pkgs))
    }
    else
        TRUE # primitive or environment
}


.pkgMethodLabel <- function(method) {
    sig <- method@target
    pkgs <- packageSlot(sig)
    if( (length(pkgs) < length(as.character(sig))) || any(!nzchar(pkgs)))
        stop("package slot missing from signature for generic ",
             sQuote(method@generic), "\n",
             "and classes ", paste(sig, collapse = ", "), "\n",
             "cannot use with duplicate class names (the package may need to be re-installed)",
             call. = FALSE, domain = NA)
    paste(pkgs, collapse = "#")
}

.resetTable <- function(table, n, signames) {
    ## protect this computation, in case it's resetting
    ## something used in the computation
    primMethods <- .allowPrimitiveMethods(FALSE)
    on.exit(.allowPrimitiveMethods(primMethods))
    ## after updating a methods table, the maximum no. of arguments in
    ## the signature increased to n.  Reassign any objects whose label
    ## does not match n classes from the defined slot
    anyLabel <- rep("ANY", n)
    anyPkg <- rep("methods", n)
    seqN <- 1L:n
    labels <- objects(envir=table, all.names = TRUE)
    for(what in labels) {
        method <- get(what, envir = table)
        if(is.primitive(method)) # stored as default ?
            newSig <- anyLabel
        else if(is(method, "MethodDefinition"))
            newSig <- method@defined
        else if(is(method, "environment")) {
            newSig <- strsplit(what, "#", fixed = TRUE)[[1]]
            .resetTable(method, n, signames)
        }
        else
            stop(gettextf("invalid object in methods table (%s), expected a method, got an object of class %s",
                          sQuote(what),
                          dQuote(class(method))),
                 domain = NA)

        if(is(method, "MethodDefinition")) {
            pkgs <- packageSlot(newSig)
            newSig <- as(ifelse(seqN > length(newSig), anyLabel, newSig), "signature")
            newSig@names <- signames
            newSig@package <-  ifelse(seqN > length(pkgs), anyPkg, pkgs)
            method@defined <- method@target <- newSig
            newLabel <- .sigLabel(newSig)
        }
        else
            newLabel <- .sigLabel(ifelse(seqN > length(newSig), anyLabel, newSig))
        remove(list=what, envir = table)
        assign(newLabel, method, envir = table)
    }
    NULL
}

### the tag associated with a method signature.
### Should perhaps use the same C code as dispatch, for consistency,
### however, that code breaks out early in the collapse loop if no match.
### This code is not used for quick matching, so efficiency less critical.
.sigLabel <- function(sig)
  paste(sig, collapse = "#")

## workhorse of selectMethod() [ -> ../Methods.R ] "
.findInheritedMethods <-
    function(classes, fdef, mtable = NULL,
             table = get(".MTable", envir = environment(fdef)),
             excluded = NULL, useInherited,
             simpleOnly = .simpleInheritanceGeneric(fdef), verbose = FALSE,
             doCache = is.environment(mtable))
{
    ## to avoid infinite recursion, and somewhat for speed, turn off S4 methods for primitives
    primMethods <- .allowPrimitiveMethods(FALSE)
    on.exit(.allowPrimitiveMethods(primMethods))
    ## classes is a list of the class(x) for each arg in generic
    ## signature, with "missing" for missing args
    if(!is.environment(table)) {
        if(is(fdef, "standardGeneric"))
          stop("invalid or unset methods table in generic function \"",
               fdef@generic,"\"", damain = NA)
        else
          stop("trying to find a methods table in a non-generic function")
    }
    hasGroup <- length(fdef@group) > 0L
    if(hasGroup)
      groupGenerics <- .getAllGroups(list(fdef))
    doExcluded <- length(excluded) > 0L
    if(verbose) {
	plist <- function(x) paste(x, collapse = ", ")
	cat(" .findInheritedMethods(): (hasGroup, doCache, doExcluded)= (",
	    plist(c("f","T")[1+c(hasGroup, doCache, doExcluded)]), ")\n",
	    if(hasGroup) paste0(" Group generics: ",
				plist(vapply(groupGenerics, slot,
					     character(1), "generic")), "\n"),
	    sep='')
    }
    nargs <- length(classes)
    if(!missing(useInherited) && length(useInherited) < nargs)
      useInherited <- rep(useInherited, length.out = nargs)
    if(hasGroup && !doExcluded) {
        ## first try for an exact match in a group generic
        ## If this matches &  is cached, it then will be treated as a non-inherited method
        ## so no further calls here should occur.
        ##
        ## doExcluded is the findNextMethod case; we don't regard group methods as
        ## inherited in the nextMethod sense, since they have the same signature
        label <- .sigLabel(classes)
        direct <- .getGroupMethods(label, groupGenerics, FALSE)
	if(length(direct)) {
	    if(doCache)
		assign(label, direct[[1L]], envir = mtable)
	    return(direct)
	}
        ## else, continue because we may want all defined methods
    }
    cl1 <- classes[[1L]]
    def <- getClass(cl1, .Force = TRUE)
    labels <-
      if(missing(useInherited) || useInherited[[1L]])
          c(cl1, .eligibleSuperClasses(def@contains, simpleOnly), "ANY")
      else cl1
    supersList <- list(labels)
    classDefs <- vector("list", nargs)
    classDefs[[1L]] <- def
    if(nargs > 1) { ## further arguments
        for(i in 2:nargs) {
            cc <- classDefs[[i]] <- getClass(classes[[i]], .Force = TRUE)
            allLabels <- if(missing(useInherited) || useInherited[[i]])
                c(cc@className, .eligibleSuperClasses(cc@contains, simpleOnly),
                  "ANY")
            else cc@className
            labels <- outerLabels(labels, allLabels)
            supersList <- c(supersList, list(allLabels))
        }
    }
    labels <- labels[-1L] # drop exact match
    labels <- unique(labels)# only needed while contains slot can have duplicates(!)
    if(verbose) {
	cat(" .fI> length(unique(method labels)) = ", length(labels))
	if(verbose >= 2) { cat(";  labels = \n") ; print(labels) }
    }
    allMethods <- objects(envir=table, all.names=TRUE)
    found <- match(labels, allMethods, 0L) > 0L
    nFound <- length(lab.found <- labels[found])
    methods <- list() # =?= vector("list", nFound) ; but fails??
    for(label in lab.found)
      methods[[label]] <- get(label, envir = table)
    if(verbose) cat(" >> found: ", nFound, "\n")
    if(hasGroup) {
        ##  add the  group methods recursively found but each time
        ## only those not already included in found.
        groupmethods <- .getGroupMethods(labels, groupGenerics, found)
        fromGroup <- c(rep(FALSE, length(methods)),
                       rep(TRUE,  length(groupmethods)))
        if(verbose) cat(" .fI> #{additional group methods}:",
                        length(groupmethods),"\n")
        methods <- c(methods, groupmethods)
    }
    else
      fromGroup <- rep(FALSE, length(methods))
    ## resolve any duplicate-class ambiguities
    if(.duplicateClassesExist()) {
        found <- integer()
        nm <- names(methods)
        for(i in seq_along(methods)) {
            m <- methods[[i]]
            if(is.environment(m)) {
                methods[[i]] <- .checkDuplicateMethodClasses(classDefs, m, nm[[i]])
                found <- c(found, i)
            }
        }
        if(length(found))
            methods <- unlist(methods, recursive = FALSE)
        if(!is.list(methods)) # reduced to a single method?
            methods <- list(methods)
    }
    if(doExcluded)
      methods <- methods[is.na(match(names(methods), as.character(excluded)))]
    ## remove default (ANY,..,ANY) if its not the only method:
    if(length(methods) > 1L) {
        defaultLabel <- paste(rep.int("ANY", nargs), collapse = "#")
        i <- match(defaultLabel, names(methods), 0L)
        if(i > 0L) {
            methods <- methods[-i]
            fromGroup <- fromGroup[-i]
        }
    }
    if(length(methods) > 1L) {
        if(verbose) cat(" .fI> length(methods) = ", length(methods),
                        " --> ambiguity\n")
        ## have ambiguity to resolve
        select <- .getBestMethods(methods, supersList, fromGroup, verbose=verbose)
        ##         --------------
        if(length(select) > 1L) {
            if(verbose) cat(" .fI> found", length(select)," best methods\n")

            target <- .sigLabel(classes)
            condAction <- getOption("ambiguousMethodSelection")
            if(is.null(condAction))
              condAction <- .ambiguousMethodMessage
            else if(!is(condAction, "function"))
              stop(gettextf("the \"ambiguousMethodSelection\" option should be a function to be called as the condition action; got an object of class %s",
                            dQuote(class(condAction))),
                   domain = NA)

            select <- withCallingHandlers(
                                          .disambiguateMethods(classes, select, fdef@generic,
                                                               methods, supersList, fromGroup,
                                                               classDefs, verbose),
                                          ambiguousMethodSelection=condAction)
        }
        methods <- methods[select]
    }
    if(simpleOnly && length(methods) == 0L) {
	## Seems to be *unused* [below, 'simpleOnly' argument was missing for years!]
	methods <- Recall(classes, fdef, mtable, table, excluded, useInherited,
			  simpleOnly, verbose, FALSE)
        if(length(methods) > 0L)
          message(gettextf("No simply inherited methods found for function %s; using non-simple method",
                           sQuote(fdef@generic)),
                  domain = NA)
    }
    if(length(methods)) {
        tlabel <- .sigLabel(classes)
        m <- methods[[1L]]
        if(is(m, "MethodDefinition"))  { # else, a primitive
            m@target <- .newSignature(classes, fdef@signature)
            ## if any of the inheritance is not simple, must insert coerce's in method body
            coerce <- .inheritedArgsExpression(m@target, m@defined, body(m))
            if(!is.null(coerce))
              body(m) <- coerce
            methods[[1L]] <- m
        }
	if(doCache) {
	    if(verbose) cat(" .fI> caching newly found methods ..\n")
	    assign(tlabel, m, envir = mtable)
	}
    }
    methods
}

.checkDuplicateMethodClasses <- function(classDefs, env, label){
    matches <- list()
    supers <- strsplit(label, "#", TRUE)[[1]]
    plabels <- strsplit(objects(env, all.names = TRUE), "#", TRUE)
    for(plabel in plabels) {
        if(.hasThisSubclass(classDefs, supers, plabel))
            matches[[plabel]] <- get(plabel, envir = env)
    }
    matches
}

.hasThisSubclass <- function(classDefs, supers, plabel) {
    for(i in seq_along(plabel)) {
        pkg <- classDefs[[i]]@package
        cl <- classDefs[[i]]@className
        si <- supers[[i]]
        pki <- plabel[[i]]
        if(identical(si, "ANY") ||
           (identical(cl, si) && identical(pkg, pki)))
            next
        cli <- getClassDef(si, package = pki)
        if(is.null(cli)) return(FALSE)
        sub <- cli@subclasses[[cl]]
        if(is.null(sub) || !identical(pkg, sub@package))
            return(FALSE)
    }
    TRUE
}

.ambiguousMethodMessage <- function(cond) {
  selected <- attr(cond, "selected")
  if(is.null(selected)) {# not properly set up, so just use the message
    message(cond$message)
  }
  else {
    possible <- attr(cond, "candidates")
    message(gettextf("Note: method with signature %s chosen for function %s,\n target signature %s.\n %s would also be valid",
                     sQuote(selected),
                     sQuote(attr(cond, "generic")),
                     sQuote(attr(cond, "target")),
		     paste0('"', possible[is.na(match(possible, selected))], '"',
			    collapse=", ")),
            domain = NA)
  }
}

.simpleInheritanceGeneric <- function(fdef) {
    identical(attr(fdef@signature, "simpleOnly"), TRUE)
}

.eligibleSuperClasses <- function(contains, simpleOnly) {
    what <- names(contains)
    if(!length(what))
      what
    else {
	eligible <-
	    sapply(contains,
		   if(simpleOnly)
		   function(x) (is.logical(x) && x) || x@simple
		   else # eliminate conditional inheritance
		   function(x) (is.logical(x) && x) || x@simple || identical(body(x@test), TRUE))
	what[eligible]
    }
}

.newSignature <- function(classes, names) {
  ## a simple version to deal with boostrapping stage, used in new() etc
    n <- min(length(classes), length(names))
  i <- seq_len(n)
    ## a corresponding set of package names
    ## <FIXME> There should be a "<unknown>" package name instead of "methods"
    ## but this requires a way to deal with that generally </FIXME>
    pkgs <- c(packageSlot(classes), rep("methods", n))[i]

  ## Simplified version ...
  .asS4(structure(as.character(classes)[i],
            class = .signatureClassName,
            names = as.character(names)[i],
            package = pkgs ))
 }

.findNextFromTable <- function(method, f, optional, envir, prev = character())
{
    fdef <- getGeneric(f)
    env <- environment(fdef)
    target <- method@target
    n <- get(".SigLength", envir = env)
    defined <- method@defined
    m <- length(defined)
    if(m > n)
        length(defined) <- n
  else if(n > m)
      ## will only really need this to be a signature when the elements
      ## have package attribute--see .sigLabel
      defined <-  new("signature", fdef, c(defined@.Data, rep("ANY", n-m)))
    excluded <- c(prev, .sigLabel(defined))
    methods <- .findInheritedMethods(defined, fdef, mtable = NULL, excluded = excluded)
    if(length(methods) == 0L) # use default method, maybe recursively.
        methods <- list(finalDefaultMethod(fdef@default)) #todo: put a label on it?
    if(length(methods) > 1L)
        warning(sprintf(ngettext(length(methods),
                                 "found %d equally good next method",
                                 "found %d equally good next methods"),
                        length(methods)),
                domain = NA)
    ## excluded slot is a list, but with methods tables, elements are just labels
    new("MethodWithNext", method, nextMethod = methods[[1L]],
        excluded = as.list(excluded))
}


  ## get the classes of the args

.InheritForDispatch <- function(classes, fdef, mtable) {
  methods <- .findInheritedMethods(classes, fdef, mtable)
  if(length(methods) == 1L)
    return(methods[[1L]]) # the method
  else if(length(methods) == 0L) {
    cnames <- paste0("\"", sapply(classes, as.character), "\"",
		     collapse = ", ")
    stop(gettextf("unable to find an inherited method for function %s for signature %s",
                  sQuote(fdef@generic),
                  sQuote(cnames)),
         domain = NA)
  }
  else
    stop("Internal error in finding inherited methods; didn't return a unique method", domain = NA)
}

.findMethodInTable <- function(signature, table, fdef = NULL)
{
    if(is(fdef, "genericFunction"))
        signature <- .matchSigLength(signature, fdef, environment(fdef), FALSE)
    label <- .sigLabel(signature)
##     allMethods <- objects(table, all.names=TRUE)
##     if(match(label, allMethods, nomatch = 0L))
    if(exists(label, envir = table, inherits = FALSE)) {
        value <- get(label, envir = table) ## else NULL
        if(is.environment(value)) {
            pkgs <- objects(value, all.names = TRUE)
            if(length(pkgs) == 1)
                value <- get(pkgs, envir = value)
            else if(length(pkgs) == 0)
                value <- NULL
            ## else, return the environment indicating multiple possibilities
        }
        value
    } # else, NULL
}

## inheritance distances:  0 for the class, 1 for immediate contains, 2 for other contains
##    and 3 for ANY
.inhDistances <- function(classDef) {
  contains <- classDef@contains
  allNames <-  unique(names(contains)) # bug allows duplicates in contains
  dist <- rep(2, length(allNames))
  for(i in seq_along(dist)) {
    ci <- contains[[i]]
    dist[[i]] <- ci@distance
  }
  dist <- c(0, dist, NA)
  names(dist) <- c(classDef@className, allNames, "ANY")
  dist
}

.leastMethodDistance <- function(methods, supersList, classDefs, fromGroup, verbose = FALSE) {
    n <- length(methods)
    dist <- rep(0, n)
    nArg <- length(classDefs)
    defClasses <- matrix("ANY", nArg, n)
    for(j in 1L:n) {
	cl <- methods[[j]]@defined@.Data
	defClasses[seq_along(cl), j] <- cl
    }
    containsDist <- lapply(classDefs, .inhDistances)
    maxDist <- max(unlist(containsDist), na.rm = TRUE) + 1
    if(verbose) { cat("** individual arguments' distances:\n"); print(containsDist) }
    ## add up the inheritance distances for each argument (row of defClasses)
    for(i in 1L:nArg) {
	ihi <- containsDist[[i]]
	ihi[is.na(ihi)] <- maxDist
	cli <- defClasses[i,]
	dist <- dist + ihi[match(cli, names(ihi))]
    }
    if(verbose) cat("** final methods' distances: (",
		    paste(formatC(dist), collapse= ", "), ")\n", sep='')
    best <- dist == min(dist)
    ## of the least distance methods, choose direct, rather than group
    ## methods, unless all the best methods are from group generics
    if(any(fromGroup[best]) && !all(fromGroup[best]))
	best <- best & !fromGroup
    (1:n)[best]
}

## currently called exactly once from .findInheritedMethods() :
.getBestMethods <- function(methods, supersList, fromGroup, verbose = FALSE) {
    n <- length(methods)      ## >= 2
    nArg <- length(supersList)## >= 1
    sigs <- matrix("ANY", nArg, n)
    for(i in 1:n) {
      sig <- methods[[i]]@defined
      if(length(sig) < nArg) { # is this still possible? --> show 'verbose'
	if(verbose) cat(sprintf(" .. method %d: length(sig) = %d < nArg = %d\n",
				i, length(sig), nArg))
	sigs[seq_along(sig), i] <- sig
      }
      else
        sigs[,i] <- sig
    }
    if(nArg < 2) { # the easy case
      return(which.min(match(sigs[1L,], supersList[[1L]])))
    }
    ## else  nArg >= 2
    best      <- rep.int(TRUE,  n)
    dominated <- rep.int(FALSE, n)
    pos <- matrix(0L, nArg, n)
    for(i in 1:nArg) {
        pos[i,] <- match(sigs[i,], supersList[[i]])
    }
    ## pairwise comparison of columns of pos.  Any way to vectorize?
    seqn <- seq_len(n)
    for(i in seqn) {
      for(j in seqn[-i]) {
        diffs <- pos[,j] - pos[,i]
	if(any(diffs < 0))  { best[i] <- FALSE; if(dominated[i]) break }
	if(all(diffs <= 0)) { dominated[i] <- TRUE; if(!best[i]) break }
      }
    }
    if(verbose)
	cat(if(any(best)) paste(" have best ones",
				paste(format(seqn[best]),collapse=","))
	    else if(any(dominated)) paste(" can eliminate dominated ones,",
				    paste(format(seqn[dominated]),collapse=",")),
	    "\n")
    ## a best method is as early in the superclasses as any other on all arguments
    ## Because the signatures are not duplicated, there can be at most one.
    if(any(best))
      seqn[best]
    ## eliminate those methods dominated by another
    else
      seqn[!dominated]
}

## currently called exactly once from .findInheritedMethods() :
.disambiguateMethods <- function(target, which, generic, methods, supersList,
                                 fromGroup, classDefs, verbose)
{
  ## save full set of possibilities for condition object
  candidates <- methods[which]
  note <- character()
  ## choose based on total generational distance
  which2 <- .leastMethodDistance(candidates, supersList, classDefs,
                                 fromGroup[which])
  if(length(which2) < length(which)) {
    note <- c(sprintf(ngettext(which2,
                               "Selecting %d method of minimum distance",
                               "Selecting %d methods of minimum distance"),
                      which2))
    which <- which[which2]
  }
  ## if some are group methods, eliminate those
  if(length(which) > 1 && any(fromGroup[which]) && !all(fromGroup[which])) {
    which <- which[!fromGroup]
    note <- c(note,  sprintf(ngettext(length(which),
                                      "Selecting %d non-group method",
                                      "Selecting %d non-group methods"),
                             length(which)))
  }
  ## prefer partially direct methods
  if(length(which) > 1) {
    direct <- sapply(methods[which], function(x, target)
                     (is(x, "MethodDefinition") && any(target == x@defined)),
                     target = target)
    if(any(direct) && !all(direct)) {
      which <- which[direct]
      note <- c(note, sprintf(ngettext(length(which),
                                       "Selecting %d partially exact-matching method",
                                       "Selecting %d partially exact-matching methods"),
                              length(which)))
    }
  }
  which <- which[[1L]]
  if(identical(as.character(generic), "coerce"))
      return(which) # as() computations not currently consistent w. selection (R 2.15.2)
  selected <- names(methods)[[which]]
  ## FIXME (?): This is not shown to the user
  msg <- sprintf(ngettext(length(candidates),
                          "Choosing method %s from %d ambiguous possibility",
                          "Choosing method %s from %d ambiguous possibilities"),
                 sQuote(selected), length(candidates))
  condObject <- simpleCondition(msg)
  ## would be nice to use an S4 class eventually
  class(condObject) <- c("ambiguousMethodSelection", class(condObject))
  attributes(condObject) <-
      c(attributes(condObject),
	list("candidates" = names(candidates),
	     "target"	  = .sigLabel(target),
	     "selected"	  = selected,
	     "generic"	  = generic,
	     "notes" = if(length(note)) paste(note, collapse ="; ") else ""))
  if(verbose) cat("   .disambiguateM*(): notes =\n\t",
		  attr(condObject, "notes"), "\n")
  signalCondition(condObject)
  which
}

# add objects to the generic function's environment that allow
# table-based dispatch of methods
.setupMethodsTables <- function(generic,
		initialize = !exists(".MTable", envir = env, inherits = FALSE))
{
    env <- environment(generic)
    if(initialize || !exists(".SigLength", envir = env, inherits = FALSE)) {
        nsig <- 1
        ## check that groups of generics agree on .SigLength; otherwise
        ## labels won't match
        for(gp in generic@group) {
            gpDef <- getGeneric(gp)
            if(is(gpDef, "genericFunction")) {
                .getMethodsTable(gpDef) # force initialization
                nsig <- max(nsig, get(".SigLength", envir = environment(gpDef)))
            }
        }
        assign(".SigLength", nsig, envir = env)
    }
    argSyms <- lapply(generic@signature, as.name)
    assign(".SigArgs", argSyms, envir = env)
    if(initialize) {
        mlist <- generic@default # from 2.11.0: method, primitive or NULL, not MethodsList
        mtable <- .mlistAddToTable(generic, mlist) # by default, adds to an empty table
        assign(".MTable", mtable, envir = env)
    }
    else ## the current .MTable
        mtable <- getMethodsForDispatch(generic)
    .resetInheritedMethods(env, mtable)
    if(is(generic, "groupGenericFunction")) {
        for(gp in generic@groupMembers) {
            gpDef <- getGeneric(gp)
            if(is(gpDef, "genericFunction"))
                .getMethodsTable(gpDef) # force initialization w. group methods
        }
    }
    NULL
}

.updateMethodsInTable <- function(generic, where, attach) {
  fenv <- environment(generic)
  reset <- identical(attach, "reset")
  if(!exists(".MTable", envir = fenv, inherits = FALSE))
    .setupMethodsTables(generic)
  mtable <- get(".MTable", envir = fenv)
  if(!reset) {
    env <- as.environment(where)
    tname <- .TableMetaName(generic@generic, generic@package)
    if(exists(tname, envir = env, inherits = FALSE)) {
      .mergeMethodsTable(generic, mtable, get(tname, envir = env), attach)
    }
    ## else used to warn, but the generic may be implicitly required
    ## by class inheritance, without any explicit methods in this package
  }
  if(length(generic@group)) {
      groups <- as.list(generic@group)
      generics <- vector("list", length(groups))
      for(i in seq_along(groups))
        generics[[i]] <- getGeneric(groups[[i]])
    .checkGroupSigLength(groups, generics)
  }
  if(is(generic, "groupGenericFunction")) {
      .checkGroupSigLength(list(generic@generic), list(generic))
      for(g in getGroupMembers(generic))
          .updateMethodsInTable(getGeneric(g), where, attach)
  }
  .resetInheritedMethods(fenv, mtable)
  mtable
}

.resetInheritedMethods <- function(fenv, mtable) {
    allObjects <- character()
    direct <- objects(mtable, all.names=TRUE)
    if(exists(".AllMTable", envir = fenv, inherits = FALSE)) {
        ## remove all inherited methods.  Note that code (e.g. setMethod) that asigns
        ## a new method to mtable is responsible for copying it to allTable as well.
        allTable <- get(".AllMTable", envir = fenv)
        allObjects <- objects(allTable, all.names=TRUE)
        remove(list= allObjects[is.na(match(allObjects, direct))], envir = allTable)
    }
    else {
        allTable <- new.env(TRUE, fenv)
        assign(".AllMTable", allTable, envir = fenv)
    }
    ## check for missing direct objects; usually a non-existent AllMTable?
    if(any(is.na(match(direct, allObjects)))) {
        direct <- objects(mtable, all.names=TRUE)
        for(what in direct)
          assign(what, get(what, envir = mtable), envir = allTable)
    }
    NULL
}

## In the following, consider separate "compute" and "print" functions/methods:
## Wish: alternative to 'classes' allow  "wild-card signature", e.g.,
##       showMethods("coerce", signature = c("dgeMatrix", "*"))
.showMethodsTable <- function(generic, includeDefs = FALSE, inherited = FALSE,
                              classes = NULL, showEmpty = TRUE, printTo = stdout())
{
    cf <- function(...) cat(file = printTo, sep = "", ...)
    sigString <- function(sig)
	paste0(names(sig), "=\"", as.character(sig), "\"", collapse = ", ")
    qs <- function(what) paste0('"', what, '"', collapse = ", ")
    doFun <- function(func, pkg) cf("Function: ", func, " (package ", pkg, ")\n")
    env <- environment(generic)
    signature <- generic@signature
    table <- get(if(inherited) ".AllMTable" else ".MTable", envir = env)
    f <- generic@generic
    p <- packageSlot(f)
    if(is.null(p)) p <- "base"
    deflt <- new("signature", generic, "ANY")
    labels <- objects(envir=table, all.names = TRUE)
    if(!is.null(classes) && length(labels)) {
	sigL <- strsplit(labels, split = "#")
	keep <- !sapply(sigL, function(x, y) all(is.na(match(x, y))), classes)
	labels <- labels[keep]
    }
    if(length(labels) == 0L) {
	if(showEmpty) {
	    doFun(f,p)
	    cf("<No methods>\n\n")
	}
	return(invisible())
    }
    ## else: non-empty methods list
    doFun(f,p)
    for(what in labels) {
	m <- get(what, envir = table)
        if(is.environment(m)) {  ## duplicate class case -- compare .findMethodInTable()
            pkgs <- objects(m)
            if(length(pkgs) == 1)
                m <- get(pkgs, envir = m)
            else if(length(pkgs) > 1)
                cf("  (", length(pkgs), " methods defined for this signature, with different packages)\n")
        }
	if( is(m, "MethodDefinition")) {
	    t <- m@target
	    if(length(t) == 0L)
		t <- deflt
	    d <- m@defined
	    if(length(d) == 0L)
		d <- deflt
	    cf(sigString(t), "\n")
	    if(!identical(t, d))
		cf("    (inherited from: ", sigString(d), ")\n")
            if(!.identC(m@generic, f) && length(m@generic) == 1L &&
               nzchar(m@generic))
		cf("    (definition from function \"", m@generic, "\")\n")
	}
	if(includeDefs && is(m, "function")) {
	    if(is(m, "MethodDefinition"))
		m <- m@.Data
	    cat(deparse(m), sep="\n", "\n", file = printTo)
	}
    }
    cat("\n", file = printTo)
}

## temporary switch for tables
useMTable <- function(onOff = NA)
  .Call(C_R_set_method_dispatch, as.logical(onOff))

## get all the group generic functions, in breadth-first order since
## direct group inheritance is closer than indirect (all existing
## groups are mutually exclusive, but multiple group membership is
## allowed)
.getAllGroups <- function(funs) {
  start <- length(funs)
  for(i in seq_along(funs)) {
    groups <- funs[[i]]@group
    funs <- c(funs, lapply(groups,
                           function(what) {
                             f <- getGeneric(what)
                             if(!is.function(f))
                               stop("failed to find expected group generic function: ",
                                    what)
                             f
                           }))
    }
    ## now the next generations recusively
    if(length(funs) > start) {
      nmore <- length(funs) - start
      more <- Recall(funs[(start+1):length(funs)])
      ## did we add any groups?
      if(length(more) > nmore)
        funs <- c(funs, more[(nmore+1):length(more)])
    }
    funs
}

.getGroupMethods <- function(labels, generics, found) {
  methods <- list()
  for(i in seq_along(generics)) {
    gen <- generics[[i]]
    if(!is(gen,"genericFunction"))
      stop(gettextf("invalid group generic function in search for inherited method (class %s)",
                    dQuote(class(gen))),
           domain = NA)
    table <- .getMethodsTable(gen)
    allMethods <- objects(envir=table, all.names = TRUE)
    ## TODO:  possible for .SigLength to differ between group &
    ## members.  Requires expanding labels to max. length
    newFound <- rep(FALSE, length(found))
    newFound[!found] <- (match(labels[!found], allMethods, 0L) > 0L)
    found <- found | newFound
    for(what in labels[newFound])
      methods[[what]] <- get(what, envir = table)
  }
  methods
}

.getMethodsTable <- function(fdef, env = environment(fdef),
                             check = TRUE, inherited = FALSE)
{
    name <- if(inherited) ".AllMTable" else ".MTable"
    if(check && !exists(name, envir = env, inherits = FALSE)) {
	.setupMethodsTables(fdef, initialize = TRUE)
	if(!exists(name, envir = env, inherits = FALSE))
	    stop("invalid methods table request")
    }
    get(name, envir = env)
}

.getGenericSigLength <- function(fdef, env = environment(fdef), check = TRUE) {
    if(check && !exists(".SigLength", envir = env, inherits = FALSE))
      .setupMethodsTables(fdef)
    get(".SigLength", envir = env)
}



.checkGroupSigLength <- function(gnames, generics = lapply(gnames, getGeneric)) {
  funs <- gnames
  recall <- FALSE
  for(i in seq_along(gnames)) {
    what <- gnames[[i]]
    fdef <- generics[[i]]
    if(!is(fdef, "groupGenericFunction")) {
      warning(gettextf("trying to check signature length of group generic '%s', but it is not a group generic", what),
              domain = NA)
      next
    }
    if(length(fdef@group))  {# push up the check one level
      gnames[[i]] <- fdef@group
      generics[[i]] <- lapply(fdef@group, getGeneric)
      recall <- TRUE
      next
    }
    funs <- c(funs, getGroupMembers(fdef, TRUE, FALSE))
  }
  if(recall)
    return(Recall(unlist(gnames, FALSE), unlist(generics, FALSE)))
  funs <- unique(funs)
  fdefs <- lapply(funs, function(x) {
    if(is.character(x) && length(x) == 1L) getGeneric(x)
    else x})
  ## now compare the sig lengths
  sigs <- rep(0,length(funs))
  for(i in seq_along(sigs)) {
    what <- funs[[i]]
    fdef <- fdefs[[i]]
    if(is.null(fdef))
      next # getGroupMembers returns NULL if  member is not defined
    if(!is(fdef, "genericFunction"))
      warning(gettextf("trying to check signature length of generic '%s', but it is not a generic function: i = %d, funs = %s, gnames = %s",
                       what,  i, paste(unlist(funs), collapse = ", "),
                       paste(as.character(gnames), collapse = ", ")),
              domain = NA)
    else {
      ev <- environment(fdef)
      if(!exists(".SigLength", envir = ev, inherits = FALSE))
        .setupMethodsTables(fdef)
      sigs[i] <- get(".SigLength", envir = ev)
    }
  }
  n <- max(sigs)
    reset <- sigs < n & sigs > 0 # all the  sigs  for defined funs & less than max.
  if(any(reset)) {
    funs <- funs[reset]
    fdefs <- fdefs[reset]
    for(fdef in fdefs) {
        .resetSigLength(fdef, n)
    }
  }
  funs
}

## a simplified outer of paste
outerLabels <- function(labels, new) {
    ## WARNING: This code incorporates the definition of .sigLabel
    ## and so must change if that does (e.g. to include package)
    n <- length(labels)
    m <- length(new)
    paste(labels[rep.int(1L:n, rep.int(m,n))], new[rep.int(1L:m,n)], sep ="#")
}


.matchSigLength <- function(sig, fdef, fenv, reset = FALSE) {
  nargs <- .getGenericSigLength(fdef, fenv, TRUE)
  n <- length(sig)
  pkgs <- packageSlot(sig)
  if(n < nargs) {
      more <- nargs - n
      pkgs <- c(pkgs, rep("methods", more))
      sig <- c(as.character(sig), rep("ANY", more))
  }
  else if(n > nargs) { #reset table?
    if(all(sig[(nargs+1):n] == "ANY"))
      length(sig) <- length(pkgs) <- nargs
    else {
      while(sig[[n]] == "ANY")
        n <- n-1
      if(reset)
        .resetSigLength(fdef, n)
      length(sig) <- length(pkgs) <- n
    }
  }
  packageSlot(sig) <- pkgs
  sig
}

.resetSigLength <- function(fdef, n) {
    fenv <- environment(fdef)
    assign(".SigLength", n, envir = fenv)
    mtable <- .getMethodsTable(fdef, fenv, check = FALSE)
    signames <- fdef@signature
    length(signames) <- n
    .resetTable(mtable, n, signames)
    .resetInheritedMethods(fenv, mtable)
}

.TableMetaName <- function(name, package)
  methodsPackageMetaName("T", paste(name, package, sep=":"))

.TableMetaPrefix <- function()
    methodsPackageMetaName("T","")

# regexp for matching table names; semi-general but assumes the
# meta pattern starts with "." and has no other special characters
.TableMetaPattern <- function()
    paste0("^[.]",substring(methodsPackageMetaName("T",""),2))

.addToMetaTable <- function(fdef, signature, definition, where, nSig) {
  return()
}

## the real version
..addToMetaTable <- function(fdef, signature, definition, where,
                             nSig = .getGenericSigLength(fdef)) {
    ## TODO:  nSig should be a slot in the table
  tname <- .TableMetaName(fdef@generic, fdef@package)
  where <- as.environment(where)
  if(exists(tname, envir =where, inherits = FALSE)) {
     table <- get(tname, envir = where)
     if(length(signature) > nSig)
       .resetTable(table, length(signature), fdef@signature[seq_along(signature)])
  }
  else {
    table <- new.env(TRUE, environment(fdef))
    assign(tname, table, envir = where)
  }
  .cacheMethodInTable(fdef, signature, definition, table)
}

## Assertion: following is unused
.assignMethodsMetaTable <- function(mlist, generic, where, overwrite = TRUE) {
    .MlistDeprecated(".assignMethodsMetaTable")
    tname <- .TableMetaName(generic@generic, generic@package)
    if(overwrite || !exists(tname, envir = where, inherits = FALSE)) {
        table <- .mlistAddToTable(generic, mlist) # asserted never to be called.
        assign(tname, table, envir = where)
    }
}

.removeMethodsMetaTable <- function(generic, where) {
    ## does not warn if none exists, on the theory that a generic may be created
    ## but no methods defined to create a table.  The use of implicitGeneric's is an example.
    tname <- .TableMetaName(generic@generic, generic@package)
    if(exists(tname, where, inherits = FALSE))
      rm(list=tname, pos = where)
}

.getGenericSigArgs <- function(fdef, env = environment(fdef), check = TRUE) {
    if(check && !exists(".SigLength", envir = env, inherits = FALSE))
      .setupMethodsTables(fdef)
    n <- get(".SigLength", envir = env)
    args <-  get(".SigArgs", envir = env)
    length(args) <- n
    args
}


## the most simple part of listFromMethods() below; not yet exported
tableNames <- function(generic, where, table) {
    fdef <- getGeneric(generic)
    if(missing(table))
	table <-
	    if(missing(where)) .getMethodsTable(fdef)
	    else get(.TableMetaName(fdef@generic, fdef@package),
                     envir = as.environment(where), inherits = FALSE)
    objects(envir=table, all.names=TRUE)
}

listFromMethods <- function(generic, where, table) {
    fdef <- getGeneric(generic)
    if(missing(table))
	table <-
	    if(missing(where)) .getMethodsTable(fdef)
	    else get(.TableMetaName(fdef@generic, fdef@package),
		     envir = as.environment(where), inherits = FALSE)
    fev <- environment(fdef)
    nSigArgs <- .getGenericSigLength(fdef, fev)
    names <- objects(envir=table, all.names=TRUE)
    methods <- lapply(names, function(x)get(x, envir = table))
    if(nSigArgs > 1) {
        n <- length(names)
        sigs <- vector("list", n)
        namesCon <- textConnection(names)
        for(i in seq_len(n))
            sigs[[i]] <- scan(namesCon, "", sep ="#", nmax = nSigArgs, quiet=TRUE)
    }
    else
      sigs <- as.list(names)
    new("LinearMethodsList", classes=sigs, methods=methods,
        arguments = .getGenericSigArgs(fdef, fev), generic = fdef)
}

.makeMlist1 <- function(arg, objects, j = 1) {
    mnames <- character(length(objects))
    for(i in seq_along(objects)) {
        what <- objects[[i]]
        if(is.primitive(what))
          sig <- "ANY"
        else
          sig <- what@defined
        mnames[[i]] <- (if(length(sig) < j) "ANY" else sig[[j]])
    }
    names(objects) <- mnames
    new("MethodsList", argument = arg, methods = objects, allMethods = objects)
}

.makeMlist2 <- function(args, objects, j = 1) {
    ## make a list according to  argument j, convert these as needed
    mlists <- list()
    for(what in objects) {
        sig <- if(!is.primitive(what)) what@defined # else NULL
        if(length(sig) <= j)
            arg1 <- arg2 <- "ANY"
        else {
            arg1 <- sig[[j]]
            arg2 <- sig[[j+1]]
        }
        x <- list(what)
        el <- mlists[[arg1, exact = TRUE]]
        mlists[[arg1]] <- (if(is.null(el)) x else c(el, x))
    }
    jNext <- j+1
    if(jNext < length(args))
      for(i in seq_along(mlists))
          mlists[[i]] <- .makeMlist2(args, mlists[[i]], jNext)
    else {
        arg2 <- as.name(args[[jNext]])
        for(i in seq_along(mlists))
          mlists[[i]] <- .makeMlist1(arg2, mlists[[i]], jNext)
    }
    new("MethodsList", argument = as.name(args[[1L]]),
        methods = mlists, allMethods = mlists)
}

.makeMlistFromTable <- function(generic, where = NULL) {
    .getAll <- function(what, table) {
        value <- list(length(what))
        for(i in seq_along(what))
          value[[i]] <- get(what[[i]], envir = table)
        value
    }
    if(is.null(where)) {
        what <- ".MTable"
        where <- environment(generic)
    }
    else {
        where <- as.environment(where)
        what <-  .TableMetaName(generic@generic, generic@package)
    }
    if(exists(what, envir = where, inherits= FALSE))
        table <- get(what, envir = where)
    else
        table <- new.env()
    value <- new("MethodsList", argument = as.name(generic@signature[[1]]))
    allNames <- objects(envir=table, all.names = TRUE)
    if(length(allNames) == 0L)
      return(value)
    argNames <- generic@signature
    ## USES THE PATTERN OF class#class#.... in the methods tables
    nargs <- nchar(unique(gsub("[^#]","", allNames)))+1
    if(length(nargs) > 1L) {
        warning("something weird:  inconsistent number of args in methods table strings:", paste(nargs,collapse = ", ")," (using the largest value)",
                domain = NA)
        nargs <- max(nargs)
    }
    length(argNames) <- nargs # the number of args used
    if(nargs == 1)
        .makeMlist1(as.name(argNames[[1L]]), .getAll(allNames, table))
    else
      .makeMlist2(argNames, .getAll(allNames, table))
 }

## assign a methods meta-data table, by default (and usually) a copy of the table
## from the generic function with the initial methods, if any.
.assignMethodsTableMetaData <- function(name, generic, where, table) {
    what <-  .TableMetaName(generic@generic, generic@package)
    if(missing(table))
          table <- .copyEnv(.getMethodsTable(generic))
    assign(what, table, envir = as.environment(where))
}

.getMethodsTableMetaData <-  function(generic, where, optional = FALSE) {
    what <-  .TableMetaName(generic@generic, generic@package)
    if(exists(what, envir = where, inherits = FALSE))
      get(what, envir = where )
    else if(optional)
      NULL
    else
      stop(gettextf("no methods table for generic %s from package %s in package %s",
                    sQuote(generic@generic),
                    sQuote(generic@package),
                    sQuote(getPackageName(where))),
           domain = NA)
}

.inheritedArgsExpression <- function(target, defined, body) {
    expr <- substitute({}, list(DUMMY = "")) # bug if you use quote({})--is overwritten!!
    args <- names(defined)
    for(i in seq_along(defined)) {
        ei <- extends(target[[i]], defined[[i]], fullInfo = TRUE)
        if(is(ei, "SClassExtension")  && !ei@simple)
          expr[[length(expr) + 1L]] <-
            substitute(ARG <- as(ARG, DEFINED, strict = FALSE),
                       list(ARG = as.name(args[[i]]),
                            DEFINED = as.character(defined[[i]])))
    }
    if(length(expr) > 1L) {
       expr[[length(expr) + 1L]] <- body
       expr
   }
    else
      NULL
}

testInheritedMethods <- function(f, signatures, test = TRUE,  virtual = FALSE,
                                 groupMethods = TRUE,  where = .GlobalEnv)
{
  getSigs <- function(fdef)
      objects(methods:::.getMethodsTable(fdef), all.names = TRUE)

  ## Function relevantClasses is defined here to set object .undefClasses
  ## in testInheritedMethods as a marker to warn about undefined subclasses
  .relevantClasses <- function(classes, excludeVirtual, where, doinheritance) {
    classDefs <- lapply(classes, getClassDef, where)
    undefs <- sapply(classDefs, is.null)
    if(any(undefs)) {
      .undefClasses <<- unique(c(.undefClasses, classes[undefs]))
      classes <- classes[!undefs]
      classDefs <- classDefs[!undefs]
    }
    if(doinheritance) {
      allSubs <- lapply(classDefs,  function(what) names(what@subclasses))
      allSubs <- unique(unlist(allSubs))
      pattern <- sapply(allSubs, .matchSubsPattern, classes, excludeVirtual)
      ## exclude virtuals
      if(excludeVirtual) {
        excl <- nzchar(pattern)
        pattern <- pattern[excl]
        allSubs <- allSubs[excl]
      }
      if(length(allSubs)>0)
        allSubs <- sapply(split(allSubs, pattern), `[[`,1)
      else
        allSubs <- character()
    }
    else
      allSubs <- character()
    ## prepend the classes themselves, as appropriate
    iAny <- match( "ANY", classes, 0)
    if(iAny > 0) {
      classes[[iAny]] <- ".Other" # non-virtual placeholder for ANY
      classDefs[[iAny]] <- getClassDef(".Other")
    }
    if(excludeVirtual)
      classes <- classes[sapply(classDefs, function(def) identical(def@virtual, FALSE))]
    unique(c(classes, allSubs))
  }
  ## end of .relevantClasses

  if(!is(f, "genericFunction"))
    f <- getGeneric(f)
  fname <- f@generic
  if(missing(signatures)) {
    mdefs <- findMethods(f)
    mnames <- names(mdefs)
    sigs <-  findMethodSignatures(methods = mdefs)
    if(groupMethods) {
      groups <- getGroup(f, recursive = TRUE)
      for(group in groups) {
        fg <- getGeneric(group)
        mg <- findMethods(fg)
        sigsg <- findMethodSignatures(methods = mg)
        newSigs <- is.na(match(names(mg), mnames))
        mg <- mg[newSigs]
        mdefs <- c(mdefs, mg[newSigs])
        sigs <- rbind(sigs, sigsg[newSigs,])
        mnames <- c(mnames, names(mg)[newSigs])
      }
    }
    if(length(sigs) == 0)
      return(new("MethodSelectionReport", generic = fname))
    ## possible selection of which args to include with inheritance
    ok <- if(fname %in% c("coerce", "coerce<-"))
	match(colnames(sigs), "from", 0) > 0 else rep.int(TRUE, ncol(sigs))
    for(j in seq_len(ncol(sigs))) {
      classesj <- unique(sigs[,j])
      .undefClasses <- character()
      subclasses <- .relevantClasses(classesj, !virtual, where, ok[[j]])
      nj <- length(subclasses)
      ##       if(nj == 0) {  ##FIXME, wrong test
      ##         warning(gettextf("No eligible subclasses for argument '%s' found, so no contribution to analysis",
      ##                          colnames(sigs)[[j]]), domain  = NA)
      ##         next
      ##       }
      if(j > 1) {
        ## replicate all the previous elements of subclasses a la outer
        subclasses <- rep(subclasses, rep.int(ncomb, nj))
        ncomb <- ncomb * nj
        sigLabels <- paste(rep(sigLabels, times = nj), subclasses, sep = "#")
      }
      else {
        sigLabels <- subclasses
        ncomb <- nj
      }

      if(length(.undefClasses)) {
        warning(gettextf("undefined classes (%s) will be ignored for argument '%s'",
                         paste0('"',unique(.undefClasses),'"', collapse=", "),
                         colnames(sigs)[[j]]), domain = NA)
        .undefClasses <- character()
      }
    } ## loop on j
    ## now split the individual labels back into signatures
    signatures <- strsplit(sigLabels, "#", fixed = TRUE)
  } ## end of missing(signatures) case
  else if(is(signatures, "matrix") && identical(typeof(signatures), "character")
       && ncol(signatures) <= length(f@signature)) {
      ## turn signatures back into a list
      siglist <- vector("list", nrow(signatures))
      for(i in seq_len(nrow(signatures)))
        siglist[[i]] <- signatures[i,]
      signatures <- siglist
  }
  else stop("argument 'signatures' must be a character matrix whose rows are method signatures")
  ambig_target <- character()
  ambig_candidates <- list()
  ambig_selected <- character()
  ambig_note <- character()
  if(test) {
    ## define a handler that accumulates the attributes from the condition object
    warninghandler <- function(cond) {
      ambig_target <<- c(ambig_target, attr(cond, "target"))
      ambig_candidates <<- c(ambig_candidates, list(attr(cond, "candidates")))
      ambig_selected <<- c(ambig_selected, attr(cond, "selected"))
      ambig_note <<- c(ambig_note, attr(cond, "note"))
    }
    ambigOpt <- options(ambiguousMethodSelection = warninghandler)
    on.exit(options(ambigOpt))
    doSelect <-  function(sig) {
      x <- selectMethod(f = f, sig, optional = TRUE)
      if(is(x, "MethodDefinition")) {
        nsig <- x@defined
        if(length(nsig) < length(sig))
          c(nsig, rep("ANY", length(sig) - length(nsig)))
        else
          nsig
      }
      else if(is.null(x))
        rep("<NONE>", length(sig))
      else # primitive
        rep("ANY", length(sig))
    }
    signatures <- lapply(signatures, doSelect)
  }
  signatures <- sapply(signatures, paste0, collapse = "#")
  names(signatures) <- sigLabels

  new("MethodSelectionReport", generic = fname, allSelections = signatures,
      target = ambig_target, selected = ambig_selected,
      candidates = ambig_candidates, note = ambig_note)
}

.matchSubsPattern <- function(what, matchto, excludeVirtual) {
  def <- getClass(what)
  if(excludeVirtual & def@virtual)
    return("")
  matches <- match(names(def@contains), matchto, 0)
  matches <- matches[matches>0]
  paste(matches, collapse=".")
}

#  File src/library/methods/R/oldClass.R
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

## assumes oldClass has been defined as a virtual class

setOldClass <- function(Classes, prototype = NULL,
                        where = topenv(parent.frame()), test = FALSE,
                        S4Class) {
    simpleCase <- is.null(prototype)
    mainClass <- Classes[[1L]]
    prevDef <- getClassDef(mainClass, where, inherits = FALSE)
    if(!missing(S4Class)) {
        if(test)
          stop("not allowed to have test==TRUE and an S4Class definition")
        if(!is(S4Class, "classRepresentation")) {
            if(is.character(S4Class)) {
                clName <- S4Class
                S4Class <- getClass(S4Class)
                if(.identC(clName, Classes[[1L]]))
                  removeClass(clName, where = where) # so Recall() will work
            }
            else
              stop(gettextf("argument 'S4Class' must be a class definition: got an object of class %s",
                            dQuote(class(S4Class))),
                   domain = NA)
        }
        if(!is.null(prototype)) {
            S4prototype <- S4Class@prototype
            ## use the explicit attributes from the supplied argument, else S4prototype
            S4Class@prototype <- .mergeAttrs(prototype, S4prototype)
        }
        ## register simple S3 class(es), including main class, if it's not defined already
        Recall(Classes, where = where)
        return(.S4OldClass(Classes[[1L]], if(length(Classes) > 1) Classes[[2L]] else "oldClass", S4Class, where, prevDef))
    }
    if(test)
        return(.setOldIs(Classes, where))
    if(!is.null(prevDef)) {
        on.exit(.restoreClass(prevDef, where))
        removeClass(mainClass, where = where) # so Recall() will work
    }
    prevClass <- "oldClass"
    S3Class <- character()  #will accumulate the S3 classes inherited
    ## The table of S3 classes, used
    ## to convert S4 objects in S3 method dispatch.
    ## TODO:  should provide an optional argument to setOldClass()
    ## to prevednt this conversion if it's not needed
    if(!exists(".S3MethodsClasses", envir = where, inherits = FALSE)) {
      S3table <- new.env()
      assign(".S3MethodsClasses", S3table, envir = where)
    }
    else S3table <- get(".S3MethodsClasses", envir = where)
    dataPartClass <- NULL
    for(cl in rev(Classes)) {
       S3Class <- c(cl, S3Class)
        if(isClass(cl, where)) {
            def <- getClass(cl, where)
            if(!extends(def, prevClass)) {
                ## maybe an object type or other valid data part
                cl1 <- .validDataPartClass(cl, where, dataPartClass)
                if(is.null(cl1))
                  stop(gettextf("inconsistent old-style class information for %s; the class is defined but does not extend %s and is not valid as the data part",
                                dQuote(cl),
                                dQuote(prevClass)),
                       domain = NA)
                else dataPartClass <- cl1
              }
            else {
              prevP <- def@prototype
              if(missing(prototype))
                prototype <- prevP # keep track of inherited prototype for use in mainClass
              prevS3Class <- attr(prevP, ".S3Class")
              if(length(prevS3Class) > length(S3Class)) #implies cl is registered S3 class
                S3Class <- prevS3Class
            }
        }
        else {
            useP <- TRUE
            if(cl != mainClass || simpleCase) {
                setClass(cl, contains = c(prevClass, "VIRTUAL"), where = where)
            }
            else if(isClass(class(prototype)))
                setClass(cl, contains = prevClass, prototype = prototype, where = where)
            else { #exceptionally, we allow an S3 object from the S3 class as prototype
                if(.class1(prototype) != mainClass)
                  stop(gettextf('the S3 class of the prototype, "%s", is undefined; only allowed when this is the S3 class being registered ("%s")', .class1(prototype), mainClass), domain = NA)
                setClass(cl, contains = prevClass, where = where)
                useP <- FALSE
            }
            def <- getClassDef(cl, where)
            if(useP) clp <- def@prototype else clp <- prototype
            attr(clp, ".S3Class") <- S3Class
            def@prototype <- .notS4(clp)
            assignClassDef(cl, def, where = where)
            ## add the class to the table of S3 classes
            assign(cl, def, envir= S3table)
        }
       prevClass <- cl
    }
    if(!is.null(prevDef)) # cancel error action
      on.exit()
}

.restoreClass <- function(def, where) {
    cl <- def@className
    message(gettextf("restoring definition of class %s", dQuote(cl)),
            domain = NA)
    if(isClass(cl, where = where))
       removeClass(cl, where = where)
    assignClassDef(cl, def, where = where)
}

.S4OldClass <- function(Class, prevClass, def,where, prevDef) {
    ## def is the S4 version of this class def'n, maybe by another class
    ## name, and may or may not already extend oldClass
    curDef <- getClassDef(Class, where) # asserted to be defined
    ## arrange to restore previous definition if there was one.  Also done in setOldClass
    ## when no S4Class argument supplied
    if(!is.null(prevDef)) {
        on.exit(.restoreClass(prevDef, where))
        removeClass(Class, where = where) # so Recall() will work
    }
    if(!identical(def@className, curDef@className))
      def <- .renameClassDef(def, curDef@className)
    ## check that any common slots will give a valid S3 object
    .validS3Extends(def, curDef)
    def@slots <- c(def@slots, curDef@slots)
    ext <- c(def@contains, curDef@contains)
    ## correct ordering & duplicate resolution: copied from .walkClassGraph
    distOrder <- sort.list(sapply(ext, function(x)x@distance))
    ext <- ext[distOrder]
    if(anyDuplicated(names(ext)))
        ext <- .resolveSuperclasses(def, ext, where)
    def@contains <- ext
    subcls <- curDef@subclasses
    if(length(subcls) > 0) {
      def@subclasses[names(subcls)]  <- subcls
    }
    proto <- def@prototype
    if(is.null(attr(proto, ".S3Class"))) { # no S3 class slot, as will usually be true
        attr(proto, ".S3Class") <- if(.identC(prevClass, "oldClass")) Class else S3Class(curDef@prototype)
        def@prototype <- proto
    }
    assignClassDef(Class, def, where = where)
    ## allow an existing superclass relation to remain (it may have a coerce method)
    ## Otherwise, create a simple transformation, which relies on consistency
    ## in the slots.
    if(!extends(def, prevClass, maybe = FALSE))
      setIs(Class, prevClass, classDef = def, where = where)
    slotsMethod <- function(object) NULL
    body(slotsMethod) <- substitute({LIST}, list(LIST = def@slots))
    setMethod("slotsFromS3", Class, slotsMethod, where = where)
    if(!is.null(prevDef)) # cancel error action
      on.exit()
}

.validS3Extends <- function(classDef1, classDef2) {
    slots2 <- classDef2@slots
    if(length(slots2) > 0) {
        n2 <- names(slots2)
        slots1 <- classDef1@slots
        n1 <- names(slots1)
        bad <- character()
        for(what in n2[match(n2, n1, 0) > 0])
          if(!extends(elNamed(slots1, what), elNamed(slots2, what))) {
              message(gettextf("slot %s: class %s should extend class %s",
                               sQuote(what),
                               dQuote(elNamed(slots1, what)),
                               dQuote(elNamed(slots2, what))),
                      domain = NA)
              bad <- c(bad, what)
          }
        if(length(bad)>0)
          stop(
               gettextf("invalid S4 class corresponding to S3 class: slots in  S4 version must extend corresponding slots in S3 version: fails for %s",
                        paste0('"', bad, '"',  collapse = ", ")),
               domain = NA)
    }
    TRUE
}

##.initS3Classes will make this generic, with a method for "oldClass"
slotsFromS3 <- function(object) {
    list()
}

utils::globalVariables("CLASS")

.oldTestFun <- function(object) CLASS %in% attr(object, "class")
.oldCoerceFun <- function(from, strict = TRUE) {
    if(strict)
        stop(gettextf("explicit coercion of old-style class (%s) is not defined", paste(class(from), collapse = ", ")), domain = NA)
    from
}
.oldReplaceFun <- function(from, to, value)
    stop(gettextf("explicit replacement not defined for as(x, \"%s\") <- value for old-style class %s",
                  to, dQuote(class(from)[1L])),
         domain = NA)

## the inheritance of these S3 classes must be decided on a per-instance
## basis.  At one time, there were classes in base/stats that had this
## property, (e.g., POSIXt, POSIX{cl}t) but apparently no longer.
## The possibility is still allowed
## for user-defined S3 classes.
.setOldIs <- function(Classes, where) {
    if(length(Classes) != 2)
        stop(gettextf("argument 'Classes' must be a vector of two classes; got an argument of length %d", length(Classes)), domain = NA)
    for(cl in Classes) {
        if(isClass(cl, where)) {
            if(!extends(cl, "oldClass"))
                warning(gettextf("inconsistent old-style class information for %s (maybe mixing old and new classes?)",
                                 dQuote(cl)), domain = NA)
        }
        else
            setClass(cl, representation("oldClass", "VIRTUAL"), where = where)
    }
    Class1 <- Classes[[1L]]
    for(cl in Classes[-1L]) {
        tfun <- .oldTestFun
        body(tfun, envir = environment(tfun)) <-
            substitute(inherits(object, CLASS), list(CLASS = cl))
        setIs(Class1, cl, test = tfun, coerce = .oldCoerceFun,
              replace = .oldReplaceFun, where = where)
    }
    NULL
}

isXS3Class <- function(classDef) {
    ".S3Class" %in% names(classDef@slots)
}

S3Class <- function(object) {
    value <- attr(object, ".S3Class")
    if(is.null(value)) {
        if(isS4(object)) {
            if(is.na(match(".Data", names(getClass(class(object))@slots))))
                stop(gettextf("'S3Class' only defined for extensions of %s or classes with a data part:  not true of class %s",
                              dQuote("oldClass"),
                              dQuote(class(object))),
                     domain = NA)
            class(getDataPart(object))
        }
        else
          class(object)
    }
    else
      value
}

.S3Class <- S3Class # alias for functions with S3Class as an argument

.addS3Class <- function(class, prototype, contains, where) {
    for(what in contains) {
        whatDef <- getClassDef(what@superClass, where = where)
        if(isXS3Class(whatDef))
          class <- c(class, attr(whatDef@prototype, ".S3Class"))
    }
    attr(prototype, ".S3Class") <- unique(class)
    prototype
}

"S3Class<-" <- function(object, value) {
    if(isS4(object)) {
        current <- attr(object, ".S3Class")
        if(is.null(current)) {
            if(is.na(match(value, .BasicClasses)))
               stop(gettextf("'S3Class' can only assign to S4 objects that extend \"oldClass\"; not true of class %s",
                             dQuote(class(object))),
                    domain = NA)
            mode(object) <- value ## may still fail, a further check would be good
        }
        else
          slot(object, ".S3Class") <- value
    }
    else
      class(object) <- value
    object
}

## rename a class definition:  needs to change if any additional occurences of class
## name are added, other than the className slot and the super/sub class names
## in the contains, subclasses slots respectively.
.renameClassDef <- function(def, className) {
    oldName <- def@className
    validObject(def) # to catch any non-SClassExtension objects
    def@className <- className
    comp <- def@contains
    for(i in seq_along(comp))
        comp[[i]]@subClass <- className
    def@contains <- comp
    comp <- def@subclasses
    for(i in seq_along(comp))
        comp[[i]]@superClass <- className
    def@subclasses <- comp
    def
}

## extends() w/o conditional inheritance:  used for S3 inheritance, method
## selection on S4 objects
..extendsForS3 <- function(Class)
    extends(Class, maybe = FALSE)
## dummy version while generating methods package
.extendsForS3 <- function(Class)
    extends(Class)
#  File src/library/methods/R/packageName.R
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

## utilities to manage package names

getPackageName <- function(where = topenv(parent.frame()), create = TRUE) {
    pkg <- ""
    hasNameSaved <- exists(".packageName", where, inherits = FALSE)
    if(hasNameSaved)
        pkg <- get(".packageName", where)
    else  if(identical(where, 1) || identical(as.environment(where), topenv(parent.frame())))
        pkg <- Sys.getenv("R_PACKAGE_NAME")
    env <- as.environment(where)
    envName <- environmentName(env)
    if(nzchar(envName)) {
        if(regexpr("package:", envName, fixed = TRUE) == 1L)
          pkg <- sub("package:","", envName, fixed = TRUE)
    }
    if(!nzchar(pkg)) { ## is still ""
        if(identical(env, .GlobalEnv))
            pkg <- ".GlobalEnv"
        else if(identical(env, .BaseNamespaceEnv))
            pkg <- "base"
        else {
            if(is.numeric(where))
                pkg <- search()[[where]]
            else if(is.environment(where)) {
                for(db in search())
                    if(identical(as.environment(db), where)) {
                        pkg <- db; break
                    }
            }
            else if(nzchar(environmentName(env)))
                pkg <- environmentName(env)
            else
                pkg <- as.character(where)
            if(identical(substr(pkg, 1L, 8L), "package:"))
                pkg <- substr(pkg, 9L, nchar(pkg, "c"))
        }
#  Problem:  the library() function should now be putting .packageName in package environments
#   but namespace makes them invisible from outside.
        ## save the package name, but .GlobalEnv is not a package name,
        ## and package base doesn't have a .packageName (yet?)
#         if(!(identical(pkg, ".GlobalEnv") || identical(pkg, "base")) ) {
#             setPackageName(pkg, env)
#             ## packages OUGHT
#             ## to be self-identifying
#              warning("The package name \"", pkg, "\" was inferred, but not found in that package")
#         }
    }
    if(!nzchar(pkg) && create) {
        pkg <- as.character(Sys.time())
        warning(gettextf("Created a package name, %s, when none found",
                         sQuote(pkg)),
                domain = NA)
        assign(pkg, env, envir = .PackageEnvironments)
        if(!(hasNameSaved || environmentIsLocked(env)))
            setPackageName(pkg, env)
    }
    pkg
}

setPackageName <- function(pkg, env)
    assign(".packageName", pkg, envir = env)

##FIXME:  rather than an attribute, the className should have a formal class
## (but there may be bootstrap problems)
packageSlot <- function(object)
    attr(object, "package")

`packageSlot<-` <- function(object, value) {
    attr(object, "package") <- value
    object
}
#  File src/library/methods/R/promptClass.R
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

promptClass <-
function (clName, filename = NULL, type = "class",
	  keywords = "classes", where = topenv(parent.frame()),
          generatorName = clName)
{
    classInSig <- function(g, where, cl) {
        ## given a generic g, is class cl in one of the method
        ## signatures for the class?
	cl %in% unique(unlist(findMethods(g, where)@signatures))
    }
    genWithClass <- function(cl, where) {
    ## given a class cl
    ## obtain list of all generics with cl in
    ## one of its signatures
	allgen <- getGenerics(where = where)
	ok <- as.logical(unlist(lapply(allgen, classInSig, cl = cl, where = where)))
	allgen[ok]
    }

    sigsList <- function (g, where)
      ## given a generic g, obtain list with one element per signature,
      ## with argument names inserted
    {
        methods <- findMethods(g, where)
	value <- methods@signatures
        args <- methods@arguments
        if(length(value)) {
            ## name the individual signature elements for output
            length(args) <- length(value[[1]]) # all sigs are same length
            value <- lapply(value, function(x){names(x) <- args; x})
        }
        value
    }
    slotClassWithSource <- function(clname) {
	clDef <- getClassDef(clname)
	extds <- names(clDef@contains)
	allslots <- getSlots(clDef) ## establishes all slots, in the right order
	for(j in rev(seq_along(extds))) {
	    i <- extds[[j]]
	    slotsi <- getSlots(getClass(i))
	    if(length(slotsi))
		allslots[names(slotsi)] <- paste0("\"", as.character(slotsi),
						  "\", from class \"", i, "\"")
	}
	slotsi <- getSlots(clDef)
	if(length(slotsi))
	    allslots[names(slotsi)] <- paste0("\"", as.character(slotsi),"\"")
	allslots
    }
    cleanPrompt <- function(object, name) {
        ## get the prompt() result and clean out the junk
        ## lines that prompt() creates
        value <- prompt(object, name = name, filename = NA)
        for(i in seq_along(value)) {
            item <- value[[i]]
            bad <- grepl("^ *%", item)
            if(any(bad))
                value[[i]] <- item[!bad]
        }
        value
    }
    pastePar <- function(x) {
        xn <- names(x)
	x <- as.character(x)
	xn <- if(length(xn) == length(x)) paste(xn, "= ") else ""
	paste0("(", paste0(xn, "\"", x, "\"", collapse = ", "), ")")
    }
    escape <- function(txt) gsub("%", "\\\\%", txt)

    if(is.null(filename))
	filename <- paste0(utils:::topicName(type, clName), ".Rd")
    if(!missing(where) && !is.na(match(clName, getClasses(where))))
      whereClass <- where
    else {
        whereClass <- find(classMetaName(clName))
        if(length(whereClass) == 0L)
            stop(gettextf("no definition of class %s found",
                          dQuote(clName)), domain = NA)
        else if(length(whereClass) > 1L) {
            if(identical(where, topenv(parent.frame()))) {
                whereClass <- whereClass[[1L]]
                warning(gettextf("multiple definitions of %s found; using the one on %s",
                                 dQuote(clName), whereClass), domain = NA)
            }
            else {
                if(exists(classMetaName(clName), where, inherits = FALSE))
                    whereClass <- where
                else
                    stop(sprintf(ngettext(length(whereClass),
                                          "no definition of class %s in the specified position, %s, definition on : %s",
                                          "no definition of class %s in the specified position, %s, definitions on : %s"),
                                 dQuote(clName), where,
                                 paste(whereClass, collapse = ", ")),
                         domain = NA)
            }
        }
    }
    fullName <- utils:::topicName("class", clName)
    clDef <- getClass(clName, where = whereClass)
    .name <- paste0("\\name{", fullName, "}")
    .type <- paste0("\\docType{", type, "}")
    .alias <- paste0("\\alias{", fullName, "}")
    .title <- sprintf("\\title{Class \\code{\"%s\"}}", clName)
    .desc <- paste0("\\description{",
                    "\n%%  ~~ A concise (1-5 lines) description of what the class is. ~~",
                    "\n}")
    slotclasses <- getSlots(clDef)
    slotnames <- names(slotclasses)
    slotclasses <- as.character(slotclasses)
    nslots <- length(slotclasses)
    clNameQ <- paste0('"', clName, '"')
    .usage <- "\\section{Objects from the Class}"
    virtualClass <- isVirtualClass(clName)
    if(virtualClass) {
	.usage <- paste0(.usage, "{A virtual Class: No objects may be created from it.}")
        generator <- NULL # regardless of what exists
    }
    else {
        if(exists(generatorName, where, inherits = FALSE))
            generator <- get(generatorName, where, inherits = FALSE)
        else
            generator <- NULL
        if(is(generator, "classGeneratorFunction")) {
            promptGenerator <- cleanPrompt(generator, generatorName)
            callString <- .makeCallString(generator, generatorName)
            .alias <- c(.alias, promptGenerator$aliases)
            ## the rest of the promptGenerator will be added later
        }
        else {
            initMethod <- unRematchDefinition(selectMethod("initialize", clName))
            argNames <- formalArgs(initMethod)
            ## but for new() the first argument is the class name
            argNames[[1L]] <- clNameQ
            callString <- .makeCallString(initMethod, "new", argNames)
        }
	.usage <-
            c(paste0(.usage,"{"),
              paste0("Objects can be created by calls of the form \\code{",
                     callString,
                     "}."),
              "%%  ~~ describe objects here ~~ ",
              "}")
    }
    .slots <- if (nslots > 0) {
	slotclasses <- slotClassWithSource(clName)
	slotnames <- names(slotclasses)
	.slots.head <- c("\\section{Slots}{", "  \\describe{")
	.slots.body <-	paste0("    \\item{\\code{", slotnames,
                               "}:}", "{Object of class \\code{",
                               slotclasses, "} ~~ }")
	.slots.tail <- c("  }","}")
	c(.slots.head,  .slots.body,	.slots.tail)
    } else character()
    .extends <- clDef@contains
## FIXME: the superclass slots should be marked as such
##       and left *optional* to be documented
    if(length(.extends)) {
	.extends <- showExtends(.extends, printTo = FALSE)
	.extends <-
	    c("\\section{Extends}{",
	      paste0("Class \\code{\"\\linkS4class{",
		    .extends$what,
		    "}\"}, ",
		    ## Add Rd markup to 'by class "CLASS"' results
		    gsub("^(by class) (\".*\")$", "\\1 \\\\code{\\2}",
			 .extends$how),
		    "."),
	      "}")
    }
    else
	.extends <- character()
    nmeths <- length(methnms <- genWithClass(clName, where = whereClass))
    .meths.head <- "\\section{Methods}{"
    .methAliases <- ""
    if (nmeths > 0) {
	.meths.body <- "  \\describe{"
	for (i in 1L:nmeths) {
	    .sig <- sigsList(methnms[i], where = whereClass)
	    for (j in seq_along(.sig)) {
		if (!all(is.na(match(.sig[[j]],clName)))) {
		    methn.i <- escape(methnms[i])
		    .meths.body <-
			c(.meths.body,
			  paste0("    \\item{",
				 methn.i, "}{\\code{signature",
				 pastePar(.sig[[j]]), "}: ... }"))

		    cur <- paste(.sig[[j]], collapse = ",")
		    .methAliases <- paste0(.methAliases, "\\alias{",
					   methn.i, ",", cur, "-method}\n")
		}
	    }
	}
	.meths.body <- c(.meths.body, "	 }")
    }
    else {
	.meths.head <- "\\section{Methods}{"
	.meths.body <- paste("No methods defined with class", clNameQ,
                             "in the signature.")
    }
    .meths.tail <- "}"
    .keywords <- paste0("\\keyword{", keywords, "}")

    Rdtxt <-
	list(name = .name,
             version = "\\Rdversion{1.1}",
	     type = .type,
	     aliases = .alias,
	     methAliases = .methAliases,
	     title = .title,
	     description = .desc,
	     "section{Objects from the Class}" = .usage,
	     "section{Slots}" = .slots,
	     "section{Extends}" = .extends,
	     "section{Methods}" =
	     c(.meths.head, .meths.body, .meths.tail),
	     references = paste("\\references{\n%%  ~~put references to the",
	     "literature/web site here~~\n}"),
	     author = "\\author{\n%%  ~~who you are~~\n}",
	     note =
	     c("\\note{\n%%  ~~further notes~~\n}",
	       "",
	       paste("%% ~Make other sections like Warning with",
		     "\\section{Warning }{....} ~"),
	       ""),
	     seealso =
	     c("\\seealso{",
	       paste("%%  ~~objects to See Also as",
		     "\\code{\\link{~~fun~~}}, ~~~"),
	       paste("%%  ~~or \\code{\\linkS4class{CLASSNAME}}",
		     "for links to other classes ~~~"),
	       "}"),
	     examples = c("\\examples{",
	     paste0("showClass(", clNameQ, ")"),
	     "}"),
	     keywords = .keywords)

    if(is(clDef, "refClassRepresentation"))
        Rdtxt <- refClassPrompt(clDef, Rdtxt, nmeths, nslots, .meths.head)
    else if(is(generator, "classGeneratorFunction")) {
        ## add in the actual usage, arguments sections, mostly to make
        ## CMD check happy
        what <-  c("usage", "arguments")
        Rdtxt[what] <- promptGenerator[what]
    }

    if(is.na(filename)) return(Rdtxt)

    cat(unlist(Rdtxt), file = filename, sep = "\n")
    .message("A shell of class documentation has been written",
             .fileDesc(filename), ".\n")
    invisible(filename)
}

## used in promptClass() above and in promptMethods() :
.fileDesc <- function(file) {
    if(is.character(file)) {
	if(nzchar(file))
	    paste(" to the file", sQuote(file))
	else
	    " to the standard output connection"
    }
    else if(inherits(file, "connection"))
	paste(" to the connection",
              sQuote(summary(file)$description))
    else "" # what, indeed?
}

refClassPrompt <- function(clDef, Rdtxt, nmeths, nslots, .meths.head) {
    ## exclude some sections that are usually irrelevant
    sections <- names(Rdtxt)
    envRefX <- paste0("{",extends("envRefClass"), "}")
    exclude <- grep("Objects from the Class", sections)
    if(nmeths < 1)
        exclude <- c(exclude, grep("Methods", sections))
    else
        .meths.head <- "\\section{Class-Based Methods}{"
    if(nslots < 2) # just the data slot, usually
        exclude <- c(exclude, grep("Slots", sections))
    Rdtxt <- Rdtxt[-exclude]
    extdsthead <- "section{Extends}" # has to be there
    extds <- Rdtxt[[extdsthead]]
    drop <- rep(FALSE, length(extds))
    for(class in envRefX) #drop the envRefClass & its superclasses
        drop <- drop | grepl(class, extds, fixed = TRUE)
    extds <- extds[!drop]
    extds <- append(extds, "\nAll reference classes extend and inherit methods from \\code{\"\\linkS4class{envRefClass}\"}.\n", length(extds)-1)
    Rdtxt[[extdsthead]] <- extds
    fieldClasses <- refClassFields(clDef)
    nfields <- length(fieldClasses)
    .fields <- if (nfields > 0) {
	fieldnames <- names(fieldClasses)
	.fields.head <- c("\\section{Fields}{", "  \\describe{")
	.fields.body <-	paste0("    \\item{\\code{", fieldnames,
                               "}:}", "{Object of class \\code{",
                               fieldClasses, "} ~~ }")
	.fields.tail <- c("  }","}")
	c(.fields.head,  .fields.body,	.fields.tail)
    } else character()
    methodDefs <- as.list(clDef@refMethods)
    nmethods <- length(methodDefs)
    if(nmethods > 0) {
        thisClassDefs <- match(sapply(methodDefs, function(x) x@refClassName), clDef@className, 0) > 0
        otherMethods <- methodDefs[!thisClassDefs]
        methodDefs <- methodDefs[thisClassDefs]
        .methods <-
            c(.meths.head, .refMethodDescription(methodDefs, fieldnames, otherMethods), "}")
    }
    else
        .methods <- character()
    c(Rdtxt,
      list("section{Fields}" = .fields,
           "section{ClassMethods}" = .methods)
      )
}

.refMethodDescription <- function(methodDefs, fieldnames, otherMethods) {
    methodnames <- names(methodDefs)
    methodargs <- sapply(methodDefs, function(x)
			 paste0("(", paste(formalArgs(x), collapse=", "), ")"))
    if(length(methodnames) > 0) {
        .methods.head <- "  \\describe{"
        .methods.body <-
            paste0("    \\item{\\code{",
                   methodnames, methodargs,
                   "}:}", "{ ~~ }")
        .methods <- c(.methods.head,  .methods.body, "  }")
    }
    else
        .methods <- character()
    methodclasses <- sapply(otherMethods,
              function(x) if(is(x, "refMethodDef")) x@refClassName else "<unknown>")
    ## don't report the standard methods from envRefClass
    superclass <- methodclasses != "envRefClass"
    otherMethods <- otherMethods[superclass]
    methodclasses <- methodclasses[superclass]
    if(length(otherMethods)) {
        methodnames <- names(otherMethods)
        methodnames <- gsub("[#].*","", methodnames)
        .methods <- c(.methods,
                      "\nThe following methods are inherited (from the corresponding class):",
                      paste0(methodnames, ' ("', methodclasses,
                             '")', collapse = ", ")
                      )
    }
    .methods
}

.makeCallString <- function (def, name = substitute(def), args = formalArgs(def))
{
##
## need this for experimentation because the function is not exported
##
    if (is.character(def)) {
	if (missing(name))
	    name <- def
	def <- getFunction(def)
    }
    if (is(def, "function"))
	paste0(name, "(", paste(args, collapse = ", "), ")")
    else ""
}
#  File src/library/methods/R/rbind.R
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

#### S4-ized  rbind() --- this is entirely parallel to ./cbind() --- KEEP IN SYNC!
###  -------------------- built by
## s/cbind/rbind/ ; s/nrow/N_COL/; s/column/row/; s/colnam/rownam/;
## s/ncol/nrow/ ; s/N_COL/ncol/; s/d[2L]/d[1L]/

rbind <- function(..., deparse.level = 1)
{
    na <- nargs() - !missing(deparse.level)
    deparse.level <- as.integer(deparse.level)
    stopifnot(0 <= deparse.level, deparse.level <= 2)

    argl <- list(...)
    ## remove trailing 'NULL's:
    while(na > 0 && is.null(argl[[na]])) { argl <- argl[-na]; na <- na - 1 }
    if(na == 0) return(NULL)
    if(na == 1) {
	if(isS4(..1)) return(rbind2(..1))
	else return(.__H__.rbind(..., deparse.level = deparse.level))
    }

    ## else :  na >= 2

    if(deparse.level) {
	symarg <- as.list(sys.call()[-1L])[1L:na] # the unevaluated arguments
	## For default 'deparse.level = 1', rbind(a, b) has to give *names*!
	Nms <- function(i) { # possibly 'deparsed' names of argument  i
	    if(is.null(r <- names(symarg[i])) || r == "") {
		if(is.symbol(r <- symarg[[i]]) || deparse.level == 2)
		    deparse(r)		# else NULL
	    } else r
	}
    }
    if(na == 2) {
	r <- ..2
	fix.na <- FALSE
    }
    else { ## na >= 3 arguments: -- RECURSION -- with care
	## determine ncol(<result>)  for e.g.,	rbind(diag(2), 1, 2)
	## only when the last two argument have *no* dim attribute:
	nrs <- unname(lapply(argl, ncol)) # of length na
	iV <- sapply(nrs, is.null)# is 'vector'
	fix.na <- identical(nrs[(na-1):na], list(NULL,NULL))
	if(fix.na) {
	    ## "fix" last argument, using 1-row `matrix' of proper ncol():
	    nr <- max(if(all(iV)) sapply(argl, length) else unlist(nrs[!iV]))
	    argl[[na]] <- rbind(rep(argl[[na]], length.out = nr),
				deparse.level = 0)
	    ## and since it's a 'matrix' now, rbind() below may not name it
	}
	## need to pass argl, the evaluated arg list to do.call();
	## OTOH, these may have lost their original 'symbols'
	if(deparse.level) {
	    if(fix.na)
		fix.na <- !is.null(Nna <- Nms(na))
	    if(!is.null(nmi <- names(argl))) iV <- iV & (nmi == "")
	    ## attach `symbols' to argl[-1L] for 'vectors'[iV]
	    ii <- if(fix.na) # need to fix later ([na] is 'matrix')
		2:(na-1) else 2:na
	    if(any(iV[ii])) {
		for(i in ii[iV[ii]])
		    if (!is.null(nmi <- Nms(i))) names(argl)[i] <- nmi
	    }
	}
	r <- do.call(rbind, c(argl[-1L], list(deparse.level=deparse.level)))
    }

    d2 <- dim(r)
    r <- rbind2(..1, r)
    if(deparse.level == 0)
	return(r)
    ism1 <- !is.null(d1 <- dim(..1)) && length(d1) == 2L
    ism2 <- !is.null(d2)	     && length(d2) == 2L && !fix.na
    if(ism1 && ism2) ## two matrices
	return(r)

    ## else -- Setting rownames correctly
    ##	       when one was not a matrix [needs some diligence!]
    Nrow <- function(x) {
	d <- dim(x); if(length(d) == 2L) d[1L] else as.integer(length(x) > 0L) }
    nn1 <- !is.null(N1 <- if((l1 <- Nrow(..1)) && !ism1) Nms(1)) # else NULL
    nn2 <- !is.null(N2 <- if(na == 2 && Nrow(..2) && !ism2) Nms(2))
    if(nn1 || nn2 || fix.na) {
	if(is.null(rownames(r)))
	    rownames(r) <- rep.int("", nrow(r))
	setN <- function(i, nams)
	    rownames(r)[i] <<- if(is.null(nams)) "" else nams
	if(nn1) setN(1,	 N1)
	if(nn2) setN(1+l1, N2)
	if(fix.na) setN(nrow(r), Nna)
    }
    r
}
#  File src/library/methods/R/refClass.R
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


## Classes to support OOP-style classes with reference-semantics for fields
## and class-based methods.
## Implementation of the R-based version of these classes (using environments)


envRefInferField <- function(self, field, thisClass, selfEnv = as.environment(self)) {
    'Install a field method into the environment of object
self from reference class thisClass.'
    fields <- thisClass@fieldPrototypes
    if(exists(field, envir = fields, inherits = FALSE)) {
        ## this allows lazy installation of fields (not currently used)
        value <- get(field, envir = fields)
    }
    else {
        methods <- thisClass@refMethods
        if(exists(field, envir = methods, inherits = FALSE)) {
            value <- get(field, envir = methods)
            ## install this method and any methods it may call
            value <- installClassMethod(value, self, field, selfEnv, thisClass)
        }
        else
            stop(gettextf("%s is not a valid field or method name for reference class %s",
                          sQuote(field),
                          dQuote(thisClass@className)),
                 domain = NA)
    }
    value
}

installClassMethod <- function(def, self, me, selfEnv, thisClass) {
    if(!is(def, "refMethodDef")) {  #should not happen? => need warning
        warning(sprintf("method %s from class %s was not processed into a class method until being installed.  Possible corruption of the methods in the class.",
                         me, thisClass@className),
                domain = NA)
        def <- makeClassMethod(def, me, thisClass@className, "", objects(thisClass@refMethods, all.names = TRUE))
        .checkFieldsInMethod(def, names(thisClass@fieldClasses))
        ## cache the analysed method definition
        assign(me, def, envir = thisClass@refMethods)
    }
    depends <- def@mayCall
    environment(def) <- selfEnv # for access to fields and methods
    assign(me, def, envir = selfEnv)
    ## process those that are not in the instance environment, now that
    ## this method has been assigned.
    done <- objects(selfEnv, all.names = TRUE)
    notDone <- depends[is.na(match(depends, done))]
    superCase <- match("callSuper", notDone, 0)
    if(superCase > 0) {
        if(nzchar(def@superClassMethod))
            notDone[[superCase]] <- def@superClassMethod
        else
            stop(gettextf("a call to superClass() is in the method %s but there is no superclass definition of this method for class %s",
                          sQuote(me),
                          dQuote(thisClass@className)),
                 domain = NA)
    }
    for(what in notDone)
        installClassMethod(get(what, envir = thisClass@refMethods), self, what, selfEnv, thisClass)
    if(superCase > 0) {
        ## provide an environment with the correct callSuper() definition,
        ## with selfEnv as its parent (can't override the definition of "callSuper"
        ## in selfEnv--there may  be other methods with a callSuper() in them
        newEnv <- new.env(FALSE, parent = selfEnv)
        assign("callSuper", get(def@superClassMethod, envir = selfEnv),
               envir = newEnv)
        environment(def) <- newEnv
        assign(me, def, envir = selfEnv)
        ## the callSuper() inside def now goes to the right method
    }
    def
   }

..hasCodeTools <- FALSE
.hasCodeTools <- function() {
    if(!identical(..hasCodeTools, TRUE)) # will be FALSE when methods is built, keep checking
        .assignOverBinding("..hasCodeTools",length(list.files(system.file(package = "codetools"))) > 0,
                           .methodsNamespace, FALSE)
    ..hasCodeTools
}

.getGlobalFuns <- function(def) {
    if(.hasCodeTools())
        codetools::findGlobals(def, merge = FALSE)$functions
    else
        unique(unlist(lapply(def, all.names)))
}

makeClassMethod <- function(def, name, Class, superClassMethod = "", allMethods) {
    depends <- .getGlobalFuns(def)
    ## find the field methods called ...
    if("usingMethods" %in% depends) { # including those declared
        declared <- .declaredMethods(def)
        ## look for invalid declared methods
        if(length(declared) && any(! declared %in% allMethods))
            warning(gettextf("methods declared in usingMethods() but not found: %s",
                paste0(declared[! declared %in% allMethods], collapse = ", ")))
        depends <- c(declared, depends)
    }
    depends <- depends[match(depends, allMethods, 0) > 0]
    new("refMethodDef", def, mayCall = depends, name = name,
        refClassName = Class, superClassMethod = superClassMethod)
}

refObjectClass <- function(object) {
    Class <- class(object)
    classDef <- getClassDef(Class)
    if(is(classDef, "refClassRepresentation"))
        classDef
    else
        stop(gettextf("%s is not a reference class",
                      dQuote(Class)),
             domain = NA)
}

envRefSetField <- function(object, field,
                           thisClass = refObjectClass(object),
                           env = as.environment(object), value) {
    fieldClass <- thisClass@fieldClasses[[field]]
    if(is.null(fieldClass))
        stop(gettextf("%s is not a field in class %s",
                      sQuote(field),
                      dQuote(thisClass@className)),
             domain = NA)
    else
        assign(field, value, envir = env)
    object
}

.initForEnvRefClass <- function(.Object, ...) {
    Class <- class(.Object)
    classDef <- getClass(Class)
    selfEnv <- new.env(TRUE, .NamespaceOrPackage(classDef@package))
    ## the parent environment will be used by field methods, to make
    ## them consistent with functions in this class's package
    .Object@.xData <- selfEnv
    ## install prototypes and active bindings
    prototypes <- classDef@fieldPrototypes
    fieldClasses <- classDef@fieldClasses
    fields <- names(fieldClasses)
    for(field in fields) {
        fp <- prototypes[[field]] # prototype or NULL
        if(is(fp, "activeBindingFunction")) {
            environment(fp) <- selfEnv
            makeActiveBinding(field, fp, selfEnv)
            if(is(fp, "defaultBindingFunction")) {
                ## ensure an initial value
                class <- fieldClasses[[field]]
                if(isVirtualClass(class))
                    value <- NULL
                else
                    value <- new(class)
                assign(.bindingMetaName(field), value, envir = selfEnv)
            }
        }
        else
            assign(field, fp, envir = selfEnv)
    }
    ## assign references to the object and to its class definition
    selfEnv$.self <- .Object
    selfEnv$.refClassDef <- classDef
    if(is.function(classDef@refMethods$initialize)) {
        .Object$initialize(...)
        ## intialize methods are allowed to change .self
        .Object <- selfEnv$.self
    }
    else {
        if(nargs() > 1) {
            .Object <-
                methods::initRefFields(.Object, classDef, selfEnv, list(...))
        }
    }
    if(is.function(classDef@refMethods$finalize))
        reg.finalizer(selfEnv, function(x) x$.self$finalize())
    lockBinding(".self", selfEnv)
    lockBinding(".refClassDef", selfEnv)
    .Object
}

## old version, for back compatibility.  Could be deleted after 2.15.0
initFieldArgs <- function(.Object, classDef, selfEnv, ...)
    initRefFields(.Object, classDef, selfEnv, list(...))

initRefFields <- function(.Object, classDef, selfEnv, args) {
    if(length(args)) {
        fieldDefs <- classDef@fieldClasses
        fieldNames <- names(fieldDefs)
        snames <- allNames(args)
        which <- nzchar(snames)
        elements <- args[which]
        supers <- args[!which]
        elNames <- names(elements)
        for(super in supers) {
            if(!is(super, "refClass")) {
                warning(gettextf("unnamed arguments to $new() must be objects from a reference class; got an object of class %s",
                                 dQuote(class(super))),
                        domain = NA)
                next
            }
            fields <- names(super$.refClassDef@fieldClasses)
            ##<FIXME> need an object$fields for the above </FIXME>
            ## assign field if it is not already specified
            fields <- fields[is.na(match(fields, elNames))]
            for(field in fields)
                elements[[field]] <- super$field(field)
            elNames <- names(elements)
        }
        ## assign the fields
        for(field in elNames)
            envRefSetField(.Object, field, classDef, selfEnv, elements[[field]])
    }
    .Object
}

.dollarForEnvRefClass <- function(x, name) {
    what <- substitute(name)
    if(is.symbol(what))
        what <- as.character(what)
    else
        what <- name
    selfEnv <- as.environment(x)
    if(exists(what, envir = selfEnv, inherits = FALSE))
        ## either a field or previously cached method
        get(what, envir = selfEnv)
    else if(is(x, "envRefClass"))
        ## infer (usually) the method, cache it and return it
        envRefInferField(x, what, getClass(class(x)), selfEnv)
    else # don't know the reference class(e.g., x is the refMethods env.)
        stop(gettextf("%s is not a valid field or method name for this class",
                      sQuote(what)),
             domain = NA)
}

.dollarGetsForEnvRefClass <- function(x, name, value) {
    what <- substitute(name)
    if(is.symbol(what))
        what <- as.character(what)
    else
        what <- name
    selfEnv <- as.environment(x)
    envRefSetField(x, what, refObjectClass(x), selfEnv, value)
    invisible(x)
}

.envRefMethods <-
    list(
         export = function(Class) {
             '
Returns the result of coercing the object to
Class.  No effect on the object itself.
'
             if(match(.refClassDef@className, Class, 0) > 0)
                 return(.self)
             classDef <- getClass(Class)
             if(is(classDef, "refClassRepresentation") &&
                !is.na(match(Class, .refClassDef@refSuperClasses))) {
                 value <- new(classDef)
                 env <- as.environment(value)
                 selfEnv <- as.environment(.self)
                 fieldClasses <- classDef@fieldClasses
                 for(field in names(fieldClasses)) {
                     current <- get(field, envir = selfEnv)
                     if(!is(current, fieldClasses[[field]]))
                         stop(gettextf("the class of field %s in the object is not compatible with the desired class %s in the target",
                                       sQuote(field),
                                       dQuote(fieldClasses[[field]])),
                              domain = NA)
                     assign(field, envir = env, current)
                 }
                 value
             }
             else if(is(classDef, "classRepresentation")) # use standard S4 as()
                 methods::as(.self, Class)
             else if(is.character(Class) && length(Class) == 1)
                 stop(gettextf("%s is not a defined class in this environment",
                               dQuote(Class)),
                      domain = NA)
             else
                 stop("invalid 'Class' argument:  should be a single string")
         },
         import =   function(value, Class = class(value)) {
             '
Imports value, replacing the part of the current object
corresponding to Class (if argument Class is missing
it is taken to be class(value)).  The Class must be one
of the reference superclasses of the current class (or
that class itself, but then you could just overrwite the object).
'
             if(!missing(Class))
                 value <- value$export(Class)
             classDef <- getClass(Class)
             if(is(classDef, "refClassRepresentation") &&
                (!is.na(match(Class, .refClassDef@refSuperClasses))
                || identical(classDef@className, .refClassDef@className))) {
                 env <- as.environment(value)
                 selfEnv <- as.environment(.self)
                 fieldClasses <- .refClassDef@fieldClasses
                 for(field in names(classDef@fieldClasses)) {
                     current <- get(field, envir = env)
                     if(!is(current, fieldClasses[[field]]))
                         stop(gettextf("the class of field %s in the object is not compatible with the desired class %s in the target",
                                       sQuote(field),
                                       dQuote(fieldClasses[[field]])),
                              domain = NA)
                     assign(field, envir = selfEnv, current)
                 }
                 invisible(.self)
             }
             else
                 stop(gettextf("%s is not one of the reference super classes for this object",
                               dQuote(Class)),
                      domain = NA)
         },
         callSuper = function(...) stop("direct calls to callSuper() are invalid:  should only be called from another method"),
         initFields = function(...) {
             if(missing(...)) .self else
             initRefFields(.self, .refClassDef, as.environment(.self), list(...))
         },
         copy = function(shallow = FALSE) {
             def <- .refClassDef
             value <- new(def)
             vEnv <- as.environment(value)
             selfEnv <- as.environment(.self)
             for(field in names(def@fieldClasses)) {
                 if(shallow)
                     assign(field, get(field, envir = selfEnv), envir = vEnv)
                 else {
                     current <- get(field, envir = selfEnv)
                     if(is(current, "envRefClass"))
                         current <- current$copy(FALSE)
                     assign(field, current, envir = vEnv)
                 }
             }
             value
         },
         getRefClass = function(Class = .refClassDef) methods::getRefClass(Class),
         getClass = function(...) if(nargs()) methods::getClass(...) else .refClassDef,
         field = function(name, value) if(missing(value)) get(name, envir = .self) else {
             if(is.na(match(name, names(.refClassDef@fieldClasses))))
                 stop(gettextf("%s is not a field in this class",
                               sQuote(name)),
                      domain = NA)
             assign(name, value, envir = .self)
         },
         trace = function(..., classMethod = FALSE) {
             ' Insert trace debugging for the specified method.  The arguments are
 the same as for the trace() function in package "base".  The first argument
 should be the name of the method to be traced, quoted or not.

 The additional argument classMethod= can be supplied as TRUE (by name only)
 in order to trace a method in a generator object (e.g., "new") rather than
 in the objects generated from that class.
'
             .TraceWithMethods(..., where = .self, classMethod = classMethod)
         },
         untrace = function(..., classMethod = FALSE) {
             ' Untrace the method given as the first argument.
'
             .TraceWithMethods(..., untrace = TRUE,  where = .self, classMethod = classMethod)
         },
         show = function() {
             cat('Reference class object of class ', classLabel(class(.self)),
        '\n', sep = "")
             fields <- names(.refClassDef@fieldClasses)
             for(fi in fields) {
                 cat('Field "', fi, '":\n', sep = "")
                 methods::show(field(fi))
             }
         },
         usingMethods = function(...) {
             ' Reference methods used by this method are named as the arguments
 either quoted or unquoted.  In the code analysis phase of installing the
 the present method, the declared methods will be included.  It is essntial
 to declare any methods used in a nonstandard way (e.g., via an apply function).
 Methods called directly do not need to be declared, but it is harmless to do so.
 $usingMethods() does nothing at run time.
'
             NULL
         }
         )

## construct a list of class methods for envRefClass
makeEnvRefMethods <- function() {
    methods <- .envRefMethods
    allMethods <- names(methods)
    for(method in allMethods) {
        methods[[method]] <- makeClassMethod(methods[[method]],
                   method, "envRefClass", "", allMethods)
    }
    methods
}

## initialize some reference classes
.InitRefClasses <- function(envir)
{
    ## class to define a reference class
    ## Should be split into an abstract class and a standard version
    ## to use environments, so other variants might use interfaces
    ## to OOP languages, and proxy objects

    setClass("refClassRepresentation",
             representation(fieldClasses = "list",
                            fieldPrototypes = "environment",
                            refMethods = "environment",
                            refSuperClasses = "character"),
             contains = "classRepresentation", where = envir)
    ## the virtual class from which all true reference clases
    ## inherit.  Its subclasses require methods
    ## for getting & setting fields and related tasks
    setClassUnion("refClass", where = envir)
    ## the union of all reference objects
    ## (including those not belonging to refClass)
    setClassUnion("refObject", c("environment", "externalptr", "name",                                "refClass"), where = envir)
    ## a class for field methods, with a slot for their dependencies,
    ## allowing installation of all required instance methods
    setClassUnion("SuperClassMethod", "character")
    ## helper classes for active binding of fields
    setClass("activeBindingFunction", contains = "function")
    setClass("defaultBindingFunction",
             representation(field = "character", className = "character"),
             contains = "activeBindingFunction")
    ## class to mark uninitialized fields
    setClass("uninitializedField",
             representation(field = "character", className = "character"))
    setClass("refMethodDef",
             representation(mayCall = "character", name = "character",
                            refClassName = "character",
                            superClassMethod = "SuperClassMethod"),
             contains = "function", where = envir)
    ## and make a traceable version of the class
    .makeTraceClass(.traceClassName("refMethodDef"), "refMethodDef", FALSE)
    setIs("refMethodDef", "SuperClassMethod", where = envir)
    setClass("envRefClass", contains = c("environment","refClass"), where =envir)
    ## bootstrap envRefClass as a refClass
    def <- new("refClassRepresentation",
               refMethods = as.environment(makeEnvRefMethods()))
    as(def, "classRepresentation") <- getClassDef("envRefClass", where = envir)
    assignClassDef("envRefClass", def, where = envir)
    setMethod("initialize", "envRefClass", methods:::.initForEnvRefClass,
              where = envir)
    ## NOTE:  "$" method requires setting in methods:::.InitStructureMethods
    setMethod("$", "envRefClass", .dollarForEnvRefClass, where = envir)
    setMethod("$<-", "envRefClass", .dollarGetsForEnvRefClass, where = envir)
    setMethod("show", "envRefClass", function(object) object$show())
    setClass("refGeneratorSlot") # a temporary virtual class to allow the next definition
    ## the refClassGenerator class
    setClass("refObjectGenerator", representation(generator ="refGeneratorSlot"),
             contains = c("classGeneratorFunction", "refClass"), where = envir)

    setMethod("$", "refObjectGenerator",
              function(x, name) eval.parent(substitute(x@generator$name)), where = envir)

    setMethod("$<-", "refObjectGenerator",
              function(x, name, value) eval.parent(substitute(x@generator$name <- value)),
              where = envir)
    ## next call is touchy:  setRefClass() uses an object of class
    ## refGeneratorSlot, but the class should have been defined before
    ## that object is created.
    setRefClass("refGeneratorSlot",
                fields = list(def = "ANY", className = "ANY"),
                methods = .GeneratorMethods, where = envir)
    setMethod("show", "refClassRepresentation",
              function(object) showRefClassDef(object), where = envir)
    setMethod("show", "refObjectGenerator",
              function(object) showRefClassDef(object$def, "Generator for class"),
              where = envir)
    setMethod("show", "refMethodDef", showClassMethod, where = envir)
    ## Now do "localRefClass"; doesn't need to be precisely here
    ## but this ensures it is not done too early or too late
    setRefClass("localRefClass", methods = .localRefMethods,
                where = envir)  # should this have contains = "VIRTUAL"?

    setMethod("$<-", "localRefClass",
              function(x, name, value) {
                  w <- parent.frame()
                  x <- .ensureLocal(x, w)
                  what <- substitute(name)
                  if (is.symbol(what))
                      what <- as.character(what)
                  else what <- name
                  selfEnv <- as.environment(x)
                  envRefSetField(x, what, refObjectClass(x), selfEnv, value)
                  invisible(x)
              } , where = envir)
}

getRefSuperClasses <- function(classes, classDefs) {
    supers <- character()
    for(i in seq_along(classes)) {
        clDef <- classDefs[[i]]
        supers <- c(supers, clDef@refSuperClasses)
    }
    unique(supers)
}

.GeneratorMethods <- list(methods =  function(...) {
    methodsEnv <- def@refMethods
    if(nargs() == 0)
        return(objects(methodsEnv, all.names = TRUE))
    if(methods:::.classDefIsLocked(def))
        stop(gettextf("the definition of class %s in package %s is locked, methods may not be redefined",
                      dQuote(def@className),
                      sQuote(def@package)),
             domain = NA)
    methodDefs <- list(...)
    ## allow either name=function, ... or a single list
    if(length(methodDefs) == 1 && is.list(methodDefs[[1]]))
        methodDefs <- methodDefs[[1]]
    mnames <- names(methodDefs)
    if(is.null(mnames) || !all(nzchar(mnames)))
        stop("arguments to methods() must be named, or one named list")
    ## look for methods to remove (new definition is NULL)
    removeThese <- sapply(methodDefs, is.null)
    if(any(removeThese)) {
        rmNames <- mnames[removeThese]
        mnames <- mnames[!removeThese]
        methodDefs <- methodDefs[!removeThese]
        remove(list = rmNames, envir = methodsEnv)
        if(length(mnames) == 0)
            return(invisible(methodsEnv))
    }
    allMethods <- as.list(methodsEnv)
    ## get a list of processed methods, plus any
    ## overriden superclass methods
    newMethods <- insertClassMethods(allMethods, className, methodDefs, names(def@fieldClasses), FALSE)
    for(what in names(newMethods))
        assign(what, newMethods[[what]], envir = methodsEnv)
    ## calls to $methods() only work in package source or
    ## as load actions.  Use the topenv() if that seems like
    ## the namespace in preparation, or the namespace if available
    env <- topenv(parent.frame()); declare <- TRUE
    if(exists(".packageName", envir = env) &&
       get(".packageName", envir = env) == def@package) {}
    else if(def@package %in% loadedNamespaces())
        env <- asNamespace(def@package)
    else
        declare <- FALSE
    if(declare)
        utils::globalVariables(names(newMethods), env)
    invisible(methodsEnv)
},

fields =  function() {
    '
Returns the named vector of classes
for the fields in this class.  Fields
defined with accessor functions have
class "activeBindingFunction".
'
    unlist(def@fieldClasses)
},
new =  function(...) {
    methods::new(def, ...)
},
  help =  function(topic) {
    '
Prints simple documentation for the method or field
specified by argument topic, which should be the name
of the method or field, quoted or not.  With no topic,
prints the definition of the class.
'
    if(missing(topic)) {
        writeLines(
c('Usage:  $help(topic) where topic is the name of a method (quoted or not)',
  paste('The definition of class', className, 'follows.')))
        methods::show(def)
    }
    else {
        if(is.name(substitute(topic)))
            topic <- as.character(substitute(topic))
        else
            topic <- as.character(topic)
        env <- def@refMethods
        if(exists(topic, envir = env)) {
            writeLines(.refMethodDoc(topic, env))
        }
        else {
            cat(gettextf("topic %s is not a method name in class %s\nThe class definition follows\n",
                         sQuote(topic),
                         dQuote(className)))
            show(def)
        }
    }
},
lock =  function(...) methods:::.lockRefFields(def, ...),
## define accessor functions, store them in the refMethods environment
## of the class definition.
accessors = function(...) {
    firstCap <- function(names) {
        firstChars <- substr(names, 1,1)
        modChars <- toupper(firstChars)
        substr(names, 1, 1) <- modChars
        list(get = paste0("get", names), set = paste0("set", names))
    }
    if(methods:::.classDefIsLocked(def))
        stop(gettextf("the definition of class %s in package %s is locked so fields may not be modified",
                      dQuote(def@className),
                      sQuote(def@package)),
             domain = NA)
    fieldNames <- c(...)
    methodNames <- firstCap(fieldNames)
    getters <- methodNames$get
    setters <- methodNames$set
    accessors <- list()
    for(i in seq_along(fieldNames)) {
        what <- fieldNames[[i]]
        field <- as.name(what)
        CLASS <- def@fieldClasses[[what]]
        if(is.null(CLASS))
            stop(gettextf("%s is not a field in class %s",
                          sQuote(what),
                          dQuote(def@className)),
                 domain = NA)
        accessors[[getters[[i]] ]] <-
                     eval(substitute(function() X, list(X = field)))
        if(CLASS == "ANY")
            accessors[[setters[[i]] ]] <-
                eval(substitute(function(value) {
                    value <- as(value, CLASS, strict = FALSE)
                    X <<- value
                    invisible(value)
                    },
                                list(X = field, CLASS = CLASS)))
        else
            accessors[[setters[[i]] ]] <-
                eval(substitute(function(value) {
                    X <<- value
                    invisible(value)
                    },
                                list(X = field)))
    }
    ## install the accessors
    methods(accessors)
    invisible(accessors)
}
)

.localRefMethods <-
    list(
         ensureLocal = function() {
             'Ensure that a shallow copy has been made of this object
to localize any further changes.  Must be called before any reference
class method modifies a field.
'
             methods:::.ensureLocal(.self, parent.frame())
         }
     )

.makeCall <- function(name, x) {
    n <- length(argls <- formals(x))
    noDeflt <- if(n > 0) sapply(argls,function(x)  !is.name(x) || nzchar(as.character(x)))
    if (n) {
        arg.names <- arg.n <- names(argls)
    }
    Call <- paste0("$", name, "(")
    for (i in seq_len(n)) {
        Call <- paste0(Call, arg.names[i], if (noDeflt[[i]]) " = "
            )
        if (i != n)
            Call <- paste0(Call, ", ")
    }
    Call <- paste0(Call, ")\n")
    Call
}


`insertFields<-` <- function(fieldList, value) {
    newNames <- names(value)
    ## check for valid overrides of existing field definitions
    hasFields <- match(newNames, names(fieldList),0) > 0
    if(any(hasFields)) {
        for(field in newNames[hasFields])
            ## the new field class must be a subclass of the old
            if(is.na(match(fieldList[[field]], c(extends(value[[field]]),"ANY"))))
                stop(gettextf("the overriding class (\"%s\") of field %s is not a subclass of the existing field definition (\"%s\")",
                              value[[field]],
                              sQuote(field),
                              fieldList[[field]]),
                     domain = NA)
    }
    fieldList[newNames] <- value
    fieldList
}

.bindingMetaName <- function(fieldName)
    paste0(".->", fieldName)

.makeActiveBinding <- function(thisField) {
    if(is(thisField, "activeBindingFunction"))
     thisField
    else
     new("activeBindingFunction", thisField)
}

.makeDefaultBinding <- function(fieldName, fieldClass, readOnly = FALSE, where) {
    metaName <- .bindingMetaName(fieldName)
    if(readOnly)
        ## write-once into the metaName object
        f <-  eval(substitute(function(value) {
            if(missing(value))
                dummyFieldName
            else {
                methods:::.setDummyField(.self, dummyField, dummyClass, thisField, TRUE, value)
                value
            }
        }, list(dummyField = metaName, thisField = fieldName,
                dummyClass = fieldClass, dummyFieldName = as.name(metaName))))
    else
        f <- eval(substitute(function(value) {
            if(missing(value))
                dummyFieldName
            else {
                methods:::.setDummyField(.self, dummyField, dummyClass, thisField, FALSE, value)
                value
            }
        }, list(dummyField = metaName, dummyClass = fieldClass,
                thisField = fieldName, dummyFieldName = as.name(metaName))))
    environment(f) <- where ## <note> Does this matter? </note>
    f <- new("defaultBindingFunction", f,
             field = fieldName, className = fieldClass)
    init <- (if(isVirtualClass(fieldClass))
                 new("uninitializedField", field = fieldName,
                     className = fieldClass)
        else new(fieldClass))
    value <- list(f, init)
    names(value) <- c(fieldName, metaName)
    value
}

.setDummyField <- function(self, metaName, fieldClass, fieldName, onceOnly, value) {
    if(is(value, fieldClass))
        value <- as(value, fieldClass, strict = FALSE) # could be more efficient?
    else
        stop(gettextf("invalid assignment for reference class field %s, should be from class %s or a subclass (was class %s)",
                       sQuote(fieldName), dQuote(fieldClass), dQuote(class(value))), call. = FALSE)
    selfEnv <- as.environment(self)
    if(onceOnly) {
        if(bindingIsLocked(metaName, selfEnv))
            stop(gettextf("invalid replacement: reference class field %s is read-only", sQuote(fieldName)),
                 call. = FALSE)
        else {
            assign(metaName, value, envir = selfEnv)
            lockBinding(metaName, selfEnv)
        }
    }
    else
       assign(metaName, value, envir = selfEnv)
}

refClassInformation <- function(Class, contains, fields, refMethods, where) {
    if(length(contains) > 0) {
        superClassDefs <- lapply(contains,
                                 function(what) {
                                     if(is(what, "classRepresentation"))
                                         what
                                     else if(is.character(what))
                                         getClass(what, where = where)
                                     else
                                         stop(gettextf("the 'contains' argument should be the names of superclasses:  got an element of class %s",
                                                       dQuote(class(what))),
                                              domain = NA)
                                 })
        missingDefs <- sapply(superClassDefs, is.null)
        if(any(missingDefs))
            stop(gettextf("no definition found for inherited class: %s",
                          paste0('"',contains[missingDefs], '"', collapse = ", ")),
                 domain = NA)
        superClasses <- unlist(lapply(superClassDefs,
                          function(def) def@className), FALSE)
        isRefSuperClass <- sapply(superClassDefs, function(def)
                              is(def, "refClassRepresentation"))
    }
    else {
        superClassDefs <- list()
        superClasses <- character()
        isRefSuperClass <- logical()
    }
    if(!any(isRefSuperClass)) {
        superClasses <- c(superClasses, "envRefClass")
        isRefSuperClass <- c(isRefSuperClass, TRUE)
        superClassDefs[["envRefClass"]] <- getClass("envRefClass", where = where)
    }
    refSuperClasses <- superClasses[isRefSuperClass]
    otherRefClasses <- getRefSuperClasses(refSuperClasses, superClassDefs[isRefSuperClass])
    refSuperClasses <- unique(c(refSuperClasses, otherRefClasses))
    ## process the field definitions.  The call from setRefClass
    ## guarantees that fields is a named list.
    fieldNames <- names(fields)
    nf <- length(fields)
    fieldClasses <- character(nf)
    names(fieldClasses) <- fieldNames
    fieldPrototypes <- list()
    for(i in seq_len(nf)) {
        thisName <- fieldNames[[i]]
        thisField <- fields[[i]]
        ## a field definition can be:
        ## 1. character string name of the class
        ## 2. a binding function
        if(is.character(thisField)) {
            if(length(thisField) != 1)
                stop(gettextf("a single class name is needed for field %s, got a character vector of length %d",
                              sQuote(thisName),
                              length(thisField)),
                     domain = NA)
            if(is.null(getClassDef(thisField, where = where)))
                stop(gettextf("class %s for field %s is not defined",
                              dQuote(thisField),
                              sQuote(thisName)),
                     domain = NA)
            fieldClasses[[i]] <- thisField
            if(thisField != "ANY")
                fieldPrototypes <- c(fieldPrototypes,
                    .makeDefaultBinding(thisName, thisField, where = where))
            else
                fieldPrototypes[[thisName]] <-
    new("uninitializedField", field = thisName,
                        className = "ANY")
        }
        else if(is.function(thisField)) {
            fieldClasses[[i]] <- "activeBindingFunction"
            fieldPrototypes[[thisName]] <-
                .makeActiveBinding(thisField)
        }
        else
            stop(gettextf("field %s was supplied as an object of class %s; must be a class name or a binding function",
                          sQuote(thisName),
                          dQuote(class(thisField))),
                 domain = NA)
    }
    ## assemble inherited information
    fc <- fp <- cm <- list(); fr <- character()
    ## assign in reverse order so nearer superclass overrides
    for(cl in rev(superClassDefs[isRefSuperClass])) {
        fcl <- cl@fieldClasses
        fpl <- as.list(cl@fieldPrototypes, all.names = TRUE) # turn env into list
        cml <- as.list(cl@refMethods, all.names = TRUE) # ditto
        insertFields(fc) <- fcl
        fp[names(fpl)] <- fpl
        cm[names(cml)] <- cml
    }
    insertFields(fc) <- fieldClasses
    fp[names(fieldPrototypes)] <- fieldPrototypes

    ## process and insert reference methods
    cm <- insertClassMethods(cm, Class, refMethods, names(fc), TRUE)
    list(superClasses = superClasses, refSuperClasses = refSuperClasses,
         fieldClasses = fc, fieldPrototypes = fp,
         refMethods = cm)
}

superClassMethodName <- function(def)
    paste(def@name, def@refClassName, sep = "#")

insertClassMethods <- function(methods, Class, value, fieldNames, returnAll) {
    ## process reference methods, return either the entire updated methods
    ## or the processed new methods in value, plus superclass versions
    theseMethods <- names(value)
    prevMethods <- names(methods) # catch refs to inherited methods as well
    allMethods <- unique(c(theseMethods, prevMethods))
    if(returnAll)
        returnMethods <- methods
    else
        returnMethods <- value
    check <- TRUE
    for(method in theseMethods) {
        prevMethod <- methods[[method]] # NULL or superClass method
        if(is.null(prevMethod)) {
            ## kludge because default version of $initialize() breaks bootstrapping of methods package
            if(identical(method, "initialize"))
                superClassMethod <- "initFields"
            else
                superClassMethod <- ""
        }
        else if(identical(prevMethod@refClassName, Class))
            superClassMethod <- prevMethod@superClassMethod
        else {
            superClassMethod <- superClassMethodName(prevMethod)
            returnMethods[[superClassMethod]] <- prevMethod
        }
        def <- makeClassMethod(value[[method]], method, Class,
                               superClassMethod, allMethods)
        check <- check && .checkFieldsInMethod(def, fieldNames, allMethods)
        returnMethods[[method]] <- def
    }
    if(is.na(check) && .methodsIsLoaded())
        message(gettextf("code for methods in class %s was not checked for suspicious field assignments (recommended package %s not available?)",
                         dQuote(Class),
                         sQuote("codetools"))
                , domain = NA)
    returnMethods
}


## refField <- function(class = "ANY", get = .stdGetField, set = .stdSetField, binding = NULL,
##                      name = "", where = topenv(parent.frame())) {
##     if(identical(set, FALSE))
##         set <- .invalidSetField
##     new("refFieldDefinition",  fieldName = name, fieldClass = class,
##         get = get, set = set, binding = binding)
##   }

setRefClass <- function(Class, fields = character(),
                        contains = character(),
                        methods = list(),
                        where = topenv(parent.frame()),
                        ...) {
    fields <- inferProperties(fields, "field")
    theseMethods <- names(methods) # non-inherited, for processing later
    ## collect the method and field definitions
    info <- refClassInformation(Class, contains, fields, methods, where)
    ## make codetools happy:
    superClasses <- refSuperClasses <- fieldClasses <- fieldPrototypes <-
        refMethods <- NULL
    ## think Python's multiple assignment operator
    for(what in c("superClasses", "refSuperClasses", "fieldClasses",
                  "fieldPrototypes", "refMethods"))
        assign(what, info[[what]])
    ## temporarily assign an ordinary class definition
    ## to allow the checks and defaults from setClass to be applied
    ## and to get the classGeneratorFunction
    ## Note:  the classGeneratorFunction has the class name, not the explicit definition
    classFun <- setClass(Class, contains = superClasses,
             where = where, ...)
    ## kludge: as.environment fails on an empty list
    asEnv <- function(x) {
        if(length(x)) as.environment(x) else new.env(FALSE)
    }
    ## now, override the class definiton with the complete definition
    classDef <- new("refClassRepresentation",
                    getClassDef(Class, where = where),
                    fieldClasses = fieldClasses,
                    refMethods = asEnv(refMethods),
                    fieldPrototypes = asEnv(fieldPrototypes),
                    refSuperClasses = refSuperClasses)
    assignClassDef(Class, classDef, where)
    generator <- new("refGeneratorSlot")
    env <- as.environment(generator)
    env$def <- classDef
    env$className <- Class
    .declareVariables(classDef, where)
    value <- new("refObjectGenerator", classFun, generator = generator)
    invisible(value)
}

getRefClass <- function(Class, where = topenv(parent.frame())) {
    if(is(Class, "refClassRepresentation")) {
        classDef <- Class
        Class <- classDef@className
    }
    else if(is.character(Class)) {
        classDef <- getClass(Class, where = where)
        if(!is(classDef, "refClassRepresentation"))
            stop(gettextf("class %s is defined but is not a reference class",
                          dQuote(Class)),
                 domain = NA)
    }
    else
        stop(gettextf("class must be a reference class representation or a character string; got an object of class %s",
                      dQuote(class(Class))),
             domain = NA)
    generator <- new("refGeneratorSlot")
    env <- as.environment(generator)
    env$className <- Class
    env$def <- classDef
    classFun <- classGeneratorFunction(Class, where)
    ## but, the package is always from the class definition, not the local environment
    classFun@package <- classDef@package
    new("refObjectGenerator", classFun, generator = generator)
}

refClassFields <- function(Class) {
    ClassDef <- getClass(Class)
    if(is(ClassDef, "refClassRepresentation"))
        ClassDef@fieldClasses
    else
        stop(gettextf("not a reference class: %s", ClassDef@name),
             domain = NA)
}

refClassMethods <- function(Class) {
    ClassDef <- getClass(Class)
    if(is(ClassDef, "refClassRepresentation"))
        value <- as.list(ClassDef@refMethods)
    else
        stop(gettextf("not a reference class: %s", ClassDef@name),
             domain = NA)
    ## possibly temporary:  return methods to pure functions
    for(i in seq_along(value))
        value[[i]] <- as(value[[i]], "function")
    value
}

showClassMethod <- function(object) {
    cl <- class(object)
    cat("Class method definition")
    if(!.identC(cl, "refMethodDef"))
        cat(sprintf(" (class %s)", dQuote(cl)))
    cat(sprintf(" for method %s()\n", object@name))
    show(as(object, "function"))
    if(length(object@mayCall))
        .printNames("Methods used: ", object@mayCall)
}

.printNames <- function(header, names, separateLine = TRUE) {
    if(separateLine)
        cat("\n",header,"\n    ")
    else
        cat(header,": ",sep="")
    cat(paste0('"', names, '"'), sep = ", ", fill = TRUE)
    cat("\n")
    }

showRefClassDef <- function(object, title = "Reference Class") {
    cat(title," \"", object@className,"\":\n", sep="")
    fields <- object@fieldClasses
    if(length(fields)) {
        printPropertiesList(fields, "Class fields")
        locked <- .getLockedFieldNames(object)
        if(length(locked))
            .printNames("Locked Fields", locked, FALSE)
    }
    else
        cat("\nNo fields defined\n")
    methods <- objects(object@refMethods, all.names = TRUE)
    if(length(methods))
        .printNames("Class Methods: ", methods)
    else
        cat ("\nNo Class Methods\n")
    supers <- object@refSuperClasses
    if(length(supers))
        .printNames("Reference Superclasses: ", supers)
}


## all.equal and identical both screw up on environments
## but a bigger change is needed to all.equal than the following
## because it also screws up on, e.g., externalptr objects
if(FALSE) {
all.equal.environment <- function(target, current, ...) {
    nt <- sort(objects(target, all.names = TRUE))
    nc <- sort(objects(current, all.names = TRUE))
    tmp <- all.equal(nt, nc, ...)
    if(!identical(tmp, TRUE))
        return(paste("Different objects in target, current:", tmp))
    if(length(nt) == 0)
        return(TRUE)
    differ <- sapply(nt, function(what) {
        tmp <- all.equal(get(what, envir = target),
                         get(what, envir = current), ...)
        if(identical(tmp, TRUE)) FALSE
        else TRUE
    })
    if(any(differ))
        paste("Objects differ: ", paste(nt[differ], collapse = ", "))
    else
        TRUE
}
}

.assignExpr <- function(e) {
    value <- list()
    value[[codetools::getAssignedVar(e)]] <- deparse(e, nlines = 1L)
    value
}

.mergeAssigns <- function(previous, new) {
    for(what in names(new)) {
        if(is.null(previous[[what]]))
            previous[[what]] <- new[[what]]
        else
            previous[[what]] <- paste(previous[[what]], new[[what]], sep="; ")
    }
    previous
}


.assignedVars <- function(e) {
    locals <- list()
    globals <- list()
    walker <- codetools::makeCodeWalker(call = function(e, w) {
        callto <- e[[1]]
        if(is.symbol(callto)) switch(as.character(callto),
               "<-" = , "=" = {
                   locals <<- .mergeAssigns(locals, .assignExpr(e))
               },
               "<<-" = {
                   globals <<- .mergeAssigns(globals, .assignExpr(e))
               })
        for (ee in as.list(e))
            if (! missing(ee)) codetools::walkCode(ee, w)
    },
    leaf = function(e, w) NULL
    )
    codetools::walkCode(e, walker)
    list(locals = locals, globals = globals)
}

.checkFieldsInMethod <- function(methodDef, fieldNames, methodNames) {
    if(!.hasCodeTools())
        return(NA)
    p0q <- function(x) paste0('"', x, '"', collapse = "; ")
    if(is(methodDef, "refMethodDef")) {
        methodName <- p0q(methodDef@name)
        className <- p0q(methodDef@refClassName)
    }
    else {
        methodName <- className <- ""
    }
    assigned <- .assignedVars(body(methodDef))
    locals <- names(assigned$locals)
    localsAreFields <- match(locals, fieldNames, 0) > 0
    if(any(localsAreFields))
        warning(gettextf("local assignment to field name will not change the field:\n    %s\n Did you mean to use \"<<-\"? ( in method %s for class %s)",
                paste(unlist(assigned$locals)[localsAreFields], collapse="; "), methodName, className),
                domain = NA)
    globals <- names(assigned$globals)
    ## check non-fields, but allow to .self (will be an
    ## error except in $initialize())
    globalsNotFields <- is.na(match(globals, c(fieldNames, ".self")))
    if(any(globalsNotFields))
        warning(gettextf("non-local assignment to non-field names (possibly misspelled?)\n    %s\n( in method %s for class %s)",
                paste(unlist(assigned$globals)[globalsNotFields], collapse="; "), methodName, className),
                domain = NA)
    globalsInMethods <- match(globals, methodNames, 0) > 0
    if(any(globalsInMethods))
        stop(gettextf("non-local assignment to method names is not allowed\n    %s\n( in method %s for class %s)",
                paste(unlist(assigned$globals)[globalsInMethods], collapse="; "), methodName, className),
                domain = NA)
    !any(localsAreFields) && !any(globalsNotFields)
}

.refMethodDoc <- function(topic, env) {
    f <- get(topic, envir = env)
    msg <- c("Call:",.makeCall(topic, f), "")
    bb <- body(f)
    ## look for self-documentation
    if(is(bb, "{") && length(bb) > 1 && is(bb[[2]], "character"))
        msg <- c(msg, bb[[2]], "")
    msg
}

## the locked fields are stored as a hidden object in the fieldPrototypes environment
## but this might change, so the .get, .set functions should be used
.lockedFieldsMetaName <- ".#lockedFields"
.getLockedFieldNames <- function(def) {
    env <- def@fieldPrototypes
    value <- env[[.lockedFieldsMetaName]]
    if(is.null(value))
        character()
    else
        value
}
.setLockedFieldNames <- function(def, value) {
    env <- def@fieldPrototypes
    env[[.lockedFieldsMetaName]] <- value
    value
}

.lockRefFields <- function(def, ...) {
    lockedFields <- .getLockedFieldNames(def)
    if(nargs()<2)
        return(lockedFields)
    fields <- c(...)
    if(is.character(fields) && all(nzchar(fields))) {}
    else
        stop("arguments must all be character string names of fields")
    if(.classDefIsLocked(def))
        stop(gettextf("the definition of class %s in package %s is locked so fields may not be modified",
                      dQuote(def@className),
                      sQuote(def@package)),
             domain = NA)
    env <- def@fieldPrototypes
    className <- def@className
    for(what in fields) {
        if(what %in% lockedFields) {
            warning(gettextf("field %s is already locked", sQuote(what)),
                    domain = NA)
            next
        }
        current <- env[[what]]
        if(is.null(current))
            stop(gettextf("%s is not a field in class %s",
                          sQuote(what),
                          dQuote(className)),
                 domain = NA)
        if(is(current, "activeBindingFunction")) {
            if(is(current, "defaultBindingFunction"))
                env[[what]] <- .makeDefaultBinding(current@field,
                    current@className, TRUE, environment(current))[[what]]
            else
                stop(gettextf("field %s of class %s has a non-default binding and cannot be locked",
                              sQuote(what),
                              dQuote(className)),
                     domain = NA)
        }
        else {
            ## capture the current prototype value with a read-only binding function
            binding <- .makeDefaultBinding(current@field,
               current@className, TRUE, environment(current))
            env[[what]] <- binding[[what]]
            metaName <- .bindingMetaName(what)
            env[[metaName]] <- current
        }
        lockedFields <- c(lockedFields, what)
    }
    .setLockedFieldNames(def, lockedFields)
    invisible(env)
}

## declare field and method names global to avoid spurious
## messages from codetools
.declareVariables <- function(def, env) {
    utils::globalVariables(c(names(def@fieldClasses), objects(def@refMethods)),
                           env)
}

.declaredMethods <- function(method) {
    methods <- character()
    if(!.hasCodeTools())
        return(methods)
    .theseMethods <- function(e, w) {
        if(length(e) < 2) character()
        else
            sapply(as.list(e)[-1], function(what)
                   methods <<- c(methods, if(is.symbol(what)) as.character(what) else if(is.character(what)) what else character()))
    }
    walker <- codetools::makeCodeWalker(
                handler = function(v, w) {
                    if(identical(v, "usingMethods"))
                        .theseMethods
                    else
                        NULL
                },
                leaf = function(e, w) NULL)
    codetools::walkCode(body(method), walker)
    unique(methods)
}

getMethodsAndAccessors <- function(Class) {
    def <- getClass(Class)
    if(!is(def, "refClassRepresentation"))
        stop(gettextf("%s is not a reference class",
             dQuote(def@className)))
    ff <- def@fieldPrototypes
    accs <- sapply(ff, function(what) is(what, "activeBindingFunction") && !is(what, "defaultBindingFunction"))
    c(as.list(def@refMethods), as.list(ff)[accs])
}

## Reference classes that guarantee to change fields only in the
## local environment.  The method for `$<-` checks that the lhs object
## has been registered in a list of local reference class objects in
## the frame where the call is evaluated.  If not, a shallow copy
## of the object's .self (environment) is made, replaces the variable
## and is registered.  The effect should be that locality of assignment
## is preserved wtihout the deep copy generated by the R evaluator
## for complex assignments that are not primitives, e.g., `@<-`

.ensureLocal <- function(object, where) {
    if(!is(object, "envRefClass"))
        stop(gettextf("Class %s is not a subclass of %s; functional semantics not defined for this class", dQuote(class(object)), dQuote("envRefClass")))
    selfEnv <- as.environment(object)
    if(exists(".localRefObjects", envir = where, inherits = FALSE)) {
        locals <- get(".localRefObjects", envir = where)
        for(i in rev(seq_along(locals)))
            if(identical(selfEnv, locals[[i]]))
                return(object)
    }
    else
        locals <- list()
    ## the object should be assigned in environment where=
    what <- NULL
    for(objName in objects(envir = where, all.names = TRUE)) {
        obj <- get(objName, envir = where)
        if(is(obj, "envRefClass") && identical(selfEnv, as.environment(obj))) {
            what <- obj
            break
        }
    }
    if(is.null(what))
        stop("Could not find local object in supplied environment")
    ## do a shallow copy and record it as local
    value <- .shallowCopy(object, selfEnv)
    locals[[length(locals)+1]] <- as.environment(value)
    assign(".localRefObjects", locals, envir = where)
    value
}

## a shallow copy of a reference object
## This code depends on knowledge of how classes extend "environment"
.shallowCopy <- function(object, selfEnv) {
    newEnv <- new.env()
    for(what in objects(envir = selfEnv, all.names = TRUE))
        assign(what, get(what, envir = selfEnv), envir = newEnv)
    attr(object, ".xData") <- newEnv
    assign(".self", object, envir = newEnv)
    object
}
#  File src/library/methods/R/show.R
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

showDefault <- function(object, oldMethods = TRUE)
{
    clDef <- getClass(cl <- class(object), .Force=TRUE)
    cl <- classLabel(cl)
    if(!is.null(clDef) && isS4(object) && is.na(match(clDef@className, .BasicClasses)) ) {
        cat("An object of class ", cl, "\n", sep="")
        slots <- slotNames(clDef)
        dataSlot <- .dataSlot(slots)
        if(length(dataSlot) > 0) {
            dataPart <- slot(object, dataSlot)
            show(dataPart)
            slots <- slots[is.na(match(slots, dataSlot))]
        }
        else if(length(slots) == 0L)
            show(unclass(object))
        for(what in slots) {
            if(identical(what, ".Data"))
                next ## should have been done above
            cat("Slot \"",what, "\":\n", sep="")
            print(slot(object, what))
            cat("\n")
        }
    }
##     else if(isS4(object) && isClass(clDef) && extends(clDef, "oldClass") &&
##             length(slotNames(clDef)) > 0) {
##         ## print the old-style object
##         cat("An object of class ", cl, "\n", sep="")
##         slots <- slotNames(clDef)
##         i <- match(".S3Class", slots)
##         if(is.na(i)) { } # but should not happen with new objects
##         else {
##             S3Class <- classLabel(object@.S3Class)
##             slots <- slots[! slots %in% names(slotsFromS3(object))]
##             if(!identical(cl, S3Class)) {
##                 if(length(S3Class) > 1)
##                   cat("  (S3 class: c(", paste0('"', S3Class, '"', collapse = ", "), "))\n", sep="")
##                 else
##                   cat("  (S3 class: \"",S3Class, "\")\n", sep = "")
##             }
##         }
##         for( cl2 in rev(extends(clDef)))
##             if(!.identC(cl2, "oldClass") && extends(cl2, "oldClass")) {
##                 print(as(object, cl2), useS4 = FALSE) # see comment NBB below
##                 break
##             }
##         for(what in slots) {
##             cat("Slot \"",what, "\":\n", sep="")
##             print(slot(object, what))
##             cat("\n")
##         }
##     }
    else
        ## NBB:  This relies on the delicate fact (as of version 1.7 at least)
        ## that print will NOT recursively call show if it gets more than one argument!
        print(object, useS4 = FALSE)
    invisible() # documented return for show().
}

.extraSlotsDone <- new.env() # any unique reference value would do

showExtraSlots <- function(object, ignore) {
    if(is(ignore, "classRepresentation"))
      ignore <- slotNames(ignore)
    else if(!is(ignore, "character"))
      stop(gettextf("invalid 'ignore' argument; should be a class definition or a character vector, got an object of class %s", dQuote(class(ignore))),
           domain = NA)
    slots <- slotNames(class(object))
    for(s in slots[is.na(match(slots, ignore))]) {
        cat("Slot ",s, ":\n", sep="")
        show(slot(object, s))
    }
    .extraSlotsDone # a signal not to call this function (again)
}

## temporary definition of show, to become the default method
## when .InitShowMethods is called
show <- function(object)
    showDefault(object, FALSE)

.InitShowMethods <- function(envir) {
    if(!isGeneric("show", envir))
        setGeneric("show", where = envir, simpleInheritanceOnly = TRUE)
    setMethod("show", "MethodDefinition",
              function(object) {
                  cl <- class(object)
		  nonStandard <-
		      if(.identC(cl, "MethodDefinition"))
			  "" else paste0(" (Class ", classLabel(cl),")")
                  cat("Method Definition",nonStandard,":\n\n", sep = "")
                  show(object@.Data)
                  mm <- methodSignatureMatrix(object)
                  cat("\nSignatures:\n")
                  print(mm)
              },
              where = envir)
    setMethod("show", "MethodWithNext",
              function(object)  {
                  callNextMethod()
                  cat("\nExcluded from nextMethod:\n")
                  print(unlist(object@excluded))
              },
              where = envir)
    setMethod("show", "genericFunction",
              function(object)  {
                  cat(class(object)," for \"", object@generic,
                      "\" defined from package \"", object@package,
                      "\"\n", sep = "")
                  if(length(object@group))
                      cat("  belonging to group(s):",
                          paste(unlist(object@group), collapse =", "), "\n")
                  if(length(object@valueClass))
                      cat("  defined with value class: \"", object@valueClass,
                          "\"\n", sep="")
                  cat("\n")
                  show(object@.Data)
                  cat("Methods may be defined for arguments: ",
                      paste(object@signature, collapse=", "), "\n",
			    "Use  showMethods(\"", object@generic,
			    "\")  for currently available ones.\n", sep="")
                  if(.simpleInheritanceGeneric(object))
                      cat("(This generic function excludes non-simple inheritance; see ?setIs)\n");
              },
              where = envir)
    setMethod("show", "classRepresentation",
              function(object){
                  if(!.identC(class(object), "classRepresentation"))
                    cat("Extended class definition (", classLabel(class(object)),
                        ")\n")
                  printClassRepresentation(object)
              },
              where = envir)

    ## a show() method for the signature class
    setMethod("show", "signature", function(object) {
        message(gettextf("An object of class %s", dQuote(class(object))),
                domain = NA)
        val <- object@.Data
        names(val) <- object@names
        callNextMethod(val)
    } ,
              where = envir)
}

.showPackage <- function(className) {
    if(is.logical(opt <- getOption("showPackageForClass")))
        opt
    else
        is.list(.Call(C_R_getClassFromCache, as.character(className), .classTable))
}
## an informative string label for a class
classLabel <- function(Class) {
    if(is.character(Class) && length(Class)) {
        className <- Class[[1L]]
        packageName <- attr(Class, "package")
        if(is.null(packageName))
            packageName <- ""
    }
    else {
        if(is(Class, "classRepresentation")) {
            className <- Class@className
            packageName <- Class@package
        }
        else stop(gettextf("invalid call to 'classLabel': expected a name or a class definition, got an object of class %s", classLabel(class(Class))), domain = NA)
    }
    if(.showPackage(className)) {
	packageName <-
	    if(identical(packageName, ".GlobalEnv"))
		" (from the global environment)"
	    else
		paste0(" (from package \"", packageName, "\")")
       paste0('"', className, '"', packageName)
   }
   else
       paste0('"', className, '"')
}
#  File src/library/methods/R/substituteDirect.R
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

substituteDirect <-
  ## subsitute the for the variables named in the second argument the corresponding
  ## objects, substituting into `object'.
  ##
  ## This function differs from the ordinary `substitute' in that it treats its first argument
  ## in the standard S way, by evaluating it.  In contrast, `substitute' does
  ## not evaluate its first argument.
  function(object, frame = parent.frame(), cleanFunction = TRUE)
{
    value <- .Call(C_do_substitute_direct, object, frame)
     if(cleanFunction && is.function(value)) {
       ## unset any local environment
       environment(value) <- .GlobalEnv
     }
    value
  }

#  File src/library/methods/R/trace.R
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

## some temporary (!) hooks to trace the tracing code
.doTraceTrace <- function(on) {
 .assignOverBinding(".traceTraceState", on,
                    environment(.doTraceTrace), FALSE)
  on
}

.traceTraceState <- FALSE

## the internal functions in the evaluator.  These are all prohibited,
## although some of them could just barely be accomodated, with some
## specially designed new definitions (not using ..., for example).
## The gain does not seem worth the inconsistencies; and "if" can
## never be traced, since it has to be used to determine if tracing is
## on.  (see .doTrace())
## The remaining invalid functions create miscellaneous bugs, maybe
## related to the use of "..." as the introduced arguments.  Aside from
## .Call, tracing them seems of marginal value.

.InvalidTracedFunctions <- c("if", "where", "for", "repeat", "(", "{",
                            "next", "break", ".Call", ".Internal", ".Primitive")

.TraceWithMethods <- function(what, tracer = NULL, exit = NULL, at =
                              numeric(), print = TRUE, signature =
                              NULL, where = .GlobalEnv, edit = FALSE,
                              from = NULL, untrace = FALSE, classMethod = FALSE) {
    if(is.function(where)) {
        ## start from the function's environment:  important for
        ## tracing from a namespace
        if(is(where, "genericFunction"))
            where <- parent.env(environment(where))
        else
            where <- environment(where)
        fromPackage <- getPackageName(where)
    }
    else fromPackage <- ""
    doEdit <- !identical(edit, FALSE)
    whereF <- NULL
    pname <- character()
    def <- NULL
    tracingWhere <- "in package"
    refCase <- isS4(where) && (is(where, "envRefClass") ||
                       is(where, "refClassRepresentation"))
    if(refCase) {
        ## some error checking
        if(!is.null(signature))
            stop("argument 'signature' is not meaningful for tracing reference methods")
        .where <- where # to avoid substituting where in the eval() below
        ## A reference class object or its class or its generator
        if(is(.where, "refGeneratorSlot") && !classMethod)
            .where <- .where$def # should now be the refClassRepresentation
        if(is(.where, "refClassRepresentation")) {
            pname <- .where@className
            .where <- .where@refMethods
            tracingWhere <- "for class"
        }
        else {
            tracingWhere <- "for object from class"
            pname <- class(.where)
        }
        ## interpret as tracing .where$what
        def <- eval(substitute(.dollarForEnvRefClass(.where, what)))
        if(!is(def, "refMethodDef")) {
            thisName <- substitute(what)
            stop(gettextf("%s is not a method for reference class %s",
                          sQuote(as.character(if(is.symbol(thisName)) thisName else what)),
                          dQuote(class(where))),
                 domain = NA)
        }
        what <- def@name
        whereF <- .where
    }
    else if(is.function(what)) {
        def <- what
        if(is(def, "genericFunction")) {
            what <- def@generic
            whereF <- .genEnv(what, where)
            pname <- def@package
        }
        else {
            fname <- substitute(what)
            if(is.name(fname)) {
                what <- as.character(fname)
                temp <- .findFunEnvAndName(what, where)
                whereF <- temp$whereF
                pname <- temp$pname
            }
            else if(is.call(fname) && identical(fname[[1L]], as.name("::"))) {
                whereF <- as.character(fname[[2L]])
                require(whereF, character.only = TRUE)
                whereF <- as.environment(paste("package", whereF, sep=":"))
                pname <-  fname[[2L]]
                what <- as.character(fname[[3L]])
            }
            else if(is.call(fname) && identical(fname[[1L]], as.name(":::"))) {
                pname <- paste(fname[[2L]], "(not-exported)")
                whereF <- loadNamespace(as.character(fname[[2L]]))
                what <- as.character(fname[[3L]])
            }
            else
                stop("argument 'what' should be the name of a function")
        }
    }
    else {
        what <- as(what, "character")
        if(length(what) != 1) {
            for(f in what) {
                if(nargs() == 1)
                    trace(f)
                else
                    Recall(f, tracer, exit, at, print, signature, where, edit, from, untrace)
            }
            return(what)
        }
        temp <- .findFunEnvAndName(what, where, signature)
        whereF <- temp$whereF
        pname <- temp$pname
    }
    if(what %in% .InvalidTracedFunctions)
        stop(gettextf("tracing the internal function %s is not allowed",
                      sQuote(what)))
    if(.traceTraceState) {
        message(".TraceWithMethods: after computing what, whereF", domain = NA)
        browser()
    }
    if(nargs() == 1)
        return(.primTrace(what)) # for back compatibility
    if(is.null(whereF)) {
        allWhere <- findFunction(what, where = where)
        if(length(allWhere)==0)
            stop(gettextf("no function definition for %s found",
                          sQuote(what)),
                 domain = NA)
        whereF <- as.environment(allWhere[[1L]])
    }
    ## detect use with no action specified (old-style R trace())
    if(is.null(tracer) && is.null(exit) && identical(edit, FALSE))
        tracer <- quote({})
    if(is.null(def))
        def <- getFunction(what, where = whereF)
    if(is(def, "traceable") && identical(edit, FALSE) && !untrace)
        def <- .untracedFunction(def)
    if(!is.null(signature)) {
        fdef <- if(is.primitive(def))  getGeneric(what, TRUE, where) else def
        def <- selectMethod(what, signature, fdef = fdef, optional = TRUE)
        if(is.null(def)) {
            warning(gettextf("cannot untrace method for %s; no method defined for this signature: %s",
                             sQuote(what),
                             paste(signature, collapse = ", ")),
                    domain = NA)
            return(def)
        }
        ## pick up signature with package slot from selectMethod
        signature <- def@target
    }
    if(untrace) {
        if(.traceTraceState) {
            message(".TraceWithMethods: untrace case", domain = NA)
            browser()
        }

        if(is.null(signature)) {
            ## ensure that the version to assign is untraced
            if(is(def, "traceable")) {
                newFun <- .untracedFunction(def)
            }
            else {
                .primUntrace(what) # to be safe--no way to know if it's traced or not
                return(what)
            }
        }
        else {
            if(is(def, "traceable"))
                newFun <- .untracedFunction(def)
            else {
                warning(gettextf("the method for %s for this signature was not being traced",
                                 sQuote(what)),
                        domain = NA)
                return(what)
            }
        }
    }
    else {
        if(!is.null(exit)) {
            if(is.function(exit)) {
                tname <- substitute(exit)
                if(is.name(tname))
                    exit <- tname
                exit <- substitute(TRACE(), list(TRACE=exit))
            }
        }
        if(!is.null(tracer)) {
            if(is.function(tracer)) {
                tname <- substitute(tracer)
                if(is.name(tname))
                    tracer <- tname
                tracer <- substitute(TRACE(), list(TRACE=tracer))
            }
        }
        original <- .untracedFunction(def)
        traceClass <- .traceClassName(class(original))
        if(is.null(getClassDef(traceClass)))
            traceClass <- .makeTraceClass(traceClass, class(original))
        if(doEdit && is.environment(edit)) {
            ## trace with the version found in the edit environment
            def <- .findNewDefForTrace(what, signature, edit, fromPackage)
            environment(def) <- environment(original)
            if(is.null(c(tracer, exit))) {
                newFun <- new(traceClass, original)
                newFun@.Data <- def
            }
            else {
                newFun <- new(traceClass, def = def, tracer = tracer, exit = exit, at = at, print = print, doEdit = FALSE)
                newFun@original <- original # left as def by initialize method
            }
            newFun@source <- edit
        }
        else
            newFun <- new(traceClass,
                      def = if(doEdit) def else original, tracer = tracer, exit = exit, at = at,
                      print = print, doEdit = edit)
    }
    global <- identical(whereF, .GlobalEnv)
    if(.traceTraceState) {
        message(".TraceWithMethods: about to assign or setMethod", domain = NA)
        browser()
    }
    if(is.null(signature)) {
        if(bindingIsLocked(what, whereF))
            .assignOverBinding(what, newFun, whereF, global)
        else
            assign(what, newFun, whereF)
        if(length(grep("[^.]+[.][^.]+", what)) > 0) { #possible S3 method
            ## check for a registered version of the object
            S3MTableName <- ".__S3MethodsTable__."
            tracedFun <- get(what, envir = whereF, inherits = TRUE)
            if(exists(S3MTableName, envir = whereF, inherits = FALSE)) {
                tbl <- get(S3MTableName, envir = whereF, inherits = FALSE)
                if(exists(what, envir = tbl, inherits = FALSE))
                    assign(what, tracedFun, envir = tbl)
            }
        }
    }
    else {
        if(untrace && is(newFun, "MethodDefinition") &&
           !identical(newFun@target, newFun@defined))
            ## we promoted an inherited method for tracing, now we have
            ## to remove that method.  Assertion is that there was no directly
            ## specified method, or else defined, target would be identical
            newFun <- NULL
        ## arrange for setMethod to put the new method in the generic
        ## but NOT to assign the methods list object (binding is ignored)
        setMethod(fdef, signature, newFun, where = baseenv())
    }
    if(!global) {
        action <- if(untrace)"Untracing" else "Tracing"
        nameSpaceCase <- FALSE
        location <- if(.identC(fromPackage, "")) {
            if(length(pname)==0  && !is.null(whereF))
                pname <- getPackageName(whereF)
            nameSpaceCase <- isNamespace(whereF) &&
            !is.na(match(pname, loadedNamespaces())) &&
            identical(whereF, getNamespace(pname))
            if(length(pname)==0)  # but not possible from getPackagename ?
                ""
            else {
                if(nameSpaceCase)
                    paste0(" in environment <namespace:",  pname, ">")
                else
                    paste0(" ", tracingWhere, " \"",  pname, "\"")
            }
        }
        else paste0(" as seen from package \"", fromPackage, "\"")
        object <- if(refCase) "reference method" else if(is.null(signature)) "function" else "specified method for function"
        object <- paste0(" ", object, " \"", what, "\" ")
        .message(action, object, location)
        if(nameSpaceCase && !untrace && exists(what, envir = .GlobalEnv)) {
            untcall<- paste("untrace(\"", what, "\", where = getNamespace(\"",
                            pname, "\"))", sep="")
            .message("Warning: Tracing only in the namespace; to untrace you will need:\n    ",untcall, "\n")
        }
    }
    what
}

.makeTracedFunction <- function(def, tracer, exit, at, print, doEdit) {
    switch(typeof(def),
           builtin = , special = {
               fBody <- substitute({.prim <- DEF; .prim(...)},
                                   list(DEF = def))
               def <- eval(function(...)NULL)
               body(def, envir = .GlobalEnv) <- fBody
               warning("making a traced version of a primitive; arguments will be treated as '...'")
           }
           )
    if(!identical(doEdit, FALSE)) {
        if(is.character(doEdit) || is.function(doEdit)) {
            editor <- doEdit
            doEdit <- TRUE
        }
        else
            editor <- getOption("editor")
    }
    ## look for a request to edit the definition
    if(doEdit) {
        if(is(def, "traceable"))
            def <- as(def, "function") # retain previous tracing if editing
        if(is(editor, "character") && !is.na(match(editor, c("emacs","xemacs")))) {
            ## cater to the usual emacs modes for editing R functions
            file <- tempfile("emacs")
            file <- sub('..$', ".R", file)
        }
        else
            file <- ""
        ## insert any requested automatic tracing expressions before editing
        if(!(is.null(tracer) && is.null(exit) && length(at)==0))
            def <- Recall(def, tracer, exit, at, print, FALSE)
        def2 <- utils::edit(def, editor = editor, file = file)
        if(!is.function(def2))
            stop(gettextf("the editing in trace() can only change the body of the function; got an object of class %s",
                          dQuote(class(def2))),
                 domain = NA)
        if(!identical(args(def), args(def2)))
            stop("the editing in trace() can only change the body of the function, not the arguments or defaults")
        fBody <- body(def2)
    }
    else {
        def <- .untracedFunction(def) # throw away earlier tracing
        fBody <- body(def)
        if(length(at) > 0) {
            if(is.null(tracer))
                stop("cannot use 'at' argument without a trace expression")
            else if(!inherits(fBody, "{"))
                stop("cannot use 'at' argument unless the function body has the form '{ ... }'")
            for(i in at) {
		fBody[[i]] <-
		    if(print)
			substitute({.doTrace(TRACE, MSG); EXPR},
                                   list(TRACE = tracer,
                                        MSG = paste("step",paste(i, collapse=",")),
                                        EXPR = fBody[[i]]))
		    else
			substitute({.doTrace(TRACE); EXPR},
                                   list(TRACE=tracer, EXPR = fBody[[i]]))
            }
        }
        else if(!is.null(tracer)){
	    fBody <-
		if(print)
		    substitute({.doTrace(TRACE, MSG); EXPR},
                                    list(TRACE = tracer, MSG = paste("on entry"), EXPR = fBody))
		else
		    substitute({.doTrace(TRACE); EXPR},
                                    list(TRACE=tracer, EXPR = fBody))
        }
        if(!is.null(exit)) {
	    exit <-
		if(print)
		    substitute(.doTrace(EXPR, MSG),
                                   list(EXPR = exit, MSG = paste("on exit")))
		else
		    substitute(.doTrace(EXPR),
                                   list(EXPR = exit))
            fBody <- substitute({on.exit(TRACE); BODY},
                                list(TRACE=exit, BODY=fBody))
        }
    }
    body(def, envir = environment(def)) <- fBody
    def
}

## return the untraced version of f
.untracedFunction <- function(f) {
    while(is(f, "traceable"))
        f <- f@original
    f
}


.InitTraceFunctions <- function(envir)  {
    setClass("traceable", representation(original = "PossibleMethod", source = "environment"), contains = "VIRTUAL",
             where = envir); clList <- "traceable"
    ## create the traceable classes
    for(cl in c("function", "MethodDefinition", "MethodWithNext", "genericFunction",
                "standardGeneric", "nonstandardGeneric", "groupGenericFunction",
                "derivedDefaultMethod")) {
        .makeTraceClass(.traceClassName(cl), cl, FALSE)
        clList <- c(clList, .traceClassName(cl))
    }
    setClass("sourceEnvironment", contains = "environment",
             representation(packageName = "character", dateCreated = "POSIXt", sourceFile = "character"),
             prototype = prototype( packageName = "", dateCreated = Sys.time(), sourceFile = ""))
    clList <- c(clList, "sourceEnvironment")
    assign(".SealedClasses", c(get(".SealedClasses", envir), clList), envir)
    setMethod("initialize", "traceable",
              function(.Object, ...) .initTraceable(.Object, ...),
              where = envir)
    if(!isGeneric("show", envir))
        setGeneric("show", where = envir, simpleInheritanceOnly = TRUE)
    setMethod("show", "traceable", .showTraceable, where = envir)
    setMethod("show", "sourceEnvironment", .showSource, where = envir)
}

## allow control over whether methods & classes are cached when assigned
## to a particular environment. defaults to TRUE
cacheOnAssign <- function(env) is.null(env$.cacheOnAssign) || env$.cacheOnAssign
setCacheOnAssign <- function(env, onOff = cacheOnAssign(env))
    env$.cacheOnAssign <- if(onOff) TRUE else FALSE

.showTraceable <- function(object) {
    if(identical(object@source, emptyenv())) {
        cat("Object with tracing code, class \"", class(object),
        "\"\nOriginal definition: \n", sep="")
        callGeneric(object@original)
        cat("\n## (to see the tracing code, look at body(object))\n")
    }
    else {
        cat("Object of class \"", class(object),
            "\", from source\n", sep = "")
        callGeneric(object@.Data)
        cat("\n## (to see original from package, look at object@original)\n")
    }
}

.initTraceable <- function(.Object, def, tracer, exit, at, print, doEdit) {
    .Object@source <- emptyenv()
    if(missing(def))
        return(.Object)
    oldClass <- class(def)
    oldClassDef <- getClass(oldClass)
    if(!is.null(oldClassDef) && length(oldClassDef@slots) > 0)
        as(.Object, oldClass) <- def # to get other slots in def
    .Object@original <- def
    if(nargs() > 2) {
        if(!is.null(elNamed(getSlots(getClass(class(def))), ".Data")))
          def <- def@.Data
        .Object@.Data <- .makeTracedFunction(def, tracer, exit, at, print, doEdit)
    }
    .Object
}

.showSource <- function(object) {
    cat("Object of class \"", class(object), "\"\n", sep = "")
    cat("Source environment created ", format(object@dateCreated), "\n")
    if(nzchar(object@packageName))
        cat("For package \"",object@packageName, "\"\n", sep = "")
    if(nzchar(object@sourceFile))
        cat("From source file \"", object@sourceFile, "\"\n", sep = "")
}

.doTracePrint <- function(msg = "") {
    call <- deparse(sys.call(sys.parent(1)))
    if(length(call)>1)
        call <- paste(call[[1L]], "....")
    cat("Tracing", call, msg, "\n")
}

.traceClassName <- function(className) {
    className[] <- paste0(className, "WithTrace")
    className
}

.assignOverBinding <- function(what, value, where, verbose = TRUE) {
    if(verbose) {
        pname <- getPackageName(where)
        msg <-
            gettextf("assigning over the binding of symbol %s in environment/package %s",
                     sQuote(what), sQuote(pname))
        message(strwrap(msg), domain = NA)
    }
    warnOpt <- options(warn= -1) # kill the obsolete warning from R_LockBinding
    on.exit(options(warnOpt))
    if(is.function(value)) {
        ## assign in the namespace for the function as well
        fenv <- environment(value)
        if(is.null(fenv)) # primitives
          fenv <- baseenv()
        if(!identical(fenv, where) && exists(what, envir = fenv, inherits = FALSE #?
                                             ) && bindingIsLocked(what, fenv)) {
            unlockBinding(what, fenv)
            assign(what, value, fenv)
            lockBinding(what, fenv)
        }
    }
    if(exists(what, envir = where, inherits = FALSE) && bindingIsLocked(what, where)) {
      unlockBinding(what, where)
      assign(what, value, where)
      lockBinding(what, where)
    }
    else
      assign(what, value, where)
}

.setMethodOverBinding <- function(what, signature, method, where, verbose = TRUE) {
    if(verbose)
        warning(gettextf("setting a method over the binding of symbol %s in environment/package %s",
                         sQuote(what),
                         sQuote(getPackageName(where))),
                domain = NA)
    if(exists(what, envir = where, inherits = FALSE)) {
        fdef <- get(what, envir = where)
        hasFunction <- is(fdef, "genericFunction")
    }

        hasFunction <- FALSE
    if(hasFunction) {
        ## find the generic in the corresponding namespace
        where2 <- findFunction(what, where = environment(fdef))[[1L]] # must find it?
        unlockBinding(what, where)
        setMethod(what, signature, method, where = where)
        lockBinding(what, where)
        ## assign in the package namespace as well
        unlockBinding(what, where2)
        setMethod(what, signature, method, where = where2)
        lockBinding(what, where2)
    }
    else {
        setMethod(what, signature, method, where = where)
    }
}

### finding the package name for a loaded namespace
.searchNamespaceNames <- function(env)
    paste("namespace", getNamespaceName(env), sep=":")

.findFunEnvAndName <- function(what, where, signature = NULL) {
    pname <- character()
    if(is.null(signature)) {
        whereF <- findFunction(what, where = where)
        if(length(whereF)>0)
            whereF <- whereF[[1L]]
        else return(list(pname = pname, whereF = baseenv()))
    } else
        whereF <- .genEnv(what, where)

    ## avoid partial matches to "names"
    if("name" %in% names(attributes(whereF)))
        pname <- gsub("^.*:", "", attr(whereF, "name"))
    else if(isNamespace(whereF))
        pname <- .searchNamespaceNames(whereF)
    list(pname = pname, whereF = whereF)
}

.makeTraceClass <- function(traceClassName, className, verbose = TRUE) {
  ## called because the traceClassName not a class
  ## first check whether it may exist but not in the same package
  if(isClass(as.character(traceClassName)))
    return(as.character(traceClassName))
  if(verbose)
    message(sprintf("Constructing traceable class %s", dQuote(traceClassName)),
            domain = NA)
  env <- .classEnv(className)
  if(environmentIsLocked(env)) {
    message(gettextf("Environment of class %s is locked; using global environment for new class",
                     dQuote(className)),
            domain = NA)
    env <- .GlobalEnv
    packageSlot(traceClassName) <- NULL
  }
  setClass(traceClassName,
                 contains = c(className, "traceable"), where = env)
  if(existsMethod("show", className, env)) # override it for traceClassName
    setMethod("show", traceClassName, .showTraceable)
  traceClassName
}


utils::globalVariables("fdef")
.dummySetMethod <- function(f, signature = character(), definition,
	     where = topenv(parent.frame()), valueClass = NULL,
	     sealed = FALSE)
{
    if(is.function(f) && is(f, "genericFunction"))
        f <- fdef@generic
    else if(is.function(f)) {
        if(is.primitive(f))
            f <- .primname(f)
        else
            stop("a function for argument 'f' must be a generic function")
    } else
        f <- switch(f, "as.double" = "as.numeric", f)
    assign(.dummyMethodName(f, signature), definition, envir = where)
}

.functionsOverriden <- c("setClass", "setClassUnion", "setGeneric", "setIs", "setMethod", "setValidity")

.setEnvForSource <- function(env) {
    doNothing <- function(x, ...)x
    ## establish some dummy definitions & a special setMethod()
    for(f in .functionsOverriden)
        assign(f, switch(f, setMethod = .dummySetMethod, doNothing),
               envir = env)
    env
}

.dummyMethodName <- function(f, signature)
    paste(c(f,signature), collapse="#")

.guessPackageName <- function(env) {
    allObjects <- objects(env, all.names = TRUE)
    allObjects <- allObjects[is.na(match(allObjects, .functionsOverriden))]
    ## counts of packaages containing objects; objects not found don't count
    possible <- sort(table(unlist(lapply(allObjects, find))), decreasing = TRUE)
    message <- ""
    if(length(possible) == 0)
        stop("none of the objects in the source code could be found:  need to attach or specify the package")
    else if(length(possible) > 1L) {
        global <- match(".GlobalEnv", names(possible), 0)
        if(global > 0) {
            possible <- possible[-global] # even if it's the most common
        }
        if(length(possible) > 1L)
            warning(gettextf("objects found in multiple packages: using %s and ignoring %s",
                             sQuote(names(possible[[1L]])),
                             paste(sQuote(names(possible[-1L])),
                                   collapse = ", ")),
                    domain = NA)
    }
    sub("package:","", names(possible[1L])) # the package name, or .GlobalEnv
}

## extract the new definitions from the source file
evalSource <- function(source, package = "", lock = TRUE, cache = FALSE) {
    if(!nzchar(package))
        envp <- .GlobalEnv # will look for the package after evaluating source
    else {
        pstring <- paste("package",package, sep=":")
        packageIsVisible <- pstring %in% search()
        if(packageIsVisible) {
            envp <- as.environment(pstring)
            envns <- tryCatch(asNamespace(package), error = function(cond) NULL)
        }
        else {
            envp <- tryCatch(asNamespace(package), error = function(cond) NULL)
            envns <- envp
        }
        if(is.null(envp))
            stop(gettextf("package %s is not attached and no namespace found for it",
                          sQuote(package)),
                 domain = NA)
    }
    env <- new("sourceEnvironment", new.env(parent = envp),
        packageName = package,
        sourceFile = (if(is.character(source)) source else ""))
    env$.packageName <- package # Fixme: should be done by an initialize method
    setCacheOnAssign(env, cache)
    if(is(source, "character"))
        for(text in source) sys.source(text, envir = env)
    else if(is(source, "connection")) sys.source(source, envir = env)
    else if(!is(source, "environment"))
        stop(gettextf("invalid 'source' argument: expected file names or a connection but got an object of class %s",
                      dQuote(class(source)[[1L]])),
             domain = NA)
    if(lock)
        lockEnvironment(env, bindings = TRUE) # no further changes allowed
    env
}

insertSource <- function(source, package = "",
                         functions = allPlainObjects(),
                         methods = (if(missing(functions)) allMethodTables() else NULL)
##                         ,classes = (if(missing(functions)) allClassDefs() else NULL)
                         , force = missing(functions) & missing(methods)
                     ){
    MPattern <- .TableMetaPattern()
    CPattern <- .ClassMetaPattern()
    allPlainObjects <- function()
        allObjects[!(grepl(MPattern, allObjects) | grepl(CPattern, allObjects) | ".cacheOnAssign" == allObjects)]
    allMethodTables <- function()
        allObjects[grepl(MPattern, allObjects)]
    allClassDefs <- function()
        allObjects[grepl(CPattern, allObjects)]
    differs <- function(f1, f2)
        !(identical(body(f1), body(f2)) && identical(args(f1), args(f2)))
    if(is.environment(source) && !nzchar(package)) {
        if(is(source, "sourceEnvironment"))
            package <- source@packageName
        else if(exists(".packageName", envir = source, inherits = FALSE))
            package <- get(".packageName", envir =source)
    }
    if(is(source, "environment"))
        env <- source
    else
        env <- evalSource(source, package, FALSE) # sourceEnvironment, unlocked
    envPackage <- getPackageName(env, FALSE)
    ## identify an environment and (if possible) namespace for the package
    envp <- parent.env(env)
    if(identical(envp, .GlobalEnv) || !nzchar(envPackage)) { # no package name in the eval, guess one
        if(!nzchar(package))
            package <- .guessPackageName(env) # use find() on objects in env
        if(identical(package, ".GlobalEnv"))
            envns <- NULL
        else {
            pname <- paste0("package:", package)
            envp <- tryCatch(as.environment(pname), error = function(cond)NULL)
            if(is.null(envp)) {
                envp <- tryCatch(as.environment(pname), error = function(cond)NULL)
                if(is.null(envp))
                    stop(gettextf(
                     "cannot find an environment corresponding to package name \'%s\"",
                     package), domain = NA)
            }
            envns <- tryCatch(asNamespace(package), error = function(cond)NULL)
        }
        if(nzchar(package))
            assign(".packageName", package, envir = env)
    }
    else {
        if(isNamespace(envp))
            envns <- envp
        else
            envns <- tryCatch(asNamespace(package), error = function(cond)NULL)
    }
    if(nzchar(envPackage) && envPackage != package)
        warning(gettextf("supplied package, %s, differs from package inferred from source, %s",
                         sQuote(package), sQuote(envPackage)),
                domain = NA)
    packageSlot(env) <- package
    ## at this point, envp is the target environment (package or other)
    ## and envns is the corresponding namespace if any, or NULL
    allObjects <- objects(envir = env, all.names = TRUE)
    ## Figure out what to trace.
    if(!missing(functions)) {
        notThere <- is.na(match(functions, allObjects))
        if(any(notThere)) {
            warning(gettextf("cannot insert these (not found in source): %s",
                    paste('"',functions[notThere],'"',
                          sep = "", collapse = ", ")),
                    domain = NA)
        }
    }
    .mnames <- allMethodTables()
    if(length(methods) > 0) {
        notThere <- sapply(methods,
         function(fname) (length(grep(fname, .mnames, fixed = TRUE)) == 0)
        )
        if(any(notThere)) {
            warning(gettextf("cannot insert methods for these functions (methods table not found in source): %s",
                    paste('"',methods[notThere],'"',
                          sep = "", collapse = ", ")),
                    domain = NA)
            methods <- methods[!notThere]
        }
        methodNames <- sapply(methods,
         function(fname) .mnames[[grep(fname, .mnames, fixed = TRUE)[[1]]]]
        )
    }
    else {
        methodNames <- .mnames
        methods <- sub(.TableMetaPrefix(), "", methodNames)
        methods <- sub(":.*","",methods)
    }
    ## if(!missing(classes)) {
    ##     .mnames <- allMethodNames()
    ##     notThere <- sapply(classes,
    ##      function(fname) length(grep(fname, .mnames, fixed = TRUE) == 0)
    ##     )
    ##     if(any(notThere)) {
    ##         warning(gettextf("Can't insert these classes (class definition not found in source): %s",
    ##                 paste('"',classes[notThere],'"',
    ##                       sep = "", collapse = ", ")),
    ##                 domain = NA)
    ##         classes <- classes[!notThere]
    ##     }
    ## }
    notTraceable <- newObjects <- objectsDone <- character()
    for(i in seq_along(functions)) {
        this <- functions[[i]]
        thisWhere <- NULL
        if(is.null(envns) ||
           exists(this, envir = envp, inherits = FALSE)) {
            envwhere <- envp
            thisWhere <- get(this, envir = envp)
        }
        else {
            envwhere <- envns
            if(is.environment(envns)  &&
               exists(this, envir = envns, inherits = FALSE))
                thisWhere <- get(this, envir = envns)
        }
        thisObj <- get(this, envir = env)
        if(is.function(thisObj) && is.function(thisWhere)
           && differs(thisObj, thisWhere)) {
            suppressMessages(
               .TraceWithMethods(this, where = envwhere, edit = env))
            objectsDone <- c(objectsDone, this)
        }
        else if(force)
            assign(this, thisObj, envir = envwhere)
        else if(!is.function(thisObj))
            notTraceable <- c(notTraceable, this)
        else if(is.null(thisWhere))
            newObjects <- c(newObjects, this)
    }
    if(length(notTraceable) > 0)
        message(gettextf("Non-function objects are not currently inserted (not traceable): %s",
                         paste(notTraceable, collapse = ", ")), domain = NA)
    if(length(newObjects) > 0)
        message(gettextf("New functions are not currently inserted (not untraceable): %s",
                         paste(newObjects, collapse = ", ")), domain = NA)
    if(length(objectsDone) > 0)
        message(gettextf("Modified functions inserted through trace(): %s",
                         paste(objectsDone, collapse = ", ")), domain = NA)
    for(i in seq_along(methods)) {
        .copyMethods(methods[[i]], methodNames[[i]], env, envp)
    }
    ## for(class in classes) {
    ##     .copyClass(class, env, envwhere)
    ## }
    ## return the environment, after cleaning up the dummy functions and
    ## adding a time stamp, if the source was parssed on this call
    if(!is.environment(source)) {
        lockEnvironment(env, bindings = TRUE) # no further changes allowed
        invisible(env)
    }
    else
        invisible(source)
}

.copyMethods <- function(f, tableName, env, envwhere) {
    differs <- function(o1, o2)
        !(is.function(o2) && # o2 can be NULL
          identical(body(o2), body(o2)) && identical(args(o1), args(o2)))
    table <- get(tableName, envir=env)
    fdef <- getGeneric(f, where = envwhere)
    if(!is(fdef, "genericFunction")) {
        message(gettextf("%s() is not a generic function in the target environment -- methods will not be inserted",
                         f), domain = NA)
        return(NULL)
    }
    curTable <- getMethodsForDispatch(fdef)
    allObjects <- objects(table, all.names = TRUE)
    methodsInserted <- character()
    if(length(allObjects) > 0) {
        for(this in allObjects) {
            def <- get(this, envir = table)
            curdef <- (if(exists(this, envir = curTable, inherits = FALSE))
                get(this, envir = curTable)
                       else NULL)
            if(differs(def, curdef)) {
                suppressMessages(
                   .TraceWithMethods(f, signature = this, where = envwhere,
                              edit = env))
                methodsInserted <- c(methodsInserted, this)
            }
        }
        if(length(methodsInserted) > 0)
            message(gettextf("Methods inserted for function %s(): %s",
                  f, paste(methodsInserted, collapse =", ")),
                  domain = NA)
    }
}

.copyClass <- function(class, env, envwhere) {
    message("Pretend we inserted class ", class, domain = NA)
}

.findNewDefForTrace <- function(what, signature, env, package) {
    if(is.null(signature)) {
        if(exists(what, envir = env, inherits = FALSE))
            newObject <- get(what, envir = env)
        else
            stop(gettextf("no definition for object %s found in tracing environment",
                          sQuote(what), source),
                 domain = NA)
    }
    else {
        ## we don't know the package for the generic (which may not
        ## be active), so we search for the string w/o package
        table <- .TableMetaName(what, "")
        allObjects <- objects(env, all.names = TRUE)
        i <- grep(table, allObjects, fixed = TRUE)
        if(length(i) == 1)
            table <- get(allObjects[[i]], envir = env)
        else if(length(i) >1) {
            table <- allObjects[[i[[1]]]]
            warning(gettextf("multiple generics match pattern, using table %s", table)
                , domain = NA)
            table <- get(table, envir = env)
        }
        else
            stop(gettextf("does not seem to be a method table for generic %s in tracing environment",
                          sQuote(what)),
                 domain = NA)
        if(exists(signature, envir = table, inherits = FALSE))
          newObject <- get(signature, envir = table)
        else
          stop(gettextf("no method in methods table for %s for signature %s",
                        sQuote(what),
                        sQuote(signature)),
               domain = NA)
    }
    newObject
}
#  File src/library/methods/R/zzz.R
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

## utils::globalVariables("...onLoad")

## Initial version of .onLoad
...onLoad  <-
  ## Initialize the methods package.
  function(libname, pkgname)
{
    where <- environment(sys.function())  # the namespace
    initMethodDispatch(where)
    ## temporary empty reference to the package's own namespace
    assign(".methodsNamespace", new.env(), envir = where)
    .Call(C_R_set_method_dispatch, TRUE)
    cat("initializing class and method definitions ...")
    ## set up default prototype (uses .Call so has be at load time)
    assign(".defaultPrototype", .Call(C_Rf_allocS4Object), envir = where)
    assign(".SealedClasses", character(), envir = where)
    .InitClassDefinition(where)
    assign("possibleExtends", .possibleExtends, envir = where)
    .InitBasicClasses(where)
    .initClassSupport(where)
    .InitMethodsListClass(where)
    .setCoerceGeneric(where)
    ## now install the non-dummy versions of some functions
    assign(".classEnv", ..classEnv, envir = where)
    assign("makeGeneric", .makeGeneric, envir = where)
    assign("newClassRepresentation", .newClassRepresentation, envir = where)
    assign("classGeneratorFunction", .classGeneratorFunction, envir = where)
    assign(".mergeClassDefSlots", ..mergeClassDefSlots, envir = where)
    assign(".addToMetaTable", ..addToMetaTable, envir = where)
    assign(".extendsForS3", ..extendsForS3, envir = where)
    .makeBasicFuns(where)
    rm(.makeGeneric, .newClassRepresentation, .possibleExtends,
       ..mergeClassDefSlots, .classGeneratorFunction, ..classEnv,
       ..addToMetaTable, ..extendsForS3, envir = where)
    .InitMethodDefinitions(where)
    .InitShowMethods(where)
    assign(".isPrototype", ..isPrototype, envir = where)
    .InitClassUnion(where)
    .InitS3Classes(where)
    .InitSpecialTypesAndClasses(where)
    .InitTraceFunctions(where)
    .InitRefClasses(where)
    ## now seal the classes defined in the package
    for(cl in get(".SealedClasses", where))
        sealClass(cl, where)
    assign("isSealedMethod", .isSealedMethod, envir = where)
    assign(".requirePackage", ..requirePackage, envir = where)
    ## initialize implicit generics for base package
    ## Note that this is done before making a non-vacuous implicitGeneric()
    ## so that non-default signatures are allowed in setGeneric()
    .initImplicitGenerics(where)
    assign("implicitGeneric", .implicitGeneric, envir = where)
    cacheMetaData(where, TRUE, searchWhere = .GlobalEnv, FALSE)
    assign(".checkRequiredGenerics", ..checkRequiredGenerics, envir = where)
    assign(".methodPackageSlots", ..methodPackageSlots, envir = where)
    rm(..isPrototype, .isSealedMethod, ..requirePackage, .implicitGeneric,
       ..checkRequiredGenerics, ..methodPackageSlots, envir = where)
    ## unlock some bindings that must be modifiable
    unlockBinding(".BasicFunsList", where)
    assign(".saveImage", TRUE, envir = where)
    cat(" done\n")

    assign(".onLoad", ..onLoad, envir = where)
    rm(...onLoad, ..onLoad, envir = where)
    dbbase <- file.path(libname, pkgname, "R", pkgname)
    ns <- asNamespace(pkgname)
    vars <- ls(envir = ns, all.names = TRUE)
    ## we need to exclude the registration vars
    vars <- grep("^C_", vars, invert = TRUE, value = TRUE)
    tools:::makeLazyLoadDB(ns, dbbase, variables = vars)
}

## avoid warnings from static analysis code by extra call
.onLoad <- function(libname, pkgname) ...onLoad(libname, pkgname)

##  .onLoad for routine use, installed by ...onLoad
..onLoad <- function(libname, pkgname)
{
    where <- environment(sys.function())  # the namespace
    initMethodDispatch(where)
    .Call(C_R_set_method_dispatch, TRUE)
    assign(".methodsNamespace", where, where)
    ## assign to baseenv also, signalling methods loaded
    assign(".methodsNamespace", where, baseenv())
    if(Sys.getenv("R_S4_BIND") == "active")
        methods:::bind_activation(TRUE)
}

.onUnload <- function(libpath)
{
    message("unloading 'methods' package ...") # see when this is called
    .isMethodsDispatchOn(FALSE)
    methods:::bind_activation(FALSE)
    library.dynam.unload("methods", libpath)
}


.onAttach <- function(libname, pkgname)
{
    env <- environment(sys.function())
    ## unlock some bindings that must be modifiable
    unlockBinding(".BasicFunsList", env)
    if(methods:::.hasS4MetaData(.GlobalEnv)) {
        result <- try(cacheMetaData(.GlobalEnv, TRUE))
        ## still attach  methods package if global env has bad objets
        if(is(result, "try-error"))
          warning("apparently bad method or class metadata in saved environment;\n",
                  "move the file or remove the class/method")
    }
}

.onDetach <- function(libpath) methods:::.onUnload(libpath)

## redefining it here, invalidates the one above:
## Why don't we unload "methods" on detach() ?
.onDetach <- function(libpath) .isMethodsDispatchOn(FALSE)

## used for .methodsIsLoaded
.saveImage <- FALSE

## want ASCII quotes, not fancy nor translated ones
.dQ <- function (x) paste0('"', x, '"')
