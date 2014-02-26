#  File src/library/grid/R/components.R
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

grid.collection <- function(..., gp=gpar(), draw=TRUE, vp=NULL) {
  .Deprecated("gTree")
  gc <- gTree(children=gList(...), gp=gp, vp=vp, cl="collection")
  if (draw)
    grid.draw(gc)
  gc
}

######################################
# AXES
######################################

# Axes are extended from the "gTree" class
# This means that the standard (e.g., draw.details)
# methods for gTrees will apply

# The children of an axis are fixed to be:

# NOTE that the `at' parameter is numeric (i.e., NOT a unit) for
# grid.xaxis and grid.yaxis.  These functions assume a unit for the `at'
# values rather than letting the user specify a unit.

validDetails.axis <- function(x) {
  if (!is.null(x$at)) {
    x$at <- as.numeric(x$at)
    if (length(x$at) < 1 ||
        !is.finite(x$at))
      stop("invalid 'at' location in 'axis'")
  }
  if (!is.logical(x$label)) {
    # labels specified
    # Can only spec labels if at is not NULL
    if (is.null(x$at))
      stop("invalid to specify axis labels when 'at' is NULL")
    # Must be either language object or string
    x$label <- as.graphicsAnnot(x$label)
    # Must be same number of labels as "at" locations
    if (length(x$label) != length(x$at))
      stop("'labels' and 'at' locations must have same length")
  }
  x$main <- as.logical(x$main)
  x
}

makeContent.xaxis <- function(x) {
    # If x$at is NULL, then we must calculate the
    # tick marks on-the-fly
    if (is.null(x$at)) {
        x$at <- grid.pretty(current.viewport()$xscale)
        # Add the new output as children
        x <- addGrob(x, make.xaxis.major(x$at, x$main))
        x <- addGrob(x, make.xaxis.ticks(x$at, x$main))
        x <- updateXlabels(x)
        # Apply any edits relevant to children
        x <- applyEdits(x, x$edits)
    }
    x
}

# NOTE that this can't be for all axes because it needs to
# call make.XAXIS.ticks and make.XAXIS.labels
editDetails.xaxis <- function(x, specs) {
  slot.names <- names(specs)
  if ("at" %in% slot.names) {
    # NOTE that grid.edit has already set x$at to the new value
    # We might set at to NULL to get ticks recalculated at redraw
    if (is.null(x$at)) {
      x <- removeGrob(x, "major", warn=FALSE)
      x <- removeGrob(x, "ticks", warn=FALSE)
      x <- removeGrob(x, "labels", warn=FALSE)
    } else {
      x <- addGrob(x, make.xaxis.major(x$at, x$main))
      x <- addGrob(x, make.xaxis.ticks(x$at, x$main))
      x <- updateXlabels(x)
    }
  }
  if ("label" %in% slot.names) {
    if (!is.null(x$at))
      x <- updateXlabels(x)
  }
  if ("main" %in% slot.names)
    if (!is.null(x$at)) {
      x <- addGrob(x, make.xaxis.major(x$at, x$main))
      x <- addGrob(x, make.xaxis.ticks(x$at, x$main))
      x <- updateXlabels(x)
    }
  x
}

make.xaxis.major <- function(at, main) {
  if (main)
    y <- c(0, 0)
  else
    y <- c(1, 1)
  linesGrob(unit(c(min(at), max(at)), "native"),
            unit(y, "npc"), name="major")
}

make.xaxis.ticks <- function(at, main) {
  if (main) {
    tick.y0 <- unit(0, "npc")
    tick.y1 <- unit(-.5, "lines")
  }
  else {
    tick.y0 <- unit(1, "npc")
    tick.y1 <- unit(1, "npc") + unit(.5, "lines")
  }
  segmentsGrob(unit(at, "native"), tick.y0,
               unit(at, "native"), tick.y1,
               name="ticks")
}

make.xaxis.labels <- function(at, label, main) {
  # FIXME:  labels only character versions of "at"
  if (main)
    label.y <- unit(-1.5, "lines")
  else
    label.y <- unit(1, "npc") + unit(1.5, "lines")
  if (is.logical(label))
    labels <- as.character(at)
  else
    labels <- label
  textGrob(labels, unit(at, "native"), label.y,
           just="centre", rot=0,
           check.overlap=TRUE, name="labels")
}

updateXlabels <- function(x) {
  if (is.logical(x$label) && !x$label)
    removeGrob(x, "labels", warn=FALSE)
  else
    addGrob(x, make.xaxis.labels(x$at, x$label, x$main))
}

xaxisGrob <- function(at=NULL, label=TRUE, main=TRUE,
                      edits=NULL,
                      name=NULL, gp=gpar(), vp=NULL) {
  grid.xaxis(at=at, label=label, main=main,
             edits=edits,
             name=name, gp=gp, draw=FALSE, vp=vp)
}

# The "main" x-axis is on the bottom when vp$origin is "bottom.*"
# and on the top when vp$origin is "top.*"
grid.xaxis <- function(at=NULL, label=TRUE, main=TRUE,
                       edits=NULL, name=NULL, gp=gpar(),
                       draw=TRUE, vp=NULL) {
  if (is.null(at)) {
    # We do not have enough information to make the ticks and labels
    major <- NULL
    ticks <- NULL
    labels <- NULL
  } else {
    major <- make.xaxis.major(at, main)
    ticks <- make.xaxis.ticks(at, main)
    if (is.logical(label) && length(label) == 0)
	stop("logical 'label' supplied of length 0")
    if (is.logical(label) && !label)
      labels <- NULL
    else
      labels <- make.xaxis.labels(at, label, main)
  }
  xg <- applyEdits(gTree(at=at, label=label, main=main,
                         children=gList(major, ticks, labels),
                         edits=edits,
                         name=name, gp=gp, vp=vp,
                         cl=c("xaxis", "axis")),
                   edits)
  if (draw)
    grid.draw(xg)
  invisible(xg)
}

makeContent.yaxis <- function(x) {
    # If x$at is NULL, then we must calculate the
    # tick marks on-the-fly
    if (is.null(x$at)) {
        x$at <- grid.pretty(current.viewport()$yscale)
        # Add the new output as children
        x <- addGrob(x, make.yaxis.major(x$at, x$main))
        x <- addGrob(x, make.yaxis.ticks(x$at, x$main))
        x <- updateYlabels(x)
        # Apply any edits relevant to children
        x <- applyEdits(x, x$edits)
    }
    x
}

editDetails.yaxis <- function(x, specs) {
  slot.names <- names(specs)
  if ("at" %in% slot.names) {
    if (is.null(x$at)) {
      x <- removeGrob(x, "major", warn=FALSE)
      x <- removeGrob(x, "ticks", warn=FALSE)
      x <- removeGrob(x, "labels", warn=FALSE)
    } else {
      x <- addGrob(x, make.yaxis.major(x$at, x$main))
      x <- addGrob(x, make.yaxis.ticks(x$at, x$main))
      x <- updateYlabels(x)
    }
  }
  if ("label" %in% slot.names) {
    if (!is.null(x$at))
      x <- updateYlabels(x)
  }
  if ("main" %in% slot.names)
    if (!is.null(x$at)) {
      x <- addGrob(x, make.yaxis.major(x$at, x$main))
      x <- addGrob(x, make.yaxis.ticks(x$at, x$main))
      x <- updateYlabels(x)
    }
  x
}

make.yaxis.major <- function(at, main) {
  if (main)
    x <- c(0, 0)
  else
    x <- c(1, 1)
  linesGrob(unit(x, "npc"), unit(c(min(at), max(at)), "native"),
            name="major")
}

make.yaxis.ticks <- function(at, main) {
  if (main) {
    tick.x0 <- unit(0, "npc")
    tick.x1 <- unit(-.5, "lines")
  }
  else {
    tick.x0 <- unit(1, "npc")
    tick.x1 <- unit(1, "npc") + unit(.5, "lines")
  }
  segmentsGrob(tick.x0, unit(at, "native"),
               tick.x1, unit(at, "native"),
               name="ticks")
}

make.yaxis.labels <- function(at, label, main) {
  if (main) {
    hjust <- "right"
    label.x <- unit(-1, "lines")
  }
  else {
    hjust <- "left"
    label.x <- unit(1, "npc") + unit(1, "lines")
  }
  just <- c(hjust, "centre")
  if (is.logical(label))
    labels <- as.character(at)
  else
    labels <- label
  textGrob(labels, label.x, unit(at, "native"),
           just=just, rot=0, check.overlap=TRUE, name="labels")
}

updateYlabels <- function(x) {
  if (is.logical(x$label) && !x$label)
    removeGrob(x, "labels", warn=FALSE)
  else
    addGrob(x, make.yaxis.labels(x$at, x$label, x$main))
}

yaxisGrob <- function(at=NULL, label=TRUE, main=TRUE,
                      edits=NULL,
                      name=NULL, gp=gpar(), vp=NULL) {
  grid.yaxis(at=at, label=label, main=main, edits=edits,
             name=name, gp=gp, draw=FALSE, vp=vp)
}

# The "main" y-axis is on the left when vp$origin is "*.left"
# and on the right when vp$origin is "*.right"
grid.yaxis <- function(at=NULL, label=TRUE, main=TRUE,
                       edits=NULL,
                       name=NULL, gp=gpar(),
                       draw=TRUE, vp=NULL) {
  if (is.null(at)) {
    # We do not have enough information to make the ticks and labels
    major <- NULL
    ticks <- NULL
    labels <- NULL
  } else {
    major <- make.yaxis.major(at, main)
    ticks <- make.yaxis.ticks(at, main)
    if (is.logical(label) && length(label) == 0)
	stop("logical 'label' supplied of length 0")
    if (is.logical(label) && !label)
      labels <- NULL
    else
      labels <- make.yaxis.labels(at, label, main)
  }
  yg <- applyEdits(gTree(at=at, label=label, main=main,
                         children=gList(major, ticks, labels),
                         edits=edits,
                         name=name, gp=gp, vp=vp,
                         cl=c("yaxis", "axis")),
                   edits)
  if (draw)
    grid.draw(yg)
  invisible(yg)
}

######################################
# Simple "side-effect" plotting functions
######################################

grid.grill <- function(h=unit(seq(0.25, 0.75, 0.25), "npc"),
                       v=unit(seq(0.25, 0.75, 0.25), "npc"),
                       default.units="npc",
                       gp=gpar(col="grey"), vp=NULL) {
  if (!is.unit(h))
    h <- unit(h, default.units)
  if (!is.unit(v))
    v <- unit(v, default.units)
  # FIXME:  Should replace for loop and call to grid.lines with call to grid.segments
  # once the latter exists
  if (!is.null(vp))
    pushViewport(vp)
  grid.segments(v, unit(0, "npc"), v, unit(1, "npc"), gp=gp)
  grid.segments(unit(0, "npc"), h, unit(1, "npc"), h, gp=gp)
  if (!is.null(vp))
    popViewport()
}

#  File src/library/grid/R/curve.R
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


###############################
# CURVE primitive
###############################

calcOrigin <- function(x1, y1, x2, y2, origin, hand) {
    # Positive origin means origin to the "right"
    # Negative origin means origin to the "left"
    xm <- (x1 + x2)/2
    ym <- (y1 + y2)/2
    dx <- x2 - x1
    dy <- y2 - y1
    slope <- dy/dx
    oslope <- -1/slope
    # The origin is a point somewhere along the line between
    # the end points, rotated by 90 (or -90) degrees
    # Two special cases:
    # If slope is non-finite then the end points lie on a vertical line, so
    # the origin lies along a horizontal line (oslope = 0)
    # If oslope is non-finite then the end points lie on a horizontal line,
    # so the origin lies along a vertical line (oslope = Inf)
    tmpox <- ifelse(!is.finite(slope),
                    xm,
                    ifelse(!is.finite(oslope),
                           xm + origin*(x2 - x1)/2,
                           xm + origin*(x2 - x1)/2))
    tmpoy <- ifelse(!is.finite(slope),
                    ym + origin*(y2 - y1)/2,
                    ifelse(!is.finite(oslope),
                           ym,
                           ym + origin*(y2 - y1)/2))
    # ALWAYS rotate by -90 about midpoint between end points
    # Actually no need for "hand" because "origin" also
    # encodes direction
    # sintheta <- switch(hand, left=-1, right=1)
    sintheta <- -1
    ox <- xm - (tmpoy - ym)*sintheta
    oy <- ym + (tmpox - xm)*sintheta

    list(x=ox, y=oy)
}

# Given ncp*ncurve vector of values, ncurve vector of start values,
# ncurve vector of end values, ncurve vector of end logicals,
# combine start or end values with original values based on logicals
interleave <- function(ncp, ncurve, val, sval, eval, e) {
    sval <- rep(sval, length.out=ncurve)
    eval <- rep(eval, length.out=ncurve)
    result <- matrix(NA, ncol=ncurve, nrow=ncp+1)
    m <- matrix(val, ncol=ncurve)
    for (i in 1L:ncurve) {
        if (e[i])
            result[,i] <- c(m[,i], eval[i])
        else
            result[,i] <- c(sval[i], m[,i])
    }
    as.numeric(result)
}

# Calculate a "square" set of end points to calculate control points from
# NOTE: end points may be vector
calcSquareControlPoints <- function(x1, y1, x2, y2,
                                    curvature, angle, ncp,
                                    debug=FALSE) {
    xm <- (x1 + x2)/2
    ym <- (y1 + y2)/2
    dx <- x2 - x1
    dy <- y2 - y1
    slope <- dy/dx
    oslope <- -1/slope

    # FIXME:  There MUST be a more compact way of calculating the
    # new end point!
    end <- (slope > 1 |
            (slope < 0 & slope > -1))
    if (curvature < 0)
        end <- !end
    startx <- ifelse(end,
                     x1,
                     ifelse(abs(slope) > 1,
                            newx <- x2 - dx,
                            newx <- x2 - sign(slope)*dy))
    starty <- ifelse(end,
                     y1,
                     ifelse(abs(slope) > 1,
                            newy <- y2 - sign(slope)*dx,
                            newy <- y2 - dy))
    endx <- ifelse(end,
                   ifelse(abs(slope) > 1,
                          newx <- x1 + dx,
                          newx <- x1 + sign(slope)*dy),
                   x2)
    endy <- ifelse(end,
                   ifelse(abs(slope) > 1,
                          newy <- y1 + sign(slope)*dx,
                          newy <- y1 + dy),
                   y2)

    cps <- calcControlPoints(startx, starty, endx, endy,
                             curvature, angle, ncp,
                             debug)

    # Intereave control points and extra "square" control points
    ncurve <- length(x1)
    cps$x <- interleave(ncp, ncurve, cps$x, startx, endx, end)
    cps$y <- interleave(ncp, ncurve, cps$y, starty, endy, end)

    list(x=cps$x, y=cps$y, end=end)
}

# Find origin of rotation
# Rotate around that origin
calcControlPoints <- function(x1, y1, x2, y2, curvature, angle, ncp,
                              debug=FALSE) {
    # Negative curvature means curve to the left
    # Positive curvature means curve to the right
    # Special case curvature = 0 (straight line) has been handled
    xm <- (x1 + x2)/2
    ym <- (y1 + y2)/2
    dx <- x2 - x1
    dy <- y2 - y1
    slope <- dy/dx
    oslope <- -1/slope

    # Calculate "corner" of region to produce control points in
    # (depends on 'angle', which MUST lie between 0 and 180)
    # Find by rotating start point by angle around mid point
    if (is.null(angle)) {
        # Calculate angle automatically
        angle <- ifelse(slope < 0,
                        2*atan(abs(slope)),
                        2*atan(1/slope))
    } else {
        angle <- angle/180*pi
    }
    sina <- sin(angle)
    cosa <- cos(angle)
    # FIXME:  special case of vertical or horizontal line ?
    cornerx <- xm + (x1 - xm)*cosa - (y1 - ym)*sina
    cornery <- ym + (y1 - ym)*cosa + (x1 - xm)*sina

    # Debugging
    if (debug) {
        grid.points(cornerx, cornery, default.units="inches",
                    pch=16, size=unit(3, "mm"),
                    gp=gpar(col="grey"))
    }

    # Calculate angle to rotate region by to align it with x/y axes
    beta <- -atan((cornery - y1)/(cornerx - x1))
    sinb <- sin(beta)
    cosb <- cos(beta)
    # Rotate end point about start point to align region with x/y axes
    newx2 <- x1 + dx*cosb - dy*sinb
    newy2 <- y1 + dy*cosb + dx*sinb

    # Calculate x-scale factor to make region "square"
    # FIXME:  special case of vertical or horizontal line ?
    scalex <- (newy2 - y1)/(newx2 - x1)
    # Scale end points to make region "square"
    newx1 <- x1*scalex
    newx2 <- newx2*scalex

    # Calculate the origin in the "square" region
    # (for rotating start point to produce control points)
    # (depends on 'curvature')
    # 'origin' calculated from 'curvature'
    ratio <- 2*(sin(atan(curvature))^2)
    origin <- curvature - curvature/ratio
    # 'hand' also calculated from 'curvature'
    if (curvature > 0)
        hand <- "right"
    else
        hand <- "left"
    oxy <- calcOrigin(newx1, y1, newx2, newy2, origin, hand)
    ox <- oxy$x
    oy <- oxy$y

    # Calculate control points
    # Direction of rotation depends on 'hand'
    dir <- switch(hand,
                  left=-1,
                  right=1)
    # Angle of rotation depends on location of origin
    maxtheta <- pi + sign(origin*dir)*2*atan(abs(origin))
    theta <- seq(0, dir*maxtheta,
                 dir*maxtheta/(ncp + 1))[c(-1, -(ncp + 2))]
    costheta <- cos(theta)
    sintheta <- sin(theta)
    # May have BOTH multiple end points AND multiple
    # control points to generate (per set of end points)
    # Generate consecutive sets of control points by performing
    # matrix multiplication
    cpx <- ox + ((newx1 - ox) %*% t(costheta)) -
        ((y1 - oy) %*% t(sintheta))
    cpy <- oy + ((y1 - oy) %*% t(costheta)) +
        ((newx1 - ox) %*% t(sintheta))

    # Reverse transformations (scaling and rotation) to
    # produce control points in the original space
    cpx <- cpx/scalex
    sinnb <- sin(-beta)
    cosnb <- cos(-beta)
    finalcpx <- x1 + (cpx - x1)*cosnb - (cpy - y1)*sinnb
    finalcpy <- y1 + (cpy - y1)*cosnb + (cpx - x1)*sinnb

    # Debugging
    if (debug) {
        ox <- ox/scalex
        fox <- x1 + (ox - x1)*cosnb - (oy - y1)*sinnb
        foy <- y1 + (oy - y1)*cosnb + (ox - x1)*sinnb
        grid.points(fox, foy, default.units="inches",
                    pch=16, size=unit(1, "mm"),
                    gp=gpar(col="grey"))
        grid.circle(fox, foy, sqrt((ox - x1)^2 + (oy - y1)^2),
                    default.units="inches",
                    gp=gpar(col="grey"))
    }

    list(x=as.numeric(t(finalcpx)), y=as.numeric(t(finalcpy)))
}

# Debugging
cbDiagram <- function(x1, y1, x2, y2, cps) {
    grid.segments(x1, y1, x2, y2,
                gp=gpar(col="grey"),
                default.units="inches")
    grid.points(x1, y1, pch=16, size=unit(1, "mm"),
                gp=gpar(col="green"),
                default.units="inches")
    grid.points(x2, y2, pch=16, size=unit(1, "mm"),
                gp=gpar(col="red"),
                default.units="inches")
    grid.points(cps$x, cps$y, pch=16, size=unit(1, "mm"),
                default.units="inches",
                gp=gpar(col="blue"))
}

straightCurve <- function(x1, y1, x2, y2, arrow, debug) {
    if (debug) {
        xm <- (x1 + x2)/2
        ym <- (y1 + y2)/2
        cbDiagram(x1, y1, x2, y2, list(x=xm, y=ym))
    }

    segmentsGrob(x1, y1, x2, y2,
                 default.units="inches",
                 arrow=arrow, name="segment")
}

# Return a gTree (even if it only has one grob as a child)
# because that is the only way to get more than one child
# to draw
calcCurveGrob <- function(x, debug) {
    x1 <- x$x1
    x2 <- x$x2
    y1 <- x$y1
    y2 <- x$y2
    curvature <- x$curvature
    angle <- x$angle
    ncp <- x$ncp
    shape <- x$shape
    square <- x$square
    squareShape <- x$squareShape
    inflect <- x$inflect
    arrow <- x$arrow
    open <- x$open

    # Calculate a set of control points based on:
    # 'curvature', ' angle', and 'ncp',
    # and the start and end point locations.

    # The origin is a point along the perpendicular bisector
    # of the line between the end points.

    # The control points are found by rotating the end points
    # about the origin.

    # Do everything in inches to make things easier.
    # Because this is within a makeContent() method,
    # the conversions will not be an
    # issue (in terms of device resizes).
    x1 <- convertX(x1, "inches", valueOnly=TRUE)
    y1 <- convertY(y1, "inches", valueOnly=TRUE)
    x2 <- convertX(x2, "inches", valueOnly=TRUE)
    y2 <- convertY(y2, "inches", valueOnly=TRUE)

    # Outlaw identical end points
    if (any(x1 == x2 & y1 == y2))
        stop("end points must not be identical")

    # Rep locations to allow multiple curves from single call
    maxn <- max(length(x1),
                length(y1),
                length(x2),
                length(y2))
    x1 <- rep(x1, length.out=maxn)
    y1 <- rep(y1, length.out=maxn)
    x2 <- rep(x2, length.out=maxn)
    y2 <- rep(y2, length.out=maxn)
    if (!is.null(arrow))
        arrow <- rep(arrow, length.out=maxn)

    if (curvature == 0) {
        children <- gList(straightCurve(x1, y1, x2, y2, arrow, debug))
    } else {
        # Treat any angle less than 1 or greater than 179 degrees
        # as a straight line
        # Takes care of some nasty limit effects as well as simplifying
        # things
        if (angle < 1 || angle > 179) {
            children <- gList(straightCurve(x1, y1, x2, y2, arrow, debug))
        } else {
            # Handle 'square' vertical and horizontal lines
            # separately
            if (square && any(x1 == x2 | y1 == y2)) {
                subset <- x1 == x2 | y1 == y2
                straightGrob <- straightCurve(x1[subset], y1[subset],
                                               x2[subset], y2[subset],
                                               arrow, debug)
                # Remove these from the curves to draw
                x1 <- x1[!subset]
                x2 <- x2[!subset]
                y1 <- y1[!subset]
                y2 <- y2[!subset]
                if (!is.null(arrow))
                    arrow <- arrow[!subset]
            } else {
                straightGrob <- NULL
            }
            ncurve <- length(x1)
            # If nothing to draw, we're done
            if (ncurve == 0) {
                children <- gList(straightGrob)
            } else {
                if (inflect) {
                    xm <- (x1 + x2)/2
                    ym <- (y1 + y2)/2
                    shape1 <- rep(rep(shape, length.out=ncp), ncurve)
                    shape2 <- rev(shape1)
                    if (square) {
                      # If 'square' then add an extra control point
                        cps1 <- calcSquareControlPoints(x1, y1, xm, ym,
                                                        curvature, angle,
                                                        ncp,
                                                        debug=debug)
                        cps2 <- calcSquareControlPoints(xm, ym, x2, y2,
                                                        -curvature, angle,
                                                        ncp,
                                                        debug=debug)
                        shape1 <- interleave(ncp, ncurve, shape1,
                                             squareShape, squareShape,
                                             cps1$end)
                        shape2 <- interleave(ncp, ncurve, shape2,
                                             squareShape, squareShape,
                                             cps2$end)
                        ncp <- ncp + 1
                    } else {
                        cps1 <- calcControlPoints(x1, y1, xm, ym,
                                                  curvature, angle, ncp,
                                                  debug=debug)
                        cps2 <- calcControlPoints(xm, ym, x2, y2,
                                                  -curvature, angle, ncp,
                                                  debug=debug)
                    }

                    if (debug) {
                        cbDiagram(x1, y1, xm, ym, cps1)
                        cbDiagram(xm, ym, x2, y2, cps2)
                    }

                    idset <- 1L:ncurve
                    splineGrob <-
                        xsplineGrob(c(x1, cps1$x, xm, cps2$x, x2),
                                    c(y1, cps1$y, ym, cps2$y, y2),
                                    id=c(idset, rep(idset, each=ncp),
                                      idset, rep(idset, each=ncp),
                                      idset),
                                    default.units="inches",
                                    shape=c(rep(0, ncurve), shape1,
                                      rep(0, ncurve), shape2,
                                      rep(0, ncurve)),
                                    arrow=arrow, open=open,
                                    name="xspline")
                    if (is.null(straightGrob)) {
                        children <- gList(splineGrob)
                    } else {
                        children <- gList(straightGrob, splineGrob)
                    }
                } else {
                    shape <- rep(rep(shape, length.out=ncp), ncurve)
                    if (square) {
                      # If 'square' then add an extra control point
                        cps <- calcSquareControlPoints(x1, y1, x2, y2,
                                                       curvature, angle,
                                                       ncp,
                                                       debug=debug)
                        shape <- interleave(ncp, ncurve, shape,
                                            squareShape, squareShape,
                                            cps$end)
                        ncp <- ncp + 1
                    } else {
                        cps <- calcControlPoints(x1, y1, x2, y2,
                                                 curvature, angle, ncp,
                                                 debug=debug)
                    }
                    if (debug) {
                        cbDiagram(x1, y1, x2, y2, cps)
                    }

                    idset <- 1L:ncurve
                    splineGrob <- xsplineGrob(c(x1, cps$x, x2),
                                              c(y1, cps$y, y2),
                                              id=c(idset,
                                                rep(idset, each=ncp), idset),
                                              default.units="inches",
                                              shape=c(rep(0, ncurve), shape,
                                                rep(0, ncurve)),
                                              arrow=arrow, open=open,
                                              name="xspline")
                    if (is.null(straightGrob)) {
                        children <- gList(splineGrob)
                    } else {
                        children <- gList(straightGrob, splineGrob)
                    }
                }
            }
        }
    }
    gTree(children=children,
          name=x$name, gp=x$gp, vp=x$vp)
}

validDetails.curve <- function(x) {
    if ((!is.unit(x$x1) || !is.unit(x$y1)) ||
        (!is.unit(x$x2) || !is.unit(x$y2)))
        stop("'x1', 'y1', 'x2', and 'y2' must be units")
    x$curvature <- as.numeric(x$curvature)
    x$angle <- x$angle %% 180
    x$ncp <- as.integer(x$ncp)
    if (x$shape < -1 || x$shape > 1)
        stop("'shape' must be between -1 and 1")
    x$square <- as.logical(x$square)
    if (x$squareShape < -1 || x$squareShape > 1)
        stop("'squareShape' must be between -1 and 1")
    x$inflect <- as.logical(x$inflect)
    if (!is.null(x$arrow) && !inherits(x$arrow, "arrow"))
        stop("'arrow' must be an arrow object or NULL")
    x$open <- as.logical(x$open)
    x
}

makeContent.curve <- function(x) {
    calcCurveGrob(x, x$debug)
}

xDetails.curve <- function(x, theta) {
    cg <- calcCurveGrob(x, FALSE)
    # Could do better here
    # (result for more than 1 child is basically to give up)
    if (length(cg$children) == 1)
        xDetails(cg$children[[1]], theta)
    else
        xDetails(cg, theta)
}

yDetails.curve <- function(x, theta) {
    cg <- calcCurveGrob(x, FALSE)
    if (length(cg$children) == 1)
        yDetails(cg$children[[1]], theta)
    else
        yDetails(cg, theta)
}

widthDetails.curve <- function(x) {
    cg <- calcCurveGrob(x, FALSE)
    if (length(cg$children) == 1)
        widthDetails(cg$children[[1]])
    else
        widthDetails(cg)
}

heightDetails.curve <- function(x) {
    cg <- calcCurveGrob(x, FALSE)
    if (length(cg$children) == 1)
        heightDetails(cg$children[[1]])
    else
        heightDetails(cg)
}

curveGrob <- function(x1, y1, x2, y2, default.units="npc",
                      curvature=1, angle=90, ncp=1,
                      shape=0.5, square=TRUE, squareShape=1,
                      inflect=FALSE, arrow=NULL, open=TRUE,
                      debug=FALSE,
                      name=NULL, gp=gpar(), vp=NULL) {
    # FIXME:  add arg checking
    # FIXME:  angle MUST be between 0 and 180
    if (!is.unit(x1))
        x1 <- unit(x1, default.units)
    if (!is.unit(y1))
        y1 <- unit(y1, default.units)
    if (!is.unit(x2))
        x2 <- unit(x2, default.units)
    if (!is.unit(y2))
        y2 <- unit(y2, default.units)
    gTree(x1=x1, y1=y1, x2=x2, y2=y2,
          curvature=curvature, angle=angle, ncp=ncp,
          shape=shape, square=square, squareShape=squareShape,
          inflect=inflect, arrow=arrow, open=open, debug=debug,
          name=name, gp=gp, vp=vp,
          cl="curve")
}

grid.curve <- function(...) {
    grid.draw(curveGrob(...))
}

# Calculate the curvature to use if you want to produce control
# points lying along the arc of a circle that spans theta degrees
# (Use ncp=8 and shape=-1 to actually produce such an arc)
arcCurvature <- function(theta) {
    # Avoid limiting cases (just draw a straight line)
    if (theta < 1 || theta > 359)
        return(0)
    angle <- 0.5*theta/180*pi
    1/sin(angle) - 1/tan(angle)
}

# Label grobs in a scene

#
#  Copyright (C) 1995-2012 The R Core Team
labelGrob <- function(grob, recurse, curdepth, depth, labelfun, ...) {
    UseMethod("labelGrob")
}

# The default grob label needs to do some calculations
# on sizes so need a drawDetails method to get the
# calculations right
drawDetails.groblabel <- function(x, ...) {
    gw <- convertWidth(grobWidth(x$grob), "inches", valueOnly=TRUE)
    gh <- convertHeight(grobHeight(x$grob), "inches", valueOnly=TRUE)
    grid.rect(grobX(x$grob, "west"), grobY(x$grob, "south"),
              unit(gw, "inches"), unit(gh, "inches"),
              just=c("left", "bottom"), gp=x$gp)
    tw <- convertWidth(stringWidth(x$grob$name), "inches", valueOnly=TRUE)
    th <- convertHeight(stringHeight(x$grob$name), "inches", valueOnly=TRUE)
    eps <- .01
    # If grob is REALLY short, draw horiz at normal cex
    if (gh < eps) {
        rot <- 0
        cex <- 1
    # If grob is REALLY thin, draw vertical at normal cex
    } else if (gw < eps) {
        rot <- 90
        cex <- 1
    } else {
        gratio <- gh/gw
        if (gratio > 1 && tw > gw) {
            rot <- 90
            wratio <- th/gw
            hratio <- tw/gh
        } else {
            rot <- 0
            wratio <- tw/gw
            hratio <- th/gh
        }
        if (wratio > 1 || hratio > 1) {
            cex <- 1/max(wratio, hratio)
        } else {
            cex <- 1
        }
    }
    if (is.null(x$gp)) {
        x$gp <- gpar(cex=cex)
    } else {
        if (is.null(x$gp$cex))
            x$gp$cex <- cex
    }
    if (is.null(x$otherArgs$rot))
        x$otherArgs$rot <- rot
    do.call("grid.text", c(list(label=x$grob$name,
                                x=grobX(x$grob, "north"),
                                y=grobY(x$grob, "west"),
                                gp=x$gp),
                           x$otherArgs))
}

grobLabel <- function(grob,
                      gp=gpar(col=rgb(1, 0, 0, .5),
                        fill=rgb(1, 0, 0, .2)),
                      ...) {
    grob(grob=grob, gp=gp, otherArgs=list(...),
         cl="groblabel")
}

labelGrob.grob <- function(grob, recurse, curdepth, depth, labelfun, ...) {
    if (is.null(depth) || curdepth %in% depth) {
        gTree(children=gList(grob,
                labelfun(grob, ...)),
              # Name new gTree same as old grob so that
              # setGrob() approach works below
              # (when 'gPath' is specified)
              name=grob$name)
    } else {
        grob
    }
}

labelGrob.gTree <- function(grob, recurse, curdepth, depth, labelfun, ...) {
    if (recurse) {
        newChildren <- do.call("gList",
                               lapply(grob$children,
                                      labelGrob,
                                      recurse, curdepth + 1, depth,
                                      labelfun, ...))
        grob <- setChildren(grob, newChildren)
    }
    if (is.null(depth) || curdepth %in% depth) {
        gTree(children=gList(grob,
                labelfun(grob, ...)),
              name=grob$name)
    } else {
        grob
    }
}

showGrob <- function(x=NULL,
                     gPath=NULL, strict=FALSE, grep=FALSE,
                     recurse=TRUE, depth=NULL,
                     labelfun=grobLabel, ...) {
    if (is.null(x)) {
        # Label all or part of current scene
        # The grid display list is NOT affected
        # To remove labels use grid.redraw()
        if (is.null(gPath)) {
            # Show the current scene
            dl <- grid.Call(L_getDisplayList)[1L:grid:::grid.Call(L_getDLindex)]
            grid.newpage(recording=FALSE)
            # -1 because first element on DL is ROOT viewport
            lapply(dl[-1],
                   function(y) {
                       # Modify the grob to add a label
                       if (is.grob(y))
                           y <- labelGrob(y, recurse, 1, depth, labelfun, ...)
                       # Draw either the original object or the modified grob
                       grid.draw(y, recording=FALSE)
                   })
        } else {
            # Only label the bit of the current scene specified by gPath
            grobToLabel <- grid.get(gPath, strict=strict, grep=grep)
            # NOTE: have to 'wrap' because otherwise the grobs in the
            # captured scene have been altered
            scene <- grid.grab(wrap=TRUE)
            modScene <- setGrob(scene, gPath,
                                labelGrob(grobToLabel, recurse, 1, depth,
                                          labelfun, ...),
                                strict=strict, grep=grep)
            grid.newpage(recording=FALSE)
            grid.draw(modScene, recording=FALSE)
        }
    } else {
        # Assume grob is not current scene so start a new page
        grid.newpage()
        grid.draw(x)
        showGrob(NULL, gPath, strict, grep, recurse, depth, labelfun, ...)
    }
    invisible()
}

#############
# Labelling viewports in a scene
#############

# FIXME:  some of this code for vpLists and vpStacks and vpTrees
# assumes that the components of a vpList or vpStack or the
# vpTree parent can ONLY be a viewport (when in fact they can
# also be a vpList, vpStack, or vpTree!)

# Label a viewport
# Get physical aspect ratio of vp to determine whether to rotate
# Shrink text to fit in vp
# (Assumes that we are currently occupying 'vp'
#  so that conversions are correct)
labelVP <- function(vp, col) {
    vw <- convertWidth(unit(1, "npc"), "inches", valueOnly=TRUE)
    vh <- convertHeight(unit(1, "npc"), "inches", valueOnly=TRUE)
    tw <- convertWidth(stringWidth(vp$name), "inches", valueOnly=TRUE)
    th <- convertHeight(stringHeight(vp$name), "inches", valueOnly=TRUE)
    eps <- .01
    # If viewport is REALLY short, draw horiz at normal cex
    if (vh < eps) {
        rot <- 0
        cex <- 1
    # If viewport is REALLY thin, draw vertical at normal cex
    } else if (vw < eps) {
        rot <- 90
        cex <- 1
    } else {
        vratio <- vh/vw
        if (vratio > 1 && tw > vw) {
            rot <- 90
            wratio <- th/vw
            hratio <- tw/vh
        } else {
            rot <- 0
            wratio <- tw/vw
            hratio <- th/vh
        }
        if (wratio > 1 || hratio > 1) {
            cex <- 1/max(wratio, hratio)
        } else {
            cex <- 1
        }
    }
    # Violate any clipping that is in effect
    pushViewport(viewport(clip="off"))
    grid.text(vp$name, rot=rot, gp=gpar(col=col, cex=cex))
    upViewport()
}

# Draw a "viewport"
drawVP <- function(vp, curDepth, depth, col, fill, label) {
    UseMethod("drawVP")
}

drawVP.viewport <- function(vp, curDepth, depth, col, fill, label) {
    if (vp$name != "ROOT" &&
        (is.null(depth) || curDepth %in% depth)) {
        pushViewport(vp)
        colIndex <- (curDepth - 1) %% length(col) + 1
        fillIndex <- (curDepth - 1) %% length(fill) + 1
        grid.rect(gp=gpar(col=col[colIndex], fill=fill[fillIndex]))
        if (label)
            labelVP(vp, col[colIndex])
        upViewport()
    }
}

drawVP.vpPath <- function(vp, curDepth, depth, col, fill, label) {
    if (is.null(depth) || curDepth %in% depth) {
        downViewport(vp)
        colIndex <- (curDepth - 1) %% length(col) + 1
        fillIndex <- (curDepth - 1) %% length(fill) + 1
        grid.rect(gp=gpar(col=col[colIndex], fill=fill[fillIndex]))
        if (label)
            labelVP(vp, col[colIndex])
        upViewport(depth(vp))
    }
}

drawVP.vpList <- function(vp, curDepth, depth, col, fill, label) {
    lapply(vp, drawVP, curDepth, depth, col, fill, label)
}

drawVP.vpStack <- function(vp, curDepth, depth, col, fill, label) {
    d <- depth(vp)
    for (i in 1:length(vp)) {
        this <- vp[[i]]
        drawVP(this, curDepth, depth, col, fill, label)
        curDepth <- curDepth + depth(this)
        pushViewport(this)
    }
    upViewport(d)
}

drawVP.vpTree <- function(vp, curDepth, depth, col, fill, label) {
    if (vp$parent$name == "ROOT") {
        lapply(vp$children, drawVP, curDepth, depth, col, fill, label)
    } else {
        pushViewport(vp$parent)
        if (is.null(depth) || curDepth %in% depth) {
            colIndex <- (curDepth - 1) %% length(col) + 1
            fillIndex <- (curDepth - 1) %% length(fill) + 1
            grid.rect(gp=gpar(col=col[colIndex], fill=fill[fillIndex]))
            if (label) {
                drawLabel <- is.null(vp$children) ||
                             (!is.null(depth) &&
                              curDepth == max(depth))
                if (drawLabel)
                    labelVP(vp$parent, col[colIndex])
            }
        }
        lapply(vp$children, drawVP, curDepth + 1, depth, col, fill, label)
        upViewport()
    }
}

# Draw all viewports in same viewport
showVP <- function(vp, newpage, cvpt, depth, col, fill,
                   label) {
    # If we've started a new page, we'll need the old
    # viewport tree to navigate within
    if (newpage) {
        pushViewport(cvpt)
        # "-1" for "ROOT"
        upViewport(depth(cvpt) - 1)
    }
    # Work off a vpTree, so convert vp if it's a vpPath
    showingPath <- inherits(vp, "vpPath")
    if (showingPath) {
        path <- vp
        downViewport(path)
        vp <- current.vpTree(all=FALSE)
        upViewport(1)
    }
    drawVP(vp, 1, depth, col, fill, label)
    if (showingPath)
        # "-1" because we went down the path then back up 1 originally
        upViewport(depth(path) - 1)
    invisible()
}

# Convert a "viewport" to a set of vpPaths
leafPaths <- function(vp) {
    UseMethod("leafPaths")
}

leafPaths.viewport <- function(vp) {
    if (vp$name == "ROOT")
        NULL
    else
        vp$name
}

leafPaths.vpList <- function(vp) {
    unlist(lapply(vp, leafPaths))
}

leafPaths.vpStack <- function(vp) {
    pathList <- lapply(vp, leafPaths)
    for (i in 1:length(pathList)) {
        if (i > 1) {
            pathList[[i]] <- paste(pathList[[i - 1]],
                                   pathList[[i]],
                                   sep=.grid.pathSep)
        }
    }
    unlist(pathList)
}

leafPaths.vpTree <- function(vp) {
    if (is.null(vp$children)) {
        if (vp$parent$name == "ROOT")
            NULL
        else
            vp$parent$name
    } else {
        pathList <- lapply(vp$children, leafPaths)
        if (vp$parent$name == "ROOT") {
            unlist(pathList)
        } else {
            paste(vp$parent$name,
                  unlist(pathList),
                  sep=.grid.pathSep)
        }
    }
}

leafPaths.vpPath <- function(vp) {
    as.character(vp)
}

# Draw a vpPath
drawPath <- function(path, depth, col, fill, label) {
    n <- depth(path)
    for (i in 1:n) {
        downViewport(path[i])
        if (is.null(depth) || i %in% depth) {
            colIndex <- (i - 1) %% length(col) + 1
            fillIndex <- (i - 1) %% length(fill) + 1
            grid.rect(gp=gpar(col=col[colIndex], fill=fill[fillIndex]))
            if (label) {
                if (is.null(depth))
                    drawLabel <- i == n
                else
                    drawLabel <- i == min(n, max(depth))
                if (drawLabel)
                    labelVP(current.viewport(), col[colIndex])
            }
        }
    }
    upViewport(n)
}

# Draw each leaf in separate viewports
# FIXME: allow control over number of rows and cols
# NOTE: this does NOT leave its viewports hanging around after
showVPmatrix <- function(vp, cvpt, depth, col, fill,
                         label, # Only the leaf viewports are labelled
                         nrow, ncol) {
    # Work off a vpPath, so convert vp if it's a "viewport"
    if (is.viewport(vp)) {
        paths <- leafPaths(vp)
    } else {
        # Should not happen
        stop("how did we get here?")
    }
    firstPath <- 0
    while (length(paths) - firstPath > 0) {
        if (firstPath > 0)
            grid.newpage()
        pushViewport(viewport(layout=grid.layout(nrow, ncol)))
        for (i in 1:nrow) {
            for (j in 1:ncol) {
                theLeaf <- firstPath + (i - 1)*nrow + j
                if (theLeaf <= length(paths)) {
                    thePath <- vpPathDirect(paths[theLeaf])
                    pushViewport(viewport(layout.pos.row=i,
                                          layout.pos.col=j))
                    grid.rect(gp=gpar(col="grey80"))
                    # We may need the old vpTree to navigate within
                    # if 'vp' is a vpStack, or something similar, that
                    # contains a vpPath
                    if (!is.null(cvpt$children)) {
                        pushViewport(cvpt$children)
                        upViewport(depth(cvpt) - 1)
                    }
                    # Now push the viewport we are showing
                    pushViewport(vp)
                    upViewport(depth(vp))
                    # Now go to the particular viewport we
                    # are going to show
                    drawPath(thePath, depth, col, fill, label)
                    # Pop our placement within the layout
                    popViewport()
                }
            }
        }
        popViewport()
        firstPath <- firstPath + nrow*ncol
    }
}

showViewport <- function(vp=NULL, recurse=TRUE, depth=NULL,
                         newpage=FALSE, leaves=FALSE,
                         col=rgb(0, 0, 1, .2), fill=rgb(0, 0, 1, .1),
                         label=TRUE, nrow=3, ncol=nrow) {
    cvpt <- current.vpTree()
    if (is.null(vp))
        vp <- cvpt
    if (newpage == FALSE && leaves == TRUE)
        stop("must start new page if showing leaves separately")
    if (newpage) {
        grid.newpage()
    }
    if (!recurse)
        depth <- 1
    if (leaves) {
        # Special case of showing vpPath (i.e., only one viewport)
        # Ignores nrow & ncol
        if (inherits(vp, "vpPath"))
            showVP(vp, TRUE, cvpt, depth, col, fill, label)
        else
            showVPmatrix(vp, cvpt, depth, col, fill, label, nrow, ncol)
    } else {
        showVP(vp, newpage, cvpt, depth, col, fill, label)
    }
    invisible()
}
#  File src/library/grid/R/edit.R
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

# All args just used as args to editGrob
gEdit <- function(...) {
  edit <- list(...)
  class(edit) <- "gEdit"
  edit
}

is.gEdit <- function(x) {
  inherits(x, "gEdit")
}

applyEdit <- function(x, edit) {
  if (is.null(edit)) {
    x
  } else {
    if (!is.gEdit(edit))
      stop("invalid 'edit' information")
    # Intended to handle whether edit has gPath spec or not
    newx <- do.call("editGrob", c(list(x), edit))
    # If edit was specified for non-existent child, newx will be NULL
    if (is.null(newx))
      x
    else
      newx
  }
}

# A list of gEdit's to apply to the same grob
gEditList <- function(...) {
  edits <- list(...)
  if (!all(sapply(edits, is.gEdit)))
    stop("'gEditList' can only contain 'gEdit' objects")
  class(edits) <- "gEditList"
  edits
}

is.gEditList <- function(x) {
  inherits(x, "gEditList")
}

applyEdits <- function(x, edits) {
  if (is.null(edits)) {
    x
  } else {
    if (is.gEdit(edits))
      applyEdit(x, edits)
    else {
      if (!inherits(edits, "gEditList"))
        stop("invalid 'edit' information")
      for (i in edits)
        x <- applyEdits(x, i)
      x
    }
  }
}

#  File src/library/grid/R/frames.R
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

################
# frame class
################
# NOTE: make framevp separate slot (rather than combining with
# normal vp slot) so that it can be edited (e.g., by grid.pack)
frameGrob <- function(layout=NULL, name=NULL, gp=gpar(), vp=NULL) {
  if (!is.null(layout))
    framevp <- viewport(layout=layout)
  else
    framevp <- NULL
  gTree(framevp=framevp, name=name, gp=gp, vp=vp,
        cl="frame")
}

# draw=TRUE will not draw anything, but will mean that
# additions to the frame are drawn
grid.frame <- function(layout=NULL, name=NULL, gp=gpar(), vp=NULL,
                       draw=TRUE) {
  fg <- frameGrob(layout=layout, name=name, gp=gp, vp=vp)
  if (draw)
    grid.draw(fg)
  invisible(fg)
}

makeContext.frame <- function(x) {
    if (!is.null(x$framevp)) {
        if (!is.null(x$vp)) {
            x$vp <- vpStack(x$vp, x$framevp)
        } else {
            x$vp <- x$framevp
        }
    }
    x
}

widthDetails.frame <- function(x) {
  if (is.null(x$framevp))
    unit(1, "null")
  else
    sum(layout.widths(viewport.layout(x$framevp)))
}

heightDetails.frame <- function(x) {
  if (is.null(x$framevp))
    unit(1, "null")
  else
    sum(layout.heights(viewport.layout(x$framevp)))
}

frameDim <- function(frame) {
  if (is.null(frame$framevp))
    rep(0, 2)
  else
    c(layout.nrow(viewport.layout(frame$framevp)),
      layout.ncol(viewport.layout(frame$framevp)))
}

################
# cellGrob class
################
cellViewport <- function(col, row, border) {
  vp <- viewport(layout.pos.col=col, layout.pos.row=row)
  if (!is.null(border))
    vp <- vpStack(vp,
                  viewport(x=border[2L],
                           y=border[1L],
                           width=unit(1, "npc") - sum(border[c(2,4)]),
                           height=unit(1, "npc") - sum(border[c(1,3)]),
                           just=c("left", "bottom")))
  vp
}

cellGrob <- function(col, row, border, grob, dynamic, vp) {
  gTree(col=col, row=row, border=border, dynamic=dynamic,
        children=gList(grob), cellvp=vp, cl="cellGrob")
}

makeContext.cellGrob <- function(x) {
    if (!is.null(x$cellvp)) {
        if (!is.null(x$vp)) {
            x$vp <- vpStack(x$vp, x$cellvp)
        } else {
            x$vp <- x$cellvp
        }
    }
    x
}

# For dynamically packed grobs, need to be able to
# recalculate cell sizes
widthDetails.cellGrob <- function(x) {
  if (x$dynamic)
    unit(1, "grobwidth", gPath(x$children[[1L]]$name))
  else
    unit(1, "grobwidth", x$children[[1L]])
}

heightDetails.cellGrob <- function(x) {
  if (x$dynamic)
    unit(1, "grobheight", gPath(x$children[[1L]]$name))
  else
    unit(1, "grobheight", x$children[[1L]])
}

################
# grid.place
################
# Place an object into an already existing cell of a frame ...
# ... for a grob on the display list
grid.place <- function(gPath, grob,
                       row=1, col=1,
                       redraw=TRUE) {
  grid.set(gPath,
           placeGrob(grid.get(gPath), grob, row, col),
           redraw)
}

# ... for a grob description
placeGrob <- function(frame, grob,
                      row=NULL, col=NULL) {
  if (!inherits(frame, "frame"))
    stop("invalid 'frame'")
  if (!is.grob(grob))
    stop("invalid 'grob'")
  dim <- frameDim(frame)
  if (is.null(row))
    row <- c(1, dim[1L])
  if (is.null(col))
    col <- c(1, dim[2L])
  if (length(row) == 1)
    row <- rep(row, 2)
  if (length(col) == 1)
    col <- rep(col, 2)
  if (min(row) < 1 || max(row) > dim[1L] ||
      min(col) < 1 || max(col) > dim[2L])
    stop("invalid 'row' and/or 'col' (no such cell in frame layout)")
  cgrob <- cellGrob(col, row, NULL, grob, FALSE,
                    cellViewport(col, row, NULL))
  addGrob(frame, cgrob)
}

################
# grid.pack
################
num.col.specs <- function(side, col, col.before, col.after) {
  4 - sum(is.null(side) || any(c("top", "bottom") %in% side),
          is.null(col), is.null(col.before), is.null(col.after))
}

# We are assuming that checking has been done so that only one
# of these specifications has been given
col.spec <- function(side, col, col.before, col.after, ncol) {
  if (!is.null(side)) {
    if (side == "left")
      col <- 1
    else if (side == "right")
      col <- ncol + 1
  }
  else if (!is.null(col.before))
    col <- col.before
  else if (!is.null(col.after))
    col <- col.after + 1
  col
}

# We are assuming that checking has been done so that only one
# of these specifications has been given
new.col <- function(side, col, col.before, col.after, ncol) {
  # Special case ncol==0 for first grob added to frame
  result <- TRUE
  if (!is.null(col)) {
    # It is an error to specify a range for col which is outside 1..ncol
    if (length(col) == 2)
      if (col[1L] < 1 || col[2L] > ncol)
        stop("'col' can only be a range of existing columns")
      else
        result <- FALSE
    # It is also an error to specify a single col outside 1..ncol+1
    else
      if (col < 1 || col > ncol + 1)
        stop("invalid 'col' specification")
      else
        result <- col == ncol+1
  }
  result
}

num.row.specs <- function(side, row, row.before, row.after) {
  4 - sum(is.null(side) || any(c("left", "right") %in% side),
          is.null(row), is.null(row.before), is.null(row.after))
}

# We are assuming that checking has been done so that only one
# of these specifications has been given
row.spec <- function(side, row, row.before, row.after, nrow) {
  if (!is.null(side)) {
    if (side == "top")
      row <- 1
    else if (side == "bottom")
      row <- nrow + 1
  }
  else if (!is.null(row.before))
    row <- row.before
  else if (!is.null(row.after))
    row <- row.after + 1
  row
}

# We are assuming that checking has been done so that only one
# of these specifications has been given
new.row <- function(side, row, row.before, row.after, nrow) {
  # Special case nrow==0 for first grob added to frame
  result <- TRUE
  if (!is.null(row)) {
    # It is an error to specify a range for row which is outside 1..nrow
    if (length(row) == 2)
      if (row[1L] < 1 || row[2L] > nrow)
        stop("'row' can only be a range of existing rows")
      else
        result <- FALSE
    # It is also an error to specify a single row outside 1..nrow+1
    else
      if (row < 1 || row > nrow + 1)
        stop("invalid 'row' specification")
      else
        result <- row == nrow+1
  }
  result
}

mod.dims <- function(dim, dims, index, new.index, nindex, force) {
  # If adding a new row/col, add the new width/height to the list
  if (new.index)
    if (index == 1)
      dims <- unit.c(dim, dims)
    else if (index == nindex)
      dims <- unit.c(dims, dim)
    else
      dims <- unit.c(dims[1L:(index-1)], dim, dims[index:nindex])
  # Otherwise, if force=TRUE, we override previous width/heights for the
  # row/col, otherotherwise, the width/height of the existing row/col
  # is the maximum of the previous width/height and the new width/height
  else {
    if (!force)
      dim <- max(dim, dims[index])
    if (index==1)
      if (nindex == 1)
        dims <- dim
      else
        dims <- unit.c(dim, dims[2:nindex])
    else if (index==nindex)
      dims <- unit.c(dims[1L:(nindex-1)], dim)
    else
      dims <- unit.c(dims[1L:(index-1)], dim, dims[(index+1):nindex])
  }
  dims
}

updateCol <- function(col, added.col) {
  old.col <- col
  # If grob$col is a range ...
  if (length(old.col) == 2) {
    if (added.col <= old.col[2L])
      col <- c(old.col[1L], old.col[2L] + 1)
  }
  else
    if (added.col <= old.col)
      col <- old.col + 1
  col
}

updateRow <- function(row, added.row) {
  old.row <- row
  # If grob$row is a range ...
  if (length(old.row) == 2) {
    if (added.row <= old.row[2L])
      row <- c(old.row[1L], old.row[2L] + 1)
  }
  else
    if (added.row <= old.row)
      row <- old.row + 1
  row
}

# FIXME:  Allow specification of respect for new row/col
# Pack a child grob within a frame grob ...
# (a special sort of editing just for frame grobs)
# ... for a grob on the display list
grid.pack <- function(gPath, grob, redraw=TRUE,
                      side=NULL,
                      row=NULL, row.before=NULL, row.after=NULL,
                      col=NULL, col.before=NULL, col.after=NULL,
                      width=NULL, height=NULL,
                      force.width=FALSE, force.height=FALSE,
                      border=NULL, dynamic=FALSE) {
  grid.set(gPath,
           packGrob(grid.get(gPath), grob, side,
                    row, row.before, row.after,
                    col, col.before, col.after,
                    width, height, force.width, force.height,
                    border),
           redraw)
}

packGrob <- function(frame, grob,
                     side=NULL,
                     row=NULL, row.before=NULL, row.after=NULL,
                     col=NULL, col.before=NULL, col.after=NULL,
                     width=NULL, height=NULL,
                     force.width=FALSE, force.height=FALSE,
                     border=NULL, dynamic=FALSE) {
  if (!inherits(frame, "frame"))
    stop("invalid 'frame'")
  if (!is.grob(grob))
    stop("invalid 'grob'")
  # col/row can be given as a range, but I only want to know
  # about the min and max
  if (!is.null(col) & length(col) > 1) {
    col <- range(col)
    col.range <- TRUE
  }
  else
    col.range <- FALSE
  if (!is.null(row) & length(row) > 1) {
    row <- range(row)
    row.range <- TRUE
  }
  else
    row.range <- FALSE

  frame.vp <- frame$framevp
  if (is.null(frame.vp))
    frame.vp <- viewport()
  lay <- viewport.layout(frame.vp)
  if (is.null(lay)) {
    ncol <- 0
    nrow <- 0
  } else {
    ncol <- layout.ncol(lay)
    nrow <- layout.nrow(lay)
  }

  # (i) Check that the specifications of the location of the grob
  # give a unique location
  ncs <- num.col.specs(side, col, col.before, col.after)
  # If user does not specify a col, assume it is all cols
  if (ncs == 0) {
    # Allow for fact that this might be first grob packed
    if (ncol > 0) {
      col <- c(1, ncol)
      col.range <- TRUE
    }
    else
      col <- 1
    ncs <- 1
  }
  if (ncs != 1)
    stop("cannot specify more than one of 'side=[\"left\", \"right\"]', 'col', 'col.before', or 'col.after'")
  nrs <- num.row.specs(side, row, row.before, row.after)
  # If user does not specify a row, assume it is all rows
  if (nrs == 0) {
    # Allow for fact that this might be first grob packed
    if (nrow > 0) {
      row <- c(1, nrow)
      row.range <- TRUE
    }
    else
      row <- 1
    nrs <- 1
  }
  if (nrs != 1)
    stop("must specify exactly one of 'side=[\"top\", \"bottom\"]', 'row', 'row.before', or 'row.after'")

  # (ii) Determine that location and check that it is valid
  new.col <- new.col(side, col, col.before, col.after, ncol)
  col <- col.spec(side, col, col.before, col.after, ncol)
  new.row <- new.row(side, row, row.before, row.after, nrow)
  row <- row.spec(side, row, row.before, row.after, nrow)

  # Wrap the child in a "cellGrob" to maintain additional info
  # (like row and col occupied in frame)
  # Need to do this here so can create widths/heights based on this cell grob
  if (!is.null(grob))
    cgrob <- cellGrob(col, row, border, grob, dynamic,
                      cellViewport(col, row, border))

  # (iii) If width and height are not given, take them from the child
  #       NOTE:  if dynamic is TRUE then use a gPath to the child
  if (is.null(width))
    if (is.null(grob))
      width <- unit(1, "null")
    else
      if (dynamic)
        width <- unit(1, "grobwidth", gPath(cgrob$name))
      else
        width <- unit(1, "grobwidth", cgrob)
  if (is.null(height))
    if (is.null(grob))
      height <- unit(1, "null")
    else
      if (dynamic)
        height <- unit(1, "grobheight", gPath(cgrob$name))
      else
        height <- unit(1, "grobheight", cgrob)
  # If there is a border, include it in the width/height
  if (!is.null(border)) {
    width <- sum(border[2L], width, border[4L])
    height <- sum(border[1L], height, border[3L])
  }

  # (iv) Update the frame.vp of the frame (possibly add new row/col,
  # possibly update existing widths/heights and respect)
  if (new.col) ncol <- ncol + 1
  if (new.row) nrow <- nrow + 1
  # If we are creating the frame.vp$layout for the first time then
  # we have to initialise the layout widths and heights
  if (is.null(lay)) {
    widths <- width
    heights <- height
  } else {
    # DO NOT modify widths/heights if the grob is being added to
    # multiple columns/rows
    if (col.range)
      widths <- layout.widths(lay)
    else
      widths <- mod.dims(width, layout.widths(lay), col, new.col, ncol,
                         force.width)
    if (row.range)
      heights <- layout.heights(lay)
    else
      heights <- mod.dims(height, layout.heights(lay), row, new.row, nrow,
                          force.height)
  }
  frame.vp$layout <- grid.layout(ncol=ncol, nrow=nrow,
                                 widths=widths, heights=heights)

  # Modify the locations (row, col) of existing children in the frame
  if (new.col || new.row) {
    for (i in childNames(frame)) {
      child <- getGrob(frame, i)
      if (new.col) {
        newcol <- updateCol(child$col, col)
        child <- editGrob(child, col=newcol,
                          cellvp=cellViewport(newcol, child$row, child$border))
      }
      if (new.row) {
        newrow <- updateRow(child$row, row)
        child <- editGrob(child, row=newrow,
                          cellvp=cellViewport(child$col, newrow, child$border))
      }
      frame <- addGrob(frame, child)
    }
  }

  # Add the new grob to the frame
  if (!is.null(grob)) {
    frame <- addGrob(frame, cgrob)
  }

  editGrob(frame, framevp=frame.vp)
}

#  File src/library/grid/R/function.R
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


validDetails.functiongrob <- function(x, ...) {
    if (x$n < 1)
        stop("invalid 'n'")
    if (!(is.character(x$range) && x$range %in% c("x", "y")))
        x$range <- as.numeric(x$range)
    if (!is.function(x$f))
        stop("invalid 'f'")
    x
}

genXY <- function(x) {
    if (is.numeric(x$range)) {
        range <- range(x$range)
    } else {
        if (x$range == "x")
            range <- current.viewport()$xscale
        else
            range <- current.viewport()$yscale
    }
    input <- seq(range[1], range[2], length.out=x$n)
    x$f(input)
}

makeContent.functiongrob <- function(x) {
    xy <- genXY(x)
    linesGrob(xy$x, xy$y, default.units=x$units,
              name=x$name, gp=x$gp, vp=x$vp)
}

xDetails.functiongrob <- function(x, theta) {
    xy <- genXY(x)
    xDetails(linesGrob(xy$x, xy$y, default.units=x$units), theta)
}

yDetails.functiongrob <- function(x, theta) {
    xy <- genXY(x)
    yDetails(linesGrob(xy$x, xy$y, default.units=x$units), theta)
}

widthDetails.functiongrob <- function(x) {
    xy <- genXY(x)
    widthDetails(linesGrob(xy$x, xy$y, default.units=x$units))
}

heightDetails.functiongrob <- function(x) {
    xy <- genXY(x)
    heightDetails(linesGrob(xy$x, xy$y, default.units=x$units))
}

functionGrob <- function(f, n=101, range="x", units="native",
                         name=NULL, gp=gpar(), vp=NULL) {
    grob(f=f, range=range, units=units, n=n,
         gp=gp, vp=vp, name=name, cl="functiongrob")
}

grid.function <- function(...) {
    grid.draw(functionGrob(...))
}

# Convenience wrappers
grid.abline <- function(intercept=0, slope=1, ...) {
    grid.function(function(x) list(x=x, y=intercept + slope*x), ...)
}

##############
# Tests
tests <- function() {

    # editing
    grid.newpage()
    pushViewport(viewport(xscale=c(0, 2*pi), yscale=c(-1, 1)))
    grid.function(function(x) list(x=x, y=sin(x)), name="fg")
    grid.edit("fg", n=10)
    grid.edit("fg", f=function(x) list(x=x, y=x))

    # x/y/width/height calculations
    grid.newpage()
    pushViewport(viewport(xscale=c(0, 2*pi), yscale=c(-2, 2)))
    grid.function(function(x) list(x=x, y=sin(x)), name="fg")
    grid.segments(0, 1,
                  grobX("fg", 135), grobY("fg", 135),
                  arrow=arrow())

}
#  File src/library/grid/R/gpar.R
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


# A "gpar" object is a list of graphics parameters
# A graphics parameter is a name-value pair

gpar <- function(...) {
  gp <- validGP(list(...))
  class(gp) <- "gpar"
  gp
}

is.gpar <- function(x) {
  inherits(x, "gpar")
}

print.gpar <- function(x, ...) {
  print(unclass(x), ...)
  invisible(x)
}

validGP <- function(gpars) {
  # Check a (non-NULL) gpar is not of length 0
  check.length <- function(gparname) {
    if (length(gpars[[gparname]]) == 0)
      stop(gettextf("'gpar' element '%s' must not be length 0", gparname),
           domain = NA)
  }
  # Check a gpar is numeric and not NULL
  numnotnull <- function(gparname) {
    if (!is.na(match(gparname, names(gpars)))) {
      if (is.null(gpars[[gparname]]))
        gpars[[gparname]] <<- NULL
      else {
        check.length(gparname)
        gpars[[gparname]] <<- as.numeric(gpars[[gparname]])
      }
    }
  }
  # fontsize, lineheight, cex, lwd should be numeric and not NULL
  numnotnull("fontsize")
  numnotnull("lineheight")
  numnotnull("cex")
  numnotnull("lwd")
  numnotnull("lex")
  # gamma defunct in 2.7.0
  if ("gamma" %in% names(gpars)) {
    warning("'gamma' 'gpar' element is defunct")
    gpars$gamma <- NULL
  }
  numnotnull("alpha")
  # col and fill are converted in C code
  # BUT still want to check length > 0
  if (!is.na(match("col", names(gpars)))) {
      if (is.null(gpars$col))
          gpars$col <- NULL
      else
          check.length("col")
  }
  if (!is.na(match("fill", names(gpars)))) {
      if (is.null(gpars$fill))
          gpars$fill <- NULL
      else
          check.length("fill")
  }
  # lty converted in C code
  # BUT still want to check for NULL and check length > 0
  if (!is.na(match("lty", names(gpars)))) {
    if (is.null(gpars$lty))
      gpars$lty <- NULL
    else
      check.length("lty")
  }
  if (!is.na(match("lineend", names(gpars)))) {
    if (is.null(gpars$lineend))
      gpars$lineend <- NULL
    else
      check.length("lineend")
  }
  if (!is.na(match("linejoin", names(gpars)))) {
    if (is.null(gpars$linejoin))
      gpars$linejoin <- NULL
    else
      check.length("linejoin")
  }
  # linemitre should be larger than 1
  numnotnull("linemitre")
  if (!is.na(match("linemitre", names(gpars)))) {
    if (any(gpars$linemitre < 1))
      stop("invalid 'linemitre' value")
  }
  # alpha should be 0 to 1
  if (!is.na(match("alpha", names(gpars)))) {
    if (any(gpars$alpha < 0 || gpars$alpha > 1))
      stop("invalid 'alpha' value")
  }
  # font should be integer and not NULL
  if (!is.na(match("font", names(gpars)))) {
    if (is.null(gpars$font))
      gpars$font <- NULL
    else {
      check.length("font")
      gpars$font <- as.integer(gpars$font)
    }
  }
  # fontfamily should be character
  if (!is.na(match("fontfamily", names(gpars)))) {
    if (is.null(gpars$fontfamily))
      gpars$fontfamily <- NULL
    else {
      check.length("fontfamily")
      gpars$fontfamily <- as.character(gpars$fontfamily)
    }
  }
  # fontface can be character or integer;  map character to integer
  # store value in font
  # Illegal to specify both font and fontface
  if (!is.na(match("fontface", names(gpars)))) {
    if (!is.na(match("font", names(gpars))))
      stop("must specify only one of 'font' and 'fontface'")
    if (is.null(gpars$fontface))
      gpars$font <- NULL
    else {
      check.length("fontface")
      if (is.numeric(gpars$fontface))
        gpars$font <- as.integer(gpars$fontface)
      else {
        temp.char <- as.character(gpars$fontface)
        temp.num <- integer(length(temp.char))
        for (i in seq_along(temp.char))
          temp.num[i] <- switch(temp.char[i],
                                plain=1,
                                italic=3,
                                oblique=3,
                                bold=2,
                                bold.italic=4,
                                symbol=5,
                                # These are Hershey variants
                                cyrillic=5,
                                cyrillic.oblique=6,
                                EUC=7,
                                stop("invalid font face"))
        gpars$font <- as.integer(temp.num)
      }
    }
  }
  gpars
}

# Method for subsetting "gpar" objects
`[.gpar` <- function(x, index, ...) {
    if (length(x) == 0)
        return(gpar())
    maxn <- do.call("max", lapply(x, length))
    newgp <- lapply(x, rep, length.out=maxn)
    newgp <- lapply(X = newgp, FUN = "[", index, ...)
    class(newgp) <- "gpar"
    newgp
}

# possible gpar names
# The order must match the GP_* values in grid.h
.grid.gpar.names <- c("fill", "col", "gamma", "lty", "lwd", "cex",
                      "fontsize", "lineheight", "font", "fontfamily",
                      "alpha", "lineend", "linejoin", "linemitre",
                      "lex",
                      # Keep fontface at the end because it is never
                      # used in C code (it gets mapped to font)
                      "fontface")

set.gpar <- function(gp) {
  if (!is.gpar(gp))
    stop("argument must be a 'gpar' object")
  temp <- grid.Call(L_getGPar)
  # gamma defunct in 2.7.0
  if ("gamma" %in% names(gp)) {
      warning("'gamma' 'gpar' element is defunct")
      gp$gamma <- NULL
  }
  # Special case "cex" (make it cumulative)
  if (match("cex", names(gp), nomatch=0L))
    tempcex <- temp$cex * gp$cex
  else
    tempcex <- temp$cex
  # Special case "alpha" (make it cumulative)
  if (match("alpha", names(gp), nomatch=0L))
    tempalpha <- temp$alpha * gp$alpha
  else
    tempalpha <- temp$alpha
  # Special case "lex" (make it cumulative)
  if (match("lex", names(gp), nomatch=0L))
    templex <- temp$lex * gp$lex
  else
    templex <- temp$lex
  # All other gpars
  temp[names(gp)] <- gp
  temp$cex <- tempcex
  temp$alpha <- tempalpha
  temp$lex <- templex
  # Do this as a .Call.graphics to get it onto the base display list
  grid.Call.graphics(L_setGPar, temp)
}

get.gpar <- function(names=NULL) {
  if (is.null(names)) {
    result <- grid.Call(L_getGPar)
    # drop gamma
    result$gamma <- NULL
  } else {
    if (!is.character(names) ||
        !all(names %in% .grid.gpar.names))
      stop("must specify only valid 'gpar' names")
    # gamma deprecated
    if ("gamma" %in% names) {
      warning("'gamma' 'gpar' element is defunct")
      names <- names[-match("gamma", names)]
    }
    result <- unclass(grid.Call(L_getGPar))[names]
  }
  class(result) <- "gpar"
  result
}

# When editing a gp slot, only update the specified gpars
# Assume gp is NULL or a gpar
# assume newgp is a gpar (and not NULL)
mod.gpar <- function(gp, newgp) {
  if (is.null(gp))
    gp <- newgp
  else
    gp[names(newgp)] <- newgp
  gp
}

#  File src/library/grid/R/grab.R
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

#########
# Generate a gTree from the current display list
#
# Or from an expression
# (recording on to a null graphics device)
#########
rootVP <- function(pvp) {
  match(pvp$name, "ROOT", nomatch=FALSE)
}

# List the children of the current vp (as a vpList)
current.vpList <- function() {
  cpvp <- grid.Call(L_currentViewport)
  if (length(ls(cpvp$children, all.names=TRUE)) == 0)
    NULL
  else
    vpListFromNode(cpvp)
}

current.vpNames <- function() {
  ls(grid.Call(L_currentViewport)$children)
}

# vp might be a viewport, or a vpList, or a vpStack, or a vpTree
vpExists <- function(vp) {
  UseMethod("vpExists")
}

vpExists.viewport <- function(vp) {
  vp$name %in% ls(.Call(L_currentViewport)$children)
}

vpExists.vpStack <- function(vp) {
  vpExists(vp[[1L]])
}

vpExists.vpList <- function(vp) {
  any(sapply(vp, vpExists, simplify=TRUE))
}

vpExists.vpTree <- function(vp) {
  vpExists(vp$parent)
}

# Handle vpPaths in a vpStack or vpTree
# Not a problem to downViewport() to a viewport that already exists
vpExists.vpPath <- function(vp) {
    FALSE
}

wrap <- function(x) {
  UseMethod("wrap")
}

wrap.default <- function(x) {
  if (!is.null(x))
    stop("invalid display list element")
  NULL
}

wrap.grob <- function(x) {
  x
}

wrap.viewport <- function(x) {
  recordGrob(pushViewport(vp), list(vp=x))
}

wrap.pop <- function(x) {
  recordGrob(popViewport(n), list(n=x))
}

wrap.up <- function(x) {
  recordGrob(upViewport(n), list(n=x))
}

wrap.vpPath <- function(x) {
  recordGrob(downViewport(path), list(path=x))
}

# Grab the display list on the current device
# ... are passed to gTree
# If warn is 0, issue no warnings
# If warn is 1, issue warnings about situations that are definitely
#   NOT captured correctly (e.g., reuse of top-level grob name)
# If warn is 2, issue warnings about situations that
#   MAY not get captured correctly (e.g., top-level downViewport)
# If wrap is TRUE, grab will wrap all pushes and grobs
#   in a gTree
grabDL <- function(warn, wrap, ...) {
  gList <- NULL
  dl.index <- grid.Call(L_getDLindex)
  if (dl.index > 1) {
    if (warn > 0) {
      names <- getNames()
      # Check for overwriting existing grob
      if (length(unique(names)) != length(names))
        warning("one of more grobs overwritten (grab WILL not be faithful; try 'wrap = TRUE')")
    }
    grid.newpage(recording=FALSE)
    # Start at 2 because first element is viewport[ROOT]
    for (i in 2:dl.index) {
      # Do all of this as a big ifelse rather than
      # dispatching to a function call per element because
      # we need to work with whole DL at times, not
      # just individual elements
      elt <- grid.Call(L_getDLelt, as.integer(i - 1))
      if (wrap)
        gList <- addToGList(wrap(elt), gList)
      else {

        ###########
        # grabGrob
        ###########
        if (inherits(elt, "grob")) {
          # Enforce grob$vp now and set grob$vp to NULL
          # Will be replaced later with full vpPath
          tempvp <- elt$vp
          if (warn > 1) {
            # Check to see if about to push a viewport
            # with existing viewport name
            if (inherits(tempvp, "viewport") &&
                vpExists(tempvp))
              warning("viewport overwritten (grab MAY not be faithful)")
          }
          if (!is.null(tempvp))
            tempdepth <- depth(tempvp)
          grid.draw(tempvp, recording=FALSE)
          # vpPath after grob$vp slot has been pushed
          # Has to be recorded here in case grob drawing
          # pushes (and does not pop) more viewports
          drawPath <- current.vpPath()
          elt$vp <- NULL
          grid.draw(elt, recording=FALSE)
          if (warn > 1) {
            # Compare new vpPath
            # If not same, the grob has pushed some viewports
            # and not popped or upped them
            pathSame <- TRUE
            if (!(is.null(drawPath) && is.null(current.vpPath()))) {
              if (is.null(drawPath))
                pathSame = FALSE
              else if (is.null(current.vpPath()))
                pathSame = FALSE
              else if (as.character(drawPath) !=
                       as.character(current.vpPath()))
                pathSame = FALSE
            }
            if (!pathSame)
              warning("grob pushed viewports and did not pop/up them (grab MAY not be faithful)")
          }
          elt$vp <- drawPath
          if (!is.null(tempvp))
            upViewport(tempdepth, recording=FALSE)
          gList <- addToGList(elt, gList)
        ###########
        # grabViewport
        ###########
        } else if (inherits(elt, "viewport")) {
          # Includes viewports, vpLists, vpTrees, and vpStacks
          # Check to see if about to push a viewport
          # with existing viewport name
          if (vpExists(elt))
            warning("viewport overwritten (grab MAY not be faithful)")
          grid.draw(elt, recording=FALSE)
        ###########
        # grabPop
        ###########
        } else if (inherits(elt, "pop")) {
          # Replace pop with up
          upViewport(elt, recording=FALSE)

        ###########
        # grabDefault
        ###########
        } else {
          grid.draw(elt, recording=FALSE)
        }
      } # matches if (wrap)
    }
    # Go to top level
    upViewport(0, recording=FALSE)
    gTree(children=gList, childrenvp=current.vpList(), ...)
  } else {
    NULL
  }
}

# expr is ignored if dev is NULL
# otherwise, it should be an expression, like postscript("myfile.ps")
grid.grab <- function(warn=2, wrap=FALSE, ...) {
  grabDL(warn, wrap, ...)
}

grid.grabExpr <- function(expr, warn=2, wrap=FALSE, ...) {
    ## Start an "offline" PDF device for this function
    pdf(file=NULL)
    on.exit(dev.off())
    ## Run the graphics code in expr
    ## Rely on lazy evaluation for correct "timing"
    eval(expr)
    ## Grab the DL on the new device
    grabDL(warn, wrap, ...)
}

#########################
# A different sort of capture ...
# Just grab the screen raster image
#########################

grid.cap <- function() {
    # This does not need recording on the display list
    grid.Call(L_cap)
}


#  File src/library/grid/R/grid.R
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


# FIXME:  all grid functions should check that .grid.started is TRUE
.grid.loaded <- FALSE

push.vp <- function(vp, recording) {
  UseMethod("push.vp")
}

push.vp.default <- function(vp, recording) {
  stop("only valid to push viewports")
}

push.vp.viewport <- function(vp, recording) {
  # Record on the display list
  if (recording)
    record(vp)
  # Create a pushedvp object for the system to keep track of
  pvp <- pushedvp(vp)
  # Store the entire set of gpar settings JUST PRIOR to push
  # We refer to this when calculating the viewport transform
  # We cannot simply rely on parent's gpar because we may be
  # being pushed from within a gTree which has enforced gpar
  # settings (i.e., the gTree$gp is enforced between this viewport
  # and the this viewport's parent$gp)
  pvp$parentgpar <- grid.Call(L_getGPar)
  # Enforce gpar settings
  set.gpar(vp$gp)
  # Store the entire set of gpar settings for this viewport
  pvp$gpar <- grid.Call(L_getGPar)
  # Pass in the pushedvp structure which will be used to store
  # things like the viewport transformation, parent-child links, ...
  grid.Call.graphics(L_setviewport, pvp, TRUE)
}

# For all but the last viewport, push the
# viewport then pop it
# For the last viewport, just push
push.vp.vpList <- function(vp, recording) {
  push.vp.parallel <- function(vp, recording) {
    push.vp(vp, recording)
    upViewport(depth(vp), recording)
  }
  if (length(vp) == 1)
    push.vp(vp[[1L]], recording)
  else {
    lapply(vp[1L:(length(vp) - 1)], push.vp.parallel, recording)
    push.vp(vp[[length(vp)]], recording)
  }
}

# Push viewports in series
push.vp.vpStack <- function(vp, recording) {
  lapply(vp, push.vp, recording)
}

# Push parent
# Children are a vpList
push.vp.vpTree <- function(vp, recording) {
  # Special case if user has saved the entire vpTree
  # parent will be the ROOT viewport, which we don't want to
  # push (grid ensures it is ALWAYS there)
  if (!(vp$parent$name %in% "ROOT"))
    push.vp(vp$parent, recording)
  push.vp(vp$children, recording)
}

# "push"ing a vpPath is just a downViewport(..., strict=TRUE)
push.vp.vpPath <- function(vp, recording) {
    downViewport(vp, strict=TRUE, recording)
}

push.viewport <- function(..., recording=TRUE) {
  .Deprecated("pushViewport")
  pushViewport(..., recording=recording)
}

pushViewport <- function(..., recording=TRUE) {
  if (missing(...))
    stop("must specify at least one viewport")
  else {
    vps <- list(...)
    lapply(vps, push.vp, recording)
  }
  invisible()
}

# Helper functions called from C
no.children <- function(children) {
  length(ls(children, all.names=TRUE)) == 0
}

child.exists <- function(name, children) {
  exists(name, envir=children, inherits=FALSE)
}

child.list <- function(children) {
  ls(children, all.names=TRUE)
}

pathMatch <- function(path, pathsofar, strict) {
  if (is.null(pathsofar))
    is.null(path)
  else {
    pattern <- paste0(if(strict) "^", path, "$")
    grepl(pattern, pathsofar)
  }
}

growPath <- function(pathsofar, name) {
  paste(pathsofar, name, sep=.grid.pathSep)
}

# Rather than pushing a new viewport, navigate down to one that has
# already been pushed
downViewport <- function(name, strict=FALSE, recording=TRUE) {
  UseMethod("downViewport")
}

# For interactive use, allow user to specify
# vpPath directly (i.e., w/o calling vpPath)
downViewport.default <- function(name, strict=FALSE, recording=TRUE) {
  name <- as.character(name)
  downViewport(vpPathDirect(name), strict, recording=recording)
}

downViewport.vpPath <- function(name, strict=FALSE, recording=TRUE) {
    if (name$n == 1)
        result <- grid.Call.graphics(L_downviewport, name$name, strict)
    else
        result <- grid.Call.graphics(L_downvppath,
                                     name$path, name$name, strict)
    # If the downViewport() fails, there is an error in C code
    # so none of the following code will be run

    # Enforce the gpar settings for the viewport
    pvp <- grid.Call(L_currentViewport)
    # Do not call set.gpar because set.gpar accumulates cex
    grid.Call.graphics(L_setGPar, pvp$gpar)
    # Record the viewport operation
    # ... including the depth navigated down
    if (recording) {
        attr(name, "depth") <- result
        record(name)
    }
    invisible(result)
}

# Similar to down.viewport() except it starts searching from the
# top-level viewport, so the result may be "up" or even "across"
# the current viewport tree
seekViewport <- function(name, recording=TRUE) {
  # up to the top-level
  upViewport(0, recording=recording)
  downViewport(name, recording=recording)
}

# Depth of the current viewport
vpDepth <- function() {
  pvp <- grid.Call(L_currentViewport)
  count <- 0
  while (!is.null(pvp$parent)) {
    pvp <- pvp$parent
    count <- count + 1
  }
  count
}

pop.viewport <- function(n=1, recording=TRUE) {
  .Deprecated("popViewport")
  popViewport(n, recording=recording)
}

popViewport <- function(n=1, recording=TRUE) {
  if (n < 0)
    stop("must pop at least one viewport")
  if (n == 0)
    n <- vpDepth()
  if (n > 0) {
    grid.Call.graphics(L_unsetviewport, as.integer(n))
    # Record on the display list
    if (recording) {
      class(n) <- "pop"
      record(n)
    }
  }
  invisible()
}

# Rather than removing the viewport from the viewport stack (tree),
# simply navigate up, leaving pushed viewports in place.
upViewport <- function(n=1, recording=TRUE) {
  if (n < 0)
    stop("must navigate up at least one viewport")
  if (n == 0) {
    n <- vpDepth()
    upPath <- current.vpPath()
  }
  if (n > 0) {
    path <- current.vpPath()
    upPath <- path[(depth(path) - n + 1):depth(path)]
    grid.Call.graphics(L_upviewport, as.integer(n))
    # Record on the display list
    if (recording) {
      class(n) <- "up"
      record(n)
    }
  }
  invisible(upPath)
}

# Return the full vpPath to the current viewport
current.vpPath <- function() {
  names <- NULL
  pvp <- grid.Call(L_currentViewport)
  while (!rootVP(pvp)) {
    names <- c(names, pvp$name)
    pvp <- pvp$parent
  }
  if (!is.null(names))
    vpPathFromVector(rev(names))
  else
    names
}

# Function to obtain the current viewport
current.viewport <- function(vp=NULL) {
  if (is.null(vp))
    # The system stores a pushedvp;  the user should only
    # ever see normal viewports, so convert.
    vpFromPushedvp(grid.Call(L_currentViewport))
  else {
    warning("the 'vp' argument is deprecated")
    vp
  }
}

vpListFromNode <- function(node) {
  childnames <- ls(node$children, all.names=TRUE)
  n <- length(childnames)
  children <- vector("list", n)
  index <- 1
  for (i in childnames) {
    children[[index]] <- vpTreeFromNode(get(i, envir=node$children))
    index <- index + 1
  }
  vpListFromList(children)
}

vpTreeFromNode <- function(node) {
  # If no children then just return viewport
  if (length(ls(node$children, all.names=TRUE)) == 0)
    vpFromPushedvp(node)
  # Otherwise return vpTree
  else
    vpTree(vpFromPushedvp(node),
           vpListFromNode(node))
}

# Obtain the current viewport tree
# Either from the current location in the tree down
# or ALL of the tree
current.vpTree <- function(all=TRUE) {
  cpvp <- grid.Call(L_currentViewport)
  moving <- all && vpDepth() > 0
  if (moving) {
    savedname <- cpvp$name
    upViewport(0, recording=FALSE)
    cpvp <- grid.Call(L_currentViewport)
  }
  tree <- vpTreeFromNode(cpvp)
  if (moving) {
    downViewport(savedname, recording=FALSE)
  }
  tree
}

current.transform <- function() {
  grid.Call(L_currentViewport)$trans
}

# Call this function if you want the graphics device erased or moved
# on to a new page.  High-level plotting functions should call this.
# NOTE however, that if you write a function which calls grid.newpage,
# you should provide an argument to allow people to turn it off
# so that they can use your function within a parent viewport
# (rather than the whole device) if they want to.
grid.newpage <- function(recording=TRUE) {
    for (fun in getHook("before.grid.newpage"))  {
        if(is.character(fun)) fun <- get(fun)
        try(fun())
    }
    # NOTE that we do NOT do grid.Call here because we have to do
    # things slightly differently if grid.newpage is the first grid operation
    # on a new device
    .Call(L_newpagerecording)
    .Call(L_newpage)
    .Call(L_initGPar)
    .Call(L_initViewportStack)
    if (recording) {
        .Call(L_initDisplayList)
        grDevices:::recordPalette()
        for (fun in getHook("grid.newpage"))  {
            if(is.character(fun)) fun <- get(fun)
            try(fun())
        }
    }
    invisible()
}

###########
# DISPLAY LIST FUNCTIONS
###########

# Keep a list of all drawing operations (since last grid.newpage()) so
# that we can redraw upon edit.

inc.display.list <- function() {
  display.list <- grid.Call(L_getDisplayList)
  dl.index <- grid.Call(L_getDLindex)
  dl.index <- dl.index + 1
  n <- length(display.list)
  # The " - 1" below is because dl.index is now stored internally
  # so is a C-style zero-based index rather than an R-style
  # 1-based index
  if (dl.index > (n - 1)) {
    temp <- display.list
    display.list <- vector("list", n + 100L)
    display.list[1L:n] <- temp
  }
  grid.Call(L_setDisplayList, display.list)
  grid.Call(L_setDLindex, as.integer(dl.index))
}

# This will either ...
#   (i) turn on AND INITIALISE the display list or ...
#   (ii) turn off AND ERASE the display list
grid.display.list <- function(on=TRUE) {
  grid.Call(L_setDLon, as.logical(on))
  if (on) {
    grid.Call(L_setDisplayList, vector("list", 100L))
    grid.Call(L_setDLindex, 0L)
  }
  else
    grid.Call(L_setDisplayList, NULL)
}

record <- function(x) {
  if (grid.Call(L_getDLon))
    UseMethod("record")
}

# When there is a pop.viewport, the number of viewports popped
# gets put on the display list
record.default <- function(x) {
  if (!is.numeric(x))
    stop("invalid object inserted on the display list")
  grid.Call(L_setDLelt, x)
  inc.display.list()
}

record.grob <- function(x) {
  grid.Call(L_setDLelt, x)
  inc.display.list()
}

record.viewport <- function(x) {
  grid.Call(L_setDLelt, x)
  inc.display.list()
}

record.vpPath <- function(x) {
  grid.Call(L_setDLelt, x)
  inc.display.list()
}

# This controls whether grid is using the graphics engine's display list
engine.display.list <- function(on=TRUE) {
  grid.Call(L_setEngineDLon, as.logical(on))
}

# Rerun the grid DL
grid.refresh <- function() {
  draw.all()
}

# Call a function on each element of the grid display list
# AND replace the element with the result
# This is blood-curdlingly dangerous for the state of the
# display list
# Two token efforts at safety are made:
#   - generate all of the new elements first THEN assign them all
#     (so if there is an error in generating any one element
#      you don't end up with a trashed display list)
#   - check that the new element is either NULL or the same
#     class as the element it is replacing
grid.DLapply <- function(FUN, ...) {
    FUN <- match.fun(FUN)
    # Traverse DL and do something to each entry
    gridDL <- grid.Call(L_getDisplayList)
    gridDLindex <- grid.Call(L_getDLindex)
    newDL <- vector("list", gridDLindex)
    for (i in 1:(gridDLindex - 1)) {
        elt <- grid.Call(L_getDLelt, i)
        newElt <- FUN(elt, ...)
        if (!(is.null(newElt) || inherits(newElt, class(elt))))
            stop("invalid modification of the display list")
        newDL[[i]] <- newElt
    }
    for (i in 1:(gridDLindex - 1)) {
        grid.Call(L_setDLindex, i)
        grid.Call(L_setDLelt, newDL[[i]])
    }
    grid.Call(L_setDLindex, gridDLindex)
}

# Wrapper for .Call and .Call.graphics
# Used to make sure that grid-specific initialisation occurs just before
# the first grid graphics output OR the first querying of grid state
# (on the current device)
# The general rule is you should use these rather than .Call or
# .Call.graphics unless you have a good reason and you know what
# you are doing -- this will be a bit of overkill, but is for safety
grid.Call <- function(fnname, ...) {
  .Call(L_gridDirty)
  .Call(fnname, ...)
}

grid.Call.graphics <- function(fnname, ...) {
  # Only record graphics operations on the graphics engine's display
  # list if the engineDLon flag is set
  engineDLon <- grid.Call(L_getEngineDLon)
  if (engineDLon) {
    # NOTE that we need a .Call.graphics("L_gridDirty") so that
    # the the first thing on the engine display list is a dirty
    # operation;  this is necessary in case the display list is
    # played on another device (e.g., via replayPlot() or dev.copy())
    .Call.graphics(L_gridDirty)
    result <- .Call.graphics(fnname, ...)
  } else {
    .Call(L_gridDirty)
    result <- .Call(fnname, ...)
  }
  result
}

# A call to recordGraphics() outside of [pre|post]drawDetails methods
# will not record the expr on the grid DL.
# If a user REALLY wants to call recordGraphics(), they should use
# grid.record() instead
drawDetails.recordedGrob <- function(x, recording) {
  eval(x$expr, x$list, getNamespace("grid"))
}

grid.record <- function(expr, list,
                        name=NULL, gp=NULL, vp=NULL) {
  grid.draw(grob(expr=substitute(expr), list=list,
                 name=name, gp=gp, vp=vp, cl="recordedGrob"))
}

recordGrob <- function(expr, list,
                       name=NULL, gp=NULL, vp=NULL) {
  grob(expr=substitute(expr), list=list,
       name=name, gp=gp, vp=vp, cl="recordedGrob")
}

# Must only generate a grob, not modify drawing context
makeContent.delayedgrob <- function(x) {
    grob <- eval(x$expr, x$list, getNamespace("grid"))
    if (is.grob(grob)) {
        children <- gList(grob)
    } else if (is.gList(grob)) {
        children <- grob
    } else {
        stop("'expr' must return a grob or gList")
    }
    x <- setChildren(x, children)
    x
}

grid.delay <- function(expr, list,
                       name=NULL, gp=NULL, vp=NULL) {
    grid.draw(gTree(expr=substitute(expr), list=list,
                    name=name, gp=gp, vp=vp, cl="delayedgrob"))
}

delayGrob <- function(expr, list,
                      name=NULL, gp=NULL, vp=NULL) {
    gTree(expr=substitute(expr), list=list,
          name=name, gp=gp, vp=vp, cl="delayedgrob")
}

#  File src/library/grid/R/grob.R
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

######################################
# Grid graphical objects
#######################################

################
# CLASS DEFN
################
# A "virtual" class "gDesc" underlies both "grob" and "gPath"

initGrobAutoName <- function() {
  index <- 0
  function(prefix="GRID", suffix="GROB") {
    index <<- index + 1
    paste(prefix, suffix, index, sep=".")
  }
}

grobAutoName <- initGrobAutoName()

# Function for user to call to get "autogenerated" grob name
grobName <- function(grob=NULL, prefix="GRID") {
    if (is.null(grob))
        grobAutoName(prefix)
    else {
        if (!is.grob(grob))
            stop("invalid 'grob' argument")
        else
            grobAutoName(prefix, class(grob)[1L])
    }
}

################
# CLASS DEFN
################
# A grob has a name, a gp, and a vp
# grob inherits from gDesc
checkvpSlot <- function(vp) {
  # vp can be a viewport, a viewport name, or a viewport path
  if (!is.null(vp))
    if (!inherits(vp, "viewport") &&
        !inherits(vp, "vpPath") &&
        !is.character(vp))
      stop("invalid 'vp' slot")
  # For interactive use, allow user to specify
  # vpPath directly (i.e., w/o calling vpPath)
  if (is.character(vp))
    vp <- vpPathDirect(vp)
  vp
}

checkNameSlot <- function(x) {
  # Supply a default name if one is not given
  if (is.null(x$name))
    grobAutoName(suffix=class(x)[1L])
  else
    as.character(x$name)
}

checkgpSlot <- function(gp) {
  # gp must be a gpar
  if (!is.null(gp))
    if (!inherits(gp, "gpar"))
      stop("invalid 'gp' slot")
}

validDetails <- function(x) {
  UseMethod("validDetails")
}

validDetails.grob <- function(x) {
  x
}

validGrob <- function(x, ...) {
  UseMethod("validGrob")
}

validGrob.grob <- function(x, ...) {
  # Validate class-specific slots
  x <- validDetails(x)
  # Validate standard grob slots
  x$name <- checkNameSlot(x)
  checkgpSlot(x$gp)
  if (!is.null(x$vp))
    x$vp <- checkvpSlot(x$vp)
  return(x)
}

# This actually creates a new class derived from grob
# and returns an instance of that new class, all in one step
grob <- function(..., name=NULL, gp=NULL, vp=NULL, cl=NULL) {
  g <- list(..., name=name, gp=gp, vp=vp)
  if (!is.null(cl) &&
      !is.character(cl))
    stop("invalid 'grob' class")
  class(g) <- c(cl, "grob", "gDesc")
  validGrob(g)
}

grid.grob <- function(list.struct, cl=NULL, draw=TRUE) {
  warning("grid.grob() is deprecated; please use grob() instead")
  g <- do.call("grob", c(list.struct, cl=cl))
  if (draw)
    grid.draw(g)
  invisible(g)
}

is.grob <- function(x) {
  inherits(x, "grob")
}

as.character.grob <- function(x, ...) {
  paste0(class(x)[1L], "[", x$name, "]")
}

print.grob <- function(x, ...) {
  cat(as.character(x), "\n")
  invisible(x)
}

################
# gPath CLASS DEFN
################
# gPath is a concatenated list of names specifying a path to a grob
# Functions for creating "paths" of viewport names

gPathFromVector <- function(names) {
  n <- length(names)
  if (n < 1L)
    stop("a 'grob' path must contain at least one 'grob' name")
  if (any(bad <- !is.character(names)))
      stop(ngettext(sum(bad), "invalid grob name", "invalid grob names"),
           domain = NA)
  path <- list(path = if (n==1) NULL else
               paste(names[1L:(n-1)], collapse = .grid.pathSep),
               name = names[n], n = n)
  class(path) <- c("gPath", "path")
  path
}

gPath <- function(...) {
  names <- c(...)
  gPathFromVector(names)
}

# Create gPath from string with embedded .grid.pathSep(s)
gPathDirect <- function(path) {
  names <- unlist(strsplit(path, .grid.pathSep))
  gPathFromVector(names)
}

################
# gList CLASS DEFN
################
# Just a list of grobs
okGListelt <- function(x) {
  is.grob(x) || is.null(x) || is.gList(x)
}

is.gList <- function(x) {
    inherits(x, "gList")
}

as.gList <- function(x) {
    if (is.null(x)) {
        result <- list()
        class(result) <- "gList"
    } else if (is.grob(x)) {
        result <- list(x)
        class(result) <- "gList"
    } else if (is.gList(x)) {
        result <- x
    } else {
        stop("unable to coerce to \"gList\"")
    }
    result
}

gList <- function(...) {
    gl <- list(...)
    if (length(gl) == 0L ||
        all(sapply(gl, okGListelt, simplify=TRUE))) {
        # Ensure gList is "flat"
        # Don't want gList containing gList ...
        if (!all(sapply(gl, is.grob)))
            gl <- do.call("c", lapply(gl, as.gList))
        class(gl) <- c("gList")
        return(gl)
    } else {
        stop("only 'grobs' allowed in \"gList\"")
    }
}

addToGList <- function(x, gList) {
  UseMethod("addToGList")
}

addToGList.default <- function(x, gList) {
  if (is.null(x))
    gList
  else
    stop("invalid element to add to \"gList\"")
}

addToGList.grob <- function(x, gList) {
  if (is.null(gList))
    gList(x)
  else {
    gList[[length(gList) + 1L]] <- x
    return(gList)
  }
}

addToGList.gList <- function(x, gList) {
  gl <- c(gList, x)
  class(gl) <- "gList"
  return(gl)
}

as.character.gList <- function(x, ...) {
  paste0("(", paste(lapply(x, as.character), collapse=", "), ")")
}

print.gList <- function(x, ...) {
  cat(as.character(x), "\n")
  invisible(x)
}

`[.gList` <- function(x, index, ...) {
    cl <- class(x)
    result <- "["(unclass(x), index, ...)
    class(result) <- cl
    result
}

################
# gTree CLASS DEFN
################
# gTree extends grob
# A gTree has additional children slot
childName <- function(x) {
  x$name
}

setChildren <- function(x, children) {
  if (!inherits(x, "gTree"))
    stop("can only set 'children' for a \"gTree\"")
  if (!is.null(children) &&
      !inherits(children, "gList"))
    stop("'children' must be a \"gList\"")
  # Thin out NULL children
  if (!is.null(children)) {
    cl <- class(children)
    children <- children[!sapply(children, is.null)]
    class(children) <- cl
  }
  if (length(children)) {
    x$children <- children
    childNames <- sapply(children, childName)
    names(x$children) <- childNames
    x$childrenOrder <- childNames
  } else {
    x$children <- gList()
    x$childrenOrder <- character()
  }
  x
}

childNames <- function(gTree) {
  if (!inherits(gTree, "gTree"))
    stop("it is only valid to get 'children' from a \"gTree\"")
  gTree$childrenOrder
}

validGrob.gTree <- function(x, childrenvp, ...) {
  # Validate class-specific slots
  x <- validDetails(x)
  # Validate standard grob slots
  x$name <- checkNameSlot(x)
  checkgpSlot(x$gp)
  if (!is.null(x$vp))
    x$vp <- checkvpSlot(x$vp)
  # Only add childrenvp here so that gTree slots can
  # be validated before childrenvp get made
  # (making of childrenvp and children likely to depend
  #  on gTree slots)
  if (!is.null(childrenvp))
    x$childrenvp <- checkvpSlot(childrenvp)
  return(x)
}

gTree <- function(..., name=NULL, gp=NULL, vp=NULL,
                  children=NULL, childrenvp=NULL,
                  cl=NULL) {
  gt <- list(..., name=name, gp=gp, vp=vp)
  if (!is.null(cl) &&
      !is.character(cl))
    stop("invalid \"gTree\" class")
  class(gt) <- c(cl, "gTree", "grob", "gDesc")
  gt <- validGrob(gt, childrenvp)
  gt <- setChildren(gt, children)
  return(gt)
}

# A basic gTree that is JUST a collection of grobs
# (simply interface to gTree)
grobTree <- function(..., name=NULL, gp=NULL, vp=NULL,
                     childrenvp=NULL, cl=NULL) {
    gTree(children=gList(...),
          name=name, gp=gp, vp=vp,
          childrenvp=childrenvp, cl=cl)
}

################
# Getting just the names of the top-level grobs on the DL
################
getName <- function(elt) {
  if (inherits(elt, "grob"))
    elt$name
  else
    ""
}

getNames <- function() {
  dl <- grid.Call(L_getDisplayList)[1L:grid.Call(L_getDLindex)]
  names <- sapply(dl, getName)
  names[nzchar(names)]
}

################
# Getting/adding/removing/editing (children of [children of ...]) a gTree
################

# NOTE:  In order to cut down on repeated code, some of these
# (i.e., all but get and set) are inefficient and call get/set
# to do their work.  If speed becomes an issue, may have to
# revert to individual support for each function with highly
# repetitive code

# Get a grob from the display list
grid.get <- function(gPath, strict=FALSE, grep=FALSE, global=FALSE,
                     allDevices=FALSE) {
  if (allDevices)
    stop("'allDevices' not yet implemented")
  if (is.character(gPath))
    gPath <- gPathDirect(gPath)
  if (!inherits(gPath, "gPath"))
    stop("invalid 'gPath'")
  if (!is.logical(grep))
    stop("invalid 'grep' value")
  grep <- rep(grep, length.out=depth(gPath))
  getDLfromGPath(gPath, strict, grep, global)
}

# Just different defaults to grid.get for convenience
# Justified by usage patterns of Hadley Wickham
grid.gget <- function(..., grep=TRUE, global=TRUE) {
    grid.get(..., grep=grep, global=global)
}

# Get a child (of a child, of a child, ...) of a grob
getGrob <- function(gTree, gPath, strict=FALSE,
                    grep=FALSE, global=FALSE) {
  if (!inherits(gTree, "gTree"))
    stop("it is only valid to get a child from a \"gTree\"")
  if (is.character(gPath))
    gPath <- gPathDirect(gPath)
  if (!inherits(gPath, "gPath"))
    stop("invalid 'gPath'")
  if (depth(gPath) == 1 && strict) {
    gTree$children[[gPath$name]]
  } else {
    if (!is.logical(grep))
      stop("invalid 'grep' value")
    grep <- rep(grep, length.out=depth(gPath))
    getGTree(gTree, NULL, gPath, strict, grep, global)
  }
}

# Set a grob on the display list
# nor is it valid to specify a global destination (i.e., no global arg)
grid.set <- function(gPath, newGrob, strict=FALSE, grep=FALSE,
                     redraw=TRUE) {
  if (is.character(gPath))
    gPath <- gPathDirect(gPath)
  if (!inherits(gPath, "gPath"))
    stop("invalid 'gPath'")
  if (!is.logical(grep))
    stop("invalid 'grep' value")
  grep <- rep(grep, length.out=depth(gPath))
  result <- setDLfromGPath(gPath, newGrob, strict, grep)
  # result$index will be non-zero if matched the gPath
  if (result$index) {
    # Get the current DL index
    dl.index <- grid.Call(L_getDLindex)
    # Destructively modify the DL elt
    grid.Call(L_setDLindex, as.integer(result$index))
    grid.Call(L_setDLelt, result$grob)
    # Reset the DL index
    grid.Call(L_setDLindex, as.integer(dl.index))
    if (redraw)
      draw.all()
  } else {
    stop("'gPath' does not specify a valid child")
  }
}

# Set a grob
# nor is it valid to specify a global destination (i.e., no global arg)
setGrob <- function(gTree, gPath, newGrob, strict=FALSE, grep=FALSE) {
  if (!inherits(gTree, "gTree"))
    stop("it is only valid to set a child of a \"gTree\"")
  if (!inherits(newGrob, "grob"))
    stop("it is only valid to set a 'grob' as child of a \"gTree\"")
  if (is.character(gPath))
    gPath <- gPathDirect(gPath)
  if (!inherits(gPath, "gPath"))
    stop("invalid 'gPath'")
  if (!is.logical(grep))
    stop("invalid 'grep' value")
  grep <- rep(grep, length.out=depth(gPath))
  if (depth(gPath) == 1 && strict) {
    # gPath must specify an existing child
    if (old.pos <- nameMatch(gPath$name, gTree$childrenOrder, grep)) {
      # newGrob name must match existing name
      if (match(gTree$childrenOrder[old.pos], newGrob$name, nomatch=0L)) {
        gTree$children[[newGrob$name]] <- newGrob
      } else {
          stop(gettextf("New 'grob' name (%s) does not match 'gPath' (%s)",
                        newGrob$name, gPath), domain = NA)
      }
    } else {
        stop("'gPath' does not specify a valid child")
    }
  } else {
    gTree <- setGTree(gTree, NULL, gPath, newGrob, strict, grep)
    if (is.null(gTree))
      stop("'gPath' does not specify a valid child")
  }
  gTree
}

# Add a grob to a grob on the display list
grid.add <- function(gPath, child, strict=FALSE,
                     grep=FALSE, global=FALSE, allDevices=FALSE,
                     redraw=TRUE) {
  if (allDevices)
    stop("'allDevices' not yet implemented")
  if (is.character(gPath))
    gPath <- gPathDirect(gPath)
  if (!inherits(gPath, "gPath"))
    stop("invalid 'gPath'")
  if (!is.logical(grep))
    stop("invalid 'grep' value")
  grep <- rep(grep, length.out=depth(gPath))
  addDLfromGPath(gPath, child, strict, grep, global, redraw)
}

# Add a grob to a gTree (or a child of a (child of a ...) gTree)
addGrob <- function(gTree, child, gPath=NULL, strict=FALSE,
                    grep=FALSE, global=FALSE, warn=TRUE) {
    if (!inherits(child, "grob"))
        stop("it is only valid to add a 'grob' to a \"gTree\"")
    if (is.null(gPath)) {
        addToGTree(gTree, child)
    } else {
        if (is.character(gPath))
            gPath <- gPathDirect(gPath)
        # Only makes sense to specify a gPath for a gTree
        if (!inherits(gTree, "gTree"))
            stop("it is only valid to add a child to a \"gTree\"")
        if (!is.logical(grep))
            stop("invalid 'grep' value")
        grep <- rep(grep, length.out=depth(gPath))
        # result will be NULL if no match
        result <- addGTree(gTree, child, NULL, gPath, strict, grep, global)
        if (is.null(result)) {
            if (warn)
                warning(gettextf("'gPath' (%s) not found",
                                 as.character(gPath)),
                        domain = NA)
            gTree
        } else {
            result
        }
    }
}

# Remove a grob (or child of ...) from the display list
grid.remove <- function(gPath, warn=TRUE, strict=FALSE,
                        grep=FALSE, global=FALSE, allDevices=FALSE,
                        redraw=TRUE) {
  if (allDevices)
    stop("'allDevices' not yet implemented")
  if (is.character(gPath))
    gPath <- gPathDirect(gPath)
  if (!inherits(gPath, "gPath"))
    stop("invalid 'gPath'")
  if (!is.logical(grep))
    stop("invalid 'grep' value")
  grep <- rep(grep, length.out=depth(gPath))
  if (depth(gPath) == 1) {
    removeNameFromDL(gPath$name, strict, grep, global, warn, redraw)
  } else {
    name <- gPath$name
    gPath <- gPathDirect(gPath$path)
    greppath <- grep[-length(grep)]
    grepname <- grep[length(grep)]
    removeDLFromGPath(gPath, name, strict, greppath, grepname,
                      global, warn, redraw)
  }
}

# Just different defaults to grid.remove for convenience
# Justified by usage patterns of Hadley Wickham
grid.gremove <- function(..., grep=TRUE, global=TRUE) {
    grid.remove(..., grep=grep, global=global)
}

# Remove a child from a (child of ...) gTree
removeGrob <- function(gTree, gPath, strict=FALSE,
                       grep=FALSE, global=FALSE, warn=TRUE) {
    if (!inherits(gTree, "gTree"))
        stop("it is only valid to remove a child from a \"gTree\"")
    if (is.character(gPath))
        gPath <- gPathDirect(gPath)
    if (!inherits(gPath, "gPath"))
        stop("invalid 'gPath'")
    if (!is.logical(grep))
        stop("invalid 'grep' value")
    grep <- rep(grep, length.out=depth(gPath))
    if (depth(gPath) == 1) {
        # result will be NULL if no match
        result <- removeName(gTree, gPath$name, strict, grep, global, warn)
    } else {
        name <- gPath$name
        gPath <- gPathDirect(gPath$path)
        greppath <- grep[-length(grep)]
        grepname <- grep[length(grep)]
        # result will be NULL if no match
        result <- removeGTree(gTree, name, NULL, gPath, strict,
                              greppath, grepname, global, warn)
    }
    if (is.null(result)) {
        if (warn)
            warning(gettextf("'gPath' (%s) not found", as.character(gPath)),
                    domain = NA)
        gTree
    } else {
        result
    }
}

# Edit a grob on the display list
grid.edit <- function(gPath, ..., strict=FALSE,
                      grep=FALSE, global=FALSE, allDevices=FALSE,
                      redraw=TRUE) {
  if (allDevices)
    stop("'allDevices' not yet implemented")
  if (is.character(gPath))
    gPath <- gPathDirect(gPath)
  if (!inherits(gPath, "gPath"))
    stop("invalid 'gPath'")
  if (!is.logical(grep))
    stop("invalid 'grep' value")
  grep <- rep(grep, length.out=depth(gPath))
  specs <- list(...)
  editDLfromGPath(gPath, specs, strict, grep, global, redraw)
}

# Just different defaults to grid.edit for convenience
# Justified by usage patterns of Hadley Wickham
grid.gedit <- function(..., grep=TRUE, global=TRUE) {
    grid.edit(..., grep=grep, global=global)
}

# Edit a (child of a ...) grob
editGrob <- function(grob, gPath=NULL, ..., strict=FALSE,
                     grep=FALSE, global=FALSE, warn=TRUE) {
    specs <- list(...)
    if (is.null(gPath)) {
        editThisGrob(grob, specs)
    } else {
        if (is.character(gPath))
            gPath <- gPathDirect(gPath)
        # Only makes sense to specify a gPath for a gTree
        if (!inherits(grob, "gTree"))
            stop("it is only valid to edit a child of a \"gTree\"")
        if (!is.logical(grep))
            stop("invalid 'grep' value")
        grep <- rep(grep, length.out=depth(gPath))
        # result will be NULL if no match
        result <- editGTree(grob, specs, NULL, gPath, strict, grep, global)
        if (is.null(result)) {
            if (warn)
                warning(gettextf("'gPath' (%s) not found",
                                 as.character(gPath)),
                        domain = NA)
            grob
        } else {
            result
        }
    }
}

#########
# Generic "hook" to allow customised action on edit
#########
editDetails <- function(x, specs) {
  UseMethod("editDetails")
}

editDetails.default <- function(x, specs) {
  # Do nothing BUT return object being edited
  x
}

editDetails.gTree <- function(x, specs) {
  # Disallow editing children or childrenOrder slots directly
  if (any(specs %in% c("children", "childrenOrder")))
    stop("it is invalid to directly edit the 'children' or 'childrenOrder' slot")
  x
}

#########
# Helper functions for getting/adding/removing/editing grobs
#
# ASSUME down here that the grep argument has been replicated
# up to the length of the gPath argument
#########

# Find a "match" between a path$name and a grob$name
nameMatch <- function(pathName, grobName, grep) {
  if (grep) {
    pos <- grep(pathName, grobName)
    (length(pos) && pos == 1)
  } else {
    match(pathName, grobName, nomatch=0L)
  }
}

# Return the position of path$name in vector of names
# Return FALSE if not found
# If grep=TRUE, the answer may be a vector!
namePos <- function(pathName, names, grep) {
  if (grep) {
    pos <- grep(pathName, names)
    if (length(pos) == 0L)
      pos <- FALSE
  } else {
    pos <- match(pathName, names, nomatch=0L)
  }
  pos
}

partialPathMatch <- function(pathsofar, path, strict=FALSE, grep) {
  if (strict) {
    if (!any(grep))
      length(grep(paste0("^", pathsofar), path)) > 0L
    else {
      pathSoFarElts <- explodePath(pathsofar)
      pathElts <- explodePath(path)
      ok <- TRUE
      npsfe <- length(pathSoFarElts)
      index <- 1
      while (ok & index <= npsfe) {
        if (grep[index])
          ok <- (grep(pathSoFarElts[index], pathElts[index]) == 1)
        else
          ok <- match(pathSoFarElts[index], pathElts[index], nomatch=0L)
        index <- index + 1
      }
      ok
    }
  } else {
    # If we're not doing strict matching then anything from a full
    # path match to absolutely no match means a partial match
    # (i.e., keep looking)
    TRUE
  }
}

fullPathMatch <- function(pathsofar, gPath, strict, grep) {
  if (is.null(pathsofar))
    match <- (depth(gPath) == 1)
  else {
    path <- gPath$path
    if (!any(grep))
      if (strict)
        match <- match(pathsofar, path, nomatch=0L)
      else
        match <- (length(grep(paste0(path, "$"), pathsofar)) > 0L)
    else {
      pathSoFarElts <- explodePath(pathsofar)
      pathElts <- explodePath(path)
      npsfe <- length(pathSoFarElts)
      npe <- length(pathElts)
      if (npe > npsfe) {
        match <- FALSE
      } else {
        match <- TRUE
        index <- 1
        if (strict) {# pathSoFar same length as gPath
        } else {# pathSoFar could be longer than gPath
          pathSoFarElts <- pathSoFarElts[(npsfe - npe + 1):npsfe]
        }
        while (match && index <= npe) {
          if (grep[index])
            match <- (length(grep(pathElts[index], pathSoFarElts[index])) > 0L)
          else
            match <- match(pathSoFarElts[index], pathElts[index], nomatch = 0L)
          index <- index + 1
        }
      }
    }
  }
  match
}

#####
##### Get support
#####

# Add a grob to a result
growResult <- function(result, x) {
  UseMethod("growResult")
}

# Should only be when result is NULL
growResult.default <- function(result, x) {
  if (!is.null(result))
    stop("invalid 'result'")
  x
}

growResult.grob <- function(result, x) {
  if (is.grob(x))
    gList(result, x)
  else
    # x should be a gList
    addToGList(result, x)
}

growResult.gList <- function(result, x) {
  addToGList(x, result)
}

# A gPath may specify the child of a gTree
# (or the child of a child of a gTree, or ...)
getGrobFromGPath <- function(grob, pathsofar, gPath, strict,
                             grep, global) {
  UseMethod("getGrobFromGPath")
}

# If it's not a grob then fail
# Handles case when traversing DL
getGrobFromGPath.default <- function(grob, pathsofar, gPath, strict,
                                     grep, global) {
  NULL
}

getGrobFromGPath.grob <- function(grob, pathsofar, gPath, strict,
                                  grep, global) {
  if (depth(gPath) > 1)
    NULL
  else {
    if (nameMatch(gPath$name, grob$name, grep))
      grob
    else
      NULL
  }
}

getGTree <- function(gTree, pathsofar, gPath, strict, grep, global) {
  # Try to find pathsofar at start of gPath
  # NOTE: may be called directly with pathsofar=NULL
  if (is.null(pathsofar) ||
      (!strict && depth(gPath) == 1) ||
      partialPathMatch(pathsofar, gPath$path, strict, grep)) {
    found <- FALSE
    index <- 1
    grob <- NULL
    # Search children for match
    while (index <= length(gTree$childrenOrder) &&
           (!found || global)) {
      childName <- gTree$childrenOrder[index]
      child <- gTree$children[[childName]]
      # Special case when strict is FALSE and depth(gPath) is 1
      # Just check for gPath$name amongst children and recurse if no match
      if (!strict && depth(gPath) == 1) {
        if (nameMatch(gPath$name, childName, grep)) {
          grob <- growResult(grob, child)
          found <- TRUE
        } else {
          if (is.null(pathsofar))
            newpathsofar <- child$name
          else
            newpathsofar <- paste0(pathsofar, .grid.pathSep, childName)
          if (!is.null(newChild <- getGrobFromGPath(child, newpathsofar,
                                                    gPath, strict,
                                                    grep, global))) {
            grob <- growResult(grob, newChild)
            found <- TRUE
          }
        }
      } else {
        # Only check for match with child if have full match with pathsofar
        # If it's a complete match, look for gPath$name amongst child
        # NOTE: may be called directly with pathsofar=NULL
        if (fullPathMatch(pathsofar, gPath, strict, grep)) {
          if (nameMatch(gPath$name, childName, grep[depth(gPath)])) {
            grob <- growResult(grob, child)
            found <- TRUE
          }
        # Otherwise recurse down child
        } else {
          # NOTE: may be called directly with pathsofar=NULL
          if (is.null(pathsofar))
            newpathsofar <- child$name
          else
            newpathsofar <- paste0(pathsofar, .grid.pathSep, childName)
          if (!is.null(newChild <- getGrobFromGPath(child, newpathsofar,
                                                    gPath, strict,
                                                    grep, global))) {
            grob <- growResult(grob, newChild)
            found <- TRUE
          }
        }
      }
      index <- index + 1
    }
    if (found)
      grob
    else
      NULL
  } else {
    NULL
  }
}

getGrobFromGPath.gTree <- function(grob, pathsofar, gPath, strict,
                                   grep, global) {
  if (depth(gPath) == 1) {
    if (nameMatch(gPath$name, grob$name, grep))
      grob
    else
      if (strict)
        NULL
      else
        getGTree(grob,
                 if (is.null(pathsofar)) grob$name else pathsofar,
                 gPath, strict, grep, global)
  } else {
    getGTree(grob,
             if (is.null(pathsofar)) grob$name else pathsofar,
             gPath, strict, grep, global)
  }
}

getDLfromGPath <- function(gPath, strict, grep, global) {
  dl.index <- grid.Call(L_getDLindex)
  result <- NULL
  index <- 1
  while (index < dl.index &&
         (is.null(result) || global)) {
    grob <- getGrobFromGPath(grid.Call(L_getDLelt,
                                       as.integer(index)),
                             NULL, gPath, strict,
                             grep, global)
    if (!is.null(grob))
      result <- growResult(result, grob)
    index <- index + 1
  }
  result
}

#####
##### Set support
#####
# A gPath may specify the child of a gTree
# (or the child of a child of a gTree, or ...)
setGrobFromGPath <- function(grob, pathsofar, gPath, newGrob, strict, grep) {
  UseMethod("setGrobFromGPath")
}

# Ignore DL elements which are not grobs
setGrobFromGPath.default <- function(grob, pathsofar, gPath, newGrob,
                                     strict, grep) {
  NULL
}

setGrobFromGPath.grob <- function(grob, pathsofar, gPath, newGrob,
                                  strict, grep) {
  if (depth(gPath) > 1)
    NULL
  else {
    if (nameMatch(gPath$name, grob$name, grep))
      if (match(grob$name, newGrob$name, nomatch=0L))
        newGrob
      else
        NULL
    else
      NULL
  }
}

# Try to match gPath in gTree children
# Return NULL if cant' find match
# Return modified gTree if can find match
setGTree <- function(gTree, pathsofar, gPath, newGrob, strict, grep) {
  # Try to find pathsofar at start of gPath
  # NOTE: may be called directly with pathsofar=NULL
  if (is.null(pathsofar) ||
      (!strict && depth(gPath) == 1) ||
      partialPathMatch(pathsofar, gPath$path, strict, grep)) {
    found <- FALSE
    index <- 1
    # Search children for match
    while (index <= length(gTree$childrenOrder) && !found) {
      childName <- gTree$childrenOrder[index]
      child <- gTree$children[[childName]]
      # Special case when strict is FALSE and depth(gPath) is 1
      # Just check for gPath$name amongst children and recurse if no match
      if (!strict && depth(gPath) == 1) {
        if (nameMatch(gPath$name, childName, grep)) {
          if (match(childName, newGrob$name, nomatch=0L)) {
            gTree$children[[newGrob$name]] <- newGrob
            found <- TRUE
          } else {
            stop("the new 'grob' must have the same name as the old 'grob'")
          }
        } else {
          if (is.null(pathsofar))
            newpathsofar <- child$name
          else
            newpathsofar <- paste0(pathsofar, .grid.pathSep, childName)
          if (!is.null(newChild <- setGrobFromGPath(child, newpathsofar,
                                                    gPath, newGrob,
                                                    strict, grep))) {
            gTree$children[[childName]] <- newChild
            found <- TRUE
          }
        }
      } else {
        # Only check for match with child if have full match with pathsofar
        # If it's a complete match, look for gPath$name amongst child
        # NOTE: may be called directly with pathsofar=NULL
        if (fullPathMatch(pathsofar, gPath, strict, grep)) {
          if (nameMatch(gPath$name, childName, grep[depth(gPath)])) {
            if (match(childName, newGrob$name, nomatch=0L)) {
                gTree$children[[newGrob$name]] <- newGrob
                found <- TRUE
            } else {
                stop("the new 'grob' must have the same name as the old 'grob'")
            }
          }
        # Otherwise recurse down child
        } else {
          # NOTE: may be called directly with pathsofar=NULL
          if (is.null(pathsofar))
            newpathsofar <- child$name
          else
            newpathsofar <- paste0(pathsofar, .grid.pathSep, childName)
          if (!is.null(newChild <- setGrobFromGPath(child, newpathsofar,
                                                    gPath, newGrob,
                                                    strict, grep))) {
            gTree$children[[childName]] <- newChild
            found <- TRUE
          }
        }
      }
      index <- index + 1
    }
    if (found)
      gTree
    else
      NULL
  } else {
    NULL
  }
}

setGrobFromGPath.gTree <- function(grob, pathsofar, gPath, newGrob,
                                   strict, grep) {
  if (depth(gPath) == 1) {
    if (nameMatch(gPath$name, grob$name, grep))
      if (match(grob$name, newGrob$name, nomatch=0L))
        newGrob
      else
        stop("the new 'grob' must have the same name as the old 'grob'")
    else
      if (strict)
        NULL
      else
        setGTree(grob,
                 if (is.null(pathsofar)) grob$name else pathsofar,
                 gPath, newGrob, strict, grep)
  } else {
    setGTree(grob,
             # Initialise pathsofar if first time through
             if (is.null(pathsofar)) grob$name else pathsofar,
             gPath, newGrob, strict, grep)
  }
}

setDLfromGPath <- function(gPath, newGrob, strict, grep) {
  dl.index <- grid.Call(L_getDLindex)
  index <- 1
  result <- list(index=0, grob=NULL)
  while (index < dl.index &&
         result$index == 0) {
    result$grob <- setGrobFromGPath(grid.Call(L_getDLelt,
                                              as.integer(index)),
                                    NULL, gPath, newGrob, strict, grep)
    if (!is.null(result$grob))
      result$index <- index
    index <- index + 1
  }
  result
}

#####
##### Edit support
#####
editThisGrob <- function(grob, specs) {
  for (i in names(specs))
    if (nzchar(i))
      # Handle gp as special case
      if (match(i, "gp", nomatch=0))
        # Handle NULL as special case
        if (is.null(specs[[i]]))
          grob[i] <- list(gp=NULL)
        else
          grob$gp <- mod.gpar(grob$gp, specs$gp)
      # If there is no slot with the argument name, just ignore that argument
      else if (match(i, names(grob), nomatch=0))
        # Handle NULL as special case
        if (is.null(specs[[i]]))
          grob[i] <- eval(substitute(list(i=NULL)))
        else
          grob[[i]] <- specs[[i]]
      else
        warning(gettextf("slot '%s' not found", i), domain = NA)
  # Check grob slots are ok before trying to do anything with them
  # in editDetails
  # grob$childrenvp may be non-NULL for a gTree
  grob <- validGrob(grob, grob$childrenvp)
  editDetails(grob, specs)
}

# A gPath may specify the child of a gTree
# (or the child of a child of a gTree, or ...)
editGrobFromGPath <- function(grob, specs, pathsofar, gPath, strict,
                              grep, global) {
  UseMethod("editGrobFromGPath")
}

# If it's not a grob then fail
# Handles case when traversing DL
editGrobFromGPath.default <- function(grob, specs,
                                      pathsofar, gPath, strict,
                                      grep, global) {
  NULL
}

editGrobFromGPath.grob <- function(grob, specs,
                                   pathsofar, gPath, strict,
                                   grep, global) {
  if (depth(gPath) > 1)
    NULL
  else {
    if (nameMatch(gPath$name, grob$name, grep))
      editThisGrob(grob, specs)
    else
      NULL
  }
}

editGTree <- function(gTree, specs, pathsofar, gPath, strict,
                      grep, global) {
  # Try to find pathsofar at start of gPath
  # NOTE: may be called directly with pathsofar=NULL
  if (is.null(pathsofar) ||
      (!strict && depth(gPath) == 1) ||
      partialPathMatch(pathsofar, gPath$path, strict, grep)) {
    found <- FALSE
    index <- 1
    # Search children for match
    while (index <= length(gTree$childrenOrder) &&
           (!found || global)) {
      childName <- gTree$childrenOrder[index]
      child <- gTree$children[[childName]]
      # Special case when strict is FALSE and depth(gPath) is 1
      # Just check for gPath$name amongst children and recurse if no match
      if (!strict && depth(gPath) == 1) {
        if (nameMatch(gPath$name, childName, grep)) {
          gTree$children[[childName]] <- editThisGrob(child, specs)
          found <- TRUE
        } else {
          if (is.null(pathsofar))
            newpathsofar <- child$name
          else
            newpathsofar <- paste0(pathsofar, .grid.pathSep, childName)
          if (!is.null(newChild <- editGrobFromGPath(child, specs,
                                                     newpathsofar,
                                                     gPath, strict,
                                                     grep, global))) {
            gTree$children[[childName]] <- newChild
            found <- TRUE
          }
        }
      } else {
        # Only check for match with child if have full match with pathsofar
        # If it's a complete match, look for gPath$name amongst child
        # NOTE: may be called directly with pathsofar=NULL
        if (fullPathMatch(pathsofar, gPath, strict, grep)) {
          if (nameMatch(gPath$name, childName, grep[depth(gPath)])) {
            gTree$children[[childName]] <- editThisGrob(child, specs)
            found <- TRUE
          }
        # Otherwise recurse down child
        } else {
          # NOTE: may be called directly with pathsofar=NULL
          if (is.null(pathsofar))
            newpathsofar <- child$name
          else
            newpathsofar <- paste0(pathsofar, .grid.pathSep, childName)
          if (!is.null(newChild <- editGrobFromGPath(child, specs,
                                                     newpathsofar,
                                                     gPath, strict,
                                                     grep, global))) {
            gTree$children[[childName]] <- newChild
            found <- TRUE
          }
        }
      }
      index <- index + 1
    }
    if (found)
      gTree
    else
      NULL
  } else {
    NULL
  }
}

editGrobFromGPath.gTree <- function(grob, specs,
                                    pathsofar, gPath, strict,
                                    grep, global) {
  if (depth(gPath) == 1) {
    if (nameMatch(gPath$name, grob$name, grep))
      editThisGrob(grob, specs)
    else
      if (strict)
        NULL
      else
        editGTree(grob, specs,
                  if (is.null(pathsofar)) grob$name else pathsofar,
                  gPath, strict, grep, global)
  } else {
    editGTree(grob, specs,
              if (is.null(pathsofar)) grob$name else pathsofar,
              gPath, strict, grep, global)
  }
}

editDLfromGPath <- function(gPath, specs, strict, grep, global, redraw) {
  dl.index <- grid.Call(L_getDLindex)
  index <- 1
  grob <- NULL
  found <- FALSE
  while (index < dl.index &&
         (is.null(grob) || global)) {
    grob <- editGrobFromGPath(grid.Call(L_getDLelt,
                                        as.integer(index)),
                              specs,
                              NULL, gPath, strict, grep, global)
    if (!is.null(grob)) {
      # Destructively modify the DL elt
      grid.Call(L_setDLindex, as.integer(index))
      grid.Call(L_setDLelt, grob)
      # Reset the DL index
      grid.Call(L_setDLindex, as.integer(dl.index))
      found <- TRUE
    }
    index <- index + 1
  }
  if (!found)
    stop(gettextf("'gPath' (%s) not found", as.character(gPath)), domain = NA)
  else if (redraw)
    draw.all()
}

#####
##### Add support
#####

# Assume that child is a grob
addToGTree <- function(gTree, child) {
  if (!inherits(gTree, "gTree"))
    stop("it is only valid to add a child to a \"gTree\"")
  gTree$children[[child$name]] <- child
  # Handle case where child name already exists (so will be overwritten)
  if (old.pos <- match(child$name, gTree$childrenOrder, nomatch=0))
    gTree$childrenOrder <- gTree$childrenOrder[-old.pos]
  gTree$childrenOrder <- c(gTree$childrenOrder, child$name)
  gTree
}

# A gPath may specify the child of a gTree
# (or the child of a child of a gTree, or ...)
addGrobFromGPath <- function(grob, child, pathsofar, gPath, strict,
                             grep, global) {
  UseMethod("addGrobFromGPath")
}

# If it's not a grob then fail
# Handles case when traversing DL
addGrobFromGPath.default <- function(grob, child,
                                     pathsofar, gPath, strict,
                                     grep, global) {
  NULL
}

# If no match then fail
# If match then error!
addGrobFromGPath.grob <- function(grob, child,
                                  pathsofar, gPath, strict,
                                  grep, global) {
  if (depth(gPath) > 1)
    NULL
  else {
    if (nameMatch(gPath$name, grob$name, grep))
      stop("it is only valid to add a child to a \"gTree\"")
    else
      NULL
  }
}

# In this function, the grob being added is called "grob"
# (in all others it is called "child"
addGTree <- function(gTree, grob, pathsofar, gPath, strict,
                     grep, global) {
  # Try to find pathsofar at start of gPath
  # NOTE: may be called directly with pathsofar=NULL
  if (is.null(pathsofar) ||
      (!strict && depth(gPath) == 1) ||
      partialPathMatch(pathsofar, gPath$path, strict, grep)) {
    found <- FALSE
    index <- 1
    # Search children for match
    while (index <= length(gTree$childrenOrder) &&
           (!found || global)) {
      childName <- gTree$childrenOrder[index]
      child <- gTree$children[[childName]]
      # Special case when strict is FALSE and depth(gPath) is 1
      # Just check for gPath$name amongst children and recurse if no match
      if (!strict && depth(gPath) == 1) {
        if (nameMatch(gPath$name, childName, grep)) {
          gTree$children[[childName]] <- addToGTree(child, grob)
          found <- TRUE
        } else {
          if (is.null(pathsofar))
            newpathsofar <- child$name
          else
            newpathsofar <- paste0(pathsofar, .grid.pathSep, childName)
          if (!is.null(newChild <- addGrobFromGPath(child, grob,
                                                    newpathsofar,
                                                    gPath, strict,
                                                    grep, global))) {
            gTree$children[[childName]] <- newChild
            found <- TRUE
          }
        }
      } else {
        # Only check for match with child if have full match with pathsofar
        # If it's a complete match, look for gPath$name amongst child
        # NOTE: may be called directly with pathsofar=NULL
        if (fullPathMatch(pathsofar, gPath, strict, grep)) {
          if (nameMatch(gPath$name, childName, grep[depth(gPath)])) {
            gTree$children[[childName]] <- addToGTree(child, grob)
            found <- TRUE
          }
        # Otherwise recurse down child
        } else {
          # NOTE: may be called directly with pathsofar=NULL
          if (is.null(pathsofar))
            newpathsofar <- child$name
          else
            newpathsofar <- paste0(pathsofar, .grid.pathSep, childName)
          if (!is.null(newChild <- addGrobFromGPath(child, grob,
                                                    newpathsofar,
                                                    gPath, strict,
                                                    grep, global))) {
            gTree$children[[childName]] <- newChild
            found <- TRUE
          }
        }
      }
      index <- index + 1
    }
    if (found)
      gTree
    else
      NULL
  } else {
    NULL
  }
}

addGrobFromGPath.gTree <- function(grob, child,
                                   pathsofar, gPath, strict,
                                   grep, global) {
  if (depth(gPath) == 1) {
    if (nameMatch(gPath$name, grob$name, grep))
      addToGTree(grob, child)
    else
      if (strict)
        NULL
      else
        addGTree(grob, child,
                 if (is.null(pathsofar)) grob$name else pathsofar,
                 gPath, strict, grep, global)
  } else {
    addGTree(grob, child,
             if (is.null(pathsofar)) grob$name else pathsofar,
             gPath, strict, grep, global)
  }
}

addDLfromGPath <- function(gPath, child, strict, grep, global, redraw) {
  dl.index <- grid.Call(L_getDLindex)
  index <- 1
  grob <- NULL
  found <- FALSE
  while (index < dl.index &&
         (is.null(grob) || global)) {
    grob <- addGrobFromGPath(grid.Call(L_getDLelt,
                                       as.integer(index)),
                             child,
                             NULL, gPath, strict, grep, global)
    if (!is.null(grob)) {
      # Destructively modify the DL elt
      grid.Call(L_setDLindex, as.integer(index))
      grid.Call(L_setDLelt, grob)
      # Reset the DL index
      grid.Call(L_setDLindex, as.integer(dl.index))
      found <- TRUE
    }
    index <- index + 1
  }
  if (!found)
    stop(gettextf("'gPath' (%s) not found", gPath), domain = NA)
  else if (redraw)
    draw.all()
}

#####
##### Remove support
#####

removeFromGTree <- function(gTree, name, grep) {
  if (!inherits(gTree, "gTree"))
    stop("it is only valid to remove a child from a \"gTree\"")
  if (grep) {
    old.pos <- grep(name, gTree$childrenOrder)
    if (length(old.pos) == 0L)
      old.pos <- 0
  } else {
    old.pos <- match(name, gTree$childrenOrder, nomatch=0)
  }
  if (old.pos > 0) {
    # name might be a regexp so use real name
    gTree$children[[gTree$childrenOrder[old.pos]]] <- NULL
    gTree$childrenOrder <- gTree$childrenOrder[-old.pos]
    gTree
  } else {
    NULL
  }
}

# A gPath may specify the child of a gTree
# (or the child of a child of a gTree, or ...)
removeGrobFromGPath <- function(grob, name, pathsofar, gPath, strict,
                                grep, grepname, global, warn) {
  UseMethod("removeGrobFromGPath")
}

# If it's not a grob then fail
# Handles case when traversing DL
removeGrobFromGPath.default <- function(grob, name,
                                        pathsofar, gPath, strict,
                                        grep, grepname, global, warn) {
  NULL
}

# ALWAYS fail
# (either no match or match but grob has no children!)
removeGrobFromGPath.grob <- function(grob, name,
                                     pathsofar, gPath, strict,
                                     grep, grepname, global, warn) {
  NULL
}

removeGTree <- function(gTree, name, pathsofar, gPath, strict,
                        grep, grepname, global, warn) {
  # Try to find pathsofar at start of gPath
  # NOTE: may be called directly with pathsofar=NULL
  if (is.null(pathsofar) ||
      (!strict && depth(gPath) == 1) ||
      partialPathMatch(pathsofar, gPath$path, strict, grep)) {
    found <- FALSE
    index <- 1
    # Search children for match
    while (index <= length(gTree$childrenOrder) &&
           (!found || global)) {
      childName <- gTree$childrenOrder[index]
      child <- gTree$children[[childName]]
      # Special case when strict is FALSE and depth(gPath) is 1
      # Just check for gPath$name amongst children and recurse if no match
      if (!strict && depth(gPath) == 1) {
        # NOTE: child has to be a gTree if we hope to find a child in it!
        if (inherits(child, "gTree") &&
            nameMatch(gPath$name, childName, grep)) {
          newchild <- removeFromGTree(child, name, grepname)
          if (!is.null(newchild)) {
            gTree$children[[childName]] <- newchild
            found <- TRUE
          }
        } else {
          if (is.null(pathsofar))
            newpathsofar <- child$name
          else
            newpathsofar <- paste0(pathsofar, .grid.pathSep, childName)
          if (!is.null(newChild <- removeGrobFromGPath(child, name,
                                                       newpathsofar,
                                                       gPath, strict,
                                                       grep, grepname,
                                                       global, warn))) {
            gTree$children[[childName]] <- newChild
            found <- TRUE
          }
        }
      } else {
        # Only check for match with child if have full match with pathsofar
        # If it's a complete match, look for gPath$name amongst child
        # NOTE: may be called directly with pathsofar=NULL
        if (fullPathMatch(pathsofar, gPath, strict, grep)) {
          # NOTE: child has to be a gTree if we hope to find a child in it!
          if (inherits(child, "gTree") &&
              nameMatch(gPath$name, childName, grep[depth(gPath)])) {
            newchild <- removeFromGTree(child, name, grepname)
            if (!is.null(newchild)) {
              gTree$children[[childName]] <- newchild
              found <- TRUE
            }
          }
        # Otherwise recurse down child
        } else {
          # NOTE: may be called directly with pathsofar=NULL
          if (is.null(pathsofar))
            newpathsofar <- child$name
          else
            newpathsofar <- paste0(pathsofar, .grid.pathSep, childName)
          if (!is.null(newChild <- removeGrobFromGPath(child, name,
                                                       newpathsofar,
                                                       gPath, strict,
                                                       grep, grepname,
                                                       global, warn))) {
            gTree$children[[childName]] <- newChild
            found <- TRUE
          }
        }
      }
      index <- index + 1
    }
    if (found)
      gTree
    else
      NULL
  } else {
    NULL
  }
}

removeGrobFromGPath.gTree <- function(grob, name,
                                      pathsofar, gPath, strict,
                                      grep, grepname, global, warn) {
  if (depth(gPath) == 1) {
    if (nameMatch(gPath$name, grob$name, grep))
      removeFromGTree(grob, name, grepname)
    else
      if (strict)
        NULL
      else
        removeGTree(grob, name,
                    if (is.null(pathsofar)) grob$name else pathsofar,
                    gPath, strict, grep, grepname, global, warn)
  } else {
    removeGTree(grob, name,
                if (is.null(pathsofar)) grob$name else pathsofar,
                gPath, strict, grep, grepname, global, warn)
  }
}

removeDLFromGPath <- function(gPath, name, strict, grep, grepname, global,
                              warn, redraw) {
  dl.index <- grid.Call(L_getDLindex)
  index <- 1
  grob <- NULL
  found <- FALSE
  while (index < dl.index &&
         (is.null(grob) || global)) {
    grob <- removeGrobFromGPath(grid.Call(L_getDLelt, as.integer(index)),
                                name,
                                NULL, gPath, strict, grep, grepname,
                                global, warn)
    if (!is.null(grob)) {
      # Destructively modify the DL elt
      grid.Call(L_setDLindex, as.integer(index))
      grid.Call(L_setDLelt, grob)
      # Reset the DL index
      grid.Call(L_setDLindex, as.integer(dl.index))
      found <- TRUE
    }
    index <- index + 1
  }
  if (!found)
    stop(gettextf("gPath (%s) not found",
                  paste(gPath, name, sep=.grid.pathSep)),
                  domain = NA)
  else if (redraw)
    draw.all()
}

#####
##### Remove NAME support
#####

# NEVER called when strict=TRUE
removeGrobFromName <- function(grob, name, grep, global, warn) {
  UseMethod("removeGrobFromName")
}

removeGrobFromName.grob <- function(grob, name, grep, global, warn) {
  NULL
}

# For a gTree, just recurse straight back to removeName
removeGrobFromName.gTree <- function(grob, name, grep, global, warn) {
    removeName(grob, name, FALSE, grep, global, warn)
}

removeName <- function(gTree, name, strict, grep, global, warn) {
  found <- FALSE
  index <- 1
  # Search children for match
  while (index <= length(gTree$childrenOrder) &&
         (!found || global)) {
    childName <- gTree$childrenOrder[index]
    child <- gTree$children[[childName]]
    # Just check child name and recurse if no match
    if (nameMatch(name, childName, grep)) {
      # name might be a regexp, so get real name
      gTree$children[[gTree$childrenOrder[index]]] <- NULL
      gTree$childrenOrder <- gTree$childrenOrder[-index]
      found <- TRUE
      # If deleted the child, do NOT increase index!
    } else if (strict) {
      NULL
      index <- index + 1
    } else {
      if (!is.null(newChild <- removeGrobFromName(child, name,
                                                  grep, global, warn))) {
        gTree$children[[childName]] <- newChild
        found <- TRUE
      }
      index <- index + 1
    }
  }
  if (found)
    gTree
  else
    NULL
}

removeNameFromDL <- function(name, strict, grep, global, warn, redraw) {
  dl.index <- grid.Call(L_getDLindex)
  index <- 1
  grob <- NULL
  found <- FALSE
  while (index < dl.index &&
         (is.null(grob) || global)) {
    grob <- grid.Call(L_getDLelt, as.integer(index))
    if (inherits(grob, "grob")) {
      # If match top-level grob, remove it from DL
      if (nameMatch(name, grob$name, grep)) {
        # Destructively modify the DL elt
        grid.Call(L_setDLindex, as.integer(index))
        grid.Call(L_setDLelt, NULL)
        # Reset the DL index
        grid.Call(L_setDLindex, as.integer(dl.index))
        found <- TRUE
      # Otherwise search down it for match
      } else {
        if (!strict) {
          grob <- removeGrobFromName(grob, name, grep, global, warn)
          if (!is.null(grob)) {
            # Destructively modify the DL elt
            grid.Call(L_setDLindex, as.integer(index))
            grid.Call(L_setDLelt, grob)
            # Reset the DL index
            grid.Call(L_setDLindex, as.integer(dl.index))
            found <- TRUE
          }
        }
      }
    } else {
      grob <- NULL
    }
    index <- index + 1
  }
  if (!found) {
    if (warn)
        stop(gettextf("gPath (%s) not found", name), domain = NA)
  } else if (redraw)
    draw.all()
}

################
# Finding a grob from a grob name
################
findgrob <- function(x, name) {
  UseMethod("findgrob")
}

findgrob.default <- function(x, name) {
  NULL
}

findgrob.grob <- function(x, name) {
  if (match(name, x$name, nomatch=0L))
    x
  else
    NULL
}

findGrobinDL <- function(name) {
  dl.index <- grid.Call(L_getDLindex)
  result <- NULL
  index <- 1
  while (index < dl.index && is.null(result)) {
    result <- findgrob(grid.Call(L_getDLelt, as.integer(index)), name)
    index <- index + 1
  }
  if (is.null(result))
    stop(gettextf("grob '%s' not found", name), domain = NA)
  result
}

findGrobinChildren <- function(name, children) {
  nc <- length(children)
  result <- NULL
  index <- 1
  while (index <= nc && is.null(result)) {
    result <- findgrob(children[[index]], name)
    index <- index + 1
  }
  if (is.null(result))
    stop(gettextf("grob '%s' not found", name), domain = NA)
  result
}

################
# grid.draw
################
# Use generic function "draw" rather than generic function "print"
# because want graphics functions to produce graphics output
# without having to be evaluated at the command-line AND without having
# to necessarily produce a single graphical object as the return value
# (i.e., so that simple procedural code can be written just for its
# side-effects).
# For example, so that the following code will draw
# a rectangle AND a line:
#   temp <- function() { grid.lines(); grid.rect() }
#   temp()
grid.draw <- function(x, recording=TRUE) {
  UseMethod("grid.draw")
}

grid.draw.default <- function(x, recording) {
  # Allow for "holes" in the DL if a grob has been removed
  if (!is.null(x))
    stop("invalid element in the display list")
}

grid.draw.viewport <- function(x, recording) {
  pushViewport(x, recording=FALSE)
}

grid.draw.vpPath <- function(x, recording) {
  # Assumes strict=FALSE, BUT in order to get onto
  # display list it must have worked => strict same as non-strict
  downViewport(x, recording=FALSE)
}

grid.draw.pop <- function(x, recording) {
  popViewport(x, recording=FALSE)
}

grid.draw.up <- function(x, recording) {
  upViewport(x, recording=FALSE)
}

pushgrobvp <- function(vp) {
  UseMethod("pushgrobvp")
}

pushgrobvp.viewport <- function(vp) {
  pushViewport(vp, recording=FALSE)
}

pushgrobvp.vpPath <- function(vp) {
  downViewport(vp, strict=TRUE, recording=FALSE)
}

popgrobvp <- function(vp) {
  UseMethod("popgrobvp")
}

popgrobvp.viewport <- function(vp) {
  # NOTE that the grob's vp may be a vpStack/List/Tree
  upViewport(depth(vp), recording=FALSE)
}

popgrobvp.vpPath <- function(vp) {
  upViewport(depth(vp), recording=FALSE)
}

preDraw <- function(x) {
  UseMethod("preDraw")
}

pushvpgp <- function(x) {
  if (!is.null(x$vp))
    pushgrobvp(x$vp)
  if (!is.null(x$gp)) {
    set.gpar(x$gp)
  }
}

makeContext <- function(x) {
    UseMethod("makeContext")
}

makeContext.default <- function(x) {
    x
}

makeContent <- function(x) {
    UseMethod("makeContent")
}

makeContent.default <- function(x) {
    x
}

preDraw.grob <- function(x) {
    # Allow customisation of x$vp
    x <- makeContext(x)
    # automatically push/pop the viewport and set/unset the gpar
    pushvpgp(x)
    preDrawDetails(x)
    x
}

preDraw.gTree <- function(x) {
    # Allow customisation of x$vp (and x$childrenvp)
    x <- makeContext(x)
    # Make this gTree the "current grob" for evaluation of
    # grobwidth/height units via gPath
    # Do this as a .Call.graphics to get it onto the base display list
    grid.Call.graphics(L_setCurrentGrob, x)
    # automatically push/pop the viewport
    pushvpgp(x)
    # Push then "up" childrenvp
    if (!is.null(x$childrenvp)) {
        # Save any x$gp gpar settings
        tempgp <- grid.Call(L_getGPar)
        pushViewport(x$childrenvp, recording=FALSE)
        upViewport(depth(x$childrenvp), recording=FALSE)
        # reset the x$gp gpar settings
        # The upViewport above may have overwritten them with
        # the previous vp$gp settings
        grid.Call.graphics(L_setGPar, tempgp)
    }
    preDrawDetails(x)
    x
}

postDraw <- function(x) {
    UseMethod("postDraw")
}

postDraw.grob <- function(x) {
    postDrawDetails(x)
    if (!is.null(x$vp))
        popgrobvp(x$vp)
}

drawGrob <- function(x) {
    # Temporarily turn off the grid DL so that
    # nested calls to drawing code do not get recorded
    dlon <- grid.Call(L_setDLon, FALSE)
    # If get error or user-interrupt, need to reset state
    # Need to turn grid DL back on (if it was on)
    on.exit(grid.Call(L_setDLon, dlon))
    # Save current gpar
    tempgpar <- grid.Call(L_getGPar)
    # If get error or user-interrupt, need to reset state
    # Need to restore current grob (gtree predraw sets current grob)
    # Need to restore gpar settings (set by gtree itself and/or its vp)
    # This does not need to be a grid.Call.graphics() because
    # we are nested within a recordGraphics()
    # Do not call set.gpar because set.gpar accumulates cex
    on.exit(grid.Call(L_setGPar, tempgpar), add=TRUE)
    # Setting up the drawing context may involve modifying the grob
    # (typically only x$vp) but the modified grob is needed for postDraw()
    x <- preDraw(x)
    # Allow customisation of x
    # (should only return a basic grob that has a drawDetails()
    #  method, otherwise nothing will be drawn)
    x <- makeContent(x)
    # Do any class-specific drawing
    drawDetails(x, recording=FALSE)
    postDraw(x)
}

grid.draw.grob <- function(x, recording=TRUE) {
    engineDLon <- grid.Call(L_getEngineDLon)
    if (engineDLon)
        recordGraphics(drawGrob(x),
                       list(x=x),
                       getNamespace("grid"))
    else
        drawGrob(x)
    if (recording)
        record(x)
    invisible()
}

drawGList <- function(x) {
    # DO NOT turn off grid DL.
    # A top-level gList does not itself go on the DL,
    # but its children do.
    # A gList which is part of some other grob (e.g., children
    # of a gTree) will be "protected" by the gTree
    # turning off the DL.
    lapply(x, grid.draw)
}

grid.draw.gList <- function(x, recording=TRUE) {
    engineDLon <- grid.Call(L_getEngineDLon)
    if (engineDLon)
        recordGraphics(drawGList(x),
                       list(x=x),
                       getNamespace("grid"))
    else
        drawGList(x)
    invisible()
}

drawGTree <- function(x) {
    # Temporarily turn off the grid DL so that
    # nested calls to drawing code do not get recorded
    dlon <- grid.Call(L_setDLon, FALSE)
    # If get error or user-interrupt, need to reset state
    # Need to turn grid DL back on (if it was on)
    on.exit(grid.Call(L_setDLon, dlon))
    # Save current grob and current gpar
    tempgrob <- grid.Call(L_getCurrentGrob)
    tempgpar <- grid.Call(L_getGPar)
    # If get error or user-interrupt, need to reset state
    # Need to restore current grob (gtree predraw sets current grob)
    # Need to restore gpar settings (set by gtree itself and/or its vp)
    # This does not need to be a grid.Call.graphics() because
    # we are nested within a recordGraphics()
    # Do not call set.gpar because set.gpar accumulates cex
    on.exit({ grid.Call(L_setGPar, tempgpar)
              grid.Call(L_setCurrentGrob, tempgrob)
            }, add=TRUE)
    # Setting up the drawing context may involve modifying the grob
    # (typically only x$vp) but the modified grob is needed for postDraw()
    x <- preDraw(x)
    # Allow customisation of x (should be confined to x$children)
    x <- makeContent(x)
    # Do any class-specific drawing
    drawDetails(x, recording=FALSE)
    # Draw all children IN THE RIGHT ORDER
    for (i in x$childrenOrder)
      grid.draw(x$children[[i]], recording=FALSE)
    postDraw(x)
}

grid.draw.gTree <- function(x, recording=TRUE) {
    engineDLon <- grid.Call(L_getEngineDLon)
    if (engineDLon)
        recordGraphics(drawGTree(x),
                       list(x=x),
                       getNamespace("grid"))
    else
        drawGTree(x)
    if (recording)
        record(x)
    invisible()
}

draw.all <- function() {
    grid.newpage(recording=FALSE)
    dl.index <- grid.Call(L_getDLindex)
    if (dl.index > 1)
        # Start at 2 because first element is viewport[ROOT]
        for (i in 2:dl.index) {
            grid.draw(grid.Call(L_getDLelt, as.integer(i - 1)),
                      recording=FALSE)
        }
}

draw.details <- function(x, recording) {
    .Deprecated("drawDetails")
    UseMethod("drawDetails")
}

preDrawDetails <- function(x) {
    UseMethod("preDrawDetails")
}

preDrawDetails.grob <- function(x) {
}

postDrawDetails <- function(x) {
    UseMethod("postDrawDetails")
}

postDrawDetails.grob <- function(x) {
}

drawDetails <- function(x, recording) {
    UseMethod("drawDetails")
}

drawDetails.grob <- function(x, recording) {
}

grid.copy <- function(grob) {
    warning("this function is redundant and will disappear in future versions",
            domain = NA)
    grob
}

################################
# Flattening a grob

force <- function(x) {
    UseMethod("force")
}

# The default action is to leave 'x' untouched
# BUT it is also necessary to enforce the drawing context
# for viewports and vpPaths
force.default <- function(x) {
    grid.draw(x, recording=FALSE)
    x
}

# This allows 'x' to be modified, but may not
# change 'x' at all
force.grob <- function(x) {
    # Copy of the original object to allow a "revert"
    originalX <- x
    # Same set up as drawGrob()
    dlon <- grid.Call(L_setDLon, FALSE)
    on.exit(grid.Call(L_setDLon, dlon))
    tempgpar <- grid.Call(L_getGPar)
    on.exit(grid.Call(L_setGPar, tempgpar), add=TRUE)
    # Same drawing context set up as drawGrob()
    # including enforcing the drawing context
    x <- preDraw(x)
    # Same drawing content set up as drawGrob() ...
    x <- makeContent(x)
    # BUT NO DRAWING
    # Same context clean up as drawGrob()
    postDraw(x)
    # If 'x' has not changed, just return original 'x'
    # Also, do not bother with saving original
    # If 'x' has changed ...
    if (!identical(x, originalX)) {
        # Store the original object to allow a "revert"
        x$.ORIGINAL <- originalX
        # Return the 'x' that would have been drawn
        # This will typically be a standard R primitive
        # (which do not have makeContext() or makeContent()
        #  methods, only drawDetails())
        # BUT ot be safe add "forcedgrob" class so that subsequent
        # draws will NOT run makeContext() or makeContent()
        # methods
        class(x) <- c("forcedgrob", class(x))
    }
    x
}

# This allows 'x' to be modified, but may not
# change 'x' at all
force.gTree <- function(x) {
    # Copy of the original object to allow a "revert"
    originalX <- x
    # Same set up as drawGTree()
    dlon <- grid.Call(L_setDLon, FALSE)
    on.exit(grid.Call(L_setDLon, dlon))
    tempgrob <- grid.Call(L_getCurrentGrob)
    tempgpar <- grid.Call(L_getGPar)
    on.exit({ grid.Call(L_setGPar, tempgpar)
              grid.Call(L_setCurrentGrob, tempgrob)
            }, add=TRUE)
    # Same drawing context set up as drawGTree(),
    # including enforcing the drawing context
    x <- preDraw(x)
    # Same drawing content set up as drawGTree() ...
    x <- makeContent(x)
    # Ensure that children are also forced
    x$children <- do.call("gList", lapply(x$children, force))
    # BUT NO DRAWING
    # Same context clean up as drawGTree()
    postDraw(x)
    # If 'x' has changed ...
    if (!identical(x, originalX)) {
        # Store the original object to allow a "revert"
        x$.ORIGINAL <- originalX
        # Return the 'x' that would have been drawn
        # This will typically be a vanilla gTree with children to draw
        # (which will not have makeContext() or makeContent() methods)
        # BUT to be safe add "forcedgrob" class so that subsequent
        # draws will NOT run makeContext() or makeContent()
        # methods
        class(x) <- c("forcedgrob", class(x))
    }
    x
}

# A "forcedgrob" does NOT modify context or content at
# drawing time
makeContext.forcedgrob <- function(x) x

makeContent.forcedgrob <- function(x) x

# NOTE that things will get much trickier if allow
# grid.force(gPath = ...)
grid.force <- function(redraw=TRUE) {
    dl.index <- grid.Call(L_getDLindex)
    if (dl.index > 1) {
        # Start at 2 because first element is viewport[ROOT]
        for (i in 2:dl.index) {
            grid.Call(L_setDLindex, as.integer(i - 1))
            grid.Call(L_setDLelt,
                      force(grid.Call(L_getDLelt, as.integer(i - 1))))
        }
        grid.Call(L_setDLindex, dl.index)
    }
    if (redraw) {
        draw.all()
    }
}

revert <- function(x) {
    UseMethod("revert")
}

revert.default <- function(x) {
    x
}

# Only need to revert "forcedgrob"s
revert.forcedgrob <- function(x) {
    x$.ORIGINAL
}

# No need for recursion for gTree because if top-level grob
# changed its children then top-level grob will have retained
# revert version of its entire self (including children)

# NOTE that things will get much trickier if allow
# grid.revert(gPath = ...)
grid.revert <- function(redraw=TRUE) {
    dl.index <- grid.Call(L_getDLindex)
    if (dl.index > 1) {
        # Start at 2 because first element is viewport[ROOT]
        for (i in 2:dl.index) {
            grid.Call(L_setDLindex, as.integer(i - 1))
            grid.Call(L_setDLelt,
                      revert(grid.Call(L_getDLelt, as.integer(i - 1))))
        }
        grid.Call(L_setDLindex, dl.index)
    }
    if (redraw) {
        draw.all()
    }
}

###############################
# Reordering grobs

# Reorder the children of a gTree
# Order may be specified as a character vector
#   Character vector MUST name existing children
# Order may be specified as a numeric vector
#   (which makes it easy to say something like
#    "make last child the first child")
#   Numeric vector MUST be within range 1:numChildren
# Only unique order values used
# Any children NOT specified by order are appended to
#   front or back of order (depending on 'front' argument)
# Order is ALWAYS back-to-front
reorderGrob <- function(x, order, back=TRUE) {
    if (!inherits(x, "gTree"))
        stop("can only reorder 'children' for a \"gTree\"")
    order <- unique(order)
    oldOrder <- x$childrenOrder
    N <- length(oldOrder)
    if (is.character(order)) {
        # Convert to numeric
        order <- match(order, x$childrenOrder)
    }
    if (is.numeric(order)) {
        if (any(!is.finite(order)) ||
            !(all(order %in% 1:N))) {
            stop("Invalid 'order'")
        }
        if (back) {
            newOrder <- c(x$childrenOrder[order],
                          x$childrenOrder[-order])
        } else {
            newOrder <- c(x$childrenOrder[-order],
                          x$childrenOrder[order])
        }
    }
    x$childrenOrder <- newOrder
    x
}

# Reorder the children of a gTree on the display list
# (identified by a gPath)
# NOTE that it is possible for this operation to produce a grob
# that no longer draws (because it relies on another grob that
# used to be drawn before it, e.g., when the width of grob "b"
# is calculated from the width of grob "a")
# Do NOT allow reordering of grobs on the display list
# (it is not even clear what should happen in terms of reordering
#  grobs mixed with viewports PLUS the potential for ending up with
#  something that will not draw is pretty high)
# IF you want to reorder the grobs on the DL, do a grid.grab()
# first and then reorder the children of the resulting gTree
grid.reorder <- function(gPath, order, back=TRUE, grep=FALSE, redraw=TRUE) {
    grob <- grid.get(gPath, grep=grep)
    grid.set(gPath, reorderGrob(grob, order, back=back),
             grep=grep, redraw=redraw)
}

#  File src/library/grid/R/highlevel.R
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

######################################
## Example applications of grid    #
######################################

grid.strip <- function(label="whatever", range.full=c(0, 1),
                   range.thumb=c(.3, .6),
                   fill="#FFBF00", thumb="#FF8000",
                   vp=NULL) {
  diff.full <- diff(range.full)
  diff.thumb <- diff(range.thumb)
  if (!is.null(vp))
    pushViewport(vp)
  grid.rect(gp=gpar(col=NULL, fill=fill))
  grid.rect((range.thumb[1L] - range.full[1L])/diff.full, 0,
            diff.thumb/diff.full, 1,
            just=c("left", "bottom"),
            gp=gpar(col=NULL, fill=thumb))
  grid.text(as.character(label))
  if (!is.null(vp))
    popViewport()
}

grid.panel <- function(x = stats::runif(10), y = stats::runif(10),
                   zrange = c(0, 1), zbin = stats::runif(2),
                   xscale = extendrange(x),
                   yscale = extendrange(y),
                   axis.left = TRUE, axis.left.label = TRUE,
                   axis.right = FALSE, axis.right.label = TRUE,
                   axis.bottom = TRUE, axis.bottom.label = TRUE,
                   axis.top = FALSE, axis.top.label = TRUE,
                   vp=NULL) {
  if (!is.null(vp))
    pushViewport(vp)
  temp.vp <- viewport(layout=grid.layout(2, 1,
                         heights=unit(c(1, 1), c("lines", "null"))))
  pushViewport(temp.vp)
  strip.vp <- viewport(layout.pos.row=1, layout.pos.col=1,
                        xscale=xscale)
  pushViewport(strip.vp)
  grid.strip(range.full=zrange, range.thumb=zbin)
  grid.rect()
  if (axis.top)
    grid.xaxis(main=FALSE, label=axis.top.label)
  popViewport()
  plot.vp <- viewport(layout.pos.row=2, layout.pos.col=1,
                       xscale=xscale, yscale=yscale)
  pushViewport(plot.vp)
  grid.grill()
  grid.points(x, y, gp=gpar(col="blue"))
  grid.rect()
  if (axis.left)
    grid.yaxis(label=axis.left.label)
  if (axis.right)
    grid.yaxis(main=FALSE, label=axis.right.label)
  if (axis.bottom)
    grid.xaxis(label=axis.bottom.label)
  popViewport(2)
  if (!is.null(vp))
    popViewport()
  invisible(list(strip.vp = strip.vp, plot.vp = plot.vp))
}

grid.multipanel <- function(x = stats::runif(90), y = stats::runif(90),
                            z = stats::runif(90),
                            nplots = 9, nrow = 5, ncol = 2,
                            newpage = TRUE, vp = NULL)
{
    if (newpage)
        grid.newpage()
    if (!is.null(vp))
        pushViewport(vp)
    stopifnot(nplots >= 1)
    if((missing(nrow) || missing(ncol)) && !missing(nplots)) {
        ## determine 'smart' default ones
        rowcol <- grDevices::n2mfrow(nplots)
        nrow <- rowcol[1L]
        ncol <- rowcol[2L]
    }
    temp.vp <- viewport(layout = grid.layout(nrow, ncol))
    pushViewport(temp.vp)
    xscale <- extendrange(x)
    yscale <- extendrange(y)
    breaks <- seq.int(min(z), max(z), length.out = nplots + 1)
    for (i in 1L:nplots) {
        col <- (i - 1) %% ncol + 1
        row <- (i - 1) %/% ncol + 1
        panel.vp <- viewport(layout.pos.row = row,
                             layout.pos.col = col)
        panelx <- x[z >= breaks[i] & z <= breaks[i+1]]
        panely <- y[z >= breaks[i] & z <= breaks[i+1]]
        grid.panel(panelx, panely, range(z), c(breaks[i], breaks[i+1]),
                   xscale, yscale,
                   axis.left = (col == 1),
                   axis.right = (col == ncol || i == nplots),
                   axis.bottom = (row == nrow),
                   axis.top = (row == 1),
                   axis.left.label = is.even(row),
                   axis.right.label = is.odd(row),
                   axis.bottom.label = is.even(col),
                   axis.top.label = is.odd(col),
                   vp = panel.vp)
    }
    grid.text("Compression Ratio", unit(.5, "npc"), unit(-4, "lines"),
              gp = gpar(fontsize = 20),
              just = "center", rot = 0)
    grid.text("NOx (micrograms/J)", unit(-4, "lines"), unit(.5, "npc"),
              gp = gpar(fontsize = 20),
              just = "centre", rot = 90)
    popViewport()
    if (!is.null(vp))
        popViewport()
}

grid.show.layout <- function(l, newpage=TRUE,
                             bg="light grey",
                             cell.border="blue", cell.fill="light blue",
                             cell.label=TRUE, label.col="blue",
                             unit.col="red", vp=NULL) {
  if (!is.layout(l))
    stop("'l' must be a layout")
  if (newpage)
    grid.newpage()
  if (!is.null(vp))
    pushViewport(vp)
  grid.rect(gp=gpar(col=NULL, fill=bg))
  vp.mid <- viewport(0.5, 0.5, 0.8, 0.8, layout=l)
  pushViewport(vp.mid)
  grid.rect(gp=gpar(fill="white"))
  gp.red <- gpar(col=unit.col)
  for (i in 1L:l$nrow)
    for (j in 1L:l$ncol) {
      vp.inner <- viewport(layout.pos.row=i, layout.pos.col=j)
      pushViewport(vp.inner)
      grid.rect(gp=gpar(col=cell.border, fill=cell.fill))
      if (cell.label)
        grid.text(paste0("(", i, ", ", j, ")"), gp=gpar(col=label.col))
      if (j==1)
        # recycle heights if necessary
        grid.text(as.character("["(l$heights, i, top=FALSE)), gp=gp.red,
              just=c("right", "centre"),
              x=unit(-.05, "inches"), y=unit(.5, "npc"), rot=0)
      if (i==l$nrow)
        # recycle widths if necessary
        grid.text(as.character("["(l$widths, j, top=FALSE)), gp=gp.red,
              just=c("centre", "top"),
              x=unit(.5, "npc"), y=unit(-.05, "inches"), rot=0)
      if (j==l$ncol)
        # recycle heights if necessary
        grid.text(as.character("["(l$heights, i, top=FALSE)), gp=gp.red,
              just=c("left", "centre"),
              x=unit(1, "npc") + unit(.05, "inches"), y=unit(.5, "npc"),
              rot=0)
      if (i==1)
        # recycle widths if necessary
        grid.text(as.character("["(l$widths, j, top=FALSE)), gp=gp.red,
              just=c("centre", "bottom"),
              x=unit(.5, "npc"), y=unit(1, "npc") + unit(.05, "inches"),
              rot=0)
      popViewport()
    }
  popViewport()
  if (!is.null(vp))
    popViewport()
  ## return the viewport used to represent the parent viewport
  invisible(vp.mid)
}

grid.show.viewport <- function(v, parent.layout=NULL, newpage=TRUE,
                               border.fill="light grey",
                               vp.col="blue", vp.fill="light blue",
                               scale.col="red",
                               vp=NULL)
{
    ## if the viewport has a non-NULL layout.pos.row or layout.pos.col
    ## AND the viewport has a parent AND the parent has a layout
    ## represent the location of the viewport in the parent's layout ...
    if ((!is.null(v$layout.pos.row) || !is.null(v$layout.pos.col)) &&
        !is.null(parent.layout)) {
        if (!is.null(vp))
            pushViewport(vp)
        vp.mid <- grid.show.layout(parent.layout,
                                   cell.border="grey", cell.fill="white",
                                   cell.label=FALSE, newpage=newpage)
        pushViewport(vp.mid)
        pushViewport(v)
        gp.red <- gpar(col=scale.col)
        grid.rect(gp=gpar(col="blue", fill="light blue"))
        at <- grid.pretty(v$xscale)
        grid.xaxis(at=c(min(at), max(at)), gp=gp.red)
        at <- grid.pretty(v$yscale)
        grid.yaxis(at=c(min(at), max(at)), gp=gp.red)
        popViewport(2)
        if (!is.null(vp))
            popViewport()
    } else {
        if (newpage)
            grid.newpage()
        if (!is.null(vp))
            pushViewport(vp)
        grid.rect(gp=gpar(col=NULL, fill=border.fill))
        ## generate a viewport within the "top" viewport (vp) to represent the
        ## parent viewport of the viewport we are "show"ing (v).
        ## This is so that annotations at the edges of the
        ## parent viewport will be at least partially visible
        vp.mid <- viewport(0.5, 0.5, 0.8, 0.8)
        pushViewport(vp.mid)
        grid.rect(gp=gpar(fill="white"))
        x <- v$x
        y <- v$y
        w <- v$width
        h <- v$height
        pushViewport(v)
        grid.rect(gp=gpar(col=vp.col, fill=vp.fill))
        ## represent the "native" scale
        gp.red <- gpar(col=scale.col)
        at <- grid.pretty(v$xscale)
        grid.xaxis(at=c(min(at), max(at)), gp=gp.red)
        at <- grid.pretty(v$yscale)
        grid.yaxis(at=c(min(at), max(at)), gp=gp.red)
        grid.text(as.character(w), gp=gp.red,
                  just=c("centre", "bottom"),
                  x=unit(.5, "npc"), y=unit(1, "npc") + unit(.05, "inches"))
        grid.text(as.character(h), gp=gp.red,
                  just=c("left", "centre"),
                  x=unit(1, "npc") + unit(.05, "inches"), y=unit(.5, "npc"))
        popViewport()
        ## annotate the location and dimensions of the viewport
        grid.lines(unit.c(x, x), unit.c(unit(0, "npc"), y),
                   gp=gpar(col=scale.col, lty="dashed"))
        grid.lines(unit.c(unit(0, "npc"), x), unit.c(y, y),
                   gp=gpar(col=scale.col, lty="dashed"))
        grid.text(as.character(x), gp=gp.red,
                  just=c("centre", "top"),
                  x=x, y=unit(-.05, "inches"))
        grid.text(as.character(y), gp=gp.red,
                  just=c("bottom"),
                  x=unit(-.05, "inches"), y=y, rot=90)
        popViewport()
        if (!is.null(vp))
            popViewport()
    }
}

## old grid.legend <-
function(pch, labels, frame=TRUE,
                        hgap=unit(0.5, "lines"), vgap=unit(0.5, "lines"),
                        default.units="lines",
                        gp=gpar(), draw=TRUE,
                        vp=NULL) {
  ## Type checking on arguments
  labels <- as.character(labels)
  nkeys <- length(labels)
  if (length(pch) != nkeys)
    stop("'pch' and 'labels' not the same length")
  if (!is.unit(hgap))
    hgap <- unit(hgap, default.units)
  if (length(hgap) != 1)
    stop("'hgap' must be single unit")
  if (!is.unit(vgap))
    vgap <- unit(vgap, default.units)
  if (length(vgap) != 1)
    stop("'vgap' must be single unit")
  gf <- grid.frame(layout=grid.layout(nkeys, 2), vp=vp, gp=gp, draw=FALSE)
  for (i in 1L:nkeys) {
    if (i==1) {
      symbol.border <- unit.c(vgap, hgap, vgap, hgap)
      text.border <- unit.c(vgap, unit(0, "npc"), vgap, hgap)
    }
    else {
      symbol.border <- unit.c(vgap, hgap, unit(0, "npc"), hgap)
      text.border <- unit.c(vgap, unit(0, "npc"), unit(0, "npc"), hgap)
    }
    grid.pack(gf, grid.points(.5, .5, pch=pch[i], draw=FALSE),
              col=1, row=i, border=symbol.border,
              width=unit(1, "lines"), height=unit(1, "lines"),
              force.width=TRUE, draw=FALSE)
    grid.pack(gf, grid.text(labels[i], x=0, y=.5, just=c("left", "centre"),
                            draw=FALSE),
              col=2, row=i, border=text.border, draw=FALSE)
  }
  if (draw)
    grid.draw(gf)
  gf
}

grid.legend <-
function(pch, labels, frame=TRUE,
                        hgap=unit(0.5, "lines"), vgap=unit(0.5, "lines"),
                        default.units="lines",
                        gp=gpar(), draw=TRUE,
                        vp=NULL) {
  ## Type checking on arguments
  labels <- as.character(labels)
  nkeys <- length(labels)
  if (length(pch) != nkeys)
    stop("'pch' and 'labels' not the same length")
  if (!is.unit(hgap))
    hgap <- unit(hgap, default.units)
  if (length(hgap) != 1)
    stop("'hgap' must be single unit")
  if (!is.unit(vgap))
    vgap <- unit(vgap, default.units)
  if (length(vgap) != 1)
    stop("'vgap' must be single unit")
  legend.layout <-
    grid.layout(nkeys, 3,
                widths=unit.c(unit(2, "lines"),
                  max(unit(rep(1, nkeys), "strwidth", as.list(labels))),
                  hgap),
                heights=unit.pmax(unit(2, "lines"),
                  vgap + unit(rep(1, nkeys), "strheight", as.list(labels))))
  fg <- frameGrob(layout=legend.layout, vp=vp, gp=gp)
  for (i in 1L:nkeys) {
    fg <- placeGrob(fg, pointsGrob(.5, .5, pch=pch[i]), col=1, row=i)
    fg <- placeGrob(fg, textGrob(labels[i], x=0, y=.5,
                                 just=c("left", "centre")),
                    col=2, row=i)
  }
  if (draw)
    grid.draw(fg)
  fg
}

## Just a wrapper for a sample series of grid commands
grid.plot.and.legend <- function() {
  grid.newpage()
  top.vp <- viewport(width=0.8, height=0.8)
  pushViewport(top.vp)
  x <- stats::runif(10)
  y1 <- stats::runif(10)
  y2 <- stats::runif(10)
  pch <- 1L:3
  labels <- c("Girls", "Boys", "Other")
  lf <- frameGrob()
  plot <- gTree(children=gList(rectGrob(),
                  pointsGrob(x, y1, pch=1),
                  pointsGrob(x, y2, pch=2),
                  xaxisGrob(),
                  yaxisGrob()))
  lf <- packGrob(lf, plot)
  lf <- packGrob(lf, grid.legend(pch, labels, draw=FALSE),
                 height=unit(1,"null"), side="right")
  grid.draw(lf)
}

#  File src/library/grid/R/interactive.R
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


grid.locator <- function(unit="native") {
  location <- c(grid.Call(L_locator), 1)
  if (is.na(location[1L]))
    invisible(NULL)
  else {
    transform <- solve(current.transform())
    location <- (location %*% transform)
    # The inverse viewport transform is from device coordinates into
    # inches relative to the current viewport
    location <- unit(location/location[3L], "inches")
    list(x=convertX(location[1L], unit),
         y=convertY(location[2L], unit))
  }
}

#  File src/library/grid/R/just.R
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

# NOTE: the order of the strings in these conversion functions must
# match the order of the enums in ../src/lattice.h
# NOTE: the result of match() is an integer, but subtracting 1 converts
# to real => have to convert back to integer for passing to C code

# If the user specifies two values, the first is horizontal
# justification and the second is vertical

# If the user specifies only one value, use the following
# conversion table to give a second default value
#
# bottom  -->  centre, bottom
# left    -->  left,   centre
# right   -->  right,  centre
# top     -->  centre, top
# centre  -->  centre, centre

valid.charjust <- function(just) {
  if (length(just) == 1) {
    # single value may be any valid just
    just <- as.integer(match(just[1L], c("left", "right", "bottom", "top",
                                        "centre", "center")) - 1)
    if (any(is.na(just)))
      stop("invalid justification")
  } else if (length(just) > 1) {
    # first value must be one of "left", "right", "centre", or "center"
    just[1L] <- as.integer(match(just[1L], c("left", "right", "bottom", "top",
                                           "centre", "center")) - 1)
    if (!(just[1L] %in% c(0, 1, 4, 5)))
      stop("invalid horizontal justification")
    # second value must be one of "bottom", "top", "centre", or "center"
    just[2L] <- as.integer(match(just[2L], c("left", "right", "bottom", "top",
                                           "centre", "center")) - 1)
    if (!(just[2L] %in% c(2, 3, 4, 5)))
      stop("invalid vertical justification")
    just <- as.integer(just)
  }
  # Extend to length 2 if necessary
  if (length(just) < 2) {
    if (length(just) == 0)
      just <- c(4, 4)
    else
      just <- switch (just[1L] + 1,
                      c(0, 4), # left
                      c(1, 4), # right
                      c(4, 2), # bottom
                      c(4, 3), # top
                      c(4, 4), # centre
                      c(4, 4)) # center
  }
  # Convert to numeric
  just <- c(switch(just[1L] + 1, 0, 1, NA, NA, 0.5, 0.5),
            switch(just[2L] + 1, NA, NA, 0, 1, 0.5, 0.5))
  # Final paranoid check
  if (any(is.na(just)))
    stop("invalid justification")
  just
}

valid.numjust <- function(just) {
  if (length(just) == 0) {
    c(0.5, 0.5)
  } else {
    if (length(just) < 2) {
      c(just, 0.5)
    } else {
      just
    }
  }
}

valid.just <- function(just) {
  if (is.character(just))
    valid.charjust(just)
  else {
    valid.numjust(as.numeric(just))
  }
}

resolveHJust <- function(just, hjust) {
  if (is.null(hjust) || length(hjust) == 0)
    valid.just(just)[1L]
  else
    hjust
}

resolveVJust <- function(just, vjust) {
  if (is.null(vjust) || length(vjust) == 0)
    valid.just(just)[2L]
  else
    vjust
}
#  File src/library/grid/R/layout.R
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


is.layout <- function(l) {
  inherits(l, "layout")
}

# FIXME:  The internal C code now does a lot of recycling of
# unit values, units, and data.  Can some/most/all of the
# recycling stuff below be removed ?
valid.layout <- function(nrow, ncol, widths, heights, respect, just) {
  nrow <- as.integer(nrow)
  ncol <- as.integer(ncol)
  # make sure we're dealing with a unit object
  if (!is.logical(respect)) {
    respect <- as.matrix(respect)
    if (!is.matrix(respect) || any(dim(respect) != c(nrow, ncol)))
      stop("'respect' must be logical or an 'nrow' by 'ncol' matrix")
    }
  if (is.matrix(respect)) {
    respect.mat <- matrix(as.integer(respect),
                          dim(respect)[1L],
                          dim(respect)[2L])
    respect <- 2
  }
  else respect.mat <- matrix(0L, nrow, ncol)

  valid.just <- valid.just(just)
  l <- list(nrow = nrow, ncol = ncol,
            widths = widths, heights = heights,
            respect = respect, valid.respect=as.integer(respect),
            respect.mat = respect.mat,
            just=just, valid.just=valid.just)
  class(l) <- "layout"
  l
}

layout.torture <- function() {
  top.vp <- viewport(y=0, height=unit(1, "npc") - unit(1.5, "lines"),
                     just=c("centre", "bottom"))
  do.label <- function(label) {
    grid.rect(y=1, height=unit(1.5, "lines"),
              just=c("center", "top"))
    grid.text(label,
              y=unit(1, "npc") - unit(1, "lines"),
              gp=gpar(font=2))
  }
  # 1 = all relative widths and heights
  grid.show.layout(grid.layout(3,2), vp=top.vp)
  do.label("All dimensions relative -- no respect")
  # (1) with full respect
  grid.show.layout(grid.layout(3,2, respect=TRUE), vp=top.vp)
  do.label("All dimensions relative -- full respect")
  # (1) with partial respect
  grid.show.layout(grid.layout(3,2,respect=matrix(c(1,0,0,0,0,0), 3L, 2L, TRUE)),
                   vp=top.vp)
  do.label("All dimensions relative -- only top-left cell respected")
  # (1) with slightly weirder partial respect
  grid.show.layout(grid.layout(3,2,respect=matrix(c(1,0,0,0,0,1), 3L, 2L, TRUE)),
                   vp=top.vp)
  do.label("All relative -- top-left, bottom-right respected")
  # 2 = combination of absolute and relative widths and heights
  grid.show.layout(grid.layout(2, 3,
                       widths=unit(c(2,4,1), c("null", "cm", "null")),
                       heights=unit(c(6,4), c("cm", "null"))),
                   vp=top.vp)
  do.label("Absolute and relative -- no respect")
  # (2) with full respect
  grid.show.layout(grid.layout(2, 3,
                       widths=unit(c(2,4,1), c("null", "cm", "null")),
                       heights=unit(c(6,4), c("cm", "null")), respect=TRUE),
                   vp=top.vp)
  do.label("Absolute and relative -- full respect")
  # (2) with partial respect
  grid.show.layout(grid.layout(2, 3,
                       widths=unit(c(2,4,1), c("null", "cm", "null")),
                       heights=unit(c(6,4), c("cm", "null")),
                       respect=matrix(c(0,0,0,0,0,1), 2L, 3L, TRUE)),
                   vp=top.vp)
  do.label("Absolute and relative -- bottom-right respected")
}

# Return the region allocated by the layout of the current viewport
layoutRegion <- function(layout.pos.row=1, layout.pos.col=1) {
  region <- grid.Call(L_layoutRegion,
                      # This conversion matches the vailidity check in
                      # valid.viewport()
                      if (is.null(layout.pos.row)) layout.pos.row
                      else as.integer(rep(layout.pos.row, length.out=2)),
                      if (is.null(layout.pos.col)) layout.pos.col
                      else as.integer(rep(layout.pos.col, length.out=2)))
  list(left=unit(region[1L], "npc"),
       bottom=unit(region[2L], "npc"),
       width=unit(region[3L], "npc"),
       height=unit(region[4L], "npc"))
}

####################
# Accessors
####################

layout.nrow <- function(lay) {
  lay$nrow
}

layout.ncol <- function(lay) {
  lay$ncol
}

layout.widths <- function(lay) {
  lay$widths
}

layout.heights <- function(lay) {
  lay$heights
}

layout.respect <- function(lay) {
  switch(lay$respect + 1,
         FALSE,
         TRUE,
         lay$respect.mat)
}

####################
# Public constructor function
####################
grid.layout <- function (nrow = 1, ncol = 1,
                         widths = unit(rep(1, ncol), "null"),
                         heights = unit(rep(1, nrow), "null"),
                         default.units = "null",
                         respect = FALSE,
                         just="centre")
{
  if (!is.unit(widths))
    widths <- unit(widths, default.units)
  if (!is.unit(heights))
    heights <- unit(heights, default.units)
  valid.layout(nrow, ncol, widths, heights, respect, just)
}

####################
# Utility Functions
####################

dim.layout <- function(x) {
    c(x$nrow, x$ncol)
}
#  File src/library/grid/R/ls.R
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


# Code for listing objects in various grid "namespaces"
# (gTrees, vpTrees, and the grid display list)

# Return a "gridListing" object,
# ... either ...
# "gridVectorListing", which is just character vector,
#     "grobListing", or "vpListing", or "vpNameListing", or
#     "vpPopListing", or "vpUpListing",
# ... or ...
# "gridListListing", which is list of "gridListing" objects,
#      "grobListListing", or "vpListListing", ...
# ... or ...
# "gridTreeListing", which is list of parent "gridVectorListing" object
#                    plus children "gridListing" object,
#      "gTreeListing", or "vpTreeListing", or "vpNameTreeListing"
#      (vpStack or vpTree produces a "vpTreeListing").
#      (vpPath [depth > 1] produces a "vpNameTreeListing").
#
# "vpListListing", and all "gridTreeListing" objects have a "depth" attribute

# The print method will print these in some format, but by having
# a separate object, others can capture the result and format the
# printing themselves.

grid.ls <- function(x=NULL, grobs=TRUE, viewports=FALSE, fullNames=FALSE,
                    recursive=TRUE, print=TRUE, flatten=TRUE, ...) {
    # If 'x' is NULL, list the grobs on the DL
    if (is.null(x)) {
        listing <- gridListDL(grobs=grobs, viewports=viewports,
                              fullNames=fullNames, recursive=recursive)
    } else {
        listing <- gridList(x, grobs=grobs, viewports=viewports,
                            fullNames=fullNames, recursive=recursive)
    }
    if (flatten) {
        listing <- flattenListing(listing)
    }
    if (is.logical(print)) {
        if (print) {
            print(listing)
        }
    } else if (is.function(print)) {
        print(listing, ...)
    } else {
        stop("invalid 'print' argument")
    }
    invisible(listing)
}

gridListDL <- function(x, grobs=TRUE, viewports=FALSE,
                       fullNames=FALSE, recursive=TRUE) {
    display.list <- grid.Call(L_getDisplayList)
    dl.index <- grid.Call(L_getDLindex)
    result <- lapply(display.list[1L:dl.index], gridList,
                     grobs=grobs, viewports=viewports,
                     fullNames=fullNames, recursive=recursive)
    names(result) <- NULL
    class(result) <- c("gridListListing", "gridListing")
    result
}

gridList <- function(x, ...) {
    UseMethod("gridList")
}

gridList.default <- function(x, grobs=TRUE, viewports=FALSE,
                             fullNames=FALSE, recursive=TRUE) {
    if (is.null(x)) {
        # This handles empty slots in the display list
        result <- character()
        class(result) <- "gridListing"
    } else {
        stop("invalid object in 'listing'")
    }
    result
}

# Grob methods
gridList.grob <- function(x, grobs=TRUE, viewports=FALSE,
                          fullNames=FALSE, recursive=TRUE) {
    if (grobs) {
        if (fullNames) {
            result <- as.character(x)
        } else {
            result <- x$name
        }
        class(result) <- c("grobListing", "gridVectorListing", "gridListing")
    } else {
        result <- character()
        class(result) <- "gridListing"
    }
    if (viewports) {
        # Call makeContext() to get x$vp at drawing time
        x <- makeContext(x)
    }
    if (viewports && !is.null(x$vp)) {
        # Bit dodgy this bit
        # Emulates an "upViewport" on the DL
        n <- depth(x$vp)
        class(n) <- "up"
        result <- list(gridList(x$vp,
                               grobs=grobs, viewports=viewports,
                               fullNames=fullNames,
                               recursive=recursive),
                       result,
                       gridList(n,
                               grobs=grobs, viewports=viewports,
                               fullNames=fullNames,
                               recursive=recursive))
        class(result) <- c("gridListListing", "gridListing")
    }
    result
}

gridList.gList <- function(x, grobs=TRUE, viewports=FALSE,
                           fullNames=FALSE, recursive=TRUE) {
    # Allow for grobs=FALSE but viewports=TRUE
    if (grobs || viewports) {
        if (length(x) == 0L) {
            result <- character()
            class(result) <- "gridListing"
        } else {
            result <- lapply(x, gridList,
                             grobs=grobs, viewports=viewports,
                             fullNames=fullNames, recursive=recursive)
            class(result) <- c("gListListing", "gridListListing",
                               "gridListing")
        }
    } else {
        result <- character()
        class(result) <- "gridListing"
    }
    result
}

gridList.gTree <- function(x, grobs=TRUE, viewports=FALSE,
                           fullNames=FALSE, recursive=TRUE) {
    if (fullNames) {
        name <- as.character(x)
    } else {
        name <- x$name
    }
    class(name) <- c("grobListing", "gridVectorListing", "gridListing")
    if (viewports) {
        # Call makeContext() to get x$vp and x$childrenvp at drawing time
        x <- makeContext(x)
    }
    if (recursive) {
        # Allow for grobs=FALSE but viewports=TRUE
        result <- gridList(x$children[x$childrenOrder],
                           grobs=grobs, viewports=viewports,
                           fullNames=fullNames, recursive=recursive)
        if (viewports && !is.null(x$childrenvp)) {
            # Bit dodgy this bit
            # Emulates an "upViewport" on the DL
            n <- depth(x$childrenvp)
            class(n) <- "up"
            result <- list(gridList(x$childrenvp,
                                    grobs=grobs, viewports=viewports,
                                    fullNames=fullNames,
                                    recursive=recursive),
                           gridList(n,
                                    grobs=grobs, viewports=viewports,
                                    fullNames=fullNames,
                                    recursive=recursive),
                           result)
            class(result) <- c("gridListListing", "gridListing")
        }
        if (grobs) {
            result <- list(parent=name,
                           children=result)
            class(result) <- c("gTreeListing", "gridTreeListing",
                               "gridListing")
        } else if (!viewports) {
            result <- character()
            class(result) <- "gridListing"
        }
    } else {
        if (grobs) {
            result <- name
        } else {
            result <- character()
            class(result) <- "gridListing"
        }
    }
    if (viewports && !is.null(x$vp)) {
        # Bit dodgy this bit
        # Emulates an "upViewport" on the DL
        n <- depth(x$vp)
        class(n) <- "up"
        result <- list(gridList(x$vp,
                                grobs=grobs, viewports=viewports,
                                fullNames=fullNames,
                                recursive=recursive),
                       result,
                       gridList(n,
                                grobs=grobs, viewports=viewports,
                                fullNames=fullNames,
                                recursive=recursive))
        class(result) <- c("gridListListing", "gridListing")
    }
    result
}

# Viewport methods
gridList.viewport <- function(x, grobs=TRUE, viewports=FALSE,
                              fullNames=FALSE, recursive=TRUE) {
    if (viewports) {
        if (fullNames) {
            result <- as.character(x)
        } else {
            result <- x$name
        }
        class(result) <- c("vpListing", "gridVectorListing", "gridListing")
    } else {
        result <- character()
        class(result) <- "gridListing"
    }
    result
}

# ... are arugments to gridList
listvpListElement <- function(x, ...) {
    n <- depth(x)
    class(n) <- "up"
    result <- list(gridList(x, ...),
                   gridList(n, ...))
    class(result) <- c("gridListListing", "gridListing")
    result
}

gridList.vpList <- function(x, grobs=TRUE, viewports=FALSE,
                            fullNames=FALSE, recursive=TRUE) {
    if (viewports) {
        if (length(x) == 0L) {
            result <- character()
            class(result) <- "gridListing"
        } else if (length(x) == 1L) {
            result <- gridList(x[[1L]],
                              grobs=grobs, viewports=viewports,
                              fullNames=fullNames,
                              recursive=recursive)
        } else {
            result <- c(lapply(x[-length(x)], listvpListElement,
                               grobs=grobs, viewports=viewports,
                               fullNames=fullNames,
                               recursive=recursive),
                        list(gridList(x[[length(x)]],
                                     grobs=grobs, viewports=viewports,
                                     fullNames=fullNames,
                                     recursive=recursive)))
            attr(result, "depth") <- depth(x[[length(x)]])
            class(result) <- c("vpListListing", "gridListListing",
                               "gridListing")
        }
    } else {
        result <- character()
        class(result) <- "gridListing"
    }
    result
}

gridList.vpStack <- function(x, grobs=TRUE, viewports=FALSE,
                             fullNames=FALSE, recursive=TRUE) {
    if (viewports) {
        if (length(x) == 0L) {
            result <- character()
            class(result) <- "gridListing"
        } else if (length(x) == 1L || !recursive) {
            result <- gridList(x[[1L]],
                               grobs=grobs, viewports=viewports,
                               fullNames=fullNames, recursive=recursive)
        } else {
            theRest <- x[-1L]
            class(theRest) <- "vpStack"
            result <- gridList(theRest,
                               grobs=grobs, viewports=viewports,
                               fullNames=fullNames,
                               recursive=recursive)
            result <- list(parent=gridList(x[[1L]],
                             grobs=grobs, viewports=viewports,
                             fullNames=fullNames,
                             recursive=recursive),
                           children=result)
            attr(result, "depth") <- depth(x)
            class(result) <- c("vpTreeListing", "gridTreeListing",
                               "gridListing")
        }
    } else {
        result <- character()
        class(result) <- "gridListing"
    }
    result
}

gridList.vpTree <- function(x, grobs=TRUE, viewports=FALSE,
                            fullNames=FALSE, recursive=TRUE) {
    if (viewports) {
        if (recursive) {
            result <- gridList(x$children,
                               grobs=grobs, viewports=viewports,
                               fullNames=fullNames, recursive=recursive)
            # Parent can only be a plain viewport
            result <- list(parent=gridList(x$parent,
                             grobs=grobs, viewports=viewports,
                             fullNames=fullNames,
                             recursive=recursive),
                           children=result)
            attr(result, "depth") <- depth(x$children) + 1
            class(result) <- c("vpTreeListing", "gridTreeListing",
                               "gridListing")
        } else {
            result <- gridList(x$parent,
                               grobs=grobs, viewports=viewports,
                               fullNames=fullNames, recursive=recursive)
        }
    } else {
        result <- character()
        class(result) <- "gridListing"
    }
    result
}

# This handles downViewports in the display list
gridList.vpPath <- function(x, grobs=TRUE, viewports=FALSE,
                            fullNames=FALSE, recursive=TRUE) {
    if (viewports) {
        # Have to account for top-level downViewports that are
        # non-strict (i.e., they could navigate down quite a long way)
        # In particular, when the vpPath navigates down more
        # levels than there are names in the vpPath
        recordedDepth <- attr(x, "depth")
        if (!is.null(recordedDepth) && recordedDepth != depth(x)) {
            # In this case, need to prepend a fake path on the front
            # so that subsequent upViewport()s will work
            x <- vpPathFromVector(c(rep("...", recordedDepth - depth(x)),
                                    explodePath(as.character(x))))
        }
        # This would be simpler if paths were kept as vectors
        # but that redesign is a bit of an undertaking
        if (depth(x) == 1) {
            if (fullNames) {
                result <- paste0("downViewport[", x$name, "]")
            } else {
                result <- x$name
            }
            class(result) <- c("vpNameListing", "gridVectorListing",
                               "gridListing")
        } else if (depth(x) == 2) {
            result <- gridList(vpPath(x$name),
                               grobs=grobs, viewports=viewports,
                               fullNames=fullNames,
                               recursive=recursive)
            result <- list(parent=gridList(vpPath(x$path),
                             grobs=grobs, viewports=viewports,
                             fullNames=fullNames,
                             recursive=recursive),
                           children=result)
            attr(result, "depth") <- depth(x)
            # Inherit updateVPDepth and updateVPPath methods
            # from vpTreeListing
            class(result) <- c("vpNameTreeListing", "vpTreeListing",
                               "gridTreeListing", "gridListing")
        } else {
            path <- explodePath(x$path)
            result <- gridList(vpPathFromVector(c(path[-1L], x$name)),
                               grobs=grobs, viewports=viewports,
                               fullNames=fullNames,
                               recursive=recursive)
            result <- list(parent=gridList(vpPath(path[1L]),
                             grobs=grobs, viewports=viewports,
                             fullNames=fullNames,
                             recursive=recursive),
                           children=result)
            attr(result, "depth") <- depth(x)
            # Inherit updateVPDepth and updateVPPath methods
            # from vpTreeListing
            class(result) <- c("vpNameTreeListing", "vpTreeListing",
                               "gridTreeListing", "gridListing")
        }
    } else {
        result <- character()
        class(result) <- "gridListing"
    }
    result
}

# This handles popViewports in the display list
gridList.pop <- function(x, grobs=TRUE, viewports=FALSE,
                         fullNames=FALSE, recursive=TRUE) {
    if (viewports) {
        result <- as.character(x)
        if (fullNames) {
            result <- paste0("popViewport[", result, "]")
        }
        class(result) <- c("vpPopListing", "gridVectorListing", "gridListing")
    } else {
        result <- character()
        class(result) <- "gridListing"
    }
    result
}

# This handles upViewports in the display list
gridList.up <- function(x, grobs=TRUE, viewports=FALSE,
                        fullNames=FALSE, recursive=TRUE) {
    if (viewports) {
        result <- as.character(x)
        if (fullNames) {
            result <- paste0("upViewport[", result, "]")
        }
        class(result) <- c("vpUpListing", "gridVectorListing", "gridListing")
    } else {
        result <- character()
        class(result) <- "gridListing"
    }
    result
}

######################
# flatten methods for gridListing objects
######################

incDepth <- function(depth, n=1) {
    depth + n
}

decrDepth <- function(depth, x) {
    n <- as.numeric(gsub("^.+\\[", "",
                         gsub("\\]$", "",
                              as.character(x))))
    depth - n
}

# updateDepth modifies depth from sibling to sibling
# (flatListing methods take care of parent to child updates of depth)
updateGDepth <- function(x, gdepth) {
    UseMethod("updateGDepth")
}

updateGDepth.default <- function(x, gdepth) {
    gdepth
}

updateVPDepth <- function(x, vpdepth) {
    UseMethod("updateVPDepth")
}

updateVPDepth.default <- function(x, vpdepth) {
    vpdepth
}

updateVPDepth.vpListing <- function(x, vpdepth) {
    incDepth(vpdepth)
}

updateVPDepth.vpNameListing <- function(x, vpdepth) {
    incDepth(vpdepth)
}

updateVPDepth.vpListListing <- function(x, vpdepth) {
    incDepth(vpdepth, attr(x, "depth"))
}

updateVPDepth.vpUpListing <- function(x, vpdepth) {
    decrDepth(vpdepth, x)
}

updateVPDepth.vpPopListing <- function(x, vpdepth) {
    decrDepth(vpdepth, x)
}

updateVPDepth.vpTreeListing <- function(x, vpdepth) {
    incDepth(vpdepth, attr(x, "depth"))
}

incPath <- function(oldpath, addition) {
    if (nchar(oldpath) > 0) {
        paste0(oldpath, .grid.pathSep, as.character(addition))
    } else {
        as.character(addition)
    }
}

decrPath <- function(oldpath, x) {
    bits <- strsplit(oldpath, .grid.pathSep)[[1L]]
    n <- as.numeric(gsub("^.+\\[", "",
                         gsub("\\]$", "",
                              as.character(x))))
    if ((m <- (length(bits) - n)) == 0L) {
        ""
    } else {
	paste(bits[seq_len(m)], collapse=.grid.pathSep)
    }
}

updateGPath <- function(x, gpath) {
    UseMethod("updateGPath")
}

updateGPath.default <- function(x, gpath) {
    gpath
}

updateVPPath <- function(x, vppath) {
    UseMethod("updateVPPath")
}

updateVPPath.default <- function(x, vppath) {
    vppath
}

updateVPPath.vpListing <- function(x, vppath) {
    incPath(vppath, x)
}

updateVPPath.vpNameListing <- function(x, vppath) {
    incPath(vppath, x)
}

updateVPPath.vpListListing <- function(x, vppath) {
    incPath(vppath, x[[length(x)]])
}

updateVPPath.vpUpListing <- function(x, vppath) {
    decrPath(vppath, x)
}

updateVPPath.vpPopListing <- function(x, vppath) {
    decrPath(vppath, x)
}

updateVPPath.vpTreeListing <- function(x, vppath) {
    incPath(vppath,
            paste0(updateVPPath(x$parent, ""), .grid.pathSep,
                   updateVPPath(x$children, "")))
}

flatListing <- function(x, gDepth=0, vpDepth=0, gPath="", vpPath="") {
    UseMethod("flatListing")
}

flatListing.gridListing <- function(x, gDepth=0, vpDepth=0,
                                    gPath="", vpPath="") {
    if (length(x)) {
        list(name=as.character(x),
             gDepth=gDepth,
             vpDepth=vpDepth,
             gPath=gPath,
             vpPath=vpPath,
             type=class(x)[1L])
    } else {
        list(name=character(),
             gDepth=numeric(),
             vpDepth=numeric(),
             gPath=character(),
             vpPath=character(),
             type=character())
    }
}

flatListing.gTreeListing <- function(x, gDepth=0, vpDepth=0,
                                     gPath="", vpPath="") {
    # Increase gDepth and gPath
    flatChildren <- flatListing(x$children, incDepth(gDepth, 1), vpDepth,
                                incPath(gPath, x$parent), vpPath)
    list(name=c(as.character(x$parent), flatChildren$name),
         gDepth=c(gDepth, flatChildren$gDepth),
         vpDepth=c(vpDepth, flatChildren$vpDepth),
         gPath=c(gPath, flatChildren$gPath),
         vpPath=c(vpPath, flatChildren$vpPath),
         type=c(class(x)[1L], flatChildren$type))
}

OLDflatListing.vpTreeListing <- function(x, gDepth=0, vpDepth=0,
                                      gPath="", vpPath="") {
    # Increase vpDepth and vpPath
    flatChildren <- flatListing(x$children, gDepth, incDepth(vpDepth, 1),
                                gPath, incPath(vpPath, x$parent))
    list(name=c(as.character(x$parent), flatChildren$name),
         gDepth=c(gDepth, flatChildren$gDepth),
         vpDepth=c(vpDepth, flatChildren$vpDepth),
         gPath=c(gPath, flatChildren$gPath),
         vpPath=c(vpPath, flatChildren$vpPath),
         type=c(class(x)[1L], flatChildren$type))
}

flatListing.vpTreeListing <- function(x, gDepth=0, vpDepth=0,
                                      gPath="", vpPath="") {
    flatParent <- flatListing(x$parent, gDepth, vpDepth,
                              gPath, vpPath)
    depth <- attr(x$parent, "depth")
    if (is.null(depth)) {
        depth <- 1
    }
    # Increase vpDepth and vpPath
    flatChildren <- flatListing(x$children, gDepth, incDepth(vpDepth, depth),
                                gPath, updateVPPath(x$parent, vpPath))
    list(name=c(flatParent$name, flatChildren$name),
         gDepth=c(flatParent$gDepth, flatChildren$gDepth),
         vpDepth=c(flatParent$vpDepth, flatChildren$vpDepth),
         gPath=c(flatParent$gPath, flatChildren$gPath),
         vpPath=c(flatParent$vpPath, flatChildren$vpPath),
         type=c(flatParent$type, flatChildren$type))
}

flatListing.vpNameTreeListing <- function(x, gDepth=0, vpDepth=0,
                                      gPath="", vpPath="") {
    # Increase vpDepth and vpPath
    flatChildren <- flatListing(x$children, gDepth, incDepth(vpDepth, 1),
                                gPath, incPath(vpPath, x$parent))
    list(name=c(as.character(x$parent), flatChildren$name),
         gDepth=c(gDepth, flatChildren$gDepth),
         vpDepth=c(vpDepth, flatChildren$vpDepth),
         gPath=c(gPath, flatChildren$gPath),
         vpPath=c(vpPath, flatChildren$vpPath),
         type=c(class(x)[1L], flatChildren$type))
}

flatListing.gridListListing <- function(x, gDepth=0, vpDepth=0,
                                        gPath="", vpPath="") {
    n <- length(x)
    listListing <- list(name=character(),
                        gDepth=numeric(),
                        vpDepth=numeric(),
                        gPath=character(),
                        vpPath=character(),
                        type=character())
    for (i in 1L:n) {
        componentListing <- flatListing(x[[i]], gDepth, vpDepth,
                                        gPath, vpPath)
        listListing$name <- c(listListing$name,
                              componentListing$name)
        listListing$gDepth <- c(listListing$gDepth,
                                componentListing$gDepth)
        listListing$vpDepth <- c(listListing$vpDepth,
                                 componentListing$vpDepth)
        listListing$gPath <- c(listListing$gPath,
                               componentListing$gPath)
        listListing$vpPath <- c(listListing$vpPath,
                                componentListing$vpPath)
        listListing$type <- c(listListing$type,
                              componentListing$type)
        gPath <- updateGPath(x[[i]], gPath)
        vpPath <- updateVPPath(x[[i]], vpPath)
        gDepth <- updateGDepth(x[[i]], gDepth)
        vpDepth <- updateVPDepth(x[[i]], vpDepth)
    }
    listListing
}

flattenListing <- function(x) {
    listing <- flatListing(x)
    class(listing) <- "flatGridListing"
    listing
}

print.flatGridListing <- function(x, ...) {
    nestedListing(x, ...)
    invisible(x)
}

######################
# Print functions for flatGridListings
######################

nestedListing <- function(x, gindent="  ", vpindent=gindent) {

    makePrefix <- function(indent, depth) {
        indents <- rep(indent, length(depth))
        indents <- mapply(rep, indents, depth)
        sapply(indents, paste, collapse="")
    }

    if (!inherits(x, "flatGridListing"))
        stop("invalid listing")
    cat(paste0(makePrefix(gindent, x$gDepth),
               makePrefix(vpindent, x$vpDepth),
               x$name),
        sep = "\n")
}

pathListing <- function(x, gvpSep=" | ", gAlign=TRUE) {

    appendToPrefix <- function(path, name) {
        emptyPath <- nchar(path) == 0
        ifelse(emptyPath,
               name,
               paste(path, name, sep = .grid.pathSep))
    }

    padPrefix <- function(path, maxLen) {
        numSpaces <- maxLen - nchar(path)
        if (length(path) == 1L) {
            paste0(path, paste(rep.int(" ", numSpaces), collapse=""))
        } else {
            padding <- rep(" ", length(path))
            padding <- mapply(rep.int, padding, numSpaces)
            paste0(path, sapply(padding, paste, collapse=""))
        }
    }

    if (!inherits(x, "flatGridListing"))
        stop("invalid 'listing'")
    vpListings <- seq_along(x$name) %in% grep("^vp", x$type)
    paths <- x$vpPath
    # Only if viewport listings
    if (sum(vpListings) > 0) {
        paths[vpListings] <- appendToPrefix(paths[vpListings],
                                            x$name[vpListings])
        # If viewports are shown, then allow extra space before grobs
        maxLen <- max(nchar(paths[vpListings]))
    }
    else
	maxLen <- max(nchar(paths))

    # Only if grob listings
    if (sum(!vpListings) > 0) {
        if (gAlign) {
            paths[!vpListings] <- padPrefix(paths[!vpListings], maxLen)
        }
        paths[!vpListings] <- paste0(paths[!vpListings],
				     gvpSep,
				     appendToPrefix(x$gPath[!vpListings],
						    x$name[!vpListings]))
    }
    cat(paths, sep = "\n")
}

grobPathListing <- function(x, ...) {
    subset <- grep("^g", x$type)
    if (length(subset)) {
        cl <- class(x)
        subListing <- lapply(x, "[", subset)
        class(subListing) <- cl
        pathListing(subListing, ...)
    }
}


# API to access detailed text metric info
#
#  Copyright (C) 1995-2012 The R Core Team

# This first function does NOT return a "unit" object
# It is just access to font metric info in the calling context
# (similar to the convert*() functions, with corresponding caveats on use)

grid.textMetric <- function(string) {

}

# It should be possible to define units like "strascent" and "strdescent"
#  File src/library/grid/R/origin.R
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

valid.origin <- function(origin) {
  origin <- as.integer(match(origin,
                             c("bottom.left", "top.left",
                               "bottom.right", "top.right")) - 1)
  if (any(is.na(origin)))
    stop("invalid 'origin'")
  origin
}

origin.left <- function(origin) {
  switch (origin,
          bottom.left = TRUE,
          bottom.right = FALSE,
          top.left = TRUE,
          top.right = FALSE)
}

origin.right <- function(origin) {
  switch (origin,
          bottom.left = FALSE,
          bottom.right = TRUE,
          top.left = FALSE,
          top.right = TRUE)
}

origin.bottom <- function(origin) {
  switch (origin,
          bottom.left = TRUE,
          bottom.right = TRUE,
          top.left = FALSE,
          top.right = FALSE)
}

origin.top <- function(origin) {
  switch (origin,
          bottom.left = FALSE,
          bottom.right = FALSE,
          top.left = TRUE,
          top.right = TRUE)
}

swap.origin.horizontal <- function(origin) {
  switch (origin,
          bottom.left = "bottom.right",
          bottom.right = "bottom.left",
          top.left = "top.right",
          top.right = "top.left")
}

swap.origin.vertical <- function(origin) {
  switch (origin,
          bottom.left = "top.left",
          bottom.right = "top.right",
          top.left = "bottom.left",
          top.right = "bottom.right")
}
#  File src/library/grid/R/primitives.R
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


# Function that creates a description of an arrow head
# to add to a line
arrow <- function(angle=30, length=unit(0.25, "inches"),
                  ends="last", type="open") {
    angle <- as.numeric(angle)
    if (!is.unit(length))
        stop("'length' must be a 'unit' object")
    ends <- as.integer(match(ends, c("first", "last", "both")))
    type <- as.integer(match(type, c("open", "closed")))
    if (any(is.na(ends)) || any(is.na(type)) ||
        length(ends) == 0 || length(type) == 0)
        stop("invalid 'ends' or 'type' argument")
    a <- list(angle=angle, length=length,
              ends=ends, type=type)
    class(a) <- "arrow"
    a
}

length.arrow <- function(x) {
    max(do.call("max", lapply(x, length)),
                length(x$length))
}

rep.arrow <- function(x, ...) {
    maxn <- length(x)
    newa <- list(angle=rep(x$angle, length.out=maxn),
                 length=rep(x$length, length.out=maxn),
                 ends=rep(x$ends, length.out=maxn),
                 type=rep(x$type, length.out=maxn))
    newa <- lapply(newa, rep, ...)
    class(newa) <- "arrow"
    newa
}

# Method for subsetting "arrow" objects
`[.arrow` <- function(x, index, ...) {
    if (length(index) == 0)
        return(NULL)
    maxn <- length(x)
    newa <- list(angle=rep(x$angle, length.out=maxn),
                 length=rep(x$length, length.out=maxn),
                 ends=rep(x$ends, length.out=maxn),
                 type=rep(x$type, length.out=maxn))
    newa <- lapply(X = newa, FUN = "[", index, ...)
    class(newa) <- "arrow"
    newa
}

######################################
# move-to and line-to primitives
######################################
validDetails.move.to <- function(x) {
  if (!is.unit(x$x) ||
      !is.unit(x$y))
    stop("'x' and 'y' must be units")
  # Make sure that x and y are of length 1
  if (length(x$x) > 1 | length(x$y) > 1)
    stop("'x' and 'y' must have length 1")
  x
}

drawDetails.move.to <- function(x, recording=TRUE) {
  grid.Call.graphics(L_moveTo, x$x, x$y)
}

moveToGrob <- function(x=0, y=0,
                       default.units="npc",
                       name=NULL, vp=NULL) {
  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)
  grob(x=x, y=y,
       name=name, vp=vp, cl="move.to")
}

grid.move.to <- function(x=0, y=0,
                         default.units="npc",
                         name=NULL, draw=TRUE, vp=NULL) {
  mtg <- moveToGrob(x=x, y=y, default.units=default.units,
                    name=name, vp=vp)
  if (draw)
    grid.draw(mtg)
  invisible(mtg)
}

validDetails.line.to <- function(x) {
  if (!is.unit(x$x) ||
      !is.unit(x$y))
    stop("'x' and 'y' must be units")
  # Make sure that x and y are of length 1
  if (length(x$x) > 1 | length(x$y) > 1)
    stop("'x' and 'y' must have length 1")
  if (!(is.null(x$arrow) || inherits(x$arrow, "arrow")))
      stop("invalid 'arrow' argument")
  x
}

drawDetails.line.to <- function(x, recording=TRUE) {
  grid.Call.graphics(L_lineTo, x$x, x$y, x$arrow)
}

lineToGrob <- function(x=1, y=1,
                       default.units="npc",
                       arrow=NULL,
                       name=NULL, gp=gpar(), vp=NULL) {
  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)
  grob(x=x, y=y, arrow=arrow,
       name=name, gp=gp, vp=vp, cl="line.to")
}

grid.line.to <- function(x=1, y=1,
                         default.units="npc",
                         arrow=NULL,
                         name=NULL, gp=gpar(), draw=TRUE, vp=NULL) {
  ltg <- lineToGrob(x=x, y=y, default.units=default.units, arrow=arrow,
                    name=name, gp=gp, vp=vp)
  if (draw)
    grid.draw(ltg)
  invisible(ltg)
}

######################################
# LINES primitive
######################################
validDetails.lines <- function(x) {
  if (!is.unit(x$x) ||
      !is.unit(x$y))
    stop("'x' and 'y' must be units")
  if (!(is.null(x$arrow) || inherits(x$arrow, "arrow")))
      stop("invalid 'arrow' argument")
  x
}

drawDetails.lines <- function(x, recording=TRUE) {
    grid.Call.graphics(L_lines, x$x, x$y,
                       list(as.integer(1L:max(length(x$x), length(x$y)))),
                       x$arrow)
}

xDetails.lines <- function(x, theta) {
    bounds <- grid.Call(L_locnBounds, x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[1L], "inches")
}

yDetails.lines <- function(x, theta) {
    bounds <- grid.Call(L_locnBounds, x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[2L], "inches")
}

widthDetails.lines <- function(x) {
  bounds <- grid.Call(L_locnBounds, x$x, x$y, 0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[3L], "inches")
}

heightDetails.lines <- function(x) {
  bounds <- grid.Call(L_locnBounds, x$x, x$y, 0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[4L], "inches")
}

linesGrob <- function(x=unit(c(0, 1), "npc"),
                      y=unit(c(0, 1), "npc"),
                      default.units="npc",
                      arrow=NULL,
                      name=NULL, gp=gpar(), vp=NULL) {
  # Allow user to specify unitless vector;  add default units
  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)
  grob(x=x, y=y,
       arrow=arrow, name=name, gp=gp, vp=vp, cl="lines")
}

grid.lines <- function(x=unit(c(0, 1), "npc"),
                       y=unit(c(0, 1), "npc"),
                       default.units="npc",
                       arrow=NULL,
                       name=NULL, gp=gpar(), draw=TRUE, vp=NULL) {
  lg <- linesGrob(x=x, y=y,
                  default.units=default.units, arrow=arrow,
                  name=name, gp=gp, vp=vp)
  if (draw)
    grid.draw(lg)
  invisible(lg)
}

######################################
# POLYLINES primitive
######################################
# Very similar to LINES primitive, but allows
# multiple polylines via 'id' and 'id.lengths' args
# as per POLYGON primitive
validDetails.polyline <- function(x) {
  if (!is.unit(x$x) ||
      !is.unit(x$y))
      stop("'x' and 'y' must be units")
  if (!is.null(x$id) && !is.null(x$id.lengths))
      stop("it is invalid to specify both 'id' and 'id.lengths'")
  if (length(x$x) != length(x$y))
      stop("'x' and 'y' must be same length")
  if (!is.null(x$id) && (length(x$id) != length(x$x)))
      stop("'x' and 'y' and 'id' must all be same length")
  if (!is.null(x$id))
      x$id <- as.integer(x$id)
  if (!is.null(x$id.lengths) && (sum(x$id.lengths) != length(x$x)))
      stop("'x' and 'y' and 'id.lengths' must specify same overall length")
  if (!is.null(x$id.lengths))
      x$id.lengths <- as.integer(x$id.lengths)
  if (!(is.null(x$arrow) || inherits(x$arrow, "arrow")))
      stop("invalid 'arrow' argument")
  x
}

drawDetails.polyline <- function(x, recording=TRUE) {
    if (is.null(x$id) && is.null(x$id.lengths))
        grid.Call.graphics(L_lines, x$x, x$y,
                           list(as.integer(seq_along(x$x))),
                           x$arrow)
    else {
        if (is.null(x$id)) {
            n <- length(x$id.lengths)
            id <- rep(1L:n, x$id.lengths)
        } else {
            n <- length(unique(x$id))
            id <- x$id
        }
        index <- split(as.integer(seq_along(x$x)), id)
        grid.Call.graphics(L_lines, x$x, x$y, index, x$arrow)
    }
}

xDetails.polyline <- function(x, theta) {
    bounds <- grid.Call(L_locnBounds, x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[1L], "inches")
}

yDetails.polyline <- function(x, theta) {
    bounds <- grid.Call(L_locnBounds, x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[2L], "inches")
}

widthDetails.polyline <- function(x) {
    bounds <- grid.Call(L_locnBounds, x$x, x$y, 0)
    if (is.null(bounds))
        unit(0, "inches")
    else
        unit(bounds[3L], "inches")
}

heightDetails.polyline <- function(x) {
    bounds <- grid.Call(L_locnBounds, x$x, x$y, 0)
    if (is.null(bounds))
        unit(0, "inches")
    else
        unit(bounds[4L], "inches")
}

polylineGrob <- function(x=unit(c(0, 1), "npc"),
                         y=unit(c(0, 1), "npc"),
                         id=NULL, id.lengths=NULL,
                         default.units="npc",
                         arrow=NULL,
                         name=NULL, gp=gpar(), vp=NULL) {
    # Allow user to specify unitless vector;  add default units
    if (!is.unit(x))
        x <- unit(x, default.units)
    if (!is.unit(y))
        y <- unit(y, default.units)
    grob(x=x, y=y, id=id, id.lengths=id.lengths,
         arrow=arrow, name=name, gp=gp, vp=vp, cl="polyline")
}

grid.polyline <- function(...) {
    grid.draw(polylineGrob(...))
}

######################################
# SEGMENTS primitive
######################################
validDetails.segments <- function(x) {
  if (!is.unit(x$x0) || !is.unit(x$x1) ||
      !is.unit(x$y0) || !is.unit(x$y1))
    stop("'x0', 'y0', 'x1', and 'y1' must be units")
  if (!(is.null(x$arrow) || inherits(x$arrow, "arrow")))
      stop("invalid 'arrow' argument")
  x
}

drawDetails.segments <- function(x, recording=TRUE) {
  grid.Call.graphics(L_segments, x$x0, x$y0, x$x1, x$y1, x$arrow)
}

segmentBounds <- function(x, theta) {
    n <- max(length(x$x0), length(x$x1),
             length(x$y0), length(x$y1))
    x0 <- rep(x$x0, length.out=n)
    x1 <- rep(x$x1, length.out=n)
    y0 <- rep(x$y0, length.out=n)
    y1 <- rep(x$y1, length.out=n)
    grid.Call(L_locnBounds, unit.c(x0, x1), unit.c(y0, y1), theta)
}

xDetails.segments <- function(x, theta) {
    bounds <- segmentBounds(x, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[1L], "inches")
}

yDetails.segments <- function(x, theta) {
    bounds <- segmentBounds(x, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[2L], "inches")
}

widthDetails.segments <- function(x) {
    bounds <- segmentBounds(x, 0)
    if (is.null(bounds))
        unit(0, "inches")
    else
        unit(bounds[3L], "inches")
}

heightDetails.segments <- function(x) {
    bounds <- segmentBounds(x, 0)
    if (is.null(bounds))
        unit(0, "inches")
    else
        unit(bounds[4L], "inches")
}

segmentsGrob <- function(x0=unit(0, "npc"), y0=unit(0, "npc"),
                         x1=unit(1, "npc"), y1=unit(1, "npc"),
                         default.units="npc",
                         arrow=NULL,
                         name=NULL, gp=gpar(), vp=NULL) {
  # Allow user to specify unitless vector;  add default units
  if (!is.unit(x0))
    x0 <- unit(x0, default.units)
  if (!is.unit(x1))
    x1 <- unit(x1, default.units)
  if (!is.unit(y0))
    y0 <- unit(y0, default.units)
  if (!is.unit(y1))
    y1 <- unit(y1, default.units)
  grob(x0=x0, y0=y0, x1=x1, y1=y1, arrow=arrow, name=name, gp=gp, vp=vp,
       cl="segments")
}

grid.segments <- function(x0=unit(0, "npc"), y0=unit(0, "npc"),
                          x1=unit(1, "npc"), y1=unit(1, "npc"),
                          default.units="npc",
                          arrow=NULL,
                          name=NULL, gp=gpar(), draw=TRUE, vp=NULL) {
  sg <- segmentsGrob(x0=x0, y0=y0, x1=x1, y1=y1,
                     default.units=default.units,
                     arrow=arrow,
                     name=name, gp=gp, vp=vp)
  if (draw)
    grid.draw(sg)
  invisible(sg)
}

######################################
# ARROWS primitive
######################################

# Superceded by 'arrow' arg to line-drawing primitives
# which contains an "arrow" object
validDetails.arrows <- function(x) {
  if ((!is.null(x$x) && !is.unit(x$x)) ||
      (!is.null(x$y) && !is.unit(x$y)))
    stop("'x' and 'y' must be units or NULL")
  if (!is.unit(x$length))
    stop("'length' must be a 'unit' object")
  x$ends <- as.integer(match(x$ends, c("first", "last", "both")))
  x$type <- as.integer(match(x$type, c("open", "closed")))
  if (any(is.na(x$ends)) || any(is.na(x$type)))
    stop("invalid 'ends' or 'type' argument")
  x
}

drawDetails.arrows <- function(x, recording=TRUE) {
  if (is.null(x$x)) { # y should be null too
    if (!is.null(x$y))
      stop("corrupt 'arrows' object")
    lineThing <- getGrob(x, childNames(x))
    # This could be done via method dispatch, but that really
    # seemed like overkill
    # OTOH, this is NOT user-extensible
    # AND the code for, e.g., "lines" is not located with
    # the other grid.lines code so changes there are unlikely
    # to propagate to here (e.g., add an id arg to grid.lines?
    if (inherits(lineThing, "line.to")) {
      x1 <- NULL
      x2 <- lineThing$x
      y1 <- NULL
      y2 <- lineThing$y
      xnm1 <- NULL
      xn <- lineThing$x
      ynm1 <- NULL
      yn <- lineThing$y
    } else if (inherits(lineThing, "lines")) {
      # x or y may be recycled
      n <- max(length(lineThing$x),
               length(lineThing$y))
      xx <- rep(lineThing$x, length.out=2)
      x1 <- xx[1L]
      x2 <- xx[2L]
      xx <- rep(lineThing$x, length.out=n)
      xnm1 <- xx[n - 1]
      xn <- xx[n]
      yy <- rep(lineThing$y, length.out=2)
      y1 <- yy[1L]
      y2 <- yy[2L]
      yy <- rep(lineThing$y, length.out=n)
      ynm1 <- yy[n - 1]
      yn <- yy[n]
    } else { # inherits(lineThing, "segments")
      x1 <- lineThing$x0
      x2 <- lineThing$x1
      xnm1 <- lineThing$x0
      xn <- lineThing$x1
      y1 <- lineThing$y0
      y2 <- lineThing$y1
      ynm1 <- lineThing$y0
      yn <- lineThing$y1
    }
  } else {
    # x or y may be recycled
    n <- max(length(x$x), length(x$y))
    xx <- rep(x$x, length.out=2)
    x1 <- xx[1L]
    x2 <- xx[2L]
    xx <- rep(x$x, length.out=n)
    xnm1 <- xx[n - 1]
    xn <- xx[n]
    yy <- rep(x$y, length.out=2)
    y1 <- yy[1L]
    y2 <- yy[2L]
    yy <- rep(x$y, length.out=n)
    ynm1 <- yy[n - 1]
    yn <- yy[n]
    grid.Call.graphics(L_lines, x$x, x$y,
                       list(as.integer(1L:n)),
                       NULL)
  }
  grid.Call.graphics(L_arrows, x1, x2, xnm1, xn, y1, y2, ynm1, yn,
                     x$angle, x$length, x$ends, x$type)
}

widthDetails.arrows <- function(x) {
  if (is.null(x$x)) { # y should be null too
    if (!is.null(x$y))
      stop("corrupt 'arrows' object")
    lineThing <- getGrob(x, childNames(x))
    widthDetails(lineThing)
  } else {
    bounds <- grid.Call(L_locnBounds, x$x, x$y, 0)
    if (is.null(bounds))
      unit(0, "inches")
    else
      unit(bounds[3L], "inches")
  }
}

heightDetails.arrows <- function(x) {
  if (is.null(x$x)) { # y should be null too
    if (!is.null(x$y))
      stop("corrupt 'arrows' object")
    lineThing <- getGrob(x, childNames(x))
    heightDetails(lineThing)
  } else {
    bounds <- grid.Call(L_locnBounds, x$x, x$y, 0)
    if (is.null(bounds))
      unit(0, "inches")
    else
      unit(bounds[4L], "inches")
  }
}

arrowsGrob <- function(x=c(0.25, 0.75), y=0.5,
                        default.units="npc",
                        grob=NULL,
                        angle=30, length=unit(0.25, "inches"),
                        ends="last", type="open",
                        name=NULL, gp=gpar(), vp=NULL) {
    warning("grid.arrows() has been deprecated;  use 'arrow' arguments to line drawing functions", domain = NA)
  if (is.null(grob)) {
    if (!is.unit(x))
      x <- unit(x, default.units)
    if (!is.unit(y))
      y <- unit(y, default.units)
  }
  # Check the grob here
  # Not in validDetails.arrows because that is for checking
  # slots of an arrows object (the grob is a child of the arrows object)
  # A possible alternative design would have a copy of the grob
  # stored in a slot of the arrows object;  then it could be checked
  # in the validDetails AND it could be edited
  if (!is.null(grob)) {
    # The grob can only be a "lines" or "segments"
    # (splines would be another candidate if they existed)
    if (!(inherits(grob, "lines") ||
          inherits(grob, "segments") ||
          inherits(grob, "line.to")))
      stop("the 'grob' argument must be a 'line.to', 'lines', or 'segments' grob")
    x <- y <- NULL
  }
  gTree(x=x, y=y, children=if (is.null(grob)) NULL else gList(grob),
       angle=as.numeric(angle), length=length,
       ends=ends, type=type,
       name=name, gp=gp, vp=vp, cl="arrows")
}

grid.arrows <- function(x=c(0.25, 0.75), y=0.5,
                        default.units="npc",
                        grob=NULL,
                        angle=30, length=unit(0.25, "inches"),
                        ends="last", type="open",
                        name=NULL, gp=gpar(), draw=TRUE, vp=NULL) {
  ag <- arrowsGrob(x=x, y=y,
                   default.units=default.units,
                   grob=grob, angle=angle, length=length,
                   ends=ends, type=type,
                   name=name, gp=gp, vp=vp)
  if (draw)
    grid.draw(ag)
  invisible(ag)
}

######################################
# POLYGON primitive
######################################

validDetails.polygon <- function(x) {
  if (!is.unit(x$x) ||
      !is.unit(x$y))
    stop("'x' and 'y' must be units")
  if (!is.null(x$id) && !is.null(x$id.lengths))
    stop("it is invalid to specify both 'id' and 'id.lengths'")
  if (length(x$x) != length(x$y))
    stop("'x' and 'y' must be same length")
  if (!is.null(x$id) && (length(x$id) != length(x$x)))
    stop("'x' and 'y' and 'id' must all be same length")
  if (!is.null(x$id))
    x$id <- as.integer(x$id)
  if (!is.null(x$id.lengths) && (sum(x$id.lengths) != length(x$x)))
    stop("'x' and 'y' and 'id.lengths' must specify same overall length")
  if (!is.null(x$id.lengths))
    x$id.lengths <- as.integer(x$id.lengths)
  x
}

drawDetails.polygon <- function(x, recording=TRUE) {
  if (is.null(x$id) && is.null(x$id.lengths))
    grid.Call.graphics(L_polygon, x$x, x$y,
                       list(as.integer(seq_along(x$x))))
  else {
    if (is.null(x$id)) {
      n <- length(x$id.lengths)
      id <- rep(1L:n, x$id.lengths)
    } else {
      n <- length(unique(x$id))
      id <- x$id
    }
    index <- split(as.integer(seq_along(x$x)), id)
    grid.Call.graphics(L_polygon, x$x, x$y, index)
  }
}

xDetails.polygon <- function(x, theta) {
    bounds <- grid.Call(L_locnBounds, x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[1L], "inches")
}

yDetails.polygon <- function(x, theta) {
    bounds <- grid.Call(L_locnBounds, x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[2L], "inches")
}

widthDetails.polygon <- function(x) {
  bounds <- grid.Call(L_locnBounds, x$x, x$y, 0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[3L], "inches")
}

heightDetails.polygon <- function(x) {
  bounds <- grid.Call(L_locnBounds, x$x, x$y, 0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[4L], "inches")
}

polygonGrob <- function(x=c(0, 0.5, 1, 0.5), y=c(0.5, 1, 0.5, 0),
                        id=NULL, id.lengths=NULL,
                        default.units="npc",
                        name=NULL, gp=gpar(), vp=NULL) {
  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)
  grob(x=x, y=y, id=id,
       id.lengths=id.lengths,
       name=name, gp=gp, vp=vp, cl="polygon")
}

grid.polygon <- function(x=c(0, 0.5, 1, 0.5), y=c(0.5, 1, 0.5, 0),
                         id=NULL, id.lengths=NULL,
                         default.units="npc",
                         name=NULL, gp=gpar(), draw=TRUE, vp=NULL) {
  pg <- polygonGrob(x=x, y=y, id=id, id.lengths=id.lengths,
                    default.units=default.units,
                    name=name, gp=gp, vp=vp)
  if (draw)
    grid.draw(pg)
  invisible(pg)
}

######################################
# PATH primitive
######################################

validDetails.pathgrob <- function(x) {
    if (!is.unit(x$x) || !is.unit(x$y))
        stop("'x' and 'y' must be units")
    if (!is.null(x$id) && !is.null(x$id.lengths))
        stop("it is invalid to specify both 'id' and 'id.lengths'")
    if (length(x$x) != length(x$y))
        stop("'x' and 'y' must be same length")
    if (!is.null(x$id) && (length(x$id) != length(x$x)))
        stop("'x' and 'y' and 'id' must all be same length")
    if (!is.null(x$id))
        x$id <- as.integer(x$id)
    if (!is.null(x$id.lengths) && (sum(x$id.lengths) != length(x$x)))
        stop("'x' and 'y' and 'id.lengths' must specify same overall length")
    if (!is.null(x$id.lengths))
        x$id.lengths <- as.integer(x$id.lengths)
    x
}

xDetails.pathgrob <- function(x, theta) {
    bounds <- grid.Call(L_locnBounds, x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[1L], "inches")
}

yDetails.pathgrob <- function(x, theta) {
    bounds <- grid.Call(L_locnBounds, x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[2L], "inches")
}

widthDetails.pathgrob <- function(x) {
    bounds <- grid.Call(L_locnBounds, x$x, x$y, 0)
    if (is.null(bounds))
        unit(0, "inches")
    else
        unit(bounds[3L], "inches")
}

heightDetails.pathgrob <- function(x) {
    bounds <- grid.Call(L_locnBounds, x$x, x$y, 0)
    if (is.null(bounds))
        unit(0, "inches")
    else
        unit(bounds[4L], "inches")
}


drawDetails.pathgrob <- function(x, recording=TRUE) {
      if (is.null(x$id) && is.null(x$id.lengths))
          grid.Call.graphics(L_polygon, x$x, x$y,
                             list(as.integer(seq_along(x$x))))
  else {
    if (is.null(x$id)) {
      n <- length(x$id.lengths)
      id <- rep(1L:n, x$id.lengths)
    } else {
      n <- length(unique(x$id))
      id <- x$id
    }
    index <- split(as.integer(seq_along(x$x)), id)
    grid.Call.graphics(L_path, x$x, x$y, index,
                       switch(x$rule, winding=1L, evenodd=0L))
  }
}

pathGrob <- function(x, y,
                     id=NULL, id.lengths=NULL,
                     rule="winding",
                     default.units="npc",
                     name=NULL, gp=gpar(), vp=NULL) {
  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)
  grob(x=x, y=y, id=id,
       id.lengths=id.lengths,
       rule=rule,
       name=name, gp=gp, vp=vp, cl="pathgrob")
}

grid.path <- function(...) {
  grid.draw(pathGrob(...))
}

######################################
# XSPLINE primitive
######################################

validDetails.xspline <- function(x) {
  if (!is.unit(x$x) ||
      !is.unit(x$y))
    stop("x and y must be units")
  if (!is.null(x$id) && !is.null(x$id.lengths))
    stop("it is invalid to specify both 'id' and 'id.lengths'")
  nx <- length(x$x)
  ny <- length(x$y)
  if (nx != ny)
    stop("'x' and 'y' must be same length")
  if (!is.null(x$id) && (length(x$id) != nx))
    stop("'x' and 'y' and 'id' must all be same length")
  if (!is.null(x$id))
    x$id <- as.integer(x$id)
  if (!is.null(x$id.lengths) && (sum(x$id.lengths) != nx))
    stop("'x' and 'y' and 'id.lengths' must specify same overall length")
  if (!is.null(x$id.lengths))
    x$id.lengths <- as.integer(x$id.lengths)
  if (!(is.null(x$arrow) || inherits(x$arrow, "arrow")))
      stop("invalid 'arrow' argument")
  if (any(x$shape < -1 || x$shape > 1))
    stop("'shape' must be between -1 and 1")
  x$open <- as.logical(x$open)
  # Force all first and last shapes to be 0 for open xsplines
  if (x$open) {
      x$shape <- rep(x$shape, length.out=nx)
      # Watch out for id or id.length!
      index <- xsplineIndex(x)
      first <- sapply(index, min)
      last <- sapply(index, max)
      x$shape[c(first, last)] <- 0
  }
  x
}

xsplineIndex <- function(x) {
  if (is.null(x$id) && is.null(x$id.lengths))
      list(as.integer(seq_along(x$x)))
  else {
    if (is.null(x$id)) {
      n <- length(x$id.lengths)
      id <- rep(1L:n, x$id.lengths)
    } else {
      n <- length(unique(x$id))
      id <- x$id
    }
    split(as.integer(seq_along(x$x)), id)
  }
}

drawDetails.xspline <- function(x, recording=TRUE) {
    grid.Call.graphics(L_xspline, x$x, x$y, x$shape, x$open, x$arrow,
                       x$repEnds, xsplineIndex(x))
}

xDetails.xspline <- function(x, theta) {
  bounds <- grid.Call(L_xsplineBounds, x$x, x$y, x$shape, x$open, x$arrow,
                      x$repEnds, xsplineIndex(x), theta)
  if (is.null(bounds))
    unit(0.5, "npc")
  else
    unit(bounds[1L], "inches")
}

yDetails.xspline <- function(x, theta) {
  bounds <- grid.Call(L_xsplineBounds, x$x, x$y, x$shape, x$open, x$arrow,
                      x$repEnds, xsplineIndex(x), theta)
  if (is.null(bounds))
    unit(0.5, "npc")
  else
    unit(bounds[2L], "inches")
}

widthDetails.xspline <- function(x) {
  bounds <- grid.Call(L_xsplineBounds, x$x, x$y, x$shape, x$open, x$arrow,
                      x$repEnds, list(as.integer(seq_along(x$x))), 0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[3L], "inches")
}

heightDetails.xspline <- function(x) {
  bounds <- grid.Call(L_xsplineBounds, x$x, x$y, x$shape, x$open, x$arrow,
                      x$repEnds, list(as.integer(seq_along(x$x))), 0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[4L], "inches")
}

xsplineGrob <- function(x=c(0, 0.5, 1, 0.5), y=c(0.5, 1, 0.5, 0),
                        id=NULL, id.lengths=NULL,
                        default.units="npc",
                        shape=0, open=TRUE, arrow=NULL, repEnds=TRUE,
                        name=NULL, gp=gpar(), vp=NULL) {
  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)
  grob(x=x, y=y, shape=shape, open=open,
       id=id, id.lengths=id.lengths, arrow=arrow, repEnds=repEnds,
       name=name, gp=gp, vp=vp, cl="xspline")
}

grid.xspline <- function(...) {
  grid.draw(xsplineGrob(...))
}

xsplinePoints <- function(x) {
    # Mimic drawGrob() to ensure x$vp and x$gp enforced
    dlon <- grid.Call(L_setDLon, FALSE)
    on.exit(grid.Call(L_setDLon, dlon))
    tempgpar <- grid.Call(L_getGPar)
    on.exit(grid.Call(L_setGPar, tempgpar), add=TRUE)
    preDraw(x)
    # Raw pts in dev coords
    devPoints <- grid.Call(L_xsplinePoints,
                           x$x, x$y, x$shape, x$open, x$arrow,
                           x$repEnds, xsplineIndex(x), 0)
    postDraw(x)
    # Convert to units in inches
    unitPoints <- lapply(devPoints,
                         function(x) {
                             names(x) <- c("x", "y")
                             x$x <- unit(x$x, "inches")
                             x$y <- unit(x$y, "inches")
                             x
                         })
    if (length(unitPoints) == 1)
        unitPoints <- unitPoints[[1]]
    unitPoints
}

######################################
# BEZIER primitive
######################################

# A bezier grob that works of a (not-100% accurate) approximation
# using X-splines

# X-Spline approx to Bezier
Ms <- 1/6*rbind(c(1, 4, 1, 0),
                c(-3, 0, 3, 0),
                c(3, -6, 3, 0),
                c(-1, 3, -3, 1))
Msinv <- solve(Ms)
# Bezier control matrix
Mb <- rbind(c(1, 0, 0, 0),
            c(-3, 3, 0, 0),
            c(3, -6, 3, 0),
            c(-1, 3, -3, 1))

splinePoints <- function(xb, yb, idIndex) {
    xs <- unlist(lapply(idIndex,
                        function(i) {
                            Msinv %*% Mb %*% xb[i]
                        }))
    ys <- unlist(lapply(idIndex,
                        function(i) {
                            Msinv %*% Mb %*% yb[i]
                        }))
    list(x=xs, y=ys)
}

splinegrob <- function(x) {
    xx <- convertX(x$x, "inches", valueOnly=TRUE)
    yy <- convertY(x$y, "inches", valueOnly=TRUE)
    sp <- splinePoints(xx, yy, xsplineIndex(x))
    xsplineGrob(sp$x, sp$y, default.units="inches",
                id=x$id, id.lengths=x$id.lengths,
                shape=1, repEnds=FALSE,
                arrow=x$arrow, name=x$name,
                gp=x$gp, vp=x$vp)
}

validDetails.beziergrob <- function(x) {
    if (!is.unit(x$x) ||
        !is.unit(x$y))
        stop("x and y must be units")
    if (!is.null(x$id) && !is.null(x$id.lengths))
        stop("it is invalid to specify both 'id' and 'id.lengths'")
    nx <- length(x$x)
    ny <- length(x$y)
    if (nx != ny)
        stop("'x' and 'y' must be same length")
    if (!is.null(x$id) && (length(x$id) != nx))
        stop("'x' and 'y' and 'id' must all be same length")
    if (!is.null(x$id))
        x$id <- as.integer(x$id)
    if (!is.null(x$id.lengths) && (sum(x$id.lengths) != nx))
        stop("'x' and 'y' and 'id.lengths' must specify same overall length")
    if (!is.null(x$id.lengths))
        x$id.lengths <- as.integer(x$id.lengths)
    if (is.null(x$id) && is.null(x$id.lengths)) {
        if (length(x$x) != 4)
            stop("must have exactly 4 control points")
    } else {
        if (is.null(x$id)) {
            id <- rep(1L:n, x$id.lengths)
        } else {
            id <- x$id
        }
        xper <- split(x$x, id)
        if (any(sapply(xper, length) != 4))
            stop("must have exactly 4 control points per Bezier curve")
    }
    if (!(is.null(x$arrow) || inherits(x$arrow, "arrow")))
        stop("invalid 'arrow' argument")
    x
}

makeContent.beziergrob <- function(x) {
    splinegrob(x)
}

xDetails.beziergrob <- function(x, theta) {
    xDetails(splinegrob(x), theta)
}

yDetails.beziergrob <- function(x, theta) {
    yDetails(splinegrob(x), theta)
}

widthDetails.beziergrob <- function(x) {
    widthDetails(splinegrob(x))
}

heightDetails.beziergrob <- function(x) {
    heightDetails(splinegrob(x))
}

bezierGrob <- function(x=c(0, 0.5, 1, 0.5), y=c(0.5, 1, 0.5, 0),
                       id=NULL, id.lengths=NULL,
                       default.units="npc", arrow=NULL,
                       name=NULL, gp=gpar(), vp=NULL) {
    if (!is.unit(x))
        x <- unit(x, default.units)
    if (!is.unit(y))
        y <- unit(y, default.units)
    grob(x=x, y=y,
         id=id, id.lengths=id.lengths, arrow=arrow,
         name=name, gp=gp, vp=vp, cl="beziergrob")
}

grid.bezier <- function(...) {
    grid.draw(bezierGrob(...))
}

bezierPoints <- function(x) {
    sg <- splinegrob(x)
    # splinegrob() does not make use of x$vp
    sg$vp <- x$vp
    xsplinePoints(sg)
}


######################################
# CIRCLE primitive
######################################

validDetails.circle <- function(x) {
  if (!is.unit(x$x) ||
      !is.unit(x$y) ||
      !is.unit(x$r))
    stop("'x', 'y', and 'r' must be units")
  x
}

drawDetails.circle <- function(x, recording=TRUE) {
  grid.Call.graphics(L_circle, x$x, x$y, x$r)
}

xDetails.circle <- function(x, theta) {
  bounds <- grid.Call(L_circleBounds, x$x, x$y, x$r, theta)
  if (is.null(bounds))
    unit(0.5, "npc")
  else
    unit(bounds[1L], "inches")
}

yDetails.circle <- function(x, theta) {
  bounds <- grid.Call(L_circleBounds, x$x, x$y, x$r, theta)
  if (is.null(bounds))
    unit(0.5, "npc")
  else
    unit(bounds[2L], "inches")
}

widthDetails.circle <- function(x) {
  bounds <- grid.Call(L_circleBounds, x$x, x$y, x$r, 0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[3L], "inches")
}

heightDetails.circle <- function(x) {
  bounds <- grid.Call(L_circleBounds, x$x, x$y, x$r, 0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[4L], "inches")
}

circleGrob <- function(x=0.5, y=0.5, r=0.5,
                       default.units="npc",
                       name=NULL, gp=gpar(), vp=NULL) {
  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)
  if (!is.unit(r))
    r <- unit(r, default.units)
  grob(x=x, y=y, r=r, name=name, gp=gp, vp=vp, cl="circle")
}

grid.circle <- function(x=0.5, y=0.5, r=0.5,
                        default.units="npc",
                        name=NULL, gp=gpar(), draw=TRUE, vp=NULL) {
  cg <- circleGrob(x=x, y=y, r=r,
                   default.units=default.units,
                   name=name, gp=gp, vp=vp)
  if (draw)
    grid.draw(cg)
  invisible(cg)
}

######################################
# RECT primitive
######################################
validDetails.rect <- function(x) {
  if (!is.unit(x$x) ||
      !is.unit(x$y) ||
      !is.unit(x$width) ||
      !is.unit(x$height))
    stop("'x', 'y', 'width', and 'height' must be units")
  valid.just(x$just)
  if (!is.null(x$hjust))
    x$hjust <- as.numeric(x$hjust)
  if (!is.null(x$vjust))
    x$vjust <- as.numeric(x$vjust)
  x
}

drawDetails.rect <- function(x, recording=TRUE) {
  grid.Call.graphics(L_rect, x$x, x$y, x$width, x$height,
                     resolveHJust(x$just, x$hjust),
                     resolveVJust(x$just, x$vjust))
}

xDetails.rect <- function(x, theta) {
  bounds <- grid.Call(L_rectBounds, x$x, x$y, x$width, x$height,
                      resolveHJust(x$just, x$hjust),
                      resolveVJust(x$just, x$vjust),
                      theta)
  if (is.null(bounds))
    unit(0.5, "npc")
  else
    unit(bounds[1L], "inches")
}

yDetails.rect <- function(x, theta) {
  bounds <- grid.Call(L_rectBounds, x$x, x$y, x$width, x$height,
                      resolveHJust(x$just, x$hjust),
                      resolveVJust(x$just, x$vjust),
                      theta)
  if (is.null(bounds))
    unit(0.5, "npc")
  else
    unit(bounds[2L], "inches")
}

widthDetails.rect <- function(x) {
  bounds <- grid.Call(L_rectBounds, x$x, x$y, x$width, x$height,
                      resolveHJust(x$just, x$hjust),
                      resolveVJust(x$just, x$vjust),
                      0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[3L], "inches")
}

heightDetails.rect <- function(x) {
  bounds <- grid.Call(L_rectBounds, x$x, x$y, x$width, x$height,
                      resolveHJust(x$just, x$hjust),
                      resolveVJust(x$just, x$vjust),
                      0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[4L], "inches")
}

rectGrob <- function(x=unit(0.5, "npc"), y=unit(0.5, "npc"),
                     width=unit(1, "npc"), height=unit(1, "npc"),
                     just="centre", hjust=NULL, vjust=NULL,
                     default.units="npc",
                     name=NULL, gp=gpar(), vp=NULL) {
  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)
  if (!is.unit(width))
    width <- unit(width, default.units)
  if (!is.unit(height))
    height <- unit(height, default.units)
  grob(x=x, y=y, width=width, height=height, just=just,
       hjust=hjust, vjust=vjust,
       name=name, gp=gp, vp=vp, cl="rect")
}

grid.rect <- function(x=unit(0.5, "npc"), y=unit(0.5, "npc"),
                      width=unit(1, "npc"), height=unit(1, "npc"),
                      just="centre", hjust=NULL, vjust=NULL,
                      default.units="npc",
                      name=NULL, gp=gpar(), draw=TRUE, vp=NULL) {
  rg <- rectGrob(x=x, y=y, width=width, height=height, just=just,
                 hjust=hjust, vjust=vjust,
                 default.units=default.units,
                 name=name, gp=gp, vp=vp)
  if (draw)
    grid.draw(rg)
  invisible(rg)
}

######################################
# RASTER primitive
######################################

validDetails.rastergrob <- function(x) {
    if (!(is.raster(x$raster) || inherits(x$raster, "nativeRaster")))
        x$raster <- as.raster(x$raster)
    if (!is.unit(x$x) ||
        !is.unit(x$y) ||
        (!is.null(x$width) && !is.unit(x$width)) ||
        (!is.null(x$height) && !is.unit(x$height)))
        stop("'x', 'y', 'width', and 'height' must be units")
    valid.just(x$just)
    if (!is.null(x$hjust))
        x$hjust <- as.numeric(x$hjust)
    if (!is.null(x$vjust))
        x$vjust <- as.numeric(x$vjust)
    x
}

resolveRasterSize <- function(x) {
    if (is.null(x$width)) {
        if (is.null(x$height)) {
            rasterRatio <- dim(x$raster)[1]/dim(x$raster)[2]
            vpWidth <- convertWidth(unit(1, "npc"), "inches", valueOnly=TRUE)
            vpHeight <- convertHeight(unit(1, "npc"), "inches", valueOnly=TRUE)
            vpRatio <- vpHeight/vpWidth
            if (rasterRatio > vpRatio) {
                x$height <- unit(vpHeight, "inches")
                x$width <- unit(vpHeight*dim(x$raster)[2]/dim(x$raster)[1],
                                "inches")
            } else {
                x$width <- unit(vpWidth, "inches")
                x$height <- unit(vpWidth*dim(x$raster)[1]/dim(x$raster)[2],
                                 "inches")
            }
        } else {
            h <- convertHeight(x$height, "inches", valueOnly=TRUE)
            x$width <- unit(h*dim(x$raster)[2]/dim(x$raster)[1],
                            "inches")
        }
    } else {
        if (is.null(x$height)) {
            w <- convertWidth(x$width, "inches", valueOnly=TRUE)
            x$height <- unit(w*dim(x$raster)[1]/dim(x$raster)[2],
                             "inches")
        }
    }
    x
}

drawDetails.rastergrob <- function(x, recording=TRUE) {
    # At this point resolve NULL width/height based on
    # image dimensions
    x <- resolveRasterSize(x)
    if (is.null(x$width)) {
        if (is.null(x$height)) {
            rasterRatio <- dim(x$raster)[1]/dim(x$raster)[2]
            vpWidth <- convertWidth(unit(1, "npc"), "inches", valueOnly=TRUE)
            vpHeight <- convertHeight(unit(1, "npc"), "inches", valueOnly=TRUE)
            vpRatio <- vpHeight/vpWidth
            if (rasterRatio > vpRatio) {
                x$height <- unit(vpHeight, "inches")
                x$width <- unit(vpHeight*dim(x$raster)[2]/dim(x$raster)[1],
                                "inches")
            } else {
                x$width <- unit(vpWidth, "inches")
                x$height <- unit(vpWidth*dim(x$raster)[1]/dim(x$raster)[2],
                                 "inches")
            }
        } else {
            h <- convertHeight(x$height, "inches", valueOnly=TRUE)
            x$width <- unit(h*dim(x$raster)[2]/dim(x$raster)[1],
                            "inches")
        }
    } else {
        if (is.null(x$height)) {
            w <- convertWidth(x$width, "inches", valueOnly=TRUE)
            x$height <- unit(w*dim(x$raster)[1]/dim(x$raster)[2],
                             "inches")
        }
    }
    grid.Call.graphics(L_raster, x$raster,
                       x$x, x$y, x$width, x$height,
                       resolveHJust(x$just, x$hjust),
                       resolveVJust(x$just, x$vjust),
                       x$interpolate)
}

xDetails.rastergrob <- function(x, theta) {
    x <- resolveRasterSize(x)
    bounds <- grid.Call(L_rectBounds, x$x, x$y, x$width, x$height,
                        resolveHJust(x$just, x$hjust),
                        resolveVJust(x$just, x$vjust),
                        theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[1L], "inches")
}

yDetails.rastergrob <- function(x, theta) {
    x <- resolveRasterSize(x)
    bounds <- grid.Call(L_rectBounds, x$x, x$y, x$width, x$height,
                        resolveHJust(x$just, x$hjust),
                        resolveVJust(x$just, x$vjust),
                        theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[2L], "inches")
}

widthDetails.rastergrob <- function(x) {
    x <- resolveRasterSize(x)
    bounds <- grid.Call(L_rectBounds, x$x, x$y, x$width, x$height,
                        resolveHJust(x$just, x$hjust),
                        resolveVJust(x$just, x$vjust),
                        0)
    if (is.null(bounds))
        unit(0, "inches")
    else
        unit(bounds[3L], "inches")
}

heightDetails.rastergrob <- function(x) {
    x <- resolveRasterSize(x)
    bounds <- grid.Call(L_rectBounds, x$x, x$y, x$width, x$height,
                        resolveHJust(x$just, x$hjust),
                        resolveVJust(x$just, x$vjust),
                        0)
    if (is.null(bounds))
        unit(0, "inches")
    else
        unit(bounds[4L], "inches")
}

rasterGrob <- function(image,
                       x=unit(0.5, "npc"), y=unit(0.5, "npc"),
                       width=NULL, height=NULL,
                       just="centre", hjust=NULL, vjust=NULL,
                       interpolate=TRUE,
                       default.units="npc",
                       name=NULL, gp=gpar(), vp=NULL) {

    if (inherits(image, "nativeRaster"))
        raster <- image
    else
        raster <- as.raster(image)
    if (!is.unit(x))
        x <- unit(x, default.units)
    if (!is.unit(y))
        y <- unit(y, default.units)
    if (!is.null(width) && !is.unit(width))
        width <- unit(width, default.units)
    if (!is.null(height) && !is.unit(height))
        height <- unit(height, default.units)
    grob(raster=raster, x=x, y=y, width=width, height=height, just=just,
         hjust=hjust, vjust=vjust, interpolate=interpolate,
         name=name, gp=gp, vp=vp, cl="rastergrob")
}

grid.raster <- function(image,
                        x=unit(0.5, "npc"), y=unit(0.5, "npc"),
                        width=NULL, height=NULL,
                        just="centre", hjust=NULL, vjust=NULL,
                        interpolate=TRUE,
                        default.units="npc",
                        name=NULL, gp=gpar(), vp=NULL) {
    rg <- rasterGrob(image,
                     x=x, y=y, width=width, height=height, just=just,
                     hjust=hjust, vjust=vjust, interpolate=interpolate,
                     default.units=default.units,
                     name=name, gp=gp, vp=vp)
    grid.draw(rg)
}

######################################
# TEXT primitive
######################################
validDetails.text <- function(x) {
  if (!is.language(x$label))
    x$label <- as.character(x$label)
  if (!is.unit(x$x) ||
      !is.unit(x$y))
    stop("'x' and 'y' must be units")
  x$rot <- as.numeric(x$rot)
  if (!all(is.finite(x$rot)) || length(x$rot) == 0)
    stop("invalid 'rot' value")
  valid.just(x$just)
  if (!is.null(x$hjust))
    x$hjust <- as.numeric(x$hjust)
  if (!is.null(x$vjust))
    x$vjust <- as.numeric(x$vjust)
  x$check.overlap <- as.logical(x$check.overlap)
  x
}

drawDetails.text <- function(x, recording=TRUE) {
  grid.Call.graphics(L_text, as.graphicsAnnot(x$label),
                     x$x, x$y,
                     resolveHJust(x$just, x$hjust),
                     resolveVJust(x$just, x$vjust),
                     x$rot, x$check.overlap)
}

xDetails.text <- function(x, theta) {
  bounds <- grid.Call(L_textBounds, as.graphicsAnnot(x$label),
                      x$x, x$y,
                      resolveHJust(x$just, x$hjust),
                      resolveVJust(x$just, x$vjust),
                      x$rot, theta)
  if (is.null(bounds))
    unit(0.5, "npc")
  else
    unit(bounds[1L], "inches")
}

yDetails.text <- function(x, theta) {
  bounds <- grid.Call(L_textBounds, as.graphicsAnnot(x$label),
                      x$x, x$y,
                      resolveHJust(x$just, x$hjust),
                      resolveVJust(x$just, x$vjust),
                      x$rot, theta)
  if (is.null(bounds))
    unit(0.5, "npc")
  else
    unit(bounds[2L], "inches")
}

widthDetails.text <- function(x) {
  bounds <- grid.Call(L_textBounds, as.graphicsAnnot(x$label),
                      x$x, x$y,
                      resolveHJust(x$just, x$hjust),
                      resolveVJust(x$just, x$vjust),
                      x$rot, 0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[3L], "inches")
}

heightDetails.text <- function(x) {
  bounds <- grid.Call(L_textBounds, as.graphicsAnnot(x$label),
                      x$x, x$y,
                      resolveHJust(x$just, x$hjust),
                      resolveVJust(x$just, x$vjust),
                      x$rot, 0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[4L], "inches")
}

ascentDetails.text <- function(x) {
    if (length(x$label) == 1) {
        metrics <- grid.Call(L_stringMetric, as.graphicsAnnot(x$label))
        unit(metrics[[1]], "inches")
    } else {
        heightDetails(x)
    }
}

descentDetails.text <- function(x) {
    if (length(x$label) == 1) {
        metrics <- grid.Call(L_stringMetric, as.graphicsAnnot(x$label))
        unit(metrics[[2]], "inches")
    } else {
        unit(0, "inches")
    }
}

textGrob <- function(label, x=unit(0.5, "npc"), y=unit(0.5, "npc"),
                     just="centre", hjust=NULL, vjust=NULL,
                     rot=0, check.overlap=FALSE,
                     default.units="npc",
                     name=NULL, gp=gpar(), vp=NULL) {
  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)
  grob(label=label, x=x, y=y, just=just, hjust=hjust, vjust=vjust,
       rot=rot, check.overlap=check.overlap,
       name=name, gp=gp, vp=vp, cl="text")
}

grid.text <- function(label, x=unit(0.5, "npc"), y=unit(0.5, "npc"),
                      just="centre", hjust=NULL, vjust=NULL,
                      rot=0, check.overlap=FALSE,
                      default.units="npc",
                      name=NULL, gp=gpar(), draw=TRUE, vp=NULL) {
  tg <- textGrob(label=label, x=x, y=y, just=just,
                 hjust=hjust, vjust=vjust, rot=rot,
                 check.overlap=check.overlap,
                 default.units=default.units,
                 name=name, gp=gp, vp=vp)
  if (draw)
    grid.draw(tg)
  invisible(tg)
}

######################################
# POINTS primitive
######################################
valid.pch <- function(pch) {
  if (length(pch) == 0L)
    stop("zero-length 'pch'")
  if (is.null(pch))
    pch <- 1L
  else if (!is.character(pch))
    pch <- as.integer(pch)
  pch
}

validDetails.points <- function(x) {
  if (!is.unit(x$x) ||
      !is.unit(x$y) ||
      !is.unit(x$size))
    stop("'x', 'y' and 'size' must be units")
  if (length(x$x) != length(x$y))
    stop("'x' and 'y' must be 'unit' objects and have the same length")
  x$pch <- valid.pch(x$pch)
  x
}

drawDetails.points <- function(x, recording=TRUE) {
  grid.Call.graphics(L_points, x$x, x$y, x$pch, x$size)
}

# FIXME:  does not take into account the size of the symbols
xDetails.points <- function(x, theta) {
    bounds <- grid.Call(L_locnBounds, x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[1L], "inches")
}

yDetails.points <- function(x, theta) {
    bounds <- grid.Call(L_locnBounds, x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[2L], "inches")
}

widthDetails.points <- function(x) {
  bounds <- grid.Call(L_locnBounds, x$x, x$y, 0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[3L], "inches")
}

heightDetails.points <- function(x) {
  bounds <- grid.Call(L_locnBounds, x$x, x$y, 0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[4L], "inches")
}

pointsGrob <- function(x=stats::runif(10),
                       y=stats::runif(10),
                       pch=1, size=unit(1, "char"),
                       default.units="native",
                       name=NULL, gp=gpar(), vp=NULL) {
  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)
  grob(x=x, y=y, pch=pch, size=size,
       name=name, gp=gp, vp=vp, cl="points")
}

grid.points <- function(x=stats::runif(10),
                        y=stats::runif(10),
                        pch=1, size=unit(1, "char"),
                        default.units="native",
                        name=NULL, gp=gpar(),
                        draw=TRUE, vp=NULL) {
  pg <- pointsGrob(x=x, y=y, pch=pch, size=size,
                   default.units=default.units,
                   name=name, gp=gp, vp=vp)
  if (draw)
    grid.draw(pg)
  invisible(pg)
}

######################################
# CLIP primitive
######################################
validDetails.clip <- function(x) {
  if (!is.unit(x$x) ||
      !is.unit(x$y) ||
      !is.unit(x$width) ||
      !is.unit(x$height))
    stop("'x', 'y', 'width', and 'height' must be units")
  if (length(x$x) > 1 || length(x$y) > 1 ||
      length(x$width) > 1 || length(x$height) > 1)
    stop("'x', 'y', 'width', and 'height' must all be units of length 1")
  valid.just(x$just)
  if (!is.null(x$hjust))
    x$hjust <- as.numeric(x$hjust)
  if (!is.null(x$vjust))
    x$vjust <- as.numeric(x$vjust)
  x
}

drawDetails.clip <- function(x, recording=TRUE) {
  grid.Call.graphics(L_clip, x$x, x$y, x$width, x$height,
                     resolveHJust(x$just, x$hjust),
                     resolveVJust(x$just, x$vjust))
}

clipGrob <- function(x=unit(0.5, "npc"), y=unit(0.5, "npc"),
                     width=unit(1, "npc"), height=unit(1, "npc"),
                     just="centre", hjust=NULL, vjust=NULL,
                     default.units="npc",
                     name=NULL, vp=NULL) {
  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)
  if (!is.unit(width))
    width <- unit(width, default.units)
  if (!is.unit(height))
    height <- unit(height, default.units)
  grob(x=x, y=y, width=width, height=height, just=just,
       hjust=hjust, vjust=vjust,
       name=name, vp=vp, cl="clip")
}

grid.clip <- function(...) {
  grid.draw(clipGrob(...))
}


######################################
# NULL primitive
######################################

validDetails.null <- function(x) {
  if (!is.unit(x$x) ||
      !is.unit(x$y))
    stop("'x' and 'y' must be units")
  if (length(x$x) > 1 || length(x$y) > 1)
    stop("'x' and 'y' must all be units of length 1")
  x
}

drawDetails.null <- function(x, recording=TRUE) {
    # Deliberate null op.
    # NOTE: nothing will go on the graphics engine DL
    # This is ok I think because these grobs are only
    # useful on the grid DL (for other grid code to query
    # their size or location).
}

xDetails.null <- function(x, theta) {
    bounds <- grid.Call(L_locnBounds, x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[1L], "inches")
}

yDetails.null <- function(x, theta) {
    bounds <- grid.Call( L_locnBounds, x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[2L], "inches")
}

# Deliberately ZERO
widthDetails.null <- function(x) {
    unit(0, "inches")
}

heightDetails.null <- function(x) {
    unit(0, "inches")
}

# A grob with GUARANTEED zero-width
# also GUARANTEED NOT to draw anything
nullGrob <- function(x=unit(0.5, "npc"), y=unit(0.5, "npc"),
                     default.units="npc",
                     name=NULL, vp=NULL) {
    if (!is.unit(x))
        x <- unit(x, default.units)
    if (!is.unit(y))
        y <- unit(y, default.units)
    grob(x=x, y=y, name=name, vp=vp, cl="null")
}

# Convenient way to get nullGrob on the grid display list
grid.null <- function(...) {
    grid.draw(nullGrob(...))
}










#  File src/library/grid/R/roundrect.R
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


# Good idea to choose r as absolute unit or "snpc"
roundrectGrob <- function(x=0.5, y=0.5, width=1, height=1,
                          default.units="npc",
                          r=unit(0.1, "snpc"),
                          just="centre",
                          name=NULL, gp=NULL, vp=NULL) {
    if (!is.unit(x))
        x <- unit(x, default.units)
    if (!is.unit(y))
        y <- unit(y, default.units)
    if (!is.unit(width))
        width <- unit(width, default.units)
    if (!is.unit(height))
        height <- unit(height, default.units)
    grob(x=x, y=y, width=width, height=height, r=r, just=just,
         name=name, gp=gp, vp=vp, cl="roundrect")
}

grid.roundrect <- function(...) {
  grid.draw(roundrectGrob(...))
}

validDetails.roundrect <- function(x) {
    if (!is.unit(x$x) ||
        !is.unit(x$y) ||
        !is.unit(x$width) ||
        !is.unit(x$height))
        stop("'x', 'y', 'width', and 'height' must be units")
    if (!is.unit(x$r))
        stop("'r' must be a 'unit' object")
    valid.just(x$just)
    # Make sure that x and y are of length 1
    if (length(x$x) != 1 | length(x$y) != 1 |
        length(x$width) != 1 | length(x$height) != 1)
        stop("'x', 'y', 'width', and 'height' must have length 1")
    x
}

makeContext.roundrect <- function(x) {
    rrvp <- viewport(x$x, x$y, x$width, x$height, just=x$just,
                     name="rrvp")
    if (!is.null(x$vp)) {
        x$vp <- vpStack(x$vp, rrvp)
    } else {
        x$vp <- rrvp
    }
    x
}

# x, y, is the real corner
roundCorner <- function(num, x, y, r) {
  n <- 10*4
  t <- seq(0, 2*pi, length.out=n)
  cost <- cos(t)
  sint <- sin(t)
  if (num == 1) {
    xc <- x + r
    yc <- y + r
    subset <- (n/2):(3*n/4)
  } else if (num == 2) {
    xc <- x + r
    yc <- y - r
    subset <- (n/4):(n/2)
  } else if (num == 3) {
    xc <- x - r
    yc <- y - r
    subset <- 1L:(n/4)
  } else if (num == 4) {
    xc <- x - r
    yc <- y + r
    subset <- (3*n/4):n
  }
  list(x=xc + (cost*r)[subset], y=yc + (sint*r)[subset])
}

rrpoints <- function(x) {
  left <- 0
  bottom <- 0
  right <- convertX(unit(1, "npc"), "inches", valueOnly=TRUE)
  top <- convertY(unit(1, "npc"), "inches", valueOnly=TRUE)
  r <- min(convertWidth(x$r, "inches", valueOnly=TRUE),
           convertHeight(x$r, "inches", valueOnly=TRUE))
  corner1 <- roundCorner(1, left, bottom, r)
  corner2 <- roundCorner(2, left, top, r)
  corner3 <- roundCorner(3, right, top, r)
  corner4 <- roundCorner(4, right, bottom, r)
  xx <- unit(c(left + r, right - r, corner4$x,
               right, right, corner3$x,
               right - r, left + r, corner2$x,
               left, left, corner1$x),
             "inches")
  yy <- unit(c(bottom, bottom, corner4$y,
               bottom + r, top - r, corner3$y,
               top, top, corner2$y,
               top - r, bottom + r, corner1$y),
             "inches")
  list(x=xx, y=yy)
}

makeContent.roundrect <- function(x) {
    boundary <- rrpoints(x)
    polygonGrob(boundary$x, boundary$y,
                name=x$name, gp=x$gp, vp=x$vp)
}

xDetails.roundrect <- function(x, theta) {
    boundary <- rrpoints(x)
    bounds <- grid.Call(L_locnBounds, boundary$x, boundary$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[1L], "inches")
}

yDetails.roundrect <- function(x, theta) {
    boundary <- rrpoints(x)
    bounds <- grid.Call(L_locnBounds, boundary$x, boundary$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[2L], "inches")
}

widthDetails.roundrect <- function(x) {
    boundary <- rrpoints(x)
    bounds <- grid.Call(L_locnBounds, boundary$x, boundary$y, 0)
    if (is.null(bounds))
        unit(0, "inches")
    else
        unit(bounds[3L], "inches")
}

heightDetails.roundrect <- function(x) {
    boundary <- rrpoints(x)
    bounds <- grid.Call(L_locnBounds, boundary$x, boundary$y, 0)
    if (is.null(bounds))
        unit(0, "inches")
    else
        unit(bounds[4L], "inches")
}



#  File src/library/grid/R/size.R
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

# These functions are used to evaluate "grobwidth" and
# "grobheight" units.
# They are usually called from within the C code
# (specifically, from within unit.c)
# It should be noted that they only give the width/height
# of the grob in the current drawing context
# (i.e., evaluating the width/height in another context
#  will not necessarily give the same result)

# The C code to evaluate "grobwidth" and "grobheight" calls
# the preDrawDetails and postDrawDetails generics before and
# after the call to width/height() to allow for complex grobs which
# construct their own viewports.

#########
# X locations on edge
#########

xDetails <- function(x, theta) {
  UseMethod("xDetails")
}

xDetails.default <- function(x, theta) {
  unit(0.5, "npc")
}

#########
# Y locations on edge
#########

yDetails <- function(x, theta) {
  UseMethod("yDetails")
}

yDetails.default <- function(x, theta) {
  unit(0.5, "npc")
}

#########
# WIDTHS
#########

# We are doing this in R code to provide generics like widthDetails
# so that users can customise the behaviour for complex grobs by
# writing their own (R code!) methods
width <- function(x) {
  widthDetails(x)
}

widthDetails <- function(x) {
  UseMethod("widthDetails", x)
}

widthDetails.default <- function(x) {
  unit(1, "null")
}

#########
# HEIGHTS
#########
height <- function(x) {
  heightDetails(x)
}

heightDetails <- function(x) {
  UseMethod("heightDetails", x)
}

heightDetails.default <- function(x) {
  unit(1, "null")
}

ascentDetails <- function(x) {
    UseMethod("ascentDetails", x)
}

ascentDetails.default <- heightDetails.default

ascentDetails.grob <- function(x) {
    heightDetails(x)
}

descentDetails <- function(x) {
    UseMethod("descentDetails", x)
}

descentDetails.default <- function(x) {
    unit(0, "inches")
}

#########
# Some functions that might be useful for determining the sizes
# of your grobs
#########

# Dimensions which depend on the parent context EITHER don't make
# sense (e.g., no good to have the parent width depend on the child's
# width unit(1, "grobwidth", <child>), which depends on the parent's
# width unit(.1, "npc"), ...) OR are slightly ambiguous
# (e.g., gf <- grid.frame(); grid.pack(gf, grid.rect(width=unit(.1, "npc")))
# makes the area allocated to the rectangle .1 of the frame area, but
# then the rectangle only occupies .1 of _that_ allocated area;  my head
# hurts !).  The first sort will actually lead to infinite loops so
# watch out for that;  the second sort I just don't want to have to deal with.
#
# On the other hand, dimensions which do not depend on the parent context
# are much easier to deal with (e.g., "inches", "cm", "lines", ...)
#
# So this function takes a unit and returns absolute values
# untouched and replaces other values with unit(1, "null")

absolute.size <- function(unit) {
  absolute.units(unit)
}

#  File src/library/grid/R/unit.R
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


# Create an object of class "unit"
# Simple units are of the form 'unit(1, "cm")' or 'unit(1L:3, "cm")' or
# 'unit(c(1,3,6), c("cm", "inch", "npc"))'
# More complicated units are of the form 'unit(1, "string", "a string")'
# or 'unit(1, "grob", a.grob)'
unit <- function(x, units, data=NULL) {
    # Used to throw error if !is.numeric(x), but this way
    # user can specify unit(NA, "npc") rather than
    # having to specify unit(as.numeric(NA), "npc")
    x <- as.numeric(x)
    units <- as.character(units)
    if (length(x) == 0 || length(units) == 0)
        stop("'x' and 'units' must have length > 0")
    valid.unit(x, units, recycle.data(data, FALSE, length(x), units))
}

valid.unit <- function(x, units, data) {
  valid.units <- valid.units(units)
  data <- valid.data(rep(units, length.out=length(x)), data)
  attr(x, "unit") <- units
  attr(x, "valid.unit") <- valid.units
  attr(x, "data") <- data
  class(x) <- "unit"
  x
}

grid.convert <- function(x, unitTo, axisFrom="x", typeFrom="location",
                         axisTo=axisFrom, typeTo=typeFrom,
                         valueOnly=FALSE) {
  .Deprecated("convertUnit")
  convertUnit(x, unitTo, axisFrom, typeFrom, axisTo, typeTo, valueOnly)
}

convertUnit <- function(x, unitTo, axisFrom="x", typeFrom="location",
                        axisTo=axisFrom, typeTo=typeFrom,
                        valueOnly=FALSE) {
  whatfrom <- match(axisFrom, c("x", "y")) - 1L +
    2L*(match(typeFrom, c("location", "dimension")) - 1L)
  whatto <- match(axisTo, c("x", "y")) - 1L +
    2L*(match(typeTo, c("location", "dimension")) - 1L)
  if (!is.unit(x))
    stop("'x' argument must be a unit object")
  if (is.na(whatfrom) || is.na(whatto))
    stop("invalid 'axis' or 'type'")
  value <- grid.Call(L_convert, x, as.integer(whatfrom),
                 as.integer(whatto), valid.units(unitTo))
  if (!valueOnly)
    unit(value, unitTo)
  else
    value
}

grid.convertX <- function(x, unitTo, valueOnly=FALSE) {
  .Deprecated("convertX")
  convertX(x, unitTo, valueOnly)
}

convertX <- function(x, unitTo, valueOnly=FALSE) {
  convertUnit(x, unitTo, "x", "location", "x", "location",
              valueOnly=valueOnly)
}

grid.convertY <- function(x, unitTo, valueOnly=FALSE) {
  .Deprecated("convertY")
  convertY(x, unitTo, valueOnly)
}

convertY <- function(x, unitTo, valueOnly=FALSE) {
  convertUnit(x, unitTo, "y", "location", "y", "location",
              valueOnly=valueOnly)
}

grid.convertWidth <- function(x, unitTo, valueOnly=FALSE) {
  .Deprecated("convertWidth")
  convertWidth(x, unitTo, valueOnly)
}

convertWidth <- function(x, unitTo, valueOnly=FALSE) {
  convertUnit(x, unitTo, "x", "dimension", "x", "dimension",
              valueOnly=valueOnly)
}

grid.convertHeight <- function(x, unitTo, valueOnly=FALSE) {
  .Deprecated("convertHeight")
  convertHeight(x, unitTo, valueOnly)
}

convertHeight <- function(x, unitTo, valueOnly=FALSE) {
  convertUnit(x, unitTo, "y", "dimension", "y", "dimension",
              valueOnly=valueOnly)
}

convertNative <- function(unit, dimension="x", type="location") {
  .Deprecated("convertUnit")
  convertUnit(unit, "native", dimension, type, dimension, type,
              valueOnly=TRUE)
}

# This is like the "convert" functions:  it evaluates units (immediately)
# in the current context
calcStringMetric <- function(text) {
    # .Call rather than .Call.graphics because it is a one-off calculation
    metric <- grid.Call(L_stringMetric, text)
    names(metric) <- c("ascent", "descent", "width")
    metric
}

# NOTE: the order of the strings in these conversion functions must
# match the order of the enums in ../src/grid.h
# AND in ../src/unit.c (see UnitTable)
# NOTE: ../src/unit.c also allows some pseudonyms (e.g., "in" for "inches")
.grid.unit.list <- c("npc", "cm", "inches", "lines",
                     "native", "null", "snpc", "mm",
                     "points", "picas", "bigpts",
                     "dida", "cicero", "scaledpts",
                     "strwidth", "strheight",
                     "strascent", "strdescent",
                     "vplayoutwidth", "vplayoutheight", "char",
                     "grobx", "groby", "grobwidth", "grobheight",
                     "grobascent", "grobdescent",
                     "mylines", "mychar", "mystrwidth", "mystrheight")

stringUnit <- function(unit) {
    unit %in% c("strwidth", "strheight", "strascent", "strdescent")
}

grobUnit <- function(unit) {
    unit %in% c("grobx", "groby", "grobwidth", "grobheight",
                "grobascent", "grobdescent")
}

dataUnit <- function(unit) {
    stringUnit(unit) | grobUnit(unit)
}

recycle.data <- function(data, data.per, max.n, units) {
    # FIRST STEP:  check that data needs to be recycled
    if (any(dataUnit(units))) {
        # VERY IMPORTANT:  Even if there is only one data specified
        # and/or only one data needed, we want this to be a LIST of
        # data values so that a single data and several data can be
        # handled equivalently
        # The test for whether it is only a single value currently
        # consists of a check for mode="character" (i.e., a single
        # string) or mode="expression" (i.e., a single expression)
        # or class="grob" (i.e., a single grob) or class="gPath"
        if (is.character(data) || is.language(data) ||
            is.grob(data) || inherits(data, "gPath"))
            data <- list(data)
        if (data.per)
            n <- max.n
        else
            n <- length(data)
        original <- data
        length(data) <- n
        if (length(original) < length(data)) {
            for (i in (length(original) + 1):length(data)) {
                data[[i]] <- original[[(i - 1) %% length(original) + 1]]
            }
        }
    }
    data
}

# Make sure that and "str*" and "grob*" units have data
valid.data <- function(units, data) {
    n <- length(units)
    str.units <- stringUnit(units)
    if (any(str.units))
        for (i in (1L:n)[str.units])
            if (!(length(data) >= i &&
                  (is.character(data[[i]]) || is.language(data[[i]]))))
                stop("no string supplied for 'strwidth/height' unit")
    # Make sure that a grob has been specified
    grob.units <- grobUnit(units)
    if (any(grob.units))
        for (i in (1L:n)[grob.units]) {
            if (!(length(data) >= i &&
                  (is.grob(data[[i]]) || inherits(data[[i]], "gPath") ||
                   is.character(data[[i]]))))
                stop("no 'grob' supplied for 'grobwidth/height' unit")
            if (is.character(data[[i]]))
                data[[i]] <- gPathDirect(data[[i]])
            if (inherits(data[[i]], "gPath"))
                if (depth(data[[i]]) > 1)
                    stop("'gPath' must have depth 1 in 'grobwidth/height' units")
        }
    # Make sure that where no data is required, the data is NULL
    if (!all(sapply(data[!(str.units | grob.units)], is.null)))
        stop("non-NULL value supplied for plain unit")
    data
}

valid.units <- function(units) {
  .Call(validUnits, units)
}

as.character.unit <- function(x, ...) {
  class(x) <- NULL
  paste0(x, attr(x, "unit"))
}

#########################
# UNIT ARITHMETIC STUFF
#########################

unit.arithmetic <- function(func.name, arg1, arg2=NULL) {
  ua <- list(fname=func.name, arg1=arg1, arg2=arg2)
  class(ua) <- c("unit.arithmetic", "unit")
  ua
}

Ops.unit <- function(e1, e2) {
  ok <- switch(.Generic, "+"=TRUE, "-"=TRUE, "*"=TRUE, FALSE)
  if (!ok)
    stop(gettextf("operator '%s' not meaningful for units", .Generic),
         domain = NA)
  if (.Generic == "*")
    # can only multiply a unit by a scalar
    if (nzchar(.Method[1L])) {
      if (nzchar(.Method[2L]))
        stop("only one operand may be a unit")
      else if (is.numeric(e2))
        # NOTE that we always put the scalar first
        # Use as.numeric to force e2 to be REAL
        unit.arithmetic(.Generic, as.numeric(e2), e1)
      else
        stop("non-unit operand must be numeric")
    } else {
      if (is.numeric(e1))
        # Use as.numeric to force e1 to be REAL
        unit.arithmetic(.Generic, as.numeric(e1), e2)
      else
        stop("non-unit operand must be numeric")
    }
  else
    # Check that both arguments are units
    if (nzchar(.Method[1L]) && nzchar(.Method[2L]))
      unit.arithmetic(.Generic, e1, e2)
    else
      stop("both operands must be units")
}

## <FIXME>
## The na.rm arg is ignored here, and the S3 groupGeneric is
## Summary(x, ...)
## </FIXME>
Summary.unit <- function(..., na.rm=FALSE) {
  # NOTE that this call to unit.c makes sure that arg1 is
  # a single unit object
  x <- unit.c(...)
  ok <- switch(.Generic, "max"=TRUE, "min"=TRUE, "sum"=TRUE, FALSE)
  if (!ok)
    stop(gettextf("'Summary' function '%s' not meaningful for units",
                  .Generic), domain = NA)
  unit.arithmetic(.Generic, x)
}

is.unit.arithmetic <- function(x) {
  inherits(x, "unit.arithmetic")
}

as.character.unit.arithmetic <- function(x, ...) {
  # bit too customised for my liking, but whatever ...
  # NOTE that paste coerces arguments to mode character hence
  # this will recurse.
  fname <- x$fname
  if (fname == "+" || fname == "-" || fname == "*")
    paste0(x$arg1, fname, x$arg2)
  else
    paste0(fname, "(", paste(x$arg1, collapse=", "), ")")
}

unit.pmax <- function(...) {

  select.i <- function(unit, i) {
    `[`(unit, i, top=FALSE)
  }

  x <- list(...)
  numargs <- length(x)
  if (numargs == 0L)
    stop("no arguments where at least one expected")
  # how long will the result be?
  maxlength <- 0L
  for (i in seq_len(numargs))
    if (length(x[[i]]) > maxlength)
      maxlength <- length(x[[i]])
  # maxlength guaranteed >= 1
  result <- max(unit.list.from.list(lapply(x, select.i, 1L)))
  if (maxlength > 1L)
      for (i in 2L:maxlength)
          result <- unit.c(result, max(unit.list.from.list(lapply(x, select.i, i))))
  result
}

unit.pmin <- function(...) {

  select.i <- function(unit, i) {
    `[`(unit, i, top=FALSE)
  }

  x <- list(...)
  numargs <- length(x)
  if (numargs == 0L)
    stop("Zero arguments where at least one expected")
  # how long will the result be?
  maxlength <- 0L
  for (i in seq_len(numargs))
    if (length(x[[i]]) > maxlength)
      maxlength <- length(x[[i]])
  # maxlength guaranteed >= 1
  result <- min(unit.list.from.list(lapply(x, select.i, 1L)))
  if (maxlength > 1L)
      for (i in 2L:maxlength)
          result <- unit.c(result, min(unit.list.from.list(lapply(x, select.i, i))))
  result
}

#########################
# UNIT LISTS
# The idea with these is to allow arbitrary combinations
# of unit objects and unit arithmetic objects
#########################

# create a unit list from a unit, unit.arithmetic, or unit.list object
unit.list <- function(unit) {
  if (is.unit.list(unit))
    unit
  else {
    l <- length(unit)
    result <- vector("list", l)
    for (i in seq_len(l))
      result[[i]] <- unit[i]
    class(result) <- c("unit.list", "unit")
    result
  }
}

is.unit.list <- function(x) {
  inherits(x, "unit.list")
}

as.character.unit.list <- function(x, ...) {
  l <- length(x)
  result <- character(l)
  for (i in seq_len(l))
    result[i] <- as.character(x[[i]])
  result
}

#########################
# These work on any sort of unit object
#########################

is.unit <- function(unit) {
  inherits(unit, "unit")
}

print.unit <- function(x, ...) {
  print(as.character(x), quote=FALSE, ...)
  invisible(x)
}

#########################
# Unit subsetting
#########################

# The idea of the "top" argument is to allow the function to
# know if it has been called from the command-line or from
# a previous (recursive) call to "[.unit" or "[.unit.arithmetic"
# this allows recycling beyond the end of the unit object
# except at the top level

# NOTE that "unit" and "data" attributes will be recycled
`[.unit` <- function(x, index, top=TRUE, ...) {
  this.length <- length(x)
  if (is.logical(index))
    index <- (1L:this.length)[index]
  # Allow for negative integer index
  if (any(index < 0)) {
      if (any(index > 0))
          stop("cannot mix signs of indices")
      else
          index <- (1L:this.length)[index]
  }
  if (top && any(index > this.length))
    stop("index out of bounds ('unit' subsetting)")
  cl <- class(x)
  units <- attr(x, "unit")
  valid.units <- attr(x, "valid.unit")
  data <- attr(x, "data")
  class(x) <- NULL
  # The line below may seem slightly odd, but it should only be
  # used to recycle values when this method is called to
  # subset an argument in a unit.arithmetic object
  x <- x[(index - 1) %% this.length + 1]
  attr(x, "unit") <- units[(index - 1) %% length(units) + 1]
  attr(x, "valid.unit") <- valid.units[(index - 1) %% length(valid.units) + 1]
  data.list <- data[(index - 1) %% length(data) + 1]
  attr(x, "data") <- data.list
  class(x) <- cl
  x
}

# NOTE that units will be recycled to the length of the largest
# of the arguments
`[.unit.arithmetic` <- function(x, index, top=TRUE, ...) {
  this.length <- length(x)
  if (is.logical(index))
    index <- (1L:this.length)[index]
  # Allow for negative integer index
  if (any(index < 0)) {
      if (any(index > 0))
          stop("cannot mix signs of indices")
      else
          index <- (1L:this.length)[index]
  }
  if (top && any(index > this.length))
    stop("index out of bounds (unit arithmetic subsetting)")

  repSummaryUnit <- function(x, n) {
      newUnits <- lapply(seq_len(n), function(z) { get(x$fname)(x$arg1) })
      class(newUnits) <- c("unit.list", "unit")
      newUnits
  }

  switch(x$fname,
         "+"=`[`(x$arg1, (index - 1) %% this.length + 1, top=FALSE) +
             `[`(x$arg2, (index - 1) %% this.length + 1, top=FALSE),
         "-"=`[`(x$arg1, (index - 1) %% this.length + 1, top=FALSE) -
             `[`(x$arg2, (index - 1) %% this.length + 1, top=FALSE),
         # Recycle multiplier if necessary
         "*"=x$arg1[(index - 1) %% length(x$arg1) + 1] *
             `[`(x$arg2, (index - 1) %% this.length + 1, top=FALSE),
         "min"=repSummaryUnit(x, length(index)),
         "max"=repSummaryUnit(x, length(index)),
         "sum"=repSummaryUnit(x, length(index)))
}

`[.unit.list` <- function(x, index, top=TRUE, ...) {
  this.length <- length(x)
  if (is.logical(index))
    index <- (1L:this.length)[index]
  # Allow for negative integer index
  if (any(index < 0)) {
      if (any(index > 0))
          stop("cannot mix signs of indices")
      else
          index <- (1L:this.length)[index]
  }
  if (top && any(index > this.length))
    stop("index out of bounds (unit list subsetting)")
  cl <- class(x)
  result <- unclass(x)[(index - 1) %% this.length + 1]
  class(result) <- cl
  result
}

# Write `[<-.unit` methods too ??

#########################
# str() method
#########################

# Should work fine on atomic units and on unit.list
# The problem arises with unit.arithmetic, which are stored as lists
# but act like vectors
# (e.g., report length greater than number of list components)
str.unit.arithmetic <- function(object, ...) {
    cat("Class 'unit.arithmetic' [1:", length(object), "] ", sep="")
    str(unclass(object), ...)
}

#########################
# "c"ombining unit objects
#########################

# NOTE that I have not written methods for c()
# because method dispatch occurs on the first argument to
# "c" so c(unit(...), ...) would come here, but c(whatever, unit(...), ...)
# would go who-knows-where.
# A particularly nasty example is:  c(1, unit(1, "npc")) which will
# produce the same result as c(1, 1)
# Same problem for trying to control c(<unit>, <unit.arithmetic>)
# versus c(<unit.arithmetic>, <unit>), etc ...

# If any arguments are unit.arithmetic or unit.list, then the result will be
# unit.list
unit.c <- function(...) {
    x <- list(...)
    if (!all(sapply(x, is.unit)))
        stop("it is invalid to combine 'unit' objects with other types")
    listUnit <- function(x) {
        inherits(x, "unit.list") ||
        inherits(x, "unit.arithmetic")
    }
    ual <- any(sapply(x, listUnit))
    if (ual)
        unit.list.from.list(x)
    else {
        values <- unlist(x)
        unitUnits <- function(x) {
            rep(attr(x, "unit"), length.out=length(x))
        }
        units <- unlist(lapply(x, unitUnits))
        unitData <- function(x) {
            data <- attr(x, "data")
            if (is.null(data))
                vector("list", length(x))
            else
                recycle.data(data, TRUE, length(x), unitUnits(x))
        }
        data <- do.call("c", lapply(x, unitData))
        unit(values, units, data=data)
    }
}

unit.list.from.list <- function(x) {
    result <- do.call("c", lapply(x, unit.list))
    class(result) <- c("unit.list", "unit")
    result
}

#########################
# rep'ing unit objects
#########################

rep.unit <- function(x, times=1, length.out=NA, each=1, ...) {
    if (length(x) == 0)
        stop("invalid 'unit' object")

    # Determine an approprite index, then call subsetting code
    repIndex <- rep(seq_along(x), times=times, length.out=length.out, each=each)
    x[repIndex, top=FALSE]
}

# Vestige from when rep() was not generic
unit.rep <- function (x, ...)
{
  warning("'unit.rep' has been deprecated in favour of a unit method for the generic rep function", domain = NA)
  rep(x, ...)
}

#########################
# Length of unit objects
#########################

length.unit <- function(x) {
  length(unclass(x))
}

length.unit.list <- function(x) {
  length(unclass(x))
}

length.unit.arithmetic <- function(x) {
  switch(x$fname,
         "+"=max(length(x$arg1), length(x$arg2)),
         "-"=max(length(x$arg1), length(x$arg2)),
         "*"=max(length(x$arg1), length(x$arg2)),
         "min" = 1L,
         "max" = 1L,
         "sum" = 1L)
}

# Vestige of when length was not generic
unit.length <- function(unit) {
   warning("'unit.length' has been deprecated in favour of a unit method for the generic length function", domain = NA)
   length(unit)
}

#########################
# Convenience functions
#########################

stringWidth <- function(string) {
    n <- length(string)
    if (is.language(string)) {
        data <- vector("list", n)
        for (i in 1L:n)
            data[[i]] <- string[i]
    } else {
        data <- as.list(as.character(string))
    }
    unit(rep(1, n), "strwidth", data=data)
}

stringHeight <- function(string) {
    n <- length(string)
    if (is.language(string)) {
        data <- vector("list", n)
        for (i in 1L:n)
            data[[i]] <- string[i]
    } else {
        data <- as.list(as.character(string))
    }
    unit(rep(1, n), "strheight", data=data)
}

stringAscent <- function(string) {
    n <- length(string)
    if (is.language(string)) {
        data <- vector("list", n)
        for (i in 1L:n)
            data[[i]] <- string[i]
    } else {
        data <- as.list(as.character(string))
    }
    unit(rep(1, n), "strascent", data=data)
}

stringDescent <- function(string) {
    n <- length(string)
    if (is.language(string)) {
        data <- vector("list", n)
        for (i in 1L:n)
            data[[i]] <- string[i]
    } else {
        data <- as.list(as.character(string))
    }
    unit(rep(1, n), "strdescent", data=data)
}

convertTheta <- function(theta) {
    if (is.character(theta))
        # Allow some aliases for common angles
        switch(theta,
               east=0,
               north=90,
               west=180,
               south=270,
               stop("invalid 'theta'"))
    else
        # Ensure theta in [0, 360)
        theta <- as.numeric(theta) %% 360
}

# grobX
grobX <- function(x, theta) {
    UseMethod("grobX", x)
}

grobX.grob <- function(x, theta) {
  unit(convertTheta(theta), "grobx", data=x)
}

grobX.gList <- function(x, theta) {
  unit(rep(convertTheta(theta), length(gList)), "grobx", data=x)
}

grobX.gPath <- function(x, theta) {
  unit(convertTheta(theta), "grobx", data=x)
}

grobX.default <- function(x, theta) {
  unit(convertTheta(theta), "grobx", data=gPathDirect(as.character(x)))
}

# grobY
grobY <- function(x, theta) {
    UseMethod("grobY", x)
}

grobY.grob <- function(x, theta) {
  unit(convertTheta(theta), "groby", data=x)
}

grobY.gList <- function(x, theta) {
  unit(rep(convertTheta(theta), length(gList)), "groby", data=x)
}

grobY.gPath <- function(x, theta) {
  unit(convertTheta(theta), "groby", data=x)
}

grobY.default <- function(x, theta) {
  unit(convertTheta(theta), "groby", data=gPathDirect(as.character(x)))
}

# grobWidth
grobWidth <- function(x) {
  UseMethod("grobWidth")
}

grobWidth.grob <- function(x) {
  unit(1, "grobwidth", data=x)
}

grobWidth.gList <- function(x) {
  unit(rep(1, length(gList)), "grobwidth", data=x)
}

grobWidth.gPath <- function(x) {
  unit(1, "grobwidth", data=x)
}

grobWidth.default <- function(x) {
  unit(1, "grobwidth", data=gPathDirect(as.character(x)))
}

# grobHeight
grobHeight <- function(x) {
  UseMethod("grobHeight")
}

grobHeight.grob <- function(x) {
  unit(1, "grobheight", data=x)
}

grobHeight.gList <- function(x) {
  unit(rep(1, length(gList)), "grobheight", data=x)
}

grobHeight.gPath <- function(x) {
  unit(1, "grobheight", data=x)
}

grobHeight.default <- function(x) {
  unit(1, "grobheight", data=gPathDirect(as.character(x)))
}

# grobAscent
grobAscent <- function(x) {
  UseMethod("grobAscent")
}

grobAscent.grob <- function(x) {
  unit(1, "grobascent", data=x)
}

grobAscent.gList <- function(x) {
  unit(rep(1, length(gList)), "grobascent", data=x)
}

grobAscent.gPath <- function(x) {
  unit(1, "grobascent", data=x)
}

grobAscent.default <- function(x) {
  unit(1, "grobascent", data=gPathDirect(as.character(x)))
}

# grobDescent
grobDescent <- function(x) {
  UseMethod("grobDescent")
}

grobDescent.grob <- function(x) {
  unit(1, "grobdescent", data=x)
}

grobDescent.gList <- function(x) {
  unit(rep(1, length(gList)), "grobdescent", data=x)
}

grobDescent.gPath <- function(x) {
  unit(1, "grobdescent", data=x)
}

grobDescent.default <- function(x) {
  unit(1, "grobdescent", data=gPathDirect(as.character(x)))
}

#########################
# Function to decide which values in a unit are "absolute" (do not depend
# on parent's drawing context or size)
#########################

# Only deals with unit of length() 1
absolute <- function(unit) {
  !is.na(match(attr(unit, "unit"),
               c("cm", "inches", "lines", "null",
                 "mm", "points", "picas", "bigpts",
                 "dida", "cicero", "scaledpts",
                 "strwidth", "strheight", "strascent", "strdescent", "char",
                 "mylines", "mychar", "mystrwidth", "mystrheight")))
}

# OLD absolute.unit
absolute.units <- function(unit) {
  UseMethod("absolute.units")
}

absolute.units.unit <- function(unit) {
  n <- length(unit)
  if (absolute(unit[1L]))
    abs.unit <- unit[1L]
  else
    abs.unit <- unit(1, "null")
  new.unit <- abs.unit
  count <- 1
  while (count < n) {
    count <- count + 1
    new.unit <- unit.c(new.unit, absolute.units(unit[count]))
  }
  new.unit
}

absolute.units.unit.list <- function(unit) {
  cl <- class(unit)
  abs.ul <- lapply(unit, absolute.units)
  class(abs.ul) <- cl
  abs.ul
}

absolute.units.unit.arithmetic <- function(unit) {
  switch(unit$fname,
         "+"=unit.arithmetic("+", absolute.units(unit$arg1),
           absolute.units(unit$arg2)),
         "-"=unit.arithmetic("-", absolute.units(unit$arg1),
           absolute.units(unit$arg2)),
         "*"=unit.arithmetic("*", unit$arg1, absolute.units(unit$arg2)),
         "min"=unit.arithmetic("min", absolute.units(unit$arg1)),
         "max"=unit.arithmetic("max", absolute.units(unit$arg1)),
         "sum"=unit.arithmetic("sum", absolute.units(unit$arg1)))
}


#  File src/library/grid/R/util.R
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


# Define a convenience function that is easy to call from C code
grid.top.level.vp <- function() {
  pushedvp(viewport(clip=TRUE, name="ROOT"))
}

# An accessor for getting at the grid global state structure
# to make debugging easier for me;  all I have to type is grid:::STATE()
STATE <- function() {
  get(".GRID.STATE", envir=.GridEvalEnv)
}

is.even <- function(x) x %% 2 == 0

is.odd <- function(x) !is.even(x)


grid.pretty <- function(range) {
  if (!is.numeric(range))
    stop("'range' must be numeric")
  .Call(L_pretty, range)
}

#  File src/library/grid/R/viewport.R
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


initvpAutoName <- function() {
  index <- 0
  function() {
    index <<- index + 1
    paste0("GRID.VP.", index)
  }
}

vpAutoName <- initvpAutoName()

# NOTE: The order of the elements in viewports and pushedvps are
# VERY IMPORTANT because the C code accesses them using constant
# indices (i.e., if you change the order here the world will end!
valid.viewport <- function(x, y, width, height, just,
                           gp, clip,
                           xscale, yscale, angle,
                           layout, layout.pos.row, layout.pos.col,
                           name) {
  if (length(x) > 1 || length(y) > 1 ||
      length(width) > 1 || length(height) > 1)
    stop("'x', 'y', 'width', and 'height' must all be units of length 1")
  if (!is.gpar(gp))
    stop("invalid 'gp' value")
  if (!is.logical(clip))
    clip <- switch(as.character(clip),
                   on=TRUE,
                   off=NA,
                   inherit=FALSE,
                   stop("invalid 'clip' value"))
  if (!is.numeric(xscale) || length(xscale) != 2 ||
      any(!is.finite(xscale)))
    stop("invalid 'xscale' in viewport")
  if (!is.numeric(yscale) || length(yscale) != 2 ||
      any(!is.finite(yscale)))
    stop("invalid 'yscale' in viewport")
  if (!is.numeric(angle) || length(angle) != 1 ||
      !is.finite(angle))
    stop("invalid 'angle' in viewport")
  if (!(is.null(layout) || is.layout(layout)))
    stop("invalid 'layout' in viewport")
  if (!is.null(layout.pos.row)) {
    layout.pos.row <- as.integer(range(layout.pos.row))
    if (any(!is.finite(layout.pos.row)))
      stop("invalid 'layout.pos.row' in viewport")
  }
  if (!is.null(layout.pos.col)) {
    layout.pos.col <- as.integer(range(layout.pos.col))
    if (any(!is.finite(layout.pos.col)))
      stop("invalid 'layout.pos.col' in viewport")
  }
  # If name is NULL then we give it a default
  # Otherwise it should be a valid R name
  if (is.null(name))
    name <- vpAutoName()
  # Put all the valid things first so that are found quicker
  vp <- list(x = x, y = y, width = width, height = height,
             justification = just,
             gp = gp,
             clip = clip,
             xscale = xscale,
             yscale = yscale,
             angle = angle,
             layout = layout,
             layout.pos.row = layout.pos.row,
             layout.pos.col = layout.pos.col,
             valid.just = valid.just(just),
             valid.pos.row = layout.pos.row,
             valid.pos.col = layout.pos.col,
             name=name)
  class(vp) <- "viewport"
  vp
}

# When a viewport is pushed, an internal copy is stored along
# with plenty of additional information relevant to the state
# at the time of being pushed (this is all used to return to this
# viewport without having to repush it)
pushedvp <- function(vp) {
  pvp <- c(vp, list(gpar = NULL,
                    trans = NULL,
                    widths = NULL,
                    heights = NULL,
                    width.cm = NULL,
                    height.cm = NULL,
                    rotation = NULL,
                    cliprect = NULL,
                    parent = NULL,
                    # Children of this pushedvp will be stored
                    # in an environment
                    children = new.env(hash=TRUE, parent=baseenv()),
                    # Initial value of 0 means that the viewport will
                    # be pushed "properly" the first time, calculating
                    # transformations, etc ...
                    devwidthcm = 0,
                    devheightcm = 0,
                    # This is down here because need to keep
                    # #defines in grid.h consistent with order here
                    parentgpar = NULL))
  class(pvp) <- c("pushedvp", class(vp))
  pvp
}

vpFromPushedvp <- function(pvp) {
  vp <- pvp[c("x", "y", "width", "height",
              "justification", "gp", "clip",
              "xscale", "yscale", "angle",
              "layout", "layout.pos.row", "layout.pos.col",
              "valid.just", "valid.pos.row", "valid.pos.col",
              "name")]
  class(vp) <- "viewport"
  vp
}

as.character.viewport <- function(x, ...) {
  paste0("viewport[", x$name, "]")
}

as.character.vpList <- function(x, ...) {
  paste0("(", paste(vapply(x, as.character, ""), collapse=", "), ")")
}

as.character.vpStack <- function(x, ...) {
  paste(vapply(x, as.character, ""), collapse="->")
}

as.character.vpTree <- function(x, ...) {
  paste(x$parent, x$children, sep="->")
}

print.viewport <- function(x, ...) {
  cat(as.character(x), "\n")
  invisible(x)
}

width.details.viewport <- function(x) {
  absolute.size(x$width)
}

height.details.viewport <- function(x) {
  absolute.size(x$height)
}

# How many "levels" in viewport object
depth <- function(vp) {
  UseMethod("depth")
}

depth.viewport <- function(vp) {
  1
}

depth.vpList <- function(vp) {
  # When pushed, the last element of the vpList is pushed last
  # so we are left whereever that leaves us
  depth(vp[[length(vp)]])
}

depth.vpStack <- function(vp) {
  # Elements in the stack may be vpStacks or vpLists or vpTrees
  # so need to sum all the depths
  sum(sapply(vp, depth, simplify=TRUE))
}

depth.vpTree <- function(vp) {
  # When pushed, the last element of the vpTree$children is
  # pushed last so we are left wherever that leaves us
  depth(vp$parent) + depth(vp$children[[length(vp$children)]])
}

depth.path <- function(path) {
  path$n
}

####################
# Accessors
####################

viewport.layout <- function(vp) {
  vp$layout
}

viewport.transform <- function(vp) {
  .Deprecated("current.transform")
}

####################
# Public Constructor
####################
viewport <- function(x = unit(0.5, "npc"),
                     y = unit(0.5, "npc"),
                     width = unit(1, "npc"),
                     height = unit(1, "npc"),
                     default.units = "npc",
                     just = "centre",
                     gp = gpar(),
                     clip = "inherit",
                     # FIXME: scales are only linear at the moment
                     xscale = c(0, 1),
                     yscale = c(0, 1),
                     angle = 0,
                     # Layout for arranging children of this viewport
                     layout = NULL,
                     # Position of this viewport in parent's layout
                     layout.pos.row = NULL,
                     layout.pos.col = NULL,
                     # This is down here to avoid breaking
                     # existing code
                     name=NULL) {
  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)
  if (!is.unit(width))
    width <- unit(width, default.units)
  if (!is.unit(height))
    height <- unit(height, default.units)
  valid.viewport(x, y, width, height, just,
                 gp, clip, xscale, yscale, angle,
                 layout, layout.pos.row, layout.pos.col, name)
}

is.viewport <- function(vp) {
  inherits(vp, "viewport")
}

#############
# Some classes derived from viewport
#############

viewportorpath <- function(x) {
    is.viewport(x) || inherits(x, "vpPath")
}

vpListFromList <- function(vps) {
  if (all(sapply(vps, viewportorpath, simplify=TRUE))) {
    class(vps) <- c("vpList", "viewport")
    vps
  } else {
    stop("only viewports allowed in 'vpList'")
  }
}

# Viewports will be pushed in parallel
vpList <- function(...) {
  vps <- list(...)
  vpListFromList(vps)
}

# Viewports will be pushed in series
vpStack <- function(...) {
  vps <- list(...)
  if (all(sapply(vps, viewportorpath, simplify=TRUE))) {
    class(vps) <- c("vpStack", "viewport")
    vps
  } else {
    stop("only viewports allowed in 'vpStack'")
  }
}

# Viewports will be pushed as a tree
vpTree <- function(parent, children) {
  if (viewportorpath(parent) && inherits(children, "vpList")) {
    tree <- list(parent=parent, children=children)
    class(tree) <- c("vpTree", "viewport")
    tree
  } else {
    stop("'parent' must be a viewport and 'children' must be a 'vpList' in 'vpTree'")
  }
}

# A function for setting all gpars for vpStack/List/Tree
# Used in size.R
setvpgpar <- function(vp) {
  UseMethod("setvpgpar")
}

setvpgpar.viewport <- function(vp) {
  if (!is.null(vp$gp))
    set.gpar(vp$gp)
}

setvpgpar.vpStack <- function(vp) {
  lapply(vp, setvpgpar)
}

setvpgpar.vpList <- function(vp) {
  setvpgpar(vp[[length(vp)]])
}

setvpgpar.vpTree <- function(vp) {
  setvpgpar(vp$parent)
  setvpgpar(vp$children)
}

#############
# Functions for creating "paths" of viewport names
#############
.grid.pathSep <- "::"

vpPathFromVector <- function(names) {
  n <- length(names)
  if (n < 1)
    stop("a viewport path must contain at least one viewport name")
  if (any(bad <- !is.character(names)))
      stop(ngettext(sum(bad),
                    "invalid viewport name",
                    "invalid viewport names"),
           domain = NA)
  path <- list(path=if (n==1) NULL else
               paste(names[seq_len(n-1L)], collapse=.grid.pathSep),
               name=names[n],
               n=n)
  class(path) <- c("vpPath", "path")
  path
}

vpPath <- function(...) {
  names <- c(...)
  vpPathFromVector(names)
}

# Create vpPath from string with embedded VpPathSep(s)
vpPathDirect <- function(path) {
  names <- unlist(strsplit(path, .grid.pathSep))
  vpPathFromVector(names)
}

as.character.path <- function(x, ...) {
  if (x$n == 1)
    x$name
  else
    paste(x$path, x$name, sep=.grid.pathSep)
}

print.path <- function(x, ...) {
  cat(as.character(x), "\n")
  invisible(x)
}

`[.vpPath` <- function(x, index, ...) {
  names <- unlist(strsplit(as.character(x), .grid.pathSep))[index]
  vpPathFromVector(names)
}

# Explode path$path
explodePath <- function(path) {
  unlist(strsplit(path, .grid.pathSep))
}


#############
# Some handy viewport functions
#############

# Create a viewport with margins given in number of lines
plotViewport <- function(margins=c(5.1, 4.1, 4.1, 2.1), ...) {
  margins <- rep(as.numeric(margins), length.out=4)
  viewport(x=unit(margins[2L], "lines"),
           width=unit(1, "npc") - unit(sum(margins[c(2,4)]), "lines"),
           y=unit(margins[1L], "lines"),
           height=unit(1, "npc") - unit(sum(margins[c(1,3)]), "lines"),
           just=c("left", "bottom"),
           ...)
}

# Create a viewport from data
# If xscale not specified then determine from x
# If yscale not specified then determine from y
dataViewport <- function(xData = NULL, yData = NULL,
                         xscale = NULL, yscale = NULL, extension = 0.05, ...)
{
    extension <- rep(extension, length.out = 2)
    if (is.null(xscale)) {
        if (is.null(xData))
            stop("must specify at least one of 'x' or 'xscale'")
        xscale <- extendrange(xData, f = extension[1L])
    }
    if (is.null(yscale)) {
        if (is.null(yData))
            stop("must specify at least one of 'y' or 'yscale'")
        yscale <- extendrange(yData, f = extension[2L])
    }
    viewport(xscale = xscale, yscale = yscale, ...)
}
#  File src/library/grid/R/zzz.R
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

## environment used for evaluation in the C code
## assigned here to protect from GC, but otherwise unused at R level
.GridEvalEnv <- new.env()

# This should be the only grid global variable(?)
# It contains the list of state structures corresponding to the
# state for each device.
# The state structures are stored in here so that they do not
# get garbage collected.
assign(".GRID.STATE", vector("list", 64L), envir = .GridEvalEnv)
## 64 comes from the maximum number of R devices allowed to be open at
## one time, see R_MaxDevices in Graphics.h.

.noGenerics <- TRUE

utils::globalVariables(c("n", "vp", "path"))

.onLoad <- function(libname, pkgname)
{
    ## want eval in C code to see unexported objects
    environment(.GridEvalEnv) <- asNamespace("grid")
    .Call(L_initGrid, .GridEvalEnv)
    .grid.loaded <<- TRUE
}

.onUnload <- function(libpath)
{
    if (.grid.loaded) {
        ## Kill all existing devices to avoid replay
        ## of display list which tries to run grid code
        ## Not very friendly to other registered graphics systems
        ## but its safety first for now
        if(length(.Devices) > 1L)
            warning("shutting down all devices when unloading 'grid' namespace",
                    call. = FALSE)
        graphics.off()
        .Call(L_killGrid)
    }
    library.dynam.unload("grid", libpath)
}

## .gridplot.hook <- function()
## {
##     pushViewport(viewport(width=unit(1, "npc") - unit(1, "lines"),
## 			  x=0, just="left"))
##     grid.text(paste("help(", ..nameEx, ")"),
## 	      x=unit(1, "npc") + unit(0.5, "lines"),
## 	      y=unit(0.8, "npc"), rot=90,
## 	      gp=gpar(col="orchid"))
## }
