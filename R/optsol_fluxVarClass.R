#  optsol_fluxVarClass.R
#  FBA and friends with R.
#
#  Copyright (C) 2010-2012 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#  
#  This file is part of SyBiL.
#
#  SyBiL is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  SyBiL is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with SyBiL.  If not, see <http://www.gnu.org/licenses/>.


# optsol_fluxVarClass


#------------------------------------------------------------------------------#
#                  definition of the class optsol_fluxVar                      #
#------------------------------------------------------------------------------#

setClass("optsol_fluxVar",
         representation(
              react    = "reactId"    # reactions to be analyzed
         ),
         contains = "optsol_simpleFBA"
)


#------------------------------------------------------------------------------#
#                              user constructor                                #
#------------------------------------------------------------------------------#

optsol_fluxVar <- function(solver, method, nprob, lpdir, ncols, nrows, objf, fld, rc) {
    if (missing(solver) ||
        missing(method) ||
        missing(nprob)  ||
        missing(lpdir)  ||
        missing(ncols)  ||
        missing(nrows)  ||
        missing(objf)   ||
        missing(fld)    ||
        missing(rc)
       ) {
        stop("Not enough arguments for creating an object of class optsol_fluxVar!")
    }

    if (fld == TRUE) {
        fldist <- fluxDistribution(0, ncols, nprob)
    }
    else {
        fldist <- fluxDistribution(NA)
    }

    new("optsol_fluxVar",
        solver       = as.character(solver),
        method       = as.character(method),
        num_of_prob  = as.integer(nprob),
        lp_num_cols  = as.integer(ncols),
        lp_num_rows  = as.integer(nrows),
        lp_obj       = numeric(nprob),
        lp_ok        = integer(nprob),
        lp_stat      = integer(nprob),
        lp_dir       = as.character(lpdir),
        obj_function = as.character(objf),
        fluxdist     = fldist,
        react        = rc
       )

}


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# react
setMethod("react", signature(object = "optsol_fluxVar"),
          function(object) {
              return(object@react)
          }
)

setReplaceMethod("react", signature = (object = "optsol_fluxVar"),
                 function(object, value) {
                     object@ctrlr <- value
                     return(object)
                 }
)


#------------------------------------------------------------------------------#
#                               other methods                                  #
#------------------------------------------------------------------------------#

###
# use barplot? or better hist allow user to set plot to false
###
#setGeneric("plotRangeVar", function(object, ...) standardGeneric("plotRangeVar"))
setMethod("plotRangeVar", signature(object = "optsol_fluxVar"),
          function(object, ...) {
              range <- (
                  (abs(maxSol(object, lp_obj)) - abs(minSol(object, lp_obj)))
                  /abs(maxSol(object, lp_obj))
              )
              #range <- (abs(maxSol(object, lp_obj)) - abs(minSol(object, lp_obj)))
              #blubber <- hist(range, breaks = 500, ...)
              blubber <- hist(range, ...)
              #barplot(range,
              #        names.arg = react_id(react(object)),
              #        las = 2,
              #        ...)
              #return(blubber)
              return(range)
          }
)


setMethod("minSol", signature(object = "optsol_fluxVar"),
          function(object, slot) {
              odds     <- seq(1, num_of_prob(object), 2)

              command  <- paste("is(",
                                deparse(substitute(slot)),
                                "(",
                                deparse(substitute(object)),
                                "))", sep = "")
              slottype <- eval(parse(text = command))

              if (is.na(match("matrix", slottype))) {
                  command <- paste(deparse(substitute(slot)),
                                   "(",
                                   deparse(substitute(object)),
                                   ")[odds]", sep = "")
              }
              else {
                  command <- paste(deparse(substitute(slot)),
                                   "(",
                                   deparse(substitute(object)),
                                   ")[,odds]", sep = "")
              }
              minimalSolutions <- eval(parse(text = command))
              return(minimalSolutions)
          }
)


setMethod("maxSol", signature(object = "optsol_fluxVar"),
          function(object, slot) {
              evans    <- seq(2, num_of_prob(object), 2)

              command  <- paste("is(",
                                deparse(substitute(slot)),
                                "(",
                                deparse(substitute(object)),
                                "))", sep = "")
              slottype <- eval(parse(text = command))

              if (is.na(match("matrix", slottype))) {
                  command <- paste(deparse(substitute(slot)),
                                   "(",
                                   deparse(substitute(object)),
                                   ")[evans]", sep = "")
              }
              else {
                  command <- paste(deparse(substitute(slot)),
                                   "(",
                                   deparse(substitute(object)),
                                   ")[,evans]", sep = "")
              }
              maximalSolutions <- eval(parse(text = command))
              return(maximalSolutions)
          }
)


#setMethod("[", signature = signature(x = "optsol_fluxVar"),
#    function(x, i, j, ..., drop = FALSE) {
#
#        if ((missing(i)) || (length(i) == 0)) {
#            return(x)
#        }
#
#        if (max(i) > length(x)) {
#            stop("subscript out of bounds")
#        }
#
#        slots <- slotNames(x)
#        
#        isO <- is(x)[1]
#        
#        newClass <- paste(isO, "(",
#                          "solver = \"", solver(x), "\"",
#                          sep = "")
#        
#
#        newClass <- paste(newClass, ", ",
#                          "nprob = length(i)-1, ",
#                          "lpdir = \"lp_dir(x)\", ",
#                          "ncols = lp_num_cols(x), ",
#                          "nrows = lp_num_rows(x), ",
#                          "objf = \"obj_function(x)\", ",
#                          "fld = ",
#        sep = "")
#
#        NC_fl <- FALSE
#        if (nfluxes(x) > 1) {
#            NC_fl <- TRUE
#            newClass <- paste(newClass, TRUE, sep = "")
#        }
#        else {
#            newClass <- paste(newClass, FALSE, sep = "")
#        }
#
#        if ("delmat" %in% slots) {
#            dimdel <- dim(delmat(x))
#            newClass <- paste(newClass, ", delrows = ", dimdel[1], ", delcols = ", dimdel[2], sep = "")
#        }
#        else {
#            NC_delmat <- NA
#        }
#
#        newClass <- paste(newClass, ")", sep = "")
#
#        newSol <- eval(parse(text = newClass))
#
#        method(newSol)    <- method(x)[i]
#        algorithm(newSol) <- algorithm(x)
#        lp_obj(newSol)    <- lp_obj(x)[i]
#        lp_ok(newSol)     <- lp_ok(x)[i]
#        lp_stat(newSol)   <- lp_stat(x)[i]
#        react_id(newSol)  <- react_id(x)
#        allGenes(newSol)  <- allGenes(x)
#        chlb(newSol)      <- chlb(x)[i]
#        chub(newSol)      <- chub(x)[i]
#        dels(newSol)      <- dels(x)[i, , drop = FALSE]
#
#        if (isTRUE(NC_fl)) {
#            fluxes(newSol) <- fluxes(x)[,i, drop = FALSE]
#        }
#
#        if ("fluxdels" %in% slots) {
#            fluxdels(newSol) <- fluxdels(x)[i]
#        }
#
#        if ("hasEffect" %in% slots) {
#            hasEffect(newSol) <- hasEffect(x)[i]
#        }
#
#        if ("delmat" %in% slots) {
#            if (all(is.na(dels(newSol)[1,]))) {
#                delmi <- dels(newSol)[-1,1]
#                delmj <- dels(newSol)[-1,2]
#            }
#            else {
#                delmi <- dels(newSol)[,1]
#                delmj <- dels(newSol)[,2]
#            }
#            
#            delmat(newSol) <- delmat(x)[
#                                        as.character(delmi),
#                                        as.character(delmj),
#                                        drop = FALSE
#                                       ]
#        }
#
#        return(newSol)
#
#
#    }
#)


setMethod("plot", signature(x = "optsol_fluxVar", "missing"),
          function(x, y,
                   ylim,
                   xlab = "",
                   ylab = "Value",
                   pch = 20,
                   col = "black",
                   collower, colupper, pchupper, pchlower,
                   dottedline = FALSE,
                   baseline = 0,
                   connect = TRUE,
                   colconnect = "black",
                   ...) {

              if (missing(ylim)) {
                  largest  <- ceiling(max(lp_obj(x)))
                  smallest <- floor(min(lp_obj(x)))
                  ylim <- c(smallest, largest)
              }
              else {
                  largest <- ylim[2]
                  smallest <- ylim[1]
              }
              
              if (missing(collower)) {
                  collower <- col
              }
              
              if (missing(colupper)) {
                  colupper <- col
              }

              if (missing(pchlower)) {
                  pchlower <- pch
              }

              if (missing(pchupper)) {
                  pchupper <- pch
              }

              num_dots <- num_of_prob(x)/2
              
              xdat <- rep(1:num_dots, 2)
              minfl <- lp_obj(x)[c(seq(1, num_of_prob(x), 2))]
              maxfl <- lp_obj(x)[c(seq(2, num_of_prob(x), 2))]
             
              
              plot(xdat, c(minfl, maxfl), xlab = xlab, ylab = ylab, ylim = ylim,
                   col = c(collower, colupper),
                   pch = c(pchlower, pchupper), ...)

              if (isTRUE(connect)) {
                  arrows(1:num_dots, minfl, 1:num_dots, maxfl,
                         length = 0, col = colconnect)
              }

              if (dottedline == TRUE) {
                  segments(c(1:num_dots), rep(smallest, num_dots),
                           c(1:num_dots), minfl,
                           lty = "dotted")
              }

              if (!is.na(baseline)) {
                  points(c(1, num_dots), c(baseline, baseline),
                         type = "s", lty = "dashed")
              }
              
          }
)
