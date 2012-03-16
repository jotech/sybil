#  optsol_robAnaClass.R
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


# optsol_robAnaClass


#------------------------------------------------------------------------------#
#                  definition of the class optsol_robAna                       #
#------------------------------------------------------------------------------#

setClass("optsol_robAna",
         representation(
              ctrlr    = "reactId",   # id of the control reaction,
              ctrlfl   = "numeric"    # fixed flux value for the control reaction
         ),
         contains = "optsol_simpleFBA"
)


#------------------------------------------------------------------------------#
#                              user constructor                                #
#------------------------------------------------------------------------------#

# optsol_robAnaClass
optsol_robAna <- function(solver, method, nprob, lpdir, ncols, nrows, objf, fld, cr, crf) {
    if (missing(solver) ||
        missing(method) ||
        missing(nprob)  ||
        missing(lpdir)  ||
        missing(ncols)  ||
        missing(nrows)  ||
        missing(objf)   ||
        missing(fld)    ||
        missing(cr)     ||
        missing(crf)
       ) {
        stop("Not enough arguments for creating an object of class optsol_robAna!")
    }

    if (fld == TRUE) {
        fldist <- fluxDistribution(0, ncols, nprob)
    }
    else {
        fldist <- fluxDistribution(NA)
    }

    new("optsol_robAna",
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
        ctrlr        = cr,
        ctrlfl       = as.numeric(crf)
       )

}


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# ctrlr
setMethod("ctrlr", signature(object = "optsol_robAna"),
          function(object) {
              return(object@ctrlr)
          }
)

setReplaceMethod("ctrlr", signature = (object = "optsol_robAna"),
                 function(object, value) {
                     object@ctrlr <- value
                     return(object)
                 }
)


# ctrlfl
setMethod("ctrlfl", signature(object = "optsol_robAna"),
          function(object) {
              return(object@ctrlfl)
          }
)

setReplaceMethod("ctrlfl", signature = (object = "optsol_robAna"),
                 function(object, value) {
                     object@ctrlfl <- value
                     return(object)
                 }
)


#------------------------------------------------------------------------------#
#                               other methods                                  #
#------------------------------------------------------------------------------#

setMethod("plot", signature(x = "optsol_robAna", y = "missing"),
          function(x, y,
                   xlab = "Control Flux",
                   ylab = "Objective Function",
                   type = "b",
                   pch = 20,
                   fillColorBg = "grey",
                   fillBg = TRUE,
                   ...) {
              plot(x@ctrlfl, x@lp_obj, type = "n", xlab = xlab, ylab = ylab)
              if (fillBg == TRUE) {
                  polygon(c(x@ctrlfl[1], x@ctrlfl, x@ctrlfl[length(x@ctrlfl)]),
                          c(min(x@lp_obj), x@lp_obj, min(x@lp_obj)),
                          col = fillColorBg, border = NA)
              }
              points(x@ctrlfl, x@lp_obj, type = type, pch = pch, ...)
          }
)
