#  optsolClass.R
#  FBA and friends with R.
#
#  Copyright (C) 2010-2012 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of sybil.
#
#  sybil is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  sybil is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sybil.  If not, see <http://www.gnu.org/licenses/>.


# optsolClass


#------------------------------------------------------------------------------#
#                            class definitions                                 #
#------------------------------------------------------------------------------#

setClass("optsol",
    representation(
        mod_id       = "character",        # model id of the original model
        solver       = "character",        # the used lp solver
        method       = "character",        # the used method
        algorithm    = "character",        # the used algorithm
        num_of_prob  = "integer",          # number of problems to solve
        lp_num_cols  = "integer",          # number of reactions
        lp_num_rows  = "integer",          # number of metabolites
        lp_obj       = "numeric",          # solution of the objective function
        lp_ok        = "integer",          # exit status of the lp solver
        lp_stat      = "integer",          # solution status
        lp_dir       = "character",        # direction
        obj_coef     = "numeric",          # objective coefficients in model
        fldind       = "integer",          # indices of fluxes
        fluxdist     = "fluxDistribution"  # the flux distribution
    ),
    contains = "VIRTUAL",
    #validity = sybil:::.validoptsol
)


#------------------------------------------------------------------------------#
#                              user constructors                               #
#------------------------------------------------------------------------------#

optsol <- function(solver) {
    if (missing(solver)) {
        stop("Creating an object of class optsol needs a valid solver!")
    }
    solver <- as.character(solver)
    new("optsol", solver = solver)
}


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# mod_id
setMethod("mod_id", signature(object = "optsol"),
          function(object) {
              return(object@mod_id)
          }
)

setReplaceMethod("mod_id", signature = (object = "optsol"),
                 function(object, value) {
                     object@mod_id <- value
                     return(object)
                 }
)


# solver
setMethod("solver", signature(object = "optsol"),
          function(object) {
              return(object@solver)
          }
)

setReplaceMethod("solver", signature = (object = "optsol"),
                 function(object, value) {
                     object@solver <- value
                     return(object)
                 }
)


# method
setMethod("method", signature(object = "optsol"),
          function(object) {
              return(object@method)
          }
)

setReplaceMethod("method", signature = (object = "optsol"),
                 function(object, value) {
                     object@method <- value
                     return(object)
                 }
)


# method
setMethod("algorithm", signature(object = "optsol"),
          function(object) {
              return(object@algorithm)
          }
)

setReplaceMethod("algorithm", signature = (object = "optsol"),
                 function(object, value) {
                     object@algorithm <- value
                     return(object)
                 }
)


# num_of_prob
setMethod("num_of_prob", signature(object = "optsol"),
          function(object) {
              return(object@num_of_prob)
          }
)

setReplaceMethod("num_of_prob", signature = (object = "optsol"),
                 function(object, value) {
                     object@num_of_prob <- value
                     return(object)
                 }
)


# lp_num_cols
setMethod("lp_num_cols", signature(object = "optsol"),
          function(object) {
              return(object@lp_num_cols)
          }
)

setReplaceMethod("lp_num_cols", signature = (object = "optsol"),
                 function(object, value) {
                     object@lp_num_cols <- value
                     return(object)
                 }
)


# lp_num_rows
setMethod("lp_num_rows", signature(object = "optsol"),
          function(object) {
              return(object@lp_num_rows)
          }
)

setReplaceMethod("lp_num_rows", signature = (object = "optsol"),
                 function(object, value) {
                     object@lp_num_rows <- value
                     return(object)
                 }
)


# lp_dir
setMethod("lp_dir", signature(object = "optsol"),
          function(object) {
              return(object@lp_dir)
          }
)

setReplaceMethod("lp_dir", signature = (object = "optsol"),
                 function(object, value) {
                     object@lp_dir <- value
                     return(object)
                 }
)


# lp_obj
setMethod("lp_obj", signature(object = "optsol"),
          function(object) {
              return(object@lp_obj)
          }
)

setReplaceMethod("lp_obj", signature = (object = "optsol"),
                 function(object, value) {
                     object@lp_obj <- value
                     return(object)
                 }
)


# lp_ok
setMethod("lp_ok", signature(object = "optsol"),
          function(object) {
              return(object@lp_ok)
          }
)

setReplaceMethod("lp_ok", signature = (object = "optsol"),
                 function(object, value) {
                     object@lp_ok <- value
                     return(object)
                 }
)


# lp_stat
setMethod("lp_stat", signature(object = "optsol"),
          function(object) {
              return(object@lp_stat)
          }
)

setReplaceMethod("lp_stat", signature = (object = "optsol"),
                 function(object, value) {
                     object@lp_stat <- value
                     return(object)
                 }
)


# objective coefficient
setMethod("obj_coef", signature(object = "optsol"),
          function(object) {
              return(object@obj_coef)
          }
)

setReplaceMethod("obj_coef", signature(object = "optsol"),
          function(object, value) {
              object@obj_coef <- value
              return(object)
          }
)


# fldind
setMethod("fldind", signature(object = "optsol"),
          function(object) {
              return(object@fldind)
          }
)

setReplaceMethod("fldind", signature = (object = "optsol"),
                 function(object, value) {
                     object@fldind <- value
                     return(object)
                 }
)


# fluxdist
setMethod("fluxdist", signature(object = "optsol"),
          function(object) {
              return(object@fluxdist)
          }
)

setReplaceMethod("fluxdist", signature = (object = "optsol"),
                 function(object, value) {
                     object@fluxdist <- value
                     return(object)
                 }
)


# fluxes
setMethod("fluxes", signature(object = "optsol"),
          function(object) {
              return(fluxes(object@fluxdist))
          }
)

setReplaceMethod("fluxes", signature = (object = "optsol"),
                 function(object, value) {
                     fluxes(object@fluxdist) <- value
                     return(object)
                 }
)


#------------------------------------------------------------------------------#
#                               other methods                                  #
#------------------------------------------------------------------------------#

# number of fluxes
setMethod("nfluxes", signature(object = "optsol"),
          function(object) {
              return(num_of_fluxes(object@fluxdist))
          }
)


# check solution status
setMethod("checkStat", signature(object = "optsol"),
          function(object) {
              ng <- checkSolStat(object@lp_stat, object@solver)
              return(ng)
          }
)


# consider using sprintf here
setMethod("show", signature(object = "optsol"),
    function(object) {
        
        cat("solver:                                  ", solver(object), "\n")
        cat("method:                                  ", method(object), "\n")
        cat("algorithm:                               ", algorithm(object), "\n")
        cat("number of variables:                     ", lp_num_cols(object), "\n")
        cat("number of constraints:                   ", lp_num_rows(object), "\n")
        cat("number of problems to solve:             ", num_of_prob(object), "\n")
        ok <- sum(lp_ok(object) == 0, na.rm = TRUE)
        cat("number of successful solution processes: ", ok, "\n")
    }
)


# length of an object of class optsol
setMethod("length", signature = signature(x = "optsol"),
          function(x) {
              return(num_of_prob(x))
          }
)


# draw a histogramm (package lattice)
setMethod("histogram", signature(x = "optsol"),
          function(x,
                   main = "",
                   xlab = "value of objective function",
                   ...) {


              histogram(lp_obj(x), main = main, xlab = xlab, ...)
              
          }
)


