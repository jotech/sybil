#  optsol_blockedReactClass.R
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


# optsol_blockedReactClass


#------------------------------------------------------------------------------#
#               definition of the class optsol_blockedReact                    #
#------------------------------------------------------------------------------#

# slot obj_function is used here for the optimized reaction
setClass("optsol_blockedReact",
         representation(
              blocked      = "logical",           # blocked reaction (yes/no)
              react        = "reactId"            # checked reaction id's
        ),
        contains = "optsol"
)


#------------------------------------------------------------------------------#
#                              user constructor                                #
#------------------------------------------------------------------------------#

optsol_blockedReact <- function(solver, nprob, ncols, nrows, fld) {
    if (missing(solver) ||
        missing(nprob)  ||
        missing(ncols)  ||
        missing(nrows)  ||
        missing(fld)
       ) {
        stop("Not enough arguments for creating an object of class optsol_blockedReact!")
    }

    if (fld == TRUE) {
        fldist <- fluxDistribution(0, (ncols * 2), nprob)
    }
    else {
        fldist <- fluxDistribution(NA)
    }

    new("optsol_blockedReact",
        solver       = as.character(solver),
        num_of_prob  = as.integer(nprob),
        lp_num_cols  = as.integer(ncols),
        lp_num_rows  = as.integer(nrows),
        lp_obj       = numeric(nprob),
        lp_ok        = integer(nprob),
        lp_stat      = integer(nprob),
        lp_dir       = "max/min",
        obj_function = "identify blocked reactions",
        fluxdist     = fldist,
        blocked      = logical(ncols)
       )

}


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# blocked
setMethod("blocked", signature(object = "optsol_blockedReact"),
          function(object) {
              return(object@blocked)
          }
)

setReplaceMethod("blocked", signature = (object = "optsol_blockedReact"),
                 function(object, value) {
                     object@blocked <- value
                     return(object)
                 }
)


# react
setMethod("react", signature(object = "optsol_blockedReact"),
          function(object) {
              return(object@react)
          }
)

setReplaceMethod("react", signature = (object = "optsol_blockedReact"),
                 function(object, value) {
                     object@react <- value
                     return(object)
                 }
)


#------------------------------------------------------------------------------#
#                               other methods                                  #
#------------------------------------------------------------------------------#

setMethod("maxSol", signature(object = "optsol_blockedReact"),
          function(object, slot) {
              odds    <- seq(1, num_of_prob(object), 2)
              command <- paste(deparse(substitute(slot)), "(", deparse(substitute(object)), ")[odds]", sep = "")
              minimalSolutions <- eval(parse(text = command))
              return(minimalSolutions)
          }
)

setMethod("minSol", signature(object = "optsol_blockedReact"),
          function(object, slot) {
              odds    <- seq(2, num_of_prob(object), 2)
              command <- paste(deparse(substitute(slot)), "(", deparse(substitute(object)), ")[odds]", sep = "")
              minimalSolutions <- eval(parse(text = command))
              return(minimalSolutions)
          }
)
