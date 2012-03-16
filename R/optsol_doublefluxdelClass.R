#  optsol_doublefluxdelClass.R
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


# optsol_doublefluxdelClass


#------------------------------------------------------------------------------#
#               definition of the class optsol_doublefluxdel                   #
#------------------------------------------------------------------------------#

setClass("optsol_doublefluxdel",
           representation(
                delmat   = "matrix", # boolean matrix showing the combination of deleted fluxes
                lethal   = "integer" # integer vector containing lethal fluxes/genes
           ),
           contains = "optsol_fluxdel"
)


#------------------------------------------------------------------------------#
#                              user constructor                                #
#------------------------------------------------------------------------------#

# optsol_doublefluxdelClass
optsol_doublefluxdel <- function(solver, nprob, lpdir, nrows, ncols, delrows, delcols, objf, fld) {
    if (missing(solver)  ||
        missing(nprob)   ||
        missing(lpdir)   ||
        missing(nrows)   ||
        missing(ncols)   ||
        missing(delrows) ||
        missing(delcols) ||
        missing(objf)    ||
        missing(fld)
       ) {
        stop("Not enough arguments for creating an object of class optsol_doublefluxdel!")
    }

    if (fld == TRUE) {
        fldist <- fluxDistribution(0, ncols, nprob + 1)
    }
    else {
        fldist <- fluxDistribution(NA)
    }

    new("optsol_doublefluxdel",
        solver       = as.character(solver),
        num_of_prob  = as.integer(nprob + 1),
        lp_num_cols  = as.integer(ncols),
        lp_num_rows  = as.integer(nrows),
        lp_obj       = numeric(nprob + 1),
        lp_ok        = integer(nprob + 1),
        lp_stat      = integer(nprob + 1),
        lp_dir       = as.character(lpdir),
        obj_function = as.character(objf),
        fluxdist     = fldist,
        dels         = matrix(NA, nrow = nprob + 1, ncol = 2),
        delmat       = matrix(FALSE, nrow = delrows, ncol = delcols)
       )

}


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# delmat
setMethod("delmat", signature(object = "optsol_doublefluxdel"),
          function(object) {
              return(object@delmat)
          }
)

setReplaceMethod("delmat", signature = (object = "optsol_doublefluxdel"),
                 function(object, value) {
                     object@delmat <- value
                     return(object)
                 }
)


# lethal
setMethod("lethal", signature(object = "optsol_doublefluxdel"),
          function(object) {
              return(object@lethal)
          }
)

setReplaceMethod("lethal", signature = (object = "optsol_doublefluxdel"),
                 function(object, value) {
                     object@lethal <- value
                     return(object)
                 }
)


#------------------------------------------------------------------------------#
#                               other methods                                  #
#------------------------------------------------------------------------------#

# print the results
setMethod("ddelres", signature(object = "optsol_doublefluxdel"),
          function(object, slot, react1, react2) {

              if (!is(react1, "numeric")) {
                  stop("react1 has to be numeric!")
              }
              if (!is(react2, "numeric")) {
                  stop("react2 has to be numeric!")
              }
              
              #a <- which(dels(object)[,1] == react1)
              #b <- which(dels(object)[,2] == react2)
              a <- lapply(react1, function(x) which(dels(object)[,1] == x))
              b <- lapply(react2, function(x) which(dels(object)[,2] == x))

              pos <- lapply(a, function(x) intersect(x, unlist(b)))
              pos <- unlist(pos)

              command <- paste("is(", deparse(substitute(slot)), "(", deparse(substitute(object)), "), 'matrix')", sep = "")
              result <- eval(parse(text = command))
              
              if (result == TRUE) {
                  command <- paste(deparse(substitute(slot)), "(", deparse(substitute(object)), ")[pos,]", sep = "")
              }
              else {
                  command <- paste(deparse(substitute(slot)), "(", deparse(substitute(object)), ")[pos]", sep = "")
              }

              result <- eval(parse(text = command))

              return(result)
          }
)
