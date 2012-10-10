#  optsol_simpleFBAClass.R
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


# optsol_simpleFBAClass


#------------------------------------------------------------------------------#
#                 definition of the class optsol_simpleFBA                     #
#------------------------------------------------------------------------------#

setClass("optsol_simpleFBA",
         representation(
              preProc  = "ppProc", # preprocessing lp result
              postProc = "ppProc"  # postprocessing lp result
        ),
        contains = "optsol"
)


#------------------------------------------------------------------------------#
#                              user constructor                                #
#------------------------------------------------------------------------------#

optsol_simpleFBA <- function(solver, nprob, lpdir, ncols, nrows, objf, fld) {
    
    .Deprecated(new = "optsol_optimizeProb", old = "optsol_simpleFBA")
    
    if (missing(solver) ||
        missing(nprob)  ||
        missing(lpdir)  ||
        missing(ncols)  ||
        missing(nrows)  ||
        missing(objf)   ||
        missing(fld)
       ) {
        stop("Not enough arguments for creating an object of class optsol_simpleFBA!")
    }

    if (fld == TRUE) {
        fldist <- fluxDistribution(0, ncols, nprob)
    }
    else {
        fldist <- fluxDistribution(NA)
    }

    new("optsol_simpleFBA",
        solver       = as.character(solver),
        num_of_prob  = as.integer(nprob),
        lp_num_cols  = as.integer(ncols),
        lp_num_rows  = as.integer(nrows),
        lp_obj       = numeric(nprob),
        lp_ok        = integer(nprob),
        lp_stat      = integer(nprob),
        lp_dir       = as.character(lpdir),
        #obj_function = as.character(objf),
        fluxdist     = fldist
       )

}


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# preProc
setMethod("preProc", signature(object = "optsol_simpleFBA"),
          function(object) {
              return(object@preProc)
          }
)

setReplaceMethod("preProc", signature = (object = "optsol_simpleFBA"),
                 function(object, value) {
                     object@preProc <- value
                     return(object)
                 }
)


# postProc
setMethod("postProc", signature(object = "optsol_simpleFBA"),
          function(object) {
              return(object@postProc)
          }
)

setReplaceMethod("postProc", signature = (object = "optsol_simpleFBA"),
                 function(object, value) {
                     object@postProc <- value
                     return(object)
                 }
)
