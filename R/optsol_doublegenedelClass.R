#  optsol_doublegenedelClass.R
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


# optsol_doublegenedelClass


#------------------------------------------------------------------------------#
#               definition of the class optsol_doublegenedel                   #
#------------------------------------------------------------------------------#

setClass("optsol_doublegenedel",
           contains = c("optsol_genedel", "optsol_doublefluxdel")
)


#------------------------------------------------------------------------------#
#                              user constructor                                #
#------------------------------------------------------------------------------#

optsol_doublegenedel <- function(solver, nprob, lpdir, nrows, ncols, delrows, delcols, objf, fld) {
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
        stop("Not enough arguments for creating an object of class optsol_doublegenedel!")
    }

    if (fld == TRUE) {
        fldist <- fluxDistribution(0, ncols, nprob + 1)
    }
    else {
        fldist <- fluxDistribution(NA)
    }

    new("optsol_doublegenedel",
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
        fluxdels     = list(nprob + 1),
        hasEffect    = logical(nprob + 1),
        delmat       = matrix(FALSE, nrow = delrows, ncol = delcols)
       )

}


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#                               other methods                                  #
#------------------------------------------------------------------------------#

setMethod("ind2id", signature = (object = "optsol_doublegenedel"),
                 function(object, slotN) {
                     out <- NULL
                     switch (slotN,
                     
                         "dels" = {
                             out <- apply(dels(object), 2,
                                          function(x) allGenes(object)[x]
                                    )
                         },
                         
                         "lethal" = {
                             out <- allGenes(object)[lethal(object)]
                         },
                         
                         "fluxdels" = {
                             out <- lapply(fluxdels(object),
                                           function(x) react_id(object)[x]
                                    )
                         },
                         
                         {
                             warning(paste("'", slotN, "' is not a valid slot!",
                                           sep = ""
                                    )
                             )
                         }
                     
                     )
                 
                     return(out)
                 }
)
