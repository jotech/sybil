#  optsol_genedelClass.R
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


# optsol_genedelClass


#------------------------------------------------------------------------------#
#                  definition of the class optsol_genedel                      #
#------------------------------------------------------------------------------#

setClass("optsol_genedel",
           representation(
                fluxdels = "list",
                hasEffect = "logical"
           ),
           contains = "optsol_fluxdel"
)


#------------------------------------------------------------------------------#
#                              user constructor                                #
#------------------------------------------------------------------------------#

optsol_genedel <- function(solver, nprob, lpdir, ncols, nrows, objf, fld, comb = 1) {
    if (missing(solver) ||
        missing(nprob)  ||
        missing(lpdir)  ||
        missing(ncols)  ||
        missing(nrows)  ||
        missing(objf)   ||
        missing(fld)
       ) {
        stop("Not enough arguments for creating an object of class optsol_genedel!")
    }

    if (fld == TRUE) {
        fldist <- fluxDistribution(0, ncols, nprob + 1)
    }
    else {
        fldist <- fluxDistribution(NA)
    }

    new("optsol_genedel",
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
        dels         = matrix(NA, nrow = nprob + 1, ncol = comb),
        fluxdels     = list(nprob + 1),
        hasEffect    = logical(nprob + 1)
       )

}


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# reactions
setMethod("fluxdels", signature(object = "optsol_genedel"),
          function(object) {
              return(object@fluxdels)
          }
)

setReplaceMethod("fluxdels", signature = (object = "optsol_genedel"),
                 function(object, value) {
                     object@fluxdels <- value
                     return(object)
                 }
)


# hasEffect
setMethod("hasEffect", signature(object = "optsol_genedel"),
          function(object) {
              return(object@hasEffect)
          }
)

setReplaceMethod("hasEffect", signature = (object = "optsol_genedel"),
                 function(object, value) {
                     object@hasEffect <- value
                     return(object)
                 }
)


#------------------------------------------------------------------------------#
#                               other methods                                  #
#------------------------------------------------------------------------------#

setMethod("ind2id", signature = (object = "optsol_genedel"),
                 function(object, slotN) {
                     out <- NULL
                     switch (slotN,
                     
                         "dels" = {
                             out <- apply(dels(object), 2,
                                          function(x) allGenes(object)[x]
                                    )
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


setMethod("deleted", signature = (object = "optsol_genedel"),
                 function(object, i) {
                     value <- fluxdels(object)[[i]]
                     return(value)
                 }
)


setMethod("cpsol", signature(object = "optsol_genedel"),
          function(object, i, z, fld = FALSE) {

                     lp_obj(object)[i+1]    <- lp_obj(object)[z]
                     lp_ok(object)[i+1]     <- lp_ok(object)[z]
                     lp_stat(object)[i+1]   <- lp_stat(object)[z]
                     fluxdels(object)[i+1]  <- fluxdels(object)[z]
                     hasEffect(object)[i+1] <- hasEffect(object)[z]
                     if (identical(TRUE, fld)) {
                         fluxes(object)[,i+1]   <- fluxes(object)[,z]
                     }
                     return(object)

          }
)
