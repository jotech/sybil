#  optsol_fluxdelClass.R
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


# optsol_fluxdelClass


#------------------------------------------------------------------------------#
#                  definition of the class optsol_fluxdel                      #
#------------------------------------------------------------------------------#

setClass("optsol_fluxdel",
           representation(
               react_id  = "character",    # reaction id's of the original model
               allGenes  = "character",    # gene id's of the original model
               chlb      = "numeric",      # lower bound of changed fluxes/genes
               chub      = "numeric",      # upper bound of changed fluxes/genes
               dels      = "matrix",       # deleted fluxes (1 column)
               algorithm = "character"
           ),
           contains = "optsol_simpleFBA"
)


#------------------------------------------------------------------------------#
#                              user constructor                                #
#------------------------------------------------------------------------------#

optsol_fluxdel <- function(solver, nprob, lpdir, ncols, nrows, objf, fld, comb = 1) {
    if (missing(solver) ||
        missing(nprob)  ||
        missing(lpdir)  ||
        missing(ncols)  ||
        missing(nrows)  ||
        missing(objf)   ||
        missing(fld)
       ) {
        stop("Not enough arguments for creating an object of class optsol_fluxdel!")
    }

    if (fld == TRUE) {
        fldist <- fluxDistribution(0, ncols, nprob + 1)
    }
    else {
        fldist <- fluxDistribution(NA)
    }

    new("optsol_fluxdel",
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
        dels         = matrix(NA, nrow = nprob + 1, ncol = comb)
       )

}


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# react_id
setMethod("react_id", signature(object = "optsol_fluxdel"),
          function(object) {
              return(object@react_id)
          }
)

setReplaceMethod("react_id", signature = (object = "optsol_fluxdel"),
                 function(object, value) {
                     object@react_id <- value
                     return(object)
                 }
)


# allGenes
setMethod("allGenes", signature(object = "optsol_fluxdel"),
          function(object) {
              return(object@allGenes)
          }
)

setReplaceMethod("allGenes", signature = (object = "optsol_fluxdel"),
                 function(object, value) {
                     object@allGenes <- value
                     return(object)
                 }
)


# chlb
setMethod("chlb", signature(object = "optsol_fluxdel"),
          function(object) {
              return(object@chlb)
          }
)

setReplaceMethod("chlb", signature = (object = "optsol_fluxdel"),
                 function(object, value) {
                     object@chlb <- value
                     return(object)
                 }
)


# chub
setMethod("chub", signature(object = "optsol_fluxdel"),
          function(object) {
              return(object@chub)
          }
)

setReplaceMethod("chub", signature = (object = "optsol_fluxdel"),
                 function(object, value) {
                     object@chub <- value
                     return(object)
                 }
)


# deleted
setMethod("dels", signature(object = "optsol_fluxdel"),
          function(object) {
              return(object@dels)
          }
)

setReplaceMethod("dels", signature = (object = "optsol_fluxdel"),
                 function(object, value) {
                     object@dels <- value
                     return(object)
                 }
)


# algorithm
setMethod("algorithm", signature(object = "optsol_fluxdel"),
          function(object) {
              return(object@algorithm)
          }
)

setReplaceMethod("algorithm", signature = (object = "optsol_fluxdel"),
                 function(object, value) {
                     object@algorithm <- value
                     return(object)
                 }
)


#------------------------------------------------------------------------------#
#                               other methods                                  #
#------------------------------------------------------------------------------#

setMethod("ind2id", signature = (object = "optsol_fluxdel"),
                 function(object, slotN) {
                     out <- NULL
                     switch (slotN,
                     
                         "dels" = {
                             out <- apply(dels(object), 2,
                                          function(x) allGenes(object)[x]
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


setMethod("deleted", signature = (object = "optsol_fluxdel"),
                 function(object, i) {
                     value <- dels(object)[i,]
                     return(value)
                 }
)


setMethod("[", "optsol_fluxdel", function(x, i, j, ..., drop = FALSE) {

        if ((missing(i)) || (length(i) == 0)) {
            return(x)
        }

        if (max(i) > length(x)) {
            stop("subscript out of bounds")
        }

        slots <- slotNames(x)
        
        isO <- is(x)[1]
        
        newClass <- paste(isO, "(",
                          "solver = \"", solver(x), "\"",
                          sep = "")
        

        newClass <- paste(newClass, ", ",
                          "nprob = length(i)-1, ",
                          "lpdir = \"lp_dir(x)\", ",
                          "ncols = lp_num_cols(x), ",
                          "nrows = lp_num_rows(x), ",
                          "objf = \"obj_function(x)\", ",
                          "fld = ",
        sep = "")

        NC_fl <- FALSE
        if (nfluxes(x) > 1) {
            NC_fl <- TRUE
            newClass <- paste(newClass, TRUE, sep = "")
        }
        else {
            newClass <- paste(newClass, FALSE, sep = "")
        }

        if ("delmat" %in% slots) {
            dimdel <- dim(delmat(x))
            newClass <- paste(newClass, ", delrows = ", dimdel[1], ", delcols = ", dimdel[2], sep = "")
        }
        else {
            NC_delmat <- NA
        }

        newClass <- paste(newClass, ")", sep = "")

        newSol <- eval(parse(text = newClass))

        method(newSol)    <- method(x)[i]
        algorithm(newSol) <- algorithm(x)
        lp_obj(newSol)    <- lp_obj(x)[i]
        lp_ok(newSol)     <- lp_ok(x)[i]
        lp_stat(newSol)   <- lp_stat(x)[i]
        react_id(newSol)  <- react_id(x)
        allGenes(newSol)  <- allGenes(x)
        chlb(newSol)      <- chlb(x)[i]
        chub(newSol)      <- chub(x)[i]
        dels(newSol)      <- dels(x)[i, , drop = FALSE]

        if (isTRUE(NC_fl)) {
            fluxes(newSol) <- fluxes(x)[,i, drop = FALSE]
        }

        if ("fluxdels" %in% slots) {
            fluxdels(newSol) <- fluxdels(x)[i]
        }

        if ("hasEffect" %in% slots) {
            hasEffect(newSol) <- hasEffect(x)[i]
        }

        if ("delmat" %in% slots) {
            if (all(is.na(dels(newSol)[1,]))) {
                delmi <- dels(newSol)[-1,1]
                delmj <- dels(newSol)[-1,2]
            }
            else {
                delmi <- dels(newSol)[,1]
                delmj <- dels(newSol)[,2]
            }
            
            delmat(newSol) <- delmat(x)[
                                        as.character(delmi),
                                        as.character(delmj),
                                        drop = FALSE
                                       ]
        }

        return(newSol)


#         if ("num_of_prob" %in% slots) {
#             NC_num_of_prob <- length(i)
#             newClass <- paste(newClass, ", nprob = ", length(i)-1, sep = "")
#         }
#         else {
#             NC_num_prob <- NA
#         }
# 
#         if ("lp_dir" %in% slots) {
#             NC_lp_dir <- lp_dir(x)
#             newClass <- paste(newClass, ", lpdir = \"", lp_dir(x), "\"", sep = "")
#         }
#         else {
#             NC_lp_dir <- NA
#         }
# 
#         if ("lp_num_cols" %in% slots) {
#             NC_lp_num_cols <- lp_num_cols(x)
#             newClass <- paste(newClass, ", ncols = ", lp_num_cols(x), sep = "")
#         }
#         else {
#             NC_lp_num_cols <- NA
#         }
# 
#         if ("lp_num_rows" %in% slots) {
#             NC_lp_num_rows <- lp_num_rows(x)
#             newClass <- paste(newClass, ", nrows = ", lp_num_rows(x), sep = "")
#         }
#         else {
#             NC_lp_num_rows <- NA
#         }
# 
#         if ("obj_function" %in% slots) {
#             NC_obj_function <- obj_function(x)
#             newClass <- paste(newClass, ", objf = \"", obj_function(x), "\"", sep = "")
#         }
#         else {
#             NC_obj_function <- NA
#         }
# 
#         if ("fluxdist" %in% slots) {
#             if (is.na(fluxes(x))) {
#                 NC_fluxdist <- FALSE
#                 newClass <- paste(newClass, ", fld = ", FALSE, sep = "")
#             }
#             else {
#                 NC_fluxdist <- TRUE
#                 newClass <- paste(newClass, ", fld = ", TRUE, sep = "")
#             }
#         }
#         else {
#             NC_fluxdist <- NA
#         }
# 
#         if ("delmat" %in% slots) {
#             NC_delmat <- delmat(x)
#             dimdel <- dim(NC_delmat)
#             newClass <- paste(newClass, ", delrows = ", dimdel[1], ", delcols = ", dimdel[2], sep = "")
#         }
#         else {
#             NC_delmat <- NA
#         }
    }
)


