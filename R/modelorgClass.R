#  modelorgClass.R
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


# modelorgClass


#------------------------------------------------------------------------------#
#                      definition of the class modelorg                        #
#------------------------------------------------------------------------------#

# The class modelorg (the data part) is inspired by the data structure used
# in the COBRA Toolbox for the same purpose.

setClass("modelorg",
    representation(
         mod_desc     = "character",
         mod_name     = "character",   # model name
         mod_id       = "character",   # model id
         mod_compart  = "character",   # vector compartments
         met_num      = "integer",     # number of metabolites
         met_id       = "character",   # vector metabolite id's
         met_name     = "character",   # vector metabolite names
         met_comp     = "integer",     # vector the metabolites compartment
         met_single   = "logical",     # metabolites appearing only once in S
         met_de       = "logical",     # dead end metabolites
         rhs          = "numeric",     # vector right hand side zeros Sv = 0
         react_num    = "integer",     # number of reactions
         react_rev    = "logical",     # vector reversibilities
         react_id     = "character",   # vector reaction id's
         react_name   = "character",   # vector reaction names
         react_single = "logical",     # reactions using metabolites appearing only once in S
         react_de     = "logical",     # reactions using dead end metabolites
         S            = "Matrix",      # matrix S
         lowbnd       = "numeric",     # vector reactions lower bounds
         uppbnd       = "numeric",     # vector reactions upper bounds
         obj_coef     = "numeric",     # vector objective coefficients
         gprRules     = "character",
         genes        = "list",
         gpr          = "character",
         allGenes     = "character",
         rxnGeneMat   = "Matrix",
         subSys       = "Matrix"

    ),
    validity = .validmodelorg
)


#------------------------------------------------------------------------------#
#                              user constructor                                #
#------------------------------------------------------------------------------#

modelorg <- function(id, name) {
    if (missing(id) || missing(name)) {
        stop("Creating an object of class model needs name and id!")
    }
    id   <- as.character(id)
    name <- as.character(name)
    obj <- new("modelorg", id = id, name = name)
    return(obj)
}


#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

setMethod(f = "initialize",
          signature = "modelorg",
          definition = function(.Object, id, name) {

              if ( (!missing(id)) || (!missing(name)) ) {
                  .Object@mod_id   <- as.character(id)
                  .Object@mod_name <- as.character(name)
              }

              return(.Object)
          }
)


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# model id
setMethod("mod_id", signature(object = "modelorg"),
          function(object) {
              return(object@mod_id)
          }
)

setReplaceMethod("mod_id", signature = (object = "modelorg"),
                 function(object, value) {
                     object@mod_id <- value
                     return(object)
                 }
)


# model name
setMethod("mod_name", signature(object = "modelorg"),
          function(object) {
              return(object@mod_name)
          }
)

setReplaceMethod("mod_name", signature = (object = "modelorg"),
                 function(object, value) {
                     object@mod_name <- value
                     return(object)
                 }
)


# model description
setMethod("mod_desc", signature(object = "modelorg"),
          function(object) {
              return(object@mod_desc)
          }
)

setReplaceMethod("mod_desc", signature = (object = "modelorg"),
                 function(object, value) {
                     object@mod_desc <- value
                     return(object)
                 }
)


# model compartments
setMethod("mod_compart", signature(object = "modelorg"),
          function(object) {
              return(object@mod_compart)
          }
)

setReplaceMethod("mod_compart", signature(object = "modelorg"),
          function(object, value) {
              object@mod_compart <- value
              return(object)
          }
)


# number of metabolites
setMethod("met_num", signature(object = "modelorg"),
          function(object) {
              return(object@met_num)
          }
)

setReplaceMethod("met_num", signature(object = "modelorg"),
          function(object, value) {
              object@met_num <- value
              return(object)
          }
)


# metabolite id's
setMethod("met_id", signature(object = "modelorg"),
          function(object) {
              return(object@met_id)
          }
)

setReplaceMethod("met_id", signature(object = "modelorg"),
          function(object, value) {
              object@met_id <- value
              return(object)
          }
)


# metabolite names
setMethod("met_name", signature(object = "modelorg"),
          function(object) {
              return(object@met_name)
          }
)

setReplaceMethod("met_name", signature(object = "modelorg"),
          function(object, value) {
              object@met_name <- value
              return(object)
          }
)


# metabolites compartments
setMethod("met_comp", signature(object = "modelorg"),
          function(object) {
              return(object@met_comp)
          }
)

setReplaceMethod("met_comp", signature(object = "modelorg"),
          function(object, value) {
              object@met_comp <- value
              return(object)
          }
)


# singletons
setMethod("met_single", signature(object = "modelorg"),
          function(object) {
              return(object@met_single)
          }
)

setReplaceMethod("met_single", signature(object = "modelorg"),
          function(object, value) {
              object@met_single <- value
              return(object)
          }
)


# dead ends
setMethod("met_de", signature(object = "modelorg"),
          function(object) {
              return(object@met_de)
          }
)

setReplaceMethod("met_de", signature(object = "modelorg"),
          function(object, value) {
              object@met_de <- value
              return(object)
          }
)


# right hand site
setMethod("rhs", signature(object = "modelorg"),
          function(object) {
              return(object@rhs)
          }
)

setReplaceMethod("rhs", signature(object = "modelorg"),
          function(object, value) {
              object@rhs <- value
              return(object)
          }
)


# number of reactions
setMethod("react_num", signature(object = "modelorg"),
          function(object) {
              return(object@react_num)
          }
)

setReplaceMethod("react_num", signature(object = "modelorg"),
          function(object, value) {
              object@react_num <- value
              return(object)
          }
)


# reversibilities
setMethod("react_rev", signature(object = "modelorg"),
          function(object) {
              return(object@react_rev)
          }
)

setReplaceMethod("react_rev", signature(object = "modelorg"),
          function(object, value) {
              object@react_rev <- value
              return(object)
          }
)


# reaction id's
setMethod("react_id", signature(object = "modelorg"),
          function(object) {
              return(object@react_id)
          }
)

setReplaceMethod("react_id", signature(object = "modelorg"),
          function(object, value) {
              object@react_id <- value
              return(object)
          }
)


# reaction names
setMethod("react_name", signature(object = "modelorg"),
          function(object) {
              return(object@react_name)
          }
)

setReplaceMethod("react_name", signature(object = "modelorg"),
          function(object, value) {
              object@react_name <- value
              return(object)
          }
)


# singletons
setMethod("react_single", signature(object = "modelorg"),
          function(object) {
              return(object@react_single)
          }
)

setReplaceMethod("react_single", signature(object = "modelorg"),
          function(object, value) {
              object@react_single <- value
              return(object)
          }
)


# dead ends
setMethod("react_de", signature(object = "modelorg"),
          function(object) {
              return(object@react_de)
          }
)

setReplaceMethod("react_de", signature(object = "modelorg"),
          function(object, value) {
              object@react_de <- value
              return(object)
          }
)


# stoichiometric matrix
setMethod("S", signature(object = "modelorg"),
          function(object) {
              return(object@S)
          }
)

setReplaceMethod("S", signature(object = "modelorg"),
          function(object, value) {
              object@S <- value
              return(object)
          }
)

# number of non zero elements in S
setMethod("Snnz", signature(object = "modelorg"),
          function(object) {
              return(length(object@S@x))
          }
)

# slot dimension of stoichiometric matrix (SparseM)
setMethod("dim", signature(x = "modelorg"),
          function(x) {
              return(dim(x@S))
          }
)


# lower bounds
setMethod("lowbnd", signature(object = "modelorg"),
          function(object) {
              return(object@lowbnd)
          }
)

setReplaceMethod("lowbnd", signature(object = "modelorg"),
          function(object, value) {
              object@lowbnd <- value
              return(object)
          }
)


# upper bounds
setMethod("uppbnd", signature(object = "modelorg"),
          function(object) {
              return(object@uppbnd)
          }
)

setReplaceMethod("uppbnd", signature(object = "modelorg"),
          function(object, value) {
              object@uppbnd <- value
              return(object)
          }
)


# objective coefficient
setMethod("obj_coef", signature(object = "modelorg"),
          function(object) {
              return(object@obj_coef)
          }
)

setReplaceMethod("obj_coef", signature(object = "modelorg"),
          function(object, value) {
              object@obj_coef <- value
              return(object)
          }
)


# gprRules
setMethod("gprRules", signature(object = "modelorg"),
          function(object) {
              return(object@gprRules)
          }
)

setReplaceMethod("gprRules", signature(object = "modelorg"),
          function(object, value) {
              object@gprRules <- value
              return(object)
          }
)


# genes
setMethod("genes", signature(object = "modelorg"),
          function(object) {
              return(object@genes)
          }
)

setReplaceMethod("genes", signature(object = "modelorg"),
          function(object, value) {
              object@genes <- value
              return(object)
          }
)


# gpr associations
setMethod("gpr", signature(object = "modelorg"),
          function(object) {
              return(object@gpr)
          }
)

setReplaceMethod("gpr", signature(object = "modelorg"),
          function(object, value) {
              object@gpr <- value
              return(object)
          }
)


# list of all genes
setMethod("allGenes", signature(object = "modelorg"),
          function(object) {
              return(object@allGenes)
          }
)

setReplaceMethod("allGenes", signature(object = "modelorg"),
          function(object, value) {
              object@allGenes <- value
              return(object)
          }
)


# reaction to gene mapping
setMethod("rxnGeneMat", signature(object = "modelorg"),
          function(object) {
              return(object@rxnGeneMat)
          }
)

setReplaceMethod("rxnGeneMat", signature(object = "modelorg"),
          function(object, value) {
              object@rxnGeneMat <- value
              return(object)
          }
)


# reaction sub systems
setMethod("subSys", signature(object = "modelorg"),
          function(object) {
              return(object@subSys)
          }
)

setReplaceMethod("subSys", signature(object = "modelorg"),
          function(object, value) {
              object@subSys <- value
              return(object)
          }
)


#------------------------------------------------------------------------------#
#                               other methods                                  #
#------------------------------------------------------------------------------#

setMethod("show", signature(object = "modelorg"),
    function(object) {
        cat("model name:            ", mod_name(object), "\n")
        cat("number of compartments ", length(mod_compart(object)), "\n")
        lapply(mod_compart(object), function(x)
        cat("                       ", x, "\n")
              )
        cat("number of reactions:   ", react_num(object), "\n")
        cat("number of metabolites: ", met_num(object), "\n")
        if (length(allGenes) > 0) {
            cat("number of unique genes:", length(allGenes(object)), "\n")
        }
        cat("objective function:    ", printObjFunc(object), "\n")
    }
)

