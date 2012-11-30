#  sysBiolAlgClass.R
#  FBA and friends with R.
#
#  Copyright (C) 2010-2012 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of sybil.
#
#  Sybil is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Sybil is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sybil.  If not, see <http://www.gnu.org/licenses/>.


#------------------------------------------------------------------------------#
#                    definition of the class sysBiolAlg                        #
#------------------------------------------------------------------------------#

setClass(Class = "sysBiolAlg",
         representation(
             problem   = "optObj",
             algorithm = "character",
             nr        = "integer",
             nc        = "integer",
             fldind    = "integer"
         ),
         contains = "VIRTUAL",
         validity = sybil:::.validsysBiolAlg
)


#------------------------------------------------------------------------------#
#                              user constructor                                #
#------------------------------------------------------------------------------#

sysBiolAlg <- function(model,
                       algorithm = SYBIL_SETTINGS("ALGORITHM"),
                       prefix = "sysBiolAlg", sep = "_",
                       ...) {

    stopifnot(is(model, "modelorg"), is(algorithm, "character"))

    prob <- tryCatch(new(paste(prefix, algorithm, sep = sep), model, ...),
                     error = function(e) e)

    if (is(prob, "simpleError")) {
        stop(prob)
    }

    return(prob)

}


#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

setMethod(f = "initialize",
          signature = "sysBiolAlg",
          definition = function(.Object,
                                solver = SYBIL_SETTINGS("SOLVER"),
                                method = SYBIL_SETTINGS("METHOD"),
                                solverParm = SYBIL_SETTINGS("SOLVER_CTRL_PARM"),
                                sbalg, pType = "lp", scaling = NULL, fi, nCols,
                                nRows, mat, ub, lb, obj, rlb, rtype,
                                lpdir = "max", rub = NULL, ctype = NULL) {

              if ( (!missing(solver)) ||
                   (!missing(method)) ||
                   (!missing(sbalg)) ) {

                  stopifnot(is(solver, "character"),
                            is(method, "character"),
                            is(sbalg,  "character"),
                            is(pType,  "character"),
                            is(fi,     "numeric"),
                            is(nCols,  "numeric"),
                            is(nRows,  "numeric"),
                            is(mat,    "Matrix"),
                            is(ub,     "numeric"),
                            is(lb,     "numeric"),
                            is(obj,    "numeric"),
                            is(rlb,    "numeric"),
                            is(rtype,  "character"),
                            is(lpdir,  "character"))


                  # ---------------------------------------------
                  # build problem object
                  # ---------------------------------------------

                  lp <- optObj(solver = solver, method = method, pType = pType)
                  lp <- initProb(lp, nrows = nRows, ncols = nCols)

                  # ---------------------------------------------
                  # set control parameters
                  # ---------------------------------------------

                  if (!any(is.na(solverParm))) {
                      setSolverParm(lp, solverParm)
                  }


                  # ---------------------------------------------
                  # load problem data
                  # ---------------------------------------------

                  loadLPprob(lp,
                             nCols = nCols,
                             nRows = nRows,
                             mat   = mat,
                             ub    = ub,
                             lb    = lb,
                             obj   = obj,
                             rlb   = rlb,
                             rtype = rtype,
                             lpdir = lpdir,
                             rub   = rub,
                             ctype = ctype)


                  # ---------------------------------------------
                  # scaling
                  # ---------------------------------------------

                  if (!is.null(scaling)) {
                      scaleProb(lp, scaling)
                  }

                  .Object@problem   <- lp
                  .Object@algorithm <- sbalg
                  .Object@nr        <- as.integer(nRows)
                  .Object@nc        <- as.integer(nCols)
                  .Object@fldind    <- as.integer(fi)

                  validObject(.Object)

              }

              return(.Object)
          }
)


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# problem
setMethod("problem", signature(object = "sysBiolAlg"),
          function(object) {
              return(object@problem)
          }
)


# algorithm
setMethod("algorithm", signature(object = "sysBiolAlg"),
          function(object) {
              return(object@algorithm)
          }
)

setReplaceMethod("algorithm", signature = (object = "sysBiolAlg"),
                 function(object, value) {
                     object@algorithm <- value
                     return(object)
                 }
)


# nr
setMethod("nr", signature(object = "sysBiolAlg"),
          function(object) {
              return(object@nr)
          }
)

setReplaceMethod("nr", signature = (object = "sysBiolAlg"),
                 function(object, value) {
                     object@nr <- value
                     return(object)
                 }
)


# nc
setMethod("nc", signature(object = "sysBiolAlg"),
          function(object) {
              return(object@nc)
          }
)

setReplaceMethod("nc", signature = (object = "sysBiolAlg"),
                 function(object, value) {
                     object@nc <- value
                     return(object)
                 }
)


# fldind
setMethod("fldind", signature(object = "sysBiolAlg"),
          function(object) {
              return(object@fldind)
          }
)

setReplaceMethod("fldind", signature = (object = "sysBiolAlg"),
                 function(object, value) {
                     object@fldind <- value
                     return(object)
                 }
)


#------------------------------------------------------------------------------#
#                               other methods                                  #
#------------------------------------------------------------------------------#

setMethod("show", signature(object = "sysBiolAlg"),
    function(object) {
        cat("Algorithm type: ", algorithm(object), "\n", sep = "")
        cat("Slot problem:\n")
        show(problem(object))
        cat("Slot fldind:\n")
        str(fldind(object))
    }
)


#------------------------------------------------------------------------------#

setMethod("optimizeProb", signature(object = "sysBiolAlg"),
    function(object, react = NULL,
             lb = NULL,
             ub = NULL,
             obj_coef = NULL,
             lpdir = NA,
             resetChanges = TRUE,
             #prCmd = NULL, poCmd = NULL,
             prCmd = NA, poCmd = NA,
             prCil = NA, poCil = NA) {


        # check the argument react
        if (is.null(react)) {
            del <- FALSE
            obj <- FALSE
        }
        else {
            # if model is of class "sysBiolAlg", react is given by a
            # preceeding function
            stopifnot(is(react, "numeric"))
    
            if ( (is.null(lb)) || (is.null(ub)) ) {
                del <- FALSE
            }
            else {
                del <- TRUE
                stopifnot(is(lb, "numeric"),
                          is(ub, "numeric"),
                          length(lb) == length(react),
                          length(ub) == length(react))
            }
    
            # check argument obj_coef
            if (is.null(obj_coef)) {
                obj <- FALSE
            }
            else {
                if ( (is(obj_coef, "numeric")) &&
                     (length(obj_coef) == length(react)) ) {
                    obj <- TRUE
                }
                else {
                    stop("argument ", sQuote("obj_coef"), "must be numeric ",
                         " and of same length as argument react")
                }
            }
        }

        # check argument lpdir
        if ( (length(lpdir) > 1L) || (is.na(lpdir)) ) {
            ld <- FALSE
        }
        else {
            ld <- TRUE
            lpdir <- ifelse(lpdir == "max", "max", "min")
        }

    
        # -------------------------------------------------------------- #
        # modifications to problem object
        # -------------------------------------------------------------- #
    
        fi    <- fldind(object)
        lpmod <- problem(object)

        if (isTRUE(del)) {
            # store default lower and upper bounds
            lowb_tmp <- getColsLowBnds(lpmod, fi[react])
            uppb_tmp <- getColsUppBnds(lpmod, fi[react])
    
            # change bounds of fluxes in react
            check <- changeColsBnds(lpmod, fi[react], lb, ub)
        }
    
        if (isTRUE(obj)) {
            # store default objective function
            obj_tmp <- getObjCoefs(lpmod, fi[react])
            
            # change objective function
            check <- changeObjCoefs(lpmod, fi[react], obj_coef)
        }

        if (isTRUE(ld)) {
            # store default optimization direction
            ld_tmp <- getObjDir(lpmod)
            
            # change objective function
            check <- setObjDir(lpmod, lpdir)
        }
    
    
        # -------------------------------------------------------------- #
        # optimization
        # -------------------------------------------------------------- #
    
        # do some kind of preprocessing
        preP <- sybil:::.ppProcessing(lpprob  = lpmod,
                                      ppCmd   = prCmd,
                                      loopvar = prCil)

        lp_ok     <- solveLp(lpmod)
        lp_stat   <- getSolStat(lpmod)
        if (is.na(lp_stat)) {
            lp_stat <- lp_ok
        }

        lp_obj    <- getObjVal(lpmod)
        lp_fluxes <- getFluxDist(lpmod)

        # do some kind of postprocessing
        postP <- sybil:::.ppProcessing(lpprob = lpmod,
                                       ppCmd = poCmd,
                                       loopvar = poCil)
    
    
        # -------------------------------------------------------------- #
        # reset modifications
        # -------------------------------------------------------------- #
    
        if (isTRUE(resetChanges)) {
            # reset the default bounds
            if (isTRUE(del)) {
                check <- changeColsBnds(lpmod, fi[react], lowb_tmp, uppb_tmp)
            }
        
            # reset the default objective function
            if (isTRUE(obj)) {
                check <- changeObjCoefs(lpmod, fi[react], obj_tmp)
            }
    
            # reset the default optimization direction
            if (isTRUE(ld)) {
                check <- setObjDir(lpmod, ld_tmp)
            }
        }


        # -------------------------------------------------------------- #
        # store solution
        # -------------------------------------------------------------- #

        return(list(ok = lp_ok,
                    obj = lp_obj,
                    stat = lp_stat,
                    fluxes = lp_fluxes,
                    preP = preP,
                    postP = postP))

    }
)


