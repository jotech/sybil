#  optObjClass.R
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


#------------------------------------------------------------------------------#
#                       definition of the class optObj                         #
#------------------------------------------------------------------------------#

setClass(Class = "lpExtPtr")       # lpSolveAPI
setClass(Class = "glpkPtr")        # glpkAPI
setClass(Class = "clpPtr")         # clpAPI
setClass(Class = "cplexPtr")       # cplexAPI

setClass(Class = "cplexPointer",
         representation(
             env = "cplexPtr",
             lp  = "cplexPtr"
         )
)

setClassUnion(name    = "pointerToProb",
              members = c("externalptr",
                          "lpExtPtr",
                          "glpkPtr",
                          "clpPtr",
                          "cplexPointer"
             )
)

setClass(Class = "optObj",
         representation(
              oobj     = "pointerToProb",
              solver   = "character",
              method   = "character",
              probType = "character"
         )
)

# derivatives
#setClass(Class = "optObj_boot", contains = "optObj")


#------------------------------------------------------------------------------#
#                              user constructor                                #
#------------------------------------------------------------------------------#

optObj <- function(solver, method, pType = "lp") {

    if (missing(solver)) {
        stop("Creating an object of class optObj needs a solver!")
    }

    if (missing(method)) {
        method <- checkDefaultMethod(solver,
                                     SYBIL_SETTINGS("METHOD"))[["met"]]
    }

    solver <- as.character(solver)
    method <- as.character(method)
    pType  <- as.character(pType)

    slv_obj <- switch(solver,
        # ----------------------- #
        "glpk" = {
            checkPackage <- require("glpkAPI")
            if(isTRUE(checkPackage)) {
                "optObj_glpk"
            }
            else {
                stop("Package glpkAPI not found.")
            }
        },
        # ----------------------- #
        "clp" = {
            checkPackage <- require("clpAPI")
            if(isTRUE(checkPackage)) {
                "optObj_clp"
            }
            else {
                stop("Package clpAPI not found.")
            }
        },
        # ----------------------- #
        "lpSolveAPI" = {
            checkPackage <- require("lpSolveAPI")
            if(isTRUE(checkPackage)) {
                method <- "lp_solve"
                "optObj_lpSolveAPI"
            }
            else {
                stop("Package lpSolveAPI not found.")
            }
        },
        # ----------------------- #
        "cplex" = {
            checkPackage <- require("cplexAPI")
            if(isTRUE(checkPackage)) {
                 "optObj_cplex"
           }
            else {
                stop("Package cplexAPI not found.")
            }
        },
        # ----------------------- #
#         "boot" = {
#             checkPackage <- require("boot")
#             if(isTRUE(checkPackage)) {
#                 "optObj_boot"
#             }
#             else {
#                 stop("Package boot not found.")
#             }
#         },
        # ----------------------- #
        {
            stop("not a valid solver")
        }
    )

    obj <- new(slv_obj, sv = solver, mt = method, pt = pType)
    return(obj)
}


#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

# contructor for class optObj
setMethod(f = "initialize",
          signature = "optObj",
          definition = function(.Object, sv, mt, pt) {

              if ( (!missing(sv)) || (!missing(mt)) || (!missing(pt)) ) {
                  
                  .Object@solver   <- as.character(sv)
                  .Object@method   <- as.character(mt)
                  .Object@probType <- as.character(pt)
                  
              }
              return(.Object)
          }
)

# contructor for class cplexPointer
setMethod(f = "initialize",
          signature = "cplexPointer",
          definition = function(.Object, en, pr) {

              if ( (!missing(en)) || (!missing(pr)) ) {
                  if ( (cplexAPI::isCPLEXenvPointer(en)) &&
                       (cplexAPI::isCPLEXprobPointer(pr)) ) {
                  
                      .Object@env <- en
                      .Object@lp  <- pr

                  }
              }
              return(.Object)
          }
)


#------------------------------------------------------------------------------#
#                                  getters                                     #
#------------------------------------------------------------------------------#

# solver
setMethod("solver", signature(object = "optObj"),
          function(object) {
              return(object@solver)
          }
)

# setReplaceMethod("solver", signature = (object = "optObj"),
#                  function(object, value) {
#                      object@solver <- value
#                      return(object)
#                  }
# )


# method
setMethod("method", signature(object = "optObj"),
          function(object) {
              return(object@method)
          }
)


# probType
setMethod("probType", signature(object = "optObj"),
          function(object) {
              return(object@probType)
          }
)


#------------------------------------------------------------------------------#
#                               other methods                                  #
#------------------------------------------------------------------------------#

# get the current dimension of the constraint matrix
setMethod("dim", "optObj",

    function(x) {

        out <- c(0, 0)
        out[1] <- getNumRows(x)
        out[2] <- getNumCols(x)

        return(out)
    }
)


#------------------------------------------------------------------------------#
# oobj

# all further methods are interface methods to the lp solvers

setMethod("delProb", signature(lp = "optObj"),

    function(lp, closeEnv = TRUE) {

        out <- FALSE

        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                glpkAPI::delProbGLPK(lp@oobj)
                out <- TRUE
            },
            # ----------------------- #
            "clp" = {
                clpAPI::delProbCLP(lp@oobj)
                out <- TRUE
            },
            # ----------------------- #
            "lpSolveAPI" = {
                #finalizeLpSolveProb(lp@oobj)
                #lpSolveAPI::delete.lp(lp@oobj)
                #lp@oobj <- NULL
                lp <- new("optObj_lpSolveAPI")
                out <- TRUE
            },
            # ----------------------- #
            "cplex" = {
                if (isTRUE(closeEnv)) {
                    cplexAPI::closeProbCPLEX(list(env = lp@oobj@env,
                                                  lp = lp@oobj@lp))
                } else {
                    cplexAPI::delProbCPLEX(lp@oobj@env, lp@oobj@lp)
                }
                out <- TRUE
            },
            # ----------------------- #
            {
                warning("not a valid instance of class optObj")
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("initProb", signature(lp = "optObj"),

    function(lp, nrows = 0, ncols = 0) {

        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                lp@oobj <- glpkAPI::initProbGLPK()
                glpkAPI::termOutGLPK(glpkAPI::GLP_OFF)
                #glpkAPI::setSimplexDefaultParmGLPK()
            },
            # ----------------------- #
            "clp" = {
                lp@oobj <- clpAPI::initProbCLP()
                clpAPI::setLogLevelCLP(lp@oobj, 0)
            },
            # ----------------------- #
            "lpSolveAPI" = {
                lp@oobj <- lpSolveAPI::make.lp(nrow = nrows, ncol = ncols)
                #reg.finalizer(lp@oobj, finalizeLpSolveProb, TRUE)
            },
            # ----------------------- #
            "cplex" = {
                #lp@oobj <- cplexAPI::openProbCPLEX()
                tmp <- cplexAPI::openProbCPLEX()
                lp@oobj <- new("cplexPointer",
                               en = tmp[["env"]],
                               pr = tmp[["lp"]])
                out <- cplexAPI::setIntParmCPLEX(lp@oobj@env,
                                                 cplexAPI::CPX_PARAM_SCRIND,
                                                 cplexAPI::CPX_OFF)
            },
            # ----------------------- #
            {
                warning("not a valid instance of class optObj")
                return(FALSE)
            }
        )

        return(lp)
    }
)


#------------------------------------------------------------------------------#

setMethod("backupProb", signature(lp = "optObj"),

    function(lp) {

        out <- FALSE
        np  <- FALSE

        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                # reset parameters!!! Parameters are reset to default when
                # doing initProbGLPK() !!!
                np <- glpkAPI::initProbGLPK()
                glpkAPI::copyProbGLPK(lp@oobj, np)
            },
            # ----------------------- #
            "clp" = {
                #fname <- paste("CLP_PROB_",
                #               format(Sys.time(), "%H%M%OS3"), sep = "")
                fname <- tempfile(pattern = "CLP_PROB_", fileext = ".tmp")
                ft <- clpAPI::saveModelCLP(lp@oobj, fname)
                if (ft != 0) {
                    stop("cannot save model")
                }
                np <- clpAPI::initProbCLP()
                clpAPI::setLogLevelCLP(np, 0)
                ft <- clpAPI::restoreModelCLP(np, fname)
                if (ft != 0) {
                    stop("cannot read model")
                }
                else {
                    unlink(fname)
                }
            },
            # ----------------------- #
            "lpSolveAPI" = {
                #fname <- paste("LPSOLVE_PROB_",
                #               format(Sys.time(), "%H%M%OS3"), sep = "")
                fname <- tempfile(pattern = "LPSOLVE_PROB_", fileext = ".tmp")
                lpSolveAPI::write.lp(lp@oobj, filename = fname, type = "lp")
                if (isTRUE(file.exists(fname))) {
                    np <- lpSolveAPI::read.lp(fname, type = "lp")
                    unlink(fname)
                }
                else {
                    stop("cannot read model")
                }
            },
            # ----------------------- #
            "cplex" = {
                np <- cplexAPI::cloneProbCPLEX(lp@oobj@env, lp@oobj@lp)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        # create new lp problem
        if (!identical(np, FALSE)) {
            out <- optObj(lp@solver, lp@method)
            if (lp@solver == "cplex") {
                out@oobj <- list(env = lp@oobj@env, lp = np)
            }
            else {
                out@oobj <- np
            }

        }

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("setSolverParm", signature(lp = "optObj"),

    function(lp, solverParm) {

        out <- FALSE
        if (!is(solverParm, "data.frame")) {
            warning(paste("Argument 'solverParm' must be of class",
                          "'data.frame'."
                    )
            )
            return(out)
        }

        if (any(is.na(solverParm))) {
            warning("Argument 'solverParm' contains 'NA' values.")
            return(out)
        }

        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                parm <- sapply(dimnames(solverParm)[[2]],
                               function(x) eval(parse(text = x)))
                val  <- solverParm[1,]
                if (lp@method == "interior") {
                    glpkAPI::setInteriorParmGLPK(parm, val)
                    out <- TRUE
                }
                else {
                    glpkAPI::setSimplexParmGLPK(parm, val)
                    out <- TRUE
                }
            },
            # ----------------------- #
            "clp" = {
            # no parameters in COIN-OR CLP yet.
            #    lp@oobj <- clpAPI::initProbCLP()
            #    clpAPI::setLogLevelCLP(lp@oobj, 0)
            },
            # ----------------------- #
            "lpSolveAPI" = {
                pname <- colnames(solverParm)
                for (i in seq(along = solverParm)) {
                    command <- paste("lpSolveAPI::lp.control(lp@oobj, ",
                                     pname[i], "='" , solverParm[[i]], "')", sep = "")
                    #print(command)
                    eval(parse(text = command))
                }
                #print(lp.control(lp@oobj))
                out <- TRUE
            },
            # ----------------------- #
            "cplex" = {
                 intdbl  <- sapply(solverParm, is.integer)
                 strparm <- sapply(solverParm, is.numeric)
                 int  <- solverParm[intdbl]
                 dbl  <- solverParm[intdbl == FALSE & strparm == TRUE]
                 char <- solverParm[strparm == FALSE]

                 if (length(int) > 0) {
                     intp <- sapply(dimnames(int)[[2]],
                                    function(x) eval(parse(text = x)))
                     intv <- int[1,]
                     for (i in seq(along = int)) {
                         out  <- cplexAPI::setIntParmCPLEX(lp@oobj@env,
                                                           intp[i], intv[i])
                     }
                 }

                 if (length(dbl) > 0) {
                     dblp <- sapply(dimnames(dbl)[[2]],
                                    function(x) eval(parse(text = x)))
                     dblv <- dbl[1,]
                     for (i in seq(along = dbl)) {
                         out  <- cplexAPI::setDblParmCPLEX(lp@oobj@env,
                                                           dblp[i], dblv[i])
                     }
                 }

                 if (length(char) > 0) {
                     charp <- sapply(dimnames(char)[[2]],
                                     function(x) eval(parse(text = x)))
                     charv <- char[1,]
                     for (i in seq(along = char)) {
                         out  <- cplexAPI::setStrParmCPLEX(lp@oobj@env,
                                                           charp[i], charv[i])
                     }
                 }

            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getSolverParm", signature(lp = "optObj"),

    function(lp) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                if (lp@method == "interior") {
                    out <- glpkAPI::getInteriorParmGLPK()
                }
                else {
                    out <- glpkAPI::getSimplexParmGLPK()
                }
            },
            # ----------------------- #
            "clp" = {
                wrong_solver_msg(lp, "getSolverParm")
                out <- FALSE
            },
            # ----------------------- #
            "lpSolveAPI" = {
                lpSolveAPI::lp.control(lp@oobj)
                out <- TRUE
            },
            # ----------------------- #
            "cplex" = {
                out <- cplexAPI::writeParmCPLEX(lp@oobj@env,
                                                "cplex_parameters.prm")
                message("Wrote the file 'cplex_parameters.prm'.")
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("setObjDir", signature(lp = "optObj"),

    function(lp, lpdir) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                dr <- ifelse(lpdir == "max", glpkAPI::GLP_MAX, glpkAPI::GLP_MIN)
                glpkAPI::setObjDirGLPK(lp@oobj, dr)
                out <- TRUE
            },
            # ----------------------- #
            "clp" = {
                dr <- ifelse(lpdir == "max", -1, 1)
                clpAPI::setObjDirCLP(lp@oobj, dr)
                out <- TRUE
            },
            # ----------------------- #
            "lpSolveAPI" = {
                lpSolveAPI::lp.control(lp@oobj, sense = lpdir)
                out <- TRUE
            },
            # ----------------------- #
            "cplex" = {
                dr <- ifelse(lpdir == "max",
                             cplexAPI::CPX_MAX,
                             cplexAPI::CPX_MIN)
                cplexAPI::setObjDirCPLEX(lp@oobj@env, lp@oobj@lp, dr)
                out <- TRUE
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getObjDir", signature(lp = "optObj"),

    function(lp) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                out <- glpkAPI::getObjDirGLPK(lp@oobj)
            },
            # ----------------------- #
            "clp" = {
                out <- clpAPI::getObjDirCLP(lp@oobj)
            },
            # ----------------------- #
            "lpSolveAPI" = {
                out <- lpSolveAPI::lp.control(lp@oobj)[["sense"]]
            },
            # ----------------------- #
            "cplex" = {
                out <- cplexAPI::getObjDirCPLEX(lp@oobj@env, lp@oobj@lp)
            },
             # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("addRows", signature(lp = "optObj"),

    function(lp, nrows) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                out <- glpkAPI::addRowsGLPK(lp@oobj, nrows)
            },
            # ----------------------- #
            "clp" = {
                # maybe we can do here something with resize()!
                wrong_solver_msg(lp, "addRows", printOut = TRUE)
            },
            # ----------------------- #
            "lpSolveAPI" = {
                ncols <- dim(lp@oobj)[2]
                out <- lpSolveAPI::resize.lp(lp@oobj, nrows, ncols)
            },
            # ----------------------- #
            "cplex" = {
                out <- cplexAPI::newRowsCPLEX(lp@oobj@env, lp@oobj@lp, nrows)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("addCols", signature(lp = "optObj"),

    function(lp, ncols) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                out <- glpkAPI::addColsGLPK(lp@oobj, ncols)
            },
            # ----------------------- #
            "clp" = {
                # maybe we can do here something with resize()!
                wrong_solver_msg(lp, "addCols", printOut = TRUE)
            },
            # ----------------------- #
            "lpSolveAPI" = {
                nrows <- dim(lp@oobj)[1]
                out <- lpSolveAPI::resize.lp(lp@oobj, nrows, ncols)
            },
            # ----------------------- #
            "cplex" = {
                out <- cplexAPI::newColsCPLEX(lp@oobj@env, lp@oobj@lp, ncols)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("addRowsCols", signature(lp = "optObj"),

    function(lp, nrows, ncols) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                outi <- glpkAPI::addRowsGLPK(lp@oobj, nrows)
                outj <- glpkAPI::addColsGLPK(lp@oobj, ncols)
                out  <- c(outi, outj)
            },
            # ----------------------- #
            "clp" = {
                clpAPI::resizeCLP(lp@oobj, nrows, ncols)
                out <- TRUE
            },
            # ----------------------- #
            "lpSolveAPI" = {
                out <- lpSolveAPI::resize.lp(lp@oobj, nrows, ncols)
            },
            # ----------------------- #
            "cplex" = {
                outi <- cplexAPI::newRowsCPLEX(lp@oobj@env, lp@oobj@lp, nrows)
                outj <- cplexAPI::newColsCPLEX(lp@oobj@env, lp@oobj@lp, ncols)
                out  <- c(outi, outj)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)

#------------------------------------------------------------------------------#

setMethod("getNumRows", signature(lp = "optObj"),

    function(lp) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                out <- glpkAPI::getNumRowsGLPK(lp@oobj)
            },
            # ----------------------- #
            "clp" = {
                out <- clpAPI::getNumRowsCLP(lp@oobj)
            },
            # ----------------------- #
            "lpSolveAPI" = {
                out <- dim(lp@oobj)[1]
            },
            # ----------------------- #
            "cplex" = {
                out <- cplexAPI::getNumRowsCPLEX(lp@oobj@env, lp@oobj@lp)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getNumCols", signature(lp = "optObj"),

    function(lp) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                out <- glpkAPI::getNumColsGLPK(lp@oobj)
            },
            # ----------------------- #
            "clp" = {
                out <- clpAPI::getNumColsCLP(lp@oobj)
            },
            # ----------------------- #
            "lpSolveAPI" = {
                out <- dim(lp@oobj)[2]
            },
            # ----------------------- #
            "cplex" = {
                out <- cplexAPI::getNumColsCPLEX(lp@oobj@env, lp@oobj@lp)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("addRowsToProb", signature(lp = "optObj"),

    # i: vector containing the new row indices (must be ascending)
    # cind: list, containing the column indices of the new nz elements
    # nzval: list, containing the new nz elements
    #
    # i, type, lb, cind and nzval must have the same length
    #
    # type can be one of the following:
    # "F" = free variable                -INF <  x <  INF
    # "L" = variable with lower bound      lb <= x <  INF
    # "U" = variable with upper bound    -INF <  x <= ub
    # "D" = double-bounded variable        lb <= x <= ub
    # "E" = fixed variable                 lb  = x  = ub
    # "R" = ranged constraint

    function(lp, i, type, lb, ub, cind, nzval) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                ord <- glpkAPI::addRowsGLPK(lp@oobj, length(i))
                gtype = integer(length(type))
                for (k in seq(along = i)) {
                    gtype[k] <- switch(type[k],
                                       "F" = { glpkAPI::GLP_FR },
                                       "L" = { glpkAPI::GLP_LO },
                                       "U" = { glpkAPI::GLP_UP },
                                       "D" = { glpkAPI::GLP_DB },
                                       "E" = { glpkAPI::GLP_FX },
                                             { glpkAPI::GLP_FX }
                    )
                    glpkAPI::setMatRowGLPK(lp@oobj, i[k],
                                           length(cind[[k]]),
                                           cind[[k]], nzval[[k]])
                }
                out <- glpkAPI::setRowsBndsGLPK(lp@oobj, i, lb, ub, gtype)
            },
            # ----------------------- #
            "clp" = {
                cst <- c(0, cumsum(unlist(lapply(cind, length))))
                out <- clpAPI::addRowsCLP(lp@oobj, length(i), lb, ub,
                                          cst, unlist(cind)-1, unlist(nzval))
            },
            # ----------------------- #
            "lpSolveAPI" = {
                for (k in seq(along = i)) {
                    ltype <- switch(type[k],
                                    "L" = { 2 },
                                    "U" = { 1 },
                                    "E" = { 3 },
                                          { 3 }
                    )
                    out <- lpSolveAPI::add.constraint(lp@oobj, nzval[[k]],
                                                      ltype, lb[k], cind[[k]])
                }
            },
            # ----------------------- #
            "cplex" = {
                cptype = character(length(type))
                for (l in seq(along = type)) {
                    cptype[l] <- switch(type[l],
                        "L" = { "G" },
                        "U" = { "L" },
                        "E" = { "E" },
                        "R" = { "R" },
                              { "E" }
                    )
                }
                beg <- c(0, cumsum(unlist(lapply(cind, length))))
                out <- cplexAPI::addRowsCPLEX(lp@oobj@env, lp@oobj@lp, 0,
                                              length(i), length(nzval), beg,
                                              unlist(cind)-1, unlist(nzval),
                                              lb, cptype)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeColsBnds", signature(lp = "optObj"),

    function(lp, j, lb, ub) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                glpkAPI::setColsBndsGLPK(lp@oobj, j, lb, ub)
                out <- TRUE
            },
            # ----------------------- #
            "clp" = {
                tmp_lb <- clpAPI::getColLowerCLP(lp@oobj)
                tmp_ub <- clpAPI::getColUpperCLP(lp@oobj)
                tmp_lb[j] <- lb
                tmp_ub[j] <- ub
                clpAPI::chgColLowerCLP(lp@oobj, tmp_lb)
                clpAPI::chgColUpperCLP(lp@oobj, tmp_ub)
                out <- TRUE
            },
            # ----------------------- #
            "lpSolveAPI" = {
                out <- lpSolveAPI::set.bounds(lp@oobj, lb, ub, j)
            },
            # ----------------------- #
            "cplex" = {
                out <- cplexAPI::chgColsBndsCPLEX(lp@oobj@env,
                                                  lp@oobj@lp, j-1, lb, ub)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeColsBndsObjCoefs", signature(lp = "optObj"),

    function(lp, j, lb, ub, obj_coef) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                glpkAPI::setColsBndsObjCoefsGLPK(lp@oobj, j, lb, ub, obj_coef)
                out <- TRUE
            },
            # ----------------------- #
            "clp" = {
                # usable only for model creation!
                clpAPI::chgColLowerCLP(lp@oobj, lb)
                clpAPI::chgColUpperCLP(lp@oobj, ub)
                clpAPI::chgObjCoefsCLP(lp@oobj, obj_coef)
                out <- TRUE
            },
            # ----------------------- #
            "lpSolveAPI" = {
                outb <- lpSolveAPI::set.bounds(lp@oobj, lb, ub, j)
                outo <- lpSolveAPI::set.objfn(lp@oobj, obj_coef, j)
                out  <- c(outb, outo)
            },
            # ----------------------- #
            "cplex" = {
                outb <- cplexAPI::chgColsBndsCPLEX(lp@oobj@env,
                                                   lp@oobj@lp, j-1, lb, ub)
                outo <- cplexAPI::chgObjCPLEX(lp@oobj@env, lp@oobj@lp,
                                              length(j), j-1, obj_coef)
                out  <- c(outb, outo)
                # usable only for model creation!
                # out <- cplexAPI::newColsCPLEX(lp@oobj@env, lp@oobj@lp,
                #                               length(j), obj_coef, lb, ub)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getColsLowBnds", signature(lp = "optObj"),

    function(lp, j) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                out <- glpkAPI::getColsLowBndsGLPK(lp@oobj, j)
            },
            # ----------------------- #
            "clp" = {
                out <- clpAPI::getColLowerCLP(lp@oobj)[j]
            },
            # ----------------------- #
            "lpSolveAPI" = {
                out <- lpSolveAPI::get.bounds(lp@oobj, j)[["lower"]]
            },
            # ----------------------- #
            "cplex" = {
                out <- cplexAPI::getLowBndsIdsCPLEX(lp@oobj@env,
                                                    lp@oobj@lp, j-1)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getColsUppBnds", signature(lp = "optObj"),

    function(lp, j) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                out <- glpkAPI::getColsUppBndsGLPK(lp@oobj, j)
            },
            # ----------------------- #
            "clp" = {
                out <- clpAPI::getColUpperCLP(lp@oobj)[j]
            },
            # ----------------------- #
            "lpSolveAPI" = {
                out <- lpSolveAPI::get.bounds(lp@oobj, j)[["upper"]]
            },
            # ----------------------- #
            "cplex" = {
                out <- cplexAPI::getUppBndsIdsCPLEX(lp@oobj@env,
                                                    lp@oobj@lp, j-1)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeRowsBnds", signature(lp = "optObj"),

    function(lp, i, lb, ub) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                glpkAPI::setRowsBndsGLPK(lp@oobj, i, lb, ub)
                out <- TRUE
            },
            # ----------------------- #
            "clp" = {
                tmp_lb <- clpAPI::getRowLowerCLP(lp@oobj)
                tmp_ub <- clpAPI::getRowUpperCLP(lp@oobj)
                tmp_lb[i] <- lb
                tmp_ub[i] <- ub
                clpAPI::chgRowLowerCLP(lp@oobj, tmp_lb)
                clpAPI::chgRowUpperCLP(lp@oobj, tmp_ub)
                out <- TRUE
            },
            # ----------------------- #
            "lpSolveAPI" = {
                out <- lpSolveAPI::set.constr.value(lp@oobj,
                                                    rhs = ub,
                                                    lhs = lb,
                                                    constraints = i)
            },
            # ----------------------- #
            "cplex" = {
                wrong_solver_msg(lp, "changeRowsBnds", printOut = TRUE)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("setRhsZero", signature(lp = "optObj"),

    function(lp) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                glpkAPI::setRhsZeroGLPK(lp@oobj)
                out <- TRUE
            },
            # ----------------------- #
            "clp" = {
                nrows <- clpAPI::getNumRowsCLP(lp@oobj)
                zeros <- rep(0, nrows)
                clpAPI::chgRowLowerCLP(lp@oobj, zeros)
                clpAPI::chgRowUpperCLP(lp@oobj, zeros)
                out <- TRUE
            },
            # ----------------------- #
            "lpSolveAPI" = {
                nrows <- dim(lp@oobj)[1]
                outb <- lpSolveAPI::set.constr.value(lp@oobj,
                                                     rhs = rep(0, nrows),
                                                     lhs = NULL,
                                                     constraints = c(1:nrows))
                outt <- lpSolveAPI::set.constr.type(lp@oobj,
                                                    rep(3, nrows), c(1:nrows))
                out <- c(outb, outt)
            },
            # ----------------------- #
            "cplex" = {
                nrows  <- cplexAPI::getNumRowsCPLEX(lp@oobj@env, lp@oobj@lp)
                zeros  <- rep(0, nrows)
                indic  <- c(0:(nrows-1))
                outb   <- cplexAPI::chgRhsCPLEX(lp@oobj@env, lp@oobj@lp,
                                                nrows, indic, zeros)
                outt   <- cplexAPI::chgSenseCPLEX(lp@oobj@env, lp@oobj@lp,
                                                  nrows, indic, rep("E", nrows))
                out <- c(outb, outt)
                # usable only for model creation!
                # ( Variable nrows has to be argument of setRhsZero()! )
                # out <- cplexAPI::newRowsCPLEX(lp@oobj@env, lp@oobj@lp, nrows,
                #                               rep(0, nrows), rep("E", nrows))
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getRowsLowBnds", signature(lp = "optObj"),

    function(lp, i) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                out <- glpkAPI::getRowsLowBndsGLPK(lp@oobj, i)
            },
            # ----------------------- #
            "clp" = {
                out <- clpAPI::getRowLowerCLP(lp@oobj)[i]
            },
            # ----------------------- #
            "lpSolveAPI" = {
                out <- lpSolveAPI::get.constr.value(lp@oobj,
                                                    side = "lhs",
                                                    constraints = i)
            },
            # ----------------------- #
            "cplex" = {
                wrong_solver_msg(lp, "getRowsLowBnds", printOut = TRUE)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getRowsUppBnds", signature(lp = "optObj"),

    function(lp, i) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                out <- glpkAPI::getRowsUppBndsGLPK(lp@oobj, i)
            },
            # ----------------------- #
            "clp" = {
                out <- clpAPI::getRowUpperCLP(lp@oobj)[i]
            },
            # ----------------------- #
            "lpSolveAPI" = {
                out <- lpSolveAPI::get.constr.value(lp@oobj,
                                                    side = "rhs",
                                                    constraints = i)
            },
            # ----------------------- #
            "cplex" = {
                wrong_solver_msg(lp, "getRowsUppBnds", printOut = TRUE)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeObjCoefs", signature(lp = "optObj"),

    function(lp, j, obj_coef) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                glpkAPI::setObjCoefsGLPK(lp@oobj, j, obj_coef)
                out <- TRUE
            },
            # ----------------------- #
            "clp" = {
                tmp_obj_coef <- clpAPI::getObjCoefsCLP(lp@oobj)
                tmp_obj_coef[j] <- obj_coef
                clpAPI::chgObjCoefsCLP(lp@oobj, tmp_obj_coef)
                out <- TRUE
            },
            # ----------------------- #
            "lpSolveAPI" = {
                out <- lpSolveAPI::set.objfn(lp@oobj, obj_coef, j)
            },
            # ----------------------- #
            "cplex" = {
                out <- cplexAPI::chgObjCPLEX(lp@oobj@env, lp@oobj@lp,
                                             length(j), j-1, obj_coef)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getObjCoefs", signature(lp = "optObj"),

    function(lp, j) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                out <- glpkAPI::getObjCoefsGLPK(lp@oobj, j)
            },
            # ----------------------- #
            "clp" = {
                out <- clpAPI::getObjCoefsCLP(lp@oobj)[j]
            },
            # ----------------------- #
            "lpSolveAPI" = {
                #wrong_solver_msg(lp, "getObjCoefs", printOut = TRUE)
                out <- numeric(length(j))
                for (i in seq(along = j)) {
                    out[i] <- lpSolveAPI::get.column(lp@oobj, j[i])$column[1]
                }
            },
            # ----------------------- #
            "cplex" = {
                if (length(j) > 1) {
                    b <- min(j) - 1
                    e <- max(j) - 1
                }
                else {
                    b <- j - 1
                    e <- j - 1
                }
                out <- cplexAPI::getObjCPLEX(lp@oobj@env, lp@oobj@lp, b, e)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("loadMatrix", signature(lp = "optObj"),

    function(lp, ...) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                glpkAPI::loadMatrixGLPK(lp@oobj, ...)
                out <- TRUE
            },
            # ----------------------- #
            "clp" = {
                clpAPI::loadMatrixCLP(lp@oobj, ...)
                out <- TRUE
            },
            # ----------------------- #
            "lpSolveAPI" = {
                out <- loadMatrixPerColumnLPSOLVE(lp@oobj, ...)
            },
            # ----------------------- #
            "cplex" = {
                out <- cplexAPI::chgCoefListCPLEX(lp@oobj@env, lp@oobj@lp, ...)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("loadProblemData", signature(lp = "optObj"),

    function(lp, model) {

        if (!is(model, "modelorg")) {
            stop("needs an object of class modelorg!")
        }

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                out <- vector(mode = "list", length = 4)
                out[[1]] <- addRowsCols(lp, met_num(model), react_num(model))
                TMPmat <- as(S(model), "TsparseMatrix")
                out[[2]] <- loadMatrix(lp,
                                       length(TMPmat@x),
                                       TMPmat@i + 1,
                                       TMPmat@j + 1,
                                       TMPmat@x)

                # add upper and lower bounds: ai < vi < bi
                out[[3]] <- changeColsBndsObjCoefs(lp,
                                                   c(1:react_num(model)),
                                                   lowbnd(model),
                                                   uppbnd(model),
                                                   obj_coef(model))

                # set the right hand side Sv = 0
                out[[4]] <- setRhsZero(lp)
            },
            # ----------------------- #
            "clp" = {
                zeros    <- rep(0, met_num(model))
                TMPmat <- as(S(model), "CsparseMatrix")
                check    <- clpAPI::loadProblemCLP(lp@oobj,
                                                   react_num(model),
                                                   met_num(model),
                                                   TMPmat@i,
                                                   TMPmat@p,
                                                   TMPmat@x,
                                                   lowbnd(model),
                                                   uppbnd(model),
                                                   obj_coef(model),
                                                   zeros,
                                                   zeros)
                out <- TRUE
            },
            # ----------------------- #
            "lpSolveAPI" = {
                out <- vector(mode = "list", length = 3)
                # add rows and columns to the problem
                out[[1]] <- loadMatrix(lp, as(S(model), "CsparseMatrix"))
                out[[2]] <- changeColsBndsObjCoefs(lp,
                                                   c(1:react_num(model)),
                                                   lowbnd(model),
                                                   uppbnd(model),
                                                   obj_coef(model))
                out[[3]] <- setRhsZero(lp)
            },
            # ----------------------- #
            "cplex" = {
                out <- vector(mode = "list", length = 3)
                TMPmat <- as(S(model), "TsparseMatrix")

                nrows <- met_num(model)

                out[[1]] <- cplexAPI::newRowsCPLEX(lp@oobj@env, lp@oobj@lp,
                                                   nrows, rep(0, nrows),
                                                   rep("E", nrows))

                out[[2]] <- cplexAPI::newColsCPLEX(lp@oobj@env, lp@oobj@lp,
                                                   react_num(model),
                                                   obj_coef(model),
                                                   lowbnd(model), uppbnd(model))

                out[[3]] <- cplexAPI::chgCoefListCPLEX(lp@oobj@env, lp@oobj@lp,
                                                       length(TMPmat@x),
                                                       TMPmat@i,
                                                       TMPmat@j,
                                                       TMPmat@x)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("loadProblemDataLM", signature(lp = "optObj"),

    function(lp, model, wtflux, nCols, nRows, COBRAflag = FALSE, lpdir = NA) {

        if (!is(model, "modelorg")) {
            stop("needs an object of class modelorg!")
        }

        out <- FALSE

        nc     <- react_num(model)
        nr     <- met_num(model)
        #lpdir  <- getObjDir(lp)
        absMAX <- SYBIL_SETTINGS("MAXIMUM")


        #  the problem: minimize
        #
        #            |      |      |                ]
        #        Swt |  0   |  0   |  = 0           ] left out if not COBRA
        #            |      |      |                ]
        #       -------------------------
        #            |      |      |
        #         0  | Sdel |  0   |  = 0
        #            |      |      |
        #       -------------------------
        #            |      |      |
        #         v1 - v2   |delta-| >= 0
        #         v2 - v1   |delta+| >= 0
        #            |      |      |
        #       -------------------------
        #       c_wt |  0   |  0   | >= c^T * v_wt  ] left out if not COBRA
        #            |      |      |
        #  lb   wt_lb|del_lb|  0   |                ] COBRA version
        #  ub   wt_ub|del_ub|10000 |                ] COBRA version
        #            |      |      |
        #  lb   v_wt |del_lb|  0   |
        #  ub   v_wt |del_ub|10000 |
        #            |      |      |
        #            |      |      |
        #  obj    0  |  0   |  1   |


        # ---------------------------------------------
        # constraint matrix
        # ---------------------------------------------

        # the initial matrix dimensions
     	LHS <- Matrix::Matrix(0, nrow = nr+2*nc, ncol = 4*nc, sparse = TRUE)

        # rows for the wild type strain
        LHS[1:nr,(nc+1):(2*nc)] <- S(model)

        # rows for the delta match matrix
		diag(LHS[(nr+1)   :(nr+2*nc),1       :(2*nc)]) <- -1
		diag(LHS[(nr+1)   :(nr+2*nc),(2*nc+1):(4*nc)]) <- 1
        diag(LHS[(nr+1)   :(nr+nc)  ,(nc+1)  :(2*nc)]) <- 1
        diag(LHS[(nr+nc+1):(nr+2*nc),1       :nc    ]) <- 1


        # contraint matrix for linearMOMA, COBRA version
        if (isTRUE(COBRAflag)) {

            # rows for the wild type strain
            LHSwt <- Matrix::Matrix(0, nrow = nr, ncol = 4*nc, sparse = TRUE)
            LHSwt[1:nr,1:nc] <- S(model)

            # fix the value of the objective function
            crow <- Matrix::Matrix(c(obj_coef(model), rep(0, 3*nc)),
                                   nrow = 1, ncol = 4*nc, sparse = TRUE)

            # the final contraint matrix
            LHS <- rBind(LHSwt, LHS, crow)

        }


        # ---------------------------------------------
        # lower and upper bounds
        # ---------------------------------------------

        if (isTRUE(COBRAflag)) {
            # Here we calculate wild type and deletion strain simultaineously,
            # so we need upper and lower bounds for both, the wild type and the
            # deletion strain.
            # All the delta's are positive.
            lower <- c(lowbnd(model), lowbnd(model), rep(0, 2*nc))
            upper <- c(uppbnd(model), uppbnd(model), rep(absMAX, 2*nc))


            rlower <- c(rep(0, nRows-1), wtflux)
            rupper <- c(rep(0, 2*nr), rep(absMAX, 2*nc), wtflux)
            #rupper <- rlower
            #rupper <- c(rep(0, 2*nr), rep(0, 2*nc), 0)
        }
        else {
            # Here, we keep the wild type flux distribution fixed.
            lower  <- c(wtflux, lowbnd(model), rep(0, 2*nc))
            upper  <- c(wtflux, uppbnd(model), rep(absMAX, 2*nc))
            rlower <- c(rep(0, nRows))
            rupper <- c(rep(0, nr), rep(absMAX, 2*nc))
        }


        # ---------------------------------------------
        # objective function
        # ---------------------------------------------

        cobj <- c(rep(0, 2*nc), rep(1, 2*nc))


        # ---------------------------------------------
        # build problem object
        # ---------------------------------------------

        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                out <- vector(mode = "list", length = 4)

                if (isTRUE(COBRAflag)) {
                    rtype <- c(rep(glpkAPI::GLP_FX, 2*nr),
                               rep(glpkAPI::GLP_LO, 2*nc))
                    if (lpdir == "max") {
                        rtype <- c(rtype, glpkAPI::GLP_LO)
                    }
                    else {
                        rtype <- c(rtype, glpkAPI::GLP_UP)
                    }
                    #rtype <- c(rep(glpkAPI::GLP_FX, 2*nr),
                    #           rep(glpkAPI::GLP_DB, 2*nc),
                    #           glpkAPI::GLP_DB)
                }
                else {
                    rtype  <- c(rep(glpkAPI::GLP_FX, nr),
                                rep(glpkAPI::GLP_LO, 2*nc))
                }

                TMPmat <- as(LHS, "TsparseMatrix")

                out[[1]] <- addRowsCols(lp, nRows, nCols)
                out[[2]] <- loadMatrix(lp,
                                       length(TMPmat@x),
                                       TMPmat@i + 1,
                                       TMPmat@j + 1,
                                       TMPmat@x)

                # add upper and lower bounds: ai < vi < bi
                out[[3]] <- changeColsBndsObjCoefs(lp,
                                                   c(1:nCols),
                                                   lower,
                                                   upper,
                                                   cobj)

                # set the right hand side Sv = b
                out[[4]] <- glpkAPI::setRowsBndsGLPK(lp@oobj,
                                                     c(1:nRows),
                                                     rlower,
                                                     rupper,
                                                     rtype)

            },
            # ----------------------- #
            "clp" = {

                TMPmat <- as(LHS, "CsparseMatrix")
                out   <- clpAPI::loadProblemCLP(lp@oobj,
                                                nCols,
                                                nRows,
                                                TMPmat@i,
                                                TMPmat@p,
                                                TMPmat@x,
                                                lower,
                                                upper,
                                                cobj,
                                                rlower,
                                                rupper)
                #out <- TRUE
            },
            # ----------------------- #
            "lpSolveAPI" = {
                out <- vector(mode = "list", length = 4)

                if (isTRUE(COBRAflag)) {
                    rtype <- c(rep(3, 2*nr), rep(2, 2*nc))
                    if (lpdir == "max") {
                        rtype <- c(rtype, 2)
                    }
                    else {
                        rtype <- c(rtype, 1)
                    }
                }
                else {
                    rtype  <- c(rep(3, nr), rep(2, 2*nc))
                }

                out[[1]] <- loadMatrix(lp, LHS)
                out[[2]] <- changeColsBndsObjCoefs(lp,
                                                   c(1:nCols),
                                                   lower,
                                                   upper,
                                                   cobj)
                out[[3]] <- lpSolveAPI::set.constr.value(lprec = lp@oobj,
                                                       rhs = rlower,
                                                       constraints = c(1:nRows))

                out[[3]] <- lpSolveAPI::set.constr.type(lp@oobj,
                                                        rtype,
                                                        c(1:nRows))

            },
            # ----------------------- #
            "cplex" = {
                out <- vector(mode = "list", length = 3)

                if (isTRUE(COBRAflag)) {
                    rtype <- c(rep("E", 2*nr), rep("G", 2*nc))
                    if (lpdir == "max") {
                        rtype <- c(rtype, "G")
                    }
                    else {
                        rtype <- c(rtype, "L")
                    }
                }
                else {
                    rtype  <- c(rep("E", nr), rep("G", 2*nc))
                }

                TMPmat <- as(LHS, "TsparseMatrix")

                out[[1]] <- cplexAPI::newRowsCPLEX(lp@oobj@env, lp@oobj@lp,
                                                   nRows, rlower, rtype)

                out[[2]] <- cplexAPI::newColsCPLEX(lp@oobj@env, lp@oobj@lp,
                                                   nCols, cobj, lower, upper)

                out[[3]] <- cplexAPI::chgCoefListCPLEX(lp@oobj@env, lp@oobj@lp,
                                                       length(TMPmat@x),
                                                       TMPmat@i,
                                                       TMPmat@j,
                                                       TMPmat@x)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("loadProblemDataMTF", signature(lp = "optObj"),

    function(lp, model, wtflux, nCols, nRows) {

        if (!is(model, "modelorg")) {
            stop("needs an object of class modelorg!")
        }

        out <- FALSE

        nc     <- react_num(model)
        nr     <- met_num(model)
        absMAX <- SYBIL_SETTINGS("MAXIMUM")


        #  the problem: minimize:
        #
        #            |      |      |
        #         S  |  0   |  0   |  = b
        #            |      |      |
        #       -------------------------
        #            |      |      |
        #         1  |  1   |  0   | >= 0
        #            |      |      |
        #       -------------------------
        #            |      |      |
        #         -1 |  0   |  1   | >= 0
        #            |      |      |
        #       -------------------------
        #       c_wt |  0   |  0   | >= c^T * v_wt
        #            |      |      |
        #  lb   wt_lb|  0   |  0   |
        #  ub   wt_ub|10000 |10000 |
        #            |      |      |
        #  obj    0  |  1   |  1   |


        # ---------------------------------------------
        # constraint matrix
        # ---------------------------------------------

        # the initial matrix dimensions
        LHS <- Matrix::Matrix(0, nrow = nr+2*nc+1, ncol = 3*nc, sparse = TRUE)

        # rows for the mutant strain
        LHS[1:nr,1:nc] <- S(model)

        # rows for the delta match matrix
		diag(LHS[(nr+1)   :(nr+nc)  ,1       :nc    ]) <- 1
		diag(LHS[(nr+1)   :(nr+nc)  ,(nc+1)  :(2*nc)]) <- 1
        diag(LHS[(nr+nc+1):(nr+2*nc),1       :nc    ]) <- -1
        diag(LHS[(nr+nc+1):(nr+2*nc),(2*nc+1):(3*nc)]) <- 1

	    # fix the value of the objective function
        LHS[(nr+2*nc+1),1:nc] <- obj_coef(model)


        # ---------------------------------------------
        # lower and upper bounds
        # ---------------------------------------------

		lower  <- c(lowbnd(model), rep(0, 2*nc))
		upper  <- c(uppbnd(model), rep(absMAX, 2*nc))
		rlower <- c(rhs(model), rep(0, 2*nc), wtflux)
		rupper <- c(rhs(model), rep(absMAX, 2*nc + 1))


        # ---------------------------------------------
        # objective function
        # ---------------------------------------------

        cobj <- c(rep(0, nc), rep(1, 2*nc))


        # ---------------------------------------------
        # build problem object
        # ---------------------------------------------

        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                out <- vector(mode = "list", length = 4)

                rtype <- c(rep(glpkAPI::GLP_FX, nr),
                           rep(glpkAPI::GLP_LO, 2*nc + 1))

                TMPmat <- as(LHS, "TsparseMatrix")

                out[[1]] <- addRowsCols(lp, nRows, nCols)
                out[[2]] <- loadMatrix(lp,
                                       length(TMPmat@x),
                                       TMPmat@i + 1,
                                       TMPmat@j + 1,
                                       TMPmat@x)

                # add upper and lower bounds: ai < vi < bi
                out[[3]] <- changeColsBndsObjCoefs(lp,
                                                   c(1:nCols),
                                                   lower,
                                                   upper,
                                                   cobj)

                # set the right hand side Sv = b
                out[[4]] <- glpkAPI::setRowsBndsGLPK(lp@oobj,
                                                     c(1:nRows),
                                                     rlower,
                                                     rupper,
                                                     rtype)
            },
            # ----------------------- #
            "clp" = {

                TMPmat <- as(LHS, "CsparseMatrix")
                out   <- clpAPI::loadProblemCLP(lp@oobj,
                                                nCols,
                                                nRows,
                                                TMPmat@i,
                                                TMPmat@p,
                                                TMPmat@x,
                                                lower,
                                                upper,
                                                cobj,
                                                rlower,
                                                rupper)
                #out <- TRUE
            },
            # ----------------------- #
            "lpSolveAPI" = {
                out <- vector(mode = "list", length = 4)

                rtype <- c(rep(3, nr), rep(2, 2*nc + 1))

                out[[1]] <- loadMatrix(lp, LHS)
                out[[2]] <- changeColsBndsObjCoefs(lp,
                                                   c(1:nCols),
                                                   lower,
                                                   upper,
                                                   cobj)
                out[[3]] <- lpSolveAPI::set.constr.value(lprec = lp@oobj,
                                                       rhs = rlower,
                                                       constraints = c(1:nRows))

                out[[3]] <- lpSolveAPI::set.constr.type(lp@oobj,
                                                        rtype,
                                                        c(1:nRows))

            },
            # ----------------------- #
            "cplex" = {
                out <- vector(mode = "list", length = 3)

                rtype <- c(rep("E", nr), rep("G", 2*nc + 1))

                TMPmat <- as(LHS, "TsparseMatrix")

                out[[1]] <- cplexAPI::newRowsCPLEX(lp@oobj@env, lp@oobj@lp,
                                                   nRows, rlower, rtype)

                out[[2]] <- cplexAPI::newColsCPLEX(lp@oobj@env, lp@oobj@lp,
                                                   nCols, cobj, lower, upper)

                out[[3]] <- cplexAPI::chgCoefListCPLEX(lp@oobj@env, lp@oobj@lp,
                                                       length(TMPmat@x),
                                                       TMPmat@i,
                                                       TMPmat@j,
                                                       TMPmat@x)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("scaleProb", signature(lp = "optObj"),

    function(lp, opt) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                # check if tryCatch works here!!
                glpkAPI::scaleProbGLPK(lp@oobj, opt)
                out <- TRUE
            },
            # ----------------------- #
            "clp" = {
                clpAPI::scaleModelCLP(lp@oobj, opt)
                out <- TRUE
            },
            # ----------------------- #
            "lpSolveAPI" = {
                invisible(lpSolveAPI::lp.control(lp@oobj, scaling = opt))
                out <- TRUE
            },
            # ----------------------- #
            "cplex" = {
                out <- cplexAPI::setIntParmCPLEX(lp@oobj@env,
                                                 cplexAPI::CPX_PARAM_REDUCE,
                                                 opt)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("solveLp", signature(lp = "optObj"),

    function(lp) {

        out          <- FALSE
        method_valid <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
#                 if (glpkAPI::bfExistsGLPK(lp@oobj) != 0) {
#                     if (glpkAPI::bfUpdatedGLPK(lp@oobj) != 0) {
#                         basis <- glpkAPI::factorizeGLPK(lp@oobj)
#                         #print(basis)
#                     }
#                 }
                switch(lp@method,
                    "interior" = {
                        out <- glpkAPI::solveInteriorGLPK(lp@oobj)
                        method_valid <- TRUE
                    },
                    "exact" = {
                        out <- glpkAPI::solveSimplexExactGLPK(lp@oobj)
                        method_valid <- TRUE
                    },
                    {
                        out <- glpkAPI::solveSimplexGLPK(lp@oobj)
                    }
                )
            },
            # ----------------------- #
            "clp" = {
                switch(lp@method,
                    "inidual" = {
                        out <- clpAPI::solveInitialDualCLP(lp@oobj)
                        method_valid <- TRUE
                    },
                    "iniprimal" = {
                        out <- clpAPI::solveInitialPrimalCLP(lp@oobj)
                        method_valid <- TRUE
                    },
                    "inibarrier" = {
                        out <- clpAPI::solveInitialBarrierCLP(lp@oobj)
                        method_valid <- TRUE
                    },
                    "inibarriernoc" = {
                        out <- clpAPI::solveInitialBarrierNoCrossCLP(lp@oobj)
                        method_valid <- TRUE
                    },
                    "dual" = {
                        out <- clpAPI::dualCLP(lp@oobj)
                        method_valid <- TRUE
                    },
                    "primal" = {
                        out <- clpAPI::primalCLP(lp@oobj)
                        method_valid <- TRUE
                    },
                    "idiot" = {
                        clpAPI::idiotCLP(lp@oobj) # idiotCLP has no return value
                        out <- 0
                        method_valid <- TRUE
                    },
                    {
                        out <- clpAPI::solveInitialCLP(lp@oobj)
                    }
                )
            },
            # ----------------------- #
            "lpSolveAPI" = {
                out <- solve(lp@oobj)
            },
            # ----------------------- #
            "cplex" = {
                switch(lp@method,
                    "primopt" = {
                        out <- cplexAPI::primoptCPLEX(lp@oobj@env, lp@oobj@lp)
                        method_valid <- TRUE
                    },
                    "dualopt" = {
                        out <- cplexAPI::dualoptCPLEX(lp@oobj@env, lp@oobj@lp)
                        method_valid <- TRUE
                    },
                    "baropt" = {
                        out <- cplexAPI::baroptCPLEX(lp@oobj@env, lp@oobj@lp)
                        method_valid <- TRUE
                    },
                    "hybbaropt" = {
                        out <- cplexAPI::hybbaroptCPLEX(lp@oobj@env, lp@oobj@lp,
                                                        method = 0)
                        method_valid <- TRUE
                    },
                    "hybnetopt" = {
                        out <- cplexAPI::hybnetoptCPLEX(lp@oobj@env, lp@oobj@lp,
                                              method = cplexAPI::CPX_ALG_PRIMAL)
                        method_valid <- TRUE
                    },
                    "siftopt" = {
                        out <- cplexAPI::siftoptCPLEX(lp@oobj@env, lp@oobj@lp)
                        method_valid <- TRUE
                    },
                    {
                        out <- cplexAPI::lpoptCPLEX(lp@oobj@env, lp@oobj@lp)
                    }
                )
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getObjVal", signature(lp = "optObj"),

    function(lp) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                if (lp@method == "interior") {
                    out <- glpkAPI::getObjValIptGLPK(lp@oobj)
                }
                else {
                    out <- glpkAPI::getObjValGLPK(lp@oobj)
                }
            },
            # ----------------------- #
            "clp" = {
                out <- clpAPI::getObjValCLP(lp@oobj)
            },
            # ----------------------- #
            "lpSolveAPI" = {
                out <- lpSolveAPI::get.objective(lp@oobj)
            },
            # ----------------------- #
            "cplex" = {
                out <- cplexAPI::getObjValCPLEX(lp@oobj@env, lp@oobj@lp)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getRedCosts", signature(lp = "optObj"),

    function(lp) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                if (lp@method == "interior") {
                    out <- glpkAPI::getColsDualIptGLPK(lp@oobj)
                }
                else {
                    out <- glpkAPI::getColsDualGLPK(lp@oobj)
                }
            },
            # ----------------------- #
            "clp" = {
                out <- clpAPI::getColDualCLP(lp@oobj)
            },
            # ----------------------- #
            "lpSolveAPI" = {
                out <- lpSolveAPI::get.dual.solution(lp@oobj)
            },
            # ----------------------- #
            "cplex" = {
                nr  <- cplexAPI::getNumRowsCPLEX(lp@oobj@env, lp@oobj@lp)
                out <- cplexAPI::getDjCPLEX(lp@oobj@env, lp@oobj@lp, 0, nr-1)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getSolStat", signature(lp = "optObj"),

    function(lp) {

        out <- NA
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                if (lp@method == "interior") {
                    out <- glpkAPI::getSolStatIptGLPK(lp@oobj)
                }
                else {
                    out <- glpkAPI::getSolStatGLPK(lp@oobj)
                }
            },
            # ----------------------- #
            "clp" = {
                out <- clpAPI::getSolStatusCLP(lp@oobj)
            },
            # ----------------------- #
            "lpSolveAPI" = {
                wrong_solver_msg(lp, "getSolStat", printOut = FALSE)
            },
            # ----------------------- #
            "cplex" = {
                out <- cplexAPI::getStatCPLEX(lp@oobj@env, lp@oobj@lp)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getFluxDist", signature(lp = "optObj"),

    function(lp) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                if (lp@method == "interior") {
                    out <- glpkAPI::getColsPrimIptGLPK(lp@oobj)
                }
                else {
                    out <- glpkAPI::getColsPrimGLPK(lp@oobj)
                }
            },
            # ----------------------- #
            "clp" = {
                out <- clpAPI::getColPrimCLP(lp@oobj)
            },
            # ----------------------- #
            "lpSolveAPI" = {
                out <- lpSolveAPI::get.variables(lp@oobj)
            },
            # ----------------------- #
            "cplex" = {
                ncols <- cplexAPI::getNumColsCPLEX(lp@oobj@env, lp@oobj@lp)
                out   <- cplexAPI::getProbVarCPLEX(lp@oobj@env, lp@oobj@lp,
                                                   0, ncols - 1)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getColPrim", signature(lp = "optObj"),

    function(lp, j) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                if (lp@method == "interior") {
                    out <- glpkAPI::getColPrimIptGLPK(lp@oobj, j)
                }
                else {
                    out <- glpkAPI::getColPrimGLPK(lp@oobj, j)
                }
            },
            # ----------------------- #
            "clp" = {
                out <- clpAPI::getColPrimCLP(lp@oobj)[j]
            },
            # ----------------------- #
            "lpSolveAPI" = {
                out <- lpSolveAPI::get.variables(lp@oobj)[j]
            },
            # ----------------------- #
            "cplex" = {
                out <- cplexAPI::getProbVarCPLEX(lp@oobj@env, lp@oobj@lp, j, j)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getNumNnz", signature(lp = "optObj"),

    function(lp) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                out <- glpkAPI::getNumNnzGLPK(lp@oobj)
            },
            # ----------------------- #
            "clp" = {
                out <- clpAPI::getNumNnzCLP(lp@oobj)
            },
            # ----------------------- #
            "lpSolveAPI" = {
                wrong_solver_msg(lp, "getNumNnz", printOut = TRUE)
            },
            # ----------------------- #
            "cplex" = {
                out <- cplexAPI::getNumNnzCPLEX(lp@oobj@env, lp@oobj@lp)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("writeProb", signature(lp = "optObj"),

    function(lp, fmt = "lp", ...) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                switch(fmt,
                   "lp"  = {
                       fl <- glpkAPI::writeLPGLPK(lp@oobj, ...)
                   },
                   "mps" = {
                       fl <- glpkAPI::writeMPSGLPK(lp@oobj, ...)
                   },
                   {
                       message("wrong format!")
                   }
                )
                out <- ifelse(fl == 0, TRUE, fl)
            },
            # ----------------------- #
            "clp" = {
                out <- clpAPI::saveModelCLP(lp@oobj, ...)
            },
            # ----------------------- #
            "lpSolveAPI" = {
                out <- lpSolveAPI::write.lp(lp@oobj, type = fmt, ...)
            },
            # ----------------------- #
            "cplex" = {
                tp  <- ifelse(is.null(fmt), NULL, toupper(fmt))
                fl  <- cplexAPI::writeProbCPLEX(lp@oobj@env, lp@oobj@lp,
                                                ftype = tp, ...)
                out <- ifelse(fl == 0, TRUE, fl)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("sensitivityAnalysis", signature(lp = "optObj"),

    function(lp, ...) {

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpk" = {
                out <- glpkAPI::printRangesGLPK(lp@oobj, ...)
                if (out == 0) {
                    message("wrote the file 'sar.txt'")
                }
                else {
                    warning("sensitivity analysis failed")
                }
            },
            # ----------------------- #
            "clp" = {
                wrong_solver_msg(lp, "sensitivityAnalysis", printOut = TRUE)
            },
            # ----------------------- #
            "lpSolveAPI" = {
                wrong_solver_msg(lp, "sensitivityAnalysis", printOut = TRUE)
            },
            # ----------------------- #
            "cplex" = {
                # number of columns and rows
                nc <- cplexAPI::getNumColsCPLEX(lp@oobj@env, lp@oobj@lp)
                nr <- cplexAPI::getNumRowsCPLEX(lp@oobj@env, lp@oobj@lp)
        
                out <- vector(mode = "list", length = 3)
                names(out) <- c("bound", "obj", "rhs")
                
                out[["bound"]] <- cplexAPI::boundSaCPLEX(lp@oobj@env,
                                                         lp@oobj@lp, 0, nc-1)
                out[["obj"]]   <- cplexAPI::objSaCPLEX(lp@oobj@env,
                                                       lp@oobj@lp, 0, nc-1)
                out[["rhs"]]   <- cplexAPI::rhsSaCPLEX(lp@oobj@env,
                                                       lp@oobj@lp, 0, nr-1)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(out)
    }
)

