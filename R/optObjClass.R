#  optObjClass.R
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
#                       definition of the class optObj                         #
#------------------------------------------------------------------------------#

setClass(Class = "optObj",
         representation(
              oobj     = "pointerToProb",
              solver   = "character",
              method   = "character",
              probType = "character"
         ),
         contains = "VIRTUAL"
)

# derivatives
#setClass(Class = "optObj_boot", contains = "optObj")


#------------------------------------------------------------------------------#
#                              user constructor                                #
#------------------------------------------------------------------------------#

optObj <- function(solver = SYBIL_SETTINGS("SOLVER"),
                   method = SYBIL_SETTINGS("METHOD"),
                   pType = "lp", prefix = "optObj", sep = "_") {

    validSoMe <- checkDefaultMethod(solver, method, pType)

    obj <- new(paste(prefix, validSoMe$sol, sep = sep),
               sv = validSoMe$sol,
               mt = validSoMe$met,
               pt = as.character(pType))

    return(obj)
}


#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

# contructor for class optObj
setMethod(f = "initialize",
          signature = "optObj",
          definition = function(.Object, sv, mt, pt) {

              if ( (!missing(sv)) && (!missing(mt)) && (!missing(pt)) ) {
                  
                  .Object@solver   <- as.character(sv)
                  .Object@method   <- as.character(mt)
                  .Object@probType <- as.character(pt)
                  
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

setMethod("show", signature(object = "optObj"),
    function(object) {
        if (length(probType(object)) > 0) {
            switch (probType(object),
                "lp" = {
                    cat("linear programming problem object\n")
                },
                "mip" = {
                    cat("mixed integer linear programming problem object\n")
                },
                "qp" = {
                    cat("continous problem object with quadratic objective\n")
                },
                {
                    cat("problem object of type ", probType(object),"\n")
                }
            )
            cat("solver:", solver(object), "\n")
            cat("method:", method(object), "\n")
            size <- tryCatch(dim(object), error = function(e) NA)
            if (any(is.na(size))) {
                cat("problem is not initialized\n")
            }
            else if (all(size == 0)) {
                cat("problem is currently empty\n")
            }
            else {
                cat("problem has", size[2],
                    ngettext(size[2], "variable", "variables"))
                cat(" and", size[1],
                    ngettext(size[1], "constraint", "constraints"), "\n")
            }
        }
        else {
            cat("empty problem object\n")
        }
    }
)


#------------------------------------------------------------------------------#
#                                 deprecated                                   #
#------------------------------------------------------------------------------#

# oobj
#------------------------------------------------------------------------------#

# all further methods are interface methods to the lp solvers

setMethod("loadMatrix", signature(lp = "optObj"),

    function(lp, ...) {

        .Deprecated(new = "loadLPprob", old = "loadMatrix")
        
        switch(lp@solver,
            # ----------------------- #
            "glpkAPI" = {
                glpkAPI::loadMatrixGLPK(lp@oobj, ...)
            },
            # ----------------------- #
            "clpAPI" = {
                clpAPI::loadMatrixCLP(lp@oobj, ...)
            },
            # ----------------------- #
            "lpSolveAPI" = {
                out <- loadMatrixPerColumnLPSOLVE(lp@oobj, ...)
            },
            # ----------------------- #
            "cplexAPI" = {
                out <- cplexAPI::chgCoefListCPLEX(lp@oobj@env, lp@oobj@lp, ...)
            },
            # ----------------------- #
            {
                wrong_type_msg(lp)
            }
        )

        return(TRUE)
    }
)


#------------------------------------------------------------------------------#

setMethod("loadProblemData", signature(lp = "optObj"),

    function(lp, model) {

        .Deprecated(new = "loadLPprob", old = "loadProblemData")
        
        if (!is(model, "modelorg")) {
            stop("needs an object of class modelorg!")
        }

        out <- FALSE
        switch(lp@solver,
            # ----------------------- #
            "glpkAPI" = {
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
            "clpAPI" = {
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
            "cplexAPI" = {
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

        .Deprecated(new = "loadLPprob", old = "loadProblemDataLM")

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
            "glpkAPI" = {
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
            "clpAPI" = {

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
            "cplexAPI" = {
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

        .Deprecated(new = "loadLPprob", old = "loadProblemDataMTF")

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
            "glpkAPI" = {
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
            "clpAPI" = {

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
            "cplexAPI" = {
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


