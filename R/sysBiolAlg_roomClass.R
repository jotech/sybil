#  sysBiolAlg_roomClass.R
#  FBA and friends with R.
#
#  Copyright (C) 2010-2012 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of sybil.
#
#  sybil is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  sybil is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sybil.  If not, see <http://www.gnu.org/licenses/>.


# Original code written by Marc Andre Daxer during his bachelor thesis:
# "Analysis of Gene Defects with Mixed Integer Linear Programming" (2011
# at the Heinrich-Heine-University Duesseldorf, Dpt. for Bioinformatics).
# He wrote the package sybilROOM.


#------------------------------------------------------------------------------#
#                 definition of the class sysBiolAlg_room                      #
#------------------------------------------------------------------------------#

setClass(Class = "sysBiolAlg_room",
         representation(
             wu  = "numeric",
             wl  = "numeric",
             fnc = "integer",
             fnr = "integer"
         ),
         contains = "sysBiolAlg"
)


#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

# contructor for class sysBiolAlg_room
setMethod(f = "initialize",
          signature = "sysBiolAlg_room",
          definition = function(.Object,
                                model,
                                wtflux,
                                delta = 0.03,
                                epsilon = 0.001,
                                LPvariant = FALSE,
                                scaling = NULL, ...) {

              if ( ! missing(model) ) {

                  if (missing(wtflux)) {
                      tmp <- sybil:::.generateWT(model, ...)
                      wtflux <- tmp$fluxes[tmp$fldind]
                  }

                  stopifnot(is(model, "modelorg"),
                            is(wtflux, "numeric"),
                            is(delta, "numeric"),
                            is(epsilon, "numeric"),
                            is(LPvariant, "logical"))
                  
                  stopifnot(length(wtflux) == react_num(model))

                  #    the problem: minimize
                  #
                  #             |       |
                  #         S   |   0   |  = 0
                  #             |       |
                  #       -------------------------
                  #       1     |       |
                  #         1   | -vwl  |  >= wl
                  #           1 |       |
                  #       -------------------------
                  #       1     |       |
                  #         1   | -vwu  |  <= wu
                  #           1 |       |
                  #       -------------------------
                  #
                  #        lb   |   0
                  #        ub   |   1
                  # ctype: C    |   B
                  # obj:   0    |   1


                  # ---------------------------------------------
                  # problem dimensions
                  # ---------------------------------------------

                  nc <- react_num(model)
                  nr <- met_num(model)

                  nCols <- (2 * nc)
                  nRows <- (nr + 2 * nc)

                  absMAX <- SYBIL_SETTINGS("MAXIMUM")


                  # ---------------------------------------------
                  # constraint matrix
                  # ---------------------------------------------

                  # ROOM-Boundaries
                  
                  # if we use the formulation as linear program, delta and
                  # epsilon are set to zero (see Shlomi et al.)
                  #if (isTRUE(LPvariant)) {
                  #    wu <- wtflux
                  #    wl <- wtflux
                  #}
                  #else {
                      wu <- wtflux + delta * abs(wtflux) + epsilon
                      wl <- wtflux - delta * abs(wtflux) - epsilon
                  #}

                  vwu <- uppbnd(model) - wu
                  vwl <- lowbnd(model) - wl


                  # the initial matrix dimensions
                  LHS <- Matrix::Matrix(0, 
                                        nrow = nRows,
                                        ncol = nCols,
                                        sparse = TRUE)

                  # rows for the flux variables
                  LHS[1:nr,1:nc] <- S(model)

                  # location of the wild type strain
                  fi <- c(1:nc)

                  # constraint-matrix for ROOM
                  diag(LHS[(nr+1)    :(nr+nc),1     :nc   ]) <-        1
                  diag(LHS[(nr+nc+1) :nRows  ,1     :nc   ]) <-        1
                  diag(LHS[(nr+1)    :(nr+nc),(nc+1):nCols]) <- vwl * -1
                  diag(LHS[(nr+nc+1) :nRows  ,(nc+1):nCols]) <- vwu * -1


                  # ---------------------------------------------
                  # columns
                  # ---------------------------------------------

                  clower <- c(lowbnd(model), rep(0, nc))
                  cupper <- c(uppbnd(model), rep(1, nc))

                  if (isTRUE(LPvariant)) {
                      ctype  <- NULL
                      pt     <- "lp"
                  }
                  else {
                      ctype  <- c(rep("C", nc),  rep("B", nc))
                      pt     <- "mip"
                  }

                  # ---------------------------------------------
                  # rows
                  # ---------------------------------------------

                  rlower <- c(rhs(model), wl, wu)
                  rupper <- c(rhs(model), rep(absMAX, nc), wu)
                  rtype  <- c(rep("E", nr), rep("L", nc), rep("U", nc))


                  # ---------------------------------------------
                  # objective function
                  # ---------------------------------------------

                  cobj <- c(rep(0, nc), rep(1, nc))


                  # ---------------------------------------------
                  # build problem object
                  # ---------------------------------------------

                  .Object <- callNextMethod(.Object,
                                            alg        = "room",
                                            pType      = pt,
                                            scaling    = scaling,
                                            fi         = fi,
                                            nCols      = nCols,
                                            nRows      = nRows,
                                            mat        = LHS,
                                            ub         = cupper,
                                            lb         = clower,
                                            obj        = cobj,
                                            rlb        = rlower,
                                            rub        = rupper,
                                            rtype      = rtype,
                                            ctype      = ctype,
                                            lpdir      = "min",
                                            ...)

                  .Object@wu  <- as.numeric(wu)
                  .Object@wl  <- as.numeric(wl)
                  .Object@fnr <- as.integer(nr)
                  .Object@fnc <- as.integer(nc)

#
#                  # ---------------------------------------------
#                  # build problem object
#                  # ---------------------------------------------
#
#                  lp <- optObj(solver = solver, method = method, pType = "mip")
#                  lp <- initProb(lp, nrows = nRows, ncols = nCols)
#
#                  # ---------------------------------------------
#                  # set control parameters
#                  # ---------------------------------------------
#
#                  if (!any(is.na(solverParm))) {
#                      setSolverParm(lp, solverParm)
#                  }
#    
#
#                  loadLPprob(lp,
#                             nCols = nCols,
#                             nRows = nRows,
#                             mat   = LHS,
#                             ub    = cupper,
#                             lb    = clower,
#                             obj   = cobj,
#                             rlb   = rlower,
#                             rub   = rupper,
#                             rtype = rtype,
#                             ctype = ctype,
#                             lpdir = "min"
#                  )
#                  
#                  if (!is.null(scaling)) {
#                      scaleProb(lp, scaling)
#                  }
#
#                  .Object@problem   <- lp
#                  .Object@algorithm <- "room"
#                  .Object@nr        <- as.integer(nRows)
#                  .Object@nc        <- as.integer(nCols)
#                  .Object@fldind    <- as.integer(fi)
#                  .Object@wu        <- as.numeric(wu)
#                  .Object@wl        <- as.numeric(wl)
#                  .Object@fnr       <- as.integer(nr)
#                  .Object@fnc       <- as.integer(nc)
#                  validObject(.Object)
                  
              }
              return(.Object)
          }
)


#------------------------------------------------------------------------------#

setMethod("optimizeProb", signature(object = "sysBiolAlg_room"),
    function(object,
             react = NA,
             lb = NA,
             ub = NA,
             prCmd = NA, poCmd = NA,
             prCil = NA, poCil = NA) {


        # check the argument react
        if (any(is.na(react))) {
            del <- FALSE
        }
        else {
            if (is(react, "numeric")) {
                del <- TRUE
            }
            else {
                stop("argument 'react' must be numeric")
            }
    
            if (any(is.na(lb))) {
                lb <- rep(0, length(react))
            }
            else {
                if (!is(lb, "numeric")) {
                    stop("argument lb must be numeric")
                }
                if (length(lb) != length(react)) {
                    stop("argument react and lb must have same length")
                }
            }
    
            if (any(is.na(ub))) {
                ub <- rep(0, length(react))
            }
            else {
                if (!is(ub, "numeric")) {
                    stop("argument ub must be numeric")
                }
                if (length(ub) != length(react)) {
                    stop("argument react and ub must have same length")
                }
            }
        }

   
        # -------------------------------------------------------------- #
        # modifications to problem object
        # -------------------------------------------------------------- #
    
        fi    <- fldind(object)
        wu    <- object@wu
        wl    <- object@wl
        lpmod <- problem(object)

        if (isTRUE(del)) {
            # store default lower and upper bounds
            lowb_tmp <- getColsLowBnds(lpmod, fi[react])
            uppb_tmp <- getColsUppBnds(lpmod, fi[react])
    
            # change bounds of fluxes in react
            check <- changeColsBnds(lpmod, fi[react], lb, ub)

            # change constraint matrix and objective coefficients
            vwu <- (ub - wu[fi[react]]) * -1
            vwl <- (lb - wl[fi[react]]) * -1

            ri <- react + object@fnr
            ci <- react + object@fnc
            for (i in seq(along = react)) {
                changeMatrixRow(lpmod, ri[i], c(react[i], ci[i]), c(1, vwl[i]))
                changeMatrixRow(lpmod,
                                ri[i]+object@fnc,
                                c(react[i], ci[i]),
                                c(1, vwu[i]))
            }
            
            #changeObjCoefs(lpmod, ci, rep(0, length(react)))
            
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
    
        # reset the default values
        if (isTRUE(del)) {
            check <- changeColsBnds(lpmod, fi[react], lowb_tmp, uppb_tmp)

            vwu <- (uppb_tmp - wu[fi[react]]) * -1
            vwl <- (lowb_tmp - wl[fi[react]]) * -1

            for (i in seq(along = react)) {
                changeMatrixRow(lpmod, ri[i], c(react[i], ci[i]), c(1, vwl[i]))
                changeMatrixRow(lpmod,
                                ri[i]+object@fnc,
                                c(react[i], ci[i]),
                                c(1, vwu[i]))
            }
            
            #changeObjCoefs(lpmod, ci, rep(1, length(react)))

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

