#  simpleFBA.R
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


# ---------------------------------------------------------------------------- #
# Function: simpleFBA
#
#
# The function simpleFBA() is in parts inspired by the functions
# optimizeCbModel() and solveCobraLP() contained in the COBRA Toolbox.


simpleFBA <- function(model,
                      react = NA,
                      lb = NA,
                      ub = NA,
                      lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
                      minTotalFlux = FALSE,
                      minDist = FALSE,
                      wtFluxes = NA,
                      solver = SYBIL_SETTINGS("SOLVER"),
                      method = SYBIL_SETTINGS("METHOD"),
                      fld = FALSE,
                      retOptSol = FALSE,
                      checkIds = TRUE,
                      prCmd = NA, poCmd = NA, prCil = NA, poCil = NA,
                      ...
                     ) {

    if ((!is(model, "modelorg")) && (!is(model, "optObj"))) {
        stop("needs an object of class modelorg or optObj!")
    }

    # check the argument react
    if (!any(is.na(react))) {
        if (is(model, "modelorg")) {
            # if model is of class "modelorg", react is given by the user
            check <- checkReactId(model, react = react)
            if (is(check, "reactId")) {
                react <- react_pos(check)
                del <- TRUE
            }
            else {
                stop("check argument react")
            }
        }
        else {
            # if model is of class "optObj", react is given by a
            # preceeding function
            if (!is(react, "numeric")) {
                stop("check argument react")
                retOptSol <- FALSE
            }
            else {
                del <- TRUE
            }
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
    else {
        del <- FALSE
    }


    if ((isTRUE(minTotalFlux)) && (isTRUE(minDist))) {
        stop("Use 'minTotalFlux' xor 'minDist'!")
    }

    if (isTRUE(minDist)) {
        # we need a wild type flux distribution
        fld <- TRUE
    }

    if ((!is.na(wtFluxes)) && (!is(wtFluxes, "numeric"))) {
        stop("argument wtFluxes must be numeric or NA")
    }

    did_fba <- FALSE


    # -------------------------------------------------------------- #
    # standard flux balance analysis
    # -------------------------------------------------------------- #

    if(any(is.na(wtFluxes))) {

        did_fba <- TRUE

        if (!is(model, "optObj")) {

            # prepare lp model for a standard FBA problem

            lpmod <- prepProbObj(model,
                                 nCols = react_num(model),
                                 nRows = met_num(model),
                                 solver = solver,
                                 method = method,
                                 lpdir = lpdir,
                                 ...
                                )
        }
        else {
            lpmod <- model
        }

        # if react is not empty and we need a wild type flux distribution,
        # we calculate it via FBA, but without deleting any flux.
        del_tmp <- ifelse(isTRUE(minDist), FALSE, del)

        if (isTRUE(del_tmp)) {
            # store default lower and upper bounds
            lowb_tmp <- getColsLowBnds(lpmod, react)
            uppb_tmp <- getColsUppBnds(lpmod, react)

            # change bounds of fluxes in react
            check <- changeColsBnds(lpmod, react, lb, ub)
            #model <- changeBounds(model, react, lb, ub, checkIds = checkIds)
        }

        # do some kind of preprocessing
        preP <- sybil:::.ppProcessing(lpprob = lpmod,
                                      ppCmd = prCmd,
                                      loopvar = prCil)

        # optimization
        lp_ok     <- solveLp(lpmod)
        #if (solver == "cplex") {
        #    a <- getdblQualCPLEX(lpmod@oobj$env, lpmod@oobj$lp, CPX_KAPPA)
        #    print(a)
        #}
        lp_obj    <- getObjVal(lpmod)
        lp_stat   <- getSolStat(lpmod)
        if (is.na(lp_stat)) {
            lp_stat <- lp_ok
        }

        lp_fluxes <- NA
        if (isTRUE(fld)) {
            lp_fluxes <- getFluxDist(lpmod)
        }

        # do some kind of postprocessing
        postP <- sybil:::.ppProcessing(lpprob = lpmod,
                                       ppCmd = poCmd,
                                       loopvar = poCil)

        # reset the default bounds
        if (isTRUE(del_tmp)) {
            check <- changeColsBnds(lpmod, react, lowb_tmp, uppb_tmp)
        }

        # clean up
        if (!is(model, "optObj")) {
            delProb(lpmod)
            remove(lpmod)
        }

    }


    # -------------------------------------------------------------- #
    # search for minimum flux distribution
    # -------------------------------------------------------------- #

    if(isTRUE(minTotalFlux)) {

        fld <- TRUE

        if (!is(model, "optObj")) {

            if (isTRUE(did_fba)) {
                check <- checkSolStat(lp_stat, solver = solver)
                if (length(check) > 0) {
                    warning("FBA solution ended not successfull!")
                    return(list(ok = lp_ok,
                                obj = lp_obj,
                                stat = lp_stat,
                                fluxes = lp_fluxes
                          )    )
                }
                else {
                    wtFluxes <- lp_obj
                }
            }

            # prepare model for a minimization of the total flux
            #  => we need a value for the objective function to be fixed.
            #     calculated by FBA (prior) or set by the user.

            nCols  <- 3*react_num(model)
            nRows  <- met_num(model) + 2*react_num(model) + 1
            lpmod <- prepProbObj(model, nCols, nRows,
                                 wtflux = wtFluxes,
                                 minTotalFluxFLAG = TRUE,
                                 solver = solver,
                                 method = method,
                                 lpdir  = "min",
                                 ...
                                )

            #if (isTRUE(del)) {
            #    react <- react + react_num(model)
            #}

        }
        else {
            lpmod <- model
        }

        if (isTRUE(del)) {
            # store default lower and upper bounds
            lowb_tmp <- getColsLowBnds(lpmod, react)
            uppb_tmp <- getColsUppBnds(lpmod, react)

            # change bounds of fluxes in react
            check <- changeColsBnds(lpmod, react, lb, ub)
            #model <- changeBounds(model, react, lb, ub, checkIds = checkIds)
        }

        # do some kind of preprocessing
        preP <- sybil:::.ppProcessing(lpprob = lpmod, ppCmd = prCmd)

        # optimization
        lp_ok     <- solveLp(lpmod)
        lp_obj    <- getObjVal(lpmod)
        lp_stat   <- getSolStat(lpmod)
        if (is.na(lp_stat)) {
            lp_stat <- lp_ok
        }

        lp_fluxesTMP <- getFluxDist(lpmod)
        lp_fluxes    <- lp_fluxesTMP[1:react_num(model)]
        #lp_fluxes    <- lp_fluxesTMP[(react_num(model)+1):(2*react_num(model))]

        # do some kind of postprocessing
        postP <- sybil:::.ppProcessing(lpprob = lpmod, ppCmd = poCmd)

        # reset the default bounds
        if (isTRUE(del)) {
            check <- changeColsBnds(lpmod, react, lowb_tmp, uppb_tmp)
        }

        # clean up
        if (!is(model, "optObj")) {
            delProb(lpmod)
            remove(lpmod)
        }

    }


    # -------------------------------------------------------------- #
    # search for minimum distance to a given flux distribution
    # -------------------------------------------------------------- #

    if(isTRUE(minDist)) {

        if (!is(model, "optObj")) {

            if (isTRUE(did_fba)) {
                check <- checkSolStat(lp_stat, solver = solver)
                if (length(check) > 0) {
                    warning("FBA solution ended not successfull!")
                    return(list(ok = lp_ok,
                                obj = lp_obj,
                                stat = lp_stat,
                                fluxes = lp_fluxes
                          )    )
                }
                else {
                    wtFluxes <- lp_fluxes
                }
            }

            # prepare model for a minimization of the total flux
            #  => we need a value for the objective function to be fixed.
            #     calculated by FBA (prior) or set by the user.

            nCols  <- 4*react_num(model)
            nRows  <- met_num(model) + 2*react_num(model)

            # problem object
            # we handle it like a 'linearMOMA' problem
            lpmod <- prepProbObj(model, nCols, nRows,
                                 wtflux = wtFluxes,
                                 MOMAflag = TRUE,
                                 solver = solver,
                                 method = method,
                                 lpdir  = lpdir,
                                 ...
                                )

            if (isTRUE(del)) {
                react <- react + react_num(model)
            }

        }
        else {
            lpmod <- model
        }

        if (isTRUE(del)) {
            # store default lower and upper bounds
            lowb_tmp <- getColsLowBnds(lpmod, react)
            uppb_tmp <- getColsUppBnds(lpmod, react)

            # change bounds of fluxes in react
            check <- changeColsBnds(lpmod, react, lb, ub)
            #model <- changeBounds(model, react, lb, ub, checkIds = checkIds)
        }

        # do some kind of preprocessing
        preP <- sybil:::.ppProcessing(lpprob = lpmod, ppCmd = prCmd)

        # optimization
        lp_ok     <- solveLp(lpmod)
        lp_obj    <- getObjVal(lpmod)
        lp_stat   <- getSolStat(lpmod)
        if (is.na(lp_stat)) {
            lp_stat <- lp_ok
        }

        lp_fluxesTMP <- getFluxDist(lpmod)
        lp_fluxes    <- lp_fluxesTMP[(react_num(model)+1):(2*react_num(model))]

        # do some kind of postprocessing
        postP <- sybil:::.ppProcessing(lpprob = lpmod, ppCmd = poCmd)

        # reset the default bounds
        if (isTRUE(del)) {
            check <- changeColsBnds(lpmod, react, lowb_tmp, uppb_tmp)
        }

        # clean up
        if (!is(model, "optObj")) {
            delProb(lpmod)
            remove(lpmod)
        }

    }


    # -------------------------------------------------------------- #
    # store solution
    # -------------------------------------------------------------- #

    optsol <- NULL

    if (!exists("lp_fluxes")) {
        stop("set either minDist or minTotalFlux to 'TRUE'")
    }

    if (isTRUE(retOptSol)) {

        # solution object
        optsol <- optsol_simpleFBA(solver = solver,
                                   nprob  = 1,
                                   lpdir  = lpdir,
                                   ncols  = react_num(model),
                                   nrows  = met_num(model),
                                   objf   = printObjFunc(model),
                                   fld    = fld
                                  )

        lp_ok(optsol)   <- lp_ok
        lp_obj(optsol)  <- lp_obj
        lp_stat(optsol) <- lp_stat
        if (!all(is.na(lp_fluxes))) {
            fluxes(optsol)[,1] <- lp_fluxes
        }

        method(optsol) <- method

        if (is(preP, "ppProc")) {
            preProc(optsol) <- preP
        }

        if (is(postP, "ppProc")) {
            postProc(optsol) <- postP
        }

        check <- validObject(optsol, test = TRUE)

        if (check != TRUE) {
            warning(paste("Validity check failed:", check, sep = "\n    "),
                    call. = FALSE
            )
        }

        checkOptSol(optsol, onlywarn = TRUE)

        #return(optsol)
    }
    else {
        optsol <- list(ok = lp_ok,
                       obj = lp_obj,
                       stat = lp_stat,
                       fluxes = lp_fluxes,
                       preP = preP,
                       postP = postP
                  )
    }

    return(optsol)
    #return(list(ok = lp_ok, obj = lp_obj, stat = lp_stat, fluxes = lp_fluxes))

}


