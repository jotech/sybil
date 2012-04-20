#  fluxVar.R
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


################################################
# Function: fluxVar
#
# Performs a flux vaiability Analysis
# 
# The function fluxVar() is inspired by the function
# fluxVariability() contained in the COBRA Toolbox.
# The algorithm is the same.


fluxVar <- function(model, react, percentage = 100,
                    tol = SYBIL_SETTINGS("TOLERANCE"),
                    lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
                    solver = SYBIL_SETTINGS("SOLVER"),
                    method = SYBIL_SETTINGS("METHOD"),
                    solverParm = SYBIL_SETTINGS("SOLVER_CTRL_PARM"),
                    fld = FALSE, verboseMode = 2, ...) {

    # check prerequisites
    if (!is(model, "modelorg")) {
      stop("needs an object of class modelorg!")
    }

    if (missing(react)) {
        react <- reactId((1:react_num(model)), react_id(model))
    }
    else {
        if (!is(react, "reactId")) {
            react <- checkReactId(model, react, needId = TRUE)
        }
    }

    # -------------------------------------------------------------- #
  
    saveOptSol <- function(saveMod, solution, j, fld) {

        lp_ok(saveMod)[j]   <- solution$ok
        lp_obj(saveMod)[j]  <- solution$obj
        #lp_obj(saveMod)[j] <- sybil:::.ceilValues(solution$lp_obj, tol = tol)
        lp_stat(saveMod)[j] <- solution$stat
      
        if (fld == TRUE) {
            fluxes(saveMod)[,j]  <- solution$flux
        }

        return(saveMod)
    }

    # -------------------------------------------------------------- #

    if (any(obj_coef(model) != 0)) {
        hasObj <- TRUE
        optimal <- simpleFBA(model, lpdir = lpdir,
                             solver = solver, method = method,
                             solverParm = solverParm, ...)
        if (optimal$ok == 0) {
            if (lpdir == "max") {
                obj <- sybil:::.floorValues(optimal$obj,
                                            tol = tol)*percentage/100
            }
            else {
                obj <- sybil:::.ceilValues(optimal$obj,
                                           tol = tol)*percentage/100
            }
        }
        else {
            stop("No optimal solution!")
        }
    }
    else {
        hasObj <- FALSE
    }
  
    lpmod <- prepProbObj(model,
                         nCols      = react_num(model),
                         nRows      = met_num(model),
                         solver     = solver,
                         method     = method,
                         lpdir      = lpdir,
                         solverParm = solverParm
             )

    # add a row to the problem
    if (hasObj == TRUE) {
        type <- ifelse(lpdir == "max", "L", "U")
        oind <- which(obj_coef(model) != 0)
        oval <- obj_coef(model)[oind]
        addRowsToProb(lpmod, met_num(model)+1,
                      type, obj, obj, list(oind), list(oval))
        changeObjCoefs(lpmod, 1:react_num(model), rep(0, react_num(model)))
    }


#------------------------------------------------------------------------------#
#                            data strucktures                                  #
#------------------------------------------------------------------------------#


    # new object for the solution
    fluxVariability <- optsol_fluxVar(solver = solver,
                                      method = method,
                                      nprob  = 2 * length(react_pos(react)),
                                      lpdir  = lpdir,
                                      ncols  = react_num(model),
                                      nrows  = met_num(model),
                                      objf   = printObjFunc(model),
                                      fld    = fld,
                                      rc     = react
                       )



    obj_values <- numeric(2*length(react_pos(react)))

    num_of_probs_half <- length(react_pos(react))

    if (verboseMode > 0) { message("calculating min/max values ...") }

    if (verboseMode == 2) { progr <- sybil:::.progressBar() }

    j <- 1
    for (i in 1:num_of_probs_half) {

        if (verboseMode == 2) {
            progr <- sybil:::.progressBar(i, num_of_probs_half, progr)
        }

        changeObjCoefs(lpmod, react_pos(react)[i], 1)

        setObjDir(lpmod, "min")
        fluxVarSol_MIN <- simpleFBA(lpmod, lpdir = "min",
                                    solver = solver,
                                    method = method, fld = fld,
                                    prCil = j, poCil = j, ...)

        setObjDir(lpmod, "max")
        fluxVarSol_MAX <- simpleFBA(lpmod, lpdir = "max",
                                    solver = solver,
                                    method = method, fld = fld,
                                    prCil = j+1, poCil = j+1, ...)

        changeObjCoefs(lpmod, react_pos(react)[i], 0)

        if (verboseMode > 2) {
            print(sprintf("%-5s %-15s %12s %12s", i,
                          substr(react_id(model)[react_pos(react)[i]], 1, 15),
                          sprintf("%.6f", fluxVarSol_MIN$obj),
                          sprintf("%.6f", fluxVarSol_MAX$obj)))
        }

        if (abs(fluxVarSol_MIN$obj) > abs(fluxVarSol_MAX$obj)) {
            fluxVariability <- saveOptSol(fluxVariability,
                                          fluxVarSol_MAX, j, fld)
            j <- j + 1
            fluxVariability <- saveOptSol(fluxVariability,
                                          fluxVarSol_MIN, j, fld)
        } else {
            fluxVariability <- saveOptSol(fluxVariability,
                                          fluxVarSol_MIN, j, fld)
            j <- j + 1
            fluxVariability <- saveOptSol(fluxVariability,
                                          fluxVarSol_MAX, j, fld)
        }

        j <- j + 1
  }
  
  return(fluxVariability)

}


 
