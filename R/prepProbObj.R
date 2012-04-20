#  prepProbObj.R
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
# Function: prepProbObj
#
# Creates a problem object for different solvers.
# 
# The function needs an model object of the class model,
# built by read????mod.
#
# Parameters:
#     model:  the model object of class model
#    solver:  the problem solver (default: SYBIL_SETTINGS("SOLVER"))
#     lpdir:  direction of optimization (default: max)
#
# Return values:
#        lp:  a problem object


prepProbObj <- function(model, nCols, nRows,
                        wtflux = NA,
						MOMAflag = FALSE,
						COBRAflag = FALSE,
						minTotalFluxFLAG = FALSE,
						minDistFLAG = FALSE,
                        scaling = NA,
                        #alg = SYBIL_SETTINGS("ALGORITHM"),
                        solver = SYBIL_SETTINGS("SOLVER"),
                        method = SYBIL_SETTINGS("METHOD"),
                        lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
                        solverParm = SYBIL_SETTINGS("SOLVER_CTRL_PARM")
                       ) {

    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }

    if (sum(MOMAflag, minTotalFluxFLAG, minDistFLAG) > 1) {
	    stop("Choose 'MOMA', 'minTotalFlux' xor 'minDist'!")
	}

    check <- TRUE
    
    # create a new empty problem object
    lp <- optObj(solver = solver, method = method)
    lp <- initProb(lp, nrows = nRows, ncols = nCols)

    # -------------------------------------------------------------- #
    # set control parameters
    # -------------------------------------------------------------- #

    if (!any(is.na(solverParm))) {
        check <- setSolverParm(lp, solverParm)
    }
    
    
    
#     if (all(is.na(solverParm))) {
#         # check weather we have some default parameters
#         #parm <- checkDefaultMethod(solver, method)
#         #solverParm <- parm$parm
#         check <- setSolverParm(lp, SYBIL_SETTINGS("SOLVER_CTRL_PARM"))
#     } else if (all(solverParm == "none")) {
#         solverParm <- NA
#     }
# 
#     if (!all(is.na(solverParm))) {
#         check <- setSolverParm(lp, solverParm)
#     }


    # -------------------------------------------------------------- #
    # load problem data
    # -------------------------------------------------------------- #

    # load problem data
    if (isTRUE(MOMAflag)) {
		#COBRAflag <- ifelse(alg == "linearMOMA_COBRA", TRUE, FALSE)
		check <- loadProblemDataLM(lp, model, wtflux, nCols, nRows, COBRAflag, lpdir = lpdir)
		check <- setObjDir(lp, "min")
	}
	else if (isTRUE(minTotalFluxFLAG)) {
		check <- loadProblemDataMTF(lp, model, wtflux, nCols, nRows)
		check <- setObjDir(lp, "min")
	}
	else {
		check <- loadProblemData(lp, model)
		check <- setObjDir(lp, lpdir)
	}
    

    # -------------------------------------------------------------- #
    # scaling
    # -------------------------------------------------------------- #

    if (!is.na(scaling)) {
        check <- scaleProb(lp, scaling)
    }
    
    check <- TRUE
    
    # return problem object
    if (isTRUE(check)) {
        return(lp)
    }
    else {
        return(check)
    }
  
}

