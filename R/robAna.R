#  robAna.R
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
# Function: robAna
#
#
# The function robAna() is inspired by the function
# robustnessAnalysis() contained in the COBRA Toolbox.
# The algorithm is the same.


robAna <- function(model,
                   ctrlreact,
                   numP = 20,
                   lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
                   solver = SYBIL_SETTINGS("SOLVER"),
                   method = SYBIL_SETTINGS("METHOD"),
                   solverParm = SYBIL_SETTINGS("SOLVER_CTRL_PARAM"),
                   fld = FALSE, verboseMode = 2, ...) {

    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }
    
    if (length(ctrlreact) != 1) {
        stop("Please enter exactly one control reaction.")
    }
    
    tmp <- checkReactId(model, ctrlreact)
    
    if (!is(tmp, "reactId")) {
        stop("Check control reaction!")
    }
    
    ctrlr <- react_pos(tmp)

    # ------------------------------------------------------------------------ #
    minmaxsol <- function(mod, mmdir) {

        tmp_sol  <- simpleFBA(mod, lpdir = mmdir)
        if (tmp_sol$ok != 0) {
            stop("Optimization for min/max solution ended not successfull!")
        }

        return(tmp_sol$obj)
    }
    # ------------------------------------------------------------------------ #


#------------------------------------------------------------------------------#
#                       minimum and maximum solution                           #
#------------------------------------------------------------------------------#

    tmpModel <- changeObjFunc(model, ctrlr)

    lpmin <- minmaxsol(mod = tmpModel, mmdir = "max")
    lpmax <- minmaxsol(mod = tmpModel, mmdir = "min")

    remove(tmpModel)


    # sequence of numP numbers between lpmin and lpmax,
    # all with the same distance
    ctrlfl <- seq(lpmin, lpmax, length.out = numP)


    # new object for the solution
    robust <- optsol_robAna(solver = solver,
                            method = method,
                            nprob  = numP,
                            lpdir  = lpdir,
                            ncols  = react_num(model),
                            nrows  = met_num(model),
                            objf   = printObjFunc(model),
                            fld    = fld,
                            cr     = tmp,
                            crf    = ctrlfl)


#------------------------------------------------------------------------------#
#                                optimization                                  #
#------------------------------------------------------------------------------#

    lpmod <- prepProbObj(model,
                         nCols      = react_num(model),
                         nRows      = met_num(model),
                         solver     = solver,
                         method     = method,
                         lpdir      = lpdir,
                         solverParm = solverParm)

    if (verboseMode > 0)  { message("running optimizations ", appendLF = FALSE) }
    if (verboseMode == 1) { message("... ", appendLF = FALSE) }
    if (verboseMode == 2) { message("") }

    for (i in 1:numP){

        if (verboseMode == 2) { .progressDots(5, i, numP) }

        optsol <- simpleFBA(lpmod,
                            react = ctrlr,
                            lb = ctrlfl[i],
                            ub = ctrlfl[i],
                            lpdir = lpdir,
                            solver = solver,
                            method = method,
                            fld = fld,
                            checkIds = FALSE,
                            prCil = i,
                            poCil = i, ...)

        if (verboseMode > 2) {
            print(sprintf("%-5s %-15s %12s", i,
                          substr(react_id(model)[ctrlr], 1, 15),
                          sprintf("%.6f", optsol$obj)))
        }

        lp_ok(robust)[i]   <- optsol$ok
        lp_obj(robust)[i]  <- optsol$obj
        lp_stat(robust)[i] <- optsol$stat

        if (fld == TRUE) {
            fluxes(robust)[,i] <- optsol$flux
        }
    }

    #ctrlfl(robust) <- abs(ctrlfl(robust))
    if (verboseMode > 0) { message("OK") }

    return(robust)
}



