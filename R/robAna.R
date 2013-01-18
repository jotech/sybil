#  robAna.R
#  FBA and friends with R.
#
#  Copyright (C) 2010-2013 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
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


################################################
# Function: robAna
#
#
# The function robAna() is inspired by the function
# robustnessAnalysis() contained in the COBRA Toolbox.
# The algorithm is the same.


robAna <- function(model, ctrlreact,
                   numP = 20, verboseMode = 2, fld = FALSE, ...) {

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
    minmaxsol <- function(mod, ...) {

        obj_coef(mod) <- integer(react_num(mod))
        lp <- sysBiolAlg(mod, algorithm = "fba", ...)

        tmp_sol_min  <- optimizeProb(lp,
                                     react = ctrlr,
                                     obj_coef = 1,
                                     lpdir = "min")
        tmp_sol_max  <- optimizeProb(lp,
                                     react = ctrlr,
                                     obj_coef = 1,
                                     lpdir = "max")

        if ( (tmp_sol_min$ok != 0) || (tmp_sol_max$ok != 0) ) {
            stop("Optimization for min/max solution ended not successfull!")
        }

        delProb(problem(lp))

        return(c(tmp_sol_min$obj, tmp_sol_max$obj))
    }
    # ------------------------------------------------------------------------ #

#------------------------------------------------------------------------------#
#                       minimum and maximum solution                           #
#------------------------------------------------------------------------------#


    mm    <- minmaxsol(model, ...)

    lpmod <- sysBiolAlg(model, algorithm = "fba", ...)

    # sequence of numP numbers between lpmin and lpmax,
    # all with the same distance
    ctrlfl <- seq(mm[1], mm[2], length.out = numP)

#------------------------------------------------------------------------------#
#                                optimization                                  #
#------------------------------------------------------------------------------#

    obj   <- numeric(numP)
    ok    <- integer(numP)
    stat  <- integer(numP)
    if (isTRUE(fld)) {
        flux <- Matrix::Matrix(0, nrow = nc(lpmod), ncol = numP)
    }
    else {
        flux <- NA
    }

    if (verboseMode > 0)  { message("running optimizations ", appendLF = FALSE) }
    if (verboseMode == 1) { message("... ", appendLF = FALSE) }
    if (verboseMode == 2) { message("") }

    for (i in 1:numP){

        if (verboseMode == 2) { sybil:::.progressDots(5, i, numP) }

        sol <- optimizeProb(lpmod,
                            react = ctrlr,
                            lb = ctrlfl[i],
                            ub = ctrlfl[i])

        if (verboseMode > 2) {
            print(sprintf("%-5s %-15s %12s", i,
                          substr(react_id(model)[ctrlr], 1, 15),
                          sprintf("%.6f", sol$obj)))
        }

        ok[i]   <- sol$ok
        obj[i]  <- sol$obj
        stat[i] <- sol$stat

        if (fld == TRUE) {
            flux[ ,i] <- sol$flux
        }
    }

    #ctrlfl(robust) <- abs(ctrlfl(robust))
    if (verboseMode > 0) { message("OK") }

    optsol <- new("optsol_robAna",
        mod_id       = mod_id(model),
        mod_key      = mod_key(model),
        solver       = solver(problem(lpmod)),
        method       = method(problem(lpmod)),
        algorithm    = algorithm(lpmod),
        num_of_prob  = as.integer(numP),
        lp_num_cols  = nc(lpmod),
        lp_num_rows  = nr(lpmod),
        lp_obj       = as.numeric(obj),
        lp_ok        = as.integer(ok),
        lp_stat      = as.integer(stat),
        lp_dir       = factor(getObjDir(problem(lpmod))),
        obj_coef     = obj_coef(model),
        obj_func     = printObjFunc(model),
        fldind       = fldind(lpmod),
        fluxdist     = fluxDistribution(flux),

        ctrlr        = tmp,
        ctrlfl       = as.numeric(ctrlfl)
    )

    return(optsol)
}



