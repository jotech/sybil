#  fluxVar.R
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


################################################
# Function: fluxVar
#
# Performs a flux vaiability Analysis
# 
# The function fluxVar() is inspired by the function
# fluxVariability() contained in the COBRA Toolbox.
# The algorithm is the same.


fluxVar <- function(model, react = c(1:react_num(model)),
                    exex = FALSE, fld = FALSE, verboseMode = 2, ...) {

    if (!is(model, "modelorg")) {
      stop("needs an object of class modelorg!")
    }

#    if (missing(react)) {
#        react <- reactId((1:react_num(model)), react_id(model))
#    }
#    else {
#        if (!is(react, "reactId")) {
#            react <- checkReactId(model, react, needId = TRUE)
#        }
#    }

    # remove exchange reactions from analysis
    if (isTRUE(exex)) {
        exchReact <- findExchReact(model)
        ex <- react_pos(exchReact)
        intReact <- 1:react_num(model)
        intReact <- intReact[-ex]
        
        if (length(intReact) < 1) {
            stop("model contains no internal reactions!")
        }
    }
    else {
        intReact <- react
    }

    creact <- checkReactId(model, intReact, needId = TRUE)
    
    if (!is(creact, "reactId")) {
        stop("check argument react")
    }

#------------------------------------------------------------------------------#
#                            data strucktures                                  #
#------------------------------------------------------------------------------#

    # problem object
    lpmod <- sysBiolAlg(model, algorithm = "fv", ...)

    # simulation solutions
    nObj  <- 2 * length(react_pos(creact))
    obj   <- numeric(nObj)
    ok    <- integer(nObj)
    stat  <- integer(nObj)
    if (isTRUE(fld)) {
        flux <- Matrix::Matrix(0, nrow = react_num(model), ncol = nObj)
    }
    else {
        flux <- NA
    }

    # reactions to optimize
    num_of_probs_half <- length(creact)
    fvR               <- react_pos(creact)
    

#------------------------------------------------------------------------------#
#                               optimizations                                  #
#------------------------------------------------------------------------------#

    if (verboseMode > 0) { message("calculating min/max values ...") }

    if (verboseMode == 2) { progr <- sybil:::.progressBar() }

    # it is a bit faster, if we do first all minimisations and
    # afterwards all maximisations

    od <- "min"
    j <- 0

    if (verboseMode > 2) {
        cat(sprintf("%-5s %-5s %-15s %12s\n",
                    "dir", "no.", "reaction id", "flux rate"))
    }

    for (i in 1:nObj) {

        if (verboseMode == 2) {
            progr <- sybil:::.progressBar(i, nObj, progr)
        }

        j <- j + 1
        if (i == (num_of_probs_half+1)) {
            od <- "max"
            j <- 1
        }

        fluxVarSol <- optimizeProb(lpmod, lpdir = od,
                                   react = fvR[j], obj_coef = 1)

        if (verboseMode > 2) {
            cat(sprintf("%-5s %-5s %-15s %12s\n", od, fvR[j],
                        substr(react_id(creact)[j], 1, 15),
                        sprintf("%.6f", fluxVarSol$obj)))
        }

        obj[i]    <- fluxVarSol$obj
        ok[i]     <- fluxVarSol$ok
        stat[i]   <- fluxVarSol$stat
        if (isTRUE(fld)) {
            flux[,i]   <- fluxVarSol$flux
        }
    }
  

#------------------------------------------------------------------------------#
#                             save the results                                 #
#------------------------------------------------------------------------------#

    optsol <- new("optsol_fluxVar",
        mod_id       = mod_id(model),
        mod_key      = mod_key(model),
        solver       = solver(problem(lpmod)),
        method       = method(problem(lpmod)),
        algorithm    = algorithm(lpmod),
        num_of_prob  = as.integer(nObj),
        lp_num_cols  = nc(lpmod),
        lp_num_rows  = nr(lpmod),
        lp_obj       = as.numeric(obj),
        lp_ok        = as.integer(ok),
        lp_stat      = as.integer(stat),
        lp_dir       = factor(c(rep("min", floor(nObj/2)),
                                rep("max", floor(nObj/2)))),
        obj_coef     = obj_coef(model),
        obj_func     = printObjFunc(model),
        fldind       = fldind(lpmod),
        fluxdist     = fluxDistribution(flux),

        react        = creact
    )

    return(optsol)

}


 
