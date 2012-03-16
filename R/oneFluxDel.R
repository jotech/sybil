#  oneFluxDel.R
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
# Function: oneFluxDel
#
# This function performs a "gene deletion analysis".
# In each iteration one gene is switched of (vi = 0)
# and the objective function will be computed.
#
# Parameters:
#      model:   an object of class model
#      react:   fluxes to delete
#      lpdir:   optimization direction
#     solver:   lp problem solver (glpk [default], lpSolve)
#     method:   simplex[default], interior
#        fld:   boolean, TRUE: get the flux distribution.
#        ...:   further arguments, passed to optimizer()
#
# Returns an object of class optsol_fluxdel.


oneFluxDel <- function(model, react, lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
                       solver = SYBIL_SETTINGS("SOLVER"),
                       method = SYBIL_SETTINGS("METHOD"),
                       fld = FALSE, ...) {

    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }

    if (missing(react)) {
        react <- reactId((1:react_num(model)), react_id(model))
    }
    else {
        if (!is(react, "reactId")) {
            react <- checkReactId(model, react)
        }
    }

    react <- sort(react_pos(react))
    num_of_prob = length(react)
  
    # new object for the solution
    optsol <- optsol_fluxdel(solver = solver,
                             nprob  = num_of_prob,
                             lpdir  = lpdir,
                             ncols  = react_num(model),
                             nrows  = met_num(model),
                             objf   = printObjFunc(model),
                             fld    = fld
                            )

    react_id(optsol) <- react_id(model)
    allGenes(optsol) <- allGenes(model)
    method(optsol) <- method

    dels(optsol)[2:(num_of_prob + 1),] <- react

    optsol <- optimizer(model = model, optsol = optsol, ...)

    return(optsol)

}

