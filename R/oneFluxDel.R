#  oneFluxDel.R
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


################################################
# Function: oneFluxDel
#
# This function performs a "gene deletion analysis".
# In each iteration one gene is switched of (vi = 0)
# and the objective function will be computed.
#
# The function oneGeneDel() is inspired by the function
# singleRxnDeletion() contained in the COBRA Toolbox.


oneFluxDel <- function(model, react, ...) {

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

    optsol <- optimizer(model = model,
                        delete = matrix(react, ncol = 1),
                        geneFlag = FALSE,
                        ...)

    return(optsol)

}

