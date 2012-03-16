#  oneGeneDel.R
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
# Function: oneGeneDel
#
# This function performs a "gene deletion analysis".
# In each iteration all reactions corresponding to a
# gene were switched of (vi = 0) and the objective
# function will be computed.
#
# The function oneGeneDel() is inspired by the function
# singleGeneDeletion() contained in the COBRA Toolbox.


oneGeneDel <- function(model, geneList, lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
                       solver = SYBIL_SETTINGS("SOLVER"),
                       method = SYBIL_SETTINGS("METHOD"),
                       fld = FALSE, ...) {

    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }

    if (missing(geneList)) {
        if (length(allGenes(model)) < 1) {
            stop("Argument 'geneList' must contain at least one gene!")
        }
        else {
            geneList <- c(1:length(allGenes(model)))
        }
    }
 
    # translate the gene List in indices of allGenes(model)
    if (is(geneList, "character")) {
       geneList <- match(geneList, allGenes(model))
       if (any(is.na(geneList))) {
           stop("check genelist!")
       }
    }

    #geneList <- sort(geneList)

    gene_num <- length(geneList)

    # new object for the solution
    optsol <- optsol_genedel(solver = solver,
                             nprob  = gene_num,
                             lpdir  = lpdir,
                             ncols  = react_num(model),
                             nrows  = met_num(model),
                             objf   = printObjFunc(model),
                             fld    = fld
                            )

    react_id(optsol) <- react_id(model)
    allGenes(optsol) <- allGenes(model)
    method(optsol)   <- method
    
    dels(optsol)[2:(gene_num + 1),] <- geneList
    fluxdels(optsol)[1] <- NA

    optsol <- optimizer(model = model, optsol = optsol, ...)

    return(optsol)

}



