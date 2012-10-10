#  geneDeletion.R
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
# Function: geneDeletion
#
# This function performs a "n gene deletion analysis".
# In each iteration m genes are switched of (vi = 0)
# and the objective function will be computed.
#
# The function geneDeletion() is inspired by the functions singleGeneDeletion()
# and doubleGeneDeletion() contained in the COBRA Toolbox.


geneDeletion <- function(model, genes, combinations = 1, ...) {

    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }

    num_genes <- length(allGenes(model))

    # delGenes containes pointers to gene id's in allGenes(model)
    if (missing(genes)) {
        delGenes <- num_genes
    }
    else {
        delGenes <- genes
    }

    # make deletion matrix
    if (!is(delGenes, "matrix")) {
        delGenes <- combn(x = delGenes, m = combinations)
    }

    if (typeof(delGenes) == "character") {
        # check the id's
        dimdg    <- dim(delGenes)
        delGenes <- match(delGenes, allGenes(model))
        if (any(is.na(delGenes))) {
            stop("some gene id's are unknown, check argument genes")
        }
        else {
            attr(delGenes, which = "dim") <- dimdg
        }
    }
    else {
        if (any(delGenes > num_genes)) {
        #if ( (max(delGenes) > num_genes) || (min(delGenes) < 0) ) {
            stop("values of genes must be in [0, length(allGenes(model))]")
        }
    }

    # number of optimizations
    num_opt <- ncol(delGenes)

#------------------------------------------------------------------------------#
#                               run optimization                               #
#------------------------------------------------------------------------------#

    optsol <- optimizer(model = model,
                        delete = t(delGenes),
                        geneFlag = TRUE,
                        ...)

    return(optsol)

}

