#  doInRound.R
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
# Function: .doInRound
#
#
#


.doInRound <- function(indRounds, maxRounds) {

    if (maxRounds < 1) {
        stop("Argument 'maxRounds' must be > 0")
    }

    if (all(!is.na(indRounds))) {
        if (is(indRounds, "character")) {

            # check if keyword 'WT' is there
            nDIR <- length(indRounds)
            if (indRounds[nDIR] == "WT") {
                DIR <- as.integer(indRounds[-nDIR])
                inclWT <- TRUE
            }
            else {
                DIR <- as.integer(indRounds)
                inclWT <- FALSE
            }

            # first value is treated as offset, if indRounds
            # contains at least two elements
            if (length(DIR) > 1) {
                offs <- as.integer(DIR[1])
                DIR <- DIR[-1]
            }
            else {
                offs <- as.integer(0)
            }

            # when we will run pre/post processing
            runPP <- lapply(abs(DIR), function(x) seq((x+offs),
                                                      maxRounds, by = x
                                                  )
                     )
            if (isTRUE(inclWT)) {
                runPP <- c(runPP, as.integer(1))
            }
            runPP <- sort(unique(unlist(runPP)))
            runPP <- runPP[runPP > 0]
        }
        else if (is(indRounds, "numeric")) {
            runPP <- sort(as.integer(indRounds[indRounds > 0]))
        }
        else {
            warning("Argument 'indRounds' must be numeric or character")
            runPP <- c(1:maxRounds)
        }
    }
    else {
        runPP <- c(1:maxRounds)
    }

    return(runPP)

}

