#  checkOptSol.R
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
# Function: checkOptSol
#
# Checks if optimization returns "success"
#
# Parameters:
#     optsol:   Solution of an optimization.


checkOptSol <- function(optsol, onlywarn = FALSE) {

    if (!is(optsol, "optsol")) {
        stop("needs an object of class optsol!")
    }

    lp_check <- FALSE

    if (isTRUE(onlywarn)) {
        if (sum(lp_ok(optsol) != 0) != 0) {
            msg <- paste("some optimizations did not end successful",
                         "check results with checkOptSol()", sep = "\n")
            warning(msg, call. = FALSE)
        }
        else {
            lp_check <- TRUE
        }
    }
    else {
		lp_check <- checksol()

		if (solver(optsol) == "cplexAPI") {
			TEMP_ENV <- cplexAPI::openEnvCPLEX()
		}
		else {
			TEMP_ENV <- NULL
		}

		num_of_prob(lp_check) <- length(optsol)
		
		ec <- table(lp_ok(optsol))
		sc <- table(lp_stat(optsol))

		exit_code(lp_check)    <- as.integer(rownames(ec))
		exit_num(lp_check)     <- as.integer(ec)
		exit_meaning(lp_check) <- mapply(getMeanReturn, exit_code(lp_check),
									   MoreArgs = list(solver = solver(optsol)))

		status_code(lp_check)    <- as.integer(rownames(sc))
		status_num(lp_check)     <- as.integer(sc)
		status_meaning(lp_check) <- mapply(getMeanStatus, status_code(lp_check),
										MoreArgs = list(solver = solver(optsol),
										env = TEMP_ENV))

		if (!is.null(TEMP_ENV)) {
			cplexAPI::closeEnvCPLEX(TEMP_ENV)
		}
    }

    return(lp_check)
}

