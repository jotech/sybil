#  optObj_basicfunc.R
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


#------------------------------------------------------------------------------#

# print a warning
wrong_type_msg <- function(lp) {
    warning(paste("Slot oobj of", lp@solver,
                  "is not a valid pointer to a valid solver."
            )
    )
}


wrong_solver_msg <- function(lp, method, printOut = TRUE) {
    if (isTRUE(printOut)) {
        warning(paste("Requested method", method,
                      "is not available for solver", lp@solver
                )
        )
    }
}


#------------------------------------------------------------------------------#

# which optimizations did not end successful
checkSolStat <- function(stat, solver = SYBIL_SETTINGS("SOLVER")) {
    out <- FALSE
    switch(solver,
        "glpk" = {
            out <- which(stat != 5)
        },
        "clp" = {
            out <- which(stat != 0)
        },
        "lpSolveAPI" = {
            out <- which(stat != 0)
        },
        "cplex" = {
            out <- which(stat != 1)
        },
        {
            warning("not a valid solver")
        }
    )
    return(out)
}


#------------------------------------------------------------------------------#

getMeanReturn <- function(code, solver = SYBIL_SETTINGS("SOLVER")) {
    out <- FALSE
    switch(solver,
        "glpk" = {
            out <- return_codeGLPK(code)
        },
        "clp" = {
            out <- return_codeCLP(code)
        },
        "lpSolveAPI" = {
            out <- return_codeLPSOLVE(code)
        },
        "cplex" = {
            out <- return_codeCPLEX(code)
        },
        {
            warning("not a valid solver")
        }
    )
    return(out)
}


#------------------------------------------------------------------------------#

getMeanStatus <- function(code,
                          solver = SYBIL_SETTINGS("SOLVER"), env = NULL) {
    out <- FALSE
    switch(solver,
        "glpk" = {
            out <- status_codeGLPK(code)
        },
        "clp" = {
            out <- status_codeCLP(code)
        },
        "lpSolveAPI" = {
            out <- "see return code"
        },
        "cplex" = {
            out <- status_codeCPLEX(env, code)
        },
        {
            warning("not a valid solver")
        }
    )
    return(out)
}


#TIMEOUT <- "timeout"
#INFINITE <- "infinite"
