#  checkDefaultMethod.R
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


checkDefaultMethod <- function(solver, method) {

    # -------------------------------------------------------------- #
    # validate solver

    val_solver_ind <- match(solver, .SYBILenv$solvers)
    if (is.na(val_solver_ind)) {
        val_solver <- .SYBILenv$solvers[1]
    }
    else {
        val_solver <- .SYBILenv$solvers[val_solver_ind]
    }


    # -------------------------------------------------------------- #
    # validate method

    val_method_ind <- match(method, .SYBILenv$solverMethods[[val_solver]])
    if (is.na(val_method_ind)) {
        val_method <- .SYBILenv$solverMethods[[val_solver]][1]
    }
    else {
        val_method <- .SYBILenv$solverMethods[[val_solver]][val_method_ind]
    }


    # -------------------------------------------------------------- #

    ctrl_parm <- .SYBILenv$solverCtrlParm[[val_solver]][[val_method]]

    if (is.null(ctrl_parm)) {
        ctrl_parm <- as.data.frame(NA)
    }

    return(list(
                sol  = val_solver,
                met  = val_method,
                parm = ctrl_parm
           )
    )

}

