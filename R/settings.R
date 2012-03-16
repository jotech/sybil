#  settings.R
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



SYBIL_SETTINGS <- function(param, value) {

    if ( (missing(param)) && (missing(value)) ) {
       return(.SYBILenv$settings)
    }

    if (missing(value)) {
        return(.SYBILenv$settings[[param]])
    }
    
    if ( (length(param) != 1) ||
         ( (length(value) != 1) && (! (param == "SOLVER_CTRL_PARAM") ) ) ) {
        stop("arguments 'param' and 'value' must have a length of 1")
    }
    
    switch(param,
    
        "SOLVER" = {
            chmet <- checkDefaultMethod(value, NA)
            .SYBILenv$settings[["SOLVER"]]            <- chmet$sol
            .SYBILenv$settings[["METHOD"]]            <- chmet$met
            .SYBILenv$settings[["SOLVER_CTRL_PARAM"]] <- chmet$param
        },
    
        "METHOD" = {
            chmet <- checkDefaultMethod(SYBIL_SETTINGS("SOLVER"), value)
            .SYBILenv$settings[["SOLVER"]]            <- chmet$sol
            .SYBILenv$settings[["METHOD"]]            <- chmet$met
            .SYBILenv$settings[["SOLVER_CTRL_PARAM"]] <- chmet$param
        },
    
        "TOLERANCE" = {
            .SYBILenv$settings[["TOLERANCE"]] <- as.numeric(value)
        },
    
        "MAXIMUM" = {
            .SYBILenv$settings[["MAXIMUM"]] <- as.numeric(value)
        },
    
        "ALGORITHM" = {
            if ( (value == "FBA")            ||
                 (value == "linearMOMA")     ||
                 (value == "linearMOMA_COBRA") ) {
                .SYBILenv$settings[["ALGORITHM"]] <- as.character(value)
            }
            else {
                stop("ALGORITHM can be either 'FBA', ",
                     "'linearMOMA' or 'linearMOMA_COBRA'")
            }
        },
    
        "OPT_DIRECTION" = {
            if ( (value == "max") || (value == "min") ) {
                .SYBILenv$settings[["OPT_DIRECTION"]] <- as.character(value)
            }
            else {
                stop("OPT_DIRECTION can be either 'max' or 'min'")
            }
        },
    
        "PATH_TO_MODEL" = {
            if (file.exists(value)) {
                .SYBILenv$settings[["PATH_TO_MODEL"]] <- as.character(value)
            }
            else {
                stop("directory ", sQuote(value), " does not exist")
            }
        },
    
        "SOLVER_CTRL_PARAM" = {
            if (is.data.frame(value)) {
                .SYBILenv$settings[["SOLVER_CTRL_PARAM"]] <- value
            }
            else {
                stop("SOLVER_CTRL_PARAM must be data.frame")
            }
        },
    
        {
            stop("unknown parameter: ", sQuote(param))
        }
    )
}

