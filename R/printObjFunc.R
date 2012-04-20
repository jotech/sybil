#  printObjFunc.R
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
# Function: printObjFunc
#
#
#
#


printObjFunc <- function(model) {

    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }

    cInd <- which(obj_coef(model) != 0)
  
    # check if there is an objective function
    if (length(cInd) == 0) {
        of <- "no objective function"
    }
    else {
        of <- paste(paste(obj_coef(model)[cInd],
                          react_id(model)[cInd],
                          sep = " * "), collapse = " + ")
    }

    return(of)  
  
}

