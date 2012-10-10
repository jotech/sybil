#  changeObjFunc.R
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
# Function: changeObjFunc
#
#
#
#

changeObjFunc <- function(lpmodel, react, obj_coef, checkIds = TRUE) {

  if (!is(lpmodel, "modelorg")) {
      stop("needs an object of class modelorg!")
  }

  if (missing(obj_coef)) {
      obj_coef <- c(rep(1, length(react)))
  }
      
  if (length(react) != length(obj_coef)) {
      stop("react and obj_coef must have the same length!")
  }

  if (checkIds == TRUE) {
      checkedIds <- checkReactId(lpmodel, react)
  }
  else {
      if (!is(react, "numeric") && (!is(react, "integer"))) {
          stop("react has to be integer or numeric")
      }
      checkedIds <- reactId(react, character(length(react)))
  }
  
  #print(is(check))
  if (!is(checkedIds, "reactId")) {
      stop("Check your reaction Id's")
  }
 
#  if (is.logical(check) && (check == FALSE)) {
#      stop("Check your reaction Id's")
#  }

  # set all objective coefficients to zero
  obj_coef(lpmodel) <- numeric(react_num(lpmodel))

  obj_coef(lpmodel)[react_pos(checkedIds)] <- obj_coef

  return(lpmodel)

}

