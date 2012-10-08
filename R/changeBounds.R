#  changeBounds.R
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
# Function: changeBounds
#
#
#
#

changeBounds <- function(model, react, lb, ub) {

  if (!is(model, "modelorg")) {
      stop("needs an object of class modelorg!")
  }
  
  checkedIds <- checkReactId(model, react)
  if (!is(checkedIds, "reactId")) {
      stop("argument react is wrong")
  }

  if (missing(lb)) {
      lbnd <- rep(0, length(react))
  }
  else if (length(lb) == 1) {
      lbnd <- rep(lb, length(react))
  }
  else {
      lbnd <- lb
  }

  if (missing(ub)) {
      ubnd <- rep(0, length(react))
  }
  else if (length(ub) == 1) {
      ubnd <- rep(ub, length(react))
  }
  else {
      ubnd <- lb
  }

  if ((length(checkedIds) != length(lbnd)) || (length(checkedIds) != length(ubnd))) {
      stop("react, lb and ub must have the same length!")
  }
 
  lowbnd(model)[react_pos(checkedIds)] <- lbnd
  uppbnd(model)[react_pos(checkedIds)] <- ubnd

  return(model)

}

