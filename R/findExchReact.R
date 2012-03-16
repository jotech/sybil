#  findExchReact.R
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
# Function: findExchReact
#
# The function findExchReact() is inspired by the function
# findExcRxns() contained in the COBRA Toolbox.
# The algorithm is the same.


findExchReact <- function(model) {

  if ( (!is(model, "modelorg")) &&
       !( (is(model, "Matrix")) || (is(model, "matrix")) ) ) {
      stop("needs an object of class modelorg, Matrix or matrix!")
  }

  if (is(model, "modelorg")) {
      St <- S(model)
  }
  else {
      St <- model
  }

  # columns with only one entry
  oneEntry <- apply(St, 2, function(x) sum(x != 0) == 1)

  # exchange reactions -- with a -1 or 1
  exchangeReact <- apply(St[ , oneEntry, drop = FALSE], 2, function(x) (sum(x == 1) == 1) | (sum(x == -1) == 1))

  # vector with the reaction id's of the exchange reactions
  ex <- c(1 : dim(St)[2])[oneEntry[exchangeReact]]


  # uptake reactions
  up <- NA
  if (is(model, "modelorg_irrev")) {
      up <- apply(St[,oneEntry], 2, function(x) (sum(x == 1) == 1))
  }

  if ((is(model, "modelorg")) && (!is(model, "modelorg_irrev"))) {
      up <- lowbnd(model)[ex] < 0
  }

  if (is(model, "modelorg")) {
      react <- list(
                    exchange = reactId(ex, react_id(model)[ex]),
                    uptake   = ex[up]
                    )
  }
  else {
      react <- oneEntry[exchangeReact]
  }

  return(react)

}
