#  checkReactId.R
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
# Function: checkReactId
#
# 
# 
#

checkReactId <- function(model, react, needId = FALSE) {

  if (is(react, "reactId")) {
      return(react)
  }

  if (!is(model, "modelorg")) {
    stop("needs an object of class modelorg!")
  }
  
#------------------------------------------------------------------------------#
#                           if "react" is numeric                              #
#------------------------------------------------------------------------------#

# If react is numeric (or integer), we only need to check, weather no element of
# react is larger than the number of reactions and if all elements are positive.

  if (is.numeric(react) || is.integer(react)) {
      num_react <- react_num(model)

      if ((max(react) > num_react) || (!all(react > 0))) {
          warning("Check argument react!")
          return(NA)
      }

      if (needId == FALSE) {
          checkedIds <- reactId(react, character(length(react)))
          return(checkedIds)
      }

      checkedIds <- reactId(react, react_id(model)[react])
      return(checkedIds)

  }

#------------------------------------------------------------------------------#
#                          if "react" is character                             #
#------------------------------------------------------------------------------#

  reaction_ids <- react_id(model)

  # find the vector indices in the model of react
  pos <- match(react, reaction_ids, nomatch = 0)
  if (all(pos != 0)) {
      checkedIds <- reactId(pos, react)
      return(checkedIds)
      #return(list(id = react, pos = pos))
  }

  # if we cannot find some entries in react, we try grep
  else {
      #print(pos)
      # only those, we could not find above
      pos_null <- which(pos == 0)

      # the grep
      react_null <- sapply(pos_null, function(x) grep(react[x], reaction_ids, ignore.case = TRUE, value = TRUE))
      
      #react  <- sapply(react, function(x) grep(x, Dmodel@react_id, ignore.case = TRUE, value = TRUE))

      # check weather all results are unique
      #print(is(react_null))
      #print(react_null)
      len    <- sapply(react_null, length)
      n_uniq <- which(len != 1)

      #print(len)

      if ((length(n_uniq) != 0) || (is(react_null, "matrix"))) {
          warning("Some reaction id's were not found or are ambiguous.", call. = FALSE)
          #print(n_uniq)
          #return(FALSE)
          return(NA)
      }

      pos_react <- match(react_null, reaction_ids)
      #print(pos_null)
      #print(pos_react)
      #print(react_null)

      pos[pos_null]   <- pos_react
      react[pos_null] <- react_null
      #print(pos)
      #print(react)
      #pos <- match(react, model@react_id)

      checkedIds <- reactId(pos, react)
      return(checkedIds)
      #return(list(id = react, pos = pos))

  }

                                      
}
