#  reactIdClass.R
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


# reactIdClass


#------------------------------------------------------------------------------#
#                     definition of the class reactId                          #
#------------------------------------------------------------------------------#

setClass("reactId",
         representation(
              react_pos = "integer",
              react_id  = "character"
         ),
         validity = sybil:::.validreactId
)


#------------------------------------------------------------------------------#
#                              user constructor                                #
#------------------------------------------------------------------------------#

reactId <- function(pos, react) {
    if (missing(pos) || missing(react)) {
        stop("Creating an object of class reactId needs pos and react!")
    }
    pos   <- as.integer(pos)
    react <- as.character(react)
    obj   <- new("reactId", react_pos = pos, react_id = react)
    return(obj)
}


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# position
setMethod("react_pos", signature(object = "reactId"),
          function(object) {
              return(object@react_pos)
          }
)

setReplaceMethod("react_pos", signature = (object = "reactId"),
                 function(object, value) {
                     object@react_pos <- value
                     return(object)
                 }
)


# react_id
setMethod("react_id", signature(object = "reactId"),
          function(object) {
              return(object@react_id)
          }
)

setReplaceMethod("react_id", signature = (object = "reactId"),
                 function(object, value) {
                     object@react_id <- value
                     return(object)
                 }
)


#------------------------------------------------------------------------------#
#                               other methods                                  #
#------------------------------------------------------------------------------#

setMethod("show", signature(object = "reactId"),
    function(object) {
        posid <- cbind(react_pos(object), react_id(object))
        colnames(posid) <- c("position", "reaction_id")
        show(posid)
    }
)


# length of an object of class reactId
setMethod("length", signature = signature(x = "reactId"),
          function(x) {
              return(length(react_pos(x)))
          }
)
