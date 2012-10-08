#  optsol_optimizeProbClass.R
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


# optsol_optimizeProbClass


#------------------------------------------------------------------------------#
#               definition of the class optsol_optimizeProb                    #
#------------------------------------------------------------------------------#

setClass("optsol_optimizeProb",
         representation(
              preProc  = "ppProc", # preprocessing lp result
              postProc = "ppProc"  # postprocessing lp result
        ),
        contains = "optsol"
)


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# preProc
setMethod("preProc", signature(object = "optsol_optimizeProb"),
          function(object) {
              return(object@preProc)
          }
)

setReplaceMethod("preProc", signature = (object = "optsol_optimizeProb"),
                 function(object, value) {
                     object@preProc <- value
                     return(object)
                 }
)


# postProc
setMethod("postProc", signature(object = "optsol_optimizeProb"),
          function(object) {
              return(object@postProc)
          }
)

setReplaceMethod("postProc", signature = (object = "optsol_optimizeProb"),
                 function(object, value) {
                     object@postProc <- value
                     return(object)
                 }
)
