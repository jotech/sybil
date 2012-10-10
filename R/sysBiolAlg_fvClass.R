#  sysBiolAlg_fvClass.R
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


#------------------------------------------------------------------------------#
#                    definition of the class sysBiolAlg_fv                     #
#------------------------------------------------------------------------------#

setClass(Class = "sysBiolAlg_fv",
         contains = "sysBiolAlg"
)


#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

# contructor for class sysBiolAlg_fv
setMethod(f = "initialize",
          signature = "sysBiolAlg_fv",
          definition = function(.Object,
                                model,
                                percentage = 100,
                                Zopt = NULL,
                                fixObjVal = TRUE,
                                tol = SYBIL_SETTINGS("TOLERANCE"),
                                lpdir = SYBIL_SETTINGS("OPT_DIRECTION"),
                                scaling = NULL, ...) {

              if ( ! missing(model) ) {

                  stopifnot(is(model, "modelorg"),
                            (is.null(Zopt) || is(Zopt, "numeric")),
                            is(tol, "numeric"),
                            is(percentage, "numeric"),
                            is(lpdir, "character"))
                  
                  # problem dimensions
                  nCols <- react_num(model)
                  nRows <- met_num(model)

                  .Object <- callNextMethod(.Object,
                                            alg        = "fv",
                                            pType      = "lp",
                                            scaling    = scaling,
                                            fi         = 1:nCols,
                                            nCols      = nCols,
                                            nRows      = nRows,
                                            mat        = S(model),
                                            ub         = uppbnd(model),
                                            lb         = lowbnd(model),
                                            obj        = rep(0, nCols),
                                            rlb        = rep(0, nRows),
                                            rtype      = rep("E", nRows),
                                            lpdir      = lpdir,
                                            rub        = NULL,
                                            ctype      = NULL,
                                            ...)

                  # objective value
                  if ( (isTRUE(fixObjVal)) && (any(obj_coef(model) != 0)) ) {
                      if (is.null(Zopt)) {
                          optimal <- optimizeProb(model,
                                                  algorithm = "fba",
                                                  lpdir = lpdir,
                                                  scaling = scaling, ...)

                          if (optimal$ok == 0) {
                              if (lpdir == "max") {
                                  obj <- sybil:::.floorValues(optimal$obj,
                                                       tol = tol)*percentage/100
                              }
                              else {
                                  obj <- sybil:::.ceilValues(optimal$obj,
                                                       tol = tol)*percentage/100
                              }
                          }
                          else {
                              stop("No optimal solution!")
                          }
                      }
                      else {
                          obj <- Zopt
                      }

                      # add a row to the problem
                      type <- ifelse(lpdir == "max", "L", "U")
                      oind <- which(obj_coef(model) != 0)
                      oval <- obj_coef(model)[oind]
                      addRowsToProb(problem(.Object),
                                    met_num(model)+1,
                                    type, obj, obj,
                                    list(oind), list(oval))
                      .Object@nr <- .Object@nr + 1L

                  }
              }
              return(.Object)
          }
)


#------------------------------------------------------------------------------#
