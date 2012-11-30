#  optimizer.R
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
# Function: optimizer
#
#
#


optimizer <- function(model, lb, ub,
                      delete, geneFlag,
                      algorithm = SYBIL_SETTINGS("ALGORITHM"),
                      setToZero = FALSE,
                      checkOptSolObj = FALSE,
                      rebuildModel = FALSE,
                      fld = "none",
                      prCmd = NA, poCmd = NA,
                      prDIR = NULL, poDIR = NULL,
                      verboseMode = 2,
                      ...) {


    stopifnot(length(fld) == 1)

    #--------------------------------------------------------------------------#
    # verboseMode

    on.exit(expr = {
        if (exists("logObj")) {
            logClose(logObj) <- NA
        }
    } )

    # start logging
    logObj <- sybilLog(filename = "",
                       loglevel = -1,
                       verblevel = verboseMode)


    #--------------------------------------------------------------------------#
    # check arguments

    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }

    if (!is(delete, "matrix")) {
        stop("needs an object of class matrix!")
    }

    if (missing(lb)) {
        lb <- rep(0, nrow(delete))
    }

    if (missing(ub)) {
        ub <- rep(0, nrow(delete))
    }

    if (!checkAlgorithm(alg = algorithm, fkt = "optimizer")) {
        stop(sQuote(algorithm), " is not a valid algorithm")
    }

    if (isTRUE(fld)) {
        fld <- "all"
    }
    if (identical(fld, FALSE)) {
        fld <- "none"
    }

#    if ((isTRUE(rebuildModel)) && (isTRUE(copyModel))) {
#        stop("use 'rebuildModel' xor 'copyModel'")
#    }


    #--------------------------------------------------------------------------#
    # prepare problem object
    #--------------------------------------------------------------------------#

    lpmod <- sysBiolAlg(model, algorithm = algorithm, ...)


    #--------------------------------------------------------------------------#
    # data structures for simulation results
    #--------------------------------------------------------------------------#

    nObj  <- nrow(delete)
    obj   <- numeric(nObj)
    mobj  <- numeric(nObj)
    ok    <- integer(nObj)
    stat  <- integer(nObj)
    flux <- switch(fld,
        "all" = {
            Matrix::Matrix(0, nrow = nc(lpmod), ncol = nObj)
        },
        "fluxes" = {
            Matrix::Matrix(0, nrow = length(fldind(lpmod)), ncol = nObj)
        },
        {
            NA
        }
    )
#     if (isTRUE(fld)) {
#         flux <- Matrix::Matrix(0, nrow = nc(lpmod), ncol = nObj)
#     }
#     else {
#         flux <- NA
#     }

    if (isTRUE(geneFlag)) {
        heff  <- logical(nObj)
        fdels <- vector("list", nObj)
    }

    runPrPl  <- logical(nObj)
    runPoPl  <- logical(nObj)
    runPrPcn <- 1
    runPoPcn <- 1

    if (all(!is.na(prCmd))) {
        do_pr  <- TRUE
        prPcmd <- NULL
        runPrP <- sybil:::.doInRound(prDIR, nObj)
        prPpa  <- vector(mode = "list", length = length(runPrP))
        runPrPl[runPrP] <- TRUE
    }
    else {
        do_pr <- FALSE
    }
    if (all(!is.na(poCmd))) {
        do_po  <- TRUE
        poPcmd <- NULL
        runPoP <- sybil:::.doInRound(poDIR, nObj)
        poPpa  <- vector(mode = "list", length = length(runPoP))
        runPoPl[runPoP] <- TRUE
    }
    else {
        do_po <- FALSE
    }


#------------------------------------------------------------------------------#
#                    flux/gene deletions with linearMOMA                       #
#------------------------------------------------------------------------------#

    message("calculating ", nObj, " optimizations ... ", appendLF = FALSE)
    if (verboseMode > 1) { cat("\n") }
    if (verboseMode == 2) {
        progr <- sybil:::.progressBar()
        #progr <- txtProgressBar(min = 2, max = nObj, initial = 2, style = 3)
    }

    logComment(logObj, "compute mutant strains")
    logOptimizationTH(logObj)


#------------------------------------------------------------------------------#
#                        flux/gene deletions with FBA                          #
#------------------------------------------------------------------------------#

    fi <- fldind(lpmod)

    for (i in 1:nObj) {

        if (verboseMode == 2) {
            progr <- sybil:::.progressBar(i, nObj, progr)
            #setTxtProgressBar(progr, i)
        }

        if (isTRUE(geneFlag)) {
            # get the reactions for gene i
            tmp_del <- geneDel(model, delete[i, ])

            if (is.null(tmp_del)) {
                fdels[[i]] <- NA
            }
            else {
                # deletion of gene i has an effect (tmp_del is not empty)
                heff[i]    <- TRUE
                fdels[[i]] <- tmp_del
            }
        }
        else {
            tmp_del <- delete[i, ]
        }

        # pre/post processing
        if (isTRUE(runPrPl[i])) {
            prCmd_tmp <- prCmd
            did_pr    <- TRUE
        }
        else {
            prCmd_tmp <- NA
            did_pr    <- FALSE
        }

        if (isTRUE(runPoPl[i])) {
            poCmd_tmp <- poCmd
            did_po    <- TRUE
        }
        else {
            poCmd_tmp <- NA
            did_po    <- FALSE
        }

        # solution i
        if (isTRUE(rebuildModel)) {
            sol <- optimizeProb(model,
                                retOptSol = FALSE,
                                react = tmp_del,
                                lb = rep(lb[i], length(tmp_del)),
                                ub = rep(ub[i], length(tmp_del)),
                                algorithm = algorithm,
                                MoreArgs = list(lpdir = getObjDir(problem(lpmod)),
                                                prCmd = prCmd_tmp,
                                                poCmd = poCmd_tmp,
                                                prCil = runPrPcn,
                                                poCil = runPoPcn),
                                ...)
        }
        else {
            sol <- optimizeProb(lpmod,
                                react = tmp_del,
                                lb = rep(lb[i], length(tmp_del)),
                                ub = rep(ub[i], length(tmp_del)),
                                prCmd = prCmd_tmp, poCmd = poCmd_tmp,
                                prCil = runPrPcn, poCil = runPoPcn)
        }


        #obj[i]  <- sum(obj_coef(model) * sol$fluxes[fi])
        #mobj[i] <- crossprod(obj_coef(model), sol$fluxes[fi])
        #obj[i]  <- sol$obj
        ok[i]   <- sol$ok
        stat[i] <- sol$stat
        if (fld == "none") {
            obj[i] <- crossprod(obj_coef(model), sol$fluxes[fi])
        }
        else {
            obj[i] <- sol$obj
            if (fld == "fluxes") {
                flux[,i] <- sol$fluxes[fi]
            }
            else {
                flux[,i] <- sol$fluxes
            }
        }
#         if (isTRUE(fld)) {
#             flux[,i] <- sol$fluxes
#             obj[i]  <- sol$obj
#         }
#         else {
#             obj[i] <- crossprod(obj_coef(model), sol$fluxes[fi])
#         }

        # pre/post processing
        if (isTRUE(did_pr)) {
            if ( (runPrPcn == 1) && (is.null(prPcmd)) ) {
                prPcmd  <- cmd(sol$preP)
            }
            prPpa[[runPrPcn]] <- pa(sol$preP)
            runPrPcn <- runPrPcn+1
            did_pr <- FALSE
        }
        if (isTRUE(did_po)) {
            if ( (runPoPcn == 1) && (is.null(poPcmd)) ) {
                poPcmd  <- cmd(sol$postP)
            }
            poPpa[[runPoPcn]] <- pa(sol$postP)
            runPoPcn <- runPoPcn+1
            did_po <- FALSE
        }

        logOptimization(logObj, sol$ok, sol$stat, obj[i], tmp_del, i)

        remove(sol)
        #close(progr)
    }

    message("OK")


#------------------------------------------------------------------------------#
#                               save the results                               #
#------------------------------------------------------------------------------#

    if (fld == "fluxes") {
        fli <- 1:length(fi)
    }
    else if(fld == "none") {
        fli <- NA
    }
    else {
        fli <- fi
    }

    if (isTRUE(geneFlag)) {
        optsol <- new("optsol_genedel",
            mod_id       = mod_id(model),
            mod_key      = mod_key(model),
            solver       = solver(problem(lpmod)),
            method       = method(problem(lpmod)),
            algorithm    = algorithm(lpmod),
            num_of_prob  = as.integer(nObj),
            lp_num_cols  = nc(lpmod),
            lp_num_rows  = nr(lpmod),
            lp_obj       = as.numeric(obj),
            lp_ok        = as.integer(ok),
            lp_stat      = as.integer(stat),
            lp_dir       = factor(getObjDir(problem(lpmod))),
            obj_coef     = obj_coef(model),
            obj_func     = printObjFunc(model),
            fldind       = as.integer(fli),
            fluxdist     = fluxDistribution(flux),
    
            chlb         = as.numeric(lb),
            chub         = as.numeric(ub),
            dels         = delete,
    
            fluxdels     = fdels,
            hasEffect    = heff
        )
    }
    else {
        optsol <- new("optsol_fluxdel",
            mod_id       = mod_id(model),
            mod_key      = mod_key(model),
            solver       = solver(problem(lpmod)),
            method       = method(problem(lpmod)),
            algorithm    = algorithm(lpmod),
            num_of_prob  = as.integer(nObj),
            lp_num_cols  = nc(lpmod),
            lp_num_rows  = nr(lpmod),
            lp_obj       = as.numeric(obj),
            lp_ok        = as.integer(ok),
            lp_stat      = as.integer(stat),
            lp_dir       = factor(getObjDir(problem(lpmod))),
            obj_coef     = obj_coef(model),
            obj_func     = printObjFunc(model),
            fldind       = as.integer(fli),
            fluxdist     = fluxDistribution(flux),
    
            chlb         = as.numeric(lb),
            chub         = as.numeric(ub),
            dels         = delete
        )
    }    


    if (isTRUE(do_pr)) {
        prAna <- ppProc(prPcmd)
        pa(prAna) <- prPpa
        ind(prAna) <- runPrP
        preProc(optsol) <- prAna
    }

    if (isTRUE(do_po)) {
        poAna <- ppProc(poPcmd)
        pa(poAna) <- poPpa
        ind(poAna) <- runPoP
        postProc(optsol) <- poAna
    }


#------------------------------------------------------------------------------#
#                           return solution object                             #
#------------------------------------------------------------------------------#


    if (isTRUE(setToZero)) {
        do_again <- checkStat(optsol)
        num_new  <- length(do_again)
        lp_obj(optsol)[do_again] <- as.numeric(0)

        message("setting ", num_new, " objective values to zero")

        for (i in seq(along = do_again)) {
            logOptimization(logObj, lp_ok(optsol)[do_again[i]],
                            lp_stat(optsol)[do_again[i]], 0,
                            deleted(optsol, do_again[i]), do_again[i])
        }
    }


    if (isTRUE(checkOptSolObj)) {
        checkOptSol(optsol, onlywarn = TRUE)
    }

    delProb(problem(lpmod))
    remove(lpmod)

    logFoot(logObj)  <- TRUE
    logClose(logObj) <- NA

    return(optsol)

}



