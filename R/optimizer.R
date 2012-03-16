#  optimizer.R
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
# Function: optimizer
#
# Parameters:
#       todo:   what to do
#     solver:   lp solver
#
# Return values:
#     on derivative of the class optsol


optimizer <- function(model, optsol, lb, ub,
                      wtFluxes = NA,
                      alg = SYBIL_SETTINGS("ALGORITHM"),
                      setToZero = FALSE,
                      checkOptSolObj = FALSE,
                      rebuildModel = FALSE,
                      copyModel = FALSE,
                      prCmd = NA, poCmd = NA,
                      prDIR = NA, poDIR = NA,
                      verboseMode = 3,
                      loglevel = -1,
                      logfile = NA,
                      logfileEnc = NA,
                      ...) {


    #--------------------------------------------------------------------------#
    # logging

    on.exit(expr = {
        if (exists("logObj")) {
            logClose(logObj) <- NA
        }
    } )

    # start logging
    logObj <- sybilLog(filename = logfile,
                       loglevel = loglevel,
                       verblevel = verboseMode,
                       logfileEnc = logfileEnc)

    logHead(logObj)
    logCall(logObj, nog = 2)
    logComment(logObj, "input/output")
    logComment(logObj, paste("\tmodel id:", mod_id(model), sep = "\t"))
    logComment(logObj, paste("\tclass:", is(optsol)[1], sep = "\t\t"))


    #--------------------------------------------------------------------------#

    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }

    if (!is(optsol, "optsol")) {
        stop("needs an object of class optsol!")
    }

    if (missing(lb)) {
        lb <- rep(0, nrow(dels(optsol)))
    }

    if (missing(ub)) {
        ub <- rep(0, nrow(dels(optsol)))
    }

    if ((isTRUE(rebuildModel)) && (alg == "linearMOMA_COBRA")) {
        msg <- paste("flag 'rebuildModel' is not implemented for COBRA version",
                     "of linearMOMA switching to 'linearMOMA'")
        logWarning(logObj, msg)
        alg <- "linearMOMA"
    }

    if ((isTRUE(rebuildModel)) && (isTRUE(copyModel))) {
        stop("use 'rebuildModel' xor 'copyModel'")
    }

    if (all(is.na(fluxes(optsol)))) {
        fld <- FALSE
    }
    else {
        fld <- TRUE
    }

#     MOMAflag <- FALSE
#     if (!is.na(pmatch(alg, c("linearMOMA", "linearMOMA_COBRA")))) {
#         MOMAflag <- TRUE
#     }


    # nc    is the number of reactions (columns) in the model
    # nr    is the number of metabolites (rows) in the model
    # nCols is the number of columns for the problem object
    # nRows is the number of rows for the problem object
    nc <- react_num(model)
    nr <- met_num(model)
    nCols <- nc
    nRows <- nr

    switch(alg,
	       "FBA" = {
		       FBAflag   <- TRUE
			   MOMAflag  <- FALSE
			   lMOMAflag <- FALSE
			   COBRAflag <- FALSE
               ALG_OK    <- TRUE
		   },

		   "linearMOMA" = {
		       FBAflag   <- FALSE
			   MOMAflag  <- TRUE
			   lMOMAflag <- TRUE
			   COBRAflag <- FALSE
               nCols <- 4*nc
               nRows <- nr + 2*nc
               ALG_OK    <- TRUE
		   },

		   "linearMOMA_COBRA" = {
		       FBAflag   <- FALSE
			   MOMAflag  <- TRUE
			   lMOMAflag <- FALSE
			   COBRAflag <- TRUE
               nCols <- 4*nc
               nRows <- 2*nr + 2*nc + 1
               ALG_OK    <- TRUE
		   },

           {   ALG_OK    <- FALSE   }
	)

    if (!isTRUE(ALG_OK)) {
        stop("argument alg must be 'FBA', 'linearMOMA' or 'linearMOMA_COBRA'")
    }

#      print(FBAflag)
#      print(MOMAflag)
#      print(lMOMAflag)
#      print(COBRAflag)
#      print(nCols)
#      print(nRows)


    #--------------------------------------------------------------------------#

    # check, weather we want to use the default method
    method <- checkDefaultMethod(solver(optsol), method(optsol)[1])
    method(optsol) <- rep(method$met, num_of_prob(optsol))

    # add used algorithm to solution object
    if(is(optsol, "optsol_fluxdel")) {
        algorithm(optsol) <- alg
    }

#     # prepare problem object for the wild type solution
#     # (algorith here is always FBA)
#     lpmod <- prepProbObj(model,
#                          nCols = nc,
#                          nRows = nr,
#                          alg = "FBA",
#                          lpdir = lp_dir(optsol),
#                          solver = solver(optsol),
#                          method = method(optsol)[1],
#                          ...
#                         )

    logComment(logObj, "verifyed arguments are")
    logComment(logObj, paste("\talgorithm:", algorithm(optsol), sep = "\t"))
    logComment(logObj, paste("\tsolver:", solver(optsol), sep = "\t\t"))
    logComment(logObj, paste("\tmethod:", method(optsol)[1], sep = "\t\t"))
    

    # tmp vectors for solutions
    lmfld <- TRUE
    nObj  <- num_of_prob(optsol)
    obj   <- numeric(nObj)
    ok    <- integer(nObj)
    stat  <- integer(nObj)
    if ((isTRUE(fld)) || (isTRUE(lMOMAflag))) {
        flux <- Matrix(0, nrow = nc, ncol = nObj)
    }
    else {
        flux <- Matrix(0, nrow = 1, ncol = 1)
    }
    if (is(optsol, "optsol_genedel")) {
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
        runPrP <- .doInRound(prDIR, nObj)
        prPpa  <- vector(mode = "list", length = length(runPrP))
        runPrPl[runPrP] <- TRUE
    }
    else {
        do_pr <- FALSE
    }
    if (all(!is.na(poCmd))) {
        do_po  <- TRUE
        poPcmd <- NULL
        runPoP <- .doInRound(poDIR, nObj)
        poPpa  <- vector(mode = "list", length = length(runPoP))
        runPoPl[runPoP] <- TRUE
    }
    else {
        do_po <- FALSE
    }


#------------------------------------------------------------------------------#
#                        calculate wild type solution                          #
#------------------------------------------------------------------------------#

    # check if we have an user defined flux distribution for MOMA
    if ((isTRUE(MOMAflag)) && (!any(is.na(wtFluxes)))) {

        if (!is(wtFluxes, "numeric")) {
            stop("Argument 'wtFluxes' must be numeric!")
        }

        if ((isTRUE(COBRAflag)) && (length(wtFluxes) > 1)) {
            stop("Argument 'wtFluxes' must contain exactly one element!")
        }

        if ((isTRUE(lMOMAflag)) && (length(wtFluxes) != nc)) {
            stop("Length of argument 'wtFluxes' must be the same as the number of reactions!")
        }

        obj[1]   <- ifelse(isTRUE(COBRAflag), wtFluxes, NA)
        flux[,1] <- ifelse(isTRUE(lMOMAflag), wtFluxes, NA)
        #lp_ok(optsol)[1]   <- NA
        #lp_stat(optsol)[1] <- NA
        ok[1] <- NA
        stat[1] <- NA

        logComment(logObj, "got wild type solution (wtFluxes)")

    }
    else {

        logComment(logObj, "reference solution")
        logOptimizationTH(logObj)

        if (isTRUE(lMOMAflag)) {
            oldFLD <- fld
            fld    <- TRUE
        }

#         sol <- simpleFBA(lpmod,
#                          lpdir = lp_dir(optsol),
#                          method = method(optsol)[1],
#                          solver = solver(optsol),
#                          fld = fld
#                         )

        # pre/post processing
        if (isTRUE(runPrPl[1])) {
            prCmd_tmp <- prCmd
            did_pr    <- TRUE
        }
        else {
            prCmd_tmp <- NA
            did_pr    <- FALSE
        }

        if (isTRUE(runPoPl[1])) {
            poCmd_tmp <- poCmd
            did_po    <- TRUE
        }
        else {
            poCmd_tmp <- NA
            did_po    <- FALSE
        }

        sol <- simpleFBA(model,
                         lpdir = lp_dir(optsol),
                         solver = solver(optsol),
                         method = method(optsol)[1],
                         fld = fld,
                         prCmd = prCmd_tmp, poCmd = poCmd_tmp,
                         prCil = runPrPcn, poCil = runPoPcn,
                         ...
                        )

        #lp_obj(optsol)[1]  <- sol$obj
        #lp_ok(optsol)[1]   <- sol$ok
        #lp_stat(optsol)[1] <- sol$stat
        obj[1]  <- sol$obj
        ok[1]   <- sol$ok
        stat[1] <- sol$stat
        if (isTRUE(fld)) {
            flux[,1] <- sol$fluxes
        }

        if (isTRUE(lMOMAflag)) {
            fld    <- oldFLD
        }

        # pre/post processing
        if (isTRUE(did_pr)) {
            prPcmd     <- cmd(sol$preP)
            prPpa[[1]] <- pa(sol$preP)
            runPrPcn   <- runPrPcn+1
            did_pr     <- FALSE
        }
        if (isTRUE(did_po)) {
            poPcmd     <- cmd(sol$postP)
            poPpa[[1]] <- pa(sol$postP)
            runPoPcn   <- runPoPcn+1
            did_po     <- FALSE
        }

        logOptimization(logObj, sol$ok, sol$stat, sol$obj, NA)

    }


#------------------------------------------------------------------------------#
#                    flux/gene deletions with linearMOMA                       #
#------------------------------------------------------------------------------#

    # fluxes/genes to delete
    delete <- dels(optsol)

    chlb(optsol) <- c(NA, lb)
    chub(optsol) <- c(NA, ub)

    logStep(logObj) <- paste("calculating",
                             num_of_prob(optsol),
                             "optimizations")
    if (verboseMode > 2) {
        cat("\n")
        progr <- .progressBar()
    }

    logOptimizationTH(logObj)

    if (isTRUE(MOMAflag)) {

        if (isTRUE(COBRAflag)) {
            #wtFluxes <- lp_obj(optsol)[1]
            wtFluxes <- obj[1]
        }
        else {
            wtFluxes <- flux[,1]
        }

        # prepare problem object for use with linearMOMA
        if (!isTRUE(rebuildModel)) {
            lpmod <- prepProbObj(model,
                                 #alg = alg,
                                 nCols = nCols,
                                 nRows = nRows,
                                 wtflux = wtFluxes,
                                 MOMAflag = MOMAflag,
                                 COBRAflag = COBRAflag,
                                 lpdir = lp_dir(optsol),
                                 solver = solver(optsol),
                                 method = method(optsol)[1],
                                 ...
                                )
        }

    }
    else {

        # prepare problem object for use with FBA
        if (!isTRUE(rebuildModel)) {
            lpmod <- prepProbObj(model,
                                 #alg = "FBA",
                                 nCols = nc,
                                 nRows = nr,
                                 lpdir = lp_dir(optsol),
                                 solver = solver(optsol),
                                 method = method(optsol)[1],
                                 ...
                                )
        }

    }


#------------------------------------------------------------------------------#
#                        flux/gene deletions with FBA                          #
#------------------------------------------------------------------------------#

    for (i in 2:nObj) {

        if (verboseMode > 2) { progr <- .progressBar(i, nObj, progr) }

        if (is(optsol, "optsol_genedel")) {
            # get the reactions for gene i
            #tmp_del <- geneDel(model, allGenes(model)[delete[i,]], checkId = TRUE)
            tmp_del <- geneDel(model, delete[i, ])

            if (any(is.na(tmp_del))) {
            #if (any(is.na(tmp_del)) && (!isTRUE(MOMAflag))) {

                # deletion of gene i has no effect (no flux has to set to zero)
                #fluxdels(optsol)[[i]] <- NA
                fdels[[i]] <- NA

                # copy reference solution
                if ((isTRUE(MOMAflag)) && (isTRUE(fld))) {
                    #lp_obj(optsol)[i] <- as.numeric(0)
                    obj[i] <- as.numeric(0)
                }
                else {
                    #lp_obj(optsol)[i] <- lp_obj(optsol)[1]
                    obj[i] <- obj[1]
                }
                #lp_ok(optsol)[i]   <- lp_ok(optsol)[1]
                #lp_stat(optsol)[i] <- lp_stat(optsol)[1]
                ok[i]   <- ok[1]
                stat[i] <- stat[1]
                if (isTRUE(fld)) {
                    flux[,i] <- flux[,1]
                }

                logOptimizationNE(logObj, tmp_del)
                
                if (isTRUE(runPrPl[i])) {
                    runPrPcn <- runPrPcn+1
                }
                if (isTRUE(runPoPl[i])) {
                    runPoPcn <- runPoPcn+1
                }

                next

            }
            else {

                # deletion of gene i has an effect (tmp_del is not empty)
                #hasEffect(optsol)[i]  <- TRUE
                #fluxdels(optsol)[[i]] <- tmp_del
                heff[i]    <- TRUE
                fdels[[i]] <- tmp_del
            }
        }
        else {
            tmp_del <- delete[i,]
        }

        if (isTRUE(MOMAflag)) {
            if (!isTRUE(rebuildModel)) {
                tmp_del <- tmp_del + nc
            }
            if (!isTRUE(fld)) {
                # if fld is TRUE, we want lp_obj to be the value
                # of the objective of linear MOMA, otherwise
                # we want it to be the value of the flux in the
                # wild type
                fld <- TRUE
                lmfld <- FALSE
            }
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
        if (!isTRUE(rebuildModel)) {
            if (isTRUE(copyModel)) {
                lpmod_backup <- backupProb(lpmod)
                lpmod_modify <- lpmod_backup
            }
            else {
                lpmod_modify <- lpmod
            }

            sol <- simpleFBA(lpmod_modify,
                             react = tmp_del,
                             lb = rep(lb[i], length(tmp_del)),
                             ub = rep(ub[i], length(tmp_del)),
                             lpdir = lp_dir(optsol),
                             solver = solver(optsol),
                             method = method(optsol)[i],
                             fld = fld,
                             prCmd = prCmd_tmp, poCmd = poCmd_tmp,
                             prCil = runPrPcn, poCil = runPoPcn
                   )


            if (isTRUE(copyModel)) {
                delProb(lpmod_modify, closeEnv = FALSE)
                remove(lpmod_modify)
            }
        }
        else {
            mD <- ifelse(isTRUE(lMOMAflag), TRUE, FALSE)

            sol <- simpleFBA(model,
                             react = tmp_del,
                             lb = rep(lb[i], length(tmp_del)),
                             ub = rep(ub[i], length(tmp_del)),
                             lpdir = lp_dir(optsol),
                             minDist = mD,
                             wtFluxes = wtFluxes,
                             solver = solver(optsol),
                             method = method(optsol)[1],
                             fld = fld,
                             prCmd = prCmd_tmp, poCmd = poCmd_tmp,
                             prCil = runPrPcn, poCil = runPoPcn,
                             ...
                   )

        }

        if ((isTRUE(MOMAflag)) && (!isTRUE(lmfld))) {
            #lp_obj(optsol)[i] <- sum(obj_coef(model) * sol$fluxes[(nc+1):(2*nc)])
            if (isTRUE(rebuildModel)) {
                obj[i] <- sum(obj_coef(model) * sol$fluxes)
            }
            else {
                obj[i] <- sum(obj_coef(model) * sol$fluxes[(nc+1):(2*nc)])
            }
            fld <- FALSE
        }
        else {
            #lp_obj(optsol)[i] <- sol$obj
            obj[i] <- sol$obj
        }
        #lp_ok(optsol)[i]   <- sol$ok
        #lp_stat(optsol)[i] <- sol$stat
        ok[i]   <- sol$ok
        stat[i] <- sol$stat
        if (isTRUE(fld)) {
            if ((isTRUE(MOMAflag)) && (!isTRUE(rebuildModel))) {
                flux[,i] <- sol$fluxes[(nc+1):(2*nc)]
                #flux[,i] <- sol$fluxes[1:nc] # wild type flux distribution
            }
            else {
                flux[,i] <- sol$fluxes
            }
        }

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

        if (isTRUE(MOMAflag)) {
            tmp_del <- tmp_del - nc
        }
        logOptimization(logObj, sol$ok, sol$stat, obj[i], tmp_del)

        remove(sol)
    }

    logStep(logObj) <- NA


#------------------------------------------------------------------------------#
#                               save the results                               #
#------------------------------------------------------------------------------#

    lp_obj(optsol)    <- as.numeric(obj)
    lp_ok(optsol)     <- as.integer(ok)
    lp_stat(optsol)   <- as.integer(stat)
    remove(obj)
    remove(ok)
    remove(stat)
    if (isTRUE(fld)) {
        fluxes(optsol) <- flux
        remove(flux)
    }
    if (is(optsol, "optsol_genedel")) {
        hasEffect(optsol) <- heff
        fluxdels(optsol)  <- fdels
        remove(heff)
        remove(fdels)
    }

    if (isTRUE(do_pr)) {
        prAna <- ppProc(prPcmd)
        pa(prAna) <- prPpa
        ind(prAna) <- runPrP
        preProc(optsol) <- prAna
        remove(prPpa)
        remove(runPrP)
    }
    if (isTRUE(do_po)) {
        poAna <- ppProc(poPcmd)
        pa(poAna) <- poPpa
        ind(poAna) <- runPoP
        postProc(optsol) <- poAna
        remove(poPpa)
        remove(runPoP)
    }


#------------------------------------------------------------------------------#
#                           return solution object                             #
#------------------------------------------------------------------------------#

    if (isTRUE(setToZero)) {
        do_again <- checkStat(optsol)
        num_new  <- length(do_again)
        lp_obj(optsol)[do_again] <- as.numeric(0)

        logStep(logObj) <- paste("setting", num_new, "objective values to zero")
        logHead(logObj)
        for (i in seq(along = do_again)) {
            logOptimization(logObj, lp_ok(optsol)[do_again[i]],
                            lp_stat(optsol)[do_again[i]], 0,
                            deleted(optsol, do_again[i]))
        }
    }


    if (isTRUE(checkOptSolObj)) {
        checkOptSol(optsol, onlywarn = TRUE)
    }

    if (!isTRUE(rebuildModel)) {
        delProb(lpmod)
        remove(lpmod)
    }

    remove(runPrPl)
    remove(runPoPl)

    logFoot(logObj)  <- TRUE
    logClose(logObj) <- NA

    return(optsol)

}



