### R code from vignette source 'sybil.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: sybil.Rnw:429-430
###################################################
library(sybil)


###################################################
### code chunk number 2: sybil.Rnw:439-440
###################################################
library(help = "sybil")


###################################################
### code chunk number 3: sybil.Rnw:444-445
###################################################
help(doubleGeneDel)


###################################################
### code chunk number 4: sybil.Rnw:449-450
###################################################
help.search("flux variability analysis")


###################################################
### code chunk number 5: sybil.Rnw:453-454 (eval = FALSE)
###################################################
## vignette("sybil")


###################################################
### code chunk number 6: sybil.Rnw:475-476
###################################################
mp  <- system.file(package = "sybil", "extdata")


###################################################
### code chunk number 7: sybil.Rnw:479-480
###################################################
mod <- readTSVmod(prefix = "Ec_core", fpath = mp, quoteChar = "\"")


###################################################
### code chunk number 8: sybil.Rnw:500-501 (eval = FALSE)
###################################################
## modelorg2tsv(mod, prefix = "Ec_core")


###################################################
### code chunk number 9: sybil.Rnw:506-507
###################################################
data(Ec_core)


###################################################
### code chunk number 10: sybil.Rnw:520-521
###################################################
ex <- findExchReact(Ec_core)


###################################################
### code chunk number 11: sybil.Rnw:524-526
###################################################
upt <- uptReact(ex)
ex[upt]


###################################################
### code chunk number 12: sybil.Rnw:532-534
###################################################
mod <- changeBounds(Ec_core, ex[c("EX_glc(e)", "EX_lac_D(e)")], lb = c(0, -10))
findExchReact(mod)


###################################################
### code chunk number 13: sybil.Rnw:553-554
###################################################
optL <- optimizeProb(Ec_core, algorithm = "fba", retOptSol = FALSE)


###################################################
### code chunk number 14: sybil.Rnw:559-560
###################################################
opt <- optimizeProb(Ec_core, algorithm = "fba", retOptSol = TRUE)


###################################################
### code chunk number 15: sybil.Rnw:566-567
###################################################
lp_obj(opt)


###################################################
### code chunk number 16: sybil.Rnw:571-572
###################################################
checkOptSol(opt)


###################################################
### code chunk number 17: sybil.Rnw:586-587
###################################################
fba <- optimizeProb(Ec_core, algorithm = "fba")


###################################################
### code chunk number 18: sybil.Rnw:590-591
###################################################
mod_obj(fba)


###################################################
### code chunk number 19: sybil.Rnw:595-596
###################################################
mtf <- optimizeProb(Ec_core, algorithm = "mtf", wtobj = mod_obj(fba))


###################################################
### code chunk number 20: sybil.Rnw:599-600
###################################################
lp_obj(mtf)


###################################################
### code chunk number 21: sybil.Rnw:605-607
###################################################
nvar(fluxdist(fba))
nvar(fluxdist(mtf))


###################################################
### code chunk number 22: sybil.Rnw:612-614
###################################################
help("sysBiolAlg_fba-class")
help("sysBiolAlg_mtf-class")


###################################################
### code chunk number 23: sybil.Rnw:617-619
###################################################
?fba
?mtf


###################################################
### code chunk number 24: sybil.Rnw:624-626
###################################################
fl <- getFluxDist(mtf)
length(fl)


###################################################
### code chunk number 25: sybil.Rnw:631-633
###################################################
fd <- getFluxDist(mtf, ex)
getNetFlux(fd)


###################################################
### code chunk number 26: sybil.Rnw:638-639
###################################################
mod_obj(mtf)


###################################################
### code chunk number 27: sybil.Rnw:650-651
###################################################
ko <- optimizeProb(Ec_core, gene = "b2276", lb = 0, ub = 0)


###################################################
### code chunk number 28: sybil.Rnw:672-674
###################################################
ko <- optimizeProb(Ec_core, gene = "b2276", lb = 0, ub = 0,
                   algorithm = "lmoma", wtflux = getFluxDist(mtf))


###################################################
### code chunk number 29: sybil.Rnw:690-693 (eval = FALSE)
###################################################
## ko <- optimizeProb(Ec_core, gene = "b2276", lb = 0, ub = 0,
##                    algorithm = "room", wtflux = getFluxDist(mtf),
##                    solverParm = list(PRESOLVE = GLP_ON))


###################################################
### code chunk number 30: sybil.Rnw:725-726
###################################################
opt <- oneGeneDel(Ec_core)


###################################################
### code chunk number 31: sybil.Rnw:739-740
###################################################
checkOptSol(opt)


###################################################
### code chunk number 32: sybil.Rnw:744-745
###################################################
histogram(opt, col = "lightgray", nint = 20)


###################################################
### code chunk number 33: sybil.Rnw:751-754
###################################################
opt <- oneGeneDel(Ec_core, algorithm = "lmoma", wtflux = getFluxDist(mtf))
checkOptSol(opt)
histogram(opt, col = "lightgray", nint = 20)


###################################################
### code chunk number 34: sybil.Rnw:760-761
###################################################
opt <- geneDeletion(Ec_core)


###################################################
### code chunk number 35: sybil.Rnw:764-766 (eval = FALSE)
###################################################
## opt2 <- geneDeletion(Ec_core, combinations = 2)
## opt3 <- geneDeletion(Ec_core, combinations = 3)


###################################################
### code chunk number 36: sybil.Rnw:782-784
###################################################
opt <- fluxVar(Ec_core, percentage = 80, verboseMode = 0)
plot(opt)


###################################################
### code chunk number 37: sybil.Rnw:812-814
###################################################
opt <- robAna(Ec_core, "EX_o2(e)", verboseMode = 0)
plot(opt)


###################################################
### code chunk number 38: sybil.Rnw:828-829
###################################################
opt <- oneGeneDel(Ec_core, algorithm = "fba", fld = "all")


###################################################
### code chunk number 39: sybil.Rnw:832-833
###################################################
sum <- summaryOptsol(opt, Ec_core)


###################################################
### code chunk number 40: sybil.Rnw:847-848
###################################################
printExchange(sum, j = c(1:50), dense = TRUE)


###################################################
### code chunk number 41: sybil.Rnw:863-867
###################################################
ref    <- optimizeProb(Ec_core)
opt    <- oneGeneDel(Ec_core)
let    <- lethal(opt, wt = mod_obj(ref))
nletid <- c(1:length(allGenes(Ec_core)))[! let] 


###################################################
### code chunk number 42: sybil.Rnw:876-877 (eval = FALSE)
###################################################
## gmat <- combn(nletid, 3)


###################################################
### code chunk number 43: sybil.Rnw:882-883 (eval = FALSE)
###################################################
## opt <- multiDel(Ec_core, nProc = 4, todo = "geneDeletion", del1 = gmat)


###################################################
### code chunk number 44: sybil.Rnw:894-895 (eval = FALSE)
###################################################
## mapply(checkOptSol, opt)


###################################################
### code chunk number 45: sybil.Rnw:906-908
###################################################
opt <- optimizeProb(Ec_core, poCmd = list("getRedCosts"))
postProc(opt)


###################################################
### code chunk number 46: sybil.Rnw:944-945 (eval = FALSE)
###################################################
## optimizeProb(Ec_core, method = "exact")


###################################################
### code chunk number 47: sybil.Rnw:948-949 (eval = FALSE)
###################################################
## optimizeProb(Ec_core, solver = "cplexAPI", method = "dualopt")


###################################################
### code chunk number 48: sybil.Rnw:972-975 (eval = FALSE)
###################################################
## opt <- oneGeneDel(Ec_core,
##                   solverParm = list(TM_LIM = 1000,
##                                     PRESOLVE = GLP_ON))


###################################################
### code chunk number 49: sybil.Rnw:988-992 (eval = FALSE)
###################################################
## opt <- optimizeProb(Ec_core,
##                     solverParm = list(CPX_PARAM_SCRIND = CPX_ON,
##                                       CPX_PARAM_EPRHS = 1E-09),
##                     solver = "cplexAPI")


###################################################
### code chunk number 50: sybil.Rnw:1011-1015 (eval = FALSE)
###################################################
## opt <- optimizeProb(Ec_core,
##                     solverParm = list(verbose = "full",
##                                       timeout = 10),
##                     solver = "lpSolveAPI")


###################################################
### code chunk number 51: sybil.Rnw:1029-1030
###################################################
help(SYBIL_SETTINGS)


###################################################
### code chunk number 52: sybil.Rnw:1052-1053 (eval = FALSE)
###################################################
## SYBIL_SETTINGS("parameter name", value)


###################################################
### code chunk number 53: sybil.Rnw:1058-1059 (eval = FALSE)
###################################################
## SYBIL_SETTINGS("parameter name")


###################################################
### code chunk number 54: sybil.Rnw:1064-1065 (eval = FALSE)
###################################################
## SYBIL_SETTINGS()


###################################################
### code chunk number 55: sybil.Rnw:1080-1081
###################################################
SYBIL_SETTINGS("SOLVER", "cplexAPI", loadPackage = FALSE)


###################################################
### code chunk number 56: sybil.Rnw:1087-1088
###################################################
SYBIL_SETTINGS("METHOD")


###################################################
### code chunk number 57: sybil.Rnw:1091-1092
###################################################
SYBIL_SETTINGS("SOLVER", "glpkAPI")


###################################################
### code chunk number 58: sybil.Rnw:1095-1096
###################################################
SYBIL_SETTINGS("METHOD")


###################################################
### code chunk number 59: sybil.Rnw:1131-1133
###################################################
data(Ec_core)
Ec_core


###################################################
### code chunk number 60: sybil.Rnw:1137-1138
###################################################
help("modelorg")


###################################################
### code chunk number 61: sybil.Rnw:1145-1146
###################################################
react_num(Ec_core)


###################################################
### code chunk number 62: sybil.Rnw:1149-1150
###################################################
id <- react_id(Ec_core)


###################################################
### code chunk number 63: sybil.Rnw:1153-1154
###################################################
react_id(Ec_core)[13] <- "biomass"


###################################################
### code chunk number 64: sybil.Rnw:1158-1160
###################################################
cg <- gray(0:8/8)
image(S(Ec_core), col.regions = c(cg, rev(cg)))


###################################################
### code chunk number 65: sybil.Rnw:1173-1174 (eval = FALSE)
###################################################
## mod <- readTSVmod(reactList = "reactionList.txt")


###################################################
### code chunk number 66: sybil.Rnw:1198-1199
###################################################
help("optsol")


###################################################
### code chunk number 67: sybil.Rnw:1203-1205
###################################################
os <- optimizeProb(Ec_core)
is(os)


###################################################
### code chunk number 68: sybil.Rnw:1208-1209
###################################################
lp_obj(os)


###################################################
### code chunk number 69: sybil.Rnw:1212-1213
###################################################
getFluxDist(os)


###################################################
### code chunk number 70: sybil.Rnw:1258-1259
###################################################
lp <- optObj(solver = "glpkAPI", method = "exact")


###################################################
### code chunk number 71: sybil.Rnw:1283-1284
###################################################
lp <- initProb(lp)


###################################################
### code chunk number 72: sybil.Rnw:1288-1293
###################################################
cm <- Matrix(c(0.5, 2, 1, 1), nrow = 2)
loadLPprob(lp, nCols = 2, nRows = 2, mat = cm,
           lb = c(0, 0), ub = rep(1000, 2), obj = c(1, 1),
           rlb = c(0, 0), rub = c(4.5, 9), rtype = c("U", "U"),
           lpdir = "max")


###################################################
### code chunk number 73: sybil.Rnw:1298-1299
###################################################
lp


###################################################
### code chunk number 74: sybil.Rnw:1302-1303
###################################################
status <- solveLp(lp)


###################################################
### code chunk number 75: sybil.Rnw:1306-1307
###################################################
getMeanReturn(code = status, solver = solver(lp))


###################################################
### code chunk number 76: sybil.Rnw:1310-1312
###################################################
status <- getSolStat(lp)
getMeanStatus(code = status, solver = solver(lp))


###################################################
### code chunk number 77: sybil.Rnw:1316-1318
###################################################
getObjVal(lp)
getFluxDist(lp)


###################################################
### code chunk number 78: sybil.Rnw:1321-1322
###################################################
getRedCosts(lp)


###################################################
### code chunk number 79: sybil.Rnw:1325-1327
###################################################
delProb(lp)
lp


###################################################
### code chunk number 80: sybil.Rnw:1358-1360
###################################################
ec <- sysBiolAlg(Ec_core, algorithm = "fba")
is(ec)


###################################################
### code chunk number 81: sybil.Rnw:1365-1366
###################################################
opt <- optimizeProb(ec)


###################################################
### code chunk number 82: sybil.Rnw:1373-1376
###################################################
ecr <- sysBiolAlg(Ec_core, algorithm = "room", wtflux = opt$fluxes)
is(ecr)
ecr


###################################################
### code chunk number 83: sybil.Rnw:1425-1426 (eval = FALSE)
###################################################
## promptSysBiolAlg(algorithm = "foo")


