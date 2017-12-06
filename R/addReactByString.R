library(stringr)


addReactByString <- function(model, 
                             ids, 
                             react_str, 
                             lb=0, 
                             ub=SYBIL_SETTINGS("MAXIMUM"), 
                             obj=0, 
                             subSystem=NA, 
                             gprAssoc=NA, 
                             reactName=NA, 
                             metName=NA, 
                             metComp=NA) {
  
  # ------------------------------------------------------------------------ #
  # check arguments
  # ------------------------------------------------------------------------ #
  
  if (!is(model, "modelorg")) {
    stop("needs an object of class modelorg!")
  }
  
  stopifnot(checkVersion(model))
  
  Nids <- length(ids)
  if (length(react_str) != Nids){
    stop("all arguments have to be provided for each reaction")
  }

  mets  <- list()
  Scoef <- list()
  rev   <- vector(mode="logical", length=Nids)
  for( i in 1:Nids){
    str    <- react_str[i]
    dirtmp <- list("<==>"="BOTH", "<-->"="BOTH", "<=>"="BOTH", "<->"="BOTH",
                   "==>"="FORW", "-->"="FORW", "=>"="FORW", "->"="FORW",
                   "<=="="BACK", "<--"="BACK", "<="="BACK", "<-"="BACK") # TODO detection may cause problem, order of list is important (e.g. <== before <=)
    dir    <- str_extract(str, paste(names(dirtmp), collapse="|"))
    rev[i] <- unname(dirtmp[dir]=="BOTH")
    
    split   <- unlist(str_split(str, dir))
    edu_str <- str_strip(split[1])
    pro_str <- str_strip(split[2])
    
    edu     <- get_reactants(edu_str)
    pro     <- get_reactants(pro_str)
    
    Scoef    <- c(Scoef, list(c(-1 * as.numeric(edu[2,]), as.numeric(pro[2,]))) )
    mets     <- c(mets,  list(unname(c(edu[1,], pro[1,]))) )
  }
  
  lb <- ifelse(rev==TRUE, -1000, 0) # always the case?
  model <- addMultiReact(model=model, ids=ids, mets=mets, Scoefs=Scoef, reversible=rev, lb=lb, ub=ub, obj=obj, 
                subSystem=subSystem, gprAssoc=gprAssoc, reactName=reactName, metName=metName, metComp=metComp)
  
  return(model)
}

get_reactants <- function(str){
  split <- unlist(str_split(str, "\\+"))
  stoich<- sapply(split, function(s){
    s <- str_strip(s)
    scoef <- str_strip(str_extract(s, "\\b\\(?[0-9,.]+\\)?\\b"))
    if(is.na(scoef)) {
      scoef    <- 1.0
      met      <- str_strip(s)
    } else {
      met <- str_strip(gsub(scoef,"",s))
      met <- str_strip(gsub("\\(\\)","",met)) # TODO: needed numbers are given as '(1)' and not as '1'
    }
    c(met, scoef)
  })
  return(stoich)
}

str_strip <- function(str){
  return(gsub("^\\s+","",gsub("\\s+$","",str)))
}