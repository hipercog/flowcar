readCCSdb <- function(fname, table = c("blob", "run", "step"), tag = ""){
  require(RSQLite)
  
  t <- c("blob", "run", "step")
  retval <- list()
  
  sqlite.driver <- dbDriver("SQLite")
  connection <- dbConnect(sqlite.driver, dbname = fname)
  
  # get tables named blob, run, step
  tables <- dbListTables(connection)
  if (any(t[1] == table | table == 1)){
    blob <- dbReadTable(connection, tables[1], row.names=NULL)
    if (length(tag) > 0)
      blob$provenance <- tag
    retval <- c(retval, list(tbl_blob=blob))
  }
  if (any(t[2] == table | table == 2)){
    run <- dbReadTable(connection, tables[2], row.names=NULL)
    if (length(tag) > 0)
      run$provenance <- tag
    retval <- c(retval, list(tbl_run=run))
  }
  if (any(t[3] == table | table == 3)){
    step <- dbReadTable(connection, tables[3], row.names=NULL)
    if (length(tag) > 0)
      step$provenance <- tag
    retval <- c(retval, list(tbl_step=step))
  }
  
  dbDisconnect(connection)
  
  return(retval)
}



#' Priority queue
#' @param p order of column, i.e., output of order (x)
#' @return Lexical closure which contains the following O(1) methods:
#'  $first(), $second(), $last(), $last2(), and $remove(i) as documented
#'  in Korpela et al. (2014), appendix A.
#' @export
PQ <- function(p) {
  n <- length(p)
  idx <- 1+order(p)  # index at which item i can be found
  p <- c(NA,p,NA)      # head is at p[1] and tail at p[n+2]
  nxt <- c(2:(n+2),NA) # pointer to next
  prv <- c(NA,1:(n+1)) # pointer to previous
  first <- function() p[nxt[1]]
  second <- function() p[nxt[nxt[1]]]
  last <- function() p[prv[n+2]]
  last2 <- function() p[prv[prv[n+2]]]
  history <- rep(NA,n)
  pos <- maxpos <- 0
  remove <- function(i) {
    pos <<- maxpos <<- pos+1 # update position
    history[pos] <<- i       # add this to history
    j <- idx[i]
    prv[nxt[j]] <<- prv[j] # previous of the next is previous of the current
    nxt[prv[j]] <<- nxt[j] # next of the previous is next of the current
    pos
  }
  ## remove/unremove previously removed items. This is not really needed
  ## but costs nothing to include here at this stage...
  goto <- function(newpos) {
    if(newpos<pos) { # go backward
      for(i in pos:(newpos+1)) {
        j <- idx[history[i]]
        prv[nxt[j]] <<- j # unremove
        nxt[prv[j]] <<- j # unremove
      }
    } else if(newpos>pos) { # go forward
      for(i in pos:(newpos-1)) {
        j <- idx[history[i]]
        prv[nxt[j]] <<- prv[j] # remove
        nxt[prv[j]] <<- nxt[j] # remove
      }
    }
    pos <<- newpos
    pos
  }
  show <- function() list(idx=idx,p=p,nxt=nxt,prv=prv,pos=pos,
                          history=if(maxpos==0) c() else history[1:maxpos]) # for debugging
  list(first=first,second=second,last=last,last2=last2,remove=remove,goto=goto,show=show)
}

#' Finds MWE efficiently using separate validation data as documented in Korpela et al.
#' (2014). The algorithm has time complexity of O(m*n*log(n)), where m is the length of
#' time series and n is the rows in the training data or validation data, whichever is
#' larger.
#' @param max_k integer, maximum value of k searched
#' @param max_alpha real number, fraction of validation data samples out of MWE
#' @param data_tr nXm matrix, training data
#' @param data_va n'Xm matrix, validation data
#' @return Returns a list that contains number of time series removed k, final alpha
#'  (should be the largest alpha which is at most max_alpha), and lower and upper
#'  bounds.
#' @export
find_mwe <- function(max_k,max_alpha,data_tr,data_va) {
  m <- dim(data_tr)[2]
  # use priority queues for both training and validation set. This finds it efficient
  # to find the largest and smallest non-removed curves.
  q_tr <- apply(data_tr,2,function(x) PQ(order(x)))
  q_va <- apply(data_va,2,function(x) PQ(order(x)))
  lo0 <- up0 <- NULL
  k <- alpha <- -1
  removed <- 0
  while(k<max_k && removed/dim(data_va)[1]<=max_alpha) {
    idx <- matrix(c(rep(1:m,2),
                    sapply(q_tr,function(q) q$first()), sapply(q_tr,function(q) q$last()),
                    sapply(q_tr,function(q) q$second()),sapply(q_tr,function(q) q$last2())),2*m,3)
    lo <- data_tr[idx[  1:m ,c(2,1)]] # current lower envelope
    up <- data_tr[idx[-(1:m),c(2,1)]] # current upper envelope
    for(i in 1:m) { # column/time index i
      j <- q_va[[i]]$first() # index of the smallest item in validation set
      while(data_va[j,i]<lo[i] && removed/dim(data_va)[1]<=max_alpha) {
        ## remove if below lower limit
        sapply(q_va,function(q) q$remove(j))
        j <- q_va[[i]]$first()
        removed <- removed+1
      }
      j <- q_va[[i]]$last() # index of the largest item in validation set
      while(up[i]<data_va[j,i] && removed/dim(data_va)[1]<=max_alpha) {
        ## remove if above upper limit
        sapply(q_va,function(q) q$remove(j))
        j <- q_va[[i]]$last()
        removed <- removed+1
      }
    }
    if(removed/dim(data_va)[1]<=max_alpha) {
      lo0 <- lo
      up0 <- up
      k <- k+1
      alpha <- removed/dim(data_va)[1]
      
      ## greedy MWE algorithm: remove time series which decreases the MWE most:
      a <- aggregate(x=data.frame(gain=abs(data_tr[idx[,c(2,1)]]-data_tr[idx[,c(3,1)]])),
                     by=list(i=idx[,2]),
                     FUN=sum)
      j <- a[which.max(a[,"gain"]),"i"]
      sapply(q_tr,function(q) q$remove(j))
    }
  }
  list(k=k,alpha=alpha,lo=lo0,up=up0)
}


findcurves <- function(data,mtr=5000,mva=5000,alpha=0.05) {
  n <- dim(data)[1]
  samples_tr <- t(replicate(mtr,colMeans(data[sample.int(n,replace=TRUE),]))) # training set
  samples_va <- t(replicate(mva,colMeans(data[sample.int(n,replace=TRUE),]))) # validation set
  a <- find_mwe(mtr-2,alpha,samples_tr,samples_va)
  q <- apply(samples_tr,2,function(x) quantile(x,probs=c(alpha/2,1-alpha/2)))
  data.frame(mean0=colMeans(data),lo=a$lo,up=a$up,k=a$k,alpha=a$alpha,lo0=q[1,],up0=q[2,])
}

plotlines <- function(x, tvec, ...) {
  lines(tvec,x[,"lo0"],lty="dotted",...)
  lines(tvec,x[,"up0"],lty="dotted",...)
  lines(tvec,x[,"lo"],lty="dashed",lwd=1.5,...)
  lines(tvec,x[,"up"],lty="dashed",lwd=1.5,...)
  lines(tvec,x[,"mean0"],lty="solid",lwd=1.5,...)
}

plotlines_comp <- function(x,tvec,...) {
  lines(tvec,x[,"lo0"],lty="dotted",...)
  lines(tvec,x[,"up0"],lty="dotted",...)
  lines(tvec,x[,"mean0"],lty="solid",lwd=1.5,...)
}