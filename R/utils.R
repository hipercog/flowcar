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


### ADD YOUR flowcar-LOCAL USEFUL UTILITY FUNCTIONS HERE! ----
