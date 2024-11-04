colorlist <- function(str){
  if(str == "NULL"){
    NULL
  }else{
    mynames <- lapply(strsplit(str, ";")[[1]], function(x){
      strsplit(x, "=")[[1]][1]
    })
    res <- lapply(strsplit(str,";")[[1]], function(x){
      strsplit(strsplit(x, "=")[[1]][2], ",")[[1]]
    })
    names(res) <- mynames
    res
  }
}
tobool <- function(str){ if(str == "0") FALSE else TRUE }
TF <- function(str, default){ if(str %in% c("T", "TRUE", "true", "True")) {TRUE} else if(str %in% c("F", "FALSE", "false", "False")) {FALSE} else{default} }
toc <- function(str){ if(str == "c()") c() else strsplit(str, ",")[[1]] }
numtona <- function(str){ if(str == "NA") NA else as.numeric(str)}
strnull <- function(str){ if(str == "NULL") NULL else str }
numnull <- function(str){ if(str == "NULL") NULL else as.numeric(str) }
breaksnull <- function(str){ if(str == "NULL") NULL else as.numeric(strsplit(str, ",")[[1]]) }
labelsnull <- function(str){ if(str == "NULL") NULL else strsplit(str, ",")[[1]] }

check_region <- function(data){
  absence_sample_n <- length(which(rowSums(data) == 0))
  if(absence_sample_n > 0){
    suffix <- ifelse(absence_sample_n > 3, "...", "")
    warning(paste0("Remove ", absence_sample_n, " regions (", paste0(rownames(data)[which(rowSums(data) == 0)][0:min(3,absence_sample_n)], collapse = ",") , suffix, ") which are 0(absence) in all samples."))
    data <- data[-which(rowSums(data) == 0),]
  }
  return(data)
}
