# @Author: andreas.bender@stat.uni-muenchen.de

get_last <- function(col) {
  res <- rep(FALSE, length(col))
  for (i in 2:length(col)) {
    if (col[i] != col[i - 1]) {
      res[i - 1] <- TRUE
    }
  }
  res[length(col)] <- TRUE
  res
}

get_time_difference <- function(
  prior.var, 
  later.var, 
  data           = patient,
  time.threshold = 48,
  time.unit      = "hours",
  id.var         = "CombinedID") {
  
  first.var  <- data[, prior.var]
  second.var <- data[, later.var]
  
  ## compute time difference in time.unit units
  t.diff <- difftime(first.var, second.var, units = time.unit)
  
  ## time difference bigger than time.threshold time.unit
  t.diff.flag <- t.diff > time.threshold & !is.na(t.diff)
  
  d.th <- cbind(data[t.diff.flag, c(id.var, prior.var, later.var)], 
                tdiff=t.diff[t.diff.flag])
  
  d.th
  
}

##' checks if prior.var occured earlier than later.var (on a timeline)
#' 
#' If threshold is specified throws error if instances occur, where difference
#' (later.var - prior.var) > threshold occur, otherwise threshold=infinity
#' If no error is thrown returns data set for all cases where difference > 0
#' @import checkmate

check_time_consistency <- function(
  data, 
  prior.var, 
  later.var, 
  threshold = Inf,
  id.var    = "CombinedID",
  time.unit = "hours") {
  
  library(checkmate)
  
  assert_character(id.var, len=1, min.chars=1, any.missing=FALSE)
  ## compute time difference in time.unit units
  t.diff <- assert_tdiff(data, prior.var, later.var, threshold=threshold, 
                         time.unit=time.unit)
  
  ## later time before prior time?
  t.diff.flag <- t.diff > 0 & !is.na(t.diff)
  
  d0 <- cbind(
    data[t.diff.flag, c(id.var, prior.var, later.var)], 
    tdiff=t.diff[t.diff.flag])
  
  d0
  
}

get_tdiff <- function(prior, later, units="hours") {
  
  library(checkmate)
  
  assert_class(prior, "POSIXct")
  assert_class(later, "POSIXct")
  
  difftime(later, prior, units=units)
  
}


assert_tdiff <- function(
  data, 
  prior.var, 
  later.var, 
  threshold = 48L,
  time.unit = "hours") {
  
  library(checkmate)
  
  assert_data_frame(data)
  assert_character(prior.var, len=1, min.chars=1, any.missing=FALSE)
  assert_character(later.var, len=1, min.chars=1, any.missing=FALSE)
  assert_number(threshold)
  assert_character(time.unit, len=1, min.chars=1, any.missing=FALSE)
  assert_choice(time.unit, eval(formals(difftime)$units))
  
  t.diff <- get_tdiff(data[[prior.var]], data[[later.var]], units=time.unit)
  
  if(!all(t.diff <= threshold, na.rm=TRUE)) 
    stop(paste("Differences (", later.var, "-", prior.var, ") >", threshold, "occured"))
  
  invisible(t.diff)
  
}


## collect data on current state of patient, daily and ICU data sets
snap_data <- function() {
  
  nrow.pat       <- nrow(patient)
  cid.pat        <- unique(patient$CombinedID)
  cid.length.pat <- length(patient$CombinedID)
  
  nrow.icu     <- nrow(ICU)
  icuid        <- unique(ICU$CombinedicuID)
  icuid.length <- length(icuid)
  
  nrow.daily       <- nrow(daily)
  cid.daily        <- unique(daily$CombinedID)
  cid.length.daily <- length(cid.daily)
  
  list(
    patient = list(nrow = nrow.pat,   cid = cid.pat,   n.pat = cid.length.pat),
    ICU     = list(nrow = nrow.icu,   cid = icuid,     n.icu = icuid.length),
    daily   = list(nrow = nrow.daily, cid = cid.daily, n.pat = cid.length.daily))
  
}

eval_rmlist <- function(rm.list, which="patient", what="n.pat") {
  
  sapply(rm.list, function(z) z[[which]][[what]])
  
}
