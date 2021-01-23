#' Mean with detrending
#'
#' Calculate the mean of values with optional (linear) detrending, and perform check that there's enough data
#'
#' @param X numeric
#' @param detrend TRUE or FALSE, default = FALSE
#' @param prop_req numeric value between 0 and 1, default = 0.99
#'
#' @return numeric
#' @export
#'
#' @examples
#' y = rnorm(10)
#' Mean(y)
#' Mean(y, detrend = TRUE)
#' y[1] = NA
#' Mean(y)
#' Mean(y, prop_req = 0.5)
Mean <- function(X, detrend = FALSE, prop_req = 0.99){
  N = sum(!is.na(X))
  if((N / length(X)) < prop_req){
    print("Warning: more than threshold of rolling window is NA, returning NA")
    Mean = NA
  }else{
    if(detrend == TRUE){
      x_t = 1:length(X)
      lm_x = lm(X~x_t, na.action="na.exclude")
      xd = residuals(lm_x)
    }else{
      xd = X
    }
    Mean = mean(xd, na.rm=TRUE)
  }
  return(Mean)
}


# functions for calculating sd of rolling window stats
Ar1 <-function(x, Detrend=FALSE, propReq = 0.99){
  N = sum(!is.na(x))
  if((N / length(x)) < propReq){
    print("Warning: more than threshold of rolling window is NA, returning NA")
    ac = NA
    ac_sd = NA
  }else{
    if(Detrend == TRUE){
      x_t = 1:length(x)
      lm_x = lm(x~x_t, na.action="na.exclude")
      xd = residuals(lm_x)
    }else{
      xd = x
    }
    ac = acf(xd, lag=1, na.action=na.pass, plot=FALSE)$acf[2]
  }
  return(ac)
}

Ar1_sd <-function(x, Detrend=FALSE, propReq = 0.99){
  N = sum(!is.na(x))
  if((N / length(x)) < propReq){
    # print("Warning: more than threshold of rolling window is NA, returning NA")
    ac = NA
    ac_sd = NA
  }else{
    if(Detrend == TRUE){
      x_t = 1:length(x)
      lm_x = lm(x~x_t, na.action="na.exclude")
      xd = residuals(lm_x)
    }else{
      xd = x
    }
    ac = acf(xd, lag=1, na.action=na.pass, plot=FALSE)$acf[2]
    var_ac = ((1 - ac^2)^2)/(sum(!is.na(x)) - 3)
    ac_sd = sqrt(var_ac)
  }
  return(ac_sd)
}

SD <- function(x, Detrend=FALSE, propReq = 0.99){
  N = sum(!is.na(x))
  if((N / length(x)) < propReq){
    # print("Warning: more than threshold of rolling window is NA, returning NA")
    sd.out = NA
    sd.sd = NA
  }else{
    if(Detrend == TRUE){
      x_t = 1:length(x)
      lm_x = lm(x~x_t, na.action="na.exclude")
      xd = residuals(lm_x)
    }else{
      xd = x
    }
    n.x = sum(!is.na(xd))
    sd.out = sd(xd, na.rm = TRUE)
  }
  return(sd.out)
}

SD_sd <- function(x, Detrend=FALSE, propReq = 0.99){
  N = sum(!is.na(x))
  if((N / length(x)) < propReq){
    # print("Warning: more than threshold of rolling window is NA, returning NA")
    sd.out = NA
    sd.sd = NA
  }else{
    if(Detrend == TRUE){
      x_t = 1:length(x)
      lm_x = lm(x~x_t, na.action="na.exclude")
      xd = residuals(lm_x)
    }else{
      xd = x
    }
    n.x = sum(!is.na(xd))
    sd.out = sd(xd, na.rm = TRUE)
    kurt.x = kurtosis(xd, na.rm=TRUE)
    var.s2 = (sd.out^4)*((2/(n.x-1)) + (kurt.x / n.x))
    var.sd = ( 1/(4*sd.out*sd.out))*var.s2
    sd.sd = sqrt(var.sd)
  }
  return(sd.sd)
}

# function that calculates all R.W. stats
calcRollStats <- function(Data, IDcols=c("Lake", "Year"), Stats=c("Mean", "SD", "SD_sd", "CV", "Ar1", "Ar1_sd"), varCols, varsToLog=NULL, logOffset=.1, Widths=c(21), Detrend=FALSE, minProp=0.99){
  if(!is.null(varsToLog)){
    Data[, paste0(varsToLog, "_log10")] = NA
    for(v in 1:length(varsToLog)){
      minVal = min(Data[, varsToLog[v]], na.rm=TRUE)
      if(minVal < 0){
        addOffset = abs(minVal) + logOffset
      }else{
        addOffset = 0
      }
      print(min(c(0, abs(minVal) + logOffset)))
      Data[, paste0(varsToLog[v], "_log10")] = log10(Data[, varsToLog[v]] + addOffset)
    }
    varCols = c(varCols, paste0(varsToLog, "_log10"))
  }

  lakeYears = unique(Data[, IDcols])
  statsCombinations = expand.grid(Stat=Stats, Width=Widths, Detrend=Detrend, stringsAsFactors = FALSE)
  allLY_allStats = list()
  for(i in 1:nrow(lakeYears)){
    ly_inds = Data[, IDcols[1]] == lakeYears[i, IDcols[1]] & Data[, IDcols[2]] == lakeYears[i, IDcols[2]]
    ly_statsList = list()
    for(s in 1:nrow(statsCombinations)){
      holdResults_lys = as.data.frame(rollapplyr(data=Data[ly_inds, varCols], width=statsCombinations[s, "Width"], FUN=get(statsCombinations[s, "Stat"]), Detrend=statsCombinations[s, "Detrend"], propReq=minProp, fill=NA))
      holdResults_lys$DOYtrunc = rollapplyr(Data[ly_inds, "DOYtrunc"], width=statsCombinations[s, "Width"], FUN=max, partial=TRUE)
      holdResults_lys$Stat = statsCombinations[s, "Stat"]
      holdResults_lys$Width = statsCombinations[s, "Width"]
      ly_statsList[[s]] = holdResults_lys
    }
    ly_allStats = rbind.fill(ly_statsList)
    ly_allStats$Lake = lakeYears[i, "Lake"]
    ly_allStats$Year = lakeYears[i, "Year"]
    allLY_allStats[[i]] = ly_allStats
  }
  out_stats = rbind.fill(allLY_allStats)
  return(out_stats)
}
