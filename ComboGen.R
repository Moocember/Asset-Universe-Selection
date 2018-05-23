#' @title uniWeightGen
#' @export
uniWeightGen <- function(mu,sigma,lookBack,minAssets,maxAssets){
  library(parallel)
  cl = makeCluster(detectCores())
  clusterExport(cl, "mu")
  clusterExport(cl, "sigma")
  clusterExport(cl, "lookBack")
  
  clusterEvalQ(cl,{
    library(gtools)
    library(SIT)
    library(quadprog)
    library(FastUniSelect)
    ia = create.ia(mu[lookBack+2,],sigma[,,lookBack+2])
    constraints = clean.constraint(ia,list())
  })
  
  weights = parLapply(cl,seq(lookBack+2,nrow(roc)),function(x){
    ia$expected.return = (mu[x,])
    ia$cov = (sigma[,,x])
    #pre
    optim.max.sharpe.combn.portfolios(ia,constraints)
    #post
  })
  stopCluster(cl)
  return(weights)
}

#' @title output2
#' @export
output2 = function(roc,weights){
  roc[is.na(roc)] = 0
  stratRoc = sapply(seq(2,(len(weights)-1)),function(x){
    weight = weights[[x-1]]
    weight = weight/rowSums(abs(weight),na.rm=TRUE)
    weight[is.na(weight)] <- 0
    weight%*%t(roc[x+1+nrow(roc)-len(weights)])#index addition and subtraction compensates for mlags
  })
  return(stratRoc)
}
#' @title performance
#' @description Performance of all combinations of given assets
#' @export
#' @keywords
#' @seealso \code{\link[utils]{head}}
#' @return NULL
#' @examples \dontrun{
#' }
performance = function(combo,backTestLookBack,performanceMetric,is.expanding){
  if (is.expanding) {
    if(performanceMetric == "CAGR"){
      performance = RollingMean(combo,backTestLookBack,expanding = TRUE,na_method = "ignore")
    }
    if(performanceMetric == "sharpeRatio"){
      performance = RollingSharpe(combo,rep(0,nrow(combo)),backTestLookBack,expanding = TRUE,na_method = "ignore")
    }
  } else{
    if(performanceMetric == "CAGR"){
      performance = RollingMean(combo,backTestLookBack,na_method = "ignore")
    }
    if(performanceMetric == "sharpeRatio"){
      performance = RollingSharpe(combo,rep(0,nrow(combo)),backTestLookBack,na_method = "ignore")
    }
  }
  return(performance)
}

#' @title momentumUniverseSelection
#' @description Uses universe selection techniques to choose universes
#' @export
#' @keywords
#' @seealso \code{\link[utils]{head}}
#' @return NULL
#' @examples \dontrun{
#' }
momentumUniverseSelection <- function(combo,performance,backTestLookBack,lookBack,rebalanceDays,n){
  returns = sapply(seq((backTestLookBack+lookBack+3),nrow(combo),(rebalanceDays+1)),function(x){
    rowMeans(combo[(x-rebalanceDays):x,performance[(x-rebalanceDays),] > quantile(performance[(x-rebalanceDays),],prob = 1-n ,na.rm = TRUE)],na.rm = TRUE)
  })
  return(as.vector(returns))
}




#' @title momentumBaseUniverse
#' @description Generates benchmark universe
#' @export
#' @keywords
#' @seealso \code{\link[utils]{head}}
#' @return NULL
#' @examples \dontrun{
#' }
momentumBaseUniverse <- function(combo,backTestLookBack,lookBack,rebalanceDays){
  returns = sapply(seq((backTestLookBack+lookBack+3),nrow(combo),(rebalanceDays+1)),function(x){
    combo[(x-rebalanceDays):x,ncol(combo)]
  })
  return(as.vector(returns))
}
