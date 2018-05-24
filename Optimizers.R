#*****************************************************************
#Optim Max Sharpe Combn Portfolios
#*****************************************************************

#' @title optim.max.sharpe.combn.portfolios
#' @export
optim.max.sharpe.combn.portfolios<-function(ia,
                                            constraints,
                                            a.combs=2:ia$n,
                                            w.thresh=10^-4){
  #*****************************************************************
  #Define Assets
  #*****************************************************************
  assets=ia$symbols
  w.0=as.numeric(rep(0,ia$n),assets)
  a.combs=a.combs[order(-a.combs)]#order from largest subset to smallest
  #*****************************************************************
  #Define asset combinations
  #*****************************************************************
  mat=permutations(2,a.combs[1],v=c(0,1),repeats.allowed=TRUE)
  idx=which(rowSums(mat)<a.combs[length(a.combs)])
  mat=mat[-idx,]
  colnames(mat)=assets
  asset.mat=mat
  rs.mat = rowSums(asset.mat)
  #*****************************************************************
  #Define w.mat
  #*****************************************************************
  w.mat = asset.mat
  w.mat[]=0
  #*****************************************************************
  #while
  #*****************************************************************
  nr=nrow(asset.mat)
  
  idx.in.store=rep(F,nr)
  iter=0
  while(any(!idx.in.store)){
    iter=iter+1
    #*****************************************************************
    #index largest remaining subsest
    #*****************************************************************
    idx=which.max(rs.mat)
    a.mat.idx=asset.mat[idx,]
    a = names(a.mat.idx[a.mat.idx>0])
    
    
    if(max(ia$expected.return[a])<=0){
      idx=which(rs.mat==(max(rs.mat,na.rm=T)))
      a.mat.idx=asset.mat[idx,]
      rets = ia$expected.return*t(a.mat.idx)
      rets[rets<0]=0
      rets=colSums(rets)
      idx.2=which.max(rets)
      if(rets[idx.2]<=0)
        break
      idx=idx[idx.2]
      a.mat.idx=asset.mat[idx,]
      a = names(a.mat.idx[a.mat.idx>0])
    }
    
    #*****************************************************************
    #subset ia and constraints
    #*****************************************************************
    ia.tmp = create.sub.ia(ia,a)
    constraints.tmp = clean.constraint(ia.tmp,constraints)
    #*****************************************************************
    #optimize
    #*****************************************************************
    w = optim.max.sharpe.portfolio(ia.tmp,constraints.tmp)
    w.in = w[w>w.thresh]
    #*****************************************************************
    #subset index
    #*****************************************************************
    a.in = names(w.in)
    
    idx.in = rowSums(asset.mat[,a.in,drop=F])==length(w.in)
    
    if(iter>1) idx.in[idx.in.store&idx.in]=F#avoid repeats
    #*****************************************************************
    #fill in weights
    #*****************************************************************
    w.mat[idx.in,a.in]=matrix(w.in,length(which(idx.in)),length(w.in),byrow=T)
    rs.mat[idx.in]=NA#=index[!idx.in]
    idx.in.store[idx.in]=T
  }
  #out=list(Asset.Matrix=asset.mat,Weight=w.mat)
  return(w.mat)
}


#ia = create.ia(mu[311,],sigma[,,311])
#constraints = clean.constraint(ia,list())
#optim.max.sharpe.combn.portfolios(ia,constraints)

#' @title sharpeWeights
#' @export
sharpeWeights <- function(ia,constraints){
  tot.assets = ia$n
  tickers = ia$symbols
  excluded = is.na(ia$expected.return)
  if (any(!excluded)){
    ia$cov = ia$cov[!excluded,!excluded]
    ia$expected.return = ia$expected.return[!excluded]
    ia$n = len(ia$expected.return)
    ia$symbols = ia$symbols[!excluded]
    
    if (ia$n == tot.assets){weight = optim.max.sharpe.portfolio(ia,constraints)
    } else if (ia$n==1){weight = as.integer(as.logical(!excluded))
    } else{weight = matrix(0,nrow=1,ncol=tot.assets)
    weight[-which(excluded)] = optim.max.sharpe.portfolio(ia,constraints)
    }
    
  } else {
    weight = rep(0,tot.assets)
  }
  names(weight)=(tickers)
  return(weight)
}


#' @title create.ia
#' @export
create.ia = function(eR,covar){
  ia = list()
  ia$symbols = names(eR)
  ia$n = len(eR)
  ia$expected.return = (eR)
  ia$cov = (covar)
  return(ia)
}

#*****************************************************************
#sub ia
#*****************************************************************
#' @title create.sub.ia
#' @export
create.sub.ia<-function(ia, assets){
  ia$symbols=assets
  ia$n = length(assets)
  
  ia$risk = ia$risk[assets]
  ia$correlation =ia$correlation[assets,assets]
  ia$cov = ia$cov[assets,assets]
  ia$expected.return = ia$expected.return[assets]
  return(ia)
}


#*****************************************************************
#Clean Constraint
#*****************************************************************
#' @title clean.constraint
#' @export
clean.constraint<-function(ia,constraints){
  #assets
  a=ia$symbols
  
  #index constraints
  constraints$n=length(a)
  
  #lower and upper bounds
  constraints$lb = constraints$lb[a]
  constraints$ub = constraints$ub[a]
  
  #constraint matrix
  a.meq=c(colnames(constraints$A)[constraints$meq],a)
  a.meq=which(colnames(constraints$A)%in%a.meq)
  constraints$A=constraints$A[a,a.meq]
  constraints$b=constraints$b[a.meq]
  
  return(constraints)
  
  
}
#' @title optim.max.sharpe.portfolio
#' @export
optim.max.sharpe.portfolio = function(ia.tmp,constraints.tmp){
  port = max.sharpe.portfolio()(ia.tmp,constraints.tmp)
  names(port)= ia.tmp$symbols
  return(port)
}

#' @title clean.ia
#' @export
clean.ia = function(eR,covar){
  excluded = is.na(eR)
  ia = list()
  ia$expected.return = matrix(eR)[!excluded]
  ia$cov = covar[!excluded,!excluded]
  ia$symbols = names(eR)[!excluded]
  ia$n = len(ia$expected.return)
  names(ia$expected.return)=ia$symbols
  return(list(ia,excluded))
}


