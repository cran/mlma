# data clean #start line #292
data.org<-function(x, levelx=1, levely=1, m, l1=NULL,l2=NULL, c1=NULL, c1r=rep(1,length(c1)), #levelx is the level of x
                   c2=NULL, c2r=rep(1,length(c2)), f01y=NULL, f10y=NULL,                  #level is the level of observations
                   f02ky=NULL, f20ky=NULL, f01km1=NULL, f01km2=NULL, f10km=NULL,          #weight is the level 1 weight of cases
                   level=1:length(x), weight=rep(1,length(x)))                                        #weight2 is the level 2 weight of cases, weight2=rep(1,length(unique(level[!is.na(level)])))   
{ns.dev<-function (x, df = NULL, knots = NULL, qnots=NULL,intercept = FALSE, Boundary.knots = range(x),derivs1=0) 
{
  nx <- names(x)
  x <- as.vector(x)
  nax <- is.na(x)
  if (nas <- any(nax)) 
    x <- x[!nax]
  if (!missing(Boundary.knots)) {
    Boundary.knots <- sort(Boundary.knots)
    outside <- (ol <- x < Boundary.knots[1L]) | (or <- x > 
                                                   Boundary.knots[2L])
  }
  else outside <- FALSE
  if (!is.null(df) && is.null(knots)) {
    nIknots <- df - 1L - intercept
    if (nIknots < 0L) {
      nIknots <- 0L
      warning(gettextf("'df' was too small; have used %d", 
                       1L + intercept), domain = NA)
    }
    knots <- if (nIknots > 0L) {
      knots <- seq.int(0, 1, length.out = nIknots + 2L)[-c(1L, 
                                                           nIknots + 2L)]
      stats::quantile(x[!outside], knots)
    }
  }
  else {if(is.null(df) && is.null(knots) && !is.null(qnots))
    knots<-quantile(x[!outside], qnots)
  nIknots <- length(knots)}
  Aknots <- sort(c(rep(Boundary.knots, 4L), knots))
  if (any(outside)) {
    basis <- array(0, c(length(x), nIknots + 4L))
    if (any(ol)) {
      k.pivot <- Boundary.knots[1L]
      xl <- cbind(1, x[ol] - k.pivot)
      tt <- splineDesign(Aknots, rep(k.pivot, 2L), 4, c(0, 
                                                        1),derivs=rep(derivs1,2L))
      basis[ol, ] <- xl %*% tt
    }
    if (any(or)) {
      k.pivot <- Boundary.knots[2L]
      xr <- cbind(1, x[or] - k.pivot)
      tt <- splineDesign(Aknots, rep(k.pivot, 2L), 4, c(0,1),derivs=rep(derivs1,2L))
      basis[or, ] <- xr %*% tt
    }
    if (any(inside <- !outside)) 
      basis[inside, ] <- splineDesign(Aknots, x[inside], 
                                      4,derivs=rep(derivs1,length(x[inside])))
  }
  else basis <- splineDesign(Aknots, x, 4,derivs=rep(derivs1,length(x)))
  const <- splineDesign(Aknots, Boundary.knots, 4, c(2, 2),derivs=rep(derivs1,length(Boundary.knots)))
  if (!intercept) {
    const <- const[, -1, drop = FALSE]
    basis <- basis[, -1, drop = FALSE]
  }
  qr.const <- qr(t(const))
  basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1L:2L), 
                                                     drop = FALSE])
  n.col <- ncol(basis)
  if (nas) {
    nmat <- matrix(NA, length(nax), n.col)
    nmat[!nax, ] <- basis
    basis <- nmat
  }
  dimnames(basis) <- list(nx, 1L:n.col)
  a <- list(degree = 3L, knots = if (is.null(knots)) numeric() else knots, 
            Boundary.knots = Boundary.knots, intercept = intercept)
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("ns", "basis", "matrix")
  basis
}



bs.dev<-function (x, df = NULL, knots = NULL, degree = 3, intercept = FALSE, 
                  Boundary.knots = range(x),derivs1=0) 
{
  nx <- names(x)
  x <- as.vector(x)
  nax <- is.na(x)
  if (nas <- any(nax)) 
    x <- x[!nax]
  if (!missing(Boundary.knots)) {
    Boundary.knots <- sort(Boundary.knots)
    outside <- (ol <- x < Boundary.knots[1L]) | (or <- x > 
                                                   Boundary.knots[2L])
  }
  else outside <- FALSE
  ord <- 1L + (degree <- as.integer(degree))
  if (ord <= 1) 
    stop("'degree' must be integer >= 1")
  if (!is.null(df) && is.null(knots)) {
    nIknots <- df - ord + (1L - intercept)
    if (nIknots < 0L) {
      nIknots <- 0L
      warning(gettextf("'df' was too small; have used %d", 
                       ord - (1L - intercept)), domain = NA)
    }
    knots <- if (nIknots > 0L) {
      knots <- seq.int(from = 0, to = 1, length.out = nIknots + 
                         2L)[-c(1L, nIknots + 2L)]
      stats::quantile(x[!outside], knots)
    }
  }
  Aknots <- sort(c(rep(Boundary.knots, ord), knots))
  if (any(outside)) {
    warning("some 'x' values beyond boundary knots may cause ill-conditioned bases")
    derivs <- 0:degree
    scalef <- gamma(1L:ord)
    basis <- array(0, c(length(x), length(Aknots) - degree - 
                          1L))
    if (any(ol)) {
      k.pivot <- Boundary.knots[1L]
      xl <- cbind(1, outer(x[ol] - k.pivot, 1L:degree, 
                           "^"))
      tt <- splineDesign(Aknots, rep(k.pivot, ord), ord, 
                         derivs+derivs1)
      basis[ol, ] <- xl %*% (tt/scalef)
    }
    if (any(or)) {
      k.pivot <- Boundary.knots[2L]
      xr <- cbind(1, outer(x[or] - k.pivot, 1L:degree, 
                           "^"))
      tt <- splineDesign(Aknots, rep(k.pivot, ord), ord, 
                         derivs+derivs1)
      basis[or, ] <- xr %*% (tt/scalef)
    }
    if (any(inside <- !outside)) 
      basis[inside, ] <- splineDesign(Aknots, x[inside], 
                                      ord,derivs=rep(derivs1,length(x[inside])))
  }
  else basis <- splineDesign(Aknots, x, ord, derivs=rep(derivs1,length(x)))
  if (!intercept) 
    basis <- basis[, -1L, drop = FALSE]
  n.col <- ncol(basis)
  if (nas) {
    nmat <- matrix(NA, length(nax), n.col)
    nmat[!nax, ] <- basis
    basis <- nmat
  }
  dimnames(basis) <- list(nx, 1L:n.col)
  a <- list(degree = degree, knots = if (is.null(knots)) numeric(0L) else knots, 
            Boundary.knots = Boundary.knots, intercept = intercept)
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("bs", "basis", "matrix")
  basis
}

  x2fx<-function(x,func) #x is the list of original numerical vector, func is a vector of character functions. 
{ # eg.  func <- c("x","x+1","x+2","x+3","log(x)")
  func.list <- list()
  #test.data <- matrix(data=rep(x,length(func)),length(x),length(func))
  #test.data <- data.frame(test.data)
  result<-NULL
  for(i in 1:length(func)){
    func.list[[i]] <- function(x){}
    body(func.list[[i]]) <- parse(text=func[i])
  }
  #result <- mapply(do.call,func.list,lapply(test.data,list))
  col_fun<-NULL
  z<-1
  for (i in 1:length(func.list))
  {res<-as.matrix(func.list[[i]](x))
   result<-cbind(result,res)
   col_fun<-cbind(col_fun,c(z,z+ncol(res)-1))
   z<-z+ncol(res)}
  list(values=as.matrix(result),col_fun=as.matrix(col_fun))
}

one2two<-function(x,l1,weight=rep(1,length(x))) #l1 is a vector that distribute x to different level 1 groups
{x1<-rep(NA,length(x))
 l1.1<-l1[!is.na(l1)]
 x.1<-x[!is.na(l1)]
 weight.1<-weight[!is.na(l1)]
 weight.1<-ifelse(is.na(weight.1),0,weight.1)   #missing weight will not be counted when calculate the level 1 variable
 x1.1<-rep(0,length(x.1))
 for (i in unique(l1.1))
   x1.1[l1.1==i]<-weighted.mean(x.1[l1.1==i],weight.1[l1.1==i],na.rm=TRUE)
 x1[!is.na(l1)]<-x1.1
 cbind(x1,x-x1)
}

two <- function(x, level, weight = rep(1, nrow(as.matrix(x)))) {
  x <- as.matrix(x)
  levels <- unique(level[!is.na(level)])
  x2 <- matrix(NA, length(levels), dim(x)[2])
  for (i in 1:length(levels)) {
    cho <- (level == levels[i])
    if (sum(cho) > 0) {
      if(sum(cho)==1)
        x2[i, ] <- x[level == levels[i], ]
      else {
        temp <- as.matrix(x[level == levels[i], ])
        weight1 <- weight[level == levels[i]]
        x2[i, ] <- apply(temp, 2, weighted.mean, weight1, 
                         na.rm = TRUE)}
    }
  }
  colnames(x2) <- colnames(x)
  x2
}

x2fdx<-function(x,func)  #x is the list of original numerical vector, func is a vector of character functions. 
{ fdx<-NULL              # eg.  func <- c("x","x+1","x+2","x+3","log(x)")
for(i in 1:length(func)){
  if (length(grep("ifelse",func[i]))>0)
  {str<-unlist(strsplit(func[i],","))
  fun1<-D(parse(text=str[2]), "x")
  fun2<-D(parse(text=unlist(strsplit(str[3],")"))),"x")
  x1<-eval(fun1)
  x2<-eval(fun2)
  if(length(x1)==1)
    x1<-rep(x1,length(x))
  if(length(x2)==1)
    x2<-rep(x2,length(x))
  fun3<-paste(str[1],"x1,x2)",sep=",")
  fdx<-cbind(fdx,eval(parse(text=fun3)))
  }
  else if(length(grep("ns",func[i]))>0)
  {temp<-paste("ns.dev",substring(func[i],3,nchar(func[i])-1),",derivs1=1)",sep="")
   fdx<-cbind(fdx,eval(parse(text=temp)))}
  else if(length(grep("bs",func[i]))>0)
  {temp<-paste("bs.dev",substring(func[i],3,nchar(func[i])-1),",derivs1=1)",sep="")
   fdx<-cbind(fdx,eval(parse(text=temp)))}
  else{
    dx2x <- D(parse(text=func[i]), "x") 
    temp<-eval(dx2x)
    if(length(temp)==1)
      fdx<-cbind(fdx,rep(temp,length(x)))
    else fdx<-cbind(fdx,temp)}
}
as.matrix(fdx)
}

order_char<-function(char1,char2)  #find the position of char2 in char1
{a<-1:length(char1)
 b<-NULL
 for (i in 1:length(char2))
   b<-c(b,a[char1==char2[i]])
 b
}

cattobin<-function(m1y,m1y.der,m1,m,cat1,cat2=rep(1,length(cat1)), level=rep(1,dim(m)[1]),weight=rep(1,dim(m)[1]))
{ ad1<-function(vec)
{vec1<-vec[-1]
 vec1[vec[1]]<-1
 vec1
}

m<-as.matrix(m)
if(is.null(m1y))
{m1<-list(cat1)
 dim1<-c(dim(m)[1],0)}
else{
  dim1<-dim(m1y)
  m1[[1]]<-c(m1[[1]],cat1)}

m12y<-NULL
m12y.der<-NULL
m12<-list(cat1)
g<-dim1[2]
ntemp<-colnames(m)[cat1]
j<-1
p<-0
for (i in cat1)
{a<-m[,i]
 d<-rep(0,dim1[1])
 b<-sort(unique(a[a!=cat2[j]]))
 l<-1
 for (k in b)
 {d[a==k]<-l
  l<-l+1}
 d[a==cat2[j]]<-l
 f<-matrix(0,dim1[1],l-1) 
 colnames(f)<-paste(ntemp[j],b,sep=".")  ##
 hi<-d[d!=l & !is.na(d)]
 f[d!=l & !is.na(d),]<-t(apply(cbind(hi,f[d!=l & !is.na(d),]),1,ad1))
 f[is.na(d),]<-NA
 z<-apply(f,2,one2two,level,weight)
 m1y<-cbind(m1y,f)
 m1y.der<-cbind(m1y.der,f)
 m1<-append(m1,list((g+1):(g+l-1)))
 temp3<-as.matrix(z[1:dim1[1],])
 colnames(temp3)<-paste(ntemp[j],"12",b,sep=".")  ##
 m12y<-cbind(m12y,temp3)
 m12y.der<-cbind(m12y.der,matrix(1,nrow(temp3),ncol(temp3))) ## for changes per unit percent
 m12<-append(m12,list((p+1):(p+l-1)))
 g<-g+length(b)
 p<-p+length(b)
 j<-j+1
}
colnames(m12y.der)<-colnames(m12y)
list(m1y=m1y,m1y.der=m1y.der,m1=m1,m12y=m12y,m12y.der=m12y.der,m12=m12)
}

n<-length(x)
mnames<-colnames(as.matrix(m))
if(is.character(l1))
  l1<-unlist(sapply(l1,grep,mnames))
if(is.character(c1))
  c1<-unlist(sapply(c1,grep,mnames))
if(is.character(l2))
  l2<-unlist(sapply(l2,grep,mnames))
if(is.character(c2))
  c2<-unlist(sapply(c2,grep,mnames))

if(is.character(f02ky[[1]]))
  f02ky[[1]]<-unlist(sapply(f02ky[[1]],grep,mnames))
if(is.character(f20ky[[1]]))
  f20ky[[1]]<-unlist(sapply(f20ky[[1]],grep,mnames))
if(is.character(f01km1[[1]]))
  f01km1[[1]]<-unlist(sapply(f01km1[[1]],grep,mnames))
if(is.character(f01km2[[1]]))
  f01km2[[1]]<-unlist(sapply(f01km2[[1]],grep,mnames))
if(is.character(f10km[[1]]))
  f10km[[1]]<-unlist(sapply(f10km[[1]],grep,mnames))

if(levelx==2)
{l1x<-NULL               #create level 1 and level 2 x variables that are used to predict y
 if(is.null(f01y))
 {x1<-x
  l2x<-1
  x1.der<-rep(1,n)}
 else
 {x1<-x2fx(x,f01y)$values
  x1.der<-x2fdx(x,f01y)
  x1<-as.matrix(x1)
  l2x<-1:ncol(x1)
  colnames(x1)<-paste("x",1:ncol(x1),sep=".")
  x1.der<-as.matrix(x1.der)
  colnames(x1.der)<-paste("x",1:ncol(x1),sep=".")}

 if(!is.null(l2))        #create level 2 m variables to explain y
 {m2y<-as.matrix(m[,l2])
  colnames(m2y)<-colnames(m)[l2]
  m2y.der<-matrix(1,n,length(l2))
  colnames(m2y.der)<-colnames(m)[l2]
  m2<-as.list(c(1,1:length(l2)))
  m2[[1]]<-l2
  if(!is.null(f02ky)) 
    for (i in 2:length(f02ky))
    {a<-f02ky[[1]][i-1]
     if(sum(l2==a)>0)
     {b<-(1:length(l2))[l2==a]
      d<-as.matrix(x2fx(m[,a],f02ky[[i]])$values)
      d.der<-as.matrix(x2fdx(m[,a],f02ky[[i]]))
      colnames(d)<-paste(mnames[a],1:ncol(d),sep=".")
      colnames(d.der)<-paste(mnames[a],1:ncol(d.der),sep=".")
      if(dim(as.matrix(d))[2]==1)
      {m2y[,b]<-d
       m2y.der[,b]<-d.der}
      else {
       m2y[,b]<-d[,1]
       m2y.der[,b]<-d.der[,1]
       m2[[b+1]]<-c(m2[[b+1]],(dim(as.matrix(m2y))[2]+1):(dim(as.matrix(m2y))[2]+ncol(d)-1))
       temp.m2<-c(colnames(m2y),colnames(d)[-1])
       m2y<-cbind(m2y,d[,-1])
       colnames(m2y)<-temp.m2
       m2y.der<-cbind(m2y.der,d.der[,-1])
       colnames(m2y.der)<-temp.m2
      }}
    }
  m2yd<-dim(as.matrix(m2y))[2]}
 else
 {m2y<-NULL
  m2y.der<-NULL
  m2yd<-0
  m2<-NULL
 }
 
 if(!is.null(c2))        #binarize level 2 categorical variables to exaplain y
 {temp<-cattobin(m2y,m2y.der,m2,m,c2,c2r)
  m2y<-temp$m1y
  m2y.der<-temp$m1y.der
  m2<-temp$m1
  m2yd<-dim(as.matrix(m2y))[2]}
 
 if(!is.null(l1))         #create level 1 and corresponding level 2 m variables to explain y
 {temp<-apply(as.matrix(m[,l1]), 2, one2two, level, weight)
  m1y<-as.matrix(temp[(n+1):(2*n),])
  m1y.der<-matrix(1,n,ncol(m1y))
  colnames(m1y)<-mnames[l1]
  colnames(m1y.der)<-mnames[l1]
  m12y<-as.matrix(temp[1:n,])
  colnames(m12y)<-paste(mnames[l1],"12",sep=".")
  m12y.der<-matrix(1,n,ncol(temp))
  colnames(m12y.der)<-paste(mnames[l1],"12",sep=".")
  m1<-as.list(c(1,1:length(l1)))
  m12<-as.list(c(1,(m2yd+1):(m2yd+length(l1))))
  m1[[1]]<-l1
  m12[[1]]<-l1
  #
  if(!is.null(f02ky)) 
    for (i in 2:length(f02ky))
    {a<-f02ky[[1]][i-1]
    if(sum(l1==a)>0)
    {b<-(1:length(l1))[l1==a]
     tt<-one2two(m[,a],level,weight)
     d<-as.matrix(x2fx(tt[,1],f02ky[[i]])$values)
     d.der<-as.matrix(x2fdx(tt[,1],f02ky[[i]]))
     colnames(d)<-paste(mnames[a],1:ncol(d),sep=".")
     colnames(d.der)<-paste(mnames[a],1:ncol(d),sep=".")
     if(dim(as.matrix(d))[2]==1)
     {m12y[,b]<-d
      m12y.der[,b]<-d.der}
    else {
      m12y[,b]<-d[,1]
      m12y.der[,b]<-d.der[,1]
      m12[[b+1]]<-c(m12[[b+1]],(dim(as.matrix(m12y))[2]+1):(dim(as.matrix(m12y))[2]+ncol(d)-1))
      temp.m12<-c(colnames(m12y),colnames(d)[-1])
      m12y<-cbind(m12y,d[,-1])
      colnames(m12y)<-temp.m12
      m12y.der<-cbind(m12y.der,d.der[,-1])
      colnames(m12y.der)<-temp.m12
    }}
    } 
  # added in case the transformation f02ky is on the aggregated level 1 mediator
  if(!is.null(f20ky)) 
    for (i in 2:length(f20ky))
    {a<-f20ky[[1]][i-1]
     b<-(1:length(l1))[l1==a]
     tt<-one2two(m[,a],level,weight)
     d<-as.matrix(x2fx(tt[,2],f20ky[[i]])$values)
     d.der<-as.matrix(x2fdx(tt[,2],f20ky[[i]]))
     colnames(d)<-paste(mnames[a],1:ncol(d),sep=".")
     colnames(d.der)<-paste(mnames[a],1:ncol(d),sep=".")
     if(ncol(d)==1)
     {m1y[,b]<-d
      m1y.der[,b]<-d.der
       }
     else {
       temp.namem1y<-c(colnames(m1y),colnames(d)[-1])
       m1y[,b]<-d[,1]
       m1y.der[,b]<-d.der[,1]
       m1[[b+1]]<-c(m1[[b+1]],(dim(as.matrix(m1y))[2]+1):(dim(as.matrix(m1y))[2]+ncol(d)-1))
       m1y<-cbind(m1y,d[,-1])
       m1y.der<-cbind(m1y.der,d.der[,-1])
       colnames(m1y)<-temp.namem1y
       colnames(m1y.der)<-temp.namem1y
     }
    }
  m2y<-cbind(m2y,m12y)
  m2y.der<-cbind(m2y.der,m12y.der)
 }
 else 
 {m1y<-NULL
  m1y.der<-NULL
  m1<-NULL
  m12<-list(NULL)}
 
 if(!is.null(c1))                  #binarize level 1 categorical variables and corresponding level 2 variables to explain y 
 {temp<-cattobin(m1y,m1y.der,m1,m,c1,c1r,level,weight)
  m1y<-temp$m1y
  m1y.der<-temp$m1y.der
  m1<-temp$m1
  dim2y<-ifelse(is.null(dim(m2y)),0,dim(m2y)[2])
  m2y<-cbind(m2y,temp$m12y)
  m2y.der<-cbind(m2y.der,temp$m12y.der)
  for (i in 2:length(temp$m12))
    m12<-append(m12,list(temp$m12[[i]]+dim2y))
  m12[[1]]<-c(m12[[1]],temp$m12[[1]])  
 }

 lc2<-c(l2,c2)
 if(!is.null(lc2))                  #create x variables to explain level 2 mediatiors
 {temp<-cbind(x,m[,lc2])
  temp2<-two(temp,level)
  xm2<-as.matrix(temp2[,1])
  colnames(xm2)<-"x.2"
  xm2.der<-as.matrix(rep(1,length(temp2[,1])))
  colnames(xm2.der)<-"x.2"
  m.2<-as.matrix(temp2[,-1])
  colnames(m.2)<-mnames[lc2]
  fm22<-as.list(rep(1,length(lc2)+1))
  fm22[[1]]<-lc2
  if(!is.null(f01km2))
  {allfun<-f01km2[[2]]
   if (length(f01km2)>2)
     for(i in 3:length(f01km2))
       allfun<-c(allfun,f01km2[[i]])
   unifun<-unique(allfun)
   unifun1<-unifun[unifun!="x"]
   unifun2<-c("x",unifun1)
   d_d<-x2fx(xm2,unifun1)
   d.der<-as.matrix(x2fdx(xm2,unifun1))
   d<-as.matrix(d_d$values)
   temp.xm2<-c(colnames(xm2),paste("x.2",1:ncol(d),sep="."))
   xm2.der<-cbind(xm2.der,d.der)
   xm2<-cbind(xm2,d)
   colnames(xm2)<-temp.xm2
   colnames(xm2.der)<-temp.xm2
   col_fun<-cbind(c(1,1),d_d$col_fun+1)
   for(i in 2:length(f01km2))
    {ttemp<-order_char(unifun2,f01km2[[i]])
     ttemp1<-NULL
     for (j in 1:length(ttemp))
       ttemp1<-c(ttemp1,col_fun[1,ttemp[j]]:col_fun[2,ttemp[j]])
     fm22[[(1:length(lc2))[fm22[[1]]==f01km2[[1]][i-1]]+1]]<-ttemp1}
  }}
 else
 {m.2<-NULL
  xm2<-NULL
  xm2.der<-NULL
  fm22<-NULL}
 
 fm11<-NULL                        #fm11 is the list of level 1 x to explain level 1 m
 lc1<-c(l1,c1)
 if(!is.null(lc1))                  #create x variables to explain level 1 mediatiors
 {xm1<-x
  xm1.der<-rep(1,n)
  fm12<-as.list(rep(1,length(lc1)+1))
  fm12[[1]]<-lc1
  if(!is.null(f01km1))
  {allfun<-f01km1[[2]]
   if (length(f01km1)>2)
     for(i in 3:length(f01km1))
       allfun<-c(allfun,f01km1[[i]])
   unifun<-unique(allfun)
   unifun1<-unifun[unifun!="x"]
   unifun2<-c("x",unifun1)
   d_d<-x2fx(x,unifun1)
   xm1<-cbind(x,d_d$values)
   d<-as.matrix(d_d$values)
   xm1.der<-cbind(xm1.der,x2fdx(x,unifun1))
   colnames(xm1)<-paste("x2",1:(1+ncol(d)),sep=".")
   colnames(xm1.der)<-paste("x2",1:(1+ncol(d)),sep=".")
   col_fun<-cbind(c(1,1),d_d$col_fun)
   for(i in 2:length(f01km1))
    {ttemp<-order_char(unifun2,f01km1[[i]])
     ttemp1<-NULL
     for (j in 1:length(ttemp))
       ttemp1<-c(ttemp1,col_fun[1,ttemp[j]]:col_fun[2,ttemp[j]])
     fm12[[(1:length(lc1))[fm12[[1]]==f01km1[[1]][i-1]]+1]]<-ttemp1}
  }}
 else
 {xm1<-NULL
  xm1.der<-NULL
  fm12<-NULL}
}
else                          # to deal with level 1 x
{x2<-one2two(x,level,weight)               #create level 1 and level 2 x variables that are used to predict y
 if(!is.null(f01y))
 {x1<-as.matrix(x2fx(x2[,1],f01y)$values)
  x1.der<-x2fdx(x2[,1],f01y)
  l2x<-1:ncol(x1)
  colnames(x1)<-paste("x.j",1:ncol(x1),sep=".")
  colnames(x1.der)<-paste("x.j",1:ncol(x1),sep=".")}
 else
 {x1<-as.matrix(x2[,1])
  x1.der<-matrix(1,n,1)
  colnames(x1)<-"x.j"
  colnames(x1.der)<-"x.j"
  l2x<-1}
 if(!is.null(f10y))
 {x3<-as.matrix(x2fx(x2[,2],f10y)$values)
  x3.der<-x2fdx(x2[,2],f10y)
  l1x<-(length(l2x)+1):(length(l2x)+ncol(x3))
  temp.name<-c(colnames(x1),paste("xij",1:ncol(x3),sep="."))
  x1<-cbind(x1,x3)
  x1.der<-cbind(x1.der,x3.der)
  colnames(x1)<-temp.name
  colnames(x1.der)<-temp.name}
 else
 {l1x<-length(l2x)+1
  temp.name<-c(colnames(x1),"xij")
  x1<-cbind(x1,x2[,2])
  colnames(x1)<-temp.name
  x1.der<-cbind(x1.der,rep(1,n))
  colnames(x1.der)<-temp.name}
 
#if(!(is.null(l2) & is.null(c2)))  #this part has changed to allow level 2 mediators when predictor is level 1
#  stop("x is a level 1 variable, m cannot be level 2: l2 and c2 should be empty!")
 
if(!is.null(l2))        #create level 2 m variables to explain y
{m2y<-as.matrix(m[,l2])
 colnames(m2y)<-colnames(m)[l2]
 m2y.der<-matrix(1,n,length(l2))
 colnames(m2y.der)<-colnames(m)[l2]
 m2<-as.list(c(1,1:length(l2)))
 m2[[1]]<-l2
 if(!is.null(f02ky)) 
  for (i in 2:length(f02ky))
  {a<-f02ky[[1]][i-1]
   if(sum(l2==a)>0)
   {b<-(1:length(l2))[l2==a]
    d_d<-x2fx(m[,a],f02ky[[i]])
    d<-as.matrix(d_d$values)
    d.der<-as.matrix(x2fdx(m[,a],f02ky[[i]]))
    colnames(d)<-paste(mnames[a],1:ncol(d),sep=".")
    colnames(d.der)<-paste(mnames[a],1:ncol(d),sep=".")
    if(dim(as.matrix(d))[2]==1)
     {m2y[,b]<-d
      m2y.der[,b]<-d.der}
   else {
    m2y[,b]<-d[,1]
    m2y.der[,b]<-d.der[,1]
    m2[[b+1]]<-c(m2[[b+1]],(dim(as.matrix(m2y))[2]+1):(dim(as.matrix(m2y))[2]+ncol(d)-1))
    temp.m2<-c(colnames(m2y),colnames(d)[-1])
    m2y<-cbind(m2y,d[,-1])
    colnames(m2y)<-temp.m2
    m2y.der<-cbind(m2y.der,d.der[,-1])
    colnames(m2y.der)<-temp.m2
  }}
  }
 m2yd<-dim(as.matrix(m2y))[2]}
else
{m2y<-NULL
 m2y.der<-NULL
 m2yd<-0
 m2<-NULL
}

if(!is.null(c2))        #binarize level 2 categorical variables to exaplain y
{temp<-cattobin(m2y,m2y.der,m2,m,c2,c2r)
 m2y<-temp$m1y
 m2y.der<-temp$m1y.der
 m2<-temp$m1
 m2yd<-dim(as.matrix(m2y))[2]}





 if(!is.null(l1))         #create level 1 and corresponding level 2 m variables to explain y
 {temp<-apply(as.matrix(m[,l1]), 2, one2two, level, weight)
  m1y<-as.matrix(temp[(n+1):(2*n),])
  m1y.der<-matrix(1,n,ncol(m1y))
  colnames(m1y)<-mnames[l1]
  colnames(m1y.der)<-mnames[l1]
  m2.names<-colnames(m2y)
  m2y<-cbind(m2y,as.matrix(temp[1:n,]))
  m2y.der<-cbind(m2y.der,matrix(1,n,ncol(temp)))
  colnames(m2y)<-c(m2.names,paste(mnames[l1],"12",sep="."))
  colnames(m2y.der)<-c(m2.names,paste(mnames[l1],"12",sep="."))
  m1<-as.list(c(1,1:length(l1)))
  m12<-as.list(c(1,(m2yd+1):(m2yd+length(l1))))
  m1[[1]]<-l1
  m12[[1]]<-l1
  #
  if(!is.null(f02ky)) 
    for (i in 2:length(f02ky))
    {a<-f02ky[[1]][i-1]
    if(sum(l1==a)>0)
    {b<-(1:length(l1))[l1==a]
     tt<-one2two(m[,a],level,weight)
     d_d<-x2fx(tt[,1],f02ky[[i]])
     d<-as.matrix(d_d$vaules)
     d.der<-as.matrix(x2fdx(tt[,1],f02ky[[i]]))
     colnames(d)<-paste(mnames[a],1:ncol(d),sep=".")
     colnames(d.der)<-paste(mnames[a],1:ncol(d),sep=".")
    if(dim(as.matrix(d))[2]==1)
    {m2y[,m2yd+b]<-d
     m2y.der[,m2yd+b]<-d.der}
    else {
      m2y[,m2yd+b]<-d[,1]
      m2y.der[,m2yd+b]<-d.der[,1]
      m12[[b+1]]<-c(m12[[b+1]],(dim(as.matrix(m2y))[2]+1):(dim(as.matrix(m2y))[2]+ncol(d)-1))
      temp.m2<-c(colnames(m2y),colnames(d)[-1])
      m2y<-cbind(m2y,d[,-1])
      colnames(m2y)<-temp.m2
      m2y.der<-cbind(m2y.der,d.der[,-1])
      colnames(m2y.der)<-temp.m2
    }}
    } 
  # added in case the transformation f02ky is on the aggregated level 1 mediator
  
  if(!is.null(f20ky)) 
    for (i in 2:length(f20ky))
    {a<-f20ky[[1]][i-1]
     b<-(1:length(l1))[l1==a]
     tt<-one2two(m[,a],level,weight)
     d_d<-x2fx(tt[,2],f20ky[[i]])
     d<-as.matrix(d_d$values)
     d.der<-x2fdx(tt[,2],f20ky[[i]])
     colnames(d)<-paste(mnames[a],1:ncol(d),sep=".")
     colnames(d.der)<-paste(mnames[a],1:ncol(d),sep=".")
     if(ncol(d)==1)
     {m1y[,b]<-d
      m1y.der[,b]<-d.der
      }
     else {
       m1y[,b]<-d[,1]
       m1y.der[,b]<-d.der[,1]
       tempname.m1y<-c(colnames(m1y),colnames(d)[-1])
       m1[[b+1]]<-c(m1[[b+1]],(dim(as.matrix(m1y))[2]+1):(dim(as.matrix(m1y))[2]+ncol(d)-1))
       m1y<-cbind(m1y,d[,-1])
       m1y.der<-cbind(m1y.der,d.der[,-1])
       colnames(m1y)<-tempname.m1y
       colnames(m1y.der)<-tempname.m1y
     }
    }
 }
 else 
 {m1y<-NULL
  m1y.der<-NULL
  m1<-NULL
  m12<-NULL
  #m2y<-NULL   
  #m2y.der<-NULL
  }
 
 if(!is.null(c1))                  #binarize level 1 categorical variables and corresponding level 2 variables to explain y 
 {temp<-cattobin(m1y,m1y.der,m1,m,c1,c1r,level,weight)
  m1y<-temp$m1y
  m1y.der<-temp$m1y.der
  m1<-temp$m1
  dim2y<-ifelse(is.null(m2y),0,dim(as.matrix(m2y))[2])
  if(is.null(m2y))
    m12<-temp$m12
  else
  {for (i in 2:length(temp$m12))
    m12<-append(m12,list(temp$m12[[i]]+dim2y))
   m12[[1]]<-c(m12[[1]],temp$m12[[1]])}
  m2y<-cbind(m2y,temp$m12y)
  m2y.der<-cbind(m2y.der,temp$m12y.der)
 }
 

lc2<-c(l2,c2)                      #added: aggregated level 1 x to explain level 2 m
if(!is.null(lc2))                  #create x variables to explain level 2 mediatiors
{temp<-cbind(x,m[,lc2])
 temp2<-two(temp,level)
 xm2<-as.matrix(temp2[,1])
 colnames(xm2)<-"x.2"
 xm2.der<-as.matrix(rep(1,length(temp2[,1])))
 colnames(xm2.der)<-"x.2"
 m.2<-as.matrix(temp2[,-1])
 colnames(m.2)<-mnames[lc2]
 fm22<-as.list(rep(1,length(lc2)+1))
 fm22[[1]]<-lc2
 if(!is.null(f01km2))
 {allfun<-f01km2[[2]]
  if (length(f01km2)>2)
   for(i in 3:length(f01km2))
    allfun<-c(allfun,f01km2[[i]])
  unifun<-unique(allfun)
  unifun1<-unifun[unifun!="x"]
  unifun2<-c("x",unifun1)
  xm2.der<-cbind(xm2.der,x2fdx(xm2,unifun1))
  d_d<-x2fx(xm2,unifun1)
  d<-as.matrix(d_d$values)
  temp.xm2<-c(colnames(xm2),paste("x.2",1:ncol(d),sep="."))
  xm2<-cbind(xm2,d)
  colnames(xm2)<-temp.xm2
  colnames(xm2.der)<-temp.xm2
  col_fun<-cbind(c(1,1),d_d$col_fun+1)
  for(i in 2:length(f01km2))
   {ttemp<-order_char(unifun2,f01km2[[i]])
    ttemp1<-NULL
    for (j in 1:length(ttemp))
      ttemp1<-c(ttemp1,col_fun[1,ttemp[j]]:col_fun[2,ttemp[j]])
    fm22[[(1:length(lc2))[fm22[[1]]==f01km2[[1]][i-1]]+1]]<-ttemp1}
}}
else
{m.2<-NULL
 xm2<-NULL
 xm2.der<-NULL
 fm22<-NULL}




 lc1<-c(l1,c1)
 if(!is.null(lc1))                  
 {fm12<-as.list(rep(1,length(lc1)+1))  #create level 2 x variables to explain level 1 mediatiors
  fm12[[1]]<-lc1
  if(!is.null(f01km1))
  {allfun<-f01km1[[2]]
   if (length(f01km1)>2)
     for(i in 3:length(f01km1))
       allfun<-c(allfun,f01km1[[i]])
   unifun<-unique(allfun)
   unifun1<-unifun[unifun!="x"]
   unifun2<-c("x",unifun1)
   d_d<-x2fx(x2[,1],unifun1)
   d<-as.matrix(d_d$values)
   xm1<-cbind(x2[,1],d)
   xm1.der<-cbind(rep(1,n),x2fdx(x2[,1],unifun1))
   colnames(xm1)<-c("x.2",paste("x.2",1:ncol(d),sep="."))
   colnames(xm1.der)<-c("x.2",paste("x.2",1:ncol(d),sep="."))
   col_fun<-cbind(c(1,1),d_d$col_fun+1)
   for(i in 2:length(f01km1))
    {ttemp<-order_char(unifun2,f01km1[[i]])
     ttemp1<-NULL
     for (j in 1:length(ttemp))
      ttemp1<-c(ttemp1,col_fun[1,ttemp[j]]:col_fun[2,ttemp[j]])
     fm12[[(1:length(lc1))[fm12[[1]]==f01km1[[1]][i-1]]+1]]<-ttemp1}
  }
  else
  {xm1<-as.matrix(x2[,1])
   colnames(xm1)<-"x.2"
   xm1.der<-matrix(1,n,1)
   colnames(xm1.der)<-"x.2"}
  dimxm1<-dim(as.matrix(xm1))
  fm11<-as.list(rep(dimxm1[2]+1,length(lc1)+1))  #create level 1 x variables to explain level 1 mediatiors
  fm11[[1]]<-lc1
  if(!is.null(f10km))
  {allfun<-f10km[[2]]
   if (length(f10km)>2)
     for(i in 3:length(f10km))
       allfun<-c(allfun,f10km[[i]])
   unifun<-unique(allfun)
   unifun1<-unifun[unifun!="x"]
   unifun2<-c("x",unifun1)
   d_d<-x2fx(x2[,2],unifun1)
   temp2<-as.matrix(d_d$values)
   #temp3<-apply(temp2,2,one2two,level,weight)
   temp2.der<-x2fdx(x2[,2],unifun1)
   #temp3.der<-apply(temp2.der,2,one2two,level,weight)
   #for(i in 2:length(f10km))
   #   fm12[[(1:length(lc1))[fm12[[1]]==f10km[[1]][i-1]]+1]]<-c(fm12[[(1:length(lc1))[fm12[[1]]==f10km[[1]][i-1]]+1]],
   #                                                           order_char(unifun2,f10km[[i]])+dimxm1[2])
   #dimxm1<-c(dimxm1[1],dimxm1[2]+dim(temp3)[2])
   col_fun<-cbind(c(1,1),d_d$col_fun+1)
   for(i in 2:length(f10km))
    {ttemp<-order_char(unifun2,f10km[[i]])
     ttemp1<-NULL
     for (j in 1:length(ttemp))
      ttemp1<-c(ttemp1,col_fun[1,ttemp[j]]:col_fun[2,ttemp[j]])
     fm11[[(1:length(lc1))[fm11[[1]]==f10km[[1]][i-1]]+1]]<-ttemp1+dimxm1[2]}
   temp.xm1<-c(colnames(xm1),"x1",paste("x1",1:ncol(temp2),sep="."))
   xm1<-cbind(xm1,x2[,2],temp2)
   xm1.der<-cbind(xm1.der,rep(1,n),temp2.der)
   colnames(xm1)<-temp.xm1
   colnames(xm1.der)<-temp.xm1
  }
  else
  {xm1<-cbind(xm1,x1=x2[,2])
   xm1.der<-cbind(xm1.der,rep(1,n))
   colnames(xm1.der)<-colnames(xm1)}
 }
 else
 {xm1<-NULL
  xm1.der<-NULL
  fm12<-NULL
  fm11<-NULL}
}

if(!is.null(x1))
{x1<-as.matrix(x1)
if(is.null(colnames(x1)))
colnames(x1)<-"x1"
x1.der<-as.matrix(x1.der)
colnames(x1.der)<-colnames(x1)}
if(!is.null(xm1))
{xm1<-as.matrix(xm1)
if(is.null(colnames(xm1)))
colnames(xm1)<-"xm1"
xm1.der<-as.matrix(xm1.der)
colnames(xm1.der)<-colnames(xm1)}

list(x1=x1, x1.der=x1.der, l1x=l1x, l2x=l2x, m1y=m1y, m1y.der=m1y.der,
     m1=m1, m2y=m2y, m2y.der=m2y.der, m2=m2, m12=m12, xm1=xm1, xm1.der=xm1.der,
     fm11=fm11, fm12=fm12, m.2=m.2, xm2=xm2, xm2.der=xm2.der, fm22=fm22,
     parameter=list(levelx=levelx, levely=levely, mnames=mnames,l1 = l1, l2 = l2,  
                    c1 = c1, c1r = c1r, c2 = c2, c2r = c2r, f01y = f01y, level=level,
                    f10y = f10y, f02ky = f02ky, f20ky = f20ky, f01km1 = f01km1, 
                    f01km2 = f01km2, f10km = f10km, level, weight = weight))
}

#multilevel mediation analysis
mlma<-function(y, biny=FALSE, data1=NULL, x, levelx=1, levely=1, m, l1=NULL,l2=NULL, c1=NULL, #levelx is the level of x
               c1r=rep(1,length(c1)), c2=NULL, c2r=rep(1,length(c2)), level=1:length(y),  
               weight=rep(1,length(x)), random="(1|level)", random.m1=NULL,intercept=TRUE, 
               covariates=NULL,cy1=NULL, cy2=NULL, cm=NULL,joint=NULL,org.data=FALSE,
               f01y=NULL, f10y=NULL, f02ky=NULL, f20ky=NULL, f01km1=NULL, f01km2=NULL, f10km=NULL)                               
  
{two<-function(x, level, weight=rep(1,nrow(as.matrix(x))))
{x<-as.matrix(x)
levels<-unique(level[!is.na(level)])
x2<-matrix(NA,length(levels),dim(x)[2])
for(i in 1:length(levels))
{cho<-(level==levels[i])
if(sum(cho)>0)
{if(sum(cho)==1)
  x2[i,]<-x[level==levels[i],]
else
{temp<-as.matrix(x[level==levels[i],])
weight1<-weight[level==levels[i]]
x2[i,]<-apply(temp,2,weighted.mean,weight1,na.rm=TRUE)}}}
colnames(x2)<-colnames(x)
x2
}

getformula<-function(expl,random="(1|level)",intercept=TRUE)
{temp.name<-colnames(expl)
 formula<-temp.name[1]
 if(ncol(expl)>1)
   for (i in 2:ncol(expl))
     formula<-paste(formula,temp.name[i],sep="+")
 formula<-ifelse(intercept,formula,paste(formula,"1",sep="-"))
 formula<-paste("y~",formula)
 if (!is.null(random))
   formula<-paste(formula,random,sep="+")
 formula
}

#n<-length(y)
#tt2<-cbind(data1$x1[,data1$l1x],covariates[,cy2])
#colnames(tt2)<-c(colnames(data1$x1)[data1$l1x],colnames(covariates)[cy2])
#tt1<-cbind(data1$x1[,data1$l2x],covariates[,cy1])
#colnames(tt1)<-c(colnames(data1$x1)[data1$l2x],colnames(covariates)[cy1])
#expl<-cbind(tt1,data1$m2y,tt2,data1$m1y)   

n<-length(y)
mnames<-colnames(as.matrix(m))
if(is.character(l1))
  l1<-unlist(sapply(l1,grep,mnames))
if(is.character(c1))
  c1<-unlist(sapply(c1,grep,mnames))
if(is.character(l2))
  l2<-unlist(sapply(l2,grep,mnames))
if(is.character(c2))
  c2<-unlist(sapply(c2,grep,mnames))

#organize the data if it has not been done
if(org.data)
  data1<-data.org(x, levelx, levely, m, l1,l2, c1, c1r,c2, c2r, f01y, f10y,
                  f02ky, f20ky, f01km1, f01km2, f10km, level, weight)
      
if(!is.null(data1$l1x))  #l1x could be NULL when x is level 2
 {tt2<-as.matrix(data1$x1[,data1$l1x])
  colnames(tt2)<-colnames(data1$x1)[data1$l1x]}
else
  tt2<-NULL             
tt1<-as.matrix(data1$x1[,data1$l2x])
colnames(tt1)<-colnames(data1$x1)[data1$l2x]
if(!is.null(c(cy1,cy2)))           ###added covariates to explain y
 {cova<-as.matrix(covariates[,c(cy1,cy2)])
  colnames(cova)<-colnames(covariates)[c(cy1,cy2)]}
else
  cova<-NULL
expl<-cbind(tt1,data1$m2y,tt2,data1$m1y,cova) 

temp.data<-cbind(y=y,level=level,data1$x1,data1$m1y,data1$m2y, covariates)
colnames(temp.data)<-c("y","level",colnames(data1$x1),colnames(data1$m1y),
                       colnames(data1$m2y), colnames(covariates))
if(levely==1)
 {frml<-getformula(expl,random,intercept)
  if(biny)
    f1<-glmer(frml,data=data.frame(temp.data),family=binomial(link="logit"))
  else
    f1<-lmer(frml,data=data.frame(temp.data))}
else
  {temp.data2<-two(temp.data, level, weight)
   frml<-getformula(expl,random=NULL,intercept)
   if(biny)
    f1<-glm(frml,data=data.frame(temp.data2),family=binomial(link="logit"))
   else
    f1<-lm(frml,data=data.frame(temp.data2))}

lc1<-c(l1,c1)
lc2<-c(l2,c2)

if(intercept)
  coef.f1<-summary(f1)$coefficient[-1,1] 
else 
  coef.f1<-summary(f1)$coefficient[,1]
len<-c(length(data1$l2x),ifelse(is.null(data1$m2y),0,ncol(data1$m2y)),
       ifelse(is.null(data1$l1x),0,length(data1$l1x)), 
       ifelse(is.null(data1$m1y),0,ncol(data1$m1y)))
if (len[1]>1)                       #calculate level 2 direct effect
  DE2<-data1$x1.der[,data1$l2x]%*%coef.f1[1:length(data1$l2x)]
else
  DE2<-data1$x1.der[,data1$l2x]*coef.f1[1]
if (len[2]==0)
{ie2_1<-NULL
 ie2_list<-NULL}
else 
{ie2_1<-NULL
 ie12_1<-NULL
 ie2_list<-NULL
 z<-0
 if (!is.null(data1$m2))
 {if(!is.null(l2))
 {for (k in 1:length(l2))
 {if (length(data1$m2[[k+1]])==1)
   ie2_1<-cbind(ie2_1,data1$m2y.der[,data1$m2[[k+1]]]*coef.f1[len[1]+data1$m2[[k+1]]])
  else
    ie2_1<-cbind(ie2_1,data1$m2y.der[,data1$m2[[k+1]]]%*%coef.f1[len[1]+data1$m2[[k+1]]])}
 colnames(ie2_1)<-colnames(m)[l2]
 ie2_list<-list(colnames(m)[l2])
 ie2_list<-append(ie2_list,1:length(l2))
 z<-length(l2)}
 if(!is.null(c2))
   for (i in 1:length(c2)) 
   {temp.name<-colnames(ie2_1)
    temp<-(1:length(data1$m2[[1]]))[data1$m2[[1]]==c2[i]]
    if (length(data1$m2[[temp + 1]])==1)
        ie2_1 <- cbind(ie2_1, data1$m2y.der[, data1$m2[[temp + 1]]] * 
                         coef.f1[len[1] + data1$m2[[temp + 1]]])
    else
       for (k in 1:length(data1$m2[[temp + 1]])) 
        ie2_1 <- cbind(ie2_1, data1$m2y.der[, data1$m2[[temp + 1]]][, k] * 
                         coef.f1[len[1] + data1$m2[[temp + 1]][k]])
    k<-length(data1$m2[[temp + 1]])
    colnames(ie2_1)<-c(temp.name,paste(colnames(m)[c2[i]],1:k,sep="."))
    ie2_list[[1]]<-c(ie2_list[[1]],colnames(m)[c2[i]])
    ie2_list<-append(ie2_list,list((z+1):(z+k)))
    z<-z+k
   }  }        
 if (!is.null(data1$m12))
 {if(!is.null(l1))
   {for (k in 1:length(l1))
   {if (length(data1$m12[[k+1]])==1)
     ie12_1<-cbind(ie12_1,data1$m2y.der[,data1$m12[[k+1]]]*coef.f1[len[1]+data1$m12[[k+1]]])
    else
      ie12_1<-cbind(ie12_1,data1$m2y.der[,data1$m12[[k+1]]]%*%coef.f1[len[1]+data1$m12[[k+1]]])}
    colnames(ie12_1)<-paste(colnames(m)[l1],"12",sep=".")}
  if(!is.null(c1))
    for (i in 1:length(c1)) 
    {temp.name<-colnames(ie12_1)
     temp<-(1:length(data1$m12[[1]]))[data1$m12[[1]]==c1[i]]
     if(length(data1$m12[[temp + 1]])==1)
       ie12_1 <- cbind(ie12_1,data1$m2y.der[, data1$m12[[temp + 1]]] * coef.f1[len[1] + data1$m12[[temp +1]]])
     else
       for (k in 1:length(data1$m12[[temp + 1]])) 
         ie12_1 <- cbind(ie12_1,data1$m2y.der[, data1$m12[[temp+1]]][,k]*coef.f1[len[1]+data1$m12[[temp +1]][k]])
     k<-length(data1$m12[[temp + 1]])
     colnames(ie12_1)<-c(temp.name,paste(colnames(m)[c1[i]],"12",1:k,sep="."))
    }          
 }
}
if (len[3]==0)
  DE1<-NULL
else if (len[3]==1)
  DE1<-data1$x1.der[,data1$l1x]*coef.f1[len[1]+len[2]+1]
else
  DE1<-data1$x1.der[,data1$l1x]%*%coef.f1[(len[1]+len[2]+1):(len[1]+len[2]+len[3])]
if (len[4]==0)
{ie1_1<-NULL
 ie1_list<-NULL}
else 
{ie1_1<-NULL
 ie1_list<-NULL
 z<-0
 if (!is.null(data1$m1))
 {if (!is.null(l1))
 {for (k in 1:length(l1))
 {if(length(data1$m1[[k+1]])==1)
   ie1_1<-cbind(ie1_1,data1$m1y.der[,data1$m1[[k+1]]]*coef.f1[len[1]+len[2]+len[3]+data1$m1[[k+1]]])
  else
    ie1_1<-cbind(ie1_1,data1$m1y.der[,data1$m1[[k+1]]]%*%coef.f1[len[1]+len[2]+len[3]+data1$m1[[k+1]]])}
 colnames(ie1_1)<-colnames(m)[l1]
 ie1_list<-list(colnames(m)[l1])
 ie1_list<-append(ie1_list,1:length(l1))
 z<-length(l1)
 }
 if(!is.null(c1))
   for (i in 1:length(c1)) 
   {temp<-(1:length(data1$m1[[1]]))[data1$m1[[1]]==c1[i]]
    temp.name<-colnames(ie1_1)
    if (length(data1$m1[[temp + 1]])==1)
      ie1_1 <- cbind(ie1_1,data1$m1y.der[, data1$m1[[temp + 1]]] * 
                       coef.f1[len[1] + len[2] + len[3] + data1$m1[[temp + 1]]])
    else
      for (k in 1:length(data1$m1[[temp + 1]])) 
        ie1_1 <- cbind(ie1_1,data1$m1y.der[, data1$m1[[temp + 1]]][, k] * 
                         coef.f1[len[1] + len[2] + len[3] + data1$m1[[temp + 1]][k]])
    k<-length(data1$m1[[temp + 1]])
    colnames(ie1_1)<-c(temp.name,paste(colnames(m)[c1[i]],1:k,sep="."))
    ie1_list[[1]]<-c(ie1_list[[1]],colnames(m)[c1[i]])
    ie1_list<-append(ie1_list,list((z+1):(z+k)))
    z<-z+k
   }  } }



if(levelx==2)    #analysis when x is a level 2 variable
{fm1<-list(NULL)            #models for x to explain level 1 mediators
 ie12_2<-NULL
 if(!is.null(l1))
 {fm1[[1]]<-l1
  for (i in 1:length(l1))
  {if(sum(cm[[1]]==l1[i])>0)
   {temp2<-(1:length(cm[[1]]))[cm[[1]]==l1[i]]+1
    temp.cov<-as.matrix(covariates[,cm[[temp2]]])
    colnames(temp.cov)<-colnames(covariates)[cm[[temp2]]]}
   else
    temp.cov<-NULL
   expl.m<-as.matrix(data1$xm1[,data1$fm12[[i+1]]])
   colnames(expl.m)<-colnames(data1$xm1)[data1$fm12[[i+1]]]
   numx<-length(data1$fm12[[i+1]]) #1:
   expl.m<-cbind(expl.m,temp.cov)
   frml.m<-getformula(expl.m,random,intercept)
   temp.data<-cbind(y=m[,l1[i]],level=level,expl.m)
   colnames(temp.data)<-c("y","level",colnames(expl.m))
   model<-lmer(frml.m,data=data.frame(temp.data))
   fm1<-append(fm1,list(model))
   coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
   if(intercept)
      coef.temp<-coef.temp[-1]
   coef.temp<-coef.temp[1:numx]
   if(length(coef.temp)==1)
     ie12_2<-cbind(ie12_2,coef.temp*data1$xm1.der[,data1$fm12[[i+1]]])
   else
     ie12_2<-cbind(ie12_2,data1$xm1.der[,data1$fm12[[i+1]]]%*%coef.temp)
  }
  j<-i+1
  colnames(ie12_2)<-colnames(m)[l1]
 }
 else
   j<-1
 if(!is.null(c1))
   for (i in 1:length(c1))
   {if(sum(cm[[1]]==c1[i])>0)
    {temp2<-(1:length(cm[[1]]))[cm[[1]]==c1[i]]+1
     temp.cov<-as.matrix(covariates[,cm[[temp2]]])
     colnames(temp.cov)<-colnames(covariates)[cm[[temp2]]]}
    else
     temp.cov<-NULL
    k1<-(1:length(data1$fm12[[1]]))[data1$fm12[[1]]==c1[i]]
    expl.m<-as.matrix(data1$xm1[,data1$fm12[[k1+1]]])
    colnames(expl.m)<-colnames(data1$xm1)[data1$fm12[[k1+1]]]
    numx<-length(data1$fm12[[k1+1]]) #1:
    expl.m<-cbind(expl.m,temp.cov)
    frml.m<-getformula(expl.m,random,intercept)
    temp<-(1:length(data1$m1[[1]]))[data1$m1[[1]]==c1[i]]
    name.temp<-colnames(ie12_2)
    if(length(data1$m1[[temp+1]])>1)######
    {temp.4<-apply(data1$m1y[,data1$m1[[temp+1]]],1,sum)
     temp.5<-ifelse(temp.4==1,F,T)}
    else 
      temp.5<-rep(T,n)
    for (k in data1$m1[[temp+1]])
    {temp.6<-(temp.5 | (data1$m1y[,k]==1))
     temp.data<-cbind(y=data1$m1y[temp.6,k],level=level[temp.6],expl.m[temp.6,]) #temp.6
     colnames(temp.data)<-c("y","level",colnames(expl.m))
     model<-glmer(frml.m,data=data.frame(temp.data),
                  family=binomial(link="logit"))
     temp.data<-cbind(level=level,expl.m) #temp.6
     colnames(temp.data)<-c("level",colnames(expl.m))
     p.temp<-predict(model,type="response",newdata=data.frame(temp.data),allow.new.levels=T)  ####################
     fm1<-append(fm1,model)
     fm1[[1]]<-c(fm1[[1]],c1[i])
     coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
     if(intercept)
        coef.temp<-coef.temp[-1]
     coef.temp<-coef.temp[1:numx]
     if(length(coef.temp)==1)
       ie12_2<-cbind(ie12_2,p.temp*(1-p.temp)*
                       (coef.temp*data1$xm1.der[,data1$fm12[[k1+1]]]))
     else
       ie12_2<-cbind(ie12_2,p.temp*(1-p.temp)*
                       (data1$xm1.der[,data1$fm12[[k1+1]]]%*%coef.temp))
    }
    colnames(ie12_2)<-c(name.temp,paste(colnames(m)[c1[i]],
                                        1:length(data1$m1[[temp+1]]),sep="."))
    j<-j+1
   }  
 fm2<-list(NULL)            #models for x to explain level 2 mediators
 ie2_2<-NULL
 if(!is.null(covariates))
   cov.2<-two(covariates,level,weight)
 else
   cov.2<-NULL
 if(!is.null(l2))
 {fm2[[1]]<-l2
  for (i in 1:length(l2))
  {if(sum(cm[[1]]==l2[i])>0)
   {temp2<-(1:length(cm[[1]]))[cm[[1]]==l2[i]]+1
    temp.cov<-as.matrix(cov.2[,cm[[temp2]]])
    colnames(temp.cov)<-colnames(covariates)[cm[[temp2]]]}
   else
    temp.cov<-NULL
   expl.m<-as.matrix(data1$xm2[,data1$fm22[[i+1]]])
   colnames(expl.m)<-colnames(data1$xm2)[data1$fm22[[i+1]]]
   numx<-ncol(expl.m)
   expl.m<-cbind(expl.m,temp.cov)
   frml.m<-getformula(expl.m,random=NULL,intercept)
   temp.data<-cbind(y=data1$m.2[,i],expl.m)
   colnames(temp.data)<-c("y",colnames(expl.m))
   model<-lm(frml.m,data=data.frame(temp.data))
   fm2<-append(fm2,list(model))
   coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
   if(intercept)
     coef.temp<-coef.temp[-1]
   coef.temp<-coef.temp[1:numx]
   if(length(coef.temp)==1)
     ie2_2<-cbind(ie2_2,coef.temp*data1$xm2.der[,data1$fm22[[i+1]]])
   else
     ie2_2<-cbind(ie2_2,data1$xm2.der[,data1$fm22[[i+1]]]%*%coef.temp)
  }
  j<-i+1
  colnames(ie2_2)<-colnames(m)[l2]
 }
 else
   j<-1
 if(!is.null(c2))    #did not test
   for (i in 1:length(c2))
   {if(sum(cm[[1]]==c2[i])>0)
    {temp2<-(1:length(cm[[1]]))[cm[[1]]==c2[i]]+1
     temp.cov<-as.matrix(cov.2[,cm[[temp2]]])
     colnames(temp.cov)<-colnames(covariates)[cm[[temp2]]]}
    else
     temp.cov<-NULL
    expl.m<-as.matrix(data1$xm2[,data1$fm22[[j+1]]])
    colnames(expl.m)<-colnames(data1$xm2)[data1$fm22[[j+1]]]
    numx<-ncol(expl.m)
    expl.m<-cbind(expl.m,temp.cov)
    frml.m<-getformula(expl.m,random=NULL,intercept)
    name.temp<-colnames(ie2_2)
    temp<-(1:length(data1$m2[[1]]))[data1$m2[[1]]==c2[i]]
    if(length(data1$m2[[temp+1]])>1)  ####
    {temp.3<-apply(data1$m2y[,data1$m2[[temp+1]]],2,two,level)
     temp.4<-apply(temp.3,1,sum)
     temp.5<-(temp.4==0)}
    else temp.5<-rep(T,length(unique(level)))
    for (k in data1$m2[[temp+1]])
    {temp.6<-(temp.5 | (two(data1$m2y[,k],level)==1))
     temp.data<-cbind(y=two(data1$m2y[temp.6,k],level=level[temp.6]),expl.m[temp.6,]) #
     colnames(temp.data)<-c("y",colnames(expl.m)) #"level",
     model<-glm(frml.m,data=data.frame(temp.data),
                family=binomial(link="logit"))
     temp.data<-expl.m
     #colnames(temp.data)<-colnames(data1$xm2)[data1$fm22[[j+1]]]
     p.temp<-predict(model,type="response",newdata=data.frame(temp.data))
     fm2<-append(fm2,list(model))
     fm2[[1]]<-c(fm2[[1]],c2[i])
     coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
     if(intercept)
       coef.temp<-coef.temp[-1]
     coef.temp<-coef.temp[1:numx]
     if(length(coef.temp)==1)
       ie2_2<-cbind(ie2_2,p.temp*(1-p.temp)*coef.temp*data1$xm2.der[,data1$fm22[[i+1]]])
     else
       ie2_2<-cbind(ie2_2,p.temp*(1-p.temp)*
                      (data1$xm2.der[,data1$fm22[[i+1]]]%*%coef.temp))
    }
    colnames(ie2_2)<-c(name.temp,paste(colnames(m)[c2[i]],
                                       1:length(data1$m2[[temp+1]]),sep="."))
    j<-j+1
   }
 ie1_2<-NULL}
else   #analysis when x is a level 1 variable
{fm1<-list(NULL)            #models for x to explain level 1 mediators
 ie1_2<-NULL
 ie12_2<-NULL
if(!is.null(l1))
 {fm1[[1]]<-l1
  for (i in 1:length(l1))
  {if(sum(cm[[1]]==l1[i])>0) #add covariates
   {temp2<-(1:length(cm[[1]]))[cm[[1]]==l1[i]]+1
    temp.cov<-as.matrix(covariates[,cm[[temp2]]])
    colnames(temp.cov)<-colnames(covariates)[cm[[temp2]]]}
   else
    temp.cov<-NULL
   expl.m<-as.matrix(data1$xm1[,c(data1$fm12[[i+1]],data1$fm11[[i+1]])])
   colnames(expl.m)<-colnames(data1$xm1)[c(data1$fm12[[i+1]],data1$fm11[[i+1]])]
   numx<-ncol(expl.m)
   expl.m<-cbind(expl.m,temp.cov)
   m.random<-"(1|level)"
   if(!is.null(random.m1))
     if(sum(random.m1[[1]]==l1[i])>0)
       m.random<-random.m1[[2]][random.m1[[1]]==l1[i]] #1st item of random.m1 is the list of l1 med, 2nd item is the random item of the same order
   frml.m<-getformula(expl.m,m.random,intercept)
   temp.data<-cbind(y=m[,l1[i]],level=level,expl.m)
   colnames(temp.data)<-c("y","level",colnames(expl.m))
   model=lmer(frml.m,data=data.frame(temp.data))
   fm1<-append(fm1,model)
   coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
   if(intercept)
     coef.temp<-coef.temp[-1]
   coef.temp<-coef.temp[1:numx]
   coef.temp12<-coef.temp[1:length(data1$fm12[[i+1]])]
   coef.temp1<-coef.temp[-(1:length(data1$fm12[[i+1]]))]
   if(length(coef.temp12)==1)
     ie12_2<-cbind(ie12_2,coef.temp12*data1$xm1.der[,data1$fm12[[i+1]]])
   else
     ie12_2<-cbind(ie12_2,data1$xm1.der[,data1$fm12[[i+1]]]%*%coef.temp12)
   if(length(coef.temp1)==1)
     ie1_2<-cbind(ie1_2,coef.temp1*data1$xm1.der[,data1$fm11[[i+1]]])
   else
     ie1_2<-cbind(ie1_2,data1$xm1.der[,data1$fm11[[i+1]]]%*%coef.temp1)
  }
  j<-i+1
  colnames(ie12_2)<-colnames(m)[l1]
  colnames(ie1_2)<-colnames(m)[l1]
 }
 else
   j<-1
 if(!is.null(c1))  
   for (i in 1:length(c1))
   {if(sum(cm[[1]]==c1[i])>0)   #add covariates
     {temp2<-(1:length(cm[[1]]))[cm[[1]]==c1[i]]+1
      temp.cov<-as.matrix(covariates[,cm[[temp2]]])
      colnames(temp.cov)<-colnames(covariates)[cm[[temp2]]]}
     else
      temp.cov<-NULL
     k1<-NULL                        #is there level 1 predictor?
     if(levelx==1)
       k1<-(1:length(data1$fm11[[1]]))[data1$fm11[[1]]==c1[i]]
     k2<-(1:length(data1$fm12[[1]]))[data1$fm12[[1]]==c1[i]]
    expl.m<-as.matrix(data1$xm1[,c(data1$fm12[[k2+1]],data1$fm11[[k1+1]])])
    colnames(expl.m)<-colnames(data1$xm1)[c(data1$fm12[[k2+1]],data1$fm11[[k1+1]])]
    numx<-ncol(expl.m)
    expl.m<-cbind(expl.m,temp.cov)
    m.random<-"(1|level)"
    if(!is.null(random.m1))
      if(sum(random.m1[[1]]==c1[i])>0)
        m.random<-random.m1[[2]][random.m1[[1]]==c1[i]] #1st item of random.m1 is the list of l1 med, 2nd item is the random item of the same order
    frml.m<-getformula(expl.m,m.random,intercept)
    name.temp<-colnames(ie12_2)
    temp<-(1:length(data1$m1[[1]]))[data1$m1[[1]]==c1[i]]
    if(length(data1$m1[[temp+1]])>1)###
    {temp.4<-apply(data1$m1y[,data1$m1[[temp+1]]],1,sum)
     temp.5<-ifelse(temp.4==1,F,T)}
    else temp.5<-rep(T,n)
    for (k in data1$m1[[temp+1]])
    {temp.6<-(temp.5 | (data1$m1y[,k]==1))
     temp.data<-cbind(y=data1$m1y[temp.6,k],level=level[temp.6],expl.m[temp.6,])  #[temp.6]
     colnames(temp.data)<-c("y","level", colnames(expl.m))
     model=glmer(frml.m,data=data.frame(temp.data),
                 family=binomial(link="logit"))
     fm1<-append(fm1,model)
     fm1[[1]]<-c(fm1[[1]],c1[i])
     temp.data<-cbind(level=level,expl.m)  #[temp.6]
     colnames(temp.data)<-c("level",colnames(expl.m))
     p.temp<-predict(model,type="response",newdata=data.frame(temp.data),allow.new.levels=T)
     coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
     if(intercept)
        coef.temp<-coef.temp[-1]
     coef.temp<-coef.temp[1:numx]
     coef.temp12<-coef.temp[1:length(data1$fm12[[k2+1]])]
     coef.temp1<-coef.temp[-(1:length(data1$fm12[[k2+1]]))]
     if(length(coef.temp12)==1)
       ie12_2<-cbind(ie12_2,coef.temp12*data1$xm1.der[,data1$fm12[[k2+1]]])
     else
       ie12_2<-cbind(ie12_2,data1$xm1.der[,data1$fm12[[k2+1]]]%*%coef.temp12)
     if(length(coef.temp1)==1)
       ie1_2<-cbind(ie1_2,p.temp*(1-p.temp)*coef.temp1*data1$xm1.der[,data1$fm11[[k1+1]]])
     else
       ie1_2<-cbind(ie1_2,p.temp*(1-p.temp)*(data1$xm1.der[,data1$fm11[[k1+1]]]%*%coef.temp1))
    }
    colnames(ie12_2)<-c(name.temp,paste(colnames(m)[c1[i]],
                                        1:length(data1$m1[[temp+1]]),sep="."))
    colnames(ie1_2)<-c(name.temp,paste(colnames(m)[c1[i]],
                                       1:length(data1$m1[[temp+1]]),sep="."))
    j<-j+1
   }

##################added aggregated level 1 x to explain level 2 mediators
 
 fm2<-list(NULL)            #models for x to explain level 2 mediators
 ie2_2<-NULL
 if(!is.null(covariates))
   cov.2<-two(covariates,level, weight)
 else 
   cov.2<-NULL
 if(!is.null(l2))
 {fm2[[1]]<-l2
 for (i in 1:length(l2))
 {if(sum(cm[[1]]==l2[i])>0)   #add covariates
  {temp2<-(1:length(cm[[1]]))[cm[[1]]==l2[i]]+1
   temp.cov<-as.matrix(cov.2[,cm[[temp2]]])
   colnames(temp.cov)<-colnames(covariates)[cm[[temp2]]]}
  else
   temp.cov<-NULL
 expl.m<-as.matrix(data1$xm2[,data1$fm22[[i+1]]])
 colnames(expl.m)<-colnames(data1$xm2)[data1$fm22[[i+1]]]
 numx<-ncol(expl.m)
 expl.m<-cbind(expl.m, temp.cov)
 frml.m<-getformula(expl.m,random=NULL,intercept)
 temp.data<-cbind(y=data1$m.2[,i],expl.m)
 colnames(temp.data)<-c("y",colnames(expl.m))
 model<-lm(frml.m,data=data.frame(temp.data))
 fm2<-append(fm2,list(model))
 coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
 if(intercept)
   coef.temp<-coef.temp[-1]
 coef.temp<-coef.temp[1:numx]
 if(length(coef.temp)==1)
   ie2_2<-cbind(ie2_2,coef.temp*data1$xm2.der[,data1$fm22[[i+1]]])
 else
   ie2_2<-cbind(ie2_2,data1$xm2.der[,data1$fm22[[i+1]]]%*%coef.temp)
 }
 j<-i+1
 colnames(ie2_2)<-colnames(m)[l2]
 }
 else
   j<-1
 if(!is.null(c2))    #did not test
   for (i in 1:length(c2))
   {if(sum(cm[[1]]==c2[i])>0)   #add covariates
    {temp2<-(1:length(cm[[1]]))[cm[[1]]==c2[i]]+1
     temp.cov<-as.matrix(cov.2[,cm[[temp2]]])
     colnames(temp.cov)<-colnames(covariates)[cm[[temp2]]]}
    else
     temp.cov<-NULL
   expl.m<-as.matrix(data1$xm2[,data1$fm22[[j+1]]])
   colnames(expl.m)<-colnames(data1$xm2)[data1$fm22[[j+1]]]
   numx<-ncol(expl.m)
   expl.m<-cbind(expl.m, temp.cov)
   frml.m<-getformula(expl.m,random=NULL,intercept)
   name.temp<-colnames(ie2_2)
   temp<-(1:length(data1$m2[[1]]))[data1$m2[[1]]==c2[i]]
   if(length(data1$m2[[temp+1]])>1)  ####
   {temp.3<-apply(data1$m2y[,data1$m2[[temp+1]]],2,two,level)
   temp.4<-apply(temp.3,1,sum)
   temp.5<-(temp.4==0)}
   else temp.5<-rep(T,length(unique(level)))
   for (k in data1$m2[[temp+1]])
   {temp.6<-(temp.5 | (two(data1$m2y[,k],level)==1))
   temp.data<-cbind(y=two(data1$m2y[temp.6,k],level[temp.6]),expl.m[temp.6,])
   colnames(temp.data)<-c("y",colnames(expl.m)) #"level",
   model<-glm(frml.m,data=data.frame(temp.data),
              family=binomial(link="logit"))
   temp.data<-expl.m
   #colnames(temp.data)<-colnames(data1$xm2)[data1$fm22[[j+1]]]
   p.temp<-predict(model,type="response",newdata=data.frame(temp.data))
   fm2<-append(fm2,list(model))
   fm2[[1]]<-c(fm2[[1]],c2[i])
   coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
   if(intercept)
     coef.temp<-coef.temp[-1]
   coef.temp<-coef.temp[1:numx]
   if(length(coef.temp)==1)
     ie2_2<-cbind(ie2_2,p.temp*(1-p.temp)*coef.temp*data1$xm2.der[,data1$fm22[[i+1]]])
   else
     ie2_2<-cbind(ie2_2,p.temp*(1-p.temp)*
                    (data1$xm2.der[,data1$fm22[[i+1]]]%*%coef.temp))
   }
   colnames(ie2_2)<-c(name.temp,paste(colnames(m)[c2[i]],
                                      1:length(data1$m2[[temp+1]]),sep="."))
   j<-j+1
   }
}

if(!is.null(ie12_1) & !is.null(ie12_2))
{ie12_1<-apply(as.matrix(ie12_1),2,two,level)
 ie12_2<-apply(as.matrix(ie12_2),2,two,level)
 ie12<-ie12_1*ie12_2
}
else ie12<-NULL
if(!is.null(ie1_1) & !is.null(ie1_2))
  ie1<-ie1_1*ie1_2
else ie1<-NULL
if(!is.null(ie2_1) & !is.null(ie2_2))
{ie2_1<-two(ie2_1,level)
 ie2<-ie2_1*ie2_2
}
else ie2<-NULL
if(!is.null(joint))  #3/2/2016, add the potential joint effect
{joint1<-(1:length(joint[[1]]))[joint[[1]]==1]
 joint2<-(1:length(joint[[1]]))[joint[[1]]==2]
 if(length(joint1)!=0)
 {temp.name<-ie1_list[[1]]
  for (i in joint1)
  {if(!is.character(joint[[i+1]]))
    joint[[i+1]]<-mnames[joint[[i+1]]]
    joint.col<-sapply(joint[[i+1]],grep,temp.name)+1
   joint.col1<-NULL
   for (j in 1:length(joint.col))
     joint.col1<-c(joint.col1,ie1_list[[joint.col[j]]])
   ie1_list[[length(ie1_list)+1]]<-joint.col1
  }
  ie1_list[[1]]<-c(ie1_list[[1]],paste("joint",joint1,sep="."))
 }
 if(length(joint2)!=0)
 {temp.name<-ie1_list[[1]]
  temp.name2<-ie2_list[[1]]
  temp.l1<-ifelse(is.null(ie12),0,ncol(ie12))
  iej_list<-list(paste("joint",joint2,sep="."))
  for(i in joint2)
   {if(!is.character(joint[[i+1]]))
     joint[[i+1]]<-mnames[joint[[i+1]]]
    joint.col1<-unlist(sapply(joint[[i+1]],grep,temp.name))+1
    joint.col2<-unlist(sapply(joint[[i+1]],grep,temp.name2))+1
    joint.col<-NULL
    if(length(joint.col1)!=0)
     for (j in 1:length(joint.col1))
      joint.col<-c(joint.col,ie1_list[[joint.col1[j]]])
    if(length(joint.col2)!=0)
     for (j in 1:length(joint.col2))
      joint.col<-c(joint.col,ie2_list[[joint.col2[j]]]+temp.l1)
    iej_list[[length(iej_list)+1]]<-joint.col}
 }
 else
   iej_list<-NULL
}
else
  iej_list<-NULL
a<-list(de1=DE1,de2=two(DE2,level),ie1=ie1,ie2=ie2,ie12=ie12,
        f1=f1,fm1=fm1,fm2=fm2, ie1_list=ie1_list, ie2_list=ie2_list,
        iej2_list=iej_list,ie12_1=ie12_1, ie12_2=ie12_2, ie1_1=ie1_1,
        ie1_2=ie1_2,ie2_1=ie2_1,ie2_2=ie2_2,x=x, x.j=two(x,level), 
        m=m, covariates=covariates, intercept=intercept, cm=cm,
        data1=data1)
class(a)<-"mlma"
return(a)
}

boot.mlma<-function(y, biny=FALSE, x, levelx=1,levely=1, m, l1=NULL,l2=NULL, c1=NULL, 
                    c1r=rep(1,length(c1)), #levelx is the level of x
                    c2=NULL, c2r=rep(1,length(c2)), f01y=NULL, f10y=NULL,                  #level is the level of observations
                    f02ky=NULL, f20ky=NULL, f01km1=NULL, f01km2=NULL, f10km=NULL,          #weight is the level 1 weight of cases
                    level=1:length(y), weight=rep(1,length(x)), random="(1|level)",
                    random.m1=NULL, intercept=TRUE, 
                    w2=rep(1,length(unique(level[!is.na(level)]))),
                    boot=100,seed=1,covariates=NULL,cy1=NULL, cy2=NULL, cm=NULL, 
                    joint=NULL, x.new=x, m.new=m, level.new=level,weight.new=weight,
                    covariates.new=covariates)
{    two <- function(x, level, weight = rep(1, nrow(as.matrix(x)))) {
  x <- as.matrix(x)
  levels <- unique(level[!is.na(level)])
  x2 <- matrix(NA, length(levels), dim(x)[2])
  for (i in 1:length(levels)) {
    cho <- (level == levels[i])
    if (sum(cho) > 0) {
      if(sum(cho)==1)
        x2[i, ] <- x[level == levels[i], ]
      else {
        temp <- as.matrix(x[level == levels[i], ])
        weight1 <- weight[level == levels[i]]
        x2[i, ] <- apply(temp, 2, weighted.mean, weight1, 
                         na.rm = TRUE)}
    }
  }
  colnames(x2) <- colnames(x)
  x2
}


sum.mlma<-function(object,levelx,de2,ie2,de1,ie1,ie12,level,weight,w2)
{two<-function(x, level, weight=rep(1,nrow(as.matrix(x))))
{x<-as.matrix(x)
levels<-unique(level[!is.na(level)])
x2<-matrix(NA,length(levels),dim(x)[2])
for(i in 1:length(levels))
{cho<-(level==levels[i])
if(sum(cho)>0)
{if(sum(cho)==1)
  x2[i,]<-x[level==levels[i],]
else
{temp<-as.matrix(x[level==levels[i],])
weight1<-weight[level==levels[i]]
x2[i,]<-apply(temp,2,weighted.mean,weight1,na.rm=TRUE)}}}
colnames(x2)<-colnames(x)
x2
}


de2<-weighted.mean(de2,w2)
ie2_temp<-ie2
if(!is.null(ie2))
{ie2<-apply(ie2,2,weighted.mean,w2)
 ie2_show<-NULL
for (i in 1:(length(object$ie2_list)-1))
  ie2_show<-c(ie2_show,sum(ie2[object$ie2_list[[i+1]]]))
names(ie2_show)<-object$ie2_list[[1]]
}
else 
  {ie2=0
   ie2_show<-NULL}
ie12_temp<-ie12
if(!is.null(ie12))
{ie12<-apply(ie12_temp,2,weighted.mean,w2)
ie12_show<-NULL
for (i in 1:(length(object$ie1_list)-1))
  ie12_show<-c(ie12_show,sum(ie12[object$ie1_list[[i+1]]]))
names(ie12_show)<-object$ie1_list[[1]]
}
else 
  {ie12=0
   ie12_show<-NULL}

if(!is.null(object$iej2_list)) #added 3/2/2016 to show level 2 joint effects
{ie_all2<-cbind(ie12_temp,ie2_temp)
ie_j2<-apply(ie_all2,2,weighted.mean,w2)
iej2_show<-NULL
for (i in 1:(length(object$iej2_list)-1))
  iej2_show<-c(iej2_show,sum(ie_j2[object$iej2_list[[i+1]]]))
names(iej2_show)<-object$iej2_list[[1]]
}
else
  iej2_show<-NULL

te2<-sum(c(de2,ie2,ie12))
if(levelx==1)
{de1<-two(de1,level,weight)
de1<-weighted.mean(de1,w2)
if(!is.null(ie1))
{ie1<-two(ie1,level,weight)
ie1<-apply(ie1,2,weighted.mean,w2)
ie1_show<-NULL
for (i in 1:(length(object$ie1_list)-1))
  ie1_show<-c(ie1_show,sum(ie1[object$ie1_list[[i+1]]]))
names(ie1_show)<-object$ie1_list[[1]]
} 
else ie1=0
te1<-sum(c(de1,ie1))
}
else
{de1<-NULL
 ie1_show<-NULL
 te1<-NULL
}
level2=c(de2=de2,ie2=ie2_show,ie12=ie12_show, iej2=iej2_show, te2=te2)
level1=c(de1=de1,ie1=ie1_show,te1=te1)
  return(list(level2=level2,level1=level1))
}

getformula<-function(expl,random="(1|level)",intercept=TRUE)
{temp.name<-colnames(expl)
formula<-temp.name[1]
if(ncol(expl)>1)
  for (i in 2:ncol(expl))
    formula<-paste(formula,temp.name[i],sep="+")
formula<-ifelse(intercept,formula,paste(formula,"1",sep="-"))
formula<-paste("y~",formula)
if (!is.null(random))
  formula<-paste(formula,random,sep="+")
formula
}

ns.dev<-function (x, df = NULL, knots = NULL, qnots=NULL,intercept = FALSE, Boundary.knots = range(x),derivs1=0) 
{
  nx <- names(x)
  x <- as.vector(x)
  nax <- is.na(x)
  if (nas <- any(nax)) 
    x <- x[!nax]
  if (!missing(Boundary.knots)) {
    Boundary.knots <- sort(Boundary.knots)
    outside <- (ol <- x < Boundary.knots[1L]) | (or <- x > 
                                                   Boundary.knots[2L])
  }
  else outside <- FALSE
  if (!is.null(df) && is.null(knots)) {
    nIknots <- df - 1L - intercept
    if (nIknots < 0L) {
      nIknots <- 0L
      warning(gettextf("'df' was too small; have used %d", 
                       1L + intercept), domain = NA)
    }
    knots <- if (nIknots > 0L) {
      knots <- seq.int(0, 1, length.out = nIknots + 2L)[-c(1L, 
                                                           nIknots + 2L)]
      stats::quantile(x[!outside], knots)
    }
  }
  else {if(is.null(df) && is.null(knots) && !is.null(qnots))
    knots<-quantile(x[!outside], qnots)
    nIknots <- length(knots)}
  Aknots <- sort(c(rep(Boundary.knots, 4L), knots))
  if (any(outside)) {
    basis <- array(0, c(length(x), nIknots + 4L))
    if (any(ol)) {
      k.pivot <- Boundary.knots[1L]
      xl <- cbind(1, x[ol] - k.pivot)
      tt <- splineDesign(Aknots, rep(k.pivot, 2L), 4, c(0, 
                                                        1),derivs=rep(derivs1,2L))
      basis[ol, ] <- xl %*% tt
    }
    if (any(or)) {
      k.pivot <- Boundary.knots[2L]
      xr <- cbind(1, x[or] - k.pivot)
      tt <- splineDesign(Aknots, rep(k.pivot, 2L), 4, c(0,1),derivs=rep(derivs1,2L))
      basis[or, ] <- xr %*% tt
    }
    if (any(inside <- !outside)) 
      basis[inside, ] <- splineDesign(Aknots, x[inside], 
                                      4,derivs=rep(derivs1,length(x[inside])))
  }
  else basis <- splineDesign(Aknots, x, 4,derivs=rep(derivs1,length(x)))
  const <- splineDesign(Aknots, Boundary.knots, 4, c(2, 2),derivs=rep(derivs1,length(Boundary.knots)))
  if (!intercept) {
    const <- const[, -1, drop = FALSE]
    basis <- basis[, -1, drop = FALSE]
  }
  qr.const <- qr(t(const))
  basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1L:2L), 
                                                     drop = FALSE])
  n.col <- ncol(basis)
  if (nas) {
    nmat <- matrix(NA, length(nax), n.col)
    nmat[!nax, ] <- basis
    basis <- nmat
  }
  dimnames(basis) <- list(nx, 1L:n.col)
  a <- list(degree = 3L, knots = if (is.null(knots)) numeric() else knots, 
            Boundary.knots = Boundary.knots, intercept = intercept)
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("ns", "basis", "matrix")
  basis
}


bs.dev<-function (x, df = NULL, knots = NULL, degree = 3, intercept = FALSE, 
                  Boundary.knots = range(x),derivs1=0) 
{
  nx <- names(x)
  x <- as.vector(x)
  nax <- is.na(x)
  if (nas <- any(nax)) 
    x <- x[!nax]
  if (!missing(Boundary.knots)) {
    Boundary.knots <- sort(Boundary.knots)
    outside <- (ol <- x < Boundary.knots[1L]) | (or <- x > 
                                                   Boundary.knots[2L])
  }
  else outside <- FALSE
  ord <- 1L + (degree <- as.integer(degree))
  if (ord <= 1) 
    stop("'degree' must be integer >= 1")
  if (!is.null(df) && is.null(knots)) {
    nIknots <- df - ord + (1L - intercept)
    if (nIknots < 0L) {
      nIknots <- 0L
      warning(gettextf("'df' was too small; have used %d", 
                       ord - (1L - intercept)), domain = NA)
    }
    knots <- if (nIknots > 0L) {
      knots <- seq.int(from = 0, to = 1, length.out = nIknots + 
                         2L)[-c(1L, nIknots + 2L)]
      stats::quantile(x[!outside], knots)
    }
  }  

  Aknots <- sort(c(rep(Boundary.knots, ord), knots))
  if (any(outside)) {
    warning("some 'x' values beyond boundary knots may cause ill-conditioned bases")
    derivs <- 0:degree
    scalef <- gamma(1L:ord)
    basis <- array(0, c(length(x), length(Aknots) - degree - 
                          1L))
    if (any(ol)) {
      k.pivot <- Boundary.knots[1L]
      xl <- cbind(1, outer(x[ol] - k.pivot, 1L:degree, 
                           "^"))
      tt <- splineDesign(Aknots, rep(k.pivot, ord), ord, 
                         derivs+derivs1)
      basis[ol, ] <- xl %*% (tt/scalef)
    }
    if (any(or)) {
      k.pivot <- Boundary.knots[2L]
      xr <- cbind(1, outer(x[or] - k.pivot, 1L:degree, 
                           "^"))
      tt <- splineDesign(Aknots, rep(k.pivot, ord), ord, 
                         derivs+derivs1)
      basis[or, ] <- xr %*% (tt/scalef)
    }
    if (any(inside <- !outside)) 
      basis[inside, ] <- splineDesign(Aknots, x[inside], 
                                      ord,derivs=rep(derivs1,length(x[inside])))
  }
  else basis <- splineDesign(Aknots, x, ord, derivs=rep(derivs1,length(x)))
  if (!intercept) 
    basis <- basis[, -1L, drop = FALSE]
  n.col <- ncol(basis)
  if (nas) {
    nmat <- matrix(NA, length(nax), n.col)
    nmat[!nax, ] <- basis
    basis <- nmat
  }
  dimnames(basis) <- list(nx, 1L:n.col)
  a <- list(degree = degree, knots = if (is.null(knots)) numeric(0L) else knots, 
            Boundary.knots = Boundary.knots, intercept = intercept)
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("bs", "basis", "matrix")
  basis
}

x2fx<-function(x,func) #x is the list of original numerical vector, func is a vector of character functions. 
{ # eg.  func <- c("x","x+1","x+2","x+3","log(x)")
  func.list <- list()
  #test.data <- matrix(data=rep(x,length(func)),length(x),length(func))
  #test.data <- data.frame(test.data)
  result<-NULL
  for(i in 1:length(func)){
    func.list[[i]] <- function(x){}
    body(func.list[[i]]) <- parse(text=func[i])
  }
  #result <- mapply(do.call,func.list,lapply(test.data,list))
  col_fun<-NULL
  z<-1
  for (i in 1:length(func.list))
  {res<-as.matrix(func.list[[i]](x))
  result<-cbind(result,res)
  col_fun<-cbind(col_fun,c(z,z+ncol(res)-1))
  z<-z+ncol(res)}
  list(values=as.matrix(result),col_fun=as.matrix(col_fun))
}

one2two<-function(x,l1,weight=rep(1,length(x))) #l1 is a vector that distribute x to different level 1 groups
{x1<-rep(NA,length(x))
l1.1<-l1[!is.na(l1)]
x.1<-x[!is.na(l1)]
weight.1<-weight[!is.na(l1)]
weight.1<-ifelse(is.na(weight.1),0,weight.1)   #missing weight will not be counted when calculate the level 1 variable
x1.1<-rep(0,length(x.1))
for (i in unique(l1.1))
  x1.1[l1.1==i]<-weighted.mean(x.1[l1.1==i],weight.1[l1.1==i],na.rm=TRUE)
x1[!is.na(l1)]<-x1.1
cbind(x1,x-x1)
}

x2fdx<-function(x,func)  #x is the list of original numerical vector, func is a vector of character functions. 
{ fdx<-NULL              # eg.  func <- c("x","x+1","x+2","x+3","log(x)")
for(i in 1:length(func)){
  if (length(grep("ifelse",func[i]))>0)
  {str<-unlist(strsplit(func[i],","))
  fun1<-D(parse(text=str[2]), "x")
  fun2<-D(parse(text=unlist(strsplit(str[3],")"))),"x")
  x1<-eval(fun1)
  x2<-eval(fun2)
  if(length(x1)==1)
    x1<-rep(x1,length(x))
  if(length(x2)==1)
    x2<-rep(x2,length(x))
  fun3<-paste(str[1],"x1,x2)",sep=",")
  fdx<-cbind(fdx,eval(parse(text=fun3)))
  }
  else if(length(grep("ns",func[i]))>0)
  {temp<-paste("ns.dev",substring(func[i],3,nchar(func[i])-1),",derivs1=1)",sep="")
  fdx<-cbind(fdx,eval(parse(text=temp)))}
  else if(length(grep("bs",func[i]))>0)
  {temp<-paste("bs.dev",substring(func[i],3,nchar(func[i])-1),",derivs1=1)",sep="")
  fdx<-cbind(fdx,eval(parse(text=temp)))}
  else{
    dx2x <- D(parse(text=func[i]), "x") 
    temp<-eval(dx2x)
    if(length(temp)==1)
      fdx<-cbind(fdx,rep(temp,length(x)))
    else fdx<-cbind(fdx,temp)}
}
as.matrix(fdx)
}

order_char<-function(char1,char2)  #find the position of char2 in char1
{a<-1:length(char1)
b<-NULL
for (i in 1:length(char2))
  b<-c(b,a[char1==char2[i]])
b
}

cattobin<-function(m1y,m1y.der,m1,m,cat1,cat2=rep(1,length(cat1)), level=rep(1,dim(m)[1]),weight=rep(1,dim(m)[1]))
{ ad1<-function(vec)
{vec1<-vec[-1]
vec1[vec[1]]<-1
vec1
}

m<-as.matrix(m)
if(is.null(m1y))
{m1<-list(cat1)
dim1<-c(dim(m)[1],0)}
else{
  dim1<-dim(m1y)
  m1[[1]]<-c(m1[[1]],cat1)}

m12y<-NULL
m12y.der<-NULL
m12<-list(cat1)
g<-dim1[2]
ntemp<-colnames(m)[cat1]
j<-1
p<-0
for (i in cat1)
{a<-m[,i]
d<-rep(0,dim1[1])
b<-sort(unique(a[a!=cat2[j]]))
l<-1
for (k in b)
{d[a==k]<-l
l<-l+1}
d[a==cat2[j]]<-l
f<-matrix(0,dim1[1],l-1) 
colnames(f)<-paste(ntemp[j],b,sep=".")  ##
hi<-d[d!=l & !is.na(d)]
f[d!=l & !is.na(d),]<-t(apply(cbind(hi,f[d!=l & !is.na(d),]),1,ad1))
f[is.na(d),]<-NA
z<-apply(f,2,one2two,level,weight)
m1y<-cbind(m1y,f)
m1y.der<-cbind(m1y.der,f)
m1<-append(m1,list((g+1):(g+l-1)))
temp3<-as.matrix(z[1:dim1[1],])
colnames(temp3)<-paste(ntemp[j],"12",b,sep=".")  ##
m12y<-cbind(m12y,temp3)
m12y.der<-cbind(m12y.der,matrix(1,nrow(temp3),ncol(temp3))) ## for changes per unit percent
m12<-append(m12,list((p+1):(p+l-1)))
g<-g+length(b)
p<-p+length(b)
j<-j+1
}
colnames(m12y.der)<-colnames(m12y)
list(m1y=m1y,m1y.der=m1y.der,m1=m1,m12y=m12y,m12y.der=m12y.der,m12=m12)
}

one<-function(x,level)  #change a level 2 x to have the length of observations
{levels<-unique(level[!is.na(level)])
x1<-rep(NA,length(level))
for (i in 1:length(levels))
  x1[level==levels[i]]=x[i]
x1
}

mnames<-colnames(as.matrix(m))
if(is.character(l1))
  l1<-unlist(sapply(l1,grep,mnames))
if(is.character(c1))
  c1<-unlist(sapply(c1,grep,mnames))
if(is.character(l2))
  l2<-unlist(sapply(l2,grep,mnames))
if(is.character(c2))
  c2<-unlist(sapply(c2,grep,mnames))

n2<-length(unique(level[!is.na(level)]))
#1a. an anlaysis on the whole data set
data1<-data.org(x, levelx, levely, m, l1, l2, c1, c1r, c2, c2r, f01y, f10y, f02ky, f20ky,
                f01km1, f01km2, f10km, level, weight)
#data2<-data.org(x.new, levelx, levely, m.new, l1, l2, c1, c1r, c2, c2r, f01y, f10y, f02ky, f20ky,
#                f01km1, f01km2, f10km, level.new, weight.new)
full<-mlma(y, biny, data1, x, levelx, levely, m, l1,l2, c1, c1r, 
           c2, c2r, level, weight, random, random.m1,intercept, 
           covariates, cy1, cy2, cm, joint)  

#1b. prepare data
m.new1<-m.new
if(is.null(m.new1))               #if the new m is to be generated
{m<-as.matrix(m)
  if(length(x.new)<=nrow(m))     # make a fake m.new
   m.new<-m[1:length(x.new),]
  else
  {a3<-length(x.new)%/%nrow(m)
   a4<-length(x.new)%%nrow(m)
   m.new<-m[rep(1:nrow(m),a3),]
   m.new<-rbind(m.new,m[1:a4,])}
 data2<-data.org(x.new, levelx, levely, m.new, l1, l2, c1, c1r, c2, c2r, f01y, f10y, f02ky, f20ky,
                 f01km1, f01km2, f10km, level.new, weight.new)
 xm1<-data2$xm1
 xm2<-data2$xm2
 xm1.der<-data2$xm1.der
 xm2.der<-data2$xm2.der
 x1.der<-data2$x1.der
}
else
{data2<-data.org(x.new, levelx, levely, m.new, l1, l2, c1, c1r, c2, c2r, f01y, f10y, f02ky, f20ky,
                 f01km1, f01km2, f10km, level.new, weight.new)
 xm1<-data2$xm1
 xm2<-data2$xm2
 xm1.der<-data2$xm1.der
 xm2.der<-data2$xm2.der
 x1.der<-data2$x1.der
 m1y.der<-data2$m1y.der
 m2y.der<-data2$m2y.der
}

#2. bootstrap sample and prepare for data
de1.boot<-NULL
de2.boot<-NULL
ie1.boot<-NULL
ie2.boot<-NULL
ie12.boot<-NULL
#xboot<-NULL
#xjboot<-NULL
sum.boot1<-NULL
sum.boot2<-NULL

#a. create x, y, m, covariates, level, x1, m1y, m2y, xm1
m<-as.matrix(m)
t1<-table(level)
t2<-sort(unique(level))
level.boot<-rep(t2,t1)
n1<-length(y)
for(l in 1:boot)
{set.seed(seed+l)
 if(length(unique(level))==n1)
   {s1<-sample(1:n1,replace=TRUE,prob=weight)
    x.boot<-x[s1]
    y.boot<-y[s1]
    m.boot<-as.matrix(as.matrix(m)[s1,])
    x1.boot<-as.matrix(as.matrix(data1$x1)[s1,])
    if(!is.null(data1$m1y))
      m1y.boot<-as.matrix(as.matrix(data1$m1y)[s1,])
    else
      m1y.boot<-NULL
    if(!is.null(data1$m2y))
      m2y.boot<-as.matrix(as.matrix(data1$m2y)[s1,])
    else
      m2y.boot<-NULL
    if(!is.null(data1$xm1))
      xm1.boot<-as.matrix(as.matrix(data1$xm1)[s1,])
    else
      xm1.boot<-NULL
    if(!is.null(covariates))
      covariates.boot<-as.matrix(as.matrix(covariates)[s1,])
    else
      covariates.boot<-NULL
    xm2.boot<-NULL
    m.2.boot<-NULL
    }
 else{
 x.boot<-NULL
 y.boot<-NULL
 m.boot<-NULL
 x1.boot<-NULL
 m2y.boot<-NULL
 m1y.boot<-NULL
 xm1.boot<-NULL
 xm2.boot<-NULL
 m.2.boot<-NULL
 weight.boot<-NULL
 covariates.boot<-NULL
 for(j in 1:length(t2))
 {temp<-level==t2[j]
  temp.2<-sample(t1[j],replace=TRUE,prob=weight[temp])
  x.boot<-c(x.boot,x[temp][temp.2])
  y.boot<-c(y.boot,y[temp][temp.2]) 
  if(ncol(m)==1)
    m.boot<-c(m.boot,m[temp][temp.2])
  else if(t1[j]==1)                       #if there is only one observation in a level
    m.boot<-rbind(m.boot,m[temp,]) 
  else
    m.boot<-rbind(m.boot,m[temp,][temp.2,]) 
  if(ncol(data1$x1)==1)
    x1.boot<-c(x1.boot,data1$x1[temp][temp.2])
  else if(t1[j]==1)
    x1.boot<-rbind(x1.boot,data1$x1[temp,]) 
  else
    x1.boot<-rbind(x1.boot,data1$x1[temp,][temp.2,]) 
  if(!is.null(data1$m1y))
  {if(ncol(data1$m1y)==1)
    m1y.boot<-c(m1y.boot,data1$m1y[temp][temp.2])
  else if(t1[j]==1)
    m1y.boot<-rbind(m1y.boot,data1$m1y[temp,]) 
  else
    m1y.boot<-rbind(m1y.boot,data1$m1y[temp,][temp.2,]) 
  }
  if(!is.null(data1$m2y))
  {if(ncol(data1$m2y)==1)
    m2y.boot<-c(m2y.boot,data1$m2y[temp][temp.2])
  else if(t1[j]==1)
    m2y.boot<-rbind(m2y.boot,data1$m2y[temp,]) 
  else
    m2y.boot<-rbind(m2y.boot,data1$m2y[temp,][temp.2,]) 
  }
  if(!is.null(data1$xm1))
  {if(ncol(data1$xm1)==1)
    xm1.boot<-c(xm1.boot,data1$xm1[temp][temp.2])
  else if(t1[j]==1)
    xm1.boot<-rbind(xm1.boot,data1$xm1[temp,]) 
  else
    xm1.boot<-rbind(xm1.boot,data1$xm1[temp,][temp.2,]) 
  }
  if(!is.null(covariates))
  {if(ncol(covariates)==1)
    covariates.boot<-c(covariates.boot,covariates[temp][temp.2])
  else if(t1[j]==1)
    covariates.boot<-rbind(covariates.boot,covariates[temp,]) 
  else
    covariates.boot<-rbind(covariates.boot,covariates[temp,][temp.2,]) 
  }
  weight.boot<-c(weight.boot,weight[temp][temp.2])
 }
 if(!is.null(dim(m.boot)))
 {m.boot<-as.matrix(m.boot)
  colnames(m.boot)<-mnames
 }
 if(!is.null(dim(covariates.boot)))
 {covariates.boot<-as.matrix(covariates.boot)
  colnames(covariates.boot)<-colnames(covariates)
 }
 if(!is.null(dim(x1.boot)))
 {x1.boot<-as.matrix(x1.boot)
  colnames(x1.boot)<-colnames(data1$x1)
 }
 if(!is.null(dim(m1y.boot)))
 {m1y.boot<-as.matrix(m1y.boot)
  colnames(m1y.boot)<-colnames(data1$m1y)
 }
 if(!is.null(dim(m2y.boot)))
 {m2y.boot<-as.matrix(m2y.boot)
  colnames(m2y.boot)<-colnames(data1$m2y)
 }
 if(!is.null(dim(xm1.boot)))
 {xm1.boot<-as.matrix(xm1.boot)
 colnames(xm1.boot)<-colnames(data1$xm1)
 }}

#b. create xm2 and m.2
 lc2<-c(l2,c2)
 if(!is.null(lc2))                  #create x variables to explain level 2 mediatiors
 {temp<-cbind(x.boot,m.boot[,lc2])
  temp2<-two(temp,level.boot)
  xm2.boot<-as.matrix(temp2[,1])
  colnames(xm2.boot)<-"x.2"
  m.2.boot<-as.matrix(temp2[,-1])
  colnames(m.2.boot)<-mnames[lc2]
  if(!is.null(f01km2))
  {allfun<-f01km2[[2]]
   if (length(f01km2)>2)
    for(i in 3:length(f01km2))
     allfun<-c(allfun,f01km2[[i]])
   unifun<-unique(allfun)
   unifun1<-unifun[unifun!="x"]
   unifun2<-c("x",unifun1)
   d_d<-x2fx(xm2.boot,unifun1)
   d<-as.matrix(d_d$values)
   temp.xm2<-c(colnames(xm2.boot),paste("x.2",1:ncol(d),sep="."))
   xm2.boot<-cbind(xm2.boot,d)
   colnames(xm2.boot)<-temp.xm2
 }}
 else
 {m.2.boot<-NULL
  xm2.boot<-NULL}
 
 #c. create the new m matrix from fm1 and fm2 by the whole training data set
 if(is.null(m.new1))
 {n<-length(x.new)
 m.new<-matrix(0,n,length(mnames))
 colnames(m.new)<-mnames
 fm1<-full$fm1
 if(!is.null(l1))        #generate level 1 continuous mediators
 {for (i in 1:length(l1))
 {if(sum(cm[[1]]==l1[i])>0)
 {temp2<-(1:length(cm[[1]]))[cm[[1]]==l1[i]]+1
 temp.cov<-as.matrix(covariates.new[,cm[[temp2]]])
 colnames(temp.cov)<-colnames(covariates.new)[cm[[temp2]]]}
   else
     temp.cov<-NULL
   k2<-(1:length(data1$fm12[[1]]))[data1$fm12[[1]]==l1[i]]
   # k1<-NULL                        #is there level 1 predictor?
   if(levelx==1)
   {k1<-(1:length(data1$fm11[[1]]))[data1$fm11[[1]]==l1[i]]
   expl.m<-as.matrix(xm1[,c(data1$fm12[[k2+1]],data1$fm11[[k1+1]])])
   colnames(expl.m)<-colnames(xm1)[c(data1$fm12[[k2+1]],data1$fm11[[k1+1]])]}
   else 
   {expl.m<-as.matrix(xm1[,data1$fm12[[k2+1]]])
   colnames(expl.m)<-colnames(xm1)[data1$fm12[[k2+1]]]}
   expl.m<-cbind(expl.m,temp.cov)
   temp.data<-cbind(level=level,expl.m)
   colnames(temp.data)<-c("level",colnames(expl.m))
   j<-(1:length(fm1[[1]]))[fm1[[1]]==l1[i]]+1
   sd.m<-rnorm(n,mean=0,sd=as.data.frame(VarCorr(fm1[[j]]))[2,5])+         #add the randomn effects
     one(rnorm(n2,mean=0,sd=as.data.frame(VarCorr(fm1[[j]]))[1,5]),level)
   m.new[,l1[i]]<-predict(fm1[[j]],newdata=data.frame(temp.data),allow.new.levels=T)+sd.m
 }
 }
 if(!is.null(c1))      #generate level 1 categorical mediators
   for (i in 1:length(c1))
   {if(sum(cm[[1]]==c1[i])>0)
   {temp2<-(1:length(cm[[1]]))[cm[[1]]==c1[i]]+1
   temp.cov<-as.matrix(covariates.new[,cm[[temp2]]])
   colnames(temp.cov)<-colnames(covariates.new)[cm[[temp2]]]}
     else
       temp.cov<-NULL
     k1<-NULL                        #is there level 1 predictor?
     if(levelx==1)
       k1<-(1:length(data1$fm11[[1]]))[data1$fm11[[1]]==c1[i]]
     k2<-(1:length(data1$fm12[[1]]))[data1$fm12[[1]]==c1[i]]
     expl.m<-as.matrix(xm1[,c(data1$fm12[[k2+1]],data1$fm11[[k1+1]])])
     colnames(expl.m)<-colnames(xm1)[c(data1$fm12[[k2+1]],data1$fm11[[k1+1]])]
     expl.m<-cbind(expl.m,temp.cov)
     temp.data<-cbind(level=level,expl.m)
     colnames(temp.data)<-c("level",colnames(expl.m))
     j<-(1:length(fm1[[1]]))[fm1[[1]]==c1[i]]+1
     if(length(j)==1)   #binary
     {tmp.logit<-predict(fm1[[j]],newdata=data.frame(temp.data),allow.new.levels=T)+
       one(rnorm(n2,mean=0,sd=as.data.frame(VarCorr(fm1[[j]]))[,5]),level)  #add the random effect
     tmp.p<-exp(tmp.logit)/(exp(tmp.logit)+1)
     tmp.m<-rbinom(n,1,prob=tmp.p)
     m.new[,c1[i]]<-ifelse(tmp.m==0,1,2)}
     else               #multi-categorical
     {tmp.logit<-rep(1,n)
     for (k in 1:length(j))
       tmp.logit<-cbind(tmp.logit,exp(predict(fm1[[j[k]]],newdata=data.frame(temp.data),allow.new.levels=T)+
                                        one(rnorm(n2,mean=0,sd=as.data.frame(VarCorr(fm1[[j[k]]]))[,5]),level))) #add the random effect
     tmp.m<-apply(tmp.logit,1,rmultinom,n=1,size=1)
     m.new[,c1[i]]<-(1:nrow(tmp.m))%*%tmp.m   #the reference group is always 1
     }
   }
 
 fm2<-full$fm2            #generate level 2 mediators
 cov.2<-NULL
 if(!is.null(covariates.new))
   cov.2<-two(covariates.new,level)
 if(!is.null(l2))          #generate level 2 mediators
   for (i in 1:length(l2))
   {if(sum(cm[[1]]==l2[i])>0)
   {temp2<-(1:length(cm[[1]]))[cm[[1]]==l2[i]]+1
   temp.cov<-as.matrix(cov.2[,cm[[temp2]]])
   colnames(temp.cov)<-colnames(covariates.new)[cm[[temp2]]]}
     else
       temp.cov<-NULL
     expl.m<-as.matrix(xm2[,data1$fm22[[i+1]]])
     colnames(expl.m)<-colnames(xm2)[data1$fm22[[i+1]]]
     expl.m<-cbind(expl.m,temp.cov)
     j<-(1:length(fm2[[1]]))[fm2[[1]]==l2[i]]+1
     tmp.m<-predict(fm2[[j]],newdata=data.frame(expl.m))+
       rnorm(n2,mean=0,sd=(summary(fm2[[j]]))$sigma)  #add the randomness
     m.new[,l2[i]]<-one(tmp.m,level) #lm was used
   }
 if(!is.null(c2))    #did not test
   for (i in 1:length(c2))
   {if(sum(cm[[1]]==c2[i])>0)
   {temp2<-(1:length(cm[[1]]))[cm[[1]]==c2[i]]+1
   temp.cov<-as.matrix(cov.2[,cm[[temp2]]])
   colnames(temp.cov)<-colnames(covariates.new)[cm[[temp2]]]}
     else
       temp.cov<-NULL
     k<-(1:length(data1$fm22[[1]]))[data1$fm22[[1]]==c2[i]]  ####
     expl.m<-as.matrix(xm2[,data1$fm22[[k+1]]])
     colnames(expl.m)<-colnames(xm2)[data1$fm22[[k+1]]]
     expl.m<-cbind(expl.m,temp.cov)
     j<-(1:length(fm2[[1]]))[fm2[[1]]==c2[i]]+1
     if(length(j)==1) #binary
     {tmp.m<-rbinom(n2,1,prob=predict(fm2[[j]],newdata=data.frame(expl.m),type="response"))
     m.new[,c2[i]]<-one(ifelse(tmp.m==0,1,2),level)}
     else               #multi-categorical
     {tmp.logit<-rep(1,n2)
     for (k in 1:length(j))
       tmp.logit<-cbind(tmp.logit,predict(fm2[[j[k]]],newdata=data.frame(expl.m)))
     tmp.m<-apply(tmp.logit,1,rmultinom,n=1,size=1)
     m.new[,c2[i]]<-one((1:nrow(tmp.m))%*%tmp.m,level)   #the reference group is always 1
     }
   }           #glm
 
 #d. calculate m1y.der and m2y.der based on the new x and new m.
 if(!is.null(l2))        #calculate m1y.der and m2y.der
 {m2y.der<-matrix(1,n,length(l2))
 colnames(m2y.der)<-colnames(m)[l2]
 if(!is.null(f02ky)) 
   for (i in 2:length(f02ky))
   {a<-f02ky[[1]][i-1]
   if(sum(l2==a)>0)
   {b<-(1:length(l2))[l2==a]
   d.der<-as.matrix(x2fdx(m.new[,a],f02ky[[i]]))
   colnames(d.der)<-paste(mnames[a],1:ncol(d.der),sep=".")
   if(dim(as.matrix(d.der))[2]==1)
     m2y.der[,b]<-d.der
   else {
     m2y.der[,b]<-d.der[,1]
     temp.m2<-c(colnames(m2y.der),colnames(d.der)[-1])
     m2y.der<-cbind(m2y.der,d.der[,-1])
     colnames(m2y.der)<-temp.m2}
   }}}
 else
   m2y.der<-NULL
 if(!is.null(c2))        #binarize level 2 categorical variables to exaplain y
 {temp<-cattobin(m1y=NULL,m2y.der,m1=NULL,m.new,c2,c2r=rep(1,length(c2)))
  m2y.der<-temp$m1y.der
 }
 if(!is.null(l1))         #create level 1 and corresponding level 2 m variables to explain y
 {temp<-apply(as.matrix(m.new[,l1]), 2, one2two, level, weight)
 m1y.der<-matrix(1,n,length(l1))
 colnames(m1y.der)<-mnames[l1]
 m12y.der<-matrix(1,n,ncol(temp))
 colnames(m12y.der)<-paste(mnames[l1],"12",sep=".")
 if(!is.null(f02ky)) 
   for (i in 2:length(f02ky))
   {a<-f02ky[[1]][i-1]
   if(sum(l1==a)>0)
   {b<-(1:length(l1))[l1==a]
   tt<-one2two(m.new[,a],level,weight)
   d.der<-as.matrix(x2fdx(tt[,1],f02ky[[i]]))
   colnames(d.der)<-paste(mnames[a],1:ncol(d.der),sep=".")
   if(dim(as.matrix(d))[2]==1)
     m12y.der[,b]<-d.der
   else {
     m12y.der[,b]<-d.der[,1]
     temp.m12<-c(colnames(m12y.der),colnames(d.der)[-1])
     m12y.der<-cbind(m12y.der,d.der[,-1])
     colnames(m12y.der)<-temp.m12}
   }}
 # added in case the transformation f02ky is on the aggregated level 1 mediator
 if(!is.null(f20ky)) 
   for (i in 2:length(f20ky))
   {a<-f20ky[[1]][i-1]
   b<-(1:length(l1))[l1==a]
   tt<-one2two(m.new[,a],level,weight)
   d.der<-as.matrix(x2fdx(tt[,2],f20ky[[i]]))
   colnames(d.der)<-paste(mnames[a],1:ncol(d.der),sep=".")
   if(ncol(d.der)==1)
     m1y.der[,b]<-d.der
   else {
     temp.namem1y<-c(colnames(m1y.der),colnames(d.der)[-1])
     m1y.der[,b]<-d.der[,1]
     m1y.der<-cbind(m1y.der,d.der[,-1])
     colnames(m1y.der)<-temp.namem1y
   }
   }
 m2y.der<-cbind(m2y.der,m12y.der)
 }
 else 
   m1y.der<-NULL
 if(!is.null(c1))                  #binarize level 1 categorical variables and corresponding level 2 variables to explain y 
 {temp<-cattobin(m1y=NULL,m1y.der=m1y.der,m1=NULL,m=m.new,cat1=c1,cat2=rep(1,length(c1)),level,weight)
  m1y.der<-temp$m1y.der
  m2y.der<-cbind(m2y.der,temp$m12y.der)
 }
}
 
# 3. calculate the mediation effects on x.new based on f1, fm1, and fm2 built on the bootstrap sample 
 
 
# data1<-data.org(x.boot, levelx, levely, m.boot, l1, l2, c1, c1r, c2, c2r, f01y, f10y, 
#                 f02ky, f20ky, f01km1, f01km2, f10km, level.boot, weight.boot)
# temp<-mlma(y.boot, biny, data1, x.boot, levelx, levely, m.boot, l1,l2, c1, c1r, 
#            c2, c2r, level.boot, weight.boot, random, random.m1,intercept, covariates, cy1, cy2, cm,joint) 
# a<-list(de1=DE1,de2=two(DE2,level),ie1=ie1,ie2=ie2,ie12=ie12,
#         f1=f1,fm1=fm1,fm2=fm2, ie1_list=ie1_list, ie2_list=ie2_list,
#         iej2_list=iej_list,ie12_1=ie12_1, ie12_2=ie12_2, ie1_1=ie1_1,
#         ie1_2=ie1_2,ie2_1=ie2_1,ie2_2=ie2_2,x=x, x.j=two(x,level), 
#         m=m, covariates=covariates, intercept=intercept, cm=cm,
#         data1=data1)

 if(!is.null(data1$l1x))  #l1x could be NULL when x is level 2
 {tt2<-as.matrix(x1.boot[,data1$l1x])
 colnames(tt2)<-colnames(x1.boot)[data1$l1x]}
 else
   tt2<-NULL             
 tt1<-as.matrix(x1.boot[,data1$l2x])
 colnames(tt1)<-colnames(x1.boot)[data1$l2x]
 if(!is.null(c(cy1,cy2)))           ###added covariates to explain y
 {cova<-as.matrix(covariates.boot[,c(cy1,cy2)])
 colnames(cova)<-colnames(covariates)[c(cy1,cy2)]}
 else
   cova<-NULL
 expl<-cbind(tt1,m2y.boot,tt2,m1y.boot,cova) 
 
 temp.data<-cbind(y=y.boot,level=level.boot,x1.boot,m1y.boot,m2y.boot, covariates.boot)
 colnames(temp.data)<-c("y","level",colnames(data1$x1),colnames(data1$m1y),
                        colnames(data1$m2y), colnames(covariates))
 if(levely==1)
 {frml<-getformula(expl,random,intercept)
 if(biny)
   f1<-glmer(frml,data=data.frame(temp.data),family=binomial(link="logit"))
 else
   f1<-lmer(frml,data=data.frame(temp.data))}
 else
 {temp.data2<-two(temp.data, level.boot, weight.boot)
 frml<-getformula(expl,random=NULL,intercept)
 if(biny)
   f1<-glm(frml,data=data.frame(temp.data2),family=binomial(link="logit"))
 else
   f1<-lm(frml,data=data.frame(temp.data2))}
 
 lc1<-c(l1,c1)
 lc2<-c(l2,c2)
 
 if(intercept)
   coef.f1<-summary(f1)$coefficient[-1,1] 
 else 
   coef.f1<-summary(f1)$coefficient[,1]
 len<-c(length(data1$l2x),ifelse(is.null(m2y.boot),0,ncol(m2y.boot)),
        ifelse(is.null(data1$l1x),0,length(data1$l1x)), 
        ifelse(is.null(m1y.boot),0,ncol(m1y.boot)))
 if (len[1]>1)                       #calculate level 2 direct effect
   DE2<-x1.der[,data1$l2x]%*%coef.f1[1:length(data1$l2x)]
 else
   DE2<-x1.der[,data1$l2x]*coef.f1[1]
 if (len[2]==0)
   ie2_1<-NULL
 else 
 {ie2_1<-NULL
 ie12_1<-NULL
 if (!is.null(data1$m2))
 {if(!is.null(l2))
   for (k in 1:length(l2))
   {if (length(data1$m2[[k+1]])==1)
     ie2_1<-cbind(ie2_1,m2y.der[,data1$m2[[k+1]]]*coef.f1[len[1]+data1$m2[[k+1]]])
   else
     ie2_1<-cbind(ie2_1,m2y.der[,data1$m2[[k+1]]]%*%coef.f1[len[1]+data1$m2[[k+1]]])}
   
   if(!is.null(c2))
     for (i in 1:length(c2)) 
     {temp<-(1:length(data1$m2[[1]]))[data1$m2[[1]]==c2[i]]
     if (length(data1$m2[[temp + 1]])==1)
       ie2_1 <- cbind(ie2_1, m2y.der[, data1$m2[[temp + 1]]] * 
                        coef.f1[len[1] + data1$m2[[temp + 1]]])
     else
       for (k in 1:length(data1$m2[[temp + 1]])) 
         ie2_1 <- cbind(ie2_1, m2y.der[, data1$m2[[temp + 1]]][, k] * 
                          coef.f1[len[1] + data1$m2[[temp + 1]][k]])
     }  }        
 if (!is.null(data1$m12))
 {if(!is.null(l1))
   for (k in 1:length(l1))
   {if (length(data1$m12[[k+1]])==1)
     ie12_1<-cbind(ie12_1,m2y.der[,data1$m12[[k+1]]]*coef.f1[len[1]+data1$m12[[k+1]]])
   else
     ie12_1<-cbind(ie12_1,m2y.der[,data1$m12[[k+1]]]%*%coef.f1[len[1]+data1$m12[[k+1]]])}
   
   if(!is.null(c1))
     for (i in 1:length(c1)) 
     {temp<-(1:length(data1$m12[[1]]))[data1$m12[[1]]==c1[i]]
     if(length(data1$m12[[temp + 1]])==1)
       ie12_1 <- cbind(ie12_1,m2y.der[, data1$m12[[temp + 1]]] * coef.f1[len[1] + data1$m12[[temp + 1]]])
     else
       for (k in 1:length(data1$m12[[temp + 1]])) 
         ie12_1 <- cbind(ie12_1,m2y.der[, data1$m12[[temp + 1]]][,k] * coef.f1[len[1] + data1$m12[[temp + 1]][k]])
     }          
 }
 }
 if (len[3]==0)
   DE1<-NULL
 else if (len[3]==1)
   DE1<-x1.der[,data1$l1x]*coef.f1[len[1]+len[2]+1]
 else
   DE1<-x1.der[,data1$l1x]%*%coef.f1[(len[1]+len[2]+1):(len[1]+len[2]+len[3])]
 if (len[4]==0)
   ie1_1<-NULL
 else 
 {ie1_1<-NULL
 if (!is.null(data1$m1))
 {if (!is.null(l1))
   for (k in 1:length(l1))
   {if(length(data1$m1[[k+1]])==1)
     ie1_1<-cbind(ie1_1,m1y.der[,data1$m1[[k+1]]]*coef.f1[len[1]+len[2]+len[3]+data1$m1[[k+1]]])
   else
     ie1_1<-cbind(ie1_1,m1y.der[,data1$m1[[k+1]]]%*%coef.f1[len[1]+len[2]+len[3]+data1$m1[[k+1]]])}
   
   if(!is.null(c1))
     for (i in 1:length(c1)) 
     {temp<-(1:length(data1$m1[[1]]))[data1$m1[[1]]==c1[i]]
     if (length(data1$m1[[temp + 1]])==1)
       ie1_1 <- cbind(ie1_1,m1y.der[, data1$m1[[temp + 1]]] * 
                        coef.f1[len[1] + len[2] + len[3] + data1$m1[[temp + 1]]])
     else
       for (k in 1:length(data1$m1[[temp + 1]])) 
         ie1_1 <- cbind(ie1_1,m1y.der[, data1$m1[[temp + 1]]][, k] * 
                          coef.f1[len[1] + len[2] + len[3] + data1$m1[[temp + 1]][k]])
     }  } }
 

if(levelx==2)    #analysis when x is a level 2 variable
 {ie12_2<-NULL
 if(!is.null(l1))
   for (i in 1:length(l1))
   {if(sum(cm[[1]]==l1[i])>0)
   {temp2<-(1:length(cm[[1]]))[cm[[1]]==l1[i]]+1
   temp.cov<-as.matrix(covariates.boot[,cm[[temp2]]])
   colnames(temp.cov)<-colnames(covariates)[cm[[temp2]]]}
     else
       temp.cov<-NULL
     expl.m<-as.matrix(xm1.boot[,data1$fm12[[i+1]]])
     colnames(expl.m)<-colnames(data1$xm1)[data1$fm12[[i+1]]]
     numx<-length(data1$fm12[[i+1]]) #1:
     expl.m<-cbind(expl.m,temp.cov)
     frml.m<-getformula(expl.m,random,intercept)
     temp.data<-cbind(y=m.boot[,l1[i]],level=level.boot,expl.m)
     colnames(temp.data)<-c("y","level",colnames(expl.m))
     model<-lmer(frml.m,data=data.frame(temp.data))
     coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
     if(intercept)
       coef.temp<-coef.temp[-1]
     coef.temp<-coef.temp[1:numx]
     if(length(coef.temp)==1)
       ie12_2<-cbind(ie12_2,coef.temp*xm1.der[,data1$fm12[[i+1]]])
     else
       ie12_2<-cbind(ie12_2,xm1.der[,data1$fm12[[i+1]]]%*%coef.temp)
   }
 
 if(!is.null(c1))
   for (i in 1:length(c1))
   {if(sum(cm[[1]]==c1[i])>0)
   {temp2<-(1:length(cm[[1]]))[cm[[1]]==c1[i]]+1
   temp.cov<-as.matrix(covariates.boot[,cm[[temp2]]])
   colnames(temp.cov)<-colnames(covariates)[cm[[temp2]]]
   temp.cov.new<-as.matrix(covariates.new[,cm[[temp2]]])
   colnames(temp.cov.new)<-colnames(covariates)[cm[[temp2]]]}
     else
     {temp.cov<-NULL
     temp.cov.new<-NULL}
     k2<-(1:length(data1$fm12[[1]]))[data1$fm12[[1]]==c1[i]]
     expl.m<-as.matrix(xm1.boot[,data1$fm12[[k2+1]]])
     colnames(expl.m)<-colnames(data1$xm1)[data1$fm12[[k2+1]]]
     expl.m.new<-as.matrix(xm1[,data1$fm12[[k2+1]]])
     colnames(expl.m.new)<-colnames(data1$xm1)[data1$fm12[[k2+1]]]
     numx<-length(data1$fm12[[k2+1]]) #1:
     expl.m<-cbind(expl.m,temp.cov)
     expl.m.new<-cbind(expl.m.new,temp.cov.new)
     frml.m<-getformula(expl.m,random,intercept)
     temp<-(1:length(data1$m1[[1]]))[data1$m1[[1]]==c1[i]]
     if(length(data1$m1[[temp+1]])>1)######
     {temp.4<-apply(m1y.boot[,data1$m1[[temp+1]]],1,sum)
     temp.5<-ifelse(temp.4==1,F,T)}
     else 
       temp.5<-rep(T,n)
     for (k in data1$m1[[temp+1]])
     {temp.6<-(temp.5 | (m1y.boot[,k]==1))
     temp.data<-cbind(y=m1y.boot[temp.6,k],level=level.boot[temp.6],expl.m[temp.6,]) #temp.6
     colnames(temp.data)<-c("y","level",colnames(expl.m))
     model<-glmer(frml.m,data=data.frame(temp.data),
                  family=binomial(link="logit"))
     temp.data<-cbind(level=level.new,expl.m.new) #temp.6
     colnames(temp.data)<-c("level",colnames(expl.m))
     p.temp<-predict(model,type="response",newdata=data.frame(temp.data),allow.new.levels=T)  ####################
     coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
     if(intercept)
       coef.temp<-coef.temp[-1]
     coef.temp<-coef.temp[1:numx]
     if(length(coef.temp)==1)
       ie12_2<-cbind(ie12_2,p.temp*(1-p.temp)*
                       (coef.temp*xm1.der[,data1$fm12[[k2+1]]]))
     else
       ie12_2<-cbind(ie12_2,p.temp*(1-p.temp)*
                       (xm1.der[,data1$fm12[[k2+1]]]%*%coef.temp))
     }
   }
 
 #models for x to explain level 2 mediators
 ie2_2<-NULL
 if(!is.null(covariates))
 {cov.2<-two(covariates.boot,level.boot,weight.boot)
 cov.2.new<-two(covariates.new,level.new,weight.new)}
 else
 {cov.2<-NULL
 cov.2.new<-NULL}
 if(!is.null(l2))
   for (i in 1:length(l2))
   {if(sum(cm[[1]]==l2[i])>0)
   {temp2<-(1:length(cm[[1]]))[cm[[1]]==l2[i]]+1
   temp.cov<-as.matrix(cov.2[,cm[[temp2]]])
   colnames(temp.cov)<-colnames(covariates)[cm[[temp2]]]}
     else
       temp.cov<-NULL
     expl.m<-as.matrix(xm2.boot[,data1$fm22[[i+1]]])
     colnames(expl.m)<-colnames(data1$xm2)[data1$fm22[[i+1]]]
     numx<-ncol(expl.m)
     expl.m<-cbind(expl.m,temp.cov)
     frml.m<-getformula(expl.m,random=NULL,intercept)
     temp.data<-cbind(y=m.2.boot[,i],expl.m)
     colnames(temp.data)<-c("y",colnames(expl.m))
     model<-lm(frml.m,data=data.frame(temp.data))
     coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
     if(intercept)
       coef.temp<-coef.temp[-1]
     coef.temp<-coef.temp[1:numx]
     if(length(coef.temp)==1)
       ie2_2<-cbind(ie2_2,coef.temp*xm2.der[,data1$fm22[[i+1]]])
     else
       ie2_2<-cbind(ie2_2,xm2.der[,data1$fm22[[i+1]]]%*%coef.temp)
   }
 
 if(!is.null(c2))    #did not test
   for (i in 1:length(c2))
   {if(sum(cm[[1]]==c2[i])>0)
   {temp2<-(1:length(cm[[1]]))[cm[[1]]==c2[i]]+1
   temp.cov<-as.matrix(cov.2[,cm[[temp2]]])
   colnames(temp.cov)<-colnames(covariates)[cm[[temp2]]]
   temp.cov.new<-as.matrix(cov.2.new[,cm[[temp2]]])
   colnames(temp.cov.new)<-colnames(covariates)[cm[[temp2]]]}
     else
     {temp.cov<-NULL
     temp.cov.new<-NULL}
     k2<-(1:length(data1$fm22[[1]]))[data1$fm22[[1]]==c2[i]]
     expl.m<-as.matrix(xm2.boot[,data1$fm22[[k2+1]]])
     colnames(expl.m)<-colnames(data1$xm2)[data1$fm22[[k2+1]]]
     expl.m.new<-as.matrix(xm2[,data1$fm22[[k2+1]]])
     colnames(expl.m.new)<-colnames(data1$xm2)[data1$fm22[[k2+1]]]
     numx<-ncol(expl.m)
     expl.m<-cbind(expl.m,temp.cov)
     expl.m.new<-cbind(expl.m.new,temp.cov.new)
     frml.m<-getformula(expl.m,random=NULL,intercept)
     temp<-(1:length(data1$m2[[1]]))[data1$m2[[1]]==c2[i]]
     if(length(data1$m2[[temp+1]])>1)  ####
     {temp.3<-apply(m2y.boot[,data1$m2[[temp+1]]],2,two,level)
     temp.4<-apply(temp.3,1,sum)
     temp.5<-(temp.4==0)}
     else temp.5<-rep(T,length(unique(level)))
     for (k in data1$m2[[temp+1]])
     {temp.6<-(temp.5 | (two(m2y.boot[,k],level)==1))
     temp.data<-cbind(y=two(m2y.boot[temp.6,k],level=level.boot[temp.6]),expl.m[temp.6,]) #
     colnames(temp.data)<-c("y",colnames(expl.m)) #"level",
     model<-glm(frml.m,data=data.frame(temp.data),
                family=binomial(link="logit"))
     temp.data<-expl.m.new
     #colnames(temp.data)<-colnames(data1$xm2)[data1$fm22[[j+1]]]
     p.temp<-predict(model,type="response",newdata=data.frame(temp.data))
     coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
     if(intercept)
       coef.temp<-coef.temp[-1]
     coef.temp<-coef.temp[1:numx]
     if(length(coef.temp)==1)
       ie2_2<-cbind(ie2_2,p.temp*(1-p.temp)*coef.temp*xm2.der[,data1$fm22[[i+1]]])
     else
       ie2_2<-cbind(ie2_2,p.temp*(1-p.temp)*
                      (xm2.der[,data1$fm22[[i+1]]]%*%coef.temp))
     }
   }
 ie1_2<-NULL}
 else   #analysis when x is a level 1 variable
 {ie1_2<-NULL
 ie12_2<-NULL
 if(!is.null(l1))
   for (i in 1:length(l1))
   {if(sum(cm[[1]]==l1[i])>0) #add covariates
   {temp2<-(1:length(cm[[1]]))[cm[[1]]==l1[i]]+1
   temp.cov<-as.matrix(covariates.boot[,cm[[temp2]]])
   colnames(temp.cov)<-colnames(covariates)[cm[[temp2]]]
   temp.cov.new<-as.matrix(covariates.new[,cm[[temp2]]])
   colnames(temp.cov.new)<-colnames(covariates)[cm[[temp2]]]}
     else
     {temp.cov<-NULL
     temp.cov.new<-NULL}
     expl.m<-as.matrix(xm1.boot[,c(data1$fm12[[i+1]],data1$fm11[[i+1]])])
     colnames(expl.m)<-colnames(data1$xm1)[c(data1$fm12[[i+1]],data1$fm11[[i+1]])]
     expl.m.new<-as.matrix(xm1[,c(data1$fm12[[i+1]],data1$fm11[[i+1]])])
     colnames(expl.m.new)<-colnames(data1$xm1)[c(data1$fm12[[i+1]],data1$fm11[[i+1]])]
     numx<-ncol(expl.m)
     expl.m<-cbind(expl.m,temp.cov)
     expl.m.new<-cbind(expl.m.new,temp.cov.new)
     m.random<-"(1|level)"
     if(!is.null(random.m1))
       if(sum(random.m1[[1]]==l1[i])>0)
         m.random<-random.m1[[2]][random.m1[[1]]==l1[i]] #1st item of random.m1 is the list of l1 med, 2nd item is the random item of the same order
     frml.m<-getformula(expl.m,m.random,intercept)
     temp.data<-cbind(y=m.boot[,l1[i]],level=level.boot,expl.m)
     colnames(temp.data)<-c("y","level",colnames(expl.m))
     model=lmer(frml.m,data=data.frame(temp.data))
     coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
     if(intercept)
       coef.temp<-coef.temp[-1]
     coef.temp<-coef.temp[1:numx]
     coef.temp12<-coef.temp[1:length(data1$fm12[[i+1]])]
     coef.temp1<-coef.temp[-(1:length(data1$fm12[[i+1]]))]
     if(length(coef.temp12)==1)
       ie12_2<-cbind(ie12_2,coef.temp12*xm1.der[,data1$fm12[[i+1]]])
     else
       ie12_2<-cbind(ie12_2,xm1.der[,data1$fm12[[i+1]]]%*%coef.temp12)
     if(length(coef.temp1)==1)
       ie1_2<-cbind(ie1_2,coef.temp1*xm1.der[,data1$fm11[[i+1]]])
     else
       ie1_2<-cbind(ie1_2,xm1.der[,data1$fm11[[i+1]]]%*%coef.temp1)
   }
 if(!is.null(c1))  
   for (i in 1:length(c1))
   {if(sum(cm[[1]]==c1[i])>0)   #add covariates
   {temp2<-(1:length(cm[[1]]))[cm[[1]]==c1[i]]+1
   temp.cov<-as.matrix(covariates.boot[,cm[[temp2]]])
   colnames(temp.cov)<-colnames(covariates)[cm[[temp2]]]    
   temp.cov.new<-as.matrix(covariates.new[,cm[[temp2]]])
   colnames(temp.cov.new)<-colnames(covariates)[cm[[temp2]]]}
     else
     {temp.cov<-NULL
     temp.cov.new<-NULL}
     k1<-(1:length(data1$fm12[[1]]))[data1$fm12[[1]]==c1[i]]
     k2<-(1:length(data1$fm11[[1]]))[data1$fm11[[1]]==c1[i]]     
     expl.m<-as.matrix(xm1.boot[,c(data1$fm12[[k1+1]],data1$fm11[[k2+1]])])
     colnames(expl.m)<-colnames(data1$xm1)[c(data1$fm12[[k1+1]],data1$fm11[[k2+1]])]
     expl.m.new<-as.matrix(xm1[,c(data1$fm12[[k1+1]],data1$fm11[[k2+1]])])
     colnames(expl.m.new)<-colnames(data1$xm1)[c(data1$fm12[[k1+1]],data1$fm11[[k2+1]])]
     numx<-ncol(expl.m)
     expl.m<-cbind(expl.m,temp.cov)
     expl.m.new<-cbind(expl.m.new,temp.cov.new)
     m.random<-"(1|level)"
     if(!is.null(random.m1))
       if(sum(random.m1[[1]]==c1[i])>0)
         m.random<-random.m1[[2]][random.m1[[1]]==c1[i]] #1st item of random.m1 is the list of l1 med, 2nd item is the random item of the same order
     frml.m<-getformula(expl.m,m.random,intercept)
     temp<-(1:length(data1$m1[[1]]))[data1$m1[[1]]==c1[i]]
     if(length(data1$m1[[temp+1]])>1)###
     {temp.4<-apply(m1y.boot[,data1$m1[[temp+1]]],1,sum)
     temp.5<-ifelse(temp.4==1,F,T)}
     else temp.5<-rep(T,n)
     for (k in data1$m1[[temp+1]])
     {temp.6<-(temp.5 | (m1y.boot[,k]==1))
     temp.data<-cbind(y=m1y.boot[temp.6,k],level=level.boot[temp.6],expl.m[temp.6,])  #[temp.6]
     colnames(temp.data)<-c("y","level", colnames(expl.m))
     model=glmer(frml.m,data=data.frame(temp.data),
                 family=binomial(link="logit"))
     temp.data<-cbind(level=level.new,expl.m.new)  #[temp.6]
     colnames(temp.data)<-c("level",colnames(expl.m))
     p.temp<-predict(model,type="response",newdata=data.frame(temp.data),allow.new.levels=T)
     coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
     if(intercept)
       coef.temp<-coef.temp[-1]
     coef.temp<-coef.temp[1:numx]
     coef.temp12<-coef.temp[1:length(data1$fm12[[k1+1]])]
     coef.temp1<-coef.temp[-(1:length(data1$fm12[[k1+1]]))]
     if(length(coef.temp12)==1)
       ie12_2<-cbind(ie12_2,coef.temp12*xm1.der[,data1$fm12[[k1+1]]])
     else
       ie12_2<-cbind(ie12_2,xm1.der[,data1$fm12[[k1+1]]]%*%coef.temp12)
     if(length(coef.temp1)==1)
       ie1_2<-cbind(ie1_2,p.temp*(1-p.temp)*coef.temp1*xm1.der[,data1$fm11[[k2+1]]])
     else
       ie1_2<-cbind(ie1_2,p.temp*(1-p.temp)*(xm1.der[,data1$fm11[[k2+1]]]%*%coef.temp1))
     }
   }
 
 ie2_2<-NULL   #models for x to explain level 2 mediators
 if(!is.null(covariates))
 {cov.2<-two(covariates.boot,level.boot, weight.boot)
 cov.2.new<-two(covariates.new,level.new, weight.new)}
 else 
 {cov.2<-NULL
 cov.2.new<-NULL}
 if(!is.null(l2))
   for (i in 1:length(l2))
   {if(sum(cm[[1]]==l2[i])>0)   #add covariates
   {temp2<-(1:length(cm[[1]]))[cm[[1]]==l2[i]]+1
   temp.cov<-as.matrix(cov.2[,cm[[temp2]]])
   colnames(temp.cov)<-colnames(covariates)[cm[[temp2]]]}
     else
       temp.cov<-NULL
     expl.m<-as.matrix(xm2.boot[,data1$fm22[[i+1]]])
     colnames(expl.m)<-colnames(data1$xm2)[data1$fm22[[i+1]]]
     numx<-ncol(expl.m)
     expl.m<-cbind(expl.m, temp.cov)
     frml.m<-getformula(expl.m,random=NULL,intercept)
     temp.data<-cbind(y=m.2.boot[,i],expl.m)
     colnames(temp.data)<-c("y",colnames(expl.m))
     model<-lm(frml.m,data=data.frame(temp.data))
     coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
     if(intercept)
       coef.temp<-coef.temp[-1]
     coef.temp<-coef.temp[1:numx]
     if(length(coef.temp)==1)
       ie2_2<-cbind(ie2_2,coef.temp*xm2.der[,data1$fm22[[i+1]]])
     else
       ie2_2<-cbind(ie2_2,xm2.der[,data1$fm22[[i+1]]]%*%coef.temp)
   }
 
 if(!is.null(c2))    #did not test
   for (i in 1:length(c2))
   {if(sum(cm[[1]]==c2[i])>0)   #add covariates
   {temp2<-(1:length(cm[[1]]))[cm[[1]]==c2[i]]+1
   temp.cov<-as.matrix(cov.2[,cm[[temp2]]])
   colnames(temp.cov)<-colnames(covariates)[cm[[temp2]]]
   temp.cov.new<-as.matrix(cov.2.new[,cm[[temp2]]])
   colnames(temp.cov.new)<-colnames(covariates)[cm[[temp2]]]}
     else
     {temp.cov<-NULL
     temp.cov.new<-NULL}
     k2<-(1:length(data1$fm22[[1]]))[data1$fm22[[1]]==c2[i]]
     expl.m<-as.matrix(xm2.boot[,data1$fm22[[k2+1]]])
     colnames(expl.m)<-colnames(data1$xm2)[data1$fm22[[k2+1]]]
     expl.m.new<-as.matrix(xm2[,data1$fm22[[k2+1]]])
     colnames(expl.m.new)<-colnames(data1$xm2)[data1$fm22[[k2+1]]]
     numx<-ncol(expl.m)
     expl.m<-cbind(expl.m, temp.cov)
     expl.m.new<-cbind(expl.m.new, temp.cov.new)
     frml.m<-getformula(expl.m,random=NULL,intercept)
     temp<-(1:length(data1$m2[[1]]))[data1$m2[[1]]==c2[i]]
     if(length(data1$m2[[temp+1]])>1)  ####
     {temp.3<-apply(m2y.boot[,data1$m2[[temp+1]]],2,two,level)
     temp.4<-apply(temp.3,1,sum)
     temp.5<-(temp.4==0)}
     else temp.5<-rep(T,length(unique(level)))
     for (k in data1$m2[[temp+1]])
     {temp.6<-(temp.5 | (two(m2y.boot[,k],level)==1))
     temp.data<-cbind(y=two(m2y.boot[temp.6,k],level.boot[temp.6]),expl.m[temp.6,])
     colnames(temp.data)<-c("y",colnames(expl.m)) #"level",
     model<-glm(frml.m,data=data.frame(temp.data),
                family=binomial(link="logit"))
     temp.data<-expl.m.new
     p.temp<-predict(model,type="response",newdata=data.frame(temp.data))
     coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
     if(intercept)
       coef.temp<-coef.temp[-1]
     coef.temp<-coef.temp[1:numx]
     if(length(coef.temp)==1)
       ie2_2<-cbind(ie2_2,p.temp*(1-p.temp)*coef.temp*xm2.der[,data1$fm22[[i+1]]])
     else
       ie2_2<-cbind(ie2_2,p.temp*(1-p.temp)*
                      (xm2.der[,data1$fm22[[i+1]]]%*%coef.temp))
     }
   }
 }
 
 if(!is.null(ie12_1) & !is.null(ie12_2))
 {ie12_1<-apply(as.matrix(ie12_1),2,two,level)
 ie12_2<-apply(as.matrix(ie12_2),2,two,level)
 ie12<-as.matrix(ie12_1*ie12_2)
 }
 else ie12<-NULL
 if(!is.null(ie1_1) & !is.null(ie1_2))
   ie1<-ie1_1*ie1_2
 else ie1<-NULL
 if(!is.null(ie2_1) & !is.null(ie2_2))
 {ie2_1<-two(ie2_1,level)
 ie2<-ie2_1*ie2_2
 }
 else ie2<-NULL
 
# 4. organize the results
 de1.boot<-cbind(de1.boot,DE1)
 de2.boot<-cbind(de2.boot,two(DE2,level.boot))
 ie1.boot<-rbind(ie1.boot,ie1)
 if(!is.null(ie2))
   ie2.boot<-rbind(ie2.boot,as.matrix(ie2))
 if(!is.null(full$ie12))
   ie12.boot<-rbind(ie12.boot,ie12)

 sum.boot<-sum.mlma(object=full,levelx,de2=two(DE2,level.boot),ie2,de1=DE1,ie1,ie12,level.boot, weight.boot,w2)
 sum.boot1<-rbind(sum.boot1,sum.boot$level1)
 sum.boot2<-rbind(sum.boot2,sum.boot$level2)
 print(l)
}
a<-list(de1=de1.boot,de2=de2.boot,ie1=ie1.boot,ie2=ie2.boot,ie12=ie12.boot, sum.boot1=sum.boot1, 
        sum.boot2=sum.boot2,full=full,levelx=levelx, level=level.new,x=x.new, weight=weight.new, m=m.new) #, xboot=xboot, xjboot=xjboot
class(a)<-"mlma.boot"
return(a)
}

print.mlma<-function(x,...,w2=rep(1,length(object$de2)))
{object<-x
two<-function(x, level, weight=rep(1,nrow(as.matrix(x))))
{x<-as.matrix(x)
levels<-unique(level[!is.na(level)])
x2<-matrix(NA,length(levels),dim(x)[2])
for(i in 1:length(levels))
{cho<-(level==levels[i])
if(sum(cho)>0)
{if(sum(cho)==1)
  x2[i,]<-x[level==levels[i],]
else
{temp<-as.matrix(x[level==levels[i],])
weight1<-weight[level==levels[i]]
x2[i,]<-apply(temp,2,weighted.mean,weight1,na.rm=TRUE)}}}
colnames(x2)<-colnames(x)
x2
}


levelx<-object$data1$parameter$levelx
de2<-weighted.mean(object$de2,w2)
cat("Level 2 Direct Effect:")
print(de2)
if(!is.null(object$ie2))
{ie2<-apply(object$ie2,2,weighted.mean,w2)
 ie2_show<-NULL
 for (i in 1:(length(object$ie2_list)-1))
   ie2_show<-c(ie2_show,sum(ie2[object$ie2_list[[i+1]]]))
 names(ie2_show)<-object$ie2_list[[1]]
 cat("\nLevel 2 indirect effects for level 2 mediators:\n")
 print(ie2_show)}
else ie2=0
if(!is.null(object$ie12))
{ie12<-apply(object$ie12,2,weighted.mean,w2)
 ie12_show<-NULL
 for (i in 1:(length(object$ie1_list)-1))
   ie12_show<-c(ie12_show,sum(ie12[object$ie1_list[[i+1]]]))
 names(ie12_show)<-object$ie1_list[[1]]
 cat("\nLevel 2 indirect effects for aggregated level 1 mediators:\n")
 print(ie12_show)}
else ie12=0

if(!is.null(object$iej2_list)) #added 3/2/2016 to show level 2 joint effects
{ie_all2<-cbind(object$ie12,object$ie2)
 ie_j2<-apply(ie_all2,2,weighted.mean,w2)
 iej2_show<-NULL
 for (i in 1:(length(object$iej2_list)-1))
  iej2_show<-c(iej2_show,sum(ie_j2[object$iej2_list[[i+1]]]))
 names(iej2_show)<-object$iej2_list[[1]]
 cat("\nLevel 2 indirect effects for joint mediators:\n")
 print(iej2_show)}

cat("\nLevel 2 total effect is:")
print(sum(c(de2,ie2,ie12)))
if(levelx==1)
{de1<-two(object$de1,object$data1$parameter$level,object$data1$parameter$weight)
 de1<-weighted.mean(de1,w2)
 cat("\nLevel 1 Direct Effect:")
 print(de1)
 if(!is.null(object$ie1))
 {ie1<-two(object$ie1,object$data1$parameter$level,object$data1$parameter$weight)
  ie1<-apply(ie1,2,weighted.mean,w2)
  ie1_show<-NULL
  for (i in 1:(length(object$ie1_list)-1))
    ie1_show<-c(ie1_show,sum(ie1[object$ie1_list[[i+1]]]))
  names(ie1_show)<-object$ie1_list[[1]]
  cat("\nLevel 1 indirect Effects:\n")
  print(ie1_show)
 } 
 else ie1=0
 cat("\nLevel 1 total effect is:")
 print(sum(c(de1,ie1)))
} 
}

summary.mlma<-function(object,...,type="III")
{f1<-Anova(object$f1,type=type)
 mnames<-colnames(object$m)
 if(!is.null(object$fm1[[1]]))
  {temp<-object$fm1
   temp[[1]]<-NULL
   fm1<-lapply(temp,Anova,type="III")
   names(fm1)<-mnames[object$fm1[[1]]]
  }
 else fm1<-NULL
 if(!is.null(object$fm2[[1]]))
 {temp<-object$fm2
  temp[[1]]<-NULL
  fm2<-lapply(temp,Anova,type="III")
  names(fm2)<-mnames[object$fm2[[1]]]
 }
 else fm2<-NULL
a<-list(f1=f1,fm1=fm1,fm2=fm2)
class(a)<-"summary.mlma"
return(a)
}

print.summary.mlma<-function(x,...)
{cat("1. Anova on the Full Model:\n")
 print(x$f1)
 if(!is.null(x$fm1))
 {cat("\n\n2. Anova on models for Level 1 mediators:\n")
   print(x$fm1)}
 if(!is.null(x$fm2))
 {cat("\n\n3. Anova on models for Level 2 mediators:\n")
   print(x$fm2)}
}


#other functions
plot.mlma<-function(x,..., var=NULL, cate=FALSE, w2=rep(1,length(object$de2)))
{object<-x
two<-function(x, level, weight=rep(1,nrow(as.matrix(x))))
{x<-as.matrix(x)
levels<-unique(level[!is.na(level)])
x2<-matrix(NA,length(levels),dim(x)[2])
for(i in 1:length(levels))
{cho<-(level==levels[i])
if(sum(cho)>0)
{if(sum(cho)==1)
  x2[i,]<-x[level==levels[i],]
else
{temp<-as.matrix(x[level==levels[i],])
weight1<-weight[level==levels[i]]
x2[i,]<-apply(temp,2,weighted.mean,weight1,na.rm=TRUE)}}}
colnames(x2)<-colnames(x)
x2
}

op <- par(no.readonly = TRUE) # the whole list of settable par's.
if(is.null(var))
{par(mfrow=c(1,1))
 levelx<-object$data1$parameter$levelx
 de2<-weighted.mean(object$de2,w2)
 if(!is.null(object$ie2))
 {ie2<-apply(object$ie2,2,weighted.mean,w2)
  ie2_show<-NULL
  for (i in 1:(length(object$ie2_list)-1))
    ie2_show<-c(ie2_show,sum(ie2[object$ie2_list[[i+1]]]))
  names(ie2_show)<-object$ie2_list[[1]]
 }
 else {ie2=NULL
       ie2_show=NULL}
 if(!is.null(object$ie12))
 {ie12<-apply(object$ie12,2,weighted.mean,w2)
  ie12_show<-NULL
  for (i in 1:(length(object$ie1_list)-1))
    ie12_show<-c(ie12_show,sum(ie12[object$ie1_list[[i+1]]]))
  names(ie12_show)<-object$ie1_list[[1]]
 }
 else {ie12=NULL
       ie12_show=NULL}
 te2<-sum(c(de2,ie2,ie12))
 re<-c(de2,ie2_show,ie12_show)/te2
 d<-order(re)
 name1<-c("DE",colnames(object$ie2),colnames(object$ie12))[d]
 barplot(re[d],horiz=TRUE,xlab="Relative Effects of Level 2 Mediation Effects",
         names=name1[d], cex.names=0.9,beside=FALSE,cex.axis=0.9,las=1,xlim=range(re),
         col = rainbow(length(d), start = 3/6, end = 4/6))
 if(object$data1$parameter$levelx==1)
 {readline("Press <return> to continue") 
  de1<-two(object$de1,object$data1$parameter$level,object$data1$parameter$weight)
  de1<-weighted.mean(de1,w2)
  if(!is.null(object$ie1))
  {ie1<-two(object$ie1,object$data1$parameter$level,object$data1$parameter$weight)
   ie1<-apply(ie1,2,weighted.mean,w2)
   ie1_show<-NULL
   for (i in 1:(length(object$ie1_list)-1))
     ie1_show<-c(ie1_show,sum(ie1[object$ie1_list[[i+1]]]))
   names(ie1_show)<-object$ie1_list[[1]]
  } 
  else {ie1=0
        ie1_show=0}
  te1<-sum(c(de1,ie1))
  re<-c(de1,ie1_show)/te1
  d<-order(re)
  name1<-c("DE",colnames(object$ie1))[d]
  barplot(re[d],horiz=TRUE,xlab="Relative Effects of Level 1 Mediation Effects",
          names=name1[d], cex.names=0.9,beside=FALSE,cex.axis=0.9,las=1,xlim=range(re),
          col = rainbow(length(d), start = 3/6, end = 4/6))
 } 
}
else{
  m<-object$m[,var]
  y<-predict(object$f1)
  a1<-(1:length(object$ie1_list[[1]]))[object$ie1_list[[1]]==var]
  a2<-(1:length(object$ie2_list[[1]]))[object$ie2_list[[1]]==var]
  l1<-2
  if(length(a1)>0)
  {if(object$data1$parameter$levelx==1)
  {x.ord<-order(object$x)
   m.ord<-order(m)
   par(mfcol=c(3,length(object$ie1_list[[a1+1]])))
   for(i in object$ie1_list[[a1+1]])
   {plot(object$x[x.ord],object$ie1[x.ord,i],xlab="x",ylab=paste("IE of",
                                                                 colnames(object$ie1)[i],sep=" "), type="l",
         main=paste("Level 1 indirect effect of", colnames(object$ie1)[i],sep=" "))
    plot(object$x[x.ord],object$ie1_2[x.ord,i],type="l", xlab="x",ylab=colnames(object$ie1)[i],
         main=paste("Differential effect of x and",colnames(object$ie1)[i],sep=" "))
    if(!cate)
      plot(m[m.ord],y[m.ord],type="l",xlab=var,ylab="y",
           main=paste("Predicted relationship between y and",var,sep=" "))
    else
      boxplot(y[object$ie1_1[,i]==0], y[object$ie1_1[,i]!=0], 
              names=c(paste("not",colnames(object$ie1)[i],sep=" "),
                      colnames(object$ie1)[i]),
              main=paste("Predicted relationship \n between y and",colnames(object$ie1)[i],sep=" "))}
   readline("Press <return> to continue") 
  }
  if(cate)
    for(i in object$ie1_list[[a1+1]])
      m<-cbind(m,ifelse(object$ie1_1[,i]==0,0,1))
  ie2_1<-as.matrix(object$ie12_1[,object$ie1_list[[a1+1]]])
  ie2_2<-as.matrix(object$ie12_2[,object$ie1_list[[a1+1]]])
  ie2<-as.matrix(object$ie12[,object$ie1_list[[a1+1]]])
  colnames(ie2_1)<-colnames(object$ie12)[object$ie1_list[[a1+1]]]
  colnames(ie2_2)<-colnames(object$ie12)[object$ie1_list[[a1+1]]]
  colnames(ie2)<-colnames(object$ie12)[object$ie1_list[[a1+1]]]
  l1<-1
  }
  else
  {     ie2_1<-as.matrix(object$ie2_1[,object$ie2_list[[a2+1]]])
        ie2_2<-as.matrix(object$ie2_2[,object$ie2_list[[a2+1]]])
        ie2<-as.matrix(object$ie2[,object$ie2_list[[a2+1]]])
        colnames(ie2_1)<-colnames(object$ie2_1)[object$ie2_list[[a2+1]]]
        colnames(ie2_2)<-colnames(object$ie2_1)[object$ie2_list[[a2+1]]]
        colnames(ie2)<-colnames(object$ie2_1)[object$ie2_list[[a2+1]]]
  }
  x.j.ord<-order(object$x.j)
  m.j<-two(m,object$data1$parameter$level)
  y.j<-two(y,object$data1$parameter$level)
  par(mfcol=c(3,ncol(ie2)))
  for(i in 1:ncol(ie2))
  {plot(object$x.j[x.j.ord],ie2[x.j.ord,i],xlab="x.j",ylab=paste("IE of",
                                                                 colnames(ie2)[i],sep=" "), type="l",
        main=paste("Level 2 indirect effect of", colnames(ie2)[i],sep=" "))
   plot(object$x.j[x.j.ord],ie2_2[x.j.ord,i],type="l", xlab="x",
        ylab=colnames(ie2)[i],
        main=paste("Differential effect of level 2 x and",colnames(ie2)[i],sep=" "))
   if(!cate)
     plot(m.j,y.j,xlab=paste("level 2",var,sep=" "),ylab="level 2 y",
          main=paste("Predicted relationship between level 2 y and",var,sep=" "))
   else if(l1==1)
     plot(m.j[,i+1],y.j,xlab=paste("level 2",var,sep=" "),ylab="level 2 y",
          main=paste("Predicted relationship \n between level 2 y and",colnames(ie2)[i],sep=" "))    
   else
     boxplot(y.j[ie2_1[,i]==0], y.j[ie2_1[,i]!=0], 
             names=c(paste("not",colnames(ie2)[i],sep=" "),
                     colnames(ie2)[i]),
             main=paste("Predicted relationship \n between level 2 y and",colnames(ie2)[i],sep=" "))
  }
}
par(op)
}


summary.mlma.boot<-function(object,...,alpha=0.05, RE=FALSE)
{summarize.boot<-function(vector,n,a1,a2,b1,b2)  #vector is stacked results from bootstrap samples
  #n is the number of elements from each bootstrap sample
{mat<-matrix(vector,n,length(vector)/n)
 all<-apply(mat,2,mean)
 c(mean=mean(all),sd=sd(all),upbd=mean(all)+b2*sd(all),
   lwbd=mean(all)+b1*sd(all), upbd.quan=quantile(all,a2),
   lwbd.quan=quantile(all,a1))
}
a1<-alpha/2
a2<-1-a1
b1<-qnorm(a1)
b2<-qnorm(a2)
direct.effect1<-NULL
direct.effect2<-NULL
indirect.effect1<-NULL
indirect.effect2<-NULL
indirect.effect12<-NULL
total.effect1<-NULL
total.effect2<-NULL
if(!is.null(object$de1))
{de1<-apply(object$de1,2,mean)
 direct.effect1<-c(est=mean(object$full$de1),mean=mean(de1),sd=sd(de1),
                   upbd=mean(de1)+b2*sd(de1),lwbd=mean(de1)+b1*sd(de1),
                   upbd.quan=quantile(de1,a2),lwbd.quan=quantile(de1,a1))}
if(!is.null(object$ie1))
{n<-nrow(object$full$ie1)
 indirect.effect1<-rbind(est=apply(object$full$ie1,2,mean),
                         apply(object$ie1,2,summarize.boot,n,a1,a2,b1,b2))
}
if(!is.null(object$de1))
{te1<-apply(cbind(c(object$de1),object$ie1),1,sum)
 total.effect1<-c(est=sum(mean(object$full$de1),ifelse(is.null(object$ie1),0,apply(object$full$ie1,2,mean))),
                  summarize.boot(te1,n,a1,a2,b1,b2))   
}
if(!is.null(object$de2))
{de2<-apply(object$de2,2,mean)
 direct.effect2<-c(est=mean(object$full$de2),mean=mean(de2),sd=sd(de2),
                   upbd=mean(de2)+b2*sd(de2),lwbd=mean(de2)+b1*sd(de2),
                   upbd.quan=quantile(de2,a2),lwbd.quan=quantile(de2,a1))}
if(!is.null(object$ie12))
{n<-nrow(object$full$ie12)
 indirect.effect12<-rbind(est=apply(object$full$ie12,2,mean),
                          apply(object$ie12,2,summarize.boot,n,a1,a2,b1,b2))}
if(!is.null(object$full$ie2))
{n<-nrow(object$full$ie2)
 indirect.effect2<-rbind(est=apply(object$full$ie2,2,mean),
                         apply(object$ie2,2,summarize.boot,n,a1,a2,b1,b2))}
if(!is.null(object$de2))
{te2<-apply(cbind(c(object$de2),object$ie12,object$ie2),1,sum)
 total.effect2<-c(est=sum(c(mean(object$full$de2),ifelse(is.null(object$ie12),0,apply(object$full$ie12,2,mean)),
                            ifelse(is.null(object$ie2),0,apply(object$full$ie2,2,mean)))),
                  summarize.boot(te2,n,a1,a2,b1,b2))   
}
if(!is.null(object$sum.boot1))
  {combine.effect1<- apply(object$sum.boot1,2,summarize.boot,1,a1,a2,b1,b2)
   temp.re1<-object$sum.boot1[,-ncol(object$sum.boot1)]/object$sum.boot1[,ncol(object$sum.boot1)]
   combine.re1<-apply(temp.re1,2,summarize.boot,1,a1,a2,b1,b2)*100}
else
   {combine.effect1<-NULL
    combine.re1<-NULL}
combine.effect2<- apply(object$sum.boot2,2,summarize.boot,1,a1,a2,b1,b2)
temp.re2<-object$sum.boot2[,-ncol(object$sum.boot2)]/object$sum.boot2[,ncol(object$sum.boot2)]
combine.re2<-apply(temp.re2,2,summarize.boot,1,a1,a2,b1,b2)*100

a<-list(total.effect1=total.effect1, total.effect2=total.effect2,
            direct.effect1=direct.effect1,direct.effect2=direct.effect2,
            indirect.effect1=indirect.effect1,indirect.effect2=indirect.effect2,
            indirect.effect12=indirect.effect12, 
        combine.effect=list(level1=combine.effect1,level2=combine.effect2),
        combine.re=list(level1=combine.re1,level2=combine.re2),
        levelx=object$levelx,RE=RE)
class(a)<-"summary.mlma.boot"
return(a)
}

print.summary.mlma.boot<-function(x,...)
{if(x$RE)
{if(x$levelx==1)
{cat("MLMA Analysis: Estimated Relative Effects at level 1 (%):\n")
  print(x$combine.re$level1)}
  cat("MLMA Analysis: Estimated Relative Effects at level 2 (%):\n")
  print(x$combine.re$level2)
}
  else
  {if(x$levelx==1)
 {cat("MLMA Analysis: Estimated Mediation Effects at level 1:\n")
  print(x$combine.effect$level1)}
 cat("MLMA Analysis: Estimated Mediation Effects at level 2:\n")
 print(x$combine.effect$level2)}
}

plot.mlma.boot<-function(x,..., var=NULL, alpha=0.05,quant=FALSE) #quantile=True to use quantile CIS
{object<-x

two<-function(x, level, weight=rep(1,nrow(as.matrix(x))))
{x<-as.matrix(x)
levels<-unique(level[!is.na(level)])
x2<-matrix(NA,length(levels),dim(x)[2])
for(i in 1:length(levels))
{cho<-(level==levels[i])
if(sum(cho)>0)
{if(sum(cho)==1)
  x2[i,]<-x[level==levels[i],]
else
{temp<-as.matrix(x[level==levels[i],])
weight1<-weight[level==levels[i]]
x2[i,]<-apply(temp,2,weighted.mean,weight1,na.rm=TRUE)}}}
colnames(x2)<-colnames(x)
x2
}


  boot.ci<-function(x,mat,cri_val) #the mat is the booted results with row be different x, and columns diff boot
                                 #cri_val is the critical value
  {x.uniq<-sort(unique(x))
   mn<-NULL
   upbd<-NULL
   lwbd<-NULL
   for (i in x.uniq)
   {sd_dev<-sd(as.vector(mat[x==i,]))
    mn1<-mean(as.vector(mat[x==i,]))
    upbd<-c(upbd,mn1+cri_val*sd_dev)
    lwbd<-c(lwbd,mn1-cri_val*sd_dev)
    mn<-c(mn,mn1)}
   return(data.frame(x=x.uniq,F=mn,L=lwbd,U=upbd))
  }
  plot_ci<-function(df,main="IE",xlab="x",ylab="IE")
  {plot(df$x, df$F, ylim = range(c(df$L,df$U)), type = "l",main=main,xlab=xlab,ylab=ylab)
    polygon(c(df$x,rev(df$x)),c(df$L,rev(df$U)),col = "grey75", border = FALSE)
    lines(df$x, df$F, lwd = 2)
    lines(df$x, df$U, col="red",lty=2)
    lines(df$x, df$L, col="red",lty=2)}
 op <- par(no.readonly = TRUE) # the whole list of settable par's.
 levelx<-object$levelx
   if(is.null(var))
 {summary.obj<-summary(object,alpha)
  if(levelx==1)
      par(mfrow=c(2,1))
  else par(mfrow=c(1,1))
  re<-summary.obj$combine.re$level2[1,]/100
  if(quant)
   {upper<-summary.obj$combine.re$level2[5,]/100
    lower<-summary.obj$combine.re$level2[6,]/100}
  else
    {lower<-summary.obj$combine.re$level2[4,]/100
     upper<-summary.obj$combine.re$level2[3,]/100}
  d<-order(re)
  bp <- barplot2(re[d], horiz = TRUE, main="Relative Effects of Level 2 Mediation Effects", 
                 names.arg=names(re)[d],plot.ci = TRUE, ci.u = upper[d], ci.l = lower[d],
                 cex.names=0.9,beside=FALSE,cex.axis=0.9,las=1,xlim=range(c(upper,lower)),
                 col = rainbow(length(d), start = 3/6, end = 4/6))
  if(levelx==1)
   {re<-summary.obj$combine.re$level1[1,]/100
    if(quant)
    {upper<-summary.obj$combine.re$level1[5,]/100
     lower<-summary.obj$combine.re$level1[6,]/100}
    else
    {lower<-summary.obj$combine.re$level1[4,]/100
     upper<-summary.obj$combine.re$level1[3,]/100}
    d<-order(re)
    bp <- barplot2(re[d], horiz = TRUE, main="Relative Effects of Level 2 Mediation Effects",
                   names.arg=names(re)[d],plot.ci = TRUE, ci.u = upper[d], ci.l = lower[d],
                   cex.names=0.9,beside=FALSE,cex.axis=0.9,las=1,xlim=range(c(upper,lower)),
                   col = rainbow(length(d), start = 3/6, end = 4/6))
  } 
  }
 else{
   a1<-(1:length(object$full$ie1_list[[1]]))[object$full$ie1_list[[1]]==var]
   a2<-(1:length(object$full$ie2_list[[1]]))[object$full$ie2_list[[1]]==var]
   x<-object$x
   x.j<-two(x,object$level,object$weight)
   cri_val<-qnorm((1+alpha)/2)
   if(length(a1)!=0 & levelx==1)
    {par(mfrow=c(2,length(object$full$ie1_list[[a1+1]])))
      for(i in object$full$ie1_list[[a1+1]])
       {temp.name<-colnames(object$full$ie1)[i]
        ie1<-boot.ci(x,matrix(object$ie1[,i],nrow=length(x)),cri_val)
        plot_ci(ie1,main=paste("Level 1 IE of", temp.name,sep=" "))
        temp.name<-colnames(object$full$ie12)[i]
        ie12<-boot.ci(x.j, matrix(object$ie12[,i],nrow=length(x.j)),cri_val)
        plot_ci(ie12,main=paste("Level 2 IE of", temp.name,sep=" "),xlab="x.j", ylab="IE12")
       }}
                    ############################################
  else if(length(a1)!=0)
   {par(mfrow=c(1,length(object$full$ie1_list[[a1+1]])))
    for(i in object$full$ie1_list[[a1+1]])
    {temp.name<-colnames(object$full$ie12)[i]
     ie12<-boot.ci(x.j, matrix(object$ie12[,i],nrow=length(x.j)),cri_val)
     plot_ci(ie12,main=paste("Level 2 IE of", temp.name,sep=" "),xlab="x.j", ylab="IE12")
    }}
  else
  {par(mfrow=c(1,length(object$full$ie2_list[[a2+1]])))
   for(i in object$full$ie2_list[[a2+1]])
   {temp.name<-colnames(object$full$ie2)[i]
    ie2<-boot.ci(x.j, matrix(object$ie2[,i],nrow=length(x.j)),cri_val)
    plot_ci(ie2,main=paste("Level 2 IE of", temp.name,sep=" "),xlab="x.j", ylab="IE12")
   }}}
par(op)
}

