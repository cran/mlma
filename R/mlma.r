# data clean #start line #292
data.org<-function(x, m, levely=1,  y=NULL, levelx=NULL, xref=NULL, 
                   l1=NULL,l2=NULL, c1=NULL, c1r=rep(1,length(c1)), #levelx is the level of x
                   c2=NULL, c2r=rep(1,length(c2)), f01y=NULL, f10y=NULL,                  #level is the level of observations
                   f02ky=NULL, f20ky=NULL, f01km1=NULL, f01km2=NULL, f10km=NULL,          #weight is the level 1 weight of cases
                   level=1:nrow(as.matrix(x)), weight=NULL)                                        #weight2 is the level 2 weight of cases, weight2=rep(1,length(unique(level[!is.na(level)])))   
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
 cbind(x1,x)
}



two<-function(x, level, weight=rep(1,nrow(as.matrix(x))))
{x2<-NULL
if(is.null(dim(x)))
  x2=cbind(x2,(aggregate(as.numeric(x)~level, na.action=na.pass, 
                         FUN=function(x) c(mean=weighted.mean(x,weight=weight,na.rm=T))))[,2])
else
  for(i in 1:ncol(x))
    x2<-cbind(x2,(aggregate(as.numeric(x[,i])~level, na.action=na.pass,
                            FUN=function(x) c(mean=weighted.mean(x,weight=weight,na.rm=T))))[,2])
  colnames(x2)<-colnames(x)
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

cattobin<-function(m1y,m1y.der,m1,m,cat1,cat2=rep(1,length(cat1)), level=rep(1,dim(m)[1]),weight=rep(1,nrow(as.matrix(m))))
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
{a<-as.factor(m[,i])
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
 m1y.der<-cbind(m1y.der,matrix(1,nrow(f),ncol(f)))
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

find_level<-function(vari, level)    #a function to identify mediators to l1, l2, c1, c2
{a<-table(level)
b<-sort(unique(level))
l=2
for (i in 1:length(b))
  if(a[i]>1)
  {temp1<-vari[level==b[i]]
  temp2<-temp1[!is.na(temp1)]
  if(length(temp2)!=0 & sum(temp2!=temp2[1])>0)
  {l=1
  break}
  }
if(l==1 & (is.factor(vari) | is.character(vari) | nlevels(as.factor(vari))==2))
  return (list(1, levels(as.factor(vari))[1]))  #c1, c1r
else if(l==1)
  return (list(2,NA))   #l1
else if(l==2 & (is.factor(vari) | is.character(vari) | nlevels(as.factor(vari))==2))
  return (list(3, levels(as.factor(vari))[1]))  #c2, c2r
else
  return(list(4,NA))
}

n1<-nrow(as.matrix(m))                                   
x=data.frame(x)
mnames<-colnames(as.matrix(m))                #change variable names to columns in m
if(is.character(l1))
  l1<-unlist(sapply(l1,grep,mnames))
if(is.character(c1))
  c1<-unlist(sapply(c1,grep,mnames))
if(is.character(l2))
  l2<-unlist(sapply(l2,grep,mnames))
if(is.character(c2))
  c2<-unlist(sapply(c2,grep,mnames))
if(is.null(c(l1,l2,c1,c2)))                   #identify the type and level of mediators
{temp.m=data.frame(m)
 temp.mnames<-mnames
  for (i in 1:ncol(temp.m))
  {temp.1<-find_level(temp.m[,i],level)
   col.num<-grep(temp.mnames[i],mnames)
   if (temp.1[[1]]==1)
     {m[,col.num]<-as.factor(m[,col.num])
      c1<-c(c1,col.num)
      c1r<-c(c1r,temp.1[[2]])}
   else if (temp.1[[1]]==2)
     l1<-c(l1,col.num)
   else if(temp.1[[1]]==3)
     {m[,col.num]<-as.factor(m[,col.num])
      c2<-c(c2,col.num)
      c2r<-c(c2r,temp.1[[2]])}
   else l2<-c(l2,col.num)
  }
}
else if (ncol(data.frame(m)[,-c(l1,l2,c1,c2)])!=0)
 {temp.m=data.frame(m[,-c(l1,l2,c1,c2)])
  temp.mnames<-mnames[-c(l1,l2,c1,c2)]
  for (i in 1:ncol(temp.m))
  {temp.1<-find_level(temp.m[,i],level)
  col.num<-grep(temp.mnames[i],mnames)
  if (temp.1[[1]]==1)
  {m[,col.num]<-as.factor(m[,col.num])
  c1<-c(c1,col.num)
  c1r<-c(c1r,temp.1[[2]])}
  else if (temp.1[[1]]==2)
    l1<-c(l1,col.num)
  else if(temp.1[[1]]==3)
  {m[,col.num]<-as.factor(m[,col.num])
  c2<-c(c2,col.num)
  c2r<-c(c2r,temp.1[[2]])}
  else l2<-c(l2,col.num)
  }
}

nx=ncol(x)   #identify levels 1 and 2 exposures
binx<-rep(F,nx)
if(is.null(levelx))
 for (i in 1:nx){
 temp.1<-find_level(x[,i],level)
 if(temp.1[[1]]<=2)
  {levelx=c(levelx,1)
   if(temp.1[[1]]==1)
   {binx[i]=T
     if(is.null(xref))
       x[,i]<-ifelse(x[,i]==temp.1[[2]],0,1)
     else
       x[,i]<-ifelse(x[,i]==xref,0,1)}
 }
else
{levelx=c(levelx,2)
 if(temp.1[[1]]==3)
  {binx[i]=T
   if(is.null(xref))
    x[,i]<-ifelse(x[,i]==temp.1[[2]],0,1)
  else
    x[,i]<-ifelse(x[,i]==xref,0,1)}
}
}
else 
 for (i in 1:nx)
 if(nlevels(as.factor(x[,i]))==2)  
 {binx[i]=T      #define the reference group for binary predictor
  if(is.null(xref))
    x[,i]<-ifelse(x[,i]==levels(as.factor(x[,i]))[1],0,1)
  else
    x[,i]<-ifelse(x[,i]==xref,0,1)
 }

if (!is.null(c(l2,c2)) & sum(levelx==2)==0) #if there are level 2 mediator(s), but no level 2 exposure(s)
 {x.temp=NULL
  x.names=names(x)
  for (i in 1:nx)
   x.temp=cbind(x.temp,one2two(x[,i],level))
  colnames(x.temp)=paste(rep(x.names,each=2),c("_agregated",""),sep="")
  x=x.temp
  binx=c(t(matrix(c(rep(F,nx),binx),nx)))
  levelx=rep(c(2,1),nx)
  nx=2*nx
}

if(is.null(levely))
{temp.1<-find_level(y,level)
if(temp.1[[1]]<=2)
 levely=1
else
 levely=2
}

if(is.null(weight))
  weight=rep(1,n1)

xnames=colnames(x)
if(is.character(f02ky[[1]]))
  f02ky[[1]]<-unlist(sapply(f02ky[[1]],grep,mnames))
if(is.character(f20ky[[1]]))
  f20ky[[1]]<-unlist(sapply(f20ky[[1]],grep,mnames))
if(is.character(f01y[[1]]))
  f01y[[1]]<-unlist(sapply(f01y[[1]],grep,xnames))
if(is.character(f10y[[1]]))
  f10y[[1]]<-unlist(sapply(f10y[[1]],grep,xnames))

lx<-matrix(0,nx,3)                       #create x to explain y
lx[,1]=levelx
x1=NULL
x1.der=NULL
k=1
x1.names=NULL
for (j in 1:nx)
 if(!(j %in% f01y[[1]]) & !(j %in% f01y[[1]]))
 {x1<-cbind(x1,x[,j])
  x1.der<-cbind(x1.der,rep(1,n1))
  lx[j,2]<-k
  lx[j,3]<-k
  k=k+1
  x1.names=c(x1.names,xnames[j])
  }
 else if (j %in% f01y[[1]])
 {p=match(j,f01y[[1]])
  x1<-cbind(x1,x2fx(x[,j],f01y[[1+p]])$values)
  if(is.null(x1.der))
    l=ncol(x1)
  else
    l=ncol(x1)-ncol(x1.der)
  x1.der<-cbind(x1.der,x2fdx(x[,j],f01y[[1+p]]))
  x1<-as.matrix(x1)
  lx[j,2]<-k
  lx[j,3]<-k+l-1
  k=k+l
  x1.names<-c(x1.names,paste(xnames[j],1:l,sep="."))
  x1.der<-as.matrix(x1.der)
 }
 else 
 {p=match(j,f10y[[1]])
  x1<-cbind(x1,x2fx(x[,j],f10y[[1+p]])$values)
  if(is.null(x1.der))
   l=ncol(x1)
  else
   l=ncol(x1)-ncol(x1.der)
  x1.der<-cbind(x1.der,x2fdx(x[,j],f10y[[1+p]]))
  x1<-as.matrix(x1)
  lx[j,2]<-k
  lx[j,3]<-k+l-1
  k=k+l
  x1.names<-c(x1.names,paste(xnames[j],1:l,sep="."))
  x1.der<-as.matrix(x1.der)
 }
 colnames(x1)=x1.names
 colnames(x1.der)=x1.names
 if(!is.null(l2))        #create level 2 m variables to explain y
 {m2y<-as.matrix(m[,l2])
  colnames(m2y)<-colnames(m)[l2]
  m2y.der<-matrix(1,n1,length(l2))
  colnames(m2y.der)<-colnames(m)[l2]
  m2<-as.list(c(1,1:length(l2)))
  m2[[1]]<-l2
  if(!is.null(f02ky)) 
    for (i in 1:length(f02ky[[1]]))
    {a<-f02ky[[1]][i]
     if(a %in% l2)
     {b<-match(a,l2)
      d<-as.matrix(x2fx(m[,a],f02ky[[i+1]])$values)
      d.der<-as.matrix(x2fdx(m[,a],f02ky[[i+1]]))
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
 {m1y<-as.matrix(as.matrix(m[,l1]))
  m1y.der<-matrix(1,n1,ncol(m1y))
  colnames(m1y)<-mnames[l1]
  colnames(m1y.der)<-mnames[l1]
  m1<-as.list(c(1,1:length(l1)))
  m1[[1]]<-l1

  if(!is.null(f20ky)) 
    for (i in 2:length(f20ky))
    {a<-f20ky[[1]][i-1]
     b<-match(a,l1)
     d<-as.matrix(x2fx(m[,a],f20ky[[i]])$values)
     d.der<-as.matrix(x2fdx(m[,a],f20ky[[i]]))
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
 }
 else 
 {m1y<-NULL
  m1y.der<-NULL
  m1<-NULL
 }
  
if(!is.null(c1))                  #binarize level 1 categorical variables and corresponding level 2 variables to explain y 
 {temp<-cattobin(m1y,m1y.der,m1,m,c1,c1r,level,weight)
  m1y<-temp$m1y
  m1y.der<-temp$m1y.der
  m1<-temp$m1
 }

 lc2<-c(l2,c2)
 f01km2.2=f01km2
 if(!is.null(lc2))                  #create x variables to explain level 2 mediators
 {xm2<-as.matrix(two(as.matrix(x[,levelx==2]),level))
  colnames(xm2)<-paste(xnames[levelx==2],"2",sep=".")
  xm2.der<-matrix(1,nrow(xm2),ncol(xm2))
  colnames(xm2.der)<-colnames(xm2)
  
  m.2<-two(as.matrix(m[,lc2]),level)
  colnames(m.2)<-mnames[lc2]
  fm22<-as.list(rep(1,length(lc2)+1))
  fm22[[1]]<-lc2
  for(i in 1:length(lc2))
    fm22[[i+1]]=1:ncol(xm2)
  
  if(!is.null(f01km2))
  {k=unique(f01km2[[1]][,2])
   for (l in k){
    temp.2=as.matrix(two(as.matrix(x[,k]),level))
    temp<-(2:length(f01km2))[f01km2[[1]][,2]==l]
    allfun=f01km2[[temp[1]]]
    if (length(temp)>1)
     for(i in 2:length(temp))
       allfun<-c(allfun,f01km2[[temp[i]]])
    unifun<-unique(allfun)
    unifun1<-unifun[unifun!="x"]
    unifun2<-c("x",unifun1)
    d_d<-x2fx(temp.2,unifun1)
    d.der<-as.matrix(x2fdx(temp.2,unifun1))
    d<-as.matrix(d_d$values)
    temp.xm2<-c(colnames(xm2),paste(paste(xnames[l],"2",sep="."),1:ncol(d),sep="."))
    xm2.der<-cbind(xm2.der,d.der)
    colnames(xm2.der)<-temp.xm2
    temp2=match(paste(xnames[l],"2",sep="."),temp.xm2)
    col_fun<-cbind(c(temp2,temp2),d_d$col_fun+ncol(xm2))
    xm2<-cbind(xm2,d)
    colnames(xm2)<-temp.xm2
    for(i in temp)
    {ttemp<-order_char(unifun2,f01km2[[i]])
     ttemp1<-NULL
     for (j in 1:length(ttemp))
       ttemp1<-c(ttemp1,col_fun[1,ttemp[j]]:col_fun[2,ttemp[j]])
     fm22[[(1:length(lc2))[fm22[[1]]==f01km2[[1]][i-1,1]]+1]]<-c(fm22[[(1:length(lc2))[fm22[[1]]==f01km2[[1]][i-1,1]]+1]][-temp2],ttemp1)
     f01km2.2[[i]]<-ttemp1
     }
  }}}
 else
 {m.2<-NULL
  xm2<-NULL
  xm2.der<-NULL
  fm22<-NULL}
 
 lc1<-c(l1,c1)
 f01km1.2=f01km1
 f10km.2=f10km
 if(!is.null(lc1))                  
 {xm1<-x
  xm1.der<-matrix(1,n1,ncol(x))
  colnames(xm1.der)<-colnames(x)
  if(sum(levelx==2)!=0){
  fm12<-as.list(rep(1,length(lc1)+1))  #create level 2 x variables to explain level 1 mediatiors
  fm12[[1]]<-lc1
  for (i in 1:length(lc1))
    fm12[[i+1]]=(1:ncol(xm1))[levelx==2]
  
  if(!is.null(f01km1))
  {k=unique(f01km1[[1]][,2])
   for (l in k){
    temp<-(2:length(f01km1))[f01km1[[1]][,2]==l]
    allfun<-f01km1[[temp[1]]]
    if (length(temp)>1)
     for(i in 2:length(temp))
      allfun<-c(allfun,f01km1[[temp[i]]])
    unifun<-unique(allfun)
    unifun1<-unifun[unifun!="x"]
    unifun2<-c("x",unifun1)
    d_d<-x2fx(x[,l],unifun1)
    d<-as.matrix(d_d$values)
    temp.xm1<-c(colnames(xm1),paste(xnames[l],1:ncol(d),sep="."))
    xm1.der<-cbind(xm1.der,x2fdx(x[,l],unifun1))
    colnames(xm1.der)<-temp.xm1
    col_fun<-cbind(c(l,l),d_d$col_fun+ncol(xm1))
    xm1<-cbind(xm1,d)
    colnames(xm1)<-temp.xm1
    for(i in temp)
    {ttemp<-order_char(unifun2,f01km1[[i]])
     ttemp1<-NULL
     for (j in 1:length(ttemp))
      ttemp1<-c(ttemp1,col_fun[1,ttemp[j]]:col_fun[2,ttemp[j]])
     fm12[[(1:length(lc1))[fm12[[1]]==f01km1[[1]][i-1,1]]+1]]<-
       c(fm12[[(1:length(lc1))[fm12[[1]]==f01km1[[1]][i-1,1]]+1]][-match(l,fm12[[(1:length(lc1))[fm12[[1]]==f01km1[[1]][i-1,1]]+1]])],ttemp1)
     f01km1.2[[i]]<-ttemp1}
   }}}
  else{
    fm12=NULL
  }
 if(sum(levelx==1)>0){
  fm11<-as.list(rep(1,length(lc1)+1))  #create level 1 x variables to explain level 1 mediatiors
  for(i in 1:length(lc1))
   fm11[[i+1]]=(1:ncol(x))[levelx==1]
 fm11[[1]]<-lc1
 if(!is.null(f10km))
 {k=unique(f10km[[1]][,2])
  for(l in k){
   temp<-(2:length(f10km))[f10km[[1]][,2]==l]
   allfun<-f10km[[temp[[1]]]]
   if (length(temp)>1)
    for(i in 2:length(temp))
     allfun<-c(allfun,f10km[[temp[i]]])
   unifun<-unique(allfun)
   unifun1<-unifun[unifun!="x"]
   unifun2<-c("x",unifun1)
   d_d<-x2fx(x[,l],unifun1)
   d<-as.matrix(d_d$values)
   temp.xm1<-c(colnames(xm1),paste(xnames[l],1:ncol(d),sep="."))
   xm1.der<-cbind(xm1.der,x2fdx(x[,l],unifun1))
   colnames(xm1.der)=temp.xm1
   col_fun<-cbind(c(l,l),d_d$col_fun+ncol(xm1))
   xm1=cbind(xm1,d)
   colnames(xm1)=temp.xm1
   for(i in temp)
   {ttemp<-order_char(unifun2,f10km[[i]])
    ttemp1<-NULL
    for (j in 1:length(ttemp))
     ttemp1<-c(ttemp1,col_fun[1,ttemp[j]]:col_fun[2,ttemp[j]])
    fm11[[(1:length(lc1))[fm11[[1]]==f10km[[1]][i-1,1]]+1]]<-
      c(fm11[[(1:length(lc1))[fm11[[1]]==f10km[[1]][i-1,1]]+1]][-match(l,fm11[[(1:length(lc1))[fm11[[1]]==f10km[[1]][i-1,1]]+1]])],ttemp1)
    f10km.2[[i]]<-ttemp1}}
 }}
 else
 {fm11=NULL}
 }
 else
 {xm1<-NULL
 xm1.der<-NULL
 fm12<-NULL
 fm11<-NULL}

 
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

list(x1=x1, x1.der=x1.der, lx=lx, m1y=m1y, m1y.der=m1y.der,
     m1=m1, m2y=m2y, m2y.der=m2y.der, m2=m2, xm1=xm1, xm1.der=xm1.der,
     fm11=fm11, fm12=fm12, m.2=m.2, xm2=xm2, xm2.der=xm2.der, fm22=fm22,binx=binx,
     f01km1.2=f01km1.2,f01km2.2=f01km2.2,f10km.2=f10km.2,
     parameter=list(levelx=levelx, levely=levely, mnames=mnames,l1 = l1, l2 = l2,  
                    c1 = c1, c1r = c1r, c2 = c2, c2r = c2r, f01y = f01y, level=level,
                    f10y = f10y, f02ky = f02ky, f20ky = f20ky, f01km1 = f01km1, 
                    f01km2 = f01km2, f10km = f10km, weight = weight,x=x,m=m))
}

#multilevel mediation analysis
mlma<-function(y, data1=NULL, x=data1$parameter$x, m=data1$parameter$m, 
               yref=NULL, xref=NULL, levelx=data1$parameter$levelx, levely=data1$parameter$levely,
               l1=data1$parameter$l1,l2=data1$parameter$l2, 
               c1=data1$parameter$c1, #levelx is the level of x
               c1r=data1$parameter$c1r, c2=data1$parameter$c2, c2r=data1$parameter$c2r, level=data1$parameter$level,  
               weight=rep(1,nrow(data.frame(x))), random="(1|level)", random.m1=NULL,intercept=TRUE, 
               family1=NULL, familym=vector("list",ncol(m)), covariates=NULL,cy1=NULL, cy2=NULL, cm=NULL,joint=NULL,
               f01y=data1$parameter$f01y, f10y=data1$parameter$f10y, f02ky=data1$parameter$f02ky, 
               f20ky=data1$parameter$f20ky, f01km1=data1$parameter$f01km1, f01km2=data1$parameter$f01km2, 
               f10km=data1$parameter$f10km,
               data2=NULL, x.new=NULL, m.new=NULL, level.new=level, weight.new=NULL,
               covariates.new=covariates,cov.mat=FALSE)                               
{
  two<-function(x, level, weight=rep(1,nrow(as.matrix(x))))
  {x2<-NULL
  if(is.null(dim(x)))
    x2=cbind(x2,(aggregate(as.numeric(x)~level, na.action=na.pass, 
                           FUN=function(x) c(mean=weighted.mean(x,weight=weight,na.rm=T))))[,2])
  else
    for(i in 1:ncol(x))
      x2<-cbind(x2,(aggregate(as.numeric(x[,i])~level, na.action=na.pass,
                              FUN=function(x) c(mean=weighted.mean(x,weight=weight,na.rm=T))))[,2])
    colnames(x2)<-colnames(x)
    x2
  }
  



one<-function(x,level)  #change a level 2 x to have the length of observations
{levels<-unique(level[!is.na(level)])
x1<-rep(NA,length(level))
for (i in 1:length(levels))
  x1[level==levels[i]]=x[i]
x1
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

order_char<-function(char1,char2)  #find the position of char2 in char1
{a<-1:length(char1)
b<-NULL
for (i in 1:length(char2))
  b<-c(b,a[char1==char2[i]])
b
}

find_level<-function(vari, level)    #a function to identify mediators to l1, l2, c1, c2
{a<-table(level)
b<-sort(unique(level))
l=2
for (i in 1:length(b))
  if(a[i]>1)
  {temp1<-vari[level==b[i]]
  temp2<-temp1[!is.na(temp1)]
  if(length(temp2)!=0 & sum(temp2!=temp2[1])>0)
  {l=1
  break}
  }
if(l==1 & (is.factor(vari) | is.character(vari) | nlevels(as.factor(vari))==2))
  return (list(1, levels(as.factor(vari))[1]))  #c1, c1r
else if(l==1)
  return (list(2,NA))   #l1
else if(l==2 & (is.factor(vari) | is.character(vari) | nlevels(as.factor(vari))==2))
  return (list(3, levels(as.factor(vari))[1]))  #c2, c2r
else
  return(list(4,NA))
}

d_link<-function(link.fun,x) #link.fun is the link function: a character string
{  x2fx<-function(x,func) #x is the list of original numerical vector, func is a vector of character functions. 
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
if(is.null(link.fun))
  result=rep(1,length(x))
else if (link.fun=="logit")
{funct="x*(1-x)"
result=x2fx(x,funct)$values}
else if(link.fun=="log")
{funct="x"
result=x2fx(x,funct)$values}
else if(link.fun--"inverse")
{funct="-x^2"
result=x2fx(x,funct)$values}
else
  result=x2fdx(x,"x")
result
}

generate.m<-function(x.new, data1,l1,c1,l2,c2,full,covariates.new,level.new,data2,cm)
{n.new<-nrow(x.new)
  n2.new<-length(unique(level.new))
  mnames=data1$parameter$mnames
  m.new<-matrix(0,n.new,length(mnames))
  colnames(m.new)<-mnames
  fm1<-full$fm1
  if(!is.null(l1))        #generate level 1 continuous mediators
  {for (i in 1:length(l1))
  {if(sum(cm[[1]]==l1[i])>0)
  {temp2<-match(l1[i],cm[[i]])+1
  temp.cov<-as.matrix(covariates.new[,cm[[temp2]]])
  colnames(temp.cov)<-colnames(covariates.new)[cm[[temp2]]]}
    else
      temp.cov<-NULL
    expl.m<-as.matrix(data2$xm1[,c(data1$fm12[[i+1]],data1$fm11[[i+1]])])
    colnames(expl.m)<-colnames(data2$xm1)[c(data1$fm12[[i+1]],data1$fm11[[i+1]])]
    expl.m<-cbind(expl.m,temp.cov)
    temp.data<-cbind(level=level.new,expl.m)
    colnames(temp.data)<-c("level",colnames(expl.m))
    j<-(1:length(fm1[[1]]))[fm1[[1]]==l1[i]]+1
    sd.m<-rnorm(n.new,mean=0,sd=as.data.frame(VarCorr(full$fm1[[j]]))[2,5])+         #add the randomn effects
      one(rnorm(n2.new,mean=0,sd=as.data.frame(VarCorr(full$fm1[[j]]))[1,5]),level.new)
    m.new[,l1[i]]<-predict(full$fm1[[j]],newdata=data.frame(temp.data),allow.new.levels=T)+sd.m
  }
  }
  if(!is.null(c1))      #generate level 1 categorical mediators
    for (i in 1:length(c1))
    {if(sum(cm[[1]]==c1[i])>0)
    {temp2<-(1:length(cm[[1]]))[cm[[1]]==c1[i]]+1  #may have more than one case?
    temp.cov<-as.matrix(covariates.new[,cm[[temp2]]])
    colnames(temp.cov)<-colnames(covariates.new)[cm[[temp2]]]}
      else
        temp.cov<-NULL
      if(is.null(data1$fm11))
        k1<-NULL
      else
        k1<-(1:length(data1$fm11[[1]]))[data1$fm11[[1]]==c1[i]]
      if(is.null(data1$fm12))
        k2<-NULL
      else
        k2<-(1:length(data1$fm12[[1]]))[data1$fm12[[1]]==c1[i]]
      expl.m<-as.matrix(data2$xm1[,c(data1$fm12[[k2+1]],data1$fm11[[k1+1]])])
      colnames(expl.m)<-colnames(data2$xm1)[c(data1$fm12[[k2+1]],data1$fm11[[k1+1]])]
      expl.m<-cbind(expl.m,temp.cov)
      temp.data<-cbind(level=level.new,expl.m)
      colnames(temp.data)<-c("level",colnames(expl.m))
      j<-(1:length(fm1[[1]]))[full$fm1[[1]]==c1[i]]+1
      if(length(j)==1)   #binary
      {tmp.logit<-predict(full$fm1[[j]],newdata=data.frame(temp.data),allow.new.levels=T)+
        one(rnorm(n2.new,mean=0,sd=as.data.frame(VarCorr(fm1[[j]]))[,5]),level.new)  #add the random effect
      tmp.p<-exp(tmp.logit)/(exp(tmp.logit)+1)
      tmp.m<-rbinom(n.new,1,prob=tmp.p)
      m.new[,c1[i]]<-ifelse(tmp.m==0,1,2)}
      else               #multi-categorical
      {tmp.logit<-rep(1,n.new)
      for (k in 1:length(j))
        tmp.logit<-cbind(tmp.logit,exp(predict(full$fm1[[j[k]]],newdata=data.frame(temp.data),allow.new.levels=T)+
                                         one(rnorm(n2.new,mean=0,sd=as.data.frame(VarCorr(fm1[[j[k]]]))[,5]),level.new))) #add the random effect
      tmp.m<-apply(tmp.logit,1,rmultinom,n=1,size=1)
      m.new[,c1[i]]<-(1:nrow(tmp.m))%*%tmp.m   #the reference group is always 1
      }
    }
 
  fm2<-full$fm2            #generate level 2 mediators
  cov.2<-NULL
  if(!is.null(covariates.new))
    cov.2<-two(covariates.new,level.new)
  if(!is.null(l2))          #generate level 2 mediators
    for (i in 1:length(l2))
    {if(sum(cm[[1]]==l2[i])>0)
    {temp2<-(1:length(cm[[1]]))[cm[[1]]==l2[i]]+1
    temp.cov<-as.matrix(cov.2[,cm[[temp2]]])
    colnames(temp.cov)<-colnames(covariates.new)[cm[[temp2]]]}
      else
        temp.cov<-NULL
      expl.m<-as.matrix(data2$xm2[,data1$fm22[[i+1]]])
      colnames(expl.m)<-colnames(data1$xm2)[data1$fm22[[i+1]]]
      expl.m<-cbind(expl.m,temp.cov)
      j<-(1:length(fm2[[1]]))[fm2[[1]]==l2[i]]+1
      tmp.m<-predict(fm2[[j]],newdata=data.frame(expl.m))+
        rnorm(n2.new,mean=0,sd=(summary(fm2[[j]]))$sigma)  #add the randomness
      m.new[,l2[i]]<-one(tmp.m,level.new) #lm was used
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
      expl.m<-as.matrix(data2$xm2[,data1$fm22[[k+1]]])
      colnames(expl.m)<-colnames(data1$xm2)[data1$fm22[[k+1]]]
      expl.m<-cbind(expl.m,temp.cov)
      j<-(1:length(fm2[[1]]))[fm2[[1]]==c2[i]]+1
      if(length(j)==1) #binary
      {tmp.m<-rbinom(n2.new,1,prob=predict(fm2[[j]],newdata=data.frame(expl.m),type="response"))
      m.new[,c2[i]]<-one(ifelse(tmp.m==0,1,2),level.new)}
      else               #multi-categorical
      {tmp.logit<-rep(1,n2.new)
      for (k in 1:length(j))
        tmp.logit<-cbind(tmp.logit,predict(fm2[[j[k]]],newdata=data.frame(expl.m)))
      tmp.m<-apply(tmp.logit,1,rmultinom,n=1,size=1)
      m.new[,c2[i]]<-one((1:nrow(tmp.m))%*%tmp.m,level)   #the reference group is always 1
      }
    }           #glm
  m.new
}
func1<-function(a,mat)
{as.vector(a%*%mat%*%a)}
func2<-function(a,mat)
{as.vector(a%*%mat*a)}
#n<-length(y)
#tt2<-cbind(data1$x1[,data1$l1x],covariates[,cy2])
#colnames(tt2)<-c(colnames(data1$x1)[data1$l1x],colnames(covariates)[cy2])
#tt1<-cbind(data1$x1[,data1$l2x],covariates[,cy1])
#colnames(tt1)<-c(colnames(data1$x1)[data1$l2x],colnames(covariates)[cy1])
#expl<-cbind(tt1,data1$m2y,tt2,data1$m1y)   

#n<-length(y)
surv=is.Surv(y)
biny<-FALSE

if(!surv){
temp.1<-find_level(y,level)
if(temp.1[[1]]<=2)
 {levely=1
  if(temp.1[[1]]==1)
   {biny=TRUE
    if(is.null(yref))
     y<-ifelse(y==temp.1[[2]],0,1)
    else
     y<-ifelse(y==yref,0,1)}
 }
else
 {levely=2
  if(temp.1[[1]]==3)
   {biny=TRUE
    if(is.null(yref))
     y<-ifelse(y==temp.1[[2]],0,1)
    else
     y<-ifelse(y==yref,0,1)}
 }
}

#organize the data if it has not been done
if(is.null(data1))
  {data1<-data.org(x=x, m=m, levely=levely, y=y, levelx=levelx, xref=xref, l1=l1,l2=l2, 
                   c1=c1, c1r=c1r, c2=c2, c2r=c2r, f01y=f01y, f10y=f10y,
                   f02ky=f02ky, f20ky=f20ky, f01km1=f01km1, f01km2=f01km2, f10km=f10km, 
                   level=level, weight=weight)
   x=data1$parameter$x
   m=data1$parameter$m
   levelx=data1$parameter$levelx
   levely=data1$parameter$levely
   l1=data1$parameter$l1
   l2=data1$parameter$l2
   c1=data1$parameter$c1
   c1r=data1$parameter$c1r
   c2=data1$parameter$c2
   c2r=data1$parameter$c2r}

mnames<-colnames(as.matrix(m))
x.names=colnames(as.matrix(x))

if(is.null(data2) & is.null(x.new))  #using the original data
{data2=data1
 x.new=x
 m.new=m
 weight.new=weight
}
else if (is.null(data2)) #when there are x.new, generate data2
{if(is.null(weight.new))
  weight.new=rep(1,nrow(x.new))
 if(is.null(m.new))
  m.new1=m[sample(1:nrow(m),size=nrow(x.new),replace=T),]
 else
  m.new1=m.new
 data2<-data.org(x=x.new, m=m.new1, levely=levely, y=NULL, levelx=levelx, xref=xref, l1=l1,l2=l2, 
                 c1=c1, c1r=c1r, c2=c2, c2r=c2r, f01y=f01y, f10y=f10y,
                 f02ky=f02ky, f20ky=f20ky, f01km1=f01km1, f01km2=f01km2, f10km=f10km, 
                 level=level.new, weight=weight.new)
 x.new=data2$parameter$x
 if(!is.null(m.new))
  m.new=data2$parameter$m
}
else #when all .new comes from data2
{x.new=data2$parameter$x
 m.new=data2$parameter$m
 weight.new=data2$parameter$weight
}

n=nrow(x.new)

if(is.null(m.new)) #if m.new is null, create one based on full models and x.new.
{full<-mlma(y=y, data1=data1, weight=weight, 
            random=random, random.m1=random.m1,intercept=intercept, family1=family1, familym=familym,
            covariates=covariates, cy1=cy1, cy2=cy2, cm=cm, joint=joint, data2=data2,x.new=x.new, 
            m.new=m.new1, level.new=level.new,weight.new=weight.new,covariates.new=covariates.new)
 m.new<-generate.m(x.new,data1,l1,c1,l2,c2,full,covariates.new,level.new,data2,cm)
 data2<-data.org(x=x.new, m=m.new, levely=levely, y=NULL, levelx=levelx, xref=xref, l1=l1,l2=l2, 
                c1=c1, c1r=c1r, c2=c2, c2r=c2r, f01y=f01y, f10y=f10y,
                f02ky=f02ky, f20ky=f20ky, f01km1=f01km1, f01km2=f01km2, f10km=f10km, 
                level=level.new, weight=weight.new)
 m.new=data2$parameter$m
}
if(sum(data1$lx[,1]==1)!=0)  #lx[,1]indicate the level of x
 {l1x=NULL
  for (i in (1:nrow(data1$lx))[data1$lx[,1]==1])
    l1x=c(l1x,data1$lx[i,2]:data1$lx[i,3])
  tt2<-as.matrix(data1$x1[,l1x])   #tt2 is the matrix with level 1 exposures
  colnames(tt2)<-colnames(data1$x1)[l1x]}
else
  tt2<-NULL

if(sum(data1$lx[,1]==2)!=0)  #ith row of lx is for the ith column of x
{l2x=NULL
 for (i in (1:nrow(data1$lx))[data1$lx[,1]==2])
  l2x=c(l2x,data1$lx[i,2]:data1$lx[i,3])
 tt1<-as.matrix(data1$x1[,l2x])     #tt1 is the matrix with level 2 exposures
 colnames(tt1)<-colnames(data1$x1)[l2x]}
else 
  tt1=NULL

if(!is.null(c(cy1,cy2)))           ###added covariates to explain y
 {cova<-as.matrix(covariates[,c(cy1,cy2)])
  colnames(cova)<-colnames(covariates)[c(cy1,cy2)]}
else
  cova<-NULL
if(is.null(cova))
  expl<-cbind(data1$x1,data1$m2y,data1$m1y)  # all variables to explain y
else
  expl<-cbind(data1$x1,data1$m2y,data1$m1y,cova)  # all variables to explain y
if(is.null(covariates))
  temp.data<-cbind(y=y,level=level,data1$x1,data1$m1y,data1$m2y)
else
  temp.data<-cbind(y=y,level=level,data1$x1,data1$m1y,data1$m2y, covariates)
if(surv)
  colnames(temp.data)<-c("time","status","level",colnames(data1$x1),colnames(data1$m1y),
                       colnames(data1$m2y), colnames(covariates))
else
  colnames(temp.data)<-c("y","level",colnames(data1$x1),colnames(data1$m1y),
                         colnames(data1$m2y), colnames(covariates))
if(levely==1)
 {frml<-getformula(expl,random,intercept)
  if(surv)
   {frml=gsub("y~","Surv(time, status)~",frml)
    f1<-coxme(as.formula(frml),data=data.frame(temp.data))}
  else if(biny & is.null(family1))
    f1<-glmer(frml,data=data.frame(temp.data),family=binomial(link="logit"))
  else if(!is.null(family1))
    f1<-glmer(frml,data=data.frame(temp.data),family=family1)
  else
    f1<-lmer(frml,data=data.frame(temp.data))}
else
  {temp.data2<-two(temp.data, level, weight)
   frml<-getformula(expl,random=NULL,intercept)
   if(surv)
   {frml=gsub("y","Surv(time, status)",frml)
    f1<-coxph(as.formula(frml),data=data.frame(temp.data2))}
   else if(biny & is.null(family1))
    f1<-glm(frml,data=data.frame(temp.data2),family=binomial(link="logit"))
   else if(!is.null(family1))
    f1<-glm(frml,data=data.frame(temp.data2),family=family1)
   else
    f1<-lm(frml,data=data.frame(temp.data2))}

lc1<-c(l1,c1)
lc2<-c(l2,c2)

if(surv)
  {coef.f1<-f1$coefficients
   cov.f1=vcov(f1)}
else if(intercept)
  {coef.f1<-summary(f1)$coefficient[-1,1] 
   cov.f1=vcov(f1)[-1,-1]}
else 
  {coef.f1<-summary(f1)$coefficient[,1]
   cov.f1=vcov(f1)}
cov.f1=as.matrix(cov.f1)
len<-c(ncol(data1$x1),ifelse(is.null(data1$m2y),0,ncol(data1$m2y)),
       ifelse(is.null(data1$m1y),0,ncol(data1$m1y)))

DE=matrix(0,n,ncol(x))   #calculate direct effects
DE.cov=DE                #calculate the variances of de
lx=data1$lx

for (i in 1:ncol(x))
 if (lx[i,2]<lx[i,3])                       
   {DE[,i]<-data2$x1.der[,lx[i,2]:lx[i,3]]%*%coef.f1[lx[i,2]:lx[i,3]]
    DE.cov[,i]=apply(data2$x1.der[,lx[i,2]:lx[i,3]],1,func1,mat=cov.f1[lx[i,2]:lx[i,3],lx[i,2]:lx[i,3]])}
 else 
   {DE[,i]<-data2$x1.der[,lx[i,2]]*coef.f1[lx[i,2]]
    DE.cov[,i]=data2$x1.der[,lx[i,2]:lx[i,3]]*cov.f1[lx[i,2]:lx[i,3],lx[i,2]:lx[i,3]]*data2$x1.der[,lx[i,2]:lx[i,3]]}

ie2_1<-NULL   #part of the IE for level-2 mediators
ie2_1.cov<-NULL #part of the variance for level 2 IE
if (len[2]>0)
 if (!is.null(data1$m2))
 {if(!is.null(l2))
   {for (k in 1:length(l2))
     if (length(data1$m2[[k+1]])==1)
      {ie2_1<-cbind(ie2_1,data2$m2y.der[,data1$m2[[k+1]]]*coef.f1[len[1]+data1$m2[[k+1]]])
       ie2_1.cov<-cbind(ie2_1.cov,data2$m2y.der[,data1$m2[[k+1]]]*cov.f1[len[1]+data1$m2[[k+1]],len[1]+data1$m2[[k+1]]]*data2$m2y.der[,data1$m2[[k+1]]])}
     else
      {ie2_1<-cbind(ie2_1,data2$m2y.der[,data1$m2[[k+1]]]%*%coef.f1[len[1]+data1$m2[[k+1]]])
       ie2_1.cov<-cbind(ie2_1.cov,apply(data2$m2y.der[,data1$m2[[k+1]]],1,func1,mat=cov.f1[len[1]+data1$m2[[k+1]],len[1]+data1$m2[[k+1]]]))}
    temp.name<-mnames[l2]}
   else {k=0
         temp.name=NULL}
   if(!is.null(c2))
     for (i in 1:length(c2))
     {if (length(data1$m2[[k+i+1]])==1)
       {ie2_1<-cbind(ie2_1,data2$m2y.der[,data1$m2[[k+i+1]]]*coef.f1[len[1]+data1$m2[[k+i+1]]])
        ie2_1.cov<-cbind(ie2_1.cov,data2$m2y.der[,data1$m2[[k+i+1]]]*cov.f1[len[1]+data1$m2[[k+i+1]],len[1]+data1$m2[[k+i+1]]]*data2$m2y.der[,data1$m2[[k+i+1]]])}
     else
       {ie2_1<-cbind(ie2_1,data2$m2y.der[,data1$m2[[k+i+1]]]%*%diag(coef.f1[len[1]+data1$m2[[k+i+1]]]))
        ie2_1.cov<-cbind(ie2_1.cov,apply(data2$m2y.der[,data1$m2[[k+i+1]]],1,func2,mat=diag(diag(cov.f1[len[1]+data1$m2[[k+i+1]],len[1]+data1$m2[[k+i+1]]]))))}
     temp.name=c(temp.name,paste(mnames[c2[i]],1:length(data1$m2[[k+i+1]]),sep="."))}
     colnames(ie2_1)<-temp.name
     colnames(ie2_1.cov)<-temp.name}

ie1_1<-NULL   #part of the IE for level-1 mediators
ie1_1.cov<-NULL # part of the variance of ie1
if (len[3]>0)
  if (!is.null(data1$m1))
  {if(!is.null(l1))
    {for (k in 1:length(l1))
      if (length(data1$m1[[k+1]])==1)
        {ie1_1<-cbind(ie1_1,data2$m1y.der[,data1$m1[[k+1]]]*coef.f1[len[1]+len[2]+data1$m1[[k+1]]])
         ie1_1.cov<-cbind(ie1_1.cov,(data2$m1y.der[,data1$m1[[k+1]]])^2*cov.f1[len[1]+len[2]+data1$m1[[k+1]],len[1]+len[2]+data1$m1[[k+1]]])}
      else
        {ie1_1<-cbind(ie1_1,data2$m1y.der[,data1$m1[[k+1]]]%*%coef.f1[len[1]+len[2]+data1$m1[[k+1]]])
         ie1_1.cov<-cbind(ie1_1.cov,apply(data2$m1y.der[,data1$m1[[k+1]]],1,func1,mat=cov.f1[len[1]+len[2]+data1$m1[[k+1]],len[1]+len[2]+data1$m1[[k+1]]]))}
      temp.name=mnames[l1]}
   else {k=0
         temp.name=NULL}
   if(!is.null(c1))
     for (i in 1:length(c1))
       {if (length(data1$m1[[k+i+1]])==1)
         {ie1_1<-cbind(ie1_1,data2$m1y.der[,data1$m1[[k+i+1]]]*coef.f1[len[1]+len[2]+data1$m1[[k+i+1]]])
          ie1_1.cov<-cbind(ie1_1.cov,(data2$m1y.der[,data1$m1[[k+i+1]]])^2*cov.f1[len[1]+len[2]+data1$m1[[k+i+1]],len[1]+len[2]+data1$m1[[k+i+1]]])}
        else
         {ie1_1<-cbind(ie1_1,data2$m1y.der[,data1$m1[[k+i+1]]]%*%diag(coef.f1[len[1]+len[2]+data1$m1[[k+i+1]]]))
          ie1_1.cov<-cbind(ie1_1.cov,t(apply(data2$m1y.der[,data1$m1[[k+i+1]]],1,func2,mat=diag(diag(cov.f1[len[1]+len[2]+data1$m1[[k+i+1]],len[1]+len[2]+data1$m1[[k+i+1]]])))))}
        temp.name=c(temp.name,paste(mnames[c1[i]],1:length(data1$m1[[k+i+1]]),sep=""))}
    colnames(ie1_1)<-temp.name
    colnames(ie1_1.cov)<-temp.name
  }
levelx=lx[,1]
nx=length(levelx)
nx1=sum(levelx==1) #number of level 1 x
nx2=sum(levelx==2)

m.2=data1$m.2
if(!is.null(m.2))    #analysis for level 2 mediators
  {fm2<-list(NULL)            #models for x to explain level 2 mediators
   coef.fm2<-list(NULL)
   if(!is.null(covariates))
    {cov.2<-two(covariates,level,weight)
     cov.2.new=two(covariates.new, level.new, weight.new)}
   else
    {cov.2<-NULL
     cov.2.new<-NULL}
  if(!is.null(l2))
  {fm2[[1]]<-l2
   coef.fm2[[1]]<-l2
   ie2_2<-array(NA,c(nrow(m.2),length(l2),nx2))
   colnames(ie2_2)=mnames[l2]
   dimnames(ie2_2)[[3]]=names(x)[levelx==2]
   ie2_2.cov<-ie2_2
   ie2_3=matrix(NA,nrow(m.2),length(l2))
   for (i in 1:length(l2))
   {if(sum(cm[[1]]==l2[i])>0)
    {temp2<-match(l2[i],cm[[1]])+1
     temp.cov<-as.matrix(cov.2[,cm[[temp2]]])
     colnames(temp.cov)<-colnames(covariates)[cm[[temp2]]]
     temp.cov<-as.matrix(cov.2.new[,cm[[temp2]]])
     colnames(temp.cov.new)<-colnames(covariates)[cm[[temp2]]]}
    else
     {temp.cov<-NULL
      temp.cov.new<-NULL}
   expl.m<-as.matrix(data1$xm2[,data1$fm22[[i+1]]])
   colnames(expl.m)<-colnames(data1$xm2)[data1$fm22[[i+1]]]
   expl.m.new<-as.matrix(data2$xm2[,data1$fm22[[i+1]]])
   colnames(expl.m.new)<-colnames(data1$xm2)[data1$fm22[[i+1]]]
   numx<-ncol(expl.m)
   expl.m<-cbind(expl.m,temp.cov)
   expl.m.new<-cbind(expl.m.new,temp.cov.new)
   frml.m<-getformula(expl.m,random=NULL,intercept)
   temp.data<-cbind(y=data1$m.2[,i],expl.m)
   colnames(temp.data)<-c("y",colnames(expl.m))
   temp.data.new<-cbind(y=data2$m.2[,i],expl.m.new)
   colnames(temp.data.new)<-c("y",colnames(expl.m))
   if(is.null(familym[[l2[i]]]))
      model<-lm(frml.m,data=data.frame(temp.data))
   else
      model<-glm(frml.m,data=data.frame(temp.data),family=familym[[l2[i]]])
   fm2<-append(fm2,list(model))
   coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
   cov.temp<-vcov(model)
   if(intercept)
     {coef.temp<-coef.temp[-1]
      cov.temp=as.matrix(cov.temp[-1,-1])}
   coef.temp<-coef.temp[1:numx]
   cov.temp=as.matrix(cov.temp[1:numx,1:numx])
   coef.fm2<-append(coef.fm2, list(coef.temp))
   temp.trans=NULL
   if(!is.null(data1$f01km2.2))
    temp.trans=(1:nrow(data1$f01km2.2[[1]]))[data1$f01km2.2[[1]][,1]==l2[i]]
   n2x=(1:nx)[levelx==2]
   for (j in 1:length(n2x)){
     if(!data1$binx[n2x[j]]){
      if(length(temp.trans)==0)
        {ie2_2[,i,j]=coef.temp[match(j,data1$fm22[[i+1]])] #if there is no transformation, the coefficient is the ie2_2
         ie2_2.cov[,i,j]=cov.temp[match(j,data1$fm22[[i+1]]),match(j,data1$fm22[[i+1]])]}
      else {
        temp.exp=data1$f01km2.2[[1]][temp.trans,2]
        temp.trans2=match(n2x[j],temp.exp)
      if(is.na(temp.trans2))
        {ie2_2[,i,j]=coef.temp[match(j,data1$fm22[[i+1]])] #if there is no transformation of the exposure variable
         ie2_2.cov[,i,j]=cov.temp[match(j,data1$fm22[[i+1]]),match(j,data1$fm22[[i+1]])]}
      else{
       temp.row=temp.trans[temp.trans2]
       if(length(data1$f01km2.2[[temp.row+1]])==1)
          {ie2_2[,i,j]<-coef.temp[match(data1$f01km2.2[[temp.row+1]],data1$fm22[[i+1]])]*
                        data2$xm2.der[,data1$f01km2.2[[temp.row+1]]]
           ie2_2.cov[,i,j]<-cov.temp[match(data1$f01km2.2[[temp.row+1]],data1$fm22[[i+1]]),match(data1$f01km2.2[[temp.row+1]],data1$fm22[[i+1]])]*
                            data2$xm2.der[,data1$f01km2.2[[temp.row+1]]]^2}
       else
          {ie2_2[,i,j]<-data2$xm2.der[,data1$f01km2.2[[temp.row+1]]]%*%
                       coef.temp[match(data1$f01km2.2[[temp.row+1]],data1$fm22[[i+1]])]
           ie2_2.cov[,i,j]<-apply(data2$xm2.der[,data1$f01km2.2[[temp.row+1]]],1,func1,
                                  mat=cov.temp[match(data1$f01km2.2[[temp.row+1]],data1$fm22[[i+1]]),match(data1$f01km2.2[[temp.row+1]],data1$fm22[[i+1]])])}
     }}}
    else{ #for binary x, use ie2_2 only, not that the mean is weighted by sample size
      if (length(data1$m2[[i+1]])==1)
        {ie2_2[,i,j]<-(mean(data2$m2y[x[,n2x[j]]==1,data1$m2[[i+1]]],na.rm=T)-mean(data2$m2y[x[,n2x[j]]==0,data1$m2[[i+1]]],na.rm=T))*
                       coef.f1[len[1]+data1$m2[[i+1]]]
         ie2_2.cov[,i,j]<-(mean(data2$m2y[x[,n2x[j]]==1,data1$m2[[i+1]]],na.rm=T)-mean(data2$m2y[x[,n2x[j]]==0,data1$m2[[i+1]]],na.rm=T))^2*
                           cov.f1[len[1]+data1$m2[[i+1]],len[1]+data1$m2[[i+1]]]}
      else
        {ie2_2[,i,j]<-(apply(data2$m2y[x[,n2x[j]]==1,data1$m2[[i+1]]],2,mean,na.rm=T)-apply(data2$m2y[x[,n2x[j]]==0,data1$m2[[i+1]]],2,mean,na.rm=T))%*%
                       coef.f1[len[1]+data1$m2[[i+1]]]
         vec.temp=apply(data2$m2y[x[,n2x[j]]==1,data1$m2[[i+1]]],2,mean,na.rm=T)-apply(data2$m2y[x[,n2x[j]]==0,data1$m2[[i+1]]],2,mean,na.rm=T)
         ie2_2.cov[,i,j]<-vec.temp%*%cov.f1[len[1]+data1$m2[[i+1]],len[1]+data1$m2[[i+1]]]%*%vec.temp}
    } 
   }
   ie2_3[,i]=d_link(link.fun=model$family$link,x=predict(model,newdata=data.frame(temp.data.new),type="response"))
   }
  j<-i
 }
 else
   {j<-0
   ie2_2<-NULL
   ie2_2.cov<-NULL
   ie2_3=NULL
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
      temp.cov.new=NULL}
     expl.m<-as.matrix(data1$xm2[,data1$fm22[[j+i+1]]])
     colnames(expl.m)<-colnames(data1$xm2)[data1$fm22[[j+i+1]]]
     expl.m.new<-as.matrix(data2$xm2[,data1$fm22[[j+i+1]]])
     colnames(expl.m.new)<-colnames(data1$xm2)[data1$fm22[[j+i+1]]]
     numx<-ncol(expl.m)
     expl.m<-cbind(expl.m,temp.cov)
     expl.m.new<-cbind(expl.m.new,temp.cov.new)
     frml.m<-getformula(expl.m,random=NULL,intercept)
    name.temp<-colnames(ie2_2)
    temp<-(1:length(data1$m2[[1]]))[data1$m2[[1]]==c2[i]]
    if(length(data1$m2[[temp+1]])>1)  #### for multicategorical level 2 mediators
    {temp.3<-apply(data1$m2y[,data1$m2[[temp+1]]],2,two,level)
     temp.4<-apply(temp.3,1,sum)
     temp.5<-(temp.4==0)}             #the reference group is temp.5
    else temp.5<-rep(T,length(unique(level)))  #for binary mediator, all observations
    for (k in data1$m2[[temp+1]])   #for each binarized categories
    {temp.6<-(temp.5 | (two(data1$m2y[,k],level)==1))
     temp.data<-cbind(y=two(data1$m2y[temp.6,k],level=level[temp.6]),expl.m[temp.6,]) #
     colnames(temp.data)<-c("y",colnames(expl.m)) #"level",
     if(is.null(familym[[c2[i]]]))
       model<-glm(frml.m,data=data.frame(temp.data),
                family=binomial(link="logit"))
     else
       model<-glm(frml.m,data=data.frame(temp.data),
                  family=familym[[c2[i]]])
     temp.data.new<-expl.m.new
     #p.temp<-predict(model,type="response",newdata=data.frame(temp.data))
     fm2<-append(fm2,list(model))
     fm2[[1]]<-c(fm2[[1]],c2[i])
     coef.fm2[[1]]<-c(coef.fm2[[1]],c2[i])
     coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
     cov.temp<-vcov(model)
     if(intercept)
       {coef.temp<-coef.temp[-1]
        cov.temp=as.matrix(cov.temp[-1,-1])}
     coef.temp<-coef.temp[1:numx]
     cov.temp<-as.matrix(cov.temp[1:numx,1:numx])
     coef.fm2=append(coef.fm2,list(coef.temp))
     temp.trans=match(c2[i],data1$f01km2.2[[1]][,1]) 
     n2x=(1:nx)[levelx==2]
     tmp.ie2_2<-NULL
     tmp.ie2_2.cov<-NULL
     
     for (z in 1:length(n2x)){    
       if(!data1$binx[n2x[z]]){
         if(is.na(temp.trans))
           {tmp.ie2_2=array(c(tmp.ie2_2,cbind(ie2_2[,,z],coef.temp[z])),
                            c(dim(cbind(ie2_2[,,z],coef.temp[z])),z))#if there is no transformation, the coefficient is the ie2_2
            tmp.ie2_2.cov=array(c(tmp.ie2_2.cov,cbind(ie2_2.cov[,,z],cov.temp[z,z])),
                                c(dim(cbind(ie2_2.cov[,,z],cov.temp[z,z])),z))}  #if there is no transformation, the variance of the coefficient is the ie2_2.cov
         else {temp.exp=data1$f01km2.2[[1]][temp.trans,2]
           temp.trans2=match(n2x[z],temp.exp)
         if(is.na(temp.trans2))
           {tmp.ie2_2=array(c(tmp.ie2_2,cbind(ie2_2[,,z],coef.temp[z])),
                            c(dim(cbind(ie2_2[,,z],coef.temp[z])),z)) #if there is no transformation of the exposure variable
            tmp.ie2_2.cov=array(c(tmp.ie2_2.cov,cbind(ie2_2.cov[,,z],cov.temp[z,z])),
                                c(dim(cbind(ie2_2.cov[,,z],cov.temp[z,z])),z)) #if there is no transformation of the exposure variable
           }
         else{temp.row=temp.trans[temp.trans2]
         if(length(data1$f01km2.2[[temp.trans2+1]])==1)
           {tmp.ie2_2<-array(c(tmp.ie2_2,cbind(ie2_2[,,z],coef.temp[match(data1$f01km2.2[[temp.row+1]],data1$fm22[[i+1]])]*
                                                data2$xm2.der[,data1$f01km2.2[[temp.row+1]]])),
                             c(dim(cbind(ie2_2[,,z],coef.temp[match(data1$f01km2.2[[temp.row+1]],data1$fm22[[i+1]])]*
                                          data2$xm2.der[,data1$f01km2.2[[temp.row+1]]])),z))
            tmp.ie2_2.cov<-array(c(tmp.ie2_2.cov,cbind(ie2_2.cov[,,z],cov.temp[match(data1$f01km2.2[[temp.row+1]],data1$fm22[[i+1]]),match(data1$f01km2.2[[temp.row+1]],data1$fm22[[i+1]])]*
                                                data2$xm2.der[,data1$f01km2.2[[temp.row+1]]]^2)),
                            c(dim(cbind(ie2_2.cov[,,z],cov.temp[match(data1$f01km2.2[[temp.row+1]],data1$fm22[[i+1]]),match(data1$f01km2.2[[temp.row+1]],data1$fm22[[i+1]])]*
                                          data2$xm2.der[,data1$f01km2.2[[temp.row+1]]]^2)),z))}
          else
            {tmp.ie2_2<-array(c(tmp.ie2_2,cbind(ie2_2[,,z],data2$xm2.der[,data1$f01km2.2[[temp.row+1]]]%*%
                                                  coef.temp[match(data1$f01km2.2[[temp.row+1]],data1$fm22[[i+1]])])),
                              c(dim(cbind(ie2_2[,,z],data2$xm2.der[,data1$f01km2.2[[temp.row+1]]]%*%
                                            coef.temp[match(data1$f01km2.2[[temp.row+1]],data1$fm22[[i+1]])])),z))
             temp1=apply(data2$xm2.der[,data1$f01km2.2[[temp.row+1]]],1,func1,
                         cov.temp[match(data1$f01km2.2[[temp.row+1]],data1$fm22[[i+1]]),match(data1$f01km2.2[[temp.row+1]],data1$fm22[[i+1]])])
             tmp.ie2_2.cov<-array(c(tmp.ie2_2.cov,cbind(ie2_2.cov[,,z],temp1)),
                             c(dim(cbind(ie2_2.cov[,,z],temp1)),z))}
         }}
       }
      else { #for binary x, use ie2_2 only, not that the mean is weighted by sample size
        if (length(data1$m2[[j+i+1]])==1)
          {tmp.ie2_2<-array(c(tmp.ie2_2,cbind(ie2_2[,,z],(mean(data2$m2y[x[,z]==1,data1$m2[[j+i+1]]],na.rm=T)-mean(data2$m2y[x[,z]==0,data1$m2[[j+i+1]]],na.rm=T))*
                             coef.f1[len[1]+data1$m2[[i+j+1]]])),
                           c(dim(cbind(ie2_2[,,z],(mean(data2$m2y[x[,z]==1,data1$m2[[j+i+1]]],na.rm=T)-mean(data2$m2y[x[,z]==0,data1$m2[[j+i+1]]],na.rm=T))*
                                         coef.f1[len[1]+data1$m2[[i+j+1]]])),z))
           temp1=(mean(data2$m2y[x[,z]==1,data1$m2[[j+i+1]]],na.rm=T)-mean(data2$m2y[x[,z]==0,data1$m2[[j+i+1]]],na.rm=T))^2*
                  cov.f1[len[1]+data1$m2[[i+j+1]],len[1]+data1$m2[[i+j+1]]]
           tmp.ie2_2.cov<-array(c(tmp.ie2_2.cov,cbind(ie2_2.cov[,,z],temp1)),
                           c(dim(cbind(ie2_2.cov[,,z],temp1)),z))}
        else
          {tmp.ie2_2<-array(c(tmp.ie2_2,cbind(ie2_2[,,z],(apply(data2$m2y[x[,z]==1,data1$m2[[j+i+1]]],2,mean,na.rm=T)-apply(data2$m2y[x[,z]==0,data1$m2[[j+i+1]]],2,mean,na.rm=T))%*%
                             coef.f1[len[1]+data1$m2[[i+j+1]]])),
                           c(dim(cbind(ie2_2[,,z],(apply(data2$m2y[x[,z]==1,data1$m2[[j+i+1]]],2,mean,na.rm=T)-apply(data2$m2y[x[,z]==0,data1$m2[[j+i+1]]],2,mean,na.rm=T))%*%
                                         coef.f1[len[1]+data1$m2[[i+j+1]]])),z))
           temp1=(apply(data2$m2y[x[,z]==1,data1$m2[[j+i+1]]],2,mean,na.rm=T)-apply(data2$m2y[x[,z]==0,data1$m2[[j+i+1]]],2,mean,na.rm=T))%*%
                  cov.f1[len[1]+data1$m2[[i+j+1]],len[1]+data1$m2[[i+j+1]]]%*%(apply(data2$m2y[x[,z]==1,data1$m2[[j+i+1]]],2,mean,na.rm=T)-apply(data2$m2y[x[,z]==0,data1$m2[[j+i+1]]],2,mean,na.rm=T))
           tmp.ie2_2.cov<-array(c(tmp.ie2_2.cov,cbind(ie2_2.cov[,,z],temp1)),
                                c(dim(cbind(ie2_2.cov[,,z],temp1)),z))}
        
      }   
     }   
     ie2_2=tmp.ie2_2
     ie2_2.cov=tmp.ie2_2.cov
     ie2_3=cbind(ie2_3,d_link(link.fun=model$family$link,x=predict(model,type="response",newdata=data.frame(temp.data.new))))
    }
    colnames(ie2_2)<-c(name.temp,paste(colnames(m)[c2[i]],
                                       1:length(data1$m2[[temp+1]]),sep="."))
    colnames(ie2_2.cov)<-c(name.temp,paste(colnames(m)[c2[i]],
                                       1:length(data1$m2[[temp+1]]),sep="."))
    colnames(ie2_3)<-c(name.temp,paste(colnames(m)[c2[i]],
                                       1:length(data1$m2[[temp+1]]),sep="."))
   }
  } 
else
{ie2_2<-NULL
 ie2_2.cov=NULL
 ie2_3=NULL
 fm2=NULL
 coef.fm2=NULL
}


fm1<-list(NULL)            #models for x to explain level 1 mediators
coef.fm1=list(NULL)
if(!is.null(c(l1,c1)))   #analysis for level 1 mediators
{
# ie1_2<-NULL
 if(!is.null(l1))
 {fm1[[1]]<-l1
  coef.fm1[[1]]=l1
  ie1_2=array(NA,c(n,length(l1),nx))
  colnames(ie1_2)=mnames[l1]
  if(nx1!=0)
  dimnames(ie1_2)[[3]]=x.names#[levelx==1]
  ie1_2.cov=ie1_2
  ie1_3=matrix(NA,n,length(l1))
  for (i in 1:length(l1))
  {if(sum(cm[[1]]==l1[i])>0) #add covariates
   {temp2<-1+match(l1[i],cm[[1]])
    temp.cov<-as.matrix(covariates[,cm[[temp2]]])
    colnames(temp.cov)<-colnames(covariates)[cm[[temp2]]]
    temp.cov.new<-as.matrix(covariates.new[,cm[[temp2]]])
    colnames(temp.cov.new)<-colnames(covariates)[cm[[temp2]]]}
   else
    {temp.cov<-NULL
     temp.cov.new<-NULL}
    expl.m<-as.matrix(data1$xm1[,c(data1$fm12[[i+1]],data1$fm11[[i+1]])])
    colnames(expl.m)<-colnames(data1$xm1)[c(data1$fm12[[i+1]],data1$fm11[[i+1]])]
    expl.m.new<-as.matrix(data2$xm1[,c(data1$fm12[[i+1]],data1$fm11[[i+1]])])
    colnames(expl.m.new)<-colnames(data1$xm1)[c(data1$fm12[[i+1]],data1$fm11[[i+1]])]
    numx<-ncol(expl.m)
    expl.m<-cbind(expl.m,temp.cov)
    expl.m.new<-cbind(expl.m.new,temp.cov.new)
    m.random<-"(1|level)"
   if(!is.null(random.m1))
     if(sum(random.m1[[1]]==l1[i])>0)
       m.random<-random.m1[[2]][random.m1[[1]]==l1[i]] #1st item of random.m1 is the list of l1 med, 2nd item is the random item of the same order
   frml.m<-getformula(expl.m,m.random,intercept)
   temp.data<-cbind(y=m[,l1[i]],level=level,expl.m)
   colnames(temp.data)<-c("y","level",colnames(expl.m))
   temp.data.new<-cbind(y=m.new[,l1[i]],level=level.new,expl.m.new)
   colnames(temp.data.new)<-c("y","level",colnames(expl.m))
   if(is.null(familym[[l1[i]]]))
      model=lmer(frml.m,data=data.frame(temp.data))
   else
      model=glmer(frml.m,data=data.frame(temp.data),family=familym[[l1[i]]])
   
   fm1<-append(fm1,model)
   coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
   cov.temp<-vcov(model)   #find the second part of IE2 var
   if(intercept)
     {coef.temp<-coef.temp[-1]
      cov.temp=as.matrix(cov.temp[-1,-1])}
   coef.temp<-coef.temp[1:numx]
   cov.temp=as.matrix(cov.temp[1:numx,1:numx])
   coef.fm1<-append(coef.fm1,list(coef.temp))
   if(!is.null(data1$fm12[[i+1]]))
     {coef.temp12<-coef.temp[1:length(data1$fm12[[i+1]])]
      coef.temp1<-coef.temp[-(1:length(data1$fm12[[i+1]]))]
      cov.temp12<-as.matrix(cov.temp[1:length(data1$fm12[[i+1]]),1:length(data1$fm12[[i+1]])])
      cov.temp1<-as.matrix(cov.temp[-(1:length(data1$fm12[[i+1]])),-(1:length(data1$fm12[[i+1]]))])}
   else
     {coef.temp1=coef.temp
      cov.temp1=cov.temp}
   
   temp.trans2=ifelse(is.null(data1$f01km1.2[[1]]),NA,(1:nrow(data1$f01km1.2[[1]]))[data1$f01km1.2[[1]][,1]==l1[i]])
   temp.trans1=ifelse(is.null(data1$f10km.2[[1]]),NA,(1:nrow(data1$f10km.2[[1]]))[data1$f10km.2[[1]][,1]==l1[i]]) 
   for (j in 1:nx){
   if(!data1$binx[j]){
     if(is.na(temp.trans2) & is.na(temp.trans1))
       {ie1_2[,i,j]=coef.temp[match(j,c(data1$fm12[[i+1]],data1$fm11[[i+1]]))] #if there is no transformation, the coefficient is the ie2_2
        ie1_2.cov[,i,j]=cov.temp[match(j,c(data1$fm12[[i+1]],data1$fm11[[i+1]])),match(j,c(data1$fm12[[i+1]],data1$fm11[[i+1]]))] }
     else if(!is.na(temp.trans1)){
       temp.exp=data1$f10km.2[[1]][temp.trans1,2]
       temp.trans1.2=match(j,temp.exp)
       if(is.na(temp.trans1.2))
         {ie1_2[,i,j]=coef.temp[match(j,c(data1$fm12[[i+1]],data1$fm11[[i+1]]))] #if there is no transformation of the exposure variable
          ie1_2.cov[,i,j]=cov.temp[match(j,c(data1$fm12[[i+1]],data1$fm11[[i+1]])),match(j,c(data1$fm12[[i+1]],data1$fm11[[i+1]]))]}
       else{
         temp.row=temp.trans1[temp.trans1.2]
         if(length(data1$f10km.2[[temp.row+1]])==1)
           {ie1_2[,i,j]<-coef.temp1[match(data1$f10km.2[[temp.row+1]],data1$fm11[[i+1]])]*
                         data2$xm1.der[,data1$f10km.2[[temp.row+1]]]
            ie1_2.cov[,i,j]<-cov.temp1[match(data1$f10km.2[[temp.row+1]],data1$fm11[[i+1]]),match(data1$f10km.2[[temp.row+1]],data1$fm11[[i+1]])]*
                             data2$xm1.der[,data1$f10km.2[[temp.row+1]]]^2}
         else
           {ie1_2[,i,j]<-data2$xm1.der[,data1$f10km.2[[temp.row+1]]]%*%
                         coef.temp1[match(data1$f10km.2[[temp.row+1]],data1$fm11[[i+1]])]
            ie1_2.cov[,i,j]<-apply(data2$xm1.der[,data1$f10km.2[[temp.row+1]]],1,func2,
                                   mat=cov.temp1[match(data1$f10km.2[[temp.row+1]],data1$fm11[[i+1]]),match(data1$f10km.2[[temp.row+1]],data1$fm11[[i+1]])])}
       }}
     else {
       temp.exp=data1$f01km1.2[[1]][temp.trans2,2]
       temp.trans2.2=match(j,temp.exp)
       if(is.na(temp.trans2.2))
         {ie1_2[,i,j]=coef.temp[match(j,c(data1$fm12[[i+1]],data1$fm11[[i+1]]))] #if there is no transformation of the exposure variable
          ie1_2.cov[,i,j]=cov.temp[match(j,c(data1$fm12[[i+1]],data1$fm11[[i+1]])),match(j,c(data1$fm12[[i+1]],data1$fm11[[i+1]]))] }
       else{
         temp.row=temp.trans2[temp.trans2.2]
         if(length(data1$f01km1.2[[temp.row+1]])==1)
           {ie1_2[,i,j]<-coef.temp12[match(data1$f01km1.2[[temp.row+1]],data1$fm12[[i+1]])]*
                         data2$xm1.der[,data1$f01km1.2[[temp.row+1]]]
            ie1_2.cov[,i,j]<-cov.temp12[match(data1$f01km1.2[[temp.row+1]],data1$fm12[[i+1]]),match(data1$f01km1.2[[temp.row+1]],data1$fm12[[i+1]])]*
                             data2$xm1.der[,data1$f01km1.2[[temp.row+1]]]^2}
         else
           {ie1_2[,i,j]<-data2$xm1.der[,data1$f01km1.2[[temp.row+1]]]%*%
                         coef.temp12[match(data1$f01km1.2[[temp.row+1]],data1$fm12[[i+1]])]
            ie1_2.cov[,i,j]<-apply(data2$xm1.der[,data1$f01km1.2[[temp.row+1]]],1,func1,
                             mat=cov.temp12[match(data1$f01km1.2[[temp.row+1]],data1$fm12[[i+1]]),match(data1$f01km1.2[[temp.row+1]],data1$fm12[[i+1]])])}
       }}
     }
   else{ #for binary x, use ie1_2 only, note that the mean is weighted by sample size for level 2 exposure
     if (length(data1$m1[[i+1]])==1)
       {ie1_2[,i,j]<-(mean(data2$m1y[x[,j]==1,data1$m1[[i+1]]],na.rm=T)-mean(data2$m1y[x[,j]==0,data1$m1[[i+1]]],na.rm=T))*
                      coef.f1[len[1]+len[2]+data1$m1[[i+1]]]
        ie1_2.cov[,i,j]<-(mean(data2$m1y[x[,j]==1,data1$m1[[i+1]]],na.rm=T)-mean(data2$m1y[x[,j]==0,data1$m1[[i+1]]],na.rm=T))^2*
                         cov.f1[len[1]+len[2]+data1$m1[[i+1]],len[1]+len[2]+data1$m1[[i+1]]]}
     else
       {ie1_2[,i,j]<-(apply(data2$m1y[x[,j]==1,data1$m1[[i+1]]],2,mean,na.rm=T)-apply(data2$m1y[x[,j]==0,data1$m1[[i+1]]],2,mean,na.rm=T))%*%
         coef.f1[len[1]+len[2]+data1$m1[[i+1]]]
        temp.vec=(apply(data2$m1y[x[,j]==1,data1$m1[[i+1]]],2,mean,na.rm=T)-apply(data2$m1y[x[,j]==0,data1$m1[[i+1]]],2,mean,na.rm=T))
        ie1_2.cov[,i,j]<-temp.vec%*%cov.f1[len[1]+len[2]+data1$m1[[i+1]],len[1]+len[2]+data1$m1[[i+1]]]%*%temp.vec}
   } 
   } 
   #browser()
   #if(surv & levely==1)
   #   ie1_3[,i]=predict(model) #newdata is not supported yet
   #else 
     ie1_3[,i]=d_link(link.fun=(summary(model))$link,x=predict(model,newdata=data.frame(temp.data.new),type="response",allow.new.levels=T))
  }
  j<-i
 }
 else
   {j<-0
    ie1_2=NULL
    ie1_2.cov=NULL
    ie1_3=NULL}

  
 if(!is.null(c1))  
   for (i in 1:length(c1))
   {if(sum(cm[[1]]==c1[i])>0)   #add covariates
     {temp2<-(1:length(cm[[1]]))[cm[[1]]==c1[i]]+1
      temp.cov<-as.matrix(covariates[,cm[[temp2]]])
      colnames(temp.cov)<-colnames(covariates)[cm[[temp2]]]
      temp.cov.new<-as.matrix(covariates.new[,cm[[temp2]]])
      colnames(temp.cov.new)<-colnames(covariates)[cm[[temp2]]]}
     else
      {temp.cov<-NULL
       temp.cov.new<-NULL}
     if(is.null(data1$fm11))
       k1<-NULL
     else
       k1<-(1:length(data1$fm11[[1]]))[data1$fm11[[1]]==c1[i]]
     if(is.null(data1$fm12))
       k2<-NULL
     else
       k2<-(1:length(data1$fm12[[1]]))[data1$fm12[[1]]==c1[i]]
     expl.m<-as.matrix(data1$xm1[,c(data1$fm12[[k2+1]],data1$fm11[[k1+1]])])
     colnames(expl.m)<-colnames(data1$xm1)[c(data1$fm12[[k2+1]],data1$fm11[[k1+1]])]
     expl.m.new<-as.matrix(data2$xm1[,c(data1$fm12[[k2+1]],data1$fm11[[k1+1]])])
     colnames(expl.m.new)<-colnames(data1$xm1)[c(data1$fm12[[k2+1]],data1$fm11[[k1+1]])]
     numx<-ncol(expl.m)
     expl.m<-cbind(expl.m,temp.cov)
     expl.m.new<-cbind(expl.m.new,temp.cov.new)
     m.random<-"(1|level)"
    if(!is.null(random.m1))
      if(sum(random.m1[[1]]==c1[i])>0)
        m.random<-random.m1[[match(c1[i],random.m1[[1]])+1]] #1st item of random.m1 is the list of l1 med, following items are the random item of the same order
    frml.m<-getformula(expl.m,m.random,intercept)
    name.temp<-dimnames(ie1_2)[[2]]
    temp<-(1:length(data1$m1[[1]]))[data1$m1[[1]]==c1[i]]
    if(length(data1$m1[[temp+1]])>1)###
    {temp.4<-apply(data1$m1y[,data1$m1[[temp+1]]],1,sum)
     temp.5<-(temp.4==0)}
    else temp.5<-rep(T,n)   #for binary mediator, all observations
    for (k in data1$m1[[temp+1]])  #for each binarized category
    {temp.6<-(temp.5 | (data1$m1y[,k]==1)) #use only the group and the reference group
     temp.data<-cbind(y=data1$m1y[temp.6,k],level=level[temp.6],expl.m[temp.6,])  #[temp.6]
     colnames(temp.data)<-c("y","level", colnames(expl.m))
     if(is.null(familym[[c1[i]]]))
       model=glmer(frml.m,data=data.frame(temp.data),
                    family=binomial(link="logit"))
     else
       model=glmer(frml.m,data=data.frame(temp.data),
                   family=familym[[c1[i]]])
     fm1<-append(fm1,model)
     fm1[[1]]<-c(fm1[[1]],c1[i])
     coef.fm1[[1]]<-c(coef.fm1[[1]],c1[i])
     temp.data.new<-cbind(level=level.new,expl.m.new)  #[temp.6]
     colnames(temp.data.new)<-c("level",colnames(expl.m))
     #p.temp<-predict(model,type="response",newdata=data.frame(temp.data),allow.new.levels=T)
     coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
     cov.temp<-vcov(model)   #find the second part of IE2
     if(intercept)
        {coef.temp<-coef.temp[-1]
         cov.temp<-as.matrix(cov.temp[-1,-1])}
     coef.temp<-coef.temp[1:numx]
     cov.temp<-as.matrix(cov.temp[1:numx,1:numx])
     coef.fm1<-append(coef.fm1,list(coef.temp))
     coef.temp12<-coef.temp[1:length(data1$fm12[[k2+1]])]
     coef.temp1<-coef.temp[-(1:length(data1$fm12[[k2+1]]))]
     cov.temp12<-as.matrix(cov.temp[1:length(data1$fm12[[k2+1]]),1:length(data1$fm12[[k2+1]])])
     cov.temp1<-as.matrix(cov.temp[-(1:length(data1$fm12[[k2+1]])),-(1:length(data1$fm12[[k2+1]]))])
     
     temp.trans1=match(c1[i],data1$f01km1.2[[1]][,1]) 
     temp.trans2=match(c1[i],data1$f10km.2[[1]][,1]) 
     tmp.ie1_2=NULL
     tmp.ie1_2.cov=NULL
     for (z in 1:nx){
       if(!data1$binx[z]){
         if(is.na(temp.trans1) & is.na(temp.trans2))
           {tmp.ie1_2=array(c(tmp.ie1_2,cbind(ie1_2[,,z],coef.temp[match(z,c(data1$fm12[[k2+1]],data1$fm11[[k1+1]]))])),
                           c(dim(cbind(ie1_2[,,z],coef.temp[match(z,c(data1$fm12[[k2+1]],data1$fm11[[k1+1]]))])),z))#if there is no transformation, the coefficient is the ie2_2
            tmp.ie1_2.cov=array(c(tmp.ie1_2.cov,cbind(ie1_2.cov[,,z],cov.temp[match(z,c(data1$fm12[[k2+1]],data1$fm11[[k1+1]])),match(z,c(data1$fm12[[k2+1]],data1$fm11[[k1+1]]))])),
                           c(dim(cbind(ie1_2.cov[,,z],cov.temp[match(z,c(data1$fm12[[k2+1]],data1$fm11[[k1+1]])),match(z,c(data1$fm12[[k2+1]],data1$fm11[[k1+1]]))])),z))#if there is no transformation, the coefficient is the ie2_2
           }
         else if(!is.na(temp.trans1)) 
           {temp.exp=data1$f01km1.2[[1]][temp.trans1,2]
            temp.trans1.2=match(z,temp.exp)
            if(is.na(temp.trans1.2))
              {tmp.ie1_2=array(c(tmp.ie1_2,cbind(ie1_2[,,z],coef.temp[match(z,c(data1$fm12[[k2+1]],data1$fm11[[k1+1]]))])),
                              c(dim(cbind(ie1_2[,,z],coef.temp[match(z,c(data1$fm12[[k2+1]],data1$fm11[[k1+1]]))])),z))#if there is no transformation of the exposure variable
               tmp.ie1_2.cov=array(c(tmp.ie1_2.cov,cbind(ie1_2.cov[,,z],cov.temp[match(z,c(data1$fm12[[k2+1]],data1$fm11[[k1+1]])),match(z,c(data1$fm12[[k2+1]],data1$fm11[[k1+1]]))])),
                              c(dim(cbind(ie1_2.cov[,,z],cov.temp[match(z,c(data1$fm12[[k2+1]],data1$fm11[[k1+1]])),match(z,c(data1$fm12[[k2+1]],data1$fm11[[k1+1]]))])),z))}
            else{temp.row=temp.trans1[temp.trans1.2]
                 if(length(data1$f01km1.2[[temp.row+1]])==1)
                   {tmp.ie1_2<-array(c(tmp.ie1_2,cbind(ie1_2[,,z],coef.temp1[match(data1$f01km1.2[[temp.row+1]],data1$fm11[[k1+1]])]*
                                     data2$xm1.der[,data1$f01km1.2[[temp.row+1]]])),
                                    c(dim(cbind(ie1_2[,,z],coef.temp1[match(data1$f01km1.2[[temp.row+1]],data1$fm11[[k1+1]])]*
                                                  data2$xm1.der[,data1$f01km1.2[[temp.row+1]]])),z))
                   tmp.ie1_2.cov<-array(c(tmp.ie1_2.cov,cbind(ie1_2.cov[,,z],cov.temp1[match(data1$f01km1.2[[temp.row+1]],data1$fm11[[k1+1]]),match(data1$f01km1.2[[temp.row+1]],data1$fm11[[k1+1]])]*
                                                        data2$xm1.der[,data1$f01km1.2[[temp.row+1]]]^2)),
                                    c(dim(cbind(ie1_2.cov[,,z],cov.temp1[match(data1$f01km1.2[[temp.row+1]],data1$fm11[[k1+1]]),match(data1$f01km1.2[[temp.row+1]],data1$fm11[[k1+1]])]*
                                                  data2$xm1.der[,data1$f01km1.2[[temp.row+1]]]^2)),z))}
                 else
                   {tmp.ie1_2<-array(c(tmp.ie1_2,cbind(ie1_2[,,z],data2$xm1.der[,data1$f01km1.2[[temp.row+1]]]%*%
                                     coef.temp1[match(data1$f01km1.2[[temp.row+1]],data1$fm11[[k1+1]])])),
                                    c(dim(cbind(ie1_2[,,z],data2$xm1.der[,data1$f01km1.2[[temp.row+1]]]%*%
                                                  coef.temp1[match(data1$f01km1.2[[temp.row+1]],data1$fm11[[k1+1]])])),z))
                    temp.vec=apply(data2$xm1.der[,data1$f01km1.2[[temp.row+1]]],1,func1,
                                   mat=cov.temp1[match(data1$f01km1.2[[temp.row+1]],data1$fm11[[k1+1]]),match(data1$f01km1.2[[temp.row+1]],data1$fm11[[k1+1]])])
                    tmp.ie1_2.cov<-array(c(tmp.ie1_2.cov,cbind(ie1_2.cov[,,z],temp.vec)),
                                    c(dim(cbind(ie1_2.cov[,,z],temp.vec)),z))}
         }}
         else  
           {temp.exp=data1$f10km.2[[1]][temp.trans2,2]
            temp.trans2.2=match(z,temp.exp)
            if(is.na(temp.trans2.2))
              {tmp.ie1_2=array(c(tmp.ie1_2,cbind(ie1_2[,,z],coef.temp[match(c1[i],c(data1$fm12[[k2+1]],data1$fm11[[k1+1]]))])),
                              c(dim(cbind(ie1_2[,,z],coef.temp[match(c1[i],c(data1$fm12[[k2+1]],data1$fm11[[k1+1]]))])),z)) #if there is no transformation of the exposure variable
               tmp.ie1_2.cov=array(c(tmp.ie1_2.cov,cbind(ie1_2.cov[,,z],cov.temp[match(c1[i],c(data1$fm12[[k2+1]],data1$fm11[[k1+1]])),match(c1[i],c(data1$fm12[[k2+1]],data1$fm11[[k1+1]]))])),
                              c(dim(cbind(ie1_2.cov[,,z],cov.temp[match(c1[i],c(data1$fm12[[k2+1]],data1$fm11[[k1+1]])),match(c1[i],c(data1$fm12[[k2+1]],data1$fm11[[k1+1]]))])),z))}
            else{temp.row=temp.trans2[temp.trans2.2]
                 if(length(data1$f10km.2[[temp.row+1]])==1)
                   {tmp.ie1_2<-array(c(tmp.ie1_2,cbind(ie1_2[,,z],coef.temp12[match(data1$f10km.2[[temp.row+1]],data1$fm12[[k2+1]])]*
                                     data2$xm1.der[,data1$f10km.2[[temp.row+1]]])),
                                    c(dim(cbind(ie1_2[,,z],coef.temp12[match(data1$f10km.2[[temp.row+1]],data1$fm12[[k2+1]])]*
                                                  data2$xm1.der[,data1$f10km.2[[temp.row+1]]])),z))
                   tmp.ie1_2.cov<-array(c(tmp.ie1_2.cov,cbind(ie1_2.cov[,,z],cov.temp12[match(data1$f10km.2[[temp.row+1]],data1$fm12[[k2+1]]),match(data1$f10km.2[[temp.row+1]],data1$fm12[[k2+1]])]*
                                                        data2$xm1.der[,data1$f10km.2[[temp.row+1]]]^2)),
                                    c(dim(cbind(ie1_2.cov[,,z],cov.temp12[match(data1$f10km.2[[temp.row+1]],data1$fm12[[k2+1]]),match(data1$f10km.2[[temp.row+1]],data1$fm12[[k2+1]])]*
                                                  data2$xm1.der[,data1$f10km.2[[temp.row+1]]])),z))}
                 else
                   {tmp.ie1_2<-array(c(tmp.ie1_2,cbind(ie1_2[,,z],data2$xm1.der[,data1$f10km.2[[temp.row+1]]]%*%
                                     coef.temp12[match(data1$f10km.2[[temp.row+1]],data1$fm12[[k2+1]])])),
                                    c(dim(cbind(ie1_2[,,z],data2$xm1.der[,data1$f10km.2[[temp.row+1]]]%*%
                                                  coef.temp12[match(data1$f10km.2[[temp.row+1]],data1$fm12[[k2+1]])])),z))
                    temp.vec=apply(data2$xm1.der[,data1$f10km.2[[temp.row+1]]],1,func1,
                                   mat=cov.temp12[match(data1$f10km.2[[temp.row+1]],data1$fm12[[k2+1]]),match(data1$f10km.2[[temp.row+1]],data1$fm12[[k2+1]])])
                    tmp.ie1_2.cov<-array(c(tmp.ie1_2.cov,cbind(ie1_2.cov[,,z],temp.vec)),
                                    c(dim(cbind(ie1_2[,,z],temp.vec)),z))}
         }}
       }
       else { #for binary x, use ie2_2 only, not that the mean is weighted by sample size
         {tmp.ie1_2<-array(c(tmp.ie1_2,cbind(ie1_2[,,z],(mean(data2$m1y[x[,z]==1,k],na.rm=T)-
                                                          mean(data2$m1y[x[,z]==0,k],na.rm=T))*
                                              coef.f1[len[1]+len[2]+k])),
                          c(dim(cbind(ie1_2[,,z],(mean(data2$m1y[x[,z]==1,k],na.rm=T)-
                                                    mean(data2$m1y[x[,z]==0,k],na.rm=T))*
                                        coef.f1[len[1]+len[2]+k])),z))
         temp.vec=(mean(data2$m1y[x[,z]==1,k],na.rm=T)-mean(data2$m1y[x[,z]==0,k],na.rm=T))^2*cov.f1[len[1]+len[2]+k,len[1]+len[2]+k]
         tmp.ie1_2.cov<-array(c(tmp.ie1_2.cov,cbind(ie1_2.cov[,,z],temp.vec)),
                          c(dim(cbind(ie1_2.cov[,,z],temp.vec)),z))}
         
       }}
     ie1_2=tmp.ie1_2
     ie1_2.cov=tmp.ie1_2.cov
     ie1_3=cbind(ie1_3,d_link(link.fun=(summary(model))$link,x=predict(model,type="response",newdata=data.frame(temp.data.new),allow.new.levels=T)))
    }
    colnames(ie1_2)<-c(name.temp,paste(colnames(m)[c1[i]],
                                       1:length(data1$m1[[temp+1]]),sep="."))
    dimnames(ie1_2)[[3]]=x.names
    colnames(ie1_2.cov)<-c(name.temp,paste(colnames(m)[c1[i]],
                                       1:length(data1$m1[[temp+1]]),sep="."))
    dimnames(ie1_2.cov)[[3]]=x.names
    colnames(ie1_3)<-c(name.temp,paste(colnames(m)[c1[i]],
                                        1:length(data1$m1[[temp+1]]),sep="."))
   }
}


if(!is.null(ie1_1) & !is.null(ie1_2))   #get the estimates for level 1 mediators
 {if(nx1!=0)
   {ie1=array(ie1_2[,,levelx==1],c(n,dim(ie1_2)[2],nx1))
    colnames(ie1)=colnames(ie1_2)
    dimnames(ie1)[[3]]=x.names[levelx==1]
    ie1.cov=ie1
    aie1=matrix(NA,nx1,ncol(ie1))
    colnames(aie1)=colnames(ie1)
    rownames(aie1)=x.names[levelx==1]
    aie1.cov=aie1
    llx1=1}
  else
  {ie1=NULL
   aie1=NULL
   ie1.cov=ie1
   aie1.cov=aie1}
  if(nx2!=0)
  {ie12=array(NA,c(dim(two(ie1_2[,,1],level)),nx2))
   colnames(ie12)=colnames(ie1_2)
   ie12.cov=ie12
   aie12=matrix(NA,nx2,ncol(ie12))
   colnames(aie12)=colnames(ie12)
   rownames(aie12)=x.names[levelx==2]
   aie12.cov=aie12
   llx2=1}
  else
  {ie12=NULL
   aie12=NULL
   ie12.cov=NULL
   aie12.cov=NULL}
  for (i in 1:nx)
   {if(levelx[i]==1){
    if (data1$binx[i])
     {ie1[,,llx1]<-ie1_2[,,i]
      ie1.cov[,,llx1]<-ie1_2.cov[,,i]}
    else
     {ie1[,,llx1]<-ie1_1*ie1_2[,,i]*ie1_3
      ie1.cov[,,llx1]<-(ie1_2[,,i]*ie1_3)^2*ie1_1.cov+(ie1_1*ie1_3)^2*ie1_2.cov[,,i]}
    aie1[llx1,]=apply(as.matrix(ie1[,,llx1]),2,mean,na.rm=T)
    aie1.cov[llx1,]=(apply(as.matrix(sqrt(ie1.cov[,,llx1])),2,mean,na.rm=T))^2 #approximation, not exact
    llx1=llx1+1
   }
    else{ 
      ie1_2.temp=two(ie1_2[,,i],level)
      ie1_1.temp=two(ie1_1,level)
      ie1_3.temp=two(ie1_3,level)
      ie1_2.cov.temp=(two(sqrt(ie1_2.cov[,,i]),level))^2
      ie1_1.cov.temp=(two(sqrt(ie1_1.cov),level))^2
      if (data1$binx[i])
        {ie12[,,llx2]<-ie1_2.temp
         ie12.cov[,,llx2]<-ie1_2.cov.temp}
      else
        {ie12[,,llx2]<-ie1_1.temp*ie1_2.temp*ie1_3.temp
         ie12.cov[,,llx2]<-(ie1_1.temp*ie1_3.temp)^2*ie1_2.cov.temp+(ie1_3.temp*ie1_2.temp)^2*ie1_1.cov.temp}
      aie12[llx2,]=apply(as.matrix(ie12[,,llx2]),2,mean,na.rm=T)
      aie12.cov[llx2,]=(apply(as.matrix(sqrt(ie12.cov[,,llx2])),2,mean,na.rm=T))^2
      llx2=llx2+1
    }
}}
else {ie1<-NULL
      aie1=NULL
      ie12=NULL
      aie12=NULL
      ie1.cov=NULL
      aie1.cov=NULL
      ie12.cov=NULL
      aie12.cov=NULL}


if(!is.null(ie2_1) & !is.null(ie2_2))   #get the estimates for level 1 mediators
{ie2=array(0,c(dim(data2$m.2),nx2))
 dimnames(ie2)=dimnames(ie2_2)
 ie2.cov=ie2
 aie2=matrix(NA,nx2,ncol(ie2))
 colnames(aie2)=colnames(ie2)
 rownames(aie2)=x.names[levelx==2]
 aie2.cov=aie2
 numx2=(1:nx)[levelx==2]  #place of level2 exposure in x
 ie2_1.temp=two(ie2_1,level)
 ie2_1.cov.temp=(two(sqrt(ie2_1.cov),level))^2
 ie2_3.temp=ie2_3
 for (i in 1:nx2)
  {if(data1$binx[numx2[i]])
     {ie2[,,i]<-ie2_2[,,i]
      ie2.cov[,,i]<-ie2_2.cov[,,i]}
   else
     {ie2[,,i]<-ie2_1.temp*ie2_2[,,i]*ie2_3.temp
      ie2.cov[,,i]<-ie2_1.cov.temp*(ie2_2[,,i]*ie2_3.temp)^2+ie2_2.cov[,,i]*(ie2_1.temp*ie2_3.temp)^2}
   aie2[i,]=apply(as.matrix(ie2[,,i]),2,mean,na.rm=T)
   aie2.cov[i,]=(apply(as.matrix(sqrt(ie2.cov[,,i])),2,mean,na.rm=T))^2
 }
}
else {ie2<-NULL
aie2=NULL
ie2.cov<-NULL
aie2.cov=NULL
}


if(nx1>0)                         #calculate the direct effect
{de1=as.matrix(DE[,levelx==1])
colnames(de1)=x.names[levelx==1]
ade1=apply(de1,2,mean)
de1.cov=as.matrix(DE.cov[,levelx==1])
colnames(de1.cov)=x.names[levelx==1]
ade1.cov=(apply(sqrt(de1.cov),2,mean))^2
}
else 
  {de1=NULL
   ade1=NULL
   de1.cov=NULL
   ade1.cov=NULL}
   
if(nx2>0)
  {de2=two(DE[,levelx==2],level)
   ade2=apply(de2,2,mean)
   de2.cov=(two(sqrt(DE.cov[,levelx==2]),level))^2
   ade2.cov=(apply(sqrt(de2.cov),2,mean))^2}
else
  {de2=NULL
   ade2=NULL
   de2.cov=NULL
   ade2.cov=NULL}
   
if(nx1>0)                         #calculate the total effect
  {te1=as.matrix(de1+apply(ie1,c(1,3),sum,na.rm=T))
   colnames(te1)=x.names[levelx==1]
   ate1=apply(te1,2,mean)}
   else 
   {te1=NULL
    ate1=NULL}
   
if(nx2>0)
 {if(is.null(ie2))
   te2=as.matrix(de2+apply(ie12,c(1,3),sum,na.rm=T))
  else
   te2=as.matrix(de2+apply(ie2,c(1,3),sum,na.rm=T)+apply(ie12,c(1,3),sum,na.rm=T))
  colnames(te2)=x.names[levelx==2]
  ate2=apply(te2,2,mean)}
else
 {te2=NULL
  ate2=NULL}
   
if(!is.null(joint))  #get the estimates for joint mediator effect
{joint1<-(1:length(joint[[1]]))[joint[[1]]==1]
 joint2<-(1:length(joint[[1]]))[joint[[1]]==2]
 if(length(joint1)!=0)
 {je1=array(0,c(n,length(joint1),nx1))
  colnames(je1)=paste("j",joint1,sep="")
  dimnames(je1)[[3]]=x.names[levelx==1]
  for (i in length(joint1))
  {posi1=order_char(colnames(ie1),mnames[joint[joint1[i]+1]])
   je1[,i,]=apply(ie1[,posi1,],c(1,3),sum,na.rm=T)
  }
  aje1=t(apply(je1,c(2,3),mean,na.rm=T))
 }
 else
 {je1=NULL
  aje1=NULL}
 if(length(joint2)!=0)
 {je2=array(0,c(nrow(data2$m.2),length(joint2),nx2))
 colnames(je2)=paste("j",joint2,sep="")
 dimnames(je2)[[3]]=x.names[levelx==2]
 for (i in 1:length(joint2))
 {if (is.character(joint[[joint2[i]+1]]))
     temp.6=joint[[joint2[i]+1]]
   else
     temp.6=mnames[joint[[joint2[i]+1]]]
  posi1=order_char(colnames(ie12),temp.6)
  posi2=order_char(colnames(ie2),temp.6)
  temp.6=dim(ie2)
  temp.7=dim(ie12)
  if(length(posi1)==0)
    je2[,i,]=apply(array(ie2[,posi2,],c(temp.6[1],length(posi2),temp.6[3])),c(1,3),sum,na.rm=T)
  else if (length(posi2)==0)
    je2[,i,]=apply(array(ie12[,posi1,],c(temp.7[1],length(posi1),temp.7[3])),c(1,3),sum,na.rm=T)
  else
    je2[,i,]=apply(array(ie12[,posi1,],c(temp.7[1],length(posi1),temp.7[3])),c(1,3),sum,na.rm=T)+
             apply(array(ie2[,posi2,],c(temp.6[1],length(posi2),temp.6[3])),c(1,3),sum,na.rm=T)
  }
 aje2=t(apply(je2,c(2,3),mean,na.rm=T))
 }
 else
 {je2=NULL
  aje2<-NULL
 }}
else
{je1=NULL
 aje1=NULL
 je2=NULL
 aje2<-NULL}

if(cov.mat)
  a<-list(de1=de1,de2=de2,ade1=ade1,ade2=ade2,te1=te1,te2=te2,ate1=ate1,ate2=ate2,
          ie1=ie1,ie2=ie2,ie12=ie12, aie1=aie1,aie2=aie2,aie12=aie12,
          f1=f1,fm1=fm1,fm2=fm2, je1=je1,je2=je2,aje1=aje1,aje2=aje2,
          ie1_1=ie1_1,ie1_2=ie1_2,ie1_3=ie1_3,ie2_1=ie2_1,ie2_2=ie2_2, ie2_3=ie2_3,
          x=x.new, x.j=two(x.new,level.new), m=m, covariates=covariates, 
          intercept=intercept, cm=cm,data1=data1,data2=data2,surv=surv,
          coef.f1=coef.f1,coef.fm1=coef.fm1,coef.fm2=coef.fm2,
          de1.cov=de1.cov,de2.cov=de2.cov,ade1.cov=ade1.cov,ade2.cov=ade2.cov,
          ie1.cov=ie1.cov,ie2.cov=ie2.cov,ie12.cov=ie12.cov, 
          aie1.cov=aie1.cov,aie2.cov=aie2.cov,aie12.cov=aie12.cov)
else
  a<-list(de1=de1,de2=de2,ade1=ade1,ade2=ade2,te1=te1,te2=te2,ate1=ate1,ate2=ate2,
          ie1=ie1,ie2=ie2,ie12=ie12, aie1=aie1,aie2=aie2,aie12=aie12,
          f1=f1,fm1=fm1,fm2=fm2, je1=je1,je2=je2,aje1=aje1,aje2=aje2,
          ie1_1=ie1_1,ie1_2=ie1_2,ie1_3=ie1_3,ie2_1=ie2_1,ie2_2=ie2_2, ie2_3=ie2_3,
          x=x.new, x.j=two(x.new,level.new), m=m, covariates=covariates, 
          intercept=intercept, cm=cm,data1=data1,data2=data2,surv=surv,
          coef.f1=coef.f1,coef.fm1=coef.fm1,coef.fm2=coef.fm2)

class(a)<-"mlma"
return(a)
}

boot.mlma<-function(y, data1=NULL,x=data1$parameter$x, m=data1$parameter$m, levelx=data1$parameter$levelx, 
                    levely=data1$parameter$levely,xref=NULL, yref=NULL, l1=data1$parameter$l1,l2=data1$parameter$l2, 
                    c1=data1$parameter$c1, #levelx is the level of x
                    c1r=data1$parameter$c1r, c2=data1$parameter$c2, c2r=data1$parameter$c2r, 
                    level=data1$parameter$level,  weight=rep(1,nrow(as.matrix(x))), 
                    random="(1|level)", random.m1=NULL,intercept=TRUE, 
                    family1=NULL, familym=vector("list",ncol(m)),
                    covariates=NULL, cy1=NULL, cy2=NULL, cm=NULL,joint=NULL,f01y=data1$parameter$f01y,
                    f10y=data1$parameter$f10y, f02ky=data1$parameter$f02ky, 
                    f20ky=data1$parameter$f20ky, f01km1=data1$parameter$f01km1, 
                    f01km2=data1$parameter$f01km2, f10km=data1$parameter$f10km,
                    data2=NULL, x.new=NULL, m.new=m, level.new=level, weight.new=NULL,
                    covariates.new=covariates,boot=100, echo=TRUE, plot.it=TRUE, cov.mat=FALSE)
{two<-function(x, level, weight=rep(1,nrow(as.matrix(x))))
{x2<-NULL
if(is.null(dim(x)))
  x2=cbind(x2,(aggregate(as.numeric(x)~level, na.action=na.pass, 
                         FUN=function(x) c(mean=weighted.mean(x,weight=weight,na.rm=T))))[,2])
else
  for(i in 1:ncol(x))
    x2<-cbind(x2,(aggregate(as.numeric(x[,i])~level, na.action=na.pass,
                            FUN=function(x) c(mean=weighted.mean(x,weight=weight,na.rm=T))))[,2])
  colnames(x2)<-colnames(x)
  x2
}

update.data.org<-function(data1,bootsample,y.boot,x.boot,
                          m.boot,weight.boot,level.boot)  #data1 is the original data.org results, bootsample is booted sample
{x2fx<-function(x,func) #x is the list of original numerical vector, func is a vector of character functions. 
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
cbind(x1,x)
}

if (ncol(data.frame(x.boot))!=ncol(data1$parameter$x)) #if there are level 2 mediator(s), but no level 2 exposure(s)
{x.temp=NULL
x.names=colnames(x.boot)
for (i in 1:ncol(x.boot))
  x.temp=cbind(x.temp,one2two(x.boot[,i],level.boot))
x.temp=data.frame(x.temp)
colnames(x.temp)=colnames(data1$parameter$x)
x.boot=x.temp
}

data3=data1
data3$parameter$x=x.boot
data3$parameter$weight=weight.boot
data3$parameter$m=m.boot
data3$parameter$level=level.boot

if(!is.null(data1$x1))
{data3$x1=as.matrix(data1$x1[bootsample,])
colnames(data3$x1)=colnames(data1$x1)}
if(!is.null(data1$x1.der))
{data3$x1.der=as.matrix(data1$x1.der[bootsample,])
colnames(data3$x1.der)=colnames(data1$x1.der)}
if(!is.null(data1$m1y))
{data3$m1y=as.matrix(data1$m1y[bootsample,])
colnames(data3$m1y)=colnames(data1$m1y)}
if(!is.null(data1$m1y.der))
{data3$m1y.der=as.matrix(data1$m1y.der[bootsample,])
colnames(data3$m1y.der)=colnames(data1$m1y.der)}
if(!is.null(data1$m2y))
{data3$m2y=as.matrix(data1$m2y[bootsample,])
colnames(data3$m2y)=colnames(data1$m2y)}
if(!is.null(data1$m2y.der))
{data3$m2y.der=as.matrix(data1$m2y.der[bootsample,])
colnames(data3$m2y.der)=colnames(data1$m2y.der)}
if(!is.null(data1$xm1))
{data3$xm1=as.matrix(data1$xm1[bootsample,])
colnames(data3$xm1)=colnames(data1$xm1)}
if(!is.null(data1$xm1.der))
{data3$xm1.der=as.matrix(data1$xm1.der[bootsample,])
colnames(data3$xm1.der)=colnames(data1$xm1.der)}

lc2<-c(data1$parameter$l2,data1$parameter$c2)
if(!is.null(lc2))                  #create x variables to explain level 2 mediators
{data3$xm2<-as.matrix(two(as.matrix(x.boot[,data1$parameter$levelx==2]),level.boot))
#colnames(data3$xm2)<- colnames(data1$xm2)
data3$xm2.der<-matrix(1,nrow(data3$xm2),ncol(data3$xm2))
#colnames(xm2.der)<-colnames(xm2)

data3$m.2<-two(as.matrix(m.boot[,lc2]),level.boot)
colnames(data3$m.2)<-data1$parameter$mnames[lc2]

if(!is.null(data1$parameter$f01km2))
{k=unique(data1$parameter$f01km2[[1]][,2])
for (l in k){
  temp.2=as.matrix(two(as.matrix(x.boot[,k]),level.boot))
  temp<-(2:length(data1$parameter$f01km2))[data1$parameter$f01km2[[1]][,2]==l]
  allfun=data1$parameter$f01km2[[temp[1]]]
  if (length(temp)>1)
    for(i in 2:length(temp))
      allfun<-c(allfun,data1$parameter$f01km2[[temp[i]]])
  unifun<-unique(allfun)
  unifun1<-unifun[unifun!="x"]
  unifun2<-c("x",unifun1)
  d_d<-x2fx(temp.2,unifun1)
  d.der<-as.matrix(x2fdx(temp.2,unifun1))
  d<-as.matrix(d_d$values)
  data3$xm2.der<-cbind(data3$xm2.der,d.der)
  colnames(data3$xm2.der)<-colnames(data1$xm2.der)
  data3$xm2<-cbind(data3$xm2,d)
  colnames(data3$xm2)<-colnames(data1$xm2)
}}}
else
{data3$m.2<-NULL
data3$xm2<-NULL
data3$xm2.der<-NULL
}
data3
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
{a<-as.factor(m[,i])
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

find_level<-function(vari, level)    #a function to identify mediators to l1, l2, c1, c2
{a<-table(level)
b<-sort(unique(level))
l=2

for (i in 1:length(b))
  if(a[i]>1)
  {temp1<-vari[level==b[i]]
  temp2<-temp1[!is.na(temp1)]
  if(length(temp2)!=0 & sum(temp2!=temp2[1])>0)
  {l=1
  break}
  }
if(l==1 & (is.factor(vari) | is.character(vari) | nlevels(as.factor(vari))==2))
  return (list(1, levels(as.factor(vari))[1]))  #c1, c1r
else if(l==1)
  return (list(2,NA))   #l1
else if(l==2 & (is.factor(vari) | is.character(vari) | nlevels(as.factor(vari))==2))
  return (list(3, levels(as.factor(vari))[1]))  #c2, c2r
else
  return(list(4,NA))
}

generate.m<-function(x.new, data1,l1,c1,l2,c2,full,covariates.new,level.new,data2,cm)
{n.new<-nrow(x.new)
 n2.new<-length(unique(level.new))
 mnames=data1$parameter$mnames
 m.new<-matrix(0,n.new,length(mnames))
 colnames(m.new)<-mnames
 fm1<-full$fm1
 if(!is.null(l1))        #generate level 1 continuous mediators
  {for (i in 1:length(l1))
  {if(sum(cm[[1]]==l1[i])>0)
   {temp2<-match(l1[i],cm[[i]])+1
    temp.cov<-as.matrix(covariates.new[,cm[[temp2]]])
    colnames(temp.cov)<-colnames(covariates.new)[cm[[temp2]]]}
   else
    temp.cov<-NULL
    expl.m<-as.matrix(data2$xm1[,c(data1$fm12[[i+1]],data1$fm11[[i+1]])])
    colnames(expl.m)<-colnames(data2$xm1)[c(data1$fm12[[i+1]],data1$fm11[[i+1]])]
    expl.m<-cbind(expl.m,temp.cov)
    temp.data<-cbind(level=level.new,expl.m)
    colnames(temp.data)<-c("level",colnames(expl.m))
    j<-(1:length(fm1[[1]]))[fm1[[1]]==l1[i]]+1
    sd.m<-rnorm(n.new,mean=0,sd=as.data.frame(VarCorr(full$fm1[[j]]))[2,5])+         #add the randomn effects
      one(rnorm(n2.new,mean=0,sd=as.data.frame(VarCorr(full$fm1[[j]]))[1,5]),level.new)
    m.new[,l1[i]]<-predict(full$fm1[[j]],newdata=data.frame(temp.data),allow.new.levels=T)+sd.m
  }
  }
  if(!is.null(c1))      #generate level 1 categorical mediators
    for (i in 1:length(c1))
    {if(sum(cm[[1]]==c1[i])>0)
    {temp2<-(1:length(cm[[1]]))[cm[[1]]==c1[i]]+1  #may have more than one case?
    temp.cov<-as.matrix(covariates.new[,cm[[temp2]]])
    colnames(temp.cov)<-colnames(covariates.new)[cm[[temp2]]]}
      else
        temp.cov<-NULL
      if(is.null(data1$fm11))
        k1<-NULL
      else
        k1<-(1:length(data1$fm11[[1]]))[data1$fm11[[1]]==c1[i]]
      if(is.null(data1$fm12))
        k2<-NULL
      else
        k2<-(1:length(data1$fm12[[1]]))[data1$fm12[[1]]==c1[i]]
      expl.m<-as.matrix(data2$xm1[,c(data1$fm12[[k2+1]],data1$fm11[[k1+1]])])
      colnames(expl.m)<-colnames(data2$xm1)[c(data1$fm12[[k2+1]],data1$fm11[[k1+1]])]
      expl.m<-cbind(expl.m,temp.cov)
      temp.data<-cbind(level=level.new,expl.m)
      colnames(temp.data)<-c("level",colnames(expl.m))
      j<-(1:length(fm1[[1]]))[full$fm1[[1]]==c1[i]]+1
      if(length(j)==1)   #binary
      {tmp.logit<-predict(full$fm1[[j]],newdata=data.frame(temp.data),allow.new.levels=T)+
        one(rnorm(n2.new,mean=0,sd=as.data.frame(VarCorr(fm1[[j]]))[,5]),level.new)  #add the random effect
      tmp.p<-exp(tmp.logit)/(exp(tmp.logit)+1)
      tmp.m<-rbinom(n.new,1,prob=tmp.p)
      m.new[,c1[i]]<-ifelse(tmp.m==0,1,2)}
      else               #multi-categorical
      {tmp.logit<-rep(1,n.new)
      for (k in 1:length(j))
        tmp.logit<-cbind(tmp.logit,exp(predict(full$fm1[[j[k]]],newdata=data.frame(temp.data),allow.new.levels=T)+
                                         one(rnorm(n2.new,mean=0,sd=as.data.frame(VarCorr(fm1[[j[k]]]))[,5]),level.new))) #add the random effect
      tmp.m<-apply(tmp.logit,1,rmultinom,n=1,size=1)
      m.new[,c1[i]]<-(1:nrow(tmp.m))%*%tmp.m   #the reference group is always 1
      }
    }
  
  fm2<-full$fm2            #generate level 2 mediators
  cov.2<-NULL
  if(!is.null(covariates.new))
    cov.2<-two(covariates.new,level.new)
  if(!is.null(l2))          #generate level 2 mediators
    for (i in 1:length(l2))
    {if(sum(cm[[1]]==l2[i])>0)
    {temp2<-(1:length(cm[[1]]))[cm[[1]]==l2[i]]+1
    temp.cov<-as.matrix(cov.2[,cm[[temp2]]])
    colnames(temp.cov)<-colnames(covariates.new)[cm[[temp2]]]}
      else
        temp.cov<-NULL
      expl.m<-as.matrix(data2$xm2[,data1$fm22[[i+1]]])
      colnames(expl.m)<-colnames(data1$xm2)[data1$fm22[[i+1]]]
      expl.m<-cbind(expl.m,temp.cov)
      j<-(1:length(fm2[[1]]))[fm2[[1]]==l2[i]]+1
      tmp.m<-predict(fm2[[j]],newdata=data.frame(expl.m))+
        rnorm(n2.new,mean=0,sd=(summary(fm2[[j]]))$sigma)  #add the randomness
      m.new[,l2[i]]<-one(tmp.m,level.new) #lm was used
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
      expl.m<-as.matrix(data2$xm2[,data1$fm22[[k+1]]])
      colnames(expl.m)<-colnames(data1$xm2)[data1$fm22[[k+1]]]
      expl.m<-cbind(expl.m,temp.cov)
      j<-(1:length(fm2[[1]]))[fm2[[1]]==c2[i]]+1
      if(length(j)==1) #binary
      {tmp.m<-rbinom(n2.new,1,prob=predict(fm2[[j]],newdata=data.frame(expl.m),type="response"))
      m.new[,c2[i]]<-one(ifelse(tmp.m==0,1,2),level.new)}
      else               #multi-categorical
      {tmp.logit<-rep(1,n2.new)
      for (k in 1:length(j))
        tmp.logit<-cbind(tmp.logit,predict(fm2[[j[k]]],newdata=data.frame(expl.m)))
      tmp.m<-apply(tmp.logit,1,rmultinom,n=1,size=1)
      m.new[,c2[i]]<-one((1:nrow(tmp.m))%*%tmp.m,level)   #the reference group is always 1
      }
    }           #glm
  m.new
}

biny<-FALSE
surv=is.Surv(y)

if(!surv){
temp.1<-find_level(y,level)
if(temp.1[[1]]<=2)
{levely=1
if(temp.1[[1]]==1)
{biny=TRUE
if(is.null(yref))
  y<-ifelse(y==temp.1[[2]],0,1)
else
  y<-ifelse(y==yref,0,1)}
}
else
{levely=2
if(temp.1[[1]]==3)
{biny=TRUE
if(is.null(yref))
  y<-ifelse(y==temp.1[[2]],0,1)
else
  y<-ifelse(y==yref,0,1)}
}
}

n2<-length(unique(level[!is.na(level.new)]))
#1a. an anlaysis on the whole data set
if(is.null(data1))
     data1<-data.org(x=x, m=m, levely=levely, y=y, levelx=levelx, xref=xref, l1=l1,l2=l2, c1=c1, 
                     c1r=c1r,c2=c2, c2r=c2r, f01y=f01y, f10y=f10y,f02ky=f02ky, f20ky=f20ky, 
                     f01km1=f01km1, f01km2=f01km2, f10km=f10km, level=level, weight=weight)
x=data1$parameter$x
m=data1$parameter$m
levelx=data1$parameter$levelx
levely=data1$parameter$levely
l1=data1$parameter$l1
l2=data1$parameter$l2
c1=data1$parameter$c1 
c1r=data1$parameter$c1r
c2=data1$parameter$c2 
c2r=data1$parameter$c2r
binx=data1$binx

mnames<-colnames(m)
x.names<-colnames(x)

#1b. prepare new data
if(is.null(data2) & is.null(x.new))
  {data2=data1
   x.new=x
   m.new=m
   weight.new=weight}
else if (is.null(data2))
  {if(is.null(weight.new))
    weight.new=rep(1,nrow(x.new))
   if(is.null(m.new))
     m.new1=m[sample(1:nrow(m),size=nrow(x.new),replace=T),]
   else
     m.new1=m.new
   data2<-data.org(x=x.new, m=m.new1, levely=levely, y=NULL, levelx=levelx, xref=xref, l1=l1,l2=l2, 
                  c1=c1, c1r=c1r, c2=c2, c2r=c2r, f01y=f01y, f10y=f10y,
                  f02ky=f02ky, f20ky=f20ky, f01km1=f01km1, f01km2=f01km2, f10km=f10km, 
                  level=level.new, weight=weight.new)
   x.new=data2$parameter$x
 #  m.new=data2$parameter$m
  }
else #when all .new comes from data2
  {x.new=data2$parameter$x
   m.new=data2$parameter$m
   weight.new=data2$parameter$weight
  }
  
full<-mlma(y=y, data1=data1, weight=weight, 
           random=random, random.m1=random.m1,intercept=intercept, family1=family1, familym=familym,
           covariates=covariates, cy1=cy1, cy2=cy2, cm=cm, joint=joint, data2=data2,x.new=x.new, 
           m.new=data2$parameter$m, level.new=level.new,weight.new=weight.new,covariates.new=covariates.new,
           cov.mat=cov.mat)

if(is.null(m.new))  #generate the new mediator from the full model if there is no m.new assigned.
 {m.new<-generate.m(x.new,data1,l1,c1,l2,c2,full,covariates.new,level.new,data2,cm)
  data2<-data.org(x=x.new, m=m.new, levely=levely, y=NULL, levelx=levelx, xref=xref, l1=l1,l2=l2, 
                c1=c1, c1r=c1r, c2=c2, c2r=c2r, f01y=f01y, f10y=f10y,
                f02ky=f02ky, f20ky=f20ky, f01km1=f01km1, f01km2=f01km2, f10km=f10km, 
                level=level.new, weight=weight.new)
  m.new=data2$parameter$m
  full<-mlma(y=y, data1=data1, weight=weight, 
             random=random, random.m1=random.m1,intercept=intercept, family1=family1, familym=familym,
             covariates=covariates, cy1=cy1, cy2=cy2, cm=cm, joint=joint, data2=data2,x.new=x.new, 
             m.new=m.new, level.new=level.new,weight.new=weight.new,covariates.new=covariates.new,cov.mat=cov.mat)}
  
#data2<-data.org(x.new, levelx, levely, m.new, l1, l2, c1, c1r, c2, c2r, f01y, f10y, f02ky, f20ky,
#                f01km1, f01km2, f10km, level.new, weight.new)
#2. bootstrap sample and prepare for data
ade1.boot<-NULL
ade2.boot<-NULL
aie1.boot<-NULL
aje1.boot<-NULL
aie12.boot<-NULL
aie2.boot<-NULL
aje2.boot<-NULL
aje12.boot<-NULL
ate1.boot<-NULL
ate2.boot<-NULL
coef.f1<-full$coef.f1
coef.fm1<-full$coef.fm1
coef.fm2<-full$coef.fm2
if(plot.it){
ie1.boot<-NULL
ie12.boot<-NULL
ie2.boot<-NULL
ie1_1.boot<-NULL
ie1_2.boot<-NULL
ie2_1.boot<-NULL
ie2_2.boot<-NULL
ie2_3.boot<-NULL
ie1_3.boot<-NULL
de1.boot<-NULL
de2.boot<-NULL
je1.boot<-NULL
je2.boot<-NULL
je12.boot<-NULL
te1.boot<-NULL
te2.boot<-NULL
}
m<-data.frame(m)
#level=droplevels(level) #remove all levels that are empty
n1=length(level)
temp1=1:n1
t3=aggregate(temp1~level, 
             FUN=function(x) c(x))
level.boot=sort(level)
for(l in 1:boot)
{#a. create bootsample
  func<-function(vec,weight)
  {if(length(vec)==1)
    temp=vec
  else
    temp=sample(vec,replace=TRUE,prob=weight[vec])
  temp
  }
  
  temp.2=unlist(sapply(t3[,2],func,weight))
  y.boot<-y[temp.2] 
  m.boot<-m[temp.2,]
  x.boot<-x[temp.2,]
  if(!is.null(covariates))
    covariates.boot<-covariates[temp.2,]
  else
    covariates.boot=NULL
  weight.boot<-weight[temp.2]

 if(!is.null(m.boot))
 {m.boot<-data.frame(m.boot)
 colnames(m.boot)<-mnames
 }
 if(!is.null(x.boot))  #remove dim()
 {x.boot<-data.frame(x.boot)
  colnames(x.boot)<-x.names
 }
 if(!is.null(covariates.boot)) #remove dim()
 {x.boot<-data.frame(x.boot)
 colnames(x.boot)<-x.names
 }

#b. create the boot data1
  data1.boot<-update.data.org(data1=data1,bootsample=temp.2,y.boot,x.boot,m.boot,
                              weight.boot,level.boot)  #data1 is the original data.org results, bootsample is booted sample
  

#c. analysis on the boot data
 full.boot<-mlma(y=y.boot, data1=data1.boot, weight=weight.boot, 
            random=random, random.m1=random.m1,intercept=intercept, 
            family1=family1, familym=familym,covariates=covariates.boot, 
            cy1=cy1, cy2=cy2, cm=cm, joint=joint, data2=data2,x.new=x.new, 
            m.new=m.new, level.new=level.new,weight.new=weight.new,
            covariates.new=covariates.new) 

#3. organize the results
 ade1.boot<-rbind(ade1.boot,full.boot$ade1)
 ade2.boot<-rbind(ade2.boot,full.boot$ade2)
 
 ate1.boot<-rbind(ate1.boot,full.boot$ate1)
 ate2.boot<-rbind(ate2.boot,full.boot$ate2)
 
 aie1.boot<-rbind(aie1.boot,full.boot$aie1)
 aie2.boot<-rbind(aie2.boot,full.boot$aie2)
 aie12.boot<-rbind(aie12.boot,full.boot$aie12)

 aje1.boot<-rbind(aje1.boot,full.boot$aje1)
 aje2.boot<-rbind(aje2.boot,full.boot$aje2)
 aje12.boot<-rbind(aje12.boot,full.boot$aje12)
 
 coef.f1=rbind(coef.f1, full.boot$coef.f1)
 if(length(coef.fm1)>1)
   for (i in 2:length(coef.fm1))
     coef.fm1[[i]]=rbind(coef.fm1[[i]], full.boot$coef.fm1[[i]])
 if(length(coef.fm2)>1)
   for (i in 2:length(coef.fm2))
     coef.fm2[[i]]=rbind(coef.fm2[[i]], full.boot$coef.fm2[[i]])

 if(plot.it){
 ie1.boot<-abind(ie1.boot,full.boot$ie1,along=1)
 ie2.boot<-abind(ie2.boot,full.boot$ie2,along=1)
 ie12.boot<-abind(ie12.boot,full.boot$ie12,along=1)
 ie1_2.boot<-abind(ie1_2.boot,full.boot$ie1_2,along=1)
 ie2_1.boot<-rbind(ie2_1.boot,full.boot$ie2_1)
 ie1_1.boot<-rbind(ie1_1.boot,full.boot$ie1_1)
 ie2_2.boot<-abind(ie2_2.boot,full.boot$ie2_2,along=1)
 ie1_3.boot<-rbind(ie1_3.boot,full.boot$ie1_3)
 ie2_3.boot<-rbind(ie2_3.boot,full.boot$ie2_3)
 de1.boot<-rbind(de1.boot,full.boot$de1)
 de2.boot<-rbind(de2.boot,full.boot$de2)
 te1.boot<-rbind(te1.boot,full.boot$te1)
 te2.boot<-rbind(te2.boot,full.boot$te2)
 je1.boot<-abind(je1.boot,full.boot$je1,along=1)
 je2.boot<-abind(je2.boot,full.boot$je2,along=1)
 je12.boot<-abind(je12.boot,full.boot$je12,along=1)
 }
 
if(echo) 
 print(l)
}

if(plot.it)
  a<-list(de1=de1.boot,de2=de2.boot,ie1=ie1.boot,ie2=ie2.boot,ie12=ie12.boot,   
          te1=te1.boot,te2=te2.boot,je1=je1.boot,je2=je2.boot,je12=je12.boot,
          ade1=ade1.boot,ade2=ade2.boot,aie1=aie1.boot,aie2=aie2.boot,
          aie12=aie12.boot,aje1=aje1.boot,aje2=aje2.boot,aje12=aje12.boot,
          ate1=ate1.boot,ate2=ate2.boot,ie1_1=ie1_1.boot,ie1_2=ie1_2.boot,
          ie1_3=ie1_3.boot,ie2_1=ie2_1.boot,ie2_2=ie2_2.boot,ie2_3=ie2_3.boot,
          full=full,levelx=levelx,level=level.new,x=x.new, weight=weight.new,  
          m=m.new, boot=boot, plot.it=plot.it,coef.f1=coef.f1,coef.fm1=coef.fm1,coef.fm2=coef.fm2) #, xboot=xboot, xjboot=xjboot
else
  a<-list(ade1=ade1.boot,ade2=ade2.boot,aie1=aie1.boot,aie2=aie2.boot,aie12=aie12.boot,
          aje1=aje1.boot,aje2=aje2.boot,aje12=aje12.boot,ate1=ate1.boot,ate2=ate2.boot,
          full=full,levelx=levelx, level=level.new,x=x.new, weight=weight.new, m=m.new,
          boot=boot, plot.it=plot.it,coef.f1=coef.f1,coef.fm1=coef.fm1,coef.fm2=coef.fm2) #, xboot=xboot, xjboot=xjboot
class(a)<-"mlma.boot"
return(a)
}

print.mlma<-function(x,...,w2=rep(1,nrow(as.matrix(object$de2))),digits=2)
{object<-x
if(!is.null(object$ate2))
 {cat("Level 2 Third Variable Effects: \n")
  temp=cbind(object$ate2,object$ade2,object$aie2,object$aie12,object$aje2)
  colnames(temp)=c("TE","DE",colnames(object$aie2),colnames(object$aie12),colnames(object$aje2))
  rownames(temp)=rownames(object$aie2)
  print(round(temp,digits))}
if(!is.null(object$ate1))
  {cat("Level 1 Third Variable Effects: \n")
   temp=cbind(object$ate1,object$ade1,object$aie1,object$aje1)
   colnames(temp)=c("TE","DE",colnames(object$aie1),colnames(object$aje1))
   rownames(temp)=rownames(object$aie1)
   print(round(temp,digits))}
}

summary.mlma<-function(object,...,type="III")
{if(!object$surv)
  f1<-Anova(object$f1,type=type)
 else
  f1<-object$f1 
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
plot.mlma<-function(x,..., var=NULL, cate=FALSE, w2=rep(1,nrow(as.matrix(object$de2))))
{object<-x
x<-object$x
if(!is.null(var) & is.character(var))   #change var to number of column in m
{temp.name<-colnames(object$m)
 var1<-grep(var,temp.name)}

two<-function(x, level, weight=rep(1,nrow(as.matrix(x))))
{x2<-NULL

if(is.null(dim(x)))
  x2=cbind(x2,(aggregate(as.numeric(x)~level, na.action=na.pass, 
                         FUN=function(x) c(mean=weighted.mean(x,weight=weight,na.rm=T))))[,2])
else
  for(i in 1:ncol(x))
    x2<-cbind(x2,(aggregate(as.numeric(x[,i])~level, na.action=na.pass,
                            FUN=function(x) c(mean=weighted.mean(x,weight=weight,na.rm=T))))[,2])
  colnames(x2)<-colnames(x)
  x2
}


marg.den<-function(x,y)
{y<-y[!is.na(x)]
x<-x[!is.na(x)]
z1<-sort(unique(x))
z2<-rep(0,length(z1))
for (i in 1:length(z1))
  z2[i]<-mean(y[x==z1[i]],na.rm=T)
cbind(z1,z2)
}

oldpar <- par(no.readonly = TRUE)  #line i
on.exit(par(oldpar)) #line i+1

if(is.null(var))
{par(mfrow=c(length(c(object$ate1,object$ate2)),1))
 if(!is.null(object$ate1))
   for(i in 1:length(object$ate1))
    {re<-c(object$ade1[i],object$aie1[i,],object$aje1[i,])/object$ate1[i]
     d<-order(re)
     name1<-c("DE",colnames(object$aie1),colnames(object$aje1))
     barplot(re[d],horiz=TRUE,xlab=paste("Relative Effects of Level 1 Mediation Effects 
             with Exposure Variable", rownames(object$aie1)[i],sep=" "),
             names=name1[d], cex.names=0.9,beside=FALSE,cex.axis=0.9,las=1,
             xlim=range(re), col = rainbow(length(d), start = 3/6, end = 4/6))}
  if(!is.null(object$ate2))
    for(i in 1:length(object$ate2))
    {re<-c(object$ade2[i],object$aie2[i,],object$aie12[i,],object$aje2[i,])/object$ate2[i]
     d<-order(re)
     name1<-c("DE",colnames(object$aie2),colnames(object$aie12),colnames(object$aje2))
     barplot(re[d],horiz=TRUE,xlab=paste("Relative Effects of Level 2 Mediation Effects 
             with Exposure Variable", rownames(object$aie2)[i],sep=" "),
             names=name1[d], cex.names=0.9,beside=FALSE,cex.axis=0.9,las=1,
             xlim=range(re), col = rainbow(length(d), start = 3/6, end = 4/6))}
}
else{
  m<-object$m[,var]
  y<-predict(object$f1)
  cate=F
  
  if(!is.factor(m) & !is.character(m) & nlevels(as.factor(m))>2)
    m.j<-two(m,object$data2$parameter$level)
  else
    {cate=T
     if(length(grep(var,colnames(object$data2$m1y)))!=0)
       m.j<-two(object$data2$m1y[,grep(var,colnames(object$data2$m1y))],object$data2$parameter$level)
     else
       m.j<-two(object$data2$m2y[,grep(var,colnames(object$data2$m2y))],object$data2$parameter$level)
    }
  if(is.null(names(object$f1)))
    y.j<-two(y,object$data2$parameter$level)
  else
    y.j<-two(y,object$data2$parameter$level[-object$f1$na.action])
  a1<-unique(c(grep(var,colnames(object$aie1)),grep(var,colnames(object$aie12))))
  a2<-grep(var,colnames(object$aie2))
  levelx=object$data1$parameter$levelx
  x1.names=colnames(object$te1)
  x2.names=colnames(object$te2)
  x.names=colnames(object$x)
  x.names.2=colnames(object$x.j)
  if(length(a1)>0)
  {if(!is.null(object$aie1))
    for (j in 1:nrow(object$aie1)){
      j1=match(rownames(object$aie1)[j],x.names)
    if(!object$data1$binx[j1])
    {par(mfcol=c(3,length(a1)))
     for(i in 1:length(a1))
     {a<-marg.den(object$x[,j1],object$ie1[,a1[i],j])
      scatter.smooth(a[,1],a[,2],xlab="x",ylab=paste("Level 1 IE of", colnames(object$aie1)[a1[i]],sep=" "),
                     main=paste("IE with exposure variable", rownames(object$aie1)[j],sep=" "))
      a<-marg.den(object$x[,j],object$ie1_2[,a1[i],j])
      suppressWarnings(plot(a,xlab="x",ylab=colnames(object$aie1)[a1[i]],type="l",
                     main=paste("Differential effect of", rownames(object$aie1)[j], "and",colnames(object$ie1)[a1[i]],sep=" ")))
      if(is.null(names(object$f1)))
        m1=m
      else
        m1=m[-object$f1$na.action]
      if(!cate)
      {a<-marg.den(m1,y)
       suppressWarnings(scatter.smooth(a[,1],a[,2],xlab=var,ylab="y",
                                      main=paste("Predicted relationship between y and",var,sep=" ")))
      }
      else
        plot(y~as.factor(m1), ylab="y",xlab=var,
             main=paste("Predicted relationship \n between y and",colnames(object$aie1)[a1[i]],sep=" ")) 
     }}
    else 
     {par(mfcol=c(2,1))
      if(!cate)
       plot(m~as.factor(object$x[,j1]),xlab=rownames(object$aie1)[j],ylab=var) 
      else
       plot(as.factor(m)~as.factor(object$x[,j1]),xlab=rownames(object$aie1)[j],ylab=var)
    if(is.null(names(object$f1)))
      m1=m
    else
      m1=m[-object$f1$na.action]
       if(!cate)
       {a<-marg.den(m1,y)
        suppressWarnings(scatter.smooth(a[,1],a[,2],xlab=var,ylab="y",
                                       main=paste("Predicted relationship between y and",var,sep=" ")))
       }
       else
        plot(y~as.factor(m1), ylab="y",xlab=var,
             main=paste("Predicted relationship \n between y and",colnames(object$aie1)[i],sep=" ")) 
     } 
   #readline("Press <return> to continue") 
   }
    if(!is.null(object$aie12))
      for (j in 1:nrow(object$aie12)){
        j1=match(rownames(object$aie12)[j],x.names)
        j2=match(rownames(object$aie12)[j],x.names.2)
        if(!object$data1$binx[j1])
        {par(mfcol=c(3,length(a1)))
             for(i in 1:length(a1))
             {a<-marg.den(object$x.j[,j2],object$ie12[,a1[i],j])
             scatter.smooth(a[,1],a[,2],xlab="x",ylab=paste("Level 2 IE of", colnames(object$aie12)[a1[i]],sep=" "),
                            main=paste("IE with exposure variable", rownames(object$aie12)[j],sep=" "))
             a<-marg.den(object$x.j[,j2],two(object$ie1_2[,a1[i],j],object$data2$parameter$level))
             suppressWarnings(plot(a,xlab="x",ylab=colnames(object$aie12)[a1[i]],type="l",
                                   main=paste("Differential effect of", rownames(object$aie12)[j], "and",colnames(object$ie12)[a1[i]],sep=" ")))
             if(!cate)
             {a<-marg.den(m.j,y.j)
             suppressWarnings(scatter.smooth(a[,1],a[,2],xlab=var,ylab="y",
                                             main=paste("Predicted relationship between y and",var,sep=" ")))
             }
             else
             {a<-marg.den(m.j[,i],y.j)
             suppressWarnings(scatter.smooth(a[,1],a[,2],xlab=var,ylab="y",
                                             main=paste("Predicted relationship between y and",colnames(m.j)[i],sep=" ")))
             }
             #plot(y.j~as.factor(m.j[,i]), ylab="y",xlab=var,
            #        main=paste("Predicted relationship \n between y and",colnames(object$aie12)[a1[i]],sep=" ")) 
             }}
        else 
        {par(mfcol=c(2,1))
          if(!cate)
            plot(m.j~as.factor(object$x.j[,j2]),xlab=rownames(object$aie12)[j],ylab=var) 
          else
            plot(as.factor(m.j)~as.factor(object$x.j[,j2]),xlab=rownames(object$aie12)[j],ylab=var)
          if(!cate)
          {a<-marg.den(m.j,y.j)
          suppressWarnings(scatter.smooth(a[,1],a[,2],xlab=var,ylab="y",
                                          main=paste("Predicted relationship between y and",var,sep=" ")))
          }
          else
            plot(y.j~as.factor(m.j), ylab="y",xlab=var,
                 main=paste("Predicted relationship \n between y and",colnames(object$aie12)[i],sep=" ")) 
        } 
        #readline("Press <return> to continue") 
        }
  }
  
  if(length(a2)>0){   
    for (j in 1:nrow(object$aie2)){
      j1=grep(rownames(object$aie2)[j],x.names)
      j2=grep(rownames(object$aie2)[j],x.names.2)
      if(!object$data1$binx[j1])
      {par(mfcol=c(3,length(a2)))
           for(i in a2)
           {a<-marg.den(object$x.j[,j2],object$ie2[,i,j])
            scatter.smooth(a[,1],a[,2],xlab="x",ylab=paste("Level 2 IE of", colnames(object$aie2)[i],sep=" "),
                          main=paste("IE with exposure variable", rownames(object$aie2)[j],sep=" "))
            a<-marg.den(object$x.j[,j2],object$ie2_2[,i,j])
            suppressWarnings(plot(a,xlab="x",ylab=colnames(object$aie2)[i],type="l",
                          main=paste("Differential effect of", rownames(object$aie2)[j], "and",colnames(object$ie2)[i],sep=" ")))
           if(!cate)
           {a<-marg.den(m.j,y.j)
            suppressWarnings(scatter.smooth(a[,1],a[,2],xlab=var,ylab="y",
                             main=paste("Predicted relationship between y and",var,sep=" ")))
           }
           else
             plot(y.j~as.factor(m.j), ylab="y",xlab=var,
                  main=paste("Predicted relationship \n between y and",colnames(object$aie2)[i],sep=" ")) 
           }}
      else 
      {par(mfcol=c(2,1))
        if(!cate)
          plot(m.j~as.factor(object$x.j[,j2]),xlab=rownames(object$aie2)[j],ylab=var) 
        else
          plot(as.factor(m.j)~as.factor(object$x.j[,j2]),xlab=rownames(object$aie2)[j],ylab=var)
        if(!cate)
        {a<-marg.den(m.j,y.j)
         suppressWarnings(scatter.smooth(a[,1],a[,2],xlab=var,ylab="y",
                          main=paste("Predicted relationship between y and",var,sep=" ")))
        }
        else
          plot(y.j~as.factor(m.j), ylab="y",xlab=var,
               main=paste("Predicted relationship \n between y and",colnames(object$aie2)[i],sep=" ")) 
      } 
      #readline("Press <return> to continue") 
      }}
}

#par(op)
}

print.mlma.boot<-function(x, ...)
{print(x$full)}

summary.mlma.boot<-function(object,...,alpha=0.05, RE=FALSE,digits=4)
{summarize.boot<-function(vector,n,a1,a2,b1,b2)  #vector is stacked results from bootstrap samples
  #n is the number of elements from each bootstrap sample
{mat<-matrix(vector,length(vector)/n,n)
all<-summarize.boot1(t(mat),a1,a2,b1,b2)
data.frame(all)
}
summarize.boot.re<-function(vector,n,te,a1,a2,b1,b2)  #vector is stacked results from bootstrap samples
  #n is the number of elements from each bootstrap sample
{mat<-t(matrix(vector,length(vector)/n,n))*(1/te)
all<-summarize.boot1(mat,a1,a2,b1,b2)
data.frame(all)
}
summarize.boot1<-function(mat,a1,a2,b1,b2){
mean.temp=apply(mat,2,mean,na.rm=T)
sd.temp=apply(mat,2,sd,na.rm=T)
cbind(mean=mean.temp,sd=sd.temp,upbd=mean.temp+b2*sd.temp,
      lwbd=mean.temp+b1*sd.temp,upbd.quan=apply(mat,2,quantile,a2,na.rm=T),
      lwbd.quan=apply(mat,2,quantile,a1,na.rm=T))
}

a1<-alpha/2
a2<-1-a1
b1<-qnorm(a1)
b2<-qnorm(a2)
boot=object$boot
direct.effect1<-NULL
direct.effect2<-NULL
indirect.effect1<-NULL
indirect.effect2<-NULL
indirect.effect12<-NULL
jdirect.effect1<-NULL
jdirect.effect2<-NULL
total.effect1<-NULL
total.effect2<-NULL
effect1=NULL
effect2=NULL
re1=NULL
re2=NULL

if(!is.null(object$ate2))
{direct.effect2<-cbind(est=object$full$ade2,summarize.boot1(object$ade2,a1,a2,b1,b2))
total.effect2<-cbind(est=object$full$ate2,summarize.boot1(object$ate2,a1,a2,b1,b2))
if(!is.null(object$aie2))
  {indirect.effect2.1<-apply(object$aie2,2,summarize.boot,boot,a1,a2,b1,b2)
  if(is.list(indirect.effect2.1))
    {for (i in 1:nrow(object$full$aie2))
       {temp=NULL
        for (j in 1:ncol(object$full$aie2))
           temp=cbind(temp,as.numeric(indirect.effect2.1[[j]][i,]))
        indirect.effect2[[i]]=temp
        indirect.effect2[[i]]=rbind(object$full$aie2[i,],indirect.effect2[[i]])

rownames(indirect.effect2[[i]])=c("est","mean","sd","upbd","lwbd","upbd.quan","lwbd.quan")
        colnames(indirect.effect2[[i]])=colnames(object$aie2)}
     names(indirect.effect2)=rownames(object$full$aie2)}
   else
   {indirect.effect2=rbind(object$full$aie2,indirect.effect2.1)
    rownames(indirect.effect2)=c("est","mean","sd","upbd","lwbd","upbd.quan","lwbd.quan")}}
if(!is.null(object$aie12))
{indirect.effect12.1<-apply(object$aie12,2,summarize.boot,boot,a1,a2,b1,b2)
if(is.list(indirect.effect12.1))
{for (i in 1:nrow(object$full$aie12))
{temp=NULL
 for (j in 1:ncol(object$full$aie12))
  temp=cbind(temp,as.numeric(indirect.effect12.1[[j]][i,]))
 indirect.effect12[[i]]=temp
indirect.effect12[[i]]=rbind(object$full$aie12[i,],indirect.effect12[[i]])
rownames(indirect.effect12[[i]])=c("est","mean","sd","upbd","lwbd","upbd.quan","lwbd.quan")
colnames(indirect.effect12[[i]])=colnames(object$aie12)}
  names(indirect.effect12)=rownames(object$full$aie12)}
else
{indirect.effect12=rbind(object$full$aie12,indirect.effect12.1)
rownames(indirect.effect12)=c("est","mean","sd","upbd","lwbd","upbd.quan","lwbd.quan")}}

if(!is.null(object$aje2))
{jdirect.effect2.1<-apply(object$aje2,2,summarize.boot,boot,a1,a2,b1,b2)
if(is.list(jdirect.effect2.1))
{for (i in 1:nrow(object$full$aje2))
{temp=NULL
 for (j in 1:ncol(object$full$aje2))
  temp=cbind(temp,as.numeric(jdirect.effect2.1[[j]][i,]))
 jdirect.effect2[[i]]=temp
jdirect.effect2[[i]]=rbind(object$full$aje2[i,],jdirect.effect2[[i]])
rownames(jdirect.effect2[[i]])=c("est","mean","sd","upbd","lwbd","upbd.quan","lwbd.quan")
colnames(jdirect.effect2[[i]])=colnames(object$aje2)}
  names(jdirect.effect2)=rownames(object$full$aje2)}
else
{jdirect.effect2=rbind(object$full$aje2,jdirect.effect2.1)
rownames(jdirect.effect2)=c("est","mean","sd","upbd","lwbd","upbd.quan","lwbd.quan")}}

effect2=NULL
if(is.list(indirect.effect2))
 {for (i in 1:ncol(object$ate2))
   effect2[[i]]=cbind(te=total.effect2[i,],de=direct.effect2[i,],indirect.effect2[[i]],
                     indirect.effect12[[i]],jdirect.effect2[[i]])
  names(effect2)=colnames(object$ate2)}
else
  effect2=cbind(te=total.effect2[1,],de=direct.effect2[1,],indirect.effect2,
                indirect.effect12,jdirect.effect2)
}

if(!is.null(object$ate1))
{direct.effect1<-cbind(est=object$full$ade1,summarize.boot1(object$ade1,a1,a2,b1,b2))
 total.effect1<-cbind(est=object$full$ate1,summarize.boot1(object$ate1,a1,a2,b1,b2))
 if(!is.null(object$aie1))
 {indirect.effect1.1<-apply(object$aie1,2,summarize.boot,boot,a1,a2,b1,b2)
 if(is.list(indirect.effect1.1))
 {for (i in 1:nrow(object$full$aie1))
 {indirect.effect1[[i]]<-as.numeric(indirect.effect1.1[[1]][i,])
  if(ncol(object$full$aie1)>1)
  for (j in 2:ncol(object$full$aie1))
   indirect.effect1[[i]]=cbind(indirect.effect1[[i]],as.numeric(indirect.effect1.1[[j]][i,]))
 if(is.null(dim(indirect.effect1[[i]])))  #if there is only one mediator
   indirect.effect1[[i]]=matrix(c(object$full$aie1[i,],indirect.effect1[[i]]),ncol=1)
 else
   indirect.effect1[[i]]=rbind(object$full$aie1[i,],indirect.effect1[[i]])
 rownames(indirect.effect1[[i]])=c("est","mean","sd","upbd","lwbd","upbd.quan","lwbd.quan")
 colnames(indirect.effect1[[i]])=colnames(object$aie1)}
   names(indirect.effect1)=rownames(object$full$aie1)}
 else
 {indirect.effect1=rbind(object$full$aie1,indirect.effect1.1)
 rownames(indirect.effect1)=c("est","mean","sd","upbd","lwbd","upbd.quan","lwbd.quan")}}

 if(!is.null(object$aje1))
 {jdirect.effect1.1<-apply(object$aje1,2,summarize.boot,boot,a1,a2,b1,b2)
 if(is.list(jdirect.effect1.1))
 {for (i in 1:nrow(object$full$aje1))
 {for (j in 1:ncol(object$full$aje1))
   jdirect.effect1[[i]]=cbind(jdirect.effect1[[i]],as.numeric(jdirect.effect1.1[[j]][i,]))
 jdirect.effect1[[i]]=rbind(object$full$aje1[i,],jdirect.effect1[[i]])
 rownames(jdirect.effect1[[i]])=c("est","mean","sd","upbd","lwbd","upbd.quan","lwbd.quan")
 colnames(jdirect.effect1[[i]])=colnames(object$aje1)}
   names(jdirect.effect1)=rownames(object$full$aje1)}
 else
 {jdirect.effect1=rbind(object$full$aje1,jdirect.effect1.1)
 rownames(jdirect.effect1)=c("est","mean","sd","upbd","lwbd","upbd.quan","lwbd.quan")}}
 
 effect1=NULL
 if(is.list(indirect.effect1))
   {for (i in 1:ncol(object$ate1))
     effect1[[i]]=cbind(te=total.effect1[i,],de=direct.effect1[i,],indirect.effect1[[i]],
                       jdirect.effect1[[i]])
    names(effect1)=colnames(object$ate1)}
 else
   effect1=cbind(te=total.effect1[1,],de=direct.effect1[1,],indirect.effect1,
                 jdirect.effect1)
}

direct.re1<-NULL
direct.re2<-NULL
indirect.re1<-NULL
indirect.re2<-NULL
indirect.re12<-NULL
jdirect.re1<-NULL
jdirect.re2<-NULL
total.re1<-NULL
total.re2<-NULL

if(!is.null(object$ate2))
{direct.re2<-cbind(est=object$full$ade2/object$full$ate2,summarize.boot1(object$ade2/object$ate2,a1,a2,b1,b2))
 if(!is.null(object$aie2))
 {indirect.re2.1<-apply(object$aie2,2,summarize.boot.re,boot,object$ate2,a1,a2,b1,b2)
 if(is.list(indirect.re2.1))
 {for (i in 1:nrow(object$full$aie2))
 {temp=NULL
   for (j in 1:ncol(object$full$aie2))
    temp=cbind(temp,as.numeric(indirect.re2.1[[j]][i,]))
   indirect.re2[[i]]=temp
  
 indirect.re2[[i]]=rbind(object$full$aie2[i,]/object$full$ate2[i],indirect.re2[[i]])
 rownames(indirect.re2[[i]])=c("est","mean","sd","upbd","lwbd","upbd.quan","lwbd.quan")
 colnames(indirect.re2[[i]])=colnames(object$aie2)}
   names(indirect.re2)=rownames(object$full$aie2)}
 else
 {indirect.re2=rbind(object$full$aie2/object$full$ate2,indirect.re2.1)
 rownames(indirect.re2)=c("est","mean","sd","upbd","lwbd","upbd.quan","lwbd.quan")}}

if(!is.null(object$aie12))
{indirect.re12.1<-apply(object$aie12,2,summarize.boot.re,boot,object$ate2,a1,a2,b1,b2)
if(is.list(indirect.re12.1))
{for (i in 1:nrow(object$full$aie12))
{temp=NULL
  for (j in 1:ncol(object$full$aie12))
   temp=cbind(temp,as.numeric(indirect.re12.1[[j]][i,]))
 indirect.re12[[i]]=temp
indirect.re12[[i]]=rbind(object$full$aie12[i,]/object$full$ate2[i],indirect.re12[[i]])
rownames(indirect.re12[[i]])=c("est","mean","sd","upbd","lwbd","upbd.quan","lwbd.quan")
colnames(indirect.re12[[i]])=colnames(object$aie12)}
  names(indirect.re12)=rownames(object$full$aie12)}
else
{indirect.re12=rbind(object$full$aie12/object$full$ate2,indirect.re12.1)
rownames(indirect.re12)=c("est","mean","sd","upbd","lwbd","upbd.quan","lwbd.quan")}}

if(!is.null(object$aje2))
{jdirect.re2.1<-apply(object$aje2,2,summarize.boot.re,boot,object$ate2,a1,a2,b1,b2)
if(is.list(jdirect.re2.1))
{for (i in 1:nrow(object$full$aje2))
{temp=NULL
 for (j in 1:ncol(object$full$aje2))
  temp=cbind(temp,as.numeric(jdirect.re2.1[[j]][i,]))
 jdirect.re2[[i]]=temp
jdirect.re2[[i]]=rbind(object$full$aje2[i,]/object$full$ate2[i],jdirect.re2[[i]])
rownames(jdirect.re2[[i]])=c("est","mean","sd","upbd","lwbd","upbd.quan","lwbd.quan")
colnames(jdirect.re2[[i]])=colnames(object$aje2)}
  names(jdirect.re2)=rownames(object$full$aje2)}
else
{jdirect.re2=rbind(object$full$aje2/object$full$ate2,jdirect.re2.1)
rownames(jdirect.re2)=c("est","mean","sd","upbd","lwbd","upbd.quan","lwbd.quan")}}

re2=NULL
if(is.list(indirect.re2))
  {for (i in 1:ncol(object$ate2))
   re2[[i]]=cbind(de=direct.re2[i,],indirect.re2[[i]],
                      indirect.re12[[i]],jdirect.re2[[i]])
   names(re2)=colnames(object$ate2)}
else
  re2=cbind(de=direct.re2[1,],indirect.re2,
                indirect.re12,jdirect.re2)
}

if(!is.null(object$ate1))
{direct.re1<-cbind(est=object$full$ade1/object$full$ate1,summarize.boot1(object$ade1/object$ate1,a1,a2,b1,b2))
if(!is.null(object$aie1))
{indirect.re1.1<-apply(object$aie1,2,summarize.boot.re,boot,object$ate1,a1,a2,b1,b2)
if(is.list(indirect.re1.1))
{for (i in 1:nrow(object$full$aie1))
{indirect.re1[[i]]=as.numeric(indirect.re1.1[[1]][i,])
 if(ncol(object$full$aie1)>1)
  for (j in 2:ncol(object$full$aie1))
   indirect.re1[[i]]=cbind(indirect.re1[[i]],as.numeric(indirect.re1.1[[j]][i,]))
  if(is.null(dim(indirect.re1[[i]])))  #if there is only one mediator
    indirect.re1[[i]]=matrix(c(object$full$aie1[i,]/object$full$ate1[i],indirect.re1[[i]]),ncol=1)
  else
    indirect.re1[[i]]=rbind(object$full$aie1[i,]/object$full$ate1[i],indirect.re1[[i]])
rownames(indirect.re1[[i]])=c("est","mean","sd","upbd","lwbd","upbd.quan","lwbd.quan")
colnames(indirect.re1[[i]])=colnames(object$aie1)}
  names(indirect.re1)=rownames(object$full$aie1)}
else
{indirect.re1=rbind(object$full$aie1/object$full$ate1,indirect.re1.1)
rownames(indirect.re1)=c("est","mean","sd","upbd","lwbd","upbd.quan","lwbd.quan")}}

if(!is.null(object$aje1))
{jdirect.re1.1<-apply(object$aje1,2,summarize.boot.re,boot,object$ate1,a1,a2,b1,b2)
if(is.list(jdirect.re1.1))
{for (i in 1:nrow(object$full$aje1))
{for (j in 1:ncol(object$full$aje1))
  jdirect.re1[[i]]=cbind(jdirect.re1[[i]],as.numeric(jdirect.re1.1[[j]][i,]))
jdirect.re1[[i]]=rbind(object$full$aje1[i,]/object$full$ate1[i],jdirect.re1[[i]])
rownames(jdirect.re1[[i]])=c("est","mean","sd","upbd","lwbd","upbd.quan","lwbd.quan")
colnames(jdirect.re1[[i]])=colnames(object$aje1)}
  names(jdirect.re1)=rownames(object$full$aje1)}
else
{jdirect.re1=rbind(object$full$aje1/object$full$ate1,jdirect.re1.1)
rownames(jdirect.re1)=c("est","mean","sd","upbd","lwbd","upbd.quan","lwbd.quan")}}

re1=NULL
if(is.list(indirect.re1))
  {for (i in 1:ncol(object$ate1))
    re1[[i]]=cbind(te=total.re1[i,],de=direct.re1[i,],indirect.re1[[i]],
                      jdirect.re1[[i]])
   names(re1)=colnames(object$ate1)}
else
  re1=cbind(te=total.re1[1,],de=direct.re1[1,],indirect.re1,
                jdirect.re1)
}
a<-list(effect1=effect1,effect2=effect2,re1=re1,re2=re2,RE=RE,digits=digits)
class(a)<-"summary.mlma.boot"
return(a)
}

print.summary.mlma.boot<-function(x,...,digits=x$digits)
{object=x
  if(object$RE)
{if(!is.null(object$re1)) 
  {cat("MLMA Analysis: Estimated Relative Effects at level 1 (%):\n")
   if(is.list(object$re1))
    print(lapply(object$re1,round,digits))
   else
     print(round(object$re1,digits))}
  if(!is.null(object$re2))
   {cat("MLMA Analysis: Estimated Relative Effects at level 2 (%):\n")
    if(is.list(object$re2))
      print(lapply(object$re2,round,digits))
    else
      print(round(object$re2,digits))}
}
  else
  {if(!is.null(object$effect1)) 
    {cat("MLMA Analysis: Estimated Effects at level 1:\n")
     if(is.list(object$effect1))
       print(lapply(object$effect1,round,digits))
     else
       print(round(object$effect1,digits))}
   if(!is.null(object$effect2)) {
     cat("MLMA Analysis: Estimated Effects at level 2:\n")
     if(is.list(object$effect2))
       print(lapply(object$effect2,round,digits))
     else
       print(round(object$effect2,digits))}
  }
}

plot.mlma.boot<-function(x,..., var=NULL, alpha=0.05,quant=FALSE, plot.it=x$plot.it) #quantile=True to use quantile CIS
{object<-x

two<-function(x, level, weight=rep(1,nrow(as.matrix(x))))
{x2<-NULL
if(is.null(dim(x)))
  x2=cbind(x2,(aggregate(as.numeric(x)~level, na.action=na.pass, 
                         FUN=function(x) c(mean=weighted.mean(x,weight=weight,na.rm=T))))[,2])
else
  for(i in 1:ncol(x))
    x2<-cbind(x2,(aggregate(as.numeric(x[,i])~level, na.action=na.pass,
                            FUN=function(x) c(mean=weighted.mean(x,weight=weight,na.rm=T))))[,2])
  colnames(x2)<-colnames(x)
  x2
}

boot.ci<-function(x,mat,cri_val)
{
mat1=as.vector(t(mat))
x1=rep(x,each=ncol(mat))
temp1=aggregate(mat1~x1, 
                FUN=function(x) c(mean=mean(x,na.rm=T),sd_dev=sd(x,na.rm=T)), simplify = TRUE, drop = TRUE)
upbd<-temp1[,2][,1]+cri_val*temp1[,2][,2]
lwbd<-temp1[,2][,1]-cri_val*temp1[,2][,2]
return(data.frame(x=temp1[,1],F=temp1[,2][,1],L=lwbd,U=upbd))
}




plot_ci<-function(df,main="IE",xlab="x",ylab="IE",cate=FALSE)
  {if(!cate){
    plot(df$x, df$F, ylim = range(c(df$L,df$U),na.rm=T), type = "l",main=main,xlab=xlab,ylab=ylab)
    polygon(c(df$x,rev(df$x)),c(df$L,rev(df$U)),col = "grey75", border = FALSE)
    lines(df$x, df$F, lwd = 2)
    lines(df$x, df$U, col="red",lty=2)
    lines(df$x, df$L, col="red",lty=2)}
   else
     barplot2(df$F,ylab=ylab, names=df$x,  plot.ci = TRUE, ci.u = df$U, 
              ci.l = df$L, cex.names=0.9,beside=FALSE,cex.axis=0.9,las=1,
              ylim=range(c(df$U,df$L)))
  }
  
oldpar <- par(no.readonly = TRUE)  #line i
on.exit(par(oldpar)) #line i+1

 if(is.null(var))
 {summary.obj<-summary(object,alpha=alpha)
 par(mfrow=c(length(c(object$full$ate1,object$full$ate2)),1))
 if(!is.null(object$full$ate1))
   {if (length(object$full$ate1)==1)
   {re<-summary.obj$re1[[1]][1,]
   if(quant)
   {upper<-summary.obj$re1[[1]][6,]
    lower<-summary.obj$re1[[1]][7,]}
   else
   {lower<-summary.obj$re1[[1]][5,]
    upper<-summary.obj$re1[[1]][4,]}
   d<-order(re)
   name1<-c("DE",colnames(object$full$aie1),colnames(object$full$aje1))
   bp<-barplot2(re[d],horiz=TRUE,xlab=paste("Relative Effects of Level 1 Mediation Effects with Exposure Variable", rownames(object$full$aie1),sep=" "),
           names=name1[d],  plot.ci = TRUE, ci.u = upper[d], 
           ci.l = lower[d], cex.names=0.9,beside=FALSE,cex.axis=0.9,las=1,
           xlim=range(c(upper,lower)), col = rainbow(length(d), start = 3/6, end = 4/6))}
   else
   for(i in 1:length(object$full$ate1))
   {re<-summary.obj$re1[[i]][1,]
   if(quant)
   {upper<-summary.obj$re1[[i]][6,]
    lower<-summary.obj$re1[[i]][7,]}
   else
   {lower<-summary.obj$re1[[i]][5,]
    upper<-summary.obj$re1[[i]][4,]}
   d<-order(re)
   name1<-c("DE",colnames(object$aie1),colnames(object$aje1))
   bp<-barplot2(re[d],horiz=TRUE,xlab=paste("Relative Effects of Level 1 Mediation Effects with Exposure Variable", rownames(object$full$aie1)[i],sep=" "),
           names=name1[d],  plot.ci = TRUE, ci.u = upper[d], 
           ci.l = lower[d], cex.names=0.9,beside=FALSE,cex.axis=0.9,las=1,
           xlim=range(c(upper,lower)), col = rainbow(length(d), start = 3/6, end = 4/6))}}
 
   if(!is.null(object$ate2))
   {if (length(object$full$ate2)==1)
   {re<-summary.obj$re2[[1]][1,]
   if(quant)
   {upper<-summary.obj$re2[[1]][6,]
   lower<-summary.obj$re2[[1]][7,]}
   else
   {lower<-summary.obj$re2[[1]][5,]
   upper<-summary.obj$re2[[1]][4,]}
   d<-order(re)
   name1<-c("DE",colnames(object$full$aie2),colnames(object$full$aie12),colnames(object$full$aje2))
   bp<-barplot2(re[d],horiz=TRUE,xlab=paste("Relative Effects of Level 2 Mediation Effects with Exposure Variable", colnames(object$ate2),sep=" "),
           names=name1[d],  plot.ci = TRUE, ci.u = upper[d], 
           ci.l = lower[d], cex.names=0.9,beside=FALSE,cex.axis=0.9,las=1,
           xlim=range(c(upper,lower)), col = rainbow(length(d), start = 3/6, end = 4/6))}
     else
       for(i in 1:length(object$full$ate2))
       {re<-summary.obj$re2[[i]][1,]
       if(quant)
       {upper<-summary.obj$re2[[i]][6,]
       lower<-summary.obj$re2[[i]][7,]}
       else
       {lower<-summary.obj$re2[[i]][5,]
       upper<-summary.obj$re2[[i]][4,]}
       d<-order(re)
       name1<-c("DE",colnames(object$full$aie2),colnames(object$full$aie12),colnames(object$full$aje2))
       bp<-barplot2(re[d],horiz=TRUE,main=paste("Relative Effects of Level 2 Mediation Effects with Exposure Variable", colnames(object$ate2)[i],sep=" "),
               names=name1[d],  plot.ci = TRUE, ci.u = upper[d], 
               ci.l = lower[d], cex.names=0.9,beside=FALSE,cex.axis=0.9,las=1,
               xlim=range(c(upper,lower)), col = rainbow(length(d), start = 3/6, end = 4/6))}}
  }
 else if (plot.it) {
   object1=object$full
   m<-object1$m[,var]
   y<-predict(object1$f1)
   if(is.null(names(object1$f1)))
     m1=m
   else
     m1=m[-object1$f1$na.action]
   cate=F
   if(!is.factor(m) & !is.character(m) & nlevels(as.factor(m))>2)
     m.j<-two(m,object1$data2$parameter$level)
   else
   {cate=T
    if(length(grep(var,colnames(object1$data2$m1y)))!=0)
      m.j<-two(object1$data2$m1y[,grep(var,colnames(object1$data2$m1y))],object1$data2$parameter$level)
    else
      m.j<-two(object1$data2$m2y[,grep(var,colnames(object1$data2$m2y))],object1$data2$parameter$level)
   }
    if(is.null(names(object1$f1)))
      y.j<-two(y,object1$data2$parameter$level)   #what if y is binary?---get the linear part
    else
      y.j<-two(y,object1$data2$parameter$level[-object1$f1$na.action])
    
   a1<-unique(c(grep(var,colnames(object1$aie1)),grep(var,colnames(object1$aie12))))
   a2<-grep(var,colnames(object1$aie2))
   levelx=object1$data1$parameter$levelx
   x1.names=colnames(object1$te1)
   x2.names=colnames(object1$te2)
   x.names=colnames(object1$x)
   x.names.2=colnames(object1$x.j)
   cri_val<-qnorm((1+alpha)/2)
   
   if(length(a1)>0)
   {if(!is.null(object1$aie1))
     for (j in 1:nrow(object1$aie1)){
       j1=match(rownames(object1$aie1)[j],x.names)
       if(!object1$data1$binx[j1])
       {par(mfcol=c(3,length(a1)))
         for(i in 1:length(a1))
         {ie1<-boot.ci(object1$x[,j1],matrix(object$ie1[,a1[i],j],ncol=object$boot),cri_val)
         plot_ci(ie1,xlab=x.names[j1],ylab=paste("Level 1 IE of", colnames(object1$aie1)[a1[i]],sep=" "),
                 main=paste("IE with exposure variable", rownames(object1$aie1)[j],sep=" "))
         ie1_2<-boot.ci(object1$x[,j1],matrix(object$ie1_2[,a1[i],j1],,ncol=object$boot),cri_val)
         suppressWarnings(plot_ci(ie1_2,xlab=x.names[j1],ylab=colnames(object$aie1)[a1[i]],
                                  main=paste("Differential effect of", rownames(object$aie1)[j], "and",colnames(object$ie1)[a1[i]],sep=" ")))
         if(!cate)
         {ie1_1<-boot.ci(m,matrix(object$ie1_1[,a1[i]],ncol=object$boot),cri_val)
         plot_ci(ie1_1,xlab=var,ylab="y",main=paste("Predicted relationship between y and",var,sep=" "))
         }
         else
           boxplot(matrix(object$ie1_1[,a1[i]],ncol=object$boot)[1,], ylab="y",
                   xlab=dimnames(object$ie1)[[2]][a1[i]], 
                   main=paste("Predicted relationship \n between y and",colnames(object$aie1)[a1[i]],sep=" ")) 
         }}
       else
         {par(mfcol=c(2,1))
           if(!cate)
             plot(m~as.factor(object$x[,j1]),xlab=rownames(object1$aie1)[j],ylab=var,
                  main=paste("Relationship between",x.names[j1],"and",var,sep=" ")) 
           else
             plot(as.factor(m)~as.factor(object$x[,j1]),xlab=rownames(object1$aie1)[j],ylab=var,
                  main=paste("Relationship between",x.names[j1],"and",var,sep=" "))
           if(!cate)
           {ie1_1<-boot.ci(m,matrix(object$ie1_1[,a1],ncol=object$boot),cri_val)
           #browser()
           plot_ci(ie1_1,xlab=var,ylab="y",main=paste("Predicted relationship between y and",var,sep=" "))
           }
           else
             plot(y~as.factor(m1), ylab="y",xlab=var,
                  main=paste("Predicted relationship \n between y and",colnames(object$aie1)[i],sep=" ")) 
         }}   
  
     if(!is.null(object1$aie12))
       for (j in 1:nrow(object1$aie12)){
         j1=match(rownames(object1$aie12)[j],x.names)
         j2=match(rownames(object1$aie12)[j],x.names.2)
         if(!object1$data1$binx[j1])
         {par(mfcol=c(3,length(a1)))
           for(i in 1:length(a1))
           {        
           ie12<-boot.ci(object1$x.j[,j2],matrix(object$ie12[,a1[i],j],ncol=object$boot),cri_val)
           plot_ci(ie12,xlab=x.names.2[j2],ylab=paste("Level 2 IE of", colnames(object1$aie12)[a1[i]],sep=" "),
                   main=paste("IE with exposure variable", rownames(object1$aie12)[j],sep=" "))
           ie12_2<-boot.ci(object1$x[,j1],matrix(object$ie1_2[,a1[i],j1],ncol=object$boot),cri_val)
           suppressWarnings(plot_ci(ie12_2,xlab=x.names.2[j2],ylab=colnames(object$aie12)[a1[i]],
                                    main=paste("Differential effect of", rownames(object$aie12)[j], "and",
                                               colnames(object$ie12)[a1[i]],sep=" ")))
           if(!cate)
           {temp.matrix<-apply(matrix(object$ie1_1[,a1[i]],ncol=object$boot),2,two,object$level)
            ie12_1<-boot.ci(m.j,temp.matrix,cri_val)
           plot_ci(ie12_1,xlab=var,ylab="y",main=paste("Predicted relationship between y and",var,sep=" "))
           }
           else
               boxplot(matrix(object$ie1_1[,a1[i]],ncol=object$boot)[1,], ylab="y",
                     xlab=dimnames(object$ie12)[[2]][a1[i]], 
                     main=paste("Predicted relationship \n between y and",colnames(object$aie12)[a1[i]],sep=" ")) 
               }}
         else
         {par(mfcol=c(2,length(a1)))
     
           if(length(a1)>1)
             {plot(m.j[,1]~as.factor(object1$x.j[,j2]),xlab=rownames(object1$aie12)[j],ylab=var,
                  main=paste("Relationship between",x.names[j1],"and",colnames(m.j)[1],sep=" ")) 
           for (jj in 2:length(a1))
             plot(m.j[,j]~as.factor(object1$x.j[,j2]),xlab=rownames(object1$aie12)[j],ylab=var,
                main=paste("Relationship between",x.names[j1],"and",colnames(m.j)[jj],sep=" ")) }
           else
             plot(m.j~as.factor(object1$x.j[,j2]),xlab=rownames(object1$aie12)[j],ylab=var,
                  main=paste("Relationship between",x.names[j1],"and",var,sep=" ")) 
           #if(!cate){
#           browser()
           for (i in 1:length(a1))
            {temp.matrix<-apply(matrix(object$ie1_1[,a1[i]],ncol=object$boot),2,two,object$level)
            ie12_1<-boot.ci(m.j[,i],temp.matrix,cri_val)
            plot_ci(ie12_1,xlab=var,ylab="y",main=paste("Predicted relationship between y and",colnames(m.j)[i],sep=" "))
           }
 #          else
#             plot(y.j~as.factor(m.j), ylab="y",xlab=var,
#                  main=paste("Predicted relationship \n between y and",colnames(object$aie12)[i],sep=" ")) 
         }}
   }
   
   if(length(a2)>0)
     for (j in 1:nrow(object1$aie2)){
       j1=match(rownames(object1$aie2)[j],x.names)
       j2=grep(rownames(object1$aie2)[j],x.names.2)
       
       if(!object1$data1$binx[j1])
       {par(mfcol=c(3,length(a2)))
         for(i in 1:length(a2))
         {ie2<-boot.ci(object1$x.j[,j2],matrix(object$ie2[,a2[i],j],ncol=object$boot),cri_val)
          plot_ci(ie2,xlab=x.names.2[j2],ylab=paste("Level 2 IE of", colnames(object1$aie2)[a2[i]],sep=" "),
                 main=paste("IE with exposure variable", rownames(object1$aie2)[j],sep=" "))
         ie2_2<-boot.ci(object1$x.j[,j2],matrix(object$ie2_2[,a2[i],j],ncol=object$boot),cri_val)
         suppressWarnings(plot_ci(ie2_2,xlab=x.names.2[j2],ylab=colnames(object$aie2)[a2[i]],
                                  main=paste("Differential effect of", rownames(object$aie2)[j], "and",
                                             colnames(object$ie2)[a2[i]],sep=" ")))
         temp.matrix=apply(matrix(object$ie2_1[,a2[i]],ncol=object$boot),2,two,object$level)
         if(!cate)
         {
          ie2_1<-boot.ci(m.j,temp.matrix,cri_val)
          plot_ci(ie2_1,xlab=var,ylab="y",main=paste("Predicted relationship between y and",var,sep=" "))
         }
         else 
           boxplot(temp.matrix[1,], ylab="Fitted Coeficients",
                   xlab=dimnames(object$ie2)[[2]][a2[i]], 
                   main=paste("Predicted coefficients for \n",colnames(object$aie2)[a2[i]], "to predict y", sep=" ")) 
         }}
       else
       {par(mfcol=c(2,1))
         if(!cate)
           plot(m.j~as.factor(object1$x.j[,j2]),xlab=rownames(object1$aie2)[j],ylab=var,
                main=paste("Relationship between",x.names.2[j2],"and",var,sep=" ")) 
         else
           plot(as.factor(m.j)~as.factor(object1$x.j[,j2]),xlab=rownames(object1$aie2)[j],ylab=var,
                main=paste("Relationship between",x.names.2[j2],"and",var,sep=" "))
         if(!cate)
         {temp.matrix=apply(matrix(object$ie2_1[,a2],ncol=object$boot),2,two,object$level)
         ie2_1<-boot.ci(m.j,temp.matrix,cri_val)
         plot_ci(ie2_1,xlab=var,ylab="y",main=paste("Predicted relationship between y and",var,sep=" "))
         }
         else
           plot(y.j~as.factor(m.j), ylab="y",xlab=var,
                main=paste("Predicted relationship \n between y and",colnames(object$aie2)[i],sep=" ")) 
       }} }  
else
{ temp=grep(var,colnames(object$m))
  a1=object$coef.fm1[[1]]%in%temp
  a2=object$coef.fm2[[1]]%in%temp
  if(sum(a1)>0)
  { par(mfcol=c(sum(a1)+1,1))
    temp=grep(var,colnames(object$coef.f1))
    boxplot(object$coef.f1[,temp],main=paste("Coefficients of", var, "at the final model"))  
    temp1=colnames(object$coef.f1)[temp]
    a1=(1:length(a1))[a1]
    for(i in 1:length(a1))
      boxplot(object$coef.fm1[[a1[i]+1]], 
              main=paste("Coefficients of exposures in predicting",temp1[i]))
  }
  if(sum(a2)>0)
  { par(mfcol=c(sum(a2)+1,1))
    temp=grep(var,colnames(object$coef.f1))
    boxplot(object$coef.f1[,temp],main=paste("Coefficients of", var, "at the final model"))  
    temp1=colnames(object$coef.f1)[temp]
    a2=(1:length(a2))[a2]
    for(i in 1:length(a2))
      boxplot(object$coef.fm2[[a2[i]+1]], 
              main=paste("Coefficients of exposures in predicting",temp1[i]))
  }
}
  #plot(object$full, var=var)

}


joint.effect<-function(object,var.list,...,alpha=0.05, echo=FALSE)
{summarize.boot<-function(vector,n,a1,a2,b1,b2)  #vector is stacked results from bootstrap samples
  #n is the number of elements from each bootstrap sample
{mat<-matrix(vector,length(vector)/n,n)
all<-summarize.boot1(t(mat),a1,a2,b1,b2)
data.frame(all)
}
summarize.boot1<-function(mat,a1,a2,b1,b2){
  mean.temp=apply(mat,2,mean,na.rm=T)
  sd.temp=apply(mat,2,sd,na.rm=T)
  cbind(mean=mean.temp,sd=sd.temp,upbd=mean.temp+b2*sd.temp,
        lwbd=mean.temp+b1*sd.temp,upbd.quan=apply(mat,2,quantile,a2,na.rm=T),
        lwbd.quan=apply(mat,2,quantile,a1,na.rm=T))
}

a1<-alpha/2
a2<-1-a1
b1<-qnorm(a1)
b2<-qnorm(a2)
boot=object$boot

aie1=object$aie1
aie2=cbind(object$aie2,object$aie12)

m1.names=colnames(aie1)
m2.names=colnames(aie2)

if(!is.character(var.list))
  var.list=colnames(object$m)[var.list]

where1=c(unlist(sapply(var.list,grep,m1.names)))
where2=c(unlist(sapply(var.list,grep,m2.names)))

indirect.effect2=NULL
indirect.effect1=NULL
re1=NULL
re2=NULL

if(length(where2)!=0)
{temp=cbind(object$full$aie2,object$full$aie12)
 if(nrow(temp)==1)
   temp=sum(temp[,where2])
 else
  temp=apply(as.matrix(cbind(object$full$aie2,object$full$aie12)[,where2]),1,sum)
 ie2.1=apply(as.matrix(aie2[,where2]),1,sum)
 indirect.effect2<-cbind(est=temp,summarize.boot(ie2.1,boot,a1,a2,b1,b2))
 rownames(indirect.effect2)=rownames(object$full$aie2)
 re2<-cbind(est=temp/object$full$ate2,summarize.boot(ie2.1/as.vector(t(object$ate2)),boot,a1,a2,b1,b2))
}
if(length(where1)!=0)
{if(nrow(object$full$aie1)==1)
  temp=sum(object$full$aie1[,where1])
 else
  temp=apply(as.matrix(object$full$aie1[,where1]),1,sum)
 ie1.1=apply(as.matrix(aie1[,where1]),1,sum)
 indirect.effect1<-cbind(est=temp,summarize.boot(ie1.1,boot,a1,a2,b1,b2))
 rownames(indirect.effect1)=rownames(object$full$aie1)
 re1<-cbind(est=temp/object$full$ate1,summarize.boot(ie1.1/as.vector(t(object$ate1)),boot,a1,a2,b1,b2))
}

if(echo){
cat("The joint indirect effect:\n")
print(rbind(indirect.effect1,indirect.effect2))
cat("The joint relative effect:\n")
print(rbind(re1,re2))
}

list(indirect.effect1=indirect.effect1, indirect.effect2=indirect.effect2,
     re1=re1,re2=re2)
}


