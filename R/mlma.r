# data clean
data.org<-function(x, levelx, m, l1=NULL,l2=NULL, c1=NULL, c1r=rep(1,length(c1)), #levelx is the level of x
                   c2=NULL, c2r=rep(1,length(c2)), f01y=NULL, f10y=NULL,                  #level is the level of observations
                   f02ky=NULL, f20ky=NULL, f01km1=NULL, f01km2=NULL, f10km=NULL,          #weight is the level 1 weight of cases
                   level, weight=rep(1,length(x)))                                        #weight2 is the level 2 weight of cases, weight2=rep(1,length(unique(level[!is.na(level)])))   
{x2fx<-function(x,func) #x is the list of original numerical vector, func is a vector of character functions. 
{ # eg.  func <- c("x","x+1","x+2","x+3","log(x)")
  func.list <- list()
  test.data <- matrix(data=rep(x,length(func)),length(x),length(func))
  test.data <- data.frame(test.data)
  for(i in 1:length(func)){
    func.list[[i]] <- function(x){}
    body(func.list[[i]]) <- parse(text=func[i])
  }
  result <- mapply(do.call,func.list,lapply(test.data,list))
  as.matrix(result)
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

two<-function(x, level, weight=rep(1,nrow(as.matrix(x))))
{x<-as.matrix(x)
 levels<-unique(level[!is.na(level)])
 x2<-matrix(NA,length(levels),dim(x)[2])
 for(i in 1:length(levels))
 {cho<-(level==levels[i])
  if(sum(cho)>0)
  {temp<-as.matrix(x[level==levels[i],])
   weight1<-weight[level==levels[i]]
   x2[i,]<-apply(temp,2,weighted.mean,weight1,na.rm=TRUE)}}
 colnames(x2)<-colnames(x)
 x2
}

x2fdx<-function(x,func)  #x is the list of original numerical vector, func is a vector of character functions. 
{ fdx<-NULL              # eg.  func <- c("x","x+1","x+2","x+3","log(x)")
  for(i in 1:length(func)){
    dx2x <- D(parse(text=func[i]), "x") 
    temp<-eval(dx2x)
    if(length(temp)==1)
      fdx<-cbind(fdx,rep(temp,length(x)))
    else fdx<-cbind(fdx,temp)
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
 temp3<-z[1:dim1[1],]
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
if(levelx==2)
{l1x<-NULL               #create level 1 and level 2 x variables that are used to predict y
 if(is.null(f01y))
 {x1<-x
  l2x<-1
  x1.der<-rep(1,n)}
 else
 {x1<-x2fx(x,f01y)
  x1.der<-x2fdx(x,f01y)
  l2x<-1:length(f01y)
  x1<-as.matrix(x1)
  colnames(x1)<-paste("x",1:length(f01y),sep=".")
  x1.der<-as.matrix(x1.der)
  colnames(x1.der)<-paste("x",1:length(f01y),sep=".")}
 
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
     b<-(1:length(l2))[l2==a]
     d<-as.matrix(x2fx(m[,a],f02ky[[i]]))
     d.der<-as.matrix(x2fdx(m[,a],f02ky[[i]]))
     colnames(d)<-paste(mnames[a],1:length(f02ky[[i]]),sep=".")
     colnames(d.der)<-paste(mnames[a],1:length(f02ky[[i]]),sep=".")
     if(dim(as.matrix(d))[2]==1)
     {m2y[,b]<-d
      m2y.der[,b]<-d.der}
     else {
       m2y[,b]<-d[,1]
       m2y.der[,b]<-d.der[,1]
       m2[[b+1]]<-c(m2[[b+1]],(dim(as.matrix(m2y))[2]+1):(dim(as.matrix(m2y))[2]+length(f02ky[[i]])-1))
       temp.m2<-c(colnames(m2y),colnames(d)[-1])
       m2y<-cbind(m2y,d[,-1])
       colnames(m2y)<-temp.m2
       m2y.der<-cbind(m2y.der,d.der[,-1])
       colnames(m2y.der)<-temp.m2
     }
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
  m2y.der<-temp$m2y.der
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
  if(!is.null(f20ky)) 
    for (i in 2:length(f20ky))
    {a<-f20ky[[1]][i-1]
     b<-(1:length(l1))[l1==a]
     d<-as.matrix(x2fx(m[,a],f20ky[[i]]))
     d.der<-as.matrix(x2fdx(m[,a],f20ky[[i]]))
     colnames(d)<-paste(mnames[a],1:length(f20ky[[i]]),sep=".")
     colnames(d.der)<-paste(mnames[a],1:length(f20ky[[i]]),sep=".")
     if(is.null(dim(d)))
     {temp<-one2two(d,level,weight)
      m1y[,b]<-temp[,2]
      m12y[,b]<-temp[,1]
      temp.der<-one2two(d.der,level,weight)
      m1y.der[,b]<-temp.der[,2]
      m12y.der[,b]<-temp.der[,1]}
     else {
       temp<-apply(as.matrix(d), 2, one2two, level, weight)
       temp.der<-apply(as.matrix(d.der), 2, one2two, level, weight)
       temp.namem1y<-c(colnames(m1y),colnames(d)[-1])
       temp.namem12y<-c(colnames(m12y),paste(mnames[a],"12", 2:length(f20ky[[i]]),sep="."))
       m12y[,b]<-temp[1:n,1]
       m1y[,b]<-temp[(n+1):(2*n),1]
       m12y.der[,b]<-temp.der[1:n,1]
       m1y.der[,b]<-temp.der[(n+1):(2*n),1]
       m1[[b+1]]<-c(m1[[b+1]],(dim(as.matrix(m1y))[2]+1):(dim(as.matrix(m1y))[2]+length(f20ky[[i]])-1))
       m1y<-cbind(m1y,temp[(n+1):(2*n),-1])
       m1y.der<-cbind(m1y.der,temp.der[(n+1):(2*n),-1])
       m12[[b+1]]<-c(m12[[b+1]],(dim(as.matrix(m12y))[2]+1+m2yd):(dim(as.matrix(m12y))[2]+length(f20ky[[i]])-1+m2yd))
       m12y<-cbind(m12y,temp[1:n,-1])
       m12y.der<-cbind(m12y.der,temp.der[1:n,-1])
       colnames(m1y)<-temp.namem1y
       colnames(m12y)<-temp.namem12y
       colnames(m1y.der)<-temp.namem1y
       colnames(m12y.der)<-temp.namem12y
     }
    }
  m2y<-cbind(m2y,m12y)
  m2y.der<-cbind(m2y.der,m12y.der)
 }
 else 
 {m1y<-NULL
  m1y.der<-NULL
  m1<-NULL
  m12<-NULL}
 
 if(!is.null(c1))                  #binarize level 1 categorical variables and corresponding level 2 variables to explain y 
 {temp<-cattobin(m1y,m1y.der,m1,m,c1,c1r,level,weight)
  m1y<-temp$m1y
  m1y.der<-temp$m1y.der
  m1<-temp$m1
  dim2y<-dim(m2y)
  m2y<-cbind(m2y,temp$m12y)
  m2y.der<-cbind(m2y.der,temp$m12y.der)
  for (i in 2:length(temp$m12))
    m12<-append(m12,list(temp$m12[[i]]+dim2y[2]))
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
   temp.xm2<-c(colnames(xm2),paste("x.2",1:length(unifun1),sep="."))
   xm2.der<-cbind(xm2.der,x2fdx(xm2,unifun1))
   xm2<-cbind(xm2,x2fx(xm2,unifun1))
   colnames(xm2)<-temp.xm2
   colnames(xm2.der)<-temp.xm2
   for(i in 2:length(f01km2))
     fm22[[(1:length(lc2))[fm22[[1]]==f01km2[[1]][i-1]]+1]]<-order_char(unifun2,f01km2[[i]])
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
   xm1<-cbind(x,x2fx(x,unifun1))
   xm1.der<-cbind(xm1.der,x2fdx(x,unifun1))
   colnames(xm1)<-paste("x2",1:length(unifun2),sep=".")
   colnames(xm1.der)<-paste("x2",1:length(unifun2),sep=".")
   for(i in 2:length(f01km1))
     fm12[[(1:length(lc1))[fm12[[1]]==f01km1[[1]][i-1]]+1]]<-order_char(unifun2,f01km1[[i]])
  }}
 else
 {xm1<-NULL
  xm1.der<-NULL
  fm12<-NULL}
}
else                          # to deal with level 1 x
{x2<-one2two(x,level,weight)               #create level 1 and level 2 x variables that are used to predict y
 if(!is.null(f01y))
 {x1<-x2fx(x2[,1],f01y)
  x1.der<-x2fdx(x2[,1],f01y)
  l2x<-1:length(f01y)
  colnames(x1)<-paste("x",1:length(f01y),sep=".")
  colnames(x1.der)<-paste("x",1:length(f01y),sep=".")}
 else
 {x1<-NULL
  x1.der<-NULL
  l2x<-NULL}
 if(!is.null(f10y))
 {x3<-x2fx(x,f10y)
  x4<-apply(as.matrix(x3),2,one2two,level,weight)
  x3.der<-x2fdx(x,f10y)
  x4.der<-apply(as.matrix(x3.der),2,one2two,level,weight)
  l2x<-c(l2x,(length(l2x)+1):(length(l2x)+length(f10y)))
  l1x<-(length(l2x)+1):(length(l2x)+length(f10y))
  temp.name<-c(colnames(x1),paste("x.j",1:length(f10y),sep="."),paste("xij",1:length(f10y),sep="."))
  x1<-cbind(x1,x4[1:n,],x4[(n+1):(2*n),])
  x1.der<-cbind(x1.der,x4.der[1:n,],x4.der[(n+1):(2*n),])
  colnames(x1)<-temp.name
  colnames(x1.der)<-temp.name}
 else
 {l2x<-c(l2x,length(l2x)+1)
  l1x<-length(l2x)+1
  temp.name<-c(colnames(x1),"x.j","xij")
  x1<-cbind(x1,x2)
  colnames(x1)<-temp.name
  x1.der<-cbind(x1.der,matrix(1,n,2))
  colnames(x1.der)<-temp.name}
 
 if(!(is.null(l2) & is.null(c2)))
   stop("x is a level 1 variable, m cannot be level 2: l2 and c2 should be empty!")
 
 m2<-NULL
 if(!is.null(l1))         #create level 1 and corresponding level 2 m variables to explain y
 {temp<-apply(as.matrix(m[,l1]), 2, one2two, level, weight)
  m1y<-as.matrix(temp[(n+1):(2*n),])
  m1y.der<-matrix(1,n,ncol(m1y))
  colnames(m1y)<-mnames[l1]
  colnames(m1y.der)<-mnames[l1]
  m2y<-temp[1:n,]
  m2y.der<-matrix(1,n,ncol(m2y))
  colnames(m2y)<-paste(mnames[l1],"12",sep=".")
  colnames(m2y.der)<-paste(mnames[l1],"12",sep=".")
  m1<-as.list(c(1,1:length(l1)))
  m12<-as.list(c(1,1:length(l1)))
  m1[[1]]<-l1
  m12[[1]]<-l1
  if(!is.null(f20ky)) 
    for (i in 2:length(f20ky))
    {a<-f20ky[[1]][i-1]
     b<-(1:length(l1))[l1==a]
     d<-x2fx(m[,a],f20ky[[i]])
     d.der<-x2fdx(m[,a],f20ky[[i]])
     colnames(d)<-paste(mnames[a],1:length(f20ky),sep=".")
     if(is.null(dim(d)))
     {temp<-one2two(d,level,weight)
      m1y[,b]<-temp[,2]
      m2y[,b]<-temp[,1]
      temp.der<-one2two(d.der,level,weight)
      m1y.der[,b]<-temp.der[,2]
      m2y.der[,b]<-temp.der[,1]}
     else {
       temp<-apply(as.matrix(d), 2, one2two, level, weight)
       m2y[,b]<-temp[1:n,1]
       m1y[,b]<-temp[(n+1):(2*n),1]
       temp.der<-apply(as.matrix(d.der), 2, one2two, level, weight)
       m2y.der[,b]<-temp.der[1:n,1]
       m1y.der[,b]<-temp.der[(n+1):(2*n),1]
       tempname.m1y<-c(colnames(m1y),colnames(d)[-1])
       tempname.m2y<-c(colnames(m2y),paste(mnames[a],"12",1:length(f20ky),sep=".")[-1])
       m1[[b+1]]<-c(m1[[b+1]],(dim(as.matrix(m1y))[2]+1):(dim(as.matrix(m1y))[2]+length(f20ky[[i]])-1))
       m1y<-cbind(m1y,temp[(n+1):(2*n),-1])
       m1y.der<-cbind(m1y.der,temp.der[(n+1):(2*n),-1])
       m12[[b+1]]<-c(m12[[b+1]],(dim(as.matrix(m2y))[2]+1):(dim(as.matrix(m2y))[2]+length(f20ky[[i]])-1))
       m2y<-cbind(m2y,temp[1:n,-1])
       m2y.der<-cbind(m2y.der,temp.der[1:n,-1])
       colnames(m1y)<-tempname.m1y
       colnames(m2y)<-tempname.m2y
       colnames(m1y.der)<-tempname.m1y
       colnames(m2y.der)<-tempname.m2y
     }
    }
 }
 else 
 {m1y<-NULL
  m1y.der<-NULL
  m1<-NULL
  m12<-NULL
  m2y<-NULL   
  m2y.der<-NULL}
 
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
 
 fm22<-NULL                        #fm11 is the list of level 1 x to explain level 1 m 
 m.2<-NULL
 xm2<-NULL
 xm2.der<-NULL
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
   xm1<-cbind(x2[,1],x2fx(x2[,1],unifun1))
   xm1.der<-cbind(rep(1,n),x2fdx(x2[,1],unifun1))
   colnames(xm1)<-c("x.2",paste("x.2",1:length(unifun1),sep="."))
   colnames(xm1.der)<-c("x.2",paste("x.2",1:length(unifun1),sep="."))
   for(i in 2:length(f01km1))
     fm12[[(1:length(lc1))[fm12[[1]]==f01km1[[1]][i-1]]+1]]<-order_char(unifun2,f01km1[[i]])
  }
  
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
   temp2<-x2fx(x2[,2],unifun1)
   temp3<-apply(temp2,2,one2two,level,weight)
   temp2.der<-x2fdx(x2[,2],unifun1)
   temp3.der<-apply(temp2.der,2,one2two,level,weight)
   for(i in 2:length(f10km))
     fm12[[(1:length(lc1))[fm12[[1]]==f10km[[1]][i-1]]+1]]<-c(fm12[[(1:length(lc1))[fm12[[1]]==f10km[[1]][i-1]]+1]],
                                                              order_char(unifun2,f10km[[i]])+dimxm1[2])
   dimxm1<-c(dimxm1[1],dimxm1[2]+dim(temp3)[2])
   for(i in 2:length(f10km))
     fm11[[(1:length(lc1))[fm11[[1]]==f10km[[1]][i-1]]+1]]<-order_char(unifun2,f10km[[i]])+dimxm1[2]
   temp.xm1<-c(colnames(xm1),"x1",paste("x.12",1:length(unifun1),sep="."),paste("x1",1:length(unifun1),sep="."))
   xm1<-cbind(xm1,x2[,2],temp3[1:dimxm1[1],],temp3[(dimxm1[1]+1):(2*dimxm1[1]),])
   xm1.der<-cbind(xm1.der,rep(1,n),temp3.der[1:dimxm1[1],],temp3.der[(dimxm1[1]+1):(2*dimxm1[1]),])
   colnames(xm1)<-temp.xm1
   colnames(xm1.der)<-temp.xm1
  }
 }
 else
 {xm1<-NULL
  xm1.der<-NULL
  fm12<-NULL
  fm11<-NULL}
}

list(x1=x1, x1.der=x1.der, l1x=l1x, l2x=l2x, m1y=m1y, m1y.der=m1y.der,
     m1=m1, m2y=m2y, m2y.der=m2y.der, m2=m2, m12=m12, xm1=xm1, xm1.der=xm1.der,
     fm11=fm11, fm12=fm12, m.2=m.2, xm2=xm2, xm2.der=xm2.der, fm22=fm22)
}

#multilevel mediation analysis
mlma<-function(y, biny=FALSE, data1, x, levelx, m, l1=NULL,l2=NULL, c1=NULL, c1r=rep(1,length(c1)), #levelx is the level of x
               c2=NULL, c2r=rep(1,length(c2)), level, weight=rep(1,length(x)), 
               random="(1|level)", random.m1=NULL,intercept=TRUE)                               
  
{two<-function(x, level, weight=rep(1,nrow(as.matrix(x))))
{x<-as.matrix(x)
 levels<-unique(level[!is.na(level)])
 x2<-matrix(NA,length(levels),dim(x)[2])
 for(i in 1:length(levels))
 {cho<-(level==levels[i])
  if(sum(cho)>0)
  {temp<-as.matrix(x[level==levels[i],])
   weight1<-weight[level==levels[i]]
   x2[i,]<-apply(temp,2,weighted.mean,weight1,na.rm=TRUE)}}
 colnames(x2)<-colnames(x)
 x2
}

getformula<-function(expl,random="(1|level)",intercept=TRUE)
{temp.name<-colnames(expl)
 formula<-temp.name[1]
 if(ncol(expl)>1)
   for (i in 2:ncol(expl))
     formula<-paste(formula,temp.name[i],sep="+")
 formula<-ifelse(intercept,paste("1",formula,sep="+"),formula)
 formula<-paste("y~",formula)
 if (!is.null(random))
   formula<-paste(formula,random,sep="+")
 formula
}

expl<-cbind(data1$x1[,data1$l2x],data1$m2y,data1$x1[,data1$l1x],data1$m1y)
frml<-getformula(expl,random,intercept)
temp.data<-cbind(y=y,level=level,data1$x1,data1$m1y,data1$m2y)
if(biny)
  f1<-glmer(frml,data=data.frame(temp.data),family=binomial(link="logit"))
else
  f1<-lmer(frml,data=data.frame(temp.data))
lc1<-c(l1,c1)
lc2<-c(l2,c2)

coef.f1<-summary(f1)$coefficient[-1,1] #cbind(data1$x1[,data1$l2x],data1$m2y,data1$x1[,data1$l1x],data1$m1y)
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
    for (k in 1:length(data1$m2[[temp+1]]))
      ie2_1<-cbind(ie2_1,data1$m2y.der[,data1$m2[[temp+1]]][,k]*coef.f1[len[1]+data1$m2[[temp+1]][k]])
    colnames(ie2_1)<-c(temp.name,paste(colnames(m)[c2[i]],1:k,sep="."))
    ie2_list[[1]]<-c(ie2_list[[1]],colnames(m)[c2[i]])
    ie2_list<-append(ie2_list,list((z+1):(z+k)))
    z<-z+k
   }  }        
 if (!is.null(data1$m12))
 {if(!is.null(l1))
   for (k in 1:length(l1))
   {if (length(data1$m12[[k+1]])==1)
     ie12_1<-cbind(ie12_1,data1$m2y.der[,data1$m12[[k+1]]]*coef.f1[len[1]+data1$m12[[k+1]]])
    else
      ie12_1<-cbind(ie12_1,data1$m2y.der[,data1$m12[[k+1]]]%*%coef.f1[len[1]+data1$m12[[k+1]]])}
  colnames(ie12_1)<-paste(colnames(m)[l1],"12",sep=".")
  if(!is.null(c1))
    for (i in 1:length(c1)) 
    {temp.name<-colnames(ie12_1)
     temp<-(1:length(data1$m12[[1]]))[data1$m12[[1]]==c1[i]]
     for (k in 1:length(data1$m12[[temp+1]]))
       ie12_1<-cbind(ie12_1,data1$m2y.der[,data1$m12[[temp+1]]][,k]*coef.f1[len[1]+data1$m12[[temp+1]][k]])
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
    for (k in 1:length(data1$m1[[temp+1]]))
      ie1_1<-cbind(ie1_1,data1$m1y.der[,data1$m1[[temp+1]]][,k]*coef.f1[len[1]+len[2]+len[3]+data1$m1[[temp+1]][k]])
    colnames(ie1_1)<-c(temp.name,paste(colnames(m)[c1[i]],1:k,sep="."))
    ie1_list[[1]]<-c(ie1_list[[1]],colnames(m)[c1[i]])
    ie1_list<-append(ie1_list,list((z+1):(z+k)))
    z<-z+k
   }  } }



if(levelx==2)    #analysis when x is a level 2 variable
{fm1<-NULL            #models for x to explain level 1 mediators
 ie12_2<-NULL
 if(!is.null(l1))
 {fm1[[1]]<-l1
  for (i in 1:length(l1))
  {expl.m<-as.matrix(data1$xm1[,data1$fm12[[i+1]]])
   colnames(expl.m)<-colnames(data1$xm1)[data1$fm12[[i+1]]]
   frml.m<-getformula(expl.m,random,intercept)
   temp.data<-cbind(y=m[,l1[i]],level=level,data1$xm1[,data1$fm12[[i+1]]])
   colnames(temp.data)<-c("y","level",colnames(data1$xm1)[data1$fm12[[i+1]]])
   model<-lmer(frml.m,data=data.frame(temp.data))
   fm1<-append(fm1,list(model))
   coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
   if(intercept)
     coef.temp<-coef.temp[-1]
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
   {expl.m<-as.matrix(data1$xm1[,data1$fm12[[j+1]]])
    colnames(expl.m)<-colnames(data1$xm1)[data1$fm12[[j+1]]]
    frml.m<-getformula(expl.m,random,intercept)
    temp<-(1:length(data1$m1[[1]]))[data1$m1[[1]]==c1[i]]
    name.temp<-colnames(ie12_2)
    for (k in data1$m1[[temp+1]])
    {temp.data<-cbind(y=data1$m1y[,k],level=level,data1$xm1[,data1$fm12[[j+1]]])
     colnames(temp.data)<-c("y","level",colnames(data1$xm1)[data1$fm12[[j+1]]])
     model<-glmer(frml.m,data=data.frame(temp.data),
                  family=binomial(link="logit"))
     p.temp<-predict(model,type="response")
     fm1<-append(fm1,list(model))
     fm1[[1]]<-c(fm1[[1]],c1[i])
     coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
     if(intercept)
       coef.temp<-coef.temp[-1]
     if(length(coef.temp)==1)
       ie12_2<-cbind(ie12_2,p.temp*(1-p.temp)*(coef.temp*data1$xm1.der[,data1$fm12[[j+1]]]))
     else
       ie12_2<-cbind(ie12_2,p.temp*(1-p.temp)*
                       (data1$xm1.der[,data1$fm12[[j+1]]]%*%coef.temp))
    }
    colnames(ie12_2)<-c(name.temp,paste(colnames(m)[c1[i]],
                                        1:length(data1$m1[[temp+1]]),sep="."))
    j<-j+1
   }  
 fm2<-NULL            #models for x to explain level 2 mediators
 if(!is.null(l2))
 {fm2[[1]]<-l2
  ie2_2<-NULL
  for (i in 1:length(l2))
  {expl.m<-as.matrix(data1$xm2[,data1$fm22[[i+1]]])
   colnames(expl.m)<-colnames(data1$xm2)[data1$fm22[[i+1]]]
   frml.m<-getformula(expl.m,random=NULL,intercept)
   temp.data<-cbind(y=data1$m.2[,i],data1$xm2[,data1$fm22[[i+1]]])
   colnames(temp.data)<-c("y",colnames(data1$xm2)[data1$fm22[[i+1]]])
   model<-lm(frml.m,data=data.frame(temp.data))
   fm2<-append(fm2,list(model))
   coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
   if(intercept)
     coef.temp<-coef.temp[-1]
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
   {expl.m<-as.matrix(data1$xm2[,data1$fm22[[j+1]]])
    colnames(expl.m)<-colnames(data1$xm2)[data1$fm22[[j+1]]]
    frml.m<-getformula(expl.m,random=NULL,intercept)
    name.temp<-colnames(ie2_2)
    temp<-(1:length(data1$m2[[1]]))[data1$m2[[1]]==c2[i]]
    for (k in data1$m2[[temp+1]])
    {temp.data<-cbind(y=two(data1$m2y[,k],level),data1$xm2[,data1$fm22[[j+1]]])
     colnames(temp.data)<-c("y",colnames(data1$xm2)[data1$fm22[[j+1]]])
     model<-glm(frml.m,data=data.frame(temp.data),
                family=binomial(link="logit"))
     p.temp<-predict(model,type="response")
     fm2<-append(fm2,list(model))
     fm2[[1]]<-c(fm2[[1]],c2[i])
     coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
     if(intercept)
       coef.temp<-coef.temp[-1]
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
{fm1<-NULL            #models for x to explain level 1 mediators
 if(!is.null(l1))
 {fm1[[1]]<-l1
  ie1_2<-NULL
  ie12_2<-NULL
  for (i in 1:length(l1))
  {expl.m<-as.matrix(data1$xm1[,c(data1$fm12[[i+1]],data1$fm11[[i+1]])])
   colnames(expl.m)<-colnames(data1$xm1)[c(data1$fm12[[i+1]],data1$fm11[[i+1]])]
   m.random<-"(1|level)"
   if(!is.null(random.m1))
     if(sum(random.m1[[1]]==l1[i])>0)
       m.random<-random.m1[[2]][random.m1[[1]]==l1[i]] #1st item of random.m1 is the list of l1 med, 2nd item is the random item of the same order
   frml.m<-getformula(expl.m,m.random,intercept)
   temp.data<-cbind(y=m[,l1[i]],level=level,data1$xm1[,c(data1$fm12[[i+1]],data1$fm11[[i+1]])])
   colnames(temp.data)<-c("y","level",colnames(data1$xm1)[c(data1$fm12[[i+1]],data1$fm11[[i+1]])])
   model=lmer(frml.m,data=data.frame(temp.data))
   fm1<-append(fm1,list(model))
   coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
   if(intercept)
     coef.temp<-coef.temp[-1]
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
   {expl.m<-as.matrix(data1$xm1[,c(data1$fm12[[j+1]],data1$fm11[[j+1]])])
    colnames(expl.m)<-colnames(data1$xm1)[c(data1$fm12[[j+1]],data1$fm11[[j+1]])]
    m.random<-"(1|level)"
    if(!is.null(random.m1))
      if(sum(random.m1[[1]]==c1[i])>0)
        m.random<-random.m1[[2]][random.m1[[1]]==c1[i]] #1st item of random.m1 is the list of l1 med, 2nd item is the random item of the same order
    frml.m<-getformula(expl.m,m.random,intercept)
    name.temp<-colnames(ie12_2)
    temp<-(1:length(data1$m1[[1]]))[data1$m1[[1]]==c1[i]]
    for (k in data1$m1[[temp+1]])
    {temp.data<-cbind(y=data1$m1y[,k],level=level,data1$xm1[,c(data1$fm12[[j+1]],data1$fm11[[j+1]])])
     colnames(temp.data)<-c("y","level",colnames(data1$xm1)[c(data1$fm12[[j+1]],data1$fm11[[j+1]])])
     model=glmer(frml.m,data=data.frame(temp.data),
                 family=binomial(link="logit"))
     fm1<-append(fm1,list(model))
     fm1[[1]]<-c(fm1[[1]],c1[i])
     p.temp<-predict(model,type="response")
     coef.temp<-summary(model)$coefficient[,1]   #find the second part of IE2
     if(intercept)
       coef.temp<-coef.temp[-1]
     coef.temp12<-coef.temp[1:length(data1$fm12[[j+1]])]
     coef.temp1<-coef.temp[-(1:length(data1$fm12[[j+1]]))]
     if(length(coef.temp12)==1)
       ie12_2<-cbind(ie12_2,coef.temp12*data1$xm1.der[,data1$fm12[[j+1]]])
     else
       ie12_2<-cbind(ie12_2,data1$xm1.der[,data1$fm12[[j+1]]]%*%coef.temp12)
     if(length(coef.temp1)==1)
       ie1_2<-cbind(ie1_2,p.temp*(1-p.temp)*coef.temp1*data1$xm1.der[,data1$fm11[[j+1]]])
     else
       ie1_2<-cbind(ie1_2,p.temp*(1-p.temp)*(data1$xm1.der[,data1$fm11[[j+1]]]%*%coef.temp1))
    }
    colnames(ie12_2)<-c(name.temp,paste(colnames(m)[c1[i]],
                                        1:length(data1$m1[[temp+1]]),sep="."))
    colnames(ie1_2)<-c(name.temp,paste(colnames(m)[c1[i]],
                                       1:length(data1$m1[[temp+1]]),sep="."))
    j<-j+1
   }
 fm2<-NULL
 ie2_2<-NULL
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
a<-list(de1=DE1,de2=two(DE2,level),ie1=ie1,ie2=ie2,ie12=ie12,
        f1=f1,fm1=fm1,fm2=fm2, ie1_list=ie1_list, ie2_list=ie2_list,
        ie12_1=ie12_1, ie12_2=ie12_2, ie1_1=ie1_1,ie1_2=ie1_2,
        ie2_1=ie2_1,ie2_2=ie2_2,x=x, x.j=two(x,level), m=m, 
        level=level, levelx=levelx, weight=weight)
class(a)<-"mlma"
return(a)
}

boot.mlma<-function(y, biny=FALSE, x, levelx, m, l1=NULL,l2=NULL, c1=NULL, 
                    c1r=rep(1,length(c1)), #levelx is the level of x
                    c2=NULL, c2r=rep(1,length(c2)), f01y=NULL, f10y=NULL,                  #level is the level of observations
                    f02ky=NULL, f20ky=NULL, f01km1=NULL, f01km2=NULL, f10km=NULL,          #weight is the level 1 weight of cases
                    level, weight=rep(1,length(x)), random="(1|level)",
                    random.m1=NULL, intercept=TRUE, boot=100,seed=1)
{two<-function(x, level, weight=rep(1,nrow(as.matrix(x))))
{x<-as.matrix(x)
 levels<-unique(level[!is.na(level)])
 x2<-matrix(NA,length(levels),dim(x)[2])
 for(i in 1:length(levels))
 {cho<-(level==levels[i])
  if(sum(cho)>0)
  {temp<-as.matrix(x[level==levels[i],])
   weight1<-weight[level==levels[i]]
   x2[i,]<-apply(temp,2,weighted.mean,weight1,na.rm=TRUE)}}
 colnames(x2)<-colnames(x)
 x2
}

data1<-data.org(x, levelx, m, l1, l2, c1, c1r, c2, c2r, f01y, f10y, f02ky, f20ky,
                f01km1, f01km2, f10km, level, weight)
full<-mlma(y, biny, data1, x, levelx, m, l1,l2, c1, c1r, 
           c2, c2r, level, weight, random, random.m1,intercept)                               
de1<-NULL
de2<-NULL
ie1<-NULL
ie2<-NULL
ie12<-NULL
xboot<-NULL
xjboot<-NULL

m<-as.matrix(m)
t1<-table(level)
t2<-sort(unique(level))
level.boot<-rep(t2,t1)
for(i in 1:boot)
{set.seed(seed+i)
 x.boot<-NULL
 y.boot<-NULL
 m.boot<-NULL
 weight.boot<-NULL
 for(j in 1:length(t2))
 {temp<-level==t2[j]
  temp.2<-sample(t1[j],replace=TRUE,prob=weight[temp])
  x.boot<-c(x.boot,x[temp][temp.2])
  y.boot<-c(y.boot,y[temp][temp.2]) 
  m.boot<-rbind(m.boot,m[temp,][temp.2,]) 
  weight.boot<-c(weight.boot,weight[temp][temp.2])
 }
 data1<-data.org(x.boot, levelx, m.boot, l1, l2, c1, c1r, c2, c2r, f01y, f10y, 
                 f02ky, f20ky, f01km1, f01km2, f10km, level.boot, weight.boot)
 temp<-mlma(y.boot, biny, data1, x.boot, levelx, m.boot, l1,l2, c1, c1r, 
            c2, c2r, level.boot, weight.boot, random, random.m1,intercept)                               
 de1<-cbind(de1,temp$de1)
 de2<-cbind(de2,temp$de2)
 ie1<-rbind(ie1,temp$ie1)
 if(!is.null(temp$ie2))
   ie2<-rbind(ie2,as.matrix(temp$ie2))
 if(!is.null(full$ie12))
   ie12<-rbind(ie12,temp$ie12)
 xboot<-c(xboot,x.boot)
 xjboot<-c(xjboot,two(x.boot,level.boot))
 print(i)
}
a<-list(de1=de1,de2=de2,ie1=ie1,ie2=ie2,ie12=ie12, full=full,
        xboot=xboot, xjboot=xjboot,levelx=levelx, level=level)
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
  {temp<-as.matrix(x[level==levels[i],])
   weight1<-weight[level==levels[i]]
   x2[i,]<-apply(temp,2,weighted.mean,weight1,na.rm=TRUE)}}
 colnames(x2)<-colnames(x)
 x2
}

levelx<-object$levelx
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
cat("\nLevel 2 total effect is:")
print(sum(c(de2,ie2,ie12)))
if(levelx==1)
{de1<-two(object$de1,object$level,object$weight)
 de1<-weighted.mean(de1,w2)
 cat("\nLevel 1 Direct Effect:")
 print(de1)
 if(!is.null(object$ie1))
 {ie1<-two(object$ie1,object$level,object$weight)
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
  {temp<-as.matrix(x[level==levels[i],])
   weight1<-weight[level==levels[i]]
   x2[i,]<-apply(temp,2,weighted.mean,weight1,na.rm=TRUE)}}
 colnames(x2)<-colnames(x)
 x2
}

op <- par(no.readonly = TRUE) # the whole list of settable par's.
if(is.null(var))
{par(mfrow=c(1,1))
 levelx<-object$levelx
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
 if(object$levelx==1)
 {readline("Press <return> to continue") 
  de1<-two(object$de1,object$level,object$weight)
  de1<-weighted.mean(de1,w2)
  if(!is.null(object$ie1))
  {ie1<-two(object$ie1,object$level,object$weight)
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
  {if(object$levelx==1)
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
  m.j<-two(m,object$level)
  y.j<-two(y,object$level)
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


summary.mlma.boot<-function(object,...,alpha=0.05)
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
return(list(total.effect1=total.effect1, total.effect2=total.effect2,
            direct.effect1=direct.effect1,direct.effect2=direct.effect2,
            indirect.effect1=indirect.effect1,indirect.effect2=indirect.effect2,
            indirect.effect12=indirect.effect12))
}


plot.mlma.boot<-function(x,..., var=NULL, alpha=0.05, cate=FALSE)
{object<-x
  boot.ci<-function(x,mat,alpha) #the first col of mat is the booted x, the rest cols
  #are boot results
{a1<-alpha/2
 a2<-1-a1
 x.boot<-mat[,1]
 vec<-mat[,-1]
 result<-matrix(NA,length(x),3)
 for (i in 1:length(x))
 {temp<-x.boot==x[i]
  result[i,1]<-mean(vec[temp])
  result[i,2:3]<-quantile(vec[temp],c(a1,a2))
 }
 result
}
op <- par(no.readonly = TRUE) # the whole list of settable par's.
a1<-(1:length(object$full$ie1_list[[1]]))[object$full$ie1_list[[1]]==var]
a2<-(1:length(object$full$ie2_list[[1]]))[object$full$ie2_list[[1]]==var]
x<-object$full$x
x.ord<-order(x)
x.j<-object$full$x.j
x.j.ord<-order(x.j)
b1<-alpha/2
b2<-1-b1
if(length(a1)!=0 & object$levelx==1)
{par(mfrow=c(2,length(object$full$ie1_list[[a1+1]])))
 for(i in object$full$ie1_list[[a1+1]])
 {temp.name<-colnames(object$ie1)[i]
  ie1<-cbind(object$xboot,object$ie1[, i])
  ie12<-matrix(object$ie12[, i],nrow=length(x.j))
  ie1<-boot.ci(x,ie1,alpha)
  plot(x[x.ord],ie1[x.ord,1],xlab="x",ylab="IE1",type="l",
       ylim=c(min(ie1[,2],na.rm=TRUE),max(ie1[,3],na.rm=TRUE)),
       main=paste("Level 1 IE of", temp.name,sep=" "))
  lines(x[x.ord],ie1[x.ord,2],col=2,lty=2)
  lines(x[x.ord],ie1[x.ord,3],col=3,lty=2)
  ie12_all<-cbind(apply(ie12,1,mean),apply(ie12,1,quantile,b1),apply(ie12,1,quantile,b2))
  plot(x.j[x.j.ord],ie12_all[x.j.ord,1],
       xlab="x.j",ylab="IE12",type="l",
       ylim=c(min(ie12_all[,2],na.rm=TRUE),max(ie12_all[,3],na.rm=TRUE)),
       main=paste("Level 2 IE of", temp.name,sep=" "))
  lines(x.j[x.j.ord],ie12_all[x.j.ord,2],col=2,lty=2)
  lines(x.j[x.j.ord],ie12_all[x.j.ord,3],col=3,lty=2)
 }}
else if(length(a1)!=0)
{par(mfrow=c(1,length(object$full$ie1_list[[a1+1]])))
 for(i in object$full$ie1_list[[a1+1]])
 {temp.name<-colnames(object$ie1)[i]
  ie12<-matrix(object$ie12[, i],nrow=length(x.j))
  ie12_all<-cbind(apply(ie12,1,mean),apply(ie12,1,quantile,b1),apply(ie12,1,quantile,b2))
  plot(x.j[x.j.ord],ie12_all[x.j.ord,1],
       xlab="x.j",ylab="IE12",type="l",
       ylim=c(min(ie12_all[,2],na.rm=TRUE),max(ie12_all[,3],na.rm=TRUE)),
       main=paste("Level 2 IE of", temp.name,sep=" "))
  lines(x.j[x.j.ord],ie12_all[x.j.ord,2],col=2,lty=2)
  lines(x.j[x.j.ord],ie12_all[x.j.ord,3],col=3,lty=2)
 }}
else
{par(mfrow=c(1,length(object$full$ie2_list[[a2+1]])))
 for(i in object$full$ie2_list[[a2+1]])
 {temp.name<-colnames(object$ie2)[i]
  ie2<-matrix(object$ie2[, i],nrow=length(x.j))
  ie2_all<-cbind(apply(ie2,1,mean),apply(ie2,1,quantile,b1),apply(ie2,1,quantile,b2))
  plot(x.j[x.j.ord],ie2_all[x.j.ord,1],
       xlab="x.j",ylab="IE2",type="l",
       ylim=c(min(ie2_all[,2],na.rm=TRUE),max(ie2_all[,3],na.rm=TRUE)),
       main=paste("Level 2 IE of", temp.name,sep=" "))
  lines(x.j[x.j.ord],ie2_all[x.j.ord,2],col=2,lty=2)
  lines(x.j[x.j.ord],ie2_all[x.j.ord,3],col=3,lty=2)
 }}
par(op)
}

