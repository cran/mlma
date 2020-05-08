\name{data.org}
\alias{data.org}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform and Organize Data for Mediation Analysis
}
\description{
To transform variables and generate data sets for mediation analysis.
}
\usage{
data.org(x, m, levely=1,  y=NULL, levelx=NULL, xref=NULL, 
         l1 = NULL, l2 = NULL, c1 = NULL, 
         c1r = rep(1, length(c1)), c2 = NULL, 
         c2r = rep(1, length(c2)), f01y = NULL, f10y = NULL, 
         f02ky = NULL, f20ky = NULL, f01km1 = NULL, 
         f01km2 = NULL, f10km = NULL,  
         level = 1:nrow(as.matrix(x)), weight = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
the data frame of exposure variable(s). If there are multiple exposures, a categorical exposure should be binarized with the reference group coded as 0. 
}
  \item{m}{
the data frame of all mediators. The program will identify the levels and types of each mediator if not specified by l1, l2, c1, or c2. A mediator is identified as categorical if the mediator is a factor, a character, or has only two unque values.
}
  \item{levely}{
the level of y (1 or 2), 1 by default. If it is NULL, levely will be decided automatically by y.
}
  \item{y}{
only needed if levely is not given.
}
  \item{levelx}{
the vector of levels of every column of x (1 or 2). If it is not given, levelx will be decided automaticaly. If there is only level 1 exposure, but there are level 2 mediator(s), the level 2 exposure will be created by aggregate all the level 1 exposure(s) to level 2.
}
  \item{xref}{
the reference group of x if it is a single binary variable. By default it will be the first level of x.
  }
  \item{l1}{
the column numbers of level 1 continuous mediators in m or the list of names of the level 1 continuous mediators.
}
  \item{l2}{
the column numbers of level 2 continuous mediators in m or the list of names of the level 2 continuous mediators.
}
  \item{c1}{
the column numbers of level 1 categorical mediators in m or the list of names of the level 1 categorical mediators.
}
  \item{c1r}{
the reference groups of categorical mediators specified by c1.
}
  \item{c2}{
the column numbers of level 2 categorical mediators in m or the list of names of the level 2 categorical mediators.
}
  \item{c2r}{
the reference groups of categorical mediators specified by c2.
}
  \item{f01y}{
the transformation function expressions on level 2 exposures (x.jm) in explaining y (eg, list(1,c("x^2","log(x)"))). The first item lists column numbers/variable names of the level 2 exposures in x, which needs to be transformed. By that order, each of the rest items of f01y list the transformation functional expressions for each exposure. The exposures not specified in the list will not be transformed in any way.
}
  \item{f10y}{
the transformation function expressions on level 1 exposure (xijm) in explaining y. (eg, list(2,c("x^2","log(x)"))). The first item lists column numbers/variable names of the level 1 exposures in x, which needs to be transformed. By that order, each of the rest items of f10y list the transformation functional expressions for each exposure. The exposures not specified in the list will not be transformed in any way.
}
  \item{f02ky}{
the transformation-function-expression list on level 2 mediators (m.jk) in explaining y (eg, list(2:3,c("log(x)","sqrt(x)"), "2*x")). The first item lists column numbers/variable names of the level 2 mediators in m, which needs to be transformed. By that order, each of the rest items of f01ky list the transformation functional expressions for each mediator. The mediators not specified in the list will not be transformed in any way.
}
  \item{f20ky}{
the transformation-function-expression list on level 1 mediators (mijk-m.jk) in explaining y.The first item lists column numbers/variable names of the level 1 mediators in m, which needs to be transformed. By that order, each of the rest items of f02ky list the transformation functional expressions for each mediator. The mediators not specified in the list will not be transformed in any way.
}
  \item{f01km1}{
the transformation-function-expression list on level 2 predictor (x.jm) in explaining the level 1 mediators. The first item is a matrix with two columns, the first column is the column numbers of the level 1 mediators in m, which should be explained by the transformed predictor(s). The second column indicates the column number of the exposure in x that will be transformed to explain the mediator identified by the 1st column of the same row. By the order of the rows of the first item, each of the rest items of f01km1 lists the transformation functional expressions for the level 2 exposure (identified by column 2) in explaining each mediator (identified by column 1). The mediators not specified in the list will be explained by the original format of the level 2 exposures. e.g.list(matrix(c(1,2,1,1),2,2), "x^2","x^2") means x[,1]^2 is used to explain both m[,1] and m[,2].
}
  \item{f01km2}{
the transformation-function-expression list on level 2 predictor (x.jm) in explaining the level 2 mediators. The first item is a matrix with two columns, the first column is the column numbers of the level 2 mediators in m, which should be explained by the transformed predictor(s). The second column indicates the column number of the exposure in x that will be transformed to explain the mediator identified by the 1st column of the same row. By the order of the rows of the first item, each of the rest items of f01km2 lists the transformation functional expressions for the level 2 exposure (identified by column 2) in explaining each mediator (identified by column 1). The mediators not specified in the list will be explained by the original format of all level 2 exposures.
}
  \item{f10km}{
the transformation-function-expression list on level 1 predictor (xijm) in explaining the level 1 mediators. The first item is a matrix with two columns, the first column is the column numbers of the level 1 mediators in m, which should be explained by the transformed level 1 exposure(s). The second column indicates the column number of the exposure in x that will be transformed to explain the mediator identified by the 1st column of the same row. By the order of the rows of the first item, each of the rest items of f10km lists the transformation functional expressions for the level 1 exposure (identified by column 2) in explaining each mediator (identified by column 1). The mediators not specified in the list will be explained by the original format of all level 1 exposures.
}
  \item{level}{
a vector that record the group number for each observation.
}
  \item{weight}{
the weight of cases in groups.
}
}

\details{
The arguments starting with "f" are used to specify the transformation functions of the predictor or mediators in explaining y, or the transformation functions of the predictor in explaining the mediators. If the name of the argument includes a "k", the transformation is on the mediators.  If the names of the arguments end with "y", the transformation is to explain the outcome. Otherwise, the transformation is on x to predict mediators (the argument ends with "m1" or "m" (for level 1 mediator), or "m2" (for level 2 mediator)). The functions corresponds to the functions in the following multilevel additive models, reading as f+subscript+superscript.  For example, f01y specifies \eqn{{f}_{01}^Y}.
\deqn{Y_{ij} = u_{0j}^Y(X_{.j}, \mathbf{M}_{.j}, \mathbf{Z}_{.j})+{\boldsymbol{\beta}_{10}^Y}^T\mathbf{f}_{10}^Y(X_{ij}-X_{.j})+\sum_{k=1}^K{\boldsymbol{\beta}_{20k}^Y}^T\mathbf{f}_{20k}^Y(M_{ijk}-M_{.jk})+{\boldsymbol{\beta}_{30}^Y}^T(\mathbf{Z}_{ij}-\mathbf{Z}_{.j})+r_{ij}^Y,}

where \deqn{u_{0j}^Y(X_{.j}, \mathbf{M}_{.j}, \mathbf{Z}_{.j}) = c_{00}^Y + {\boldsymbol{\beta}_{01}^Y}^T\mathbf{f}_{01}^Y(X_{.j}) + \sum_{k=1}^K{\boldsymbol{\beta}_{02k}^Y}^T\mathbf{f}_{02k}^Y(M_{.jk}) + {\boldsymbol{\beta}_{03}^Y}^T\mathbf{Z}_{.j} + r_{0j}^Y.}
For \eqn{k=1,\ldots,K,}
\deqn{M_{.jk} = u_{0jk}^M(X_{.j})+{\boldsymbol{\beta}_{10k}^M}^T\mathbf{f}_{10k}^M(X_{ij}-X_{.j})+r_{ijk}^M,}
\deqn{u_{0jk}^M(X_{.j}) = c_{00k}^M + {\boldsymbol{\beta}_{01k}^M}^T\mathbf{f}_{01k}^{M1}(X_{.j}) + r_{0jk}^M.}
If for some k, \eqn{M_k} is level 2 variable,
\deqn{M_{.jk} = c_{00k}^M + {\boldsymbol{\beta}_{01k}^M}^T\mathbf{f}_{01k}^{M2}(X_{.j}) + r_{0jk}^M.} 

The transformation function can be any function that is differentiable by the function deriv(), or the ifelse function with those functions. The transformation function can also be the ns() and bs() functions for natural and b spline basis.
}
\value{
The function returns a list with transformed and organized data with the following items:
  \item{x1}{the level 1 and 2 transformed predictor variable matrix in explaining y (eg, f01y(x.j) & f10y(xij)).}
  \item{lx}{a matrix with ncol(x) rows and three columns. The ith row is for the jth column of x. The first column of the ith row gives the level (1 or 2) of the ith column of x. lx[i,2]:lx[i,3] gives the column numbers of transformed ith x in x1.}
  \item{m1y}{the level 1 mediator matrix in explaining y (eg, f20ky(mijk-m.jk) & mijk or binarized mijk for categorical mediators).}
  \item{m1}{a list where the first item identify column numbers of level 1 mediators in m (ie, c(l1,c1)). For every mediator identified by m1[[1]] and by that order, each of the rest item identify the column number(s) in m1y the (transformed) value(s) of the mediator in explaining y.}
  \item{m2y}{the level 2 mediator (original or aggregated) matrix in explaining y (eg, f02ky(m.jk) & m.jk).}
  \item{m2}{a list where the first item identify column numbers of level 2 mediators in m (i.e., c(l2,c2)). For every mediator identified by m2[[1]] and by that order, each of the rest item identify the column number(s) in m2y the (transformed) value(s) of the mediator in explaining y.}
%  \item{m12}{a list where the first item identify column numbers of aggregated level 2 mediators in m (eg, aggregated f20ky). For every mediator identified by m12[[1]] and by that order, each of the rest item identify the column number(s) in m2y the aggregated value(s) of the mediator in explaining y.}
  \item{xm1}{the (transformed) level 1 and level 2 predictor(s) in explaining level 1 mediators.}
  \item{fm11}{a list where the first item identify column numbers of level 1 mediators in m. For every mediator identified by fm11[[1]] and by that order, each of the rest item identify the column number(s) in xm1 the (transformed) level 1 predictor(s) in explaining the mediator.}
  \item{fm12}{a list where the first item identify column numbers of level 1 mediators in m. For every mediator identified by fm12[[1]] and by that order, each of the rest item identify the column number(s) in xm1 the (transformed/aggregated) level 2 predictor(s) in explaining the mediator.}
  \item{m.2}{a matrix of level 2 mediators (one row for each group).}
  \item{xm2}{the (transformed/aggregated) level 2 predictor(s) in explaining level 2 mediators (one row for each group).}
  \item{fm22}{a list where the first item identify column numbers of level 2 mediators in m. For every mediator identified by fm22[[1]] and by that order, each of the rest item identify the column number(s) in xm2 the (transformed) level 2 predictor(s) in explaining the mediator.}
  \item{x1.der, m2y.der, m1y.der, xm2.der, xm1.der}{the derivative of x1, m2y, m1y, xm2, and xm1 respectively.}
  \item{f01km1.2}{the first item is f01km1[[1]]. The rest are the column numbers of f01km1[[i]] in xm1.}
  \item{f01km2.2}{the first item is f01km2[[1]]. The rest are the column numbers of f01km2[[i]] in xm2.}
  \item{f10km.2}{the first item is f10km[[1]]. The rest are the column numbers of f10km[[i]] in xm1.}
  \item{binx}{the ith item of binx is TRUE if the ith predictor is binary, otherwise is FALSE.}
  \item{parameter}{The list of all arguments.}
}
\author{
Qingzhao Yu (qyu@lsuhsc.edu), Bin Li (bli@lsu.edu).
}
\examples{
data(sim.211)
example1<-data.org(x=sim.211$x, m=sim.211$m, 
                   f01y=list(1,c("x","log(x^2)")), 
                   f02ky=list(1,c("x","x^2")),  
                   f20ky=list(2,c("x","x^3")), 
                   f01km2=list(matrix(c(1,1),1),c("x^1.2","x^2.3")),
                   f01km1=list(matrix(c(2,1),1),"sqrt(x)+3"),
                   level=sim.211$level) 

data(sim.111)
example2<-data.org(sim.111$x, m=sim.111$m, 
                   f10y=list(1,c("x^2","sqrt(x+6)")), 
                   f20ky=list(2,c("x","x^3")), 
                   f10km=list(matrix(c(2,1),1),"log(x+2)"), 
                   level=sim.111$level) #f01y=list(1,c("x","log(x^2)")),

example3<-data.org(x=cbind(sim.211$x,sim.111$x), m=sim.211$m, 
                   f01y=list(1,c("x","log(x^2)")), f02ky=list(1,c("x","x^2")), 
                   f20ky=list(2,c("x","x^3")), f01km1=list(matrix(c(2,1),1),"sqrt(x)+3"), 
                   f01km2=list(matrix(c(1,1),1),c("x^1.2","x^2.3")), level=sim.211$level)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Data Transformation }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
