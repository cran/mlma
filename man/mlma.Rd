\name{mlma}
\alias{mlma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multilevel Mediation Analysis
}
\description{
The function transforms the data set and does multilevel mediation analysis.  The total, direct, and indirect effects will be returned as the results.
}
\usage{
mlma(y, data1=NULL, x=data1$parameter$x, m=data1$parameter$m, 
               yref=NULL, xref=NULL, levelx=data1$parameter$levelx, 
               l1=data1$parameter$l1,l2=data1$parameter$l2, 
               c1=data1$parameter$c1, #levelx is the level of x
               c1r=data1$parameter$c1r, c2=data1$parameter$c2, 
               c2r=data1$parameter$c2r,level=data1$parameter$level,  
               weight=rep(1, nrow(data.frame(x))), 
               random="(1|level)", random.m1=NULL,intercept=TRUE, 
               family1=NULL, familym=vector("list",ncol(m)),
               covariates=NULL, cy1=NULL, cy2=NULL, cm=NULL,joint=NULL,
               f01y=data1$parameter$f01y,f10y=data1$parameter$f10y, 
               f02ky=data1$parameter$f02ky, f20ky=data1$parameter$f20ky, 
               f01km1=data1$parameter$f01km1,f01km2=data1$parameter$f01km2, 
               f10km=data1$parameter$f10km, data2=NULL, x.new=NULL, 
               m.new=NULL, level.new=level, weight.new=NULL,
               covariates.new=covariates)
}
%- maybe also 'usage' for other objects documented here.
\details{
The multilevel mediation is based on the following linear multilevel additive models:
\deqn{Y_{ij} = u_{0j}^Y(X_{.j}, \mathbf{M}_{.j}, \mathbf{Z}_{.j})+{\boldsymbol{\beta}_{10}^Y}^T\mathbf{f}_{10}^Y(X_{ij}-X_{.j})+\sum_{k=1}^K{\boldsymbol{\beta}_{20k}^Y}^T\mathbf{f}_{20k}^Y(M_{ijk}-M_{.jk})+{\boldsymbol{\beta}_{30}^Y}^T(\mathbf{Z}_{ij}-\mathbf{Z}_{.j})+r_{ij}^Y,}
where \deqn{u_{0j}^Y(X_{.j}, \mathbf{M}_{.j}, \mathbf{Z}_{.j}) = c_{00}^Y + {\boldsymbol{\beta}_{01}^Y}^T\mathbf{f}_{01}^Y(X_{.j}) + \sum_{k=1}^K{\boldsymbol{\beta}_{02k}^Y}^T\mathbf{f}_{02k}^Y(M_{.jk}) + {\boldsymbol{\beta}_{03}^Y}^T\mathbf{Z}_{.j} + r_{0j}^Y.}
For \eqn{k=1,\ldots,K,}
\deqn{M_{.jk} = u_{0jk}^M(X_{.j})+{\boldsymbol{\beta}_{10k}^M}^T\mathbf{f}_{10k}^M(X_{ij}-X_{.j})+r_{ijk}^M,}
\deqn{u_{0jk}^M(X_{.j}) = c_{00k}^M + {\boldsymbol{\beta}_{01k}^M}^T\mathbf{f}_{01k}^{M1}(X_{.j}) + r_{0jk}^M.}
If for some k, \eqn{M_k} is level 2 variable,
\deqn{M_{.jk} = c_{00k}^M + {\boldsymbol{\beta}_{01k}^M}^T\mathbf{f}_{01k}^{M2}(X_{.j}) + r_{0jk}^M.} 

Note that in the models, \eqn{\mathbf{f}(\cdot)=(f_1(\cdot), f_2(\cdot), \cdots, f_l(\cdot))^T} is a set of l transformation functions on \eqn{\cdot}, with the corresponding linear coefficients vector \eqn{\boldsymbol{\beta}=(\beta_1, \beta_2, \cdots, \beta_l)^T}. \eqn{\mathbf{f}} and l are known for model fitting. l may be different with \eqn{\mathbf{f}} of different sub- and super-scripts.}
\arguments{
  \item{y}{
the vector of the outcome variable.
}
  \item{data1}{The transformed and organized data set from data.org.  If the data set has not been organized, leave data1=NULL (by default), and set the transformation functions (f arguments). Otherwise, set data1 as the output from the data.org function and do not include the arguments starting with fs.}
  \item{x}{
the vector of the predictive variable.
}
  \item{m}{
the mediators. The program will identify the levels and types of each mediator if not specified by l1, l2, c1, or c2. A mediator is identified as categorical if the mediator is a factor, a character, or has only two unque values.
}
  \item{yref}{
the reference group of y if it is binary. By default it will be the first level of y.
}
  \item{xref}{
the reference group of x if it is binary. By default it will be the first level of x.
}
  \item{levelx}{
the level of x (1 or 2). If it is not given, levelx will be decided by x.
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
  \item{level}{
a vector that record the group number for each observation.
}
  \item{weight}{
the weight of cases in groups.
}
  \item{random}{
the random effect part for the full model. random = "(1|level)" by default.
}
  \item{random.m1}{
the random effect part for model explaining the mediators. If not null, 1st item of random.m1 is the list of level 1 mediators, following items are the random item of the same order. All other random effects are random = "(1|level)" if not specified here.
}
  \item{intercept}{
True if fit an intercept to models, by default.
}
  \item{covariates}{
the covariates matrix to explain the outcome, y, and/or the mediators, m.
}
  \item{family1}{
the glm family for the response variable y. If it is null, will be binomial with logic link for binary y and gaussian with identity link for continuous y.
}
  \item{familym}{
a list of length ncol(m), each item gives the glm family for the corresponding column of m. If an item is null, the family will be binomial with logic link for binary m and gaussian with identity link for continuous m.
}
  \item{cy1}{
the column numbers of covariates that are level 1 and used to explain y.
}
  \item{cy2}{
the column numbers of covariates that are level 2 and used to explain y.
}
  \item{cm}{
the column numbers of covariates that are used to explain m. cm[[1]] gives the mediators (in l1, cl, l2, or c2) that can be partially explained by covariates. Each of the rest items of the cm list shows the column number(s) in covariates that should be used to explain each mediator listed in cm[[1]] and by that order.
}
  \item{joint}{
the list of group(s) of mediators whose joint mediation effect is of interests. joint[[1]] list the levels of mediators in each group and by the order of the list. Note that if any mediator in the group is of level 2, the level of the group should be 2.
}
%  \item{org.data}{
%if is TRUE, first organize the data set and do transformations using the function \code{"\link[=data.org]{data.org}"}. In such case, need to specify the transformation function arguments.
%}
  \item{f01y, f10y, f02ky, f20ky, f01km1, f01km2, f10km}{
the transformation functions as describe in the function \code{"\link[=data.org]{data.org}"}.  Need these arguments only when org.data=T.
  }
  \item{x.new, m.new, level.new, weight.new, covariates.new}{
the settings that we want to make inferences on the mediation effects. %If m.new=NULL, generate new mediators from x.new.
}
  \item{data2}{The transformed and organized data set from data.org on the set of new x.new and m.new etc..  If the data set has not been organized, leave data2=NULL (by default).}
}
\value{
A "mlma" mode list will be returned with the following items:
  \item{de1}{direct effect(s) of level 1 exposure(s). de1 is a matrix of dimension n by nx1, where n is the number of observations, and nx1 is the number of level 1 exposures.}
  \item{de2}{direct effect(s) of level 2 exposure(s). de2 is a matrix of dimension n2 by nx2, where n2 is the number of unique levels, and nx2 is the number of level 2 exposures.}
  \item{ade1}{average direct effect(s) of level 1 exposure(s). ade1 is a vector of length nx1.}
  \item{ade2}{average direct effect(s) of level 2 exposure(s). ade2 is a vector of length nx2.}
  \item{te1}{total effect of each level 1 exposure. te1 is a matrix of dimension n by nx1, where n is the number of observations, and nx1 is the number of level 1 exposures.}
  \item{te2}{total effect of each level 2 exposure. te2 is a matrix of dimension n2 by nx2, where n2 is the number of unique levels, and nx2 is the number of level 2 exposures.}
  \item{ate1}{average total effect(s) of level 1 exposure(s). ate1 is a vector of length nx1.}
  \item{ate2}{average total effect(s) of level 2 exposure(s). ate2 is a vector of length nx2.}
  \item{ie1}{level 1 indirect effect from level 1 exposure(s) to level 1 mediator(s) on the outcome. ie1 is an array of dimension (n,nm1,nx1), where nm1 is the number of level 1 mediators.}
  \item{ie2}{level 2 indirect effect from level 2 exposure(s) to level 2 mediator(s) on the outcome. ie2 is an array of dimension (n2,nm2,nx2), where nm2 is the number of level 2 mediators.}
  \item{ie12}{level 2 indirect effect from level 2 exposure(s) to level 1 mediator(s) on the outcome. ie12 is an array of dimension (n2,nm1,nx2).}
  \item{aie1}{level 1 average indirect effect from level 1 exposure(s) to level 1 mediator(s) on the outcome. aie1 is a matrix of dimension (nm1,nx1).}
  \item{aie2}{level 2 average indirect effect from level 2 exposure(s) to level 2 mediator(s) on the outcome. aie2 is a matrix of dimension (nm2,nx2).}
  \item{aie12}{level 2 average indirect effect from level 2 exposure(s) to level 1 mediator(s) on the outcome. aie12 is a matrix of dimension (nm1,nx2).}
  \item{je1}{joint level 1 indirect effect from level 1 exposure(s) to joint level 1 mediators on the outcome. je1 is an array of dimension (n,nj1,nx1), where nj1 is the number of groups of level 1 mediators.}
  \item{je2}{joint level 2 indirect effect from level 2 exposure(s) to joint level 2 mediators on the outcome. je2 is an array of dimension (n2,nj2,nx2), where nj2 is the number of groups oflevel 2 mediators.}
  \item{je12}{joint level 2 indirect effect from level 2 exposure(s) to joint level 1 mediators on the outcome. je12 is an array of dimension (n2,nj1,nx2).}
  \item{aje1}{average joint level 1 indirect effect from level 1 exposure(s) to joint level 1 mediators on the outcome. aje1 is a matrix of dimension (nj1,nx1).}
  \item{aje2}{average joint level 2 indirect effect from level 2 exposure(s) to joint level 2 mediators on the outcome. je2 is a matrix of dimension (nj2,nx2), where nj2 is the number of groups oflevel 2 mediators.}
  \item{aje12}{average joint level 2 indirect effect from level 2 exposure(s) to joint level 1 mediators on the outcome. je12 is a matrix of dimension (nj1,nx2).}
  \item{f1}{the overall multilevel model.}
  \item{fm1}{a list, where the first item identifies the level 1 mediators, and in that order, the following items give the prediction functions of the mediators.}
  \item{fm2}{a list, where the first item identifies the level 2 mediators, and in that order, the following items give the prediction functions of the mediators.}
  %\item{ie1_list}{a list, where the first item identifies the level 1 mediators, and in that order, the following items give the column(s) of the indirect effects of the mediator in ie1.}
  %\item{ie2_list}{a list, where the first item identifies the level 2 mediators, and in that order, the following items give the column(s) of the indirect effects of the mediator in ie2.}
  %\item{iej2_list}{a list, where the first item identifies the level 2 joint mediators, and in that order, the following items give the column(s) of the indirect effects of the mediator in cbind(ie12,ie2).}
  \item{ie1_1, ie1_2, ie1_3, ie2_1, ie2_2, ie2_3}{the first, second and third part of the corresponding indirect effects.}
  \item{x}{the matrix of the new exposure variable(s) (x.new).}
  \item{x.j}{the vector of the aggregated x at the higher level by the order of unique(level.new[!is.na(level.new)]).}
  %\item{levelx}{the level of the predictive variable.}
  %\item{weight}{the weight for each observation.}
  \item{data1}{The results from data.org on the original data (x, m, etc.).}
  \item{data2}{The results from data.org on the new data (x.new, m.new, etc.).}
}
\author{
Qingzhao Yu (qyu@lsuhsc.edu), Bin Li (bli@lsu.edu).
}
\examples{
#with a level 1 exposure
data(sim.111)
data2<-data.org(sim.111$x, m=sim.111$m, 
                   f10y=list(1,c("x^2","sqrt(x+6)")), 
                   f20ky=list(2,c("x","x^3")), 
                   f10km=list(matrix(c(2,1),1),"log(x+2)"), level=sim.111$level)
temp2<-mlma(y=sim.111$y, data1=data2)

#can also do the above analysis using the following code
temp2<-mlma(y=sim.111$y, x=sim.111$x, m=sim.111$m, 
           f10y=list(1,c("x^2","sqrt(x+6)")), 
           f20ky=list(2,c("x","x^3")), 
           f10km=list(matrix(c(2,1),1),"log(x+2)"), level=sim.111$level)
                
#with a level 2 exposure
data(sim.211)
data1<-data.org(x=sim.211$x, m=sim.211$m, 
                   f01y=list(1,c("x","log(x^2)")), f02ky=list(1,c("x","x^2")),  
                   f20ky=list(2,c("x","x^3")), f01km2=list(matrix(c(1,1),1),c("x^1.2","x^2.3")),
                   f01km1=list(matrix(c(2,1),1),"sqrt(x)+3"),level=sim.211$level) 
temp1<-mlma(y=sim.211$y, data1) 

#with both level 1 and 2 exposure
data3<-data.org(x=cbind(sim.211$x,sim.111$x), m=sim.211$m, 
                   f01y=list(1,c("x","log(x^2)")), f02ky=list(1,c("x","x^2")), 
                   f20ky=list(2,c("x","x^3")), f01km1=list(matrix(c(2,1),1),"sqrt(x)+3"), 
                   f01km2=list(matrix(c(1,1),1),c("x^1.2","x^2.3")), level=sim.211$level)
temp3<-mlma(y=sim.211$y, data3) 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Multilevel Mediation Analysis}
