\name{boot.mlma}
\alias{boot.mlma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap Method for Inference on Multilevel Mediation Analysis
}
\description{
Bootstrap samples are selected from the original data set. The bootstrap sample has the same number of groups and in each group, the same number of observations as in the original data set. Based on each bootstrap sample, a multilevel mediation analysis is done and the results saved to make inferences on the total, direct and indirect effects.
}
\usage{
boot.mlma(y, data1=NULL,x=data1$parameter$x, m=data1$parameter$m, 
          levelx=data1$parameter$levelx, xref=NULL, yref=NULL, 
          l1=data1$parameter$l1,l2=data1$parameter$l2, 
          c1=data1$parameter$c1, #levelx is the level of x
          c1r=data1$parameter$c1r, c2=data1$parameter$c2, 
          c2r=data1$parameter$c2r,level=data1$parameter$level,  
          weight=rep(1,nrow(as.matrix(x))), 
          random="(1|level)", random.m1=NULL,intercept=TRUE, 
          family1=NULL, familym=vector("list",ncol(m)),
          covariates=NULL, cy1=NULL, cy2=NULL, cm=NULL,
          joint=NULL,f01y=data1$parameter$f01y,
          f10y=data1$parameter$f10y, f02ky=data1$parameter$f02ky, 
          f20ky=data1$parameter$f20ky, f01km1=data1$parameter$f01km1, 
          f01km2=data1$parameter$f01km2, f10km=data1$parameter$f10km,
          data2=NULL, x.new=NULL, m.new=m, level.new=level, 
          weight.new=NULL,covariates.new=covariates,boot=100,echo=TRUE)
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
  \item{data1}{
the results from data.org. If data1=NULL, needs to set up x, m, levelx, xref, l1, l2, c1,c2, c2r, f01y, f10y, f02ky, f20ky, f01km1, f01km2, f10km ...
}  
  \item{x}{
the vector of the predictive variable.
}
  \item{m}{
the mediators. The program will identify the levels and types of each mediator if not specified by l1, l2, c1, or c2. A mediator is identified as categorical if the mediator is a factor, a character, or has only two unque values.
}
  \item{levelx}{
the level of x (1 or 2). If it is not given, levelx will be decided by x.
}
  \item{xref}{
the reference group of x if it is binary. By default it will be the first level of x.
}
  \item{yref}{
the reference group of y if it is binary. By default it will be the first level of y.
}
  \item{l1}{
the column numbers of level 1 continuous mediators in m.
}
  \item{l2}{
the column numbers of level 2 continuous mediators in m.
}
  \item{c1}{
the column numbers of level 1 categorical mediators in m.
}
  \item{c1r}{
the reference groups of categorical mediators specified by c1.
}
  \item{c2}{
the column numbers of level 2 categorical mediators in m.
}
  \item{c2r}{
the reference groups of categorical mediators specified by c2.
}
  \item{f01y, f10y, f02ky, f20ky, f01km1, f01km2, f10km}{
the transformation functions as describe in the function \code{"\link[=data.org]{data.org}"}.  Need these arguments only when org.data=T.
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
the random effect part for model explaining the mediators. All other random effects are random = "(1|level)" if not specified here.
}
  \item{intercept}{
True if fit an intercept to models, by default.
}
  \item{family1}{
the glm family for the response variable y. If it is null, will be binomial with logic link for binary y and gaussian with identity link for continuous y.
}
  \item{familym}{
a list of length ncol(m), each item gives the glm family for the corresponding column of m. If an item is null, the family will be binomial with logic link for binary m and gaussian with identity link for continuous m.
}
%  \item{w2}{
%the weight for observations at level 2, which should be the same order as unique(level[!is.na(level)]). 
%}
  \item{boot}{
the number of bootstrapping samples.
}
  \item{covariates}{
the covariates matrix to explain the outcome, y, and/or the mediators, m.
}
  \item{cy1}{
the column numbers of covariates that are level 1 and used to explain y.
}
  \item{cy2}{
the column numbers of covariates that are level 2 and used to explain y.
}
  \item{cm}{
the column numbers of covariates that are used to explain m. cm[[1]] gives the mediators (in l1, cl, l2, or c2) that can be partially explained by covariates. Each of the rest items of the cm list shows the column number(s) in covariates that should be used to explain each mediator listed in cm[[1]] and by that order.For example, joint=list(1,c("m.2","m.4")) means find the joint effects of level 1 mediators m.2 and m.4.
}
  \item{joint}{
the list of group(s) of mediators whose joint mediation effect is of interests. joint[[1]] list the levels of mediators in each group and by the order of the list. Note that if any mediator in the group is of level 2, the level of the group should be 2.
}
  \item{data2}{
the results from data.org on x.new, covariates.new, m.new, level.new, and weight.new. If data1=NULL, needs to rerun data.org on the new data. If x.new is also NULL, set data2=data1. The new data is the same as the old data.
}  
  \item{x.new, covariates.new, m.new, level.new, weight.new}{
the settings that we want to make inferences on the mediation effects. If m.new=NULL, generate new mediators from x.new.
}
  \item{echo}{
print a i when finish the ith bootstrap estimation if echo=T.
}
}
\value{
Return a "mlma.boot" mode list, which include the following items:
  \item{de1}{direct effect(s) of level 1 exposure(s). de1 is a matrix of dimension n times boot by nx1, where n is the number of observations, and nx1 is the number of level 1 exposures. boot is the number of bootstrap samples.}
  \item{de2}{direct effect(s) of level 2 exposure(s). de2 is a matrix of dimension n2 times boot by nx2, where n2 is the number of unique levels, and nx2 is the number of level 2 exposures.}
  \item{ade1}{average direct effect(s) of level 1 exposure(s). ade1 is a matrix of dimension boot by nx1.}
  \item{ade2}{average direct effect(s) of level 2 exposure(s). ade2 is a matric of dimension boot by nx2.}
  \item{te1}{total effect of each level 1 exposure. te1 is a matrix of dimension n times boot by nx1, where n is the number of observations, and nx1 is the number of level 1 exposures.}
  \item{te2}{total effect of each level 2 exposure. te2 is a matrix of dimension n2 times boot by nx2, where n2 is the number of unique levels, and nx2 is the number of level 2 exposures.}
  \item{ate1}{average total effect(s) of level 1 exposure(s). ate1 is a matrix of dimension boot by nx1.}
  \item{ate2}{average total effect(s) of level 2 exposure(s). ate2 is a matrix of dimension boot by nx2.}
  \item{ie1}{level 1 indirect effect from level 1 exposure(s) to level 1 mediator(s) on the outcome. ie1 is an array of dimension (n*boot,nm1,nx1), where nm1 is the number of level 1 mediators.}
  \item{ie2}{level 2 indirect effect from level 2 exposure(s) to level 2 mediator(s) on the outcome. ie2 is an array of dimension (n2*boot,nm2,nx2), where nm2 is the number of level 2 mediators.}
  \item{ie12}{level 2 indirect effect from level 2 exposure(s) to level 1 mediator(s) on the outcome. ie12 is an array of dimension (n2*boot,nm1,nx2).}
  \item{aie1}{level 1 average indirect effect from level 1 exposure(s) to level 1 mediator(s) on the outcome. aie1 is a matrix of dimension (nm1*boot,nx1).}
  \item{aie2}{level 2 average indirect effect from level 2 exposure(s) to level 2 mediator(s) on the outcome. aie2 is a matrix of dimension (nm2*boot,nx2).}
  \item{aie12}{level 2 average indirect effect from level 2 exposure(s) to level 1 mediator(s) on the outcome. aie12 is a matrix of dimension (nm1*boot,nx2).}
  \item{je1}{joint level 1 indirect effect from level 1 exposure(s) to joint level 1 mediators on the outcome. je1 is an array of dimension (n*boot,nj1,nx1), where nj1 is the number of groups of level 1 mediators.}
  \item{je2}{joint level 2 indirect effect from level 2 exposure(s) to joint level 2 mediators on the outcome. je2 is an array of dimension (n2*boot,nj2,nx2), where nj2 is the number of groups oflevel 2 mediators.}
  \item{je12}{joint level 2 indirect effect from level 2 exposure(s) to joint level 1 mediators on the outcome. je12 is an array of dimension (n2*boot,nj1,nx2).}
  \item{aje1}{average joint level 1 indirect effect from level 1 exposure(s) to joint level 1 mediators on the outcome. aje1 is a matrix of dimension (nj1*boot,nx1).}
  \item{aje2}{average joint level 2 indirect effect from level 2 exposure(s) to joint level 2 mediators on the outcome. je2 is a matrix of dimension (nj2*boot,nx2), where nj2 is the number of groups oflevel 2 mediators.}
  \item{aje12}{average joint level 2 indirect effect from level 2 exposure(s) to joint level 1 mediators on the outcome. je12 is a matrix of dimension (nj1*boot,nx2).}
  \item{full}{an "mlma" results using the original data set.}
  \item{xboot}{a n*boot vector of the level 1 predictors in all boot bootstrap samples.}
  \item{xjboot}{a g*boot vector of the (aggregated) level 2 predictors in all boot bootstrap samples.}
  \item{levelx}{inherited from the same argument.}
  \item{level}{inherited from the same argument.}
}
\author{
Qingzhao Yu (qyu@lsuhsc.edu), Bin Li (bli@lsu.edu).
}
%\Warning{
%Bootstrap sample may result in "Error in fn(x, ...) : Downdated VtV is not positive %definite", set seed can help avoid such mistakes.
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(sim.111)
sim.111$m[,3]=as.factor(sim.111$m[,3])
data2<-data.org(ifelse(sim.111$x>1.9,1,0), m=sim.111$m, 
                f10y=list(1,c("x^2","sqrt(x+6)")), 
                f20ky=list(2,c("x","x^3")), 
                f10km=list(matrix(c(2,1),1),"log(x+2)"), level=sim.111$level)
temp2.boot<-boot.mlma(y=sim.111$y, data1=data2,boot=2)

\donttest{
#can also do the above analysis using the following code
temp2.boot<-boot.mlma(y=ifelse(sim.111$y>4.5,1,0), x=sim.111$x, m=sim.111$m, 
            f10y=list(1,c("x^2","sqrt(x+6)")), 
            f20ky=list(2,c("x","x^3")), 
            f10km=list(matrix(c(2,1),1),"log(x+2)"), level=sim.111$level,boot=2)

#with a level 2 exposure
data(sim.211)
data1<-data.org(x=ifelse(sim.211$x>2.2,1,0), m=sim.211$m, 
                f02ky=list(1,c("x","x^2")),  
                f20ky=list(2,c("x","x^3")), f01km2=list(matrix(c(1,1),1),c("x^1.2","x^2.3")),
                f01km1=list(matrix(c(2,1),1),"sqrt(x)+3"),level=sim.211$level) 
temp1.boot<-boot.mlma(y=sim.211$y, data1=data1,boot=2) 

#with both level 1 and 2 exposure
data3<-data.org(x=cbind(sim.211$x,sim.111$x), m=sim.211$m, 
                f20ky=list(2,c("x","x^3")), f01km1=list(matrix(c(2,1),1),"sqrt(x)+3"), 
                f01km2=list(matrix(c(1,1),1),c("x^1.2","x^2.3")), level=sim.211$level)
temp3.boot<-boot.mlma(y=sim.211$y, data1=data3,boot=2) }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Multilevel Mediation Analysis}
\keyword{ ~Inferences on MLMA}% __ONLY ONE__ keyword per line
