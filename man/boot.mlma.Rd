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
boot.mlma(y, biny = FALSE, x, levelx, levely, m, l1 = NULL, l2 = NULL, c1 = NULL, 
          c1r = rep(1, length(c1)), c2 = NULL, c2r = rep(1, length(c2)), 
          f01y = NULL, f10y = NULL, f02ky = NULL, f20ky = NULL, f01km1 = NULL, 
          f01km2 = NULL, f10km = NULL, level, weight = rep(1, length(x)), 
          random = "(1|level)", random.m1 = NULL, intercept = TRUE,
          w2=rep(1,length(unique(level[!is.na(level)]))),
          boot = 100, seed = 1, covariates=NULL,cy1=NULL,cy2=NULL,cm=NULL,
          joint=NULL,x.new=x, m.new=m, level.new=level,weight.new=weight,
          covariates.new=covariates)
}
%- maybe also 'usage' for other objects documented here.
\details{
The multilevel mediation is based on the following linear multilevel additive models:
\deqn{Y_{ij} = u_{0j}^Y(X_{.j}, \mathbf{M}_{.j}, \mathbf{Z}_{.j})+u_{1j}^Y(X_{ij}-X_{.j})+\sum_{k=1}^Ku_{2jk}^Y(M_{ijk}-M_{.jk})+u_{3j}^Y(\mathbf{Z}_{ij}-\mathbf{Z}_{.j})+r_{ij}^Y}

where \deqn{u_{0j}^Y(X_{.j}, \mathbf{M}_{.j}, \mathbf{Z}_{.j}) = c_{00}^Y + {\boldsymbol{\beta}_{01}^Y}^T\mathbf{f}_{01}^Y(X_{.j}) + \sum_{k=1}^K{\boldsymbol{\beta}_{02k}^Y}^T\mathbf{f}_{02k}^Y(M_{.jk}) + {\boldsymbol{\beta}_{03}^Y}^T\mathbf{Z}_{.j} + r_{0j}^Y}
\deqn{u_{1j}^Y(X_{ij}-X_{.j}) = {\boldsymbol{\beta}_{10}^Y}^T\mathbf{f}_{10}^Y(X_{ij}-X_{.j}) + r_{1j}^Y}
\deqn{u_{2jk}^Y(M_{ijk}-M_{.jk}) = {\boldsymbol{\beta}_{20k}^Y}^T\mathbf{f}_{20k}^Y(M_{ijk}-M_{.jk}) + r_{2jk}^Y, \mbox{~~~~~~}k=1,\ldots,K}
\deqn{u_{3j}^Y(\mathbf{Z}_{ij}-\mathbf{Z}_{.j}) = {\boldsymbol{\beta}_{30}^Y}^T(\mathbf{Z}_{ij}-\mathbf{Z}_{.j}) + r_{3j}^Y}
For \eqn{k=1,\ldots,K,}
\deqn{M_{.jk} = c_{00k}^M+u_{1jk}^M(X_{ij}-X_{.j})+r_{ijk}^M}
\deqn{u_{0jk}^M(X_{.j}) = c_{00k}^M + {\boldsymbol{\beta}_{01k}^M}^T\mathbf{f}_{01k}^{M1}(X_{.j}) + r_{0jk}^M}
\deqn{u_{1jk}^M(X_{ij}-X_{.j}) = {\boldsymbol{\beta}_{10k}^M}^T\mathbf{f}_{10k}^M(X_{ij}-X_{.j}) + r_{1jk}^M}
If for some k, \eqn{M_k} is level 2 variable,
\deqn{M_{.jk} = c_{00k}^M + {\boldsymbol{\beta}_{01k}^M}^T\mathbf{f}_{01k}^{M2}(X_{.j}) + r_{0jk}^M} 

Note that in the models, \eqn{\mathbf{f}(\cdot)=(f_1(\cdot), f_2(\cdot), \cdots, f_l(\cdot))^T} is a set of l transformation functions on \eqn{\cdot}, with the corresponding linear coefficients vector \eqn{\boldsymbol{\beta}=(\beta_1, \beta_2, \cdots, \beta_l)^T}. \eqn{\mathbf{f}} and l are known for model fitting. l may be different with \eqn{\mathbf{f}} of different sub- and super-scripts.}
\arguments{
  \item{y}{
the vector of the outcome variable.
}
  \item{biny}{
True if the outcome is binary, otherwise False. Default is False.
}
  \item{x}{
the vector of the predictive variable.
}
  \item{levelx}{
the level of x (1 or 2), 1 by default.
}
  \item{levely}{
the level of y (1 or 2), 1 by default.
}
  \item{m}{
the matrix or vector of mediators.
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
  \item{f01y}{
the transformation function expressions on level 2 predictive variable (x.j) in explaining y (eg, c("x^2","log(x)")).
}
  \item{f10y}{
the transformation function expressions on level 1 predictive variable (xij-x.j) in explaining y.
}
  \item{f02ky}{
the transformation-function-expression list on level 2 mediators (m.jk) in explaining y (eg, list(2:3,c("log(x)","sqrt(x)"), "2*x")). The first item lists column numbers of the level 2 mediators in m, which needs to be transformed. By that order, each of the rest items of f01ky list the transformation functional expressions for each mediator. The mediators not specified in the list will not be transformed in any way.
}
  \item{f20ky}{
the transformation-function-expression list on level 1 mediators (mijk-m.jk) in explaining y.The first item lists column numbers of the level 1 mediators in m, which needs to be transformed. By that order, each of the rest items of f02ky list the transformation functional expressions for each mediator. The mediators not specified in the list will not be transformed in any way.
}
  \item{f01km1}{
the transformation-function-expression list on level 2 predictor (x.j) in explaining the level 1 mediators. The first item lists column numbers of the level 1 mediators in m, which should be explained by the transformed predictor(s). By that order, each of the rest items of f01km1 lists the transformation functional expressions for the (aggregated) level 2 predictor in explaining each mediator. The mediators not specified in the list will be explained by the original format of the level 2 predictor only.
}
  \item{f01km2}{
the transformation-function-expression list on level 2 predictor (x.j) in explaining the level 2 mediators. The first item lists column numbers of the level 2 mediators in m, which should be explained by the transformed level 2 predictor(s). By that order, each of the rest items of f01km2 lists the transformation functional expressions for the predictor in explaining each mediator. The mediators not specified in the list will be explained by the original format of the (aggregated) level 2 predictor only.
}
  \item{f10km}{
the transformation-function-expression list on level 1 predictor (xij-x.j) in explaining the level 1 mediators. The first item lists column numbers of the level 1 mediators in m, which should be explained by the transformed level 1 predictor(s). By that order, each of the rest items of f10km lists the transformation functional expressions for the level 1 predictor in explaining each mediator. The mediators not specified in the list will be explained by the original format of the level 1 predictor only.
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
  \item{w2}{
the weight for observations at level 2, which should be the same order as unique(level[!is.na(level)]). 
}
  \item{boot}{
the number of bootstrapping samples.
}
  \item{seed}{
set seed, default is 1.
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
  \item{x.new, covariates.new, m.new, level.new, weight.new}{
the settings that we want to make inferences on the mediation effects. If m.new=NULL, generate new mediators from x.new.
}
}
\value{
Return a "mlma.boot" mode list, which include the following items:
  \item{de1}{an n by boot matrix, where each column is the level 1 direct effects from one bootstrap sampling. n is the number of observations in the original data.}
  \item{de2}{an g by boot matrix, where each column is the level 2 direct effects from one bootstrap sampling. g is the number of groups in the original data.}  
  \item{ie1}{an v1 by n*boot matrix, where each column is the level 1 indirect effects from the boot bootstrap samples for one level 1 mediator. v1 is the number of level 1 mediators.}
  \item{ie1}{an v2 by g*boot matrix, where each column is the level 2 indirect effects from the boot bootstrap samples for one level 2 mediator. v2 is the number of level 2 mediators.}
    \item{ie12}{an v1 by g*boot matrix, where each column is the aggregated level 2 indirect effects from the boot bootstrap samples for one level 1 mediator. v1 is the number of level 1 mediators.}
    \item{sum.boot1}{summary results of level 1 mediation effects from bootstrap sample.}
    \item{sum.boot2}{summary results of level 2 mediation effects from bootstrap sample.}
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
temp111<-boot.mlma(y=sim.111$y, biny=FALSE, sim.111$x, levelx=1, 
                   m=sim.111$m, l1=1:2, c1=3,c1r=1,
                   f01y=c("x","log(x^2)"), f10y=c("x^2","sqrt(x+6)"), 
                   f20ky=list(2,c("x","x^3")), f01km1=list(2,"sqrt(x)+3"), 
                   f10km=list(2,"log(x+2)"), level=sim.111$level, boot=2)
  
data(sim.211)
temp211<-boot.mlma(y=sim.211$y, biny=FALSE, x=sim.211$x, levelx=2, m=sim.211$m, 
                   l1=2,l2=1, c1=3,c1r=1, f01y=c("x","log(x^2)"), 
                   f02ky=list(1,c("x","x^2")), f20ky=list(2,c("x","x^3")),
                   f01km1=list(2,"sqrt(x)+3"), f01km2=list(1,c("x^1.2","x^2.3")), 
                   level=sim.211$level, boot=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Multilevel Mediation Analysis}
\keyword{ ~Inferences on MLMA}% __ONLY ONE__ keyword per line
