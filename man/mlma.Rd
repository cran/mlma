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
mlma(y, biny = FALSE, data1, x, levelx, levely, m, l1 = NULL, l2 = NULL, c1 = NULL, 
     c1r = rep(1, length(c1)), c2 = NULL, c2r = rep(1, length(c2)), level, 
     weight = rep(1, length(x)), random = "(1|level)", random.m1 = NULL, intercept = TRUE,
     covariates=NULL, cy1=NULL, cy2=NULL, cm=NULL, joint=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
the vector of the outcome variable.
}
  \item{biny}{
True if the outcome is binary, otherwise False. Default is False.
}
  \item{data1}{The transformed and organized data set from data.org}
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
the random effect part for model explaining the mediators. All other random effects are random = "(1|level)" if not specified here.
}
  \item{intercept}{
True if fit an intercept to models, by default.
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
the column numbers of covariates that are used to explain m. cm[[1]] gives the mediators (in l1, cl, l2, or c2) that can be partially explained by covariates. Each of the rest items of the cm list shows the column number(s) in covariates that should be used to explain each mediator listed in cm[[1]] and by that order.
}
  \item{joint}{
the list of group(s) of mediators whose joint mediation effect is of interests. joint[[1]] list the levels of mediators in each group and by the order of the list. Note that if any mediator in the group is of level 2, the level of the group should be 2.
}
}
\value{
A "mlma" mode list will be returned with the following items:
  \item{de1}{level 1 direct effect from predictor.}
  \item{de2}{level 2 direct effect from predictor.}
  \item{ie1}{level 1 indirect effect from level 1 mediator.}
  \item{ie12}{level 2 indirect effect from level 1 mediator.}
  \item{ie2}{level 2 indirect effect from level 2 mediator.}
  \item{f1}{the overall multilevel model.}
  \item{fm1}{a list, where the first item identifies the level 1 mediators, and in that order, the following items give the prediction functions of the mediators.}
  \item{fm2}{a list, where the first item identifies the level 2 mediators, and in that order, the following items give the prediction functions of the mediators.}
  \item{ie1_list}{a list, where the first item identifies the level 1 mediators, and in that order, the following items give the column(s) of the indirect effects of the mediator in ie1.}
  \item{ie2_list}{a list, where the first item identifies the level 2 mediators, and in that order, the following items give the column(s) of the indirect effects of the mediator in ie2.}
  \item{iej2_list}{a list, where the first item identifies the level 2 joint mediators, and in that order, the following items give the column(s) of the indirect effects of the mediator in cbind(ie12,ie2).}
  \item{ie12_1,ie12_2, ie1_1, ie1_2, ie2_1, ie2_2}{the first and second part of the corresponding indirect effects.}
  \item{x}{the vector of the predictive variable.}
  \item{x.j}{the vector of the aggregated variable at the higher level by the order of unique(level[!is.na(level)]).}
  \item{levelx}{the level of the predictive variable.}
  \item{weight}{the weight for each observation.}
}
\author{
Qingzhao Yu (qyu@lsuhsc.edu), Bin Li (bli@lsu.edu).
}
\examples{
data(sim.111)
data1<-data.org(x=sim.111$x, levelx=1, m=sim.111$m, l1=1:2, c1=3, c1r=1, 
                f01y=c("x","log(x^2)"), f10y=c("x^2","sqrt(x+6)"), 
                f20ky=list(2,c("x","x^3")), f01km1=list(2,"sqrt(x)+3"), 
                f10km=list(2,"log(x+2)"), level=sim.111$level)
temp<-mlma(y=sim.111$y, biny=FALSE, data1=data1, x=sim.111$x, levelx=1, m=sim.111$m, 
                l1=1:2,c1=3, c1r=1,level=sim.111$level)

data(sim.211)
data1<-data.org(x=sim.211$x, levelx=2, m=sim.211$m, l1=2,l2=1, c1=3, c1r=1,
                   f01y=c("x","log(x^2)"), f02ky=list(1,c("x","x^2")), 
                   f20ky=list(2,c("x","x^3")), f01km1=list(2,"sqrt(x)+3"), 
                   f01km2=list(1,c("x^1.2","x^2.3")), level=sim.211$level)
temp<-mlma(y=sim.211$y, biny=FALSE, data1, x=sim.211$x, levelx=2, m=sim.211$m, 
           l1=2, l2=1, c1=3, c1r=1,level=sim.211$level)                               
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Multilevel Mediation Analysis}
