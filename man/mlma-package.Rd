\name{mlma-package}
\alias{mlma-package}
\docType{package}
\title{
Multilevel Mediation Analysis
}
\description{
The package is used to do mediation analysis with generalized multilevel models.
}
\details{
\code{"\link[=data.org]{data.org}"} is used to transform the variables and organize the predictor, mediators and outcome into the format that are ready to be used for multilevel mediation analysis. \code{"\link[=mlma]{mlma}"} is for multilevel mediation analysis on the original data set. \code{"\link[=boot.mlma]{boot.mlma}"} is a combined function that organized data set, do multilevel mediation analysis on original data sets and bootstrapping samples.

The multilevel mediation is based on the following linear multilevel additive models:
\deqn{Y_{ij} = u_{0j}^Y(X_{.j}, \mathbf{M}_{.j}, \mathbf{Z}_{.j})+{\boldsymbol{\beta}_{10}^Y}^T\mathbf{f}_{10}^Y(X_{ij}-X_{.j})+\sum_{k=1}^K{\boldsymbol{\beta}_{20k}^Y}^T\mathbf{f}_{20k}^Y(M_{ijk}-M_{.jk})+{\boldsymbol{\beta}_{30}^Y}^T(\mathbf{Z}_{ij}-\mathbf{Z}_{.j})+r_{ij}^Y,}
where \deqn{u_{0j}^Y(X_{.j}, \mathbf{M}_{.j}, \mathbf{Z}_{.j}) = c_{00}^Y + {\boldsymbol{\beta}_{01}^Y}^T\mathbf{f}_{01}^Y(X_{.j}) + \sum_{k=1}^K{\boldsymbol{\beta}_{02k}^Y}^T\mathbf{f}_{02k}^Y(M_{.jk}) + {\boldsymbol{\beta}_{03}^Y}^T\mathbf{Z}_{.j} + r_{0j}^Y.}
For \eqn{k=1,\ldots,K,}
\deqn{M_{.jk} = u_{0jk}^M(X_{.j})+{\boldsymbol{\beta}_{10k}^M}^T\mathbf{f}_{10k}^M(X_{ij}-X_{.j})+r_{ijk}^M,}
\deqn{u_{0jk}^M(X_{.j}) = c_{00k}^M + {\boldsymbol{\beta}_{01k}^M}^T\mathbf{f}_{01k}^{M1}(X_{.j}) + r_{0jk}^M.}
If for some k, \eqn{M_k} is level 2 variable,
\deqn{M_{.jk} = c_{00k}^M + {\boldsymbol{\beta}_{01k}^M}^T\mathbf{f}_{01k}^{M2}(X_{.j}) + r_{0jk}^M.} 

Note that in the models, \eqn{\mathbf{f}(\cdot)=(f_1(\cdot), f_2(\cdot), \cdots, f_l(\cdot))^T} is a set of l transformation functions on \eqn{\cdot}, with the corresponding linear coefficients vector \eqn{\boldsymbol{\beta}=(\beta_1, \beta_2, \cdots, \beta_l)^T}. \eqn{\mathbf{f}} and l are known for model fitting. l may be different with \eqn{\mathbf{f}} of different sub- and super-scripts.}
\author{
Qingzhao Yu  \email{qyu@lsuhsc.edu}, Bin Li \email{bli@lsu.edu}

Maintainer: Qingzhao Yu  \email{qyu@lsuhsc.edu}
}
\references{
Yu, Q. and Li, B., (2020). <doi:10.1371/journal.pone.0241072>. "Third-Variable Effect Analysis with Multilevel Additive Models," PLoS ONE 15(10): e0241072.

Yu, Q., Yu, M., Zou, J., Wu, X., Gomez, SL, Li, B. (2021). <doi:10.1177/26320843211061292>. "Multilevel Mediation Analysis on Time-to-Event Outcomes - Exploring racial/ethnic Disparities in Breast Cancer Survival in California," Research Methods in Medicine & Health Sciences. 

Yu, Q. and Li, B., 2022. Statistical Methods for Mediation, Confounding and Moderation Analysis Using R and SAS. Chapman and Hall/CRC. ISBN 9780367365479.
}