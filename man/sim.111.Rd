\name{sim.111}
\alias{sim.111}
\docType{data}
\title{
Simulated 1-1-1 Data set 
}
\description{
A simulated data set, where both predictor and outcome are level 1 variables.
}
\usage{data("sim.111")}
\format{
  The data set contains 10 groups, each group has 30 observations. The format is list, where there are four elements:
  
   x: the level 1 continuous predictor.
   
   y: the level 1 continuous outcome.
   
   m: the matrix of mediators, where there are three level 1 mediators, where m.2 and m.3 are continuous, and m.4 is categorical with 3 levels.
   
   level: the group level for each observation.
}
\examples{
data(sim.111)
}
\keyword{datasets}
