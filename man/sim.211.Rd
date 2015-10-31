\name{sim.211}
\alias{sim.211}
\docType{data}
\title{
Simulated 2-1-1 Data
}
\description{
A simulated data set, where the predictor is a level 2 and the outcome is a level 1 variable.
}
\usage{data("sim.211")}
\format{
  The data set contains 10 groups, each group has 30 observations. The format is list, where there are four elements:
  
   x: the level 1 continuous predictor.
   
   y: the level 1 continuous outcome.
   
   m: the matrix of mediators, where there are one level 2 mediator, m.1, and two level 1 mediators, m.3 and m.4. m.4 is categorical with 3 levels.
   
   level: the group level for each observation.
}
\examples{
data(sim.211)
}
\keyword{datasets}
