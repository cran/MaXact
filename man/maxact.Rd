\name{Exact Max CATT Test}
\alias{MaXact}
\alias{maxact.test}
\alias{catt.test}
\title{MAX CATT test in case-control genetic association analysis }
\description{
  Perform exact MAX3 test, MAX2 test and Cochran-Armitage trend test for
  one-locus genetic association analysis.  MAX3 is the maximum of the
  Cochran-Armitage trend test statistic for dominant, recessive and
  additive models, and MAX2 is the maximum of dominant and recessive
  models.
  
  It also provides approximated p-value of a MAX test with the normal
  approximation method.
}
\usage{
maxact.test(data, max3 = TRUE, exact=TRUE, alternative = "two.sided")
catt.test(data, theta, exact=TRUE, alternative = "two.sided")
}
\arguments{
  \item{data}{ 2x3 contingency table of data values, with rows
    representing case/control, and columns representing the three
    genotypes }
  \item{max3}{ a logical indicating whether MAX3 or MAX2 test statistic
    should be used }
  \item{theta}{ parameter of the Cochran-Armitage trend test } 
  \item{exact}{ a logical indicating whether an exact p-value should be
  computed }
  \item{alternative}{ a character string describing the alternative
    hypothesis, which must be one of '"two.sided"', '"greater"' or
    '"less"'}
}
\value{
  A list with class '"htest"' containing the following components:
  
  \item{data.name}{a character string giving the names of the data.}
  \item{p.value}{the p-value of the test.}
  \item{method}{a character string describing the method used.}
  \item{statistic}{the value of the test statistic.}
  \item{alternative}{ a character string describing the alternative hypothesis.}
}
\note{
  The normal approximation method(exact=FALSE) is provided for
  comparison purpose, and is not recommended for data analysis. The
  exact method is faster than the normal approximation method in most of
  cases, and it is fast enough for most of applications.
}
\references{ Tian, J., Xu, C., Zhan, H., Yang, Y. (2009) "Exact MAX
tests in case-control association analysis (Manuscript)." }
\examples{
data.sladek <- matrix(c(129, 326, 229, 198, 325, 143), 2, 3, byrow=TRUE)

## MAX3 test (exact p-value)
maxact.test(data.sladek) 

## MAX2 test (exact p-value):
maxact.test(data.sladek, max3=FALSE, exact=TRUE)


## MAX test (p-value calculated from the normal approximation method), it
## is not recommends, see NOTE for more information
 maxact.test(data.sladek, max3=TRUE, exact=FALSE)
 maxact.test(data.sladek, max3=FALSE, exact=FALSE)

## Cochran-Armitage trend test, assuming we know the model 
 catt.test(data.sladek, theta=0)    #recessive model
 catt.test(data.sladek, theta=0.5)  #additive model
 catt.test(data.sladek, theta=1)    #dominant model

#the  alternative hypothesis is less
 catt.test(data.sladek, theta=0, alternative="less") 


}
\author{ Jianan Tian \email{jiant17@gmail.com}, Chenliang Xu \email{lucky.random@gmail.com} }
\keyword{ htest }
