\name{Exact Max CATT Test}
\alias{MaXact}
\alias{maxact.test}
\alias{catt.test}
\title{MAX CATT test in case-control genetic association analysis }
\description{
  Perform exact MAX3 test, MAX2 test or Cochran-Armitage trend test for
  one-locus genetic association analysis and trend test for dominant,
  recessive and additive models.  It can also calculate approximated
  p-value of a MAX test using the normal approximation method. MAX3 is
  the maximized chi-squares of Cochran-Armitage trend tests for
  dominant, recessive and additive models. MAX2 is the maximization of
  the former two chi-squares.
}
\usage{
maxact.test(data, max3 = TRUE, exact=TRUE, alternative = "two.sided")
catt.test(data, theta, exact=TRUE, alternative = "two.sided")
}
\arguments{
  \item{data}{ 2x3 contingency table(rows are case/control, columns are
    the three genotypes) }
  \item{max3}{ a logical. If 'TRUE'(the default) the MAX3 is returned;
    if 'FALSE' MAX2 is computed }
  \item{theta}{ parameter theta of Cochran-Armitage trend test } 
  \item{exact}{ a logical indicating if an exact method should be used
    (default is 'TRUE') }
  \item{alternative}{ indicates the alternative hypothesis and must be
    one of '"two.sided"', '"greater"' or '"less"'}
}
\value{
  A list with class '"htest"' containing the following components:
  
  \item{data.name}{a character string giving the names of the data.}
  \item{p.value}{the p-value of the test.}
  \item{method}{a character string describing the method used.}
  \item{statistic}{the value of the MAX test statistic.}
  \item{alternative}{ a character string describing the alternative hypothesis.}
}
\note{
  The normal approximation method(exact=FALSE) is not recommended,
  unless you have really good reason for that. It is not exact, and it
  would not be faster than the exact method for normal sample size (such
  as thousands).
  }
\details{
  
}
\references{ Tian, J., Xu, C., Zhan, H., Yang, Y. (2009) "Exact MAX
tests in case-control association analysis (Manascript)." }
\seealso{
}
\examples{
data.sladek<-matrix(c(129,326,229,198,325,143),2,3,byrow=TRUE)

## MAX3 test (exact p-value)
maxact.test(data.sladek) 

## MAX2 test (exact p-value):
maxact.test(data.sladek,max3=FALSE,exact=TRUE)


## MAX test (p-value calculated from the normal approximation method), it
## is not recommends, see NOTE for more information
 maxact.test(data.sladek,max3=TRUE,exact=FALSE)
 maxact.test(data.sladek,max3=FALSE,exact=FALSE)

## Cochran-Armitage trend test, assuming we know the model 
 catt.test(data.sladek, theta=0)    #recessive model
 catt.test(data.sladek, theta=0.5)  #additive model
 catt.test(data.sladek, theta=1)    #dominant model

#the  alternative hypothesis is less
 catt.test(data.sladek, theta=0, alternative="less") 


}
\author{ Jianan Tian \email{jiant17@gmail.com}, Chenliang Xu \email{lucky.random@gmail.com} }
\keyword{ htest }
