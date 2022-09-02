#' Bootstrap test for equality of two (2) population variances
#'
#'
#' Testing equality of two (2) population variances against the alternative that both variances are not equal.
#'
#'
#' @param x1 first sample vector of data or observations  
#' @param x2 second sample vector of data or observations 
#' @param a significance level alpha
#' @param B number of bootstrap samples. At least 500 is recommended.
#'
#'
#'
#' @return list consisting of a non-numeric decision whether to reject the null hypothesis or not, the significance level, the number of bootstrap samples used, and the bootstrap P-value  calculated using the Euclidean distance.  
#'
#'
#'
#' @references
#'
#' Cahoy, DO (2010), \emph{A Bootstrap Test For Equality Of Variances,} Computational Statistics & Data Analysis, 54(10), 2306-2316.
#' <doi:10.1016/j.csda.2010.04.012>
#'
#'
#'
#' @examples
#'
#'
#' x1=sqrt(10)*runif(8, -sqrt(3), sqrt(3) )
#' x2=sqrt(1)*runif(8, -sqrt(3), sqrt(3) )
#' equa2vartest(x1,x2,0.05, 1000)
#'
#'
#'
#' x1=sqrt(1)*rexp(8)
#' x2=sqrt(1)*rexp(8)
#' equa2vartest(x1,x2,0.01, 1000)
#'
#'
#'
#' @import stats
#'
#'
#' @export
equa2vartest=function(x1, x2, a, B){
  k<-1
  x1<-na.omit(x1)
  x2<-na.omit(x2)
  y<-matrix(0,nrow=B,ncol=k+1)
  vary1<-matrix(0, nrow=B)
  vary2<-matrix(0, nrow=B)
  
  e1<-x1-mean(x1)
  e2<-x2-mean(x2)
  n1<-length(x1)
  n2<-length(x2)
  N<-(n1 + n2)
  nu2<-sum( c(  ( x1 -mean(x1) )^2 , ( x2-mean(x2) )^2  ) )/N
  ss<-c(e1,e2)
  b3<- N*sum(ss^4)/ (sum( ss^2)^2)
  vlogs21c <-(1/(n1-1))*( b3-(n1-3)/n1 )
  vlogs22c <-(1/(n2-1))*( b3-(n2-3)/n2 )
  vary1c<- (1-2/(k+1) )*vlogs21c + (1/(k+1)^2)*(vlogs21c+vlogs22c)
  vary2c<- (1-2/(k+1) )*vlogs22c + (1/(k+1)^2)*(vlogs21c+vlogs22c)
  s1<-var(e1)
  s2<-var(e2)
  yy1<-( log(s1/(s1*s2)^(1/(k+1))) )
  yy2<-( log(s2/(s1*s2)^(1/(k+1))) )
  for(j in 1:B){
    ss11<-sample(x1, n1, replace=TRUE)
    ss12<-sample(x2, n2, replace=TRUE)
    b<- N*sum( c(ss11-mean(ss11),ss12-mean(ss12) )^4 )/ ( ( sum ( c( ss11-mean(ss11),ss12-mean(ss12)  )^2) )^2)
    vlogs21<-(1/(n1-1))*( b-(n1-3)/n1 )
    vlogs22<-(1/(n2-1))*( b-(n2-3)/n2 )
    vary1[j]<- (1-2/(k+1) )*vlogs21 + (1/(k+1)^2)*(vlogs21+vlogs22)
    vary2[j]<- (1-2/(k+1) )*vlogs22 + (1/(k+1)^2)*(vlogs21+vlogs22)
    s1a<-var(ss11)
    s2a<-var(ss12)
    y[j,1]<-( log(s1a/(s1a*s2a)^(1/(k+1)))   )
    y[j,2]<-( log(s2a/(s1a*s2a)^(1/(k+1)))   )
  }
  y22<-cbind(y,sqrt(vary1),sqrt(vary2) )
  y<-y22[complete.cases(y22),]
  y<-(y[,1:(k+1)])/y[, ((k+1)+1):(2*(k+1))]
  avey<-colMeans(y)
  ycntr<-(y-matrix(rep(avey,length(y[,1])), nrow=length(y[,1]), byrow=T) )
  y<-ycntr[,1]
  y<- y[abs(y)<Inf]
  B2<-length(y)
  cval<- sort( abs(y), decreasing=T)[ceiling(B2*a) + 1]
  ts2<- c( (yy1/sqrt(vary1c)),(yy2/sqrt(vary2c)) )
  ind <- ifelse( sum( ifelse( ( abs(ts2) > rep( cval ,(k+1) ) ), 1, 0) ) >0,1,0 )
  dum=ifelse( (abs(ycntr)>abs(ts2))==c("TRUE", "TRUE"), 1,0)
  pval<-round((sum(apply(dum,1,sum)==2)+1)/B2,3)
  return(list(ifelse(ind==0,"Decision: Fail to Reject the Null", "Decision: Reject the Null"), Alpha=a, NumberOfBootSamples=B, BootPvalue=pval ) )
}
