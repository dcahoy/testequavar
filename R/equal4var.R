#' Bootstrap test for equality of four (4) population variances
#'
#'
#' Testing equality of four (4) population variances against the alternative that all variances are not equal.
#'
#'
#' @param x1 first sample vector of data or observations
#' @param x2 second sample vector of data or observations
#' @param x3 third sample vector of data or observations
#' @param x4 fourth sample vector of data or observations
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
#' x1=sqrt(10)*runif(10, -sqrt(3), sqrt(3) )
#' x2=sqrt(1)*runif(10, -sqrt(3), sqrt(3) )
#' x3=sqrt(1)*runif(10, -sqrt(3), sqrt(3) )
#' x4=sqrt(1)*runif(10, -sqrt(3), sqrt(3) )
#' equa4vartest(x1,x2,x3, x4, a=0.05, B=500)
#'
#'
#'
#' equa4vartest(rexp(10) ,rexp(10) ,rexp(10) , rexp(10),  a=0.01, B=1000)
#'
#'
#'
#' @import stats
#'
#'
#' @export
equa4vartest=function(x1, x2, x3, x4, a, B){
  k=3
  x1<-na.omit(x1)
  x2<-na.omit(x2)
  x3<-na.omit(x3)
  x4<-na.omit(x4) 
  y<-matrix(0,nrow=B,ncol=k+1)
  vary1<-matrix(0, nrow=B)
  vary2<-matrix(0, nrow=B)
  vary3<-matrix(0, nrow=B)
  vary4<-matrix(0, nrow=B)
  e1<-x1-mean(x1)
  e2<-x2-mean(x2)
  e3<-x3-mean(x3)
  e4<-x4-mean(x4)
  n1<-length(x1)
  n2<-length(x2)
  n3<-length(x3)
  n4<-length(x4)
  N<-(n1 + n2 + n3 + n4)
  nu2<-sum( c(  ( x1 -mean(x1) )^2 , ( x2-mean(x2) )^2 , ( x3-mean(x3) )^2, ( x4 -mean(x4) )^2   ) )/N
  ss<-c(e1,e2,e3, e4)
  b3<- N*sum(ss^4)/ (sum( ss^2)^2)
  vlogs21c <-(1/(n1-1))*( b3-(n1-3)/n1 )
  vlogs22c<-(1/(n2-1))*( b3-(n2-3)/n2 )
  vlogs23c<-(1/(n3-1))*( b3-(n3-3)/n3 )
  vlogs24c<-(1/(n4-1))*( b3-(n4-3)/n4 )
  vary1c <- (1-2/(k+1) )*vlogs21c + (1/(k+1)^2)*(vlogs21c+vlogs22c+vlogs23c+vlogs24c)
  vary2c <- (1-2/(k+1) )*vlogs22c + (1/(k+1)^2)*(vlogs21c+vlogs22c+vlogs23c+vlogs24c)
  vary3c <- (1-2/(k+1) )*vlogs23c + (1/(k+1)^2)*(vlogs21c+vlogs22c+vlogs23c+vlogs24c)
  vary4c<- (1-2/(k+1) )*vlogs24c + (1/(k+1)^2)*(vlogs21c+vlogs22c+vlogs23c+vlogs24c)
  s1<-var(e1)
  s2<-var(e2)
  s3<-var(e3)
  s4<-var(e4)
  yy1<-( log(s1/(s1*s2*s3*s4)^(1/(k+1))) )
  yy2<-( log(s2/(s1*s2*s3*s4)^(1/(k+1))) )
  yy3<-( log(s3/(s1*s2*s3*s4)^(1/(k+1))) )
  yy4<-( log(s4/(s1*s2*s3*s4)^(1/(k+1))) )
  for(j in 1:B){
    ss11<-sample(x1, n1, replace=TRUE)
    ss12<-sample(x2, n2, replace=TRUE)
    ss13<-sample(x3, n3, replace=TRUE)
    ss14<-sample(x4, n4, replace=TRUE)
    b <- N*sum( c(ss11-mean(ss11),ss12-mean(ss12),ss13-mean(ss13), ss14-mean(ss14) )^4 )/ ( ( sum (c(ss11-mean(ss11),ss12-mean(ss12),ss13-mean(ss13), ss14-mean(ss14))^2) )^2)
    vlogs21<-(1/(n1-1))*( b-(n1-3)/n1 )
    vlogs22<-(1/(n2-1))*( b-(n2-3)/n2 )
    vlogs23<-(1/(n3-1))*( b-(n3-3)/n3 )
    vlogs24<-(1/(n4-1))*( b-(n4-3)/n4 )
    vary1[j]<- (1-2/(k+1) )*vlogs21 + (1/(k+1)^2)*(vlogs21+vlogs22+vlogs23+vlogs24)
    vary2[j]<- (1-2/(k+1) )*vlogs22 + (1/(k+1)^2)*(vlogs21+vlogs22+vlogs23+vlogs24)
    vary3[j]<- (1-2/(k+1) )*vlogs23 + (1/(k+1)^2)*(vlogs21+vlogs22+vlogs23+vlogs24)
    vary4[j]<- (1-2/(k+1) )*vlogs24 + (1/(k+1)^2)*(vlogs21+vlogs22+vlogs23+vlogs24)
    s1a<-var(ss11)
    s2a<-var(ss12)
    s3a<-var(ss13)
    s4a<-var(ss14)
    y[j,1]<-( log(s1a/(s1a*s2a*s3a*s4a)^(1/(k+1)))   )
    y[j,2]<-( log(s2a/(s1a*s2a*s3a*s4a)^(1/(k+1)))   )
    y[j,3]<-( log(s3a/(s1a*s2a*s3a*s4a)^(1/(k+1)))   )
    y[j,4]<-( log(s4a/(s1a*s2a*s3a*s4a)^(1/(k+1)))   )
  }
  y22<-cbind(y,sqrt(vary1),sqrt(vary2),sqrt(vary3),sqrt(vary4)  )

  y<-y22[complete.cases(y22),]
  B1<-dim(y)[1]
  y<-(y[,1:(k+1)])/y[, ((k+1)+1):(2*(k+1))]
  avey<-colMeans(y)
  y<-(y-matrix(rep(avey,length(y[,1])), nrow=length(y[,1]), byrow=T) )
  ys<-sort(abs(y), decreasing=T)
  kr<-as.integer(round( (k+1)*.3*B1))
  for(r in 1:kr){
    if( as.integer(round(sum(ifelse( apply( ifelse( (  ( y <= matrix( rep (ys[r],(k+1)*B1), nrow=B1)  ) & ( y >= matrix( rep (-ys[r],(k+1)*B1), nrow=B1) )  ) , 1, 0),1,sum) ==(k+1), 1,0)) ))==as.integer(floor(B1*(1-a)))) break
  }
  cval<- ys[r]
  ts2<- c( (yy1/sqrt(vary1c)),(yy2/sqrt(vary2c)),(yy3/sqrt(vary3c)),(yy4/sqrt(vary4c)) )
  ind<- ifelse( sum( ifelse( ( abs(ts2) > rep( cval ,(k+1) ) ), 1, 0)    ) >0,1,0 )
  euclid <- function(a) sqrt(sum((a)^2))
  dum=apply(y,1,euclid) 
  dum2=euclid(ts2)
  pval<-round((sum(dum>dum2)+1)/B1,3)
  return(list(ifelse(ind==0,"Decision: Fail to Reject the Null", "Decision: Reject the Null"),  Alpha=a, NumberOfBootSamples=B, BootPvalue=pval ) )
}

