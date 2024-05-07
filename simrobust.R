library(xts)
install.packages("rmgarch")

#library(qrmdata)
library(rmgarch)
library(nlme)
data <- read.csv("dowjonesusr1.csv") #use min 500 datapoints
a=dim(data)
n=a[2]
simreturns=matrix(0,100*50,n)
m=1
for(val in seq(1,946,63))
{
a=val
b=val+499
data1<-data[a:b,]
## Load some real data
#data("FTSE")
#data("SMI")

#INDEXES <- merge(FTSE, SMI, all = FALSE)
#plot.zoo(INDEXES)

## Compute returns
#FTSE.X <- diff(log(FTSE))[-1]
#SMI.X <- diff(log(SMI))[-1]
#INDEXES.X <- merge(FTSE.X, SMI.X, all = FALSE)
#plot.zoo(INDEXES.X)

## Take 4 years of data
#data <- INDEXES.X['2006-01-01/2009-12-31']
#pairs(as.zoo(data))
#dim(data)

xspec = ugarchspec(mean.model = list(armaOrder = c(1, 0)), variance.model = list(garchOrder = c(1,1), model = 'gjrGARCH'), distribution.model = 'norm')
#fit=ugarchfit(xspec,data=data)
#likelihood(fit)
uspec = multispec(replicate(n, xspec))
spec1 = dccspec(uspec = uspec, dccOrder = c(1, 0), distribution = 'mvnorm')

#fit1 = dccfit(spec1, data = data, fit.control = list(eval.se = TRUE))
fit1 = dccfit(spec1, data = data1)
likelihood(fit1)
c=100*(m-1)+1
d=100*m
for(k in c:d)
{
aa=dccsim(fit1, n.sim = 500)
simualted_data=fitted(aa)
mean1=colMeans(simualted_data)
simreturns[k,]=mean1
}
m=m+1
print(m)
}
write.csv(simreturns,"dowjonesusfinalsimout63.csv")

