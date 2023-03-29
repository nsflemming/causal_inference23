#### unbiased 


N <- 1000
x <- runif(N, 0, 1)
y0 <- x + .25*rnorm(N) #potential outcome, control
y1 <- x*2 + .25*rnorm(N) #potential outcome treated
#units have correlated potential outcomes, b/c outcomes related to covars
## almost always the case, treatment vs. control doesn't entirely determine outcome
# but errors are uncorrelated
index <- 1:N
potcov <- data.frame(index, y1,y0,x)

plot(x, y0, col="blue", pch=19, 
     cex=.25, ylim=range(c(y1,y0)),
     main="Potential outcomes")
points(x, y1, col="red", pch=19, cex=.25)



####

#### correlated potential outcomes


N <- 1000
x <- runif(N, 0, 1)
y0 <- x + .25*rnorm(N)
y1 <- x*2 + y0+ .25*rnorm(N)
# y1 is function of y0, so error term from y0 is going into y1
## errors correlated
index <- 1:N
potcov <- data.frame(index, y1,y0,x)

plot(x, y0, col="blue", pch=19, 
     cex=.25, ylim=range(c(y1,y0)),
     main="Potential outcomes")
points(x, y1, col="red", pch=19, cex=.25)



### strong covariance between covariate and unit-level potential outcomes

N <- 1000
x <- runif(N, 0, 2)
y0 <- .5*x - .25*(x-mean(x))^2 +.25*rnorm(N)
y1 <- 2*y0 + .25 + (x-mean(x))^2 + .25*rnorm(N)
#non-linearity in contribution of x to y's
index <- 1:N
potcov <- data.frame(index, y1,y0,x)

plot(x, y0, col="blue", pch=19, 
     cex=.25, ylim=range(c(y1,y0)),
     main="Potential outcomes")
points(x, y1, col="red", pch=19, cex=.25)



plot(x, y1-y0,  pch=19, 
     cex=.25, 
     main="Potential outcomes")


plot((x-mean(x))^2, y1-y0,  pch=19, 
     cex=.25, 
     main="Potential outcomes")

######

#########BOOTSTRAP example for covaraite adjustment




n.vec <- c(20,40,80,160,320)

nsim <- 2000

# vectors to store results
sate.hat.noadj <- sate.hat.simpadj <- sate.hat.intadj <- sate <- matrix(NA, ncol=length(n.vec), nrow=nsim)

for(j in 1:length(n.vec)){
  n <- n.vec[j] #iterate through sample sizes
  
  for(i in 1:nsim){
    
    # Simple random sample without replacement
    potcov.s <- potcov[potcov$index %in% sample(potcov$index, n),] #potential outcomes
    
    sate[i,j] <- mean(potcov.s$y1 - potcov.s$y0) #sample avg treatment effect = mean of diff in potential outcomes
    
    # Complete random assignment of treatment
    n1 <- floor(.33*n)  # you can play around with treatment allocation.
    # this is set to an imbalanced design.
    potcov.s$D <- 0
    potcov.s$D[potcov.s$index %in% sample(potcov.s$index, n1)] <- 1 #assign treatment
    
    potcov.s$Y <- with(potcov.s, D*y1 + (1-D)*y0) #calculate outcome based on assignment
    
    # No adjustment (difference in means)
    ols.1 <- lm(Y~D, data=potcov.s)
    
    # Simple covariate adjustment, 
    ols.2 <- lm(Y~D+x, data=potcov.s)
    
#adding coefficients increases efficiency
    sate.hat.noadj[i,j] <- coef(ols.1)["D"]
    sate.hat.simpadj[i,j] <- coef(ols.2)["D"]
  }
}

se.sate.hat.noadj <- apply(sate.hat.noadj, 2, sd) 
bias.sate.hat.noadj <- apply(sate.hat.noadj-sate, 2, mean) #calc bias by subtract true value of SATE
std.bias.sate.hat.noadj <- bias.sate.hat.noadj/se.sate.hat.noadj #divide bias by std error
  #lowered std error by adding covars
se.sate.hat.simpadj <- apply(sate.hat.simpadj, 2, sd)
bias.sate.hat.simpadj <- apply(sate.hat.simpadj-sate, 2, mean)
std.bias.sate.hat.simpadj <- bias.sate.hat.simpadj/se.sate.hat.simpadj

res.tab <- cbind(c("SE","Bias","Bias/SE","SE","Bias","Bias/SE"),
                 round(1000*rbind(   se.sate.hat.noadj,
                                     bias.sate.hat.noadj,
                                     std.bias.sate.hat.noadj,
                                     se.sate.hat.simpadj,
                                     bias.sate.hat.simpadj,
                                     std.bias.sate.hat.simpadj),0)/1000)
res.tab.out <- cbind(c("No adjustment","","",
                       "Simple adjustment","",""), res.tab)
colnames(res.tab.out) <- c("Estimator","Statistic",paste("N=",n.vec,sep=""))


res.tab.out

#bias starts small, but declines a bit as sample size increases
#standard error goes down fast as sample size increases
#simple adjustment: std errors always smaller than no adjustement; add covars = smaller std errors
## adding adjustment should increase bias, but bias/se still better, cause SE reduction much bigger than bias increase
  ### that last bit didn't work for some reason
######




########################################


##########################################Replication of Munger (2017)
#tweetment effects paper

####conceptually, what is the motivation for controlling for log.followers? 
#### is there some risk of this being a ``bad control" ? 

# num of followers each subj has. introducing bias?
  # controlling increases efficiency
#affects reaction to treatment. popular user less likely to be dissuaded?

#Is there a relationship between the treatment and this covar? Is there relationship b/t this covar & outcome
## if treatment causes covar then it might be post-treatment var, so controlling bad
## if a collider then controlling will increase bias
## neither is the case here


#### The two datasets are based on different assumptions about how to treat 
### study attrition --- subjects whose accounts were deleted/suspended/unmeasured

### "standard" drops those subjects and proceeds as normal
### "conservative" says assign each of those subjects a treatment effect of 0
  ### assume the treatment effect was 0, racism score constant pre to post

### Explain the potential bias from using the standard assumption

###### The treatment may have caused the attrition. 
##     So dropping them is conditioning on a post treatment variable (whether they stick around)?
## conservative says attrition (sample selection bias) is related to treatment
  ## say sample selection is correlated w/ potential outcomes
## standard: may be collider (being in sample is collider)


###using the code from inclass_215.R, perform CEM on this data, matching
## on the anonymity, pre-treatment outcome measure, and log.followers.

### You'll have to restrict the data to just treatment 3 and the control (0)

## estimate a standard ols model with these cem weights. 

<<<<<<< HEAD
library(tidyverse)
library(fixest)
library(haven)
library(MatchIt)

=======
>>>>>>> 9ce91ef1d1ed90daaa209797310652e296d3d1b0


# Coarsened-Exact Matching

<<<<<<< HEAD
data_t3<-data[data$treat.f==3|data$treat.f==0,]

cem_out <- matchit(
  treat.f ~ anonymity + racism.scores.pre.2mon + log.followers,
  data = data_t3,
  method = "cem", estimand = "ATT"
)

data_t3$cem_weights = cem_out$weights

feols(
  racism.scores.post.1wk ~ i(treat.f), weights = ~cem_weights,
  data = data_t3,  vcov = "hc1"
)
=======
data$treat<-ifelse(data$treat.f == 3, 1, 0)
data_binary<-filter(data, treat == 1 | treat.f == 0)

cem_out <- matchit(
  treat ~  racism.scores.pre.2mon + log.followers + anonymity,
  data = data_binary,
  method = "cem", estimand = "ATE"
)

data_binary$cem_weights = cem_out$weights

feols(
  racism.scores.post.1wk ~ treat, weights = ~cem_weights,
  data = data_binary,  vcov = "hc1")



summary(lm(
  racism.scores.post.2wk ~ treat, weights = cem_weights,
  data = data_binary
  
))
>>>>>>> 9ce91ef1d1ed90daaa209797310652e296d3d1b0

