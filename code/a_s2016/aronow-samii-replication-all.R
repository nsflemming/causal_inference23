rm(list=ls())
# change wd as needed:
setwd("C:\Users\natha\OneDrive\Documents\GitHub\causal_inference23\code\a_s2016")

library(foreign)
library(maptools)

# Jensen example

# Files needed:
# 5603 folder with world shapefile and auxiliary files
# jensen-rep.dta
# mapnames_filled.csv

world <-  readShapePoly("world_countries_boundary_file_world_2002.shp")
jensen.cc <- read.dta("jensen-rep.dta") 

X.vars <- c(	"var5",
				"market",
				"lgdppc",
				"gdpgrowt",
				"tradeofg",
				"overallb",
				"generalg",
				"country",
				"d2",
				"d3")
X.vars.f <- paste(X.vars,collapse="+")
fit.y <- lm(as.formula(paste("Fvar5~regime+", X.vars.f, sep="")), data=jensen.cc) #predict democracy
fit.d <- lm(as.formula(paste("regime~", X.vars.f, sep="")), data=jensen.cc) #predict regime with covars
summary(fit.y)
summary(fit.d) #several variables are highly predictive of the treatment
#weights are constructed based on the residuals from this partial regression

d.tilde <- as.numeric(residuals(fit.d)) #take residuals
w <- d.tilde^2 #square residuals to get weights

w1 <- tapply(w, jensen.cc$country, mean) #apply mean to weight vector (avg weight for each country across years)
mapnames <-  read.csv("mapnames_filled.csv") #names of countries in world
mapnames$weight <- 0
mapnames$weight[match(rownames(w1), mapnames$jensen)] <- as.numeric(w1) #replace weight with new weights with row names

attributes(world)$data$weight <- 0
attributes(world)$data$weight[match(mapnames$mapname,attributes(world)$data$NAME)] <- mapnames$weight

attributes(world)$data$incl <- 0
attributes(world)$data$incl[match(mapnames$mapname,attributes(world)$data$NAME)] <- as.numeric(!is.na(mapnames$jensen))





# Sample

pdf(file="jensen-nominal-map.pdf", height=5, width=8)
plot(	world, 
		col=gray(1-.75*attributes(world)$data$incl),
		border="gray", lwd=.25)
dev.off()


# Effective sample

pdf(file="jensen-effective-map.pdf", height=5, width=8)
plot(	world, 
		col=gray(1-abs(attributes(world)$data$weight)/max(abs(attributes(world)$data$weight))),
		border="gray", lwd=.25)
dev.off()


########

##in-class


#### What are the top 10 individual observations in terms of weight?
jensen.cc$weight<-w
library(dplyr)
jensen.desc<-jensen.cc %>% arrange(desc(weight))
jensen.desc[1:10,]

#### WHat is the minimum number of individual observations
####you have to delete in order to eliminate the significant result?
jensen.desc<-jensen.cc %>% arrange(desc(weight))
jensen.desc<-jensen.desc[-c(1:15),]
fit.y <- lm(as.formula(paste("Fvar5~regime+", X.vars.f, sep="")), data=jensen.desc) #predict democracy
summary(fit.y)

###generate a new covariate that is 
### uniform random in the range of the treatment variable 
newvar <- sample(0:20,1630, TRUE)

###randomly replace half of the values of this covariate with the true value 
### of the treatment variable
replace.v<-sample(0:1630, 815)
newvar[replace.v]<-jensen.cc$regime[replace.v]
#jensen.cc$indicator<-sample.int(2,1630, replace=T)

##re-run the analysis...how does the figure look?
jensen.cc$newvar<-newvar
fit.y <- lm(as.formula(paste("Fvar5~newvar+", X.vars.f, sep="")), data=jensen.cc) #predict democracy
summary(fit.y)
fit.d <- lm(as.formula(paste("newvar~", X.vars.f, sep="")), data=jensen.cc) #predict regime with covars
summary(fit.d)

d.tilde <- as.numeric(residuals(fit.d)) 
w <- d.tilde^2 

w1 <- tapply(w, jensen.cc$country, mean)
mapnames <-  read.csv("mapnames_filled.csv") 
mapnames$weight <- 0
mapnames$weight[match(rownames(w1), mapnames$jensen)] <- as.numeric(w1) #replace weight with new weights with row names

attributes(world)$data$weight <- 0
attributes(world)$data$weight[match(mapnames$mapname,attributes(world)$data$NAME)] <- mapnames$weight

attributes(world)$data$incl <- 0
attributes(world)$data$incl[match(mapnames$mapname,attributes(world)$data$NAME)] <- as.numeric(!is.na(mapnames$jensen))

# Sample

pdf(file="jensen-nominal-map.pdf", height=5, width=8)
plot(	world, 
      col=gray(1-.75*attributes(world)$data$incl),
      border="gray", lwd=.25)
dev.off()


# Effective sample

pdf(file="jensen-effective-map.pdf", height=5, width=8)
plot(	world, 
      col=gray(1-abs(attributes(world)$data$weight)/max(abs(attributes(world)$data$weight))),
      border="gray", lwd=.25)
dev.off()
