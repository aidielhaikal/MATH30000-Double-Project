rm(list= ls())
###  IMPORT LIBRARIES ### 
# library for fitting different methods
library(MASS)   # stepwise regression
library(glmnet) # LASSO, ridge, and elastic net
library(ncpen)  # SCAD
library(ncvreg) # SCAD
library(grpreg) # group LASSO
library(lars)   # LASSO using lars algorithm 
# library for plotting
library(ggcorrplot)
library(GGally)
library(plotmo)
library(psych)
# library for diagnostic checking 
library(car)


###  FUNCTIONS ###
# function for plotting legend
col.list <- c("tan1", "steelblue2", "yellow2", "violetred2", "turquoise3", 
              "slateblue2", "thistle2", "red2", "gray8", "gray49")
leg_fun <- function(fit, ...) {
  L <- length(fit$lambda)
  x <- log(fit$lambda[L])
  y <- fit$beta[, L]
  labs <- names(y)[-1]
  legend("topright",  
         inset = c(-0.3,-0.1),
         legend=labs, 
         col=col.list, 
         lty=1,
         bty = "n",
         xpd = TRUE)
}


###  DATA EXPLORATION ### 
diabetes <- read.csv("Diabetes.csv", header = TRUE)
dim(diabetes)
diabetes[,-11] <- scale(diabetes[,-11])
diab.mat <- as.matrix(diabetes)
str(diabetes)
summary(diabetes)
ggcorrplot(cor(diabetes), lab = TRUE)


### FORMAT DATA FOR FITTING ### 
set.seed(7654321)
trainn = 332; testn  = 110 
n <- trainn + testn
p <- length(names(diabetes[,-11]))
indx <- sample(1:442,size=trainn, replace = FALSE)
# in vector
trainx <- diab.mat[indx,-11]; trainy <- diab.mat[indx,11]
testx <- diab.mat[-indx,-11]; testy <- diab.mat[-indx,11]
# in data frame
dat1 <- data.frame(trainx, stringsAsFactors = FALSE); dat1$y <- trainy
dat2 <- data.frame(testx, stringsAsFactors = FALSE); dat2$y <- testy


### ORDINARY LEAST SQUARES ###
fit.ols <- lm(y ~ ., data = dat1)
summary(fit.ols) 
betas.ols <- fit.ols$coefficients[-1]
# predict
pred.ols <- predict(fit.ols, newdata = dat2[,-11], type = "response")
# test
test.mse.ols <- mean((testy - pred.ols)^2); test.mse.ols # 3218.804
# testing collinearity
vif(fit.ols)


### STEPWISE REGRESSION USING AIC ###
fit.stepwise <- step(fit.ols, direction = "both")
betas.stepwise <- fit.stepwise$coefficients[-1]
bs=betas.stepwise
bs1=numeric(10) 
nbs1=c("age" ,"sex" ,"bmi" ,"map" ,"tc" , "ldl", "hdl" ,"tch", "ltg", "glu")
names(bs1)=nbs1
for (i in 1:length(bs)){
  for (j in 1:length(nbs1)){
    if(nbs1[j]== names(bs[i])) bs1[j]=bs[i]}}
# predict
pred.step <- predict(fit.stepwise, newdata = dat2[,-11])
# test 
test.mse.step <- mean((testy - pred.step)^2); test.mse.step # 3180.197


### RIDGE REGRESSION USING 10-FOLD CV ###
fit.ridge <- glmnet(trainx, trainy, alpha=0, intercept = FALSE, standardize = TRUE)
# cross-validation
set.seed(7654321)
cv.ridge <- cv.glmnet(trainx,trainy,alpha=0,nfolds=10)
cv.ridge$lambda.min; cv.ridge$lambda.1se         # 4.86236; 54.62005
betas.opt.ridge  <- coef(fit.ridge, s = cv.ridge$lambda.min)[-1,]
betas.opt.ridge1 <- coef(fit.ridge, s = cv.ridge$lambda.1se)[-1,]
cbind(betas.opt.ridge, betas.opt.ridge1)
# predict
pred.ridge <- predict(cv.ridge, newx = testx, s = "lambda.min")
pred.ridge1 <- predict(cv.ridge, newx = testx, s = "lambda.1se")
# test
test.mse.ridge <- mean((testy - pred.ridge)^2)   # 3140.878
test.mse.ridge1 <- mean((testy - pred.ridge1)^2) # 3222.77
test.mse.ridge; test.mse.ridge1
# plot
par(mar = c(5, 5, 2, 6))
plot(fit.ridge, xvar="lambda", col=col.list, 
     ylab = "Standardized Coefficients")
leg_fun(fit.ridge)
abline(v = log(cv.ridge$lambda.min), lty = 5)
abline(v = log(cv.ridge$lambda.1se), lty = 3)
plot(cv.ridge)


### LASSO USING 10-FOLD CV VIA GLMNET ###
fit.lasso <- glmnet(trainx,trainy, alpha=1, intercept = FALSE)
# cross-validation
set.seed(7654321)
cv.lasso = cv.glmnet(trainx,trainy,alpha = 1,nfolds = 10)
cv.lasso$lambda.min; cv.lasso$lambda.1se # 0.2009133; 8.301762
betas.opt.lasso <- coef(fit.lasso, s = cv.lasso$lambda.min)[-1,]
betas.opt.lasso1 <- coef(fit.lasso, s = cv.lasso$lambda.1se)[-1,]
cbind(betas.opt.lasso, betas.opt.lasso1)
# predict
pred.lasso  <- predict(cv.lasso, newx = testx, s = "lambda.min")
pred.lasso1 <- predict(cv.lasso, newx = testx, s = "lambda.1se")
# test
test.mse.lasso <- mean((testy - pred.lasso)^2)    # 3164.053
test.mse.lasso1 <- mean((testy - pred.lasso1)^2)  # 3246.901
test.mse.lasso; test.mse.lasso1
# plot
par(mar = c(5, 5, 3, 6))
plot(fit.lasso, xvar="lambda", 
     col=col.list, ylab = "Standardized Coefficients")
leg_fun(fit.lasso)
abline(v = log(cv.lasso$lambda.min), lty = 5)
abline(v = log(cv.lasso$lambda.1se), lty = 3)
plot(cv.lasso)


### ELASTIC NET USING 10-FOLD CV ###
# optimize alpha
alpha <- seq(from = 0, to = 1, by = 0.1)
min.cv.en <- c()
for (i in 1:11) {
  set.seed(7654321)
  cv.en <- cv.glmnet(x = trainx, y = trainy, alpha = alpha[i], nfolds =10)
  min.cv.en[i] <- min(cv.en$cvm)
}
alpha.opt <- alpha[which.min(min.cv.en)]; alpha.opt # 0.1

fit.en <- glmnet(trainx,trainy, alpha=alpha.opt, intercept = FALSE)
# cross-validation
set.seed(7654321)
cv.en = cv.glmnet(x = trainx, y = trainy, alpha = alpha.opt, nfolds =10)
cv.en$lambda.min; cv.en$lambda.1se # 0.5462005; 32.74386
betas.opt.en  <- coef(fit.en, s = cv.en$lambda.min)[-1,]
betas.opt.en1 <- coef(fit.en, s = cv.en$lambda.1se)[-1,]
cbind(betas.opt.en, betas.opt.en1)
# predict
pred.en  <- predict(cv.en, newx = testx, s = "lambda.min")
pred.en1 <- predict(cv.en, newx = testx, s = "lambda.1se")
# test
test.mse.en  <- mean((testy - pred.en)^2)
test.mse.en1 <- mean((testy - pred.en1)^2)
test.mse.en; test.mse.en1 # 3175.237; 3244.327
# plot
plot(alpha, min.cv.en, 
     type = "l", 
     xlab = "alpha", 
     ylab = "Mean Cross-Validated Error")
points(alpha, y = min.cv.en, pch = 19, col ="red")


par(mar = c(5, 5, 4, 6))
plot(fit.en, xvar="lambda", 
     col=col.list, ylab = "Standardized Coefficients")
leg_fun(fit.en)
abline(v = log(cv.en$lambda.min), lty = 5)
abline(v = log(cv.en$lambda.1se), lty = 3)
plot(cv.en)


### SCAD USING 10-FOLD CV ###
fit.scad <- ncvreg(X = trainx, y = trainy, family = "gaussian", penalty = "SCAD")
# cross-validation
set.seed(7654321)
cv.scad <- cv.ncvreg(X = trainx, y = trainy, nfolds = 10)
betas.opt.scad <- coef(cv.scad, s = "lambda.min")[-1]; betas.opt.scad
s.scad <- log(cv.scad$lambda.min); s.scad # -0.7210613
# predict
pred.scad <- predict(fit.scad, X = testx, lambda = cv.scad$lambda.min)
# test
test.mse.scad <- mean((testy - pred.scad)^2); test.mse.scad # 3216.764
# plot
plot(cv.scad, log.l = TRUE, type = "cve", vertical.line = TRUE,
     ylab = "Mean-Square Error",
     xlim = c(-3, 4))

par(mar = c(5, 5, 3, 6))
plot(fit.scad, log.l = TRUE, lty = 1, col = col.list, 
     xlab = "Log Lambda", ylab =  "Standardized Coefficients",
     xlim = c(-3, 4))
abline(v = s.scad, lty = 5)
leg_fun(fit.scad)


### GROUP LASSO VIA 10-FOLD CV ###

fit.grlas<- grpreg(X = trainx, y = trainy, penalty = "grLasso")
# cross-validation
cv.grlas <- cv.grpreg(X = trainx, y = trainy, nfolds = 10, seed = 7654321)
betas.opt.grlas <- coef(cv.grlas)[-1]; betas.opt.grlas 
s.grlas <- log(cv.grlas$lambda.min); s.grlas   # -1.604882
# predict
pred.grlas <- predict(cv.grlas, X = testx, type ="response", lambda = cv.grlas$lambda.min)
# test
test.mse.grlas <- mean((testy - pred.grlas)^2); test.mse.grlas  # 3164.412
# plot 
par(mar = c(5, 5, 3, 6))
plot(fit.grlas, log.l = TRUE, col = col.list,
     xlab = "Log Lambda", ylab =  "Standardized Coefficients",
     xlim = c(-6,4))
abline(v = s.grlas, lty = 5)
leg_fun(fit.grlas)

plot(cv.grlas, log.l = TRUE, type = "cve",
     ylab = "Mean-Squared Error", 
     ylim = c(2500,7000), xlim = c(-6,4))

### PLOT: COMPARING TEST MSE ###
sel.meth <- c("Stepwise", "Ridge", "LASSO", "Elastic Net", "SCAD", "GLASSO")
test.err <- c(test.mse.step, test.mse.ridge, test.mse.lasso, 
              test.mse.en, test.mse.scad, test.mse.grlas)
barplot(test.err, space = , names.arg = sel.meth, col= "darkred",
        ylab = "Test Error", ylim = c(3000, 3250), xpd=FALSE)
abline(h = test.mse.ols, lty = 2)

### BETA COEFFICIENTS ###
betas.ols
bs1
betas.opt.ridge
betas.opt.ridge1
betas.opt.lasso
betas.opt.lasso1
betas.opt.en
betas.opt.en1
betas.opt.scad
betas.opt.grlas

### TEST MSE ###
c(test.mse.ols, test.mse.step, test.mse.scad, test.mse.grlas)
c(test.mse.ridge, test.mse.ridge1, test.mse.lasso, test.mse.lasso1, test.mse.en, test.mse.en1)
