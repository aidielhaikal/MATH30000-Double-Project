# ### DATA ###
# # need to be commented otherwise overwrite data
# trainn = 60; validationn = 60; testn = 200
# corr = 0.9; sigma = 3
# beta = c(rep(0.2,20))
# n  = trainn + validationn + testn
# p  = length(beta)
# # compound symmetry
# cov.mat <- outer(1:p, 1:p, function(x,y){0.2^(x!=y)})
# cov.mat1 <- outer(1:p, 1:p, function(x,y){0.9^(x!=y)})
# X <- mvrnorm(n, rep(0,p), cov.mat) # nxp design matrix
# y <- X%*%beta + rnorm(n,mean = 0,sd=sigma)
# 
# ### CORRELATION PLOT ###
# cowplot::plot_grid(
#   ggcorrplot(cov.mat), ggcorrplot(cov.mat1),
#   ncol = 2, labels = c("(a)","(b)"))

### FITTING FUNCTION ###
simfun <-  function(trainn = 60, validationn = 60, testn = 200, sigma = 3, 
                    corr = 0.2, beta = c(rep(0.2,20)) ){
  n  = trainn + validationn + testn
  p  = length(beta)
  # compound symmetry
  cov.mat <- outer(1:p, 1:p, function(x,y){corr^(x!=y)})
  X <- mvrnorm(n, rep(0,p), cov.mat) 
  colnames(X)<- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8",
                  "x9", "x10", "x11", "x12", "x13", "x14", "x15", "x16",
                  "x17", "x18", "x19", "x20")
  y <- X%*%beta + rnorm(n,mean = 0,sd=sigma)

  ## SPLIT DATA
  trainx <- X[1:trainn,]
  step.trainx <- data.frame(X[1:trainn,])
  validationx <- X[(trainn+1):(trainn+validationn),]
  testx <- X[(trainn+validationn+1):n,]
  trainy <- y[1:trainn,]
  validationy <- y[(trainn+1):(trainn+validationn),]
  testy <- y[(trainn+validationn+1):n,]
  
  ### LEAST SQUARES ###
  fit.ols <- lm(trainy ~ ., data = step.trainx)
  betas.ols <- fit.ols$coefficients[-1]
  # test: compute the test error 
  test.mse.ols <- mean((testy - testx%*%betas.ols)^2)
  test.mse.ols 
  
  ### STEPWISE REGRESSION ###
  fit.stepwise <- step(fit.ols, direction = "both", trace = FALSE)
  betas.stepwise <- fit.stepwise$coefficients[-1]
  
  bs=betas.stepwise
  bs1=numeric(20) # final beta
  nbs1=c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8",
         "x9", "x10", "x11", "x12", "x13", "x14", "x15", "x16",
         "x17", "x18", "x19", "x20")
  names(bs1)=nbs1
  for (i in 1:length(bs))
  {
    for (j in 1:length(nbs1))
    {
      if(nbs1[j]== names(bs[i])) bs1[j]=bs[i]
    }
  }
  # test: compute the test error 
  test.mse.step <- mean((testy - testx %*% bs1)^2)
  test.mse.step
  
  ### RIDGE REGRESSION ###
  # training: find betas for all the lambdas
  fit.ridge <- glmnet(trainx,trainy,alpha=0, intercept = FALSE)
  betas.ridge <- fit.ridge$beta
  # validation: compute validation test error for each lambda and find the minimum
  Val.err.ridge <- colMeans((validationy-validationx%*%betas.ridge)^2)
  beta.opt.ridge <- betas.ridge[, which.min(Val.err.ridge)]
  # test: compute the test error 
  test.mse.ridge <- mean((testy - testx%*%beta.opt.ridge)^2)
  test.mse.ridge
  
  ### LASSO VIA GLMNET ###
  # training: find betas for all the lambdas
  fit.lasso <- glmnet(trainx,trainy,alpha=1, intercept = FALSE)
  betas.lasso <- fit.lasso$beta
  # validation: compute validation test error for each lambda and find the minimum
  Val.err.lasso <- colMeans((validationy - validationx%*%betas.lasso)^2)
  beta.opt.lasso <- betas.lasso[, which.min(Val.err.lasso)]
  # test: compute the test error 
  test.mse.lasso <- mean((testy - testx%*%beta.opt.lasso)^2)
  test.mse.lasso
  
  ### LASSO VIA LARS ###
  fit.lars <- lars(trainx, trainy, type = "lasso", intercept = FALSE)
  betas.lars <- t(fit.lars$beta)
  # validation: compute validation test error for each lambda and find the minimum
  Val.err.lars<- colMeans((validationy-validationx%*%betas.lars)^2)
  beta.opt.lars <- betas.lars[, which.min(Val.err.lars)]
  # test: compute the test error 
  test.mse.lars <- mean((testy - testx%*%beta.opt.lars)^2)
  test.mse.lars
  
  ### ELASTIC NET ###
  # training: find betas and alpha for all the lambdas
  alpha_to_try = seq(0,1, length.out = 11)
  alpha.test.mse = c()
  beta.mat = matrix(nrow = 11, ncol = p)
  for (i in 1:11) {
    fit.en <- glmnet(trainx,trainy,alpha=alpha_to_try[i], intercept = FALSE)
    betas.en <- fit.en$beta
    # validation: compute validation test error for each lambda and alpha
    # and find the minimum
    Val.err.en <- colMeans((validationy-validationx%*%betas.en)^2)
    beta.try.en <- betas.en[,which.min(Val.err.en)]
    # test: compute the test error 
    test.mse.en <- mean((testy - testx%*%beta.try.en)^2)
    alpha.test.mse[i] = test.mse.en
    beta.mat[i,] = beta.try.en
  }
  test.mse.en <- min(alpha.test.mse)
  alpha.opt <- alpha_to_try[which.min(alpha.test.mse)]
  beta.opt.en <- beta.mat[which.min(alpha.test.mse), ]
  
  ### SCAD ###
  # training: find betas for all the lambdas
  fit.scad <- ncpen(trainy,trainx, x.standardize = TRUE, intercept = FALSE)
  betas.scad <- fit.scad$beta
  # validation: compute validation test error for each lambda and find the minimum
  Val.err.scad <- colMeans((validationy- validationx %*% betas.scad)^2)
  beta.opt.scad <- betas.scad[, which.min(Val.err.scad)]
  # test: compute the test error 
  test.mse.scad <- mean((testy - testx%*%beta.opt.scad)^2)
  test.mse.scad
  
  ### GROUP LASSO ###
  # training: find betas for all the lambdas
  fit.glasso <- grpreg(trainx,trainy, penalty = "grLasso")
  betas.glasso <- fit.glasso$beta[-1,]
  # validation: compute validation test error for each lambda and find the minimum
  Val.err.glasso <- colMeans((validationy- validationx%*%betas.glasso)^2)
  beta.opt.glasso <- betas.glasso[, which.min(Val.err.glasso)]
  # test: compute the test error 
  test.mse.glasso <- mean((testy - testx%*%beta.opt.glasso)^2)
  test.mse.glasso
  
  list(trainx = trainx, beta = beta,
       betas.ols = betas.ols, test.mse.ols = test.mse.ols,
       betas.step = bs1, test.mse.step = test.mse.step,
       beta.opt.ridge = beta.opt.ridge, test.mse.ridge = test.mse.ridge,
       beta.opt.lasso = beta.opt.lasso, test.mse.lasso = test.mse.lasso,
       beta.opt.lars = beta.opt.lars, test.mse.lars = test.mse.lars,
       alpha.opt = alpha.opt, beta.opt.en = beta.opt.en, test.mse.en = test.mse.en,
       beta.opt.scad = beta.opt.scad, test.mse.scad = test.mse.scad,
       beta.opt.glasso = beta.opt.glasso, test.mse.glasso = test.mse.glasso
  )
}


### MEASURES OF ACCURACY & ASSESSING VARIABLE SELECTION ###
true.coef = c()
false.coef = c(1:20)

ols.test.mse = c()
ols.nonzero.coef = c()
ols.true.zero = c()
ols.false.zero = c()
ols.false.nonzero = c()

step.test.mse = c()
step.nonzero.coef = c()
step.true.zero = c()
step.false.zero = c()
step.false.nonzero = c()

ridge.test.mse = c()
ridge.nonzero.coef = c()
ridge.true.zero = c()
ridge.false.zero = c()
ridge.false.nonzero = c()

lasso.test.mse = c()
lasso.nonzero.coef = c()
lasso.true.zero = c()
lasso.false.zero = c()
lasso.false.nonzero = c()

lars.test.mse = c() 
lars.nonzero.coef = c()
lars.true.zero = c()
lars.false.zero = c()
lars.false.nonzero = c()

en.test.mse = c()
en.nonzero.coef = c()
en.alpha.opt = c()
en.true.zero = c()
en.false.zero = c()
en.false.nonzero = c()

scad.test.mse = c()
scad.nonzero.coef = c()
scad.true.zero = c()
scad.false.zero = c()
scad.false.nonzero = c()

glasso.test.mse = c()
glasso.nonzero.coef = c()
glasso.true.zero = c()
glasso.false.zero = c()
glasso.false.nonzero = c()

snr = c()

for (i in 1:100) {
  data <- simfun(trainn = 60, validationn = 60, testn = 200, sigma = 3, 
                 corr = 0.2, beta = c(rep(0.2,20)) );
  snr[i] = sqrt(sum((data$trainx %*% data$beta)**2))/3**2
  
  ols.test.mse[i] = data$test.mse.ols
  ols.nonzero.coef[i] = sum(data$betas.ols != 0)
  step.test.mse[i] = data$test.mse.step
  step.nonzero.coef[i] = sum(data$betas.step != 0)
  ridge.test.mse[i] = data$test.mse.ridge
  ridge.nonzero.coef[i] = sum(data$beta.opt.ridge != 0)
  lasso.test.mse[i] = data$test.mse.lasso
  lasso.nonzero.coef[i] = sum(data$beta.opt.lasso != 0)
  lars.test.mse[i] = data$test.mse.lars
  lars.nonzero.coef[i] = sum(data$beta.opt.lars != 0)
  en.test.mse[i] = data$test.mse.en
  en.nonzero.coef[i] = sum(data$beta.opt.en != 0)
  en.alpha.opt[i] = data$alpha.opt
  scad.test.mse[i] = data$test.mse.scad
  scad.nonzero.coef[i] = sum(data$beta.opt.scad != 0)
  glasso.test.mse[i] = data$test.mse.glasso
  glasso.nonzero.coef[i] = sum(data$beta.opt.glasso != 0)
  
  ols.sum.false = c()
  for (j in false.coef) {ols.sum.false = append(ols.sum.false,(data$betas.ols[j] ==0))}
  ols.false.zero[i] = sum(ols.sum.false)
  ols.sum.true = c()
  for (j in true.coef) {ols.sum.true= append(ols.sum.true,(data$betas.ols[j] ==0))}
  ols.true.zero[i] = sum(ols.sum.true)
  ols.sum.false.nonzero = c()
  for (j in true.coef) {ols.sum.false.nonzero = append(ols.sum.false.nonzero, (data$betas.ols[j] != 0))}
  ols.false.nonzero[i] = sum(ols.sum.false.nonzero)
  
  step.sum.false = c()
  for (j in false.coef) {step.sum.false = append(step.sum.false,(data$betas.step[j] ==0))}
  step.false.zero[i] = sum(step.sum.false)
  step.sum.true = c()
  for (j in true.coef) {step.sum.true= append(step.sum.true,(data$betas.step[j] ==0))}
  step.true.zero[i] = sum(step.sum.true)
  step.sum.false.nonzero = c()
  for (j in true.coef) {step.sum.false.nonzero = append(step.sum.false.nonzero, (data$betas.step[j] != 0))}
  step.false.nonzero[i] = sum(step.sum.false.nonzero)
  
  ridge.sum.false = c()
  for (j in false.coef) {ridge.sum.false = append(ridge.sum.false,(data$beta.opt.ridge[j] ==0))}
  ridge.false.zero[i] = sum(ridge.sum.false)
  ridge.sum.true = c()
  for (j in true.coef) {ridge.sum.true= append(ridge.sum.true,(data$beta.opt.ridge[j] ==0))}
  ridge.true.zero[i] = sum(ridge.sum.true)
  ridge.sum.false.nonzero = c()
  for (j in true.coef) {ridge.sum.false.nonzero = append(ridge.sum.false.nonzero, (data$beta.opt.ridge[j] != 0))}
  ridge.false.nonzero[i] = sum(ridge.sum.false.nonzero)
  
  lasso.sum.false = c()
  for (j in false.coef) {lasso.sum.false = append(lasso.sum.false,(data$beta.opt.lasso[j] ==0))}
  lasso.false.zero[i] = sum(lasso.sum.false)
  lasso.sum.true = c()
  for (j in true.coef) {lasso.sum.true= append(lasso.sum.true,(data$beta.opt.lasso[j] ==0))}
  lasso.true.zero[i] = sum(lasso.sum.true)
  lasso.sum.false.nonzero = c()
  for (j in true.coef) {lasso.sum.false.nonzero = append(lasso.sum.false.nonzero, (data$beta.opt.lasso[j] != 0))}
  lasso.false.nonzero[i] = sum(lasso.sum.false.nonzero)
  
  lars.sum.false = c()
  for (j in false.coef) {lars.sum.false = append(lars.sum.false,(data$beta.opt.lars[j] ==0))}
  lars.false.zero[i] = sum(lars.sum.false)
  lars.sum.true = c()
  for (j in true.coef) {lars.sum.true= append(lars.sum.true,(data$beta.opt.lars[j] ==0))}
  lars.true.zero[i] = sum(lars.sum.true)
  lars.sum.false.nonzero = c()
  for (j in true.coef) {lars.sum.false.nonzero = append(lars.sum.false.nonzero, (data$beta.opt.lars[j] != 0))}
  lars.false.nonzero[i] = sum(lars.sum.false.nonzero)
  
  en.sum.false = c()
  for (j in false.coef) {en.sum.false = append(en.sum.false,(data$beta.opt.en[j] ==0))}
  en.false.zero[i] = sum(en.sum.false)
  en.sum.true = c()
  for (j in true.coef) {en.sum.true= append(en.sum.true,(data$beta.opt.en[j] ==0))}
  en.true.zero[i] = sum(en.sum.true)
  en.sum.false.nonzero = c()
  for (j in true.coef) {en.sum.false.nonzero = append(en.sum.false.nonzero, (data$beta.opt.en[j] != 0))}
  en.false.nonzero[i] = sum(en.sum.false.nonzero)

  scad.sum.false = c()
  for (j in false.coef) {scad.sum.false = append(scad.sum.false,(data$beta.opt.scad[j] ==0))}
  scad.false.zero[i] = sum(scad.sum.false)
  scad.sum.true = c()
  for (j in true.coef) {scad.sum.true= append(scad.sum.true,(data$beta.opt.scad[j] ==0))}
  scad.true.zero[i] = sum(scad.sum.true)
  scad.sum.false.nonzero = c()
  for (j in true.coef) {scad.sum.false.nonzero = append(scad.sum.false.nonzero, (data$beta.opt.scad[j] != 0))}
  scad.false.nonzero[i] = sum(scad.sum.false.nonzero)
  
  glasso.sum.false = c()
  for (j in false.coef) {glasso.sum.false = append(glasso.sum.false,(data$beta.opt.glasso[j] ==0))}
  glasso.false.zero[i] = sum(glasso.sum.false)
  glasso.sum.true = c()
  for (j in true.coef) {glasso.sum.true= append(glasso.sum.true,(data$beta.opt.glasso[j] ==0))}
  glasso.true.zero[i] = sum(glasso.sum.true)
  glasso.sum.false.nonzero = c() 
  for (j in true.coef) {glasso.sum.false.nonzero = append(glasso.sum.false.nonzero, (data$beta.opt.glasso[j] != 0))}
  glasso.false.nonzero[i] = sum(glasso.sum.false.nonzero)
}

### OUTPUT TABLE ###
output = matrix(c(median(ols.test.mse), median(step.test.mse), median(ridge.test.mse), median(lasso.test.mse), median(lars.test.mse), median(en.test.mse), median(scad.test.mse), median(glasso.test.mse),
                  mean(ols.true.zero), mean(step.true.zero), mean(ridge.true.zero), mean(lasso.true.zero), mean(lars.true.zero), mean(en.true.zero), mean(scad.true.zero), mean(glasso.true.zero),
                  mean(ols.false.zero), mean(step.false.zero), mean(ridge.false.zero), mean(lasso.false.zero), mean(lars.false.zero), mean(en.false.zero), mean(scad.false.zero), mean(glasso.false.zero),
                  median(ols.nonzero.coef), median(step.nonzero.coef), median(ridge.nonzero.coef), median(lasso.nonzero.coef), median(lars.nonzero.coef), median(en.nonzero.coef), median(scad.nonzero.coef), median(glasso.nonzero.coef),
                  mean(ols.false.nonzero), mean(step.false.nonzero), mean(ridge.false.nonzero), mean(lasso.false.nonzero), mean(lars.false.nonzero), mean(en.false.nonzero), mean(scad.false.nonzero), mean(glasso.false.nonzero)),
                  ncol = 5, nrow = 8, byrow = FALSE)

rownames(output) <-  c("OLS", "Stepwise", "Ridge", "LASSO", "LASSO(lars)", "Elastic Net" , "SCAD", "group LASSO")
colnames(output) <-  c("median.MSE", "true.zero", "false.zero", "median.nonzero", "false.nonzero")

output = as.table(output); output
mean(snr); mean(en.alpha.opt)
