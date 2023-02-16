library("survival")
library("MASS")
library("Icens")
library("icensBKL")
library("icenReg")   
library("copula")
library("numDeriv")
## functions for obtaining the estimates of the copula parameter
data(tandmob)

source("model_functions.R")
source("step1.R")
source("step2.R")
source("likhood.R")

options(warn = -1)                  # ignore all warnings
odd_index_f  <- seq(1, nrow(tandmob[1:200,]), by = 1)  # sequence
eta_knife <- numeric()                        
# esimate of copula parameter for the conventional PIOL

copula0 <- "Clayton" # Gumbel, Frank

data_example0 <- tandmob[1:200,]
data_example <- data_example0
n <- nrow(data_example)
##################  data set
status1 <- rep(1, n)
status2 <- rep(1, n)
for (i in 1:n) {
  if (is.na(data_example$R14[i])) {
    status1[i] <- 0
  }
}
for (i in 1:n) {
  if (is.na(data_example$R44[i])) {
    status2[i] <- 0
  }
}

data_example$L14[is.na(data_example$L14)] <- 0
data_example$R14[is.na(data_example$R14)] <-
  max(na.omit(data_example$R14)) + 0.1
data_example$L44[is.na(data_example$L44)] <- 0
data_example$R44[is.na(data_example$R44)] <-
  max(na.omit(data_example$R44)) + 0.1
max_end = max(data_example$R14, data_example$R44)
t_left1 <- data_example$L14
t_right1 <- data_example$R14
t_left2 <- data_example$L44
t_right2 <- data_example$R44
data_example$ind <- rep(1, n)  # For using function 'ic_np'

NPMLE1 <- ic_np(cbind(L14, R14)~ind, data = data_example)
NPMLE2 <- ic_np(cbind(L44, R44)~ind, data = data_example)
t_four <- cbind(t_left1, t_right1, t_left2, t_right2)  # observation times
u0 <- NP_f(NP1=NPMLE1, NP2=NPMLE2, t_four=t_four)      

model_step2 = try(nlm(
  step2,
  p = 0.7,              # initial value of the copula parameter
  hessian = F,
  u = u0,
  status1 = status1,
  status2 = status2,
  copula = copula0
),silent = TRUE)
if ('try-error' %in% class(model_step2)){
  next}
eta_e = model_step2$estimate       # estimate of the copula parameter
tau_e= tau_copula(eta_e, copula0)  # estimate of Kendall's tau
print(c(tau_e, copula0))

### l_in
################
l_in = loglik_sum(para =  eta_e, u=u0, status1 = status1,
                  status2 = status2, copula = copula0)
### l_out
################
hess = hessian(loglik_sum,x=eta_e,u=u0, status1 = status1,
             status2 = status2, copula = copula0)

hess_inv = try(solve(hess),silent=T)
if (!(is.character(hess_inv) & !any(is.na(hess_inv)))){
  jaco = jacobian(loglik,x=eta_e,u=u0, status1 = status1,
                status2 = status2, copula = copula0)
  if (is.matrix(jaco)) V = t(jaco)%*%jaco else V = sum(jaco^2)
  ## IR staistic
  IR = sum(diag(-hess_inv%*%V))
  ## approximated PIOL
  eta_new = eta_e + (hess_inv)%*%t(jaco)
  l_out_approx = lapply(1:ncol(eta_new),function(i){
    loglik_out(para=eta_new[,i],u=u0[i,],status1 = status1[i],
               status2 = status2[i], copula = copula0)
  })
  PIOL_T_approx <- l_in-sum(unlist(l_out_approx))
}

## conventional PIOL
for (i in odd_index_f) {
  model_step2 = try(nlm(
    step2,
    p = 0.7,
    hessian = F,
    u = u0[-i,],
    status1 = status1[-i],
    status2 = status2[-i],
    copula = copula0
  ),silent = TRUE)
  if ('try-error' %in% class(model_step2)){
    next}
  eta = model_step2$estimate
  eta_knife[i] = eta
}

l_out = lapply(1:n,function(i){
  loglik_out(para=eta_knife[i],u=u0[i,],status1 = status1[i],
             status2 = status2[i],copula = copula0)
})

PIOL_T <- l_in-sum(unlist(l_out))
