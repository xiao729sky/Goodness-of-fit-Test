# status1, status2: censoring status for two margins, with 0 for right censor and 1 for interval censor 
# return -loglik

step2 <- function(para, u, status1, status2, copula) {
  
  eta<-para
  u1_left = u[,1]
  u1_right = u[,2]
  u2_left = u[,3]
  u2_right = u[,4]
  
  ### copula models ###
  
  if (copula == "Joe") {
    # Joe Copula function for joint distribution probability
    C_val_1 <- 1 - ( (1-u1_left)^(eta) + (1-u2_left)^(eta) - (1-u1_left)^(eta)*(1-u2_left)^(eta) )^(1/eta)
    C_val_2 <- 1 - ( (1-u1_left)^(eta) + (1-u2_right)^(eta) - (1-u1_left)^(eta)*(1-u2_right)^(eta) )^(1/eta)
    C_val_3 <- 1 - ( (1-u1_right)^(eta) + (1-u2_left)^(eta) - (1-u1_right)^(eta)*(1-u2_left)^(eta) )^(1/eta)
    C_val_4 <- 1 - ( (1-u1_right)^(eta) + (1-u2_right)^(eta) - (1-u1_right)^(eta)*(1-u2_right)^(eta) )^(1/eta)
  }
  
  if (copula == "Gumbel") {
    # Gumbel Copula function for joint distribution probability
    C_val_1 <- exp(-((-log(u1_left))^eta + (-log(u2_left))^eta)^(1/eta))
    C_val_2 <- exp(-((-log(u1_left))^eta + (-log(u2_right))^eta)^(1/eta))
    C_val_3 <- exp(-((-log(u1_right))^eta + (-log(u2_left))^eta)^(1/eta))
    C_val_4 <- exp(-((-log(u1_right))^eta + (-log(u2_right))^eta)^(1/eta))
  }
  
  if (copula == "Frank") {
    # Frank Copula function for joint distribution probability
    C_val_1 <- 1/((-1) * eta) * log(1 + (exp((-1) * eta*u1_left)-1)*(exp((-1) * eta*u2_left)-1)/(exp((-1) * eta)-1))
    C_val_2 <- 1/((-1) * eta) * log(1 + (exp((-1) * eta*u1_left)-1)*(exp((-1) * eta*u2_right)-1)/(exp((-1) * eta)-1))
    C_val_3 <- 1/((-1) * eta) * log(1 + (exp((-1) * eta*u1_right)-1)*(exp((-1) * eta*u2_left)-1)/(exp((-1) * eta)-1))
    C_val_4 <- 1/((-1) * eta) * log(1 + (exp((-1) * eta*u1_right)-1)*(exp((-1) * eta*u2_right)-1)/(exp((-1) * eta)-1))
  }
  
  if (copula == "Clayton") {
    # Copula function for joint distribution probability
    C_val_1 <-(u1_left^(-eta)+u2_left^(-eta)-1)^(-1/eta)
    C_val_2 <-(u1_left^(-eta)+u2_right^(-eta)-1)^(-1/eta)
    C_val_3 <-(u1_right^(-eta)+u2_left^(-eta)-1)^(-1/eta)
    C_val_4 <-(u1_right^(-eta)+u2_right^(-eta)-1)^(-1/eta)
  }
  
  if (copula == "AMH") {
    # AMH Copula function for joint distribution probability
    C_val_1 <- u1_left * u2_left /(1 - eta * (1-u1_left) * (1-u2_left))
    C_val_2 <- u1_left * u2_right /(1 - eta * (1-u1_left) * (1-u2_right))
    C_val_3 <- u1_right * u2_left /(1 - eta * (1-u1_right) * (1-u2_left))
    C_val_4 <- u1_right * u2_right /(1 - eta * (1-u1_right) * (1-u2_right))
  }
  
  # Use Copula functions to write each block of likelihood function
  term1 <- log(abs(C_val_1 - C_val_2 - C_val_3 + C_val_4))
  term1 <- ifelse((status1 == 1) & (status2 == 1), term1, 0)
  term1 <- ifelse(term1!=-Inf, term1, mean(term1[which(term1!=-Inf)]))
  
  term2 <- log(abs(C_val_1 - C_val_3))
  term2 <- ifelse((status1 == 1) & (status2 == 0), term2, 0)
  term2 <- ifelse(term2!=-Inf, term2, mean(term2[which(term2!=-Inf)]))
  
  term3 <- log(abs(C_val_1 - C_val_2))
  term3 <- ifelse((status1 == 0) & (status2 == 1), term3, 0)
  term3 <- ifelse(term3!=-Inf, term3, mean(term3[which(term3!=-Inf)]))
  
  term4 <- log(abs(C_val_1))
  term4 <- ifelse((status1 == 0) & (status2 == 0), term4, 0)
  term4 <- ifelse(term4!=-Inf, term4, mean(term4[which(term4!=-Inf)]))
  
  logL<-(-1)*sum( term1 + term2 + term3 + term4 )
  
  return(logL)
}
  