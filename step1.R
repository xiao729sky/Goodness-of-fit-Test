## NPMLE survival function

NP_f <- function(NP1, NP2, t_four) { 

  u <- matrix(0, nrow(t_four), 4)    # return NPMLE
  t <- t_four
  
  S1 <- NP1[["scurves"]][["1"]][["Tbull_ints"]][,2]
  S2 <- NP2[["scurves"]][["1"]][["Tbull_ints"]][,2]
  S1_SP <- NP1[["scurves"]][["1"]][["S_curves"]][["baseline"]]
  S2_SP <- NP2[["scurves"]][["1"]][["S_curves"]][["baseline"]]
  
  u[which(t[,1]<=S1[1]),1]<-1
  u[which(t[,3]<=S2[1]),3]<-1
  
  for (i in 1:nrow(t)){
    for (w in 1:2){
      j<-1
      while (j < length(S1)){
        if(t[i,w]>S1[j] & t[i,w]<=S1[j+1]){
          u[i,w] <- S1_SP[j]
        }
        j = j+1
      }
    }
    for (w in 3:4){
      j<-1
      while (j < length(S2)){
        if(t[i,w]>S2[j] & t[i,w]<=S2[j+1]){
          u[i,w] <- S2_SP[j]
        }
        j = j+1
      }
    }
  }
  u[,1] <- ifelse(u[,1]==0, 0.001, u[,1])
  u[,2] <- ifelse(u[,1]==u[,2], u[,2]-0.001, u[,2])
  u[,3] <- ifelse(u[,3]==0, 0.001, u[,3])
  u[,4] <- ifelse(u[,3]==u[,4], u[,4]-0.001, u[,4])
  
  return(u)
}

