
# LDSGE_solver.R ----------------------------------------------------------
# 線型の確率的動学的一般均衡モデル
# B E_t(x_{t+1}) = A x_t + C z_t
# x_t = [k_t d_t]'
# z_t = Phi z_{t-1} + epsilon_t
# xi_{t+1} = k_{t+1} - E_t(k_{t+1}) is exogenous
# E(epsilon_t) = 0
# E_t(xi_{t+1}) = 0
#
# x_t :n*1
# k_t :n_k*1
# d_t :(n-n_k)*1
# z_t :n_z*1
# A   :n*n
# B   :n*n
# C   :n*n_z
#
# を解いて，合理的期待均衡解
# k_{t+1} = K k_t + L z_t + xi_{t+1}
# d_t     = J k_t + N z_t
# の行列K, L, J, Nを計算する関数slvLDSGE(A,B,C,Phi,n_k)


# geigenパッケージの準備 ----------------------------------------------------------

# geigenをインストールしていなければ、インターネットに繋がっている状況で下の行を実行する
#install.packages("geigen")

# geigenを読み出す
library(geigen)

slvLDSGE <- function(A,B,C,Phi,n_k) {
  
  n <- nrow(A)
  n_z <- ncol(C)
  stopifnot(n==ncol(A),dim(A)==dim(B),n==nrow(C))
  
  #complexify
  A <- A+0i
  B <- B+0i
  
  #QZ decomposition
  QZ <- gqz(A, B, "S")
  
  #generalized eigenvalues
  ger <- gevalues(QZ)
  
  #test stability
  n_s <- sum(abs(ger)<1) # the number of stable eigenvalues
  if (n_s > n_k) {
    stop(message="Indeterminate!")
  } else if (n_s < n_k) {
    stop(message="no stable solution!")
  } else {
    message(paste("saddle path stable!"))
  }
  
  n_u <- n - n_s            # unstable eigenvalues
  
  #行列Mの計算
  G <- t(Conj(QZ$Q))[c((n_s+1):n), ] %*% C   # n*n_z
  S_22 <- QZ$S[c((n_s+1):n),c((n_s+1):n)]    # n_u*n_u
  T_22 <- QZ$T[c((n_s+1):n),c((n_s+1):n)]    # n_u*n_u
  
  M <- matrix(0, nrow=n_u, ncol=n_z)        # n_u*n_z
  R <- matrix(0, nrow=n_u, ncol=n_z)
  
  for (i in n_u:1) {
    if (i == n_u) {
      R[i, ] <- G[n_u, ] # 1*n_z
    } else {
      R[i, ] <- G[i, ] + S_22[i,c((i+1):n_u)] %*% M[c((i+1):n_u), ] - T_22[i,c((i+1):n_u)] %*% M[c((i+1):n_u), ] %*% Phi  # 1*n_z
    }
    M[i, ] <- R[i, ] %*% solve(T_22[i,i] %*% Phi - S_22[i,i] %*% diag(1, nrow = n_z, ncol = n_z)) # 1*n_z
  }
  
  
  # 行列J,N,K,Lを計算 ------------------------------------------------------------
  
  Z_11 <- QZ$Z[c(1:n_s),c(1:n_s)]
  Z_12 <- QZ$Z[c(1:n_s),c((n_s+1):n)]
  Z_21 <- QZ$Z[c((n_s+1):n),c(1:n_s)]
  Z_22 <- QZ$Z[c((n_s+1):n),c((n_s+1):n)]
  
  S_11 <- QZ$S[c(1:n_s),c(1:n_s)]
  S_12 <- QZ$S[c(1:n_s),c((n_s+1):n)]
  T_11 <- QZ$T[c(1:n_s),c(1:n_s)]
  T_12 <- QZ$T[c(1:n_s),c((n_s+1):n)]
  
  #realify
  J <- as.matrix(as.numeric(Z_21 %*% solve(Z_11)))
  N <- as.matrix(as.numeric((Z_22 - Z_21 %*% solve(Z_11) %*% Z_12) %*% M))
  K <- as.matrix(as.numeric(Z_11 %*% solve(T_11) %*% S_11 %*% solve(Z_11)))
  L <- as.matrix(as.numeric(-Z_11 %*% solve(T_11) %*% S_11 %*% solve(Z_11) %*% Z_12 %*% M + Z_11 %*% solve(T_11) %*% (S_12 %*% M - T_12 %*% M %*% Phi + t(Conj(QZ$Q))[c(1:n_s), ] %*% C) + Z_12 %*% M %*% Phi))
  
  solution <- list(n=n,n_z=n_z,ger=ger,J=J,N=N,K=K,L=L)
  return(solution)
}

#end
