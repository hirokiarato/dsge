
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
# の行列K, L, J, Nを計算する．


# geigenパッケージの準備 ----------------------------------------------------------

# geigenをインストールしていなければ、インターネットに繋がっている状況で下の行を実行する
#install.packages("geigen")

# geigenを読み出す
library(geigen)



# Example McCandless(2008) "The ABCs of RBCs" Chap.6 ----------------------


# Set parameter values --------------------------------------------------------------

beta = .99    #subjective discount factor
delta = .025  #capital depreciation rate
alpha = .36   #capital share
Abar = 1      #TFP
psi = 1.72    #labor disutility coefficient
phi_a = 0.95  #AR(1) coefficient of a~_t
#TFP shock: a~_t = 0.95 * a~_{t-1} + epsilon_t


# Calculate Steady State --------------------------------------------------

K_H <- ((alpha*Abar)/(1/beta - 1 + delta))^(1/(1-alpha))
C_H <- Abar*K_H^alpha - delta * K_H
H_1minusH <- (1-alpha)*Abar*K_H^alpha/(psi*C_H)
H_ss <- H_1minusH/(1+H_1minusH)
K_ss <- K_H * H_ss
C_ss <- C_H * H_ss
R_ss <- beta^(-1) - 1 + delta
Y_ss <- Abar * K_ss^alpha * H_ss^(1-alpha)


# Count the number of variables --------------------------------------

#x_t = [~k_t, y~_t, c~_t, r~_t]', k_t = [k~_t], d_t = [y~_t, c~_t, r~_t]',z_t = [a~_t].

n <- 4          # endogenous variables
n_k <- 1
n_z <- 1

# Set matrices ------------------------------------------------------------

B <- matrix(c(K_ss, 0, 0, 0,
              0, 0, 0, 0,
              0, 0, 0, 0,
              0, 0, 1, -beta*R_ss),
              nrow = 4, byrow = TRUE)

A <- matrix(c((1-delta)*K_ss, Y_ss, -C_ss, 0,
              alpha, -alpha, -(1-alpha), 0,
              1, -1, 0, 1,
              0, 0, 1, 0),
              nrow = 4, byrow = TRUE)
#complexify
B <- B+0i
A <- A+0i

C <- matrix(c(0,
              1,
              0,
              0),
              nrow = 4, byrow = TRUE)

Phi <- matrix(phi_a)

xi <- matrix(0)


# QZ decomposition
QZ <- gqz(A, B, "S")

ger <- gevalues(QZ)
ger

#安定性の確認
n_s <- sum(abs(ger) < 1)  # stable eigenvalues
n_s == n_k                # Test stability
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


# impulse response --------------------------------------------------------

t = 102

#response vector
z <- matrix(0, nrow = n_z, ncol = t)
k <- matrix(0, nrow = n_k, ncol = t)
d <- matrix(0, nrow = n - n_k, ncol = t)

#shock
z[1,2] <- matrix(.01)

for (i in 2:t) {
  d[ , i] <- J %*% k[ , i] + N %*% z[ , i]
  if (i < t) {
  k[ , i + 1] <- K * k[ , i] + L %*% z[ , i] + xi
  z[ , i + 1] <- Phi * z[ , i]
  }
}

Periods <- -1:(t-2)
labels <- c("a","k","y","c","r")
cols <- c("black","purple","green","red","blue")
ltys <- c(1,1,1,1)
lwds <- c(2,2,2,2)


plot(Periods,d[1,],type = "l", col = "green", ylim = c(-.005,.02),ylab = "",lwd=2) #y
par(new=T)
plot(Periods,d[2,],type = "l", col = "red", ylim = c(-.005,.02),ann=F,axes=F,ylab = "",lwd=2) #c
par(new=T)
plot(Periods,d[3,],type = "l", col = "blue", ylim = c(-.005,.02),ann=F,axes=F,ylab = "",lwd=2) #r
par(new=T)
plot(Periods,k[1,],type = "l", col = "purple", ylim = c(-.005,.02),ann=F,axes=F,ylab = "",lwd=2) #k
par(new=T)
plot(Periods,z[1,],type = "l", col = "black", ylim = c(-.005,.02),ann=F,axes=F,ylab = "",lwd=2) #a
abline(h=0,lty=1)
abline(v=0,lty=2)
legend("topright", legend = labels, col = cols, lty = ltys, lwd = lwds)

#end
