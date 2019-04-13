# Example McCandless(2008) "The ABCs of RBCs" Chap.6 Hansen's RBC----------------------

# Hansen's RBC modelを解いて，Impulse Responseを求めるコード．

#注意1. geigenパッケージをインストールしていなければ、インターネットに繋がっている状況で下の行を実行する
#install.packages("geigen")
#注意2. LDSGE_solver.Rを作業ファイルに入れておくこと．


# geigenパッケージの準備 ----------------------------------------------------------


#call solver
source("LDSGE_solver.R")



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

n_k <- 1        # predetermined variables

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

(sol <- slvLDSGE(A=A,B=B,C=C,Phi=Phi,n_k=n_k))


# impulse response --------------------------------------------------------

t = 102

#response vector
z <- matrix(0, nrow = sol$n_z, ncol = t)
k <- matrix(0, nrow = n_k, ncol = t)
d <- matrix(0, nrow = sol$n - n_k, ncol = t)

#shock
z[1,2] <- matrix(.01)

for (i in 2:t) {
  d[ , i] <- sol$J %*% k[ , i] + sol$N %*% z[ , i]
  if (i < t) {
  k[ , i + 1] <- sol$K %*% k[ , i] + sol$L %*% z[ , i] + xi
  z[ , i + 1] <- Phi %*% z[ , i]
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
