# The expriments here are based on four functions in
# Xinkai Zhou  287D final_project.R

set.seed(24500418)
########################################################
# Experiment 0: compare differnt stochastic methods
# Generate data
X <- rbind(matrix(rnorm(120), , 2), matrix(rnorm(120, mean = 3), , 2))
y <- matrix(c(rep(1, 60), rep(-1, 60)))
epsilon=10^-2
loss=2
gamma=0
lambda=10^-2

# Sto=1, Permutation
stochastic=1
expe_0_sto_1=sdca(X,y,lambda,epsilon,stochastic,loss)

# Save plot
pdf("expe_0_sto_1.pdf")
plot(expe_0_sto_1[4,]~expe_0_sto_1[1,],log="y",xlab="# of Epoches",ylab="Duality Gap")
lines(expe_0_sto_1[1,],expe_0_sto_1[4,],type='b',col=1)
title(main="stochastic=1, Perm")
dev.off()

# Sto=2, Random w Rep
stochastic=2
expe_0_sto_2=sdca(X,y,lambda,epsilon,stochastic,loss)

# Save Plot
pdf("expe_0_sto_2.pdf")
plot(expe_0_sto_2[4,]~expe_0_sto_2[1,],log="y",xlab="# of Epoches",ylab="Duality Gap")
lines(expe_0_sto_2[1,],expe_0_sto_2[4,],type='b',col=1)
title(main="stochastic=2, Random w/ Rep")
dev.off()

# Sto=3, Cyclic
stochastic=3
expe_0_sto_3=sdca(X,y,lambda,epsilon,stochastic,loss)

# Save Plot
pdf("expe_0_sto_3.pdf")
plot(expe_0_sto_3[4,]~expe_0_sto_3[1,],log="y",xlab="# of Epoches",ylab="Duality Gap")
lines(expe_0_sto_3[1,],expe_0_sto_3[4,],type='b',col=1)
title(main="stochastic=3, Cyclic")
dev.off()


########################################################
# Experiment 1: Convergence rate of non-smooth hinge loss
# Fix n and epsilon, vary lambda
X <- rbind(matrix(rnorm(120), , 2), matrix(rnorm(120, mean = 3), , 2))
y <- matrix(c(rep(1, 60), rep(-1, 60)))
epsilon=10^-2
stochastic=1
loss=2
gamma=0

# lambda=10^-2
lambda=10^-2
expe_1_2_lam_001=sdca(X,y,lambda,epsilon,stochastic,loss)

# Save Plot
pdf("expe_1_2_lam_001.pdf")
plot(expe_1_2_lam_001[4,]~expe_1_2_lam_001[1,],log="y",xlab="# of Epoches",ylab="duality gap")
lines(expe_1_2_lam_001[1,],expe_1_2_lam_001[4,],type='b',col=1)
title(main="expe_1_2, lambda=10^-2")
dev.off()

# lambda=10^-3
lambda=10^-3
expe_1_2_lam_0001=sdca(X,y,lambda,epsilon,stochastic,loss)

# Save Plot
pdf("expe_1_2_lam_0001.pdf")
plot(expe_1_2_lam_0001[4,]~expe_1_2_lam_0001[1,],log="y",xlab="# of Epoches",ylab="duality gap")
lines(expe_1_2_lam_0001[1,],expe_1_2_lam_0001[4,],type='b',col=1)
title(main="expe_1_2, lambda=10^-3")
dev.off()

# lambda=10^-4
lambda=10^-4
expe_1_2_lam_00001=sdca(X,y,lambda,epsilon,stochastic,loss)

# Save Plot
pdf("expe_1_2_lam_00001.pdf")
plot(expe_1_2_lam_00001[4,]~expe_1_2_lam_00001[1,],log="y",xlab="# of Epoches",ylab="duality gap")
lines(expe_1_2_lam_00001[1,],expe_1_2_lam_00001[4,],type='b',col=1)
title(main="expe_1_2, lambda=10^-4")
dev.off()

# Look at the dual solution, it should be increasing
pdf("expe_1_2_lam_00001_dual_solution.pdf")
plot(expe_1_2_lam_00001[3,]~expe_1_2_lam_00001[1,],xlab="# of Epoches",ylab="Dual Solution")
lines(expe_1_2_lam_00001[1,],expe_1_2_lam_00001[3,],type='b',col=1)
title(main="expe_1_2_dual_solution, lambda=10^-4")
dev.off()

# So the jump in duality gap should come from primal solution
pdf("expe_1_2_lam_00001_primal_solution.pdf")
plot(expe_1_2_lam_00001[2,]~expe_1_2_lam_00001[1,],xlab="# of Epoches",ylab="Primal Solution")
lines(expe_1_2_lam_00001[1,],expe_1_2_lam_00001[2,],type='b',col=1)
title(main="expe_1_2_primal_solution, lambda=10^-4")
dev.off()


########################################################
# Experiment 2: Convergence rate of smooth hinge loss 
# log scale the vertical axis, we should expect linear convergence rate
X <- rbind(matrix(rnorm(120), , 2), matrix(rnorm(120, mean = 3), , 2))
y <- matrix(c(rep(1, 60), rep(-1, 60)))
epsilon=10^-1
stochastic=1
loss=1
gamma=10^-1

lambda=10
expe_1_1_lam_001=sdca(X,y,lambda,epsilon,stochastic,loss,gamma)

pdf("expesmooth10.pdf")
plot(expe_1_1_lam_001[4,]~expe_1_1_lam_001[1,],log="y",xlab="# of Epoches",ylab="duality gap")
lines(expe_1_1_lam_001[1,],expe_1_1_lam_001[4,],type='b',col=1)
title(main="expe_1_1_smooth, lambda=10")
dev.off()



########################################################
# Experiment 3: Convergence rate of squared loss

# It takes two hours to get the result, so I'll just send you
# the saved result: sqr_loss_duality_gap.Rdata. 
# But the result is reproducible, just set loss=3, and use the
# data as above.

load("sqr_loss_duality_gap.Rdata")

pdf("expesqrlossdualgapall.pdf")
plot(sqr_loss_duality_gap[4,]~sqr_loss_duality_gap[1,],log="y",xlab="# of Epoches",ylab="Duality Gap")
lines(sqr_loss_duality_gap[1,],sqr_loss_duality_gap[4,],type='b',col=1)
title(main="sqr_loss_duality_gap, all")
dev.off()

pdf("expesqrlossprimalall.pdf")
plot(sqr_loss_duality_gap[2,]~sqr_loss_duality_gap[1,],log="y",xlab="# of Epoches",ylab="Primal Sol")
lines(sqr_loss_duality_gap[1,],sqr_loss_duality_gap[2,],type='b',col=1)
title(main="sqr_loss_primal_sol")
dev.off()

pdf("expesqrlossprimal500710.pdf")
plot(sqr_loss_duality_gap[2,500:710]~sqr_loss_duality_gap[1,500:710],log="y",xlab="# of Epoches",ylab="Primal Sol")
lines(sqr_loss_duality_gap[1,500:710],sqr_loss_duality_gap[2,500:710],type='b',col=1)
title(main="sqr_loss_primal_sol_last_200")
dev.off()

pdf("expesqrlossdualall.pdf")
plot(sqr_loss_duality_gap[3,]~sqr_loss_duality_gap[1,],log="y",xlab="# of Epoches",ylab="Dual Sol")
lines(sqr_loss_duality_gap[1,],sqr_loss_duality_gap[3,],type='b',col=1)
title(main="sqr_loss_dual_sol_")
dev.off()

pdf("expesqrlossdual500710.pdf")
plot(sqr_loss_duality_gap[3,500:710]~sqr_loss_duality_gap[1,500:710],log="y",xlab="# of Epoches",ylab="Dual Sol")
lines(sqr_loss_duality_gap[1,500:710],sqr_loss_duality_gap[3,500:710],type='b',col=1)
title(main="sqr_loss_dual_sol_last_200")
dev.off()

  # Lessons learned from Experiment 3:
  # Although SGD might be too agreesive at first, 
  # which is mentioned in the paper,
  # it's also true that SDCA is too conservative.

  # We see from the graph that although dual solution is approaching
  # the dual optimality, the approaching speed is really too slow.
  