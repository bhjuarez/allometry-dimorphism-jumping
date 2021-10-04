
library(nlme)
library(lmtest)
library(phytools)

#run allometry.R first to load and clean data, and calculate anatomical approximations

# Histograms --------------------------------------------------------------

tiff("Figure_1_hist.tiff", width = 8.5, height = 10, units = "in", res = 472)
par(mfrow = c(3,2), mar = c(5, 4, 4, 2) + 0.1)

hist(exp(mal.morph[,2])/exp(fem.morph[,2]), main="a", xlab="M/F SVL", font.lab=2)
abline(v = 1, lwd = 3)

hist(exp(mal.morph[,1])/exp(fem.morph[,1]), main="b", xlab="M/F Mass", font.lab=2)
abline(v = 1, lwd = 3)

hist(exp(mal.morph[,4])/exp(fem.morph[,4]), main="c", xlab="M/F Relative Muscle Volume", font.lab = 2)
abline(v = 1, lwd = 3)

hist(exp(mal.morph[,3])/exp(fem.morph[,3]), main="d", xlab= expression(bold("M/F Relative L"*''[com])), font.lab = 2)
abline(v = 1, lwd = 3)

hist(exp(mal.perf[,2])/exp(fem.perf[,2]), main = "e", xlab="M/F Mass-specific Energy", font.lab = 2)
abline(v = 1, lwd = 3)

hist(exp(mal.perf[,1])/exp(fem.perf[,1]), main="f", xlab="M/F Velocity", font.lab = 2)
abline(v = 1, lwd = 3)

dev.off()

#Repeat PGLS section of Allometry.R for easy comparison against PGLS w/ error code below
# PGLS ------------------------------------------------------------

#SVL
fit <- gls(M~F, correlation=corBrownian(phy=tree), data=data.frame(M=mal.morph[,2],F=fem.morph[,2]))
anova(fit)
summary(fit)  #0.96
T <- (coef(fit)[2]-1) / summary(fit)$tTable[2,2]
prob.T <- dt(T, df = (nrow(mal.perf) - 1))
T
prob.T # NOT sig diff from 1

#for r^2
gdf <- geomorph.data.frame(M=mal.morph[,2],F=fem.morph[,2])
fit <- lm.rrpp(M~F, Cov = vcv.phylo(tree), data = gdf, iter=9999)
anova(fit)

par(mfrow = c(3,2))
plot(fem.morph[,2], mal.morph[,2], pch = 16, xlab = "Log Female SVL", ylab = "Log Male SVL", main = "a", asp = 1, font.lab = 2)
abline(fit, lwd = 2, col = "black")
abline(0,1, lwd = 2, lty= 2, col = "red")

#Mass
fit <- gls(M~F, correlation=corBrownian(phy=tree), data=data.frame(M=mal.morph[,1],F=fem.morph[,1]))
anova(fit)
summary(fit) #0.95
T <- (coef(fit)[2]-1) / summary(fit)$tTable[2,2]
prob.T <- dt(T, df = (nrow(mal.perf) - 1))
T
prob.T # NOT sig diff from 1

#for r^2
gdf <- geomorph.data.frame(M=mal.morph[,1],F=fem.morph[,1])
fit <- lm.rrpp(M~F, Cov = vcv.phylo(tree), data = gdf, iter=9999)
anova(fit)

plot(fem.morph[,1], mal.morph[,1], pch = 16, xlab = "Log Female Mass", ylab = "Log Male Mass", main = "b", asp = 1, font.lab = 2)
abline(fit, lwd = 2, col = "black")
abline(0,1, lwd = 2, lty= 2, col = "red")

#Relative muscle volume
fit <- gls(M~F, correlation=corBrownian(phy=tree), data=data.frame(M=mal.morph[,4],F=fem.morph[,4]))
anova(fit)
summary(fit) #0.95
T <- (coef(fit)[2]-1) / summary(fit)$tTable[2,2]
prob.T <- dt(T, df = (nrow(mal.perf) - 1))
T
prob.T # NOT sig diff from 1

plot(fem.morph[,4], mal.morph[,4], pch = 16, xlab = "Log Female Relative Muscle Volume", ylab = "Log Male Relative Muscles Volume", main = "c", asp = 1, font.lab = 2)
abline(fit, lwd = 2, col = "black")
abline(0,1, lwd = 2, lty= 2, col = "red")

gdf <- geomorph.data.frame(M=mal.morph[,4],F=fem.morph[,4])
fit <- lm.rrpp(M~F, Cov = vcv.phylo(tree), data = gdf, iter=9999)
anova(fit)

#Relative Lcom
fit <- gls(M~F, correlation=corBrownian(phy=tree), data=data.frame(M=mal.morph[,3],F=fem.morph[,3]))
anova(fit)
summary(fit) #0.95
T <- (coef(fit)[2]-1) / summary(fit)$tTable[2,2]
prob.T <- dt(T, df = (nrow(mal.perf) - 1))
T
prob.T # sig

gdf <- geomorph.data.frame(M=mal.morph[,3],F=fem.morph[,3])
fit <- lm.rrpp(M~F, Cov = vcv.phylo(tree), data = gdf, iter=9999)
anova(fit)


plot(fem.morph[,3], mal.morph[,3], pch = 16, xlab = expression(bold("Log Female Relative L"*''[com])), ylab = expression(bold("Log Male Relative L"*''[com])), main = "d", asp = 1)
abline(fit, lwd = 2, col = "black")
abline(0,1, lwd = 2, lty= 2, col = "red")

#Mass-specific Energy
fit <- gls(M~F, correlation=corBrownian(phy=tree), data=data.frame(M=mal.perf[,2],F=fem.perf[,2]))
anova(fit)
summary(fit) #0.89
T <- (coef(fit)[2]-1) / summary(fit)$tTable[2,2]
prob.T <- dt(T, df = (nrow(mal.perf) - 1))
T
prob.T # sig

gdf <- geomorph.data.frame(M=mal.perf[,2],F=fem.perf[,2])
fit <- lm.rrpp(M~F, Cov = vcv.phylo(tree), data = gdf, iter=9999)
anova(fit)

plot(fem.perf[,2], mal.perf[,2], pch = 16, xlab = "Log Female Mass-specific Energy", ylab = "Log Male Mass-Specific Energy", main = "e", asp = 1, font.lab = 2)
abline(fit, lwd = 2, col = "black")
abline(0,1, lwd = 2, lty= 2, col = "red")

#Velocity
fit <- gls(M~F, correlation=corBrownian(phy=tree), data=data.frame(M=mal.perf[,1],F=fem.perf[,1]))
anova(fit)
summary(fit) #0.87
T <- (coef(fit)[2]-1) / summary(fit)$tTable[2,2]
prob.T <- dt(T, df = (nrow(mal.perf) - 1))
T
prob.T # sig

plot(fem.perf[,1], mal.perf[,1], pch = 16, xlab = "Log Female Velocity", ylab = "Log Male Velocity", main = "f", asp = 1, font.lab = 2)
abline(fit, lwd = 2, col = "black")
abline(0,1, lwd = 2, lty= 2, col = "red")

gdf <- geomorph.data.frame(M=mal.perf[,1],F=fem.perf[,1])
fit <- lm.rrpp(M~F, Cov = vcv.phylo(tree), data = gdf, iter=9999)
anova(fit)


# PGLS w/ error -----------------------------------------------------------

cxy <- rep(0, 146)
names(cxy) <- rownames(mal.morph)

#SVL
fit1.SVL.ols <- pgls.Ives(tree, X = fem.morph[,2], y = mal.morph[,2], Vx = cxy, Vy = cxy, Cxy = cxy, fixed.b1 = 1)
fit1.SVL.ols #converged; int = -0.1188584, lik = -21.172028

fit.SVL.ols <- pgls.Ives(tree, X = fem.morph[,2], y = mal.morph[,2], Vx = cxy, Vy = cxy, Cxy = cxy)
fit.SVL.ols #converged; slope = 0.9567683, likelihood = -20.572475

fit1.SVL <- pgls.Ives(tree, X = fem.morph[,2], y = mal.morph[,2], Vx = fem.morph[,6], Vy = mal.morph[,6], Cxy = cxy, fixed.b1 = 1)
fit1.SVL #converged; int = -0.1168917,lik = -8.402692

fit.SVL <- pgls.Ives(tree, X = fem.morph[,2], y = mal.morph[,2], Vx = fem.morph[,6], Vy = mal.morph[,6], Cxy = cxy)
fit.SVL #converged; slope = 0.94948, lik  = -7.108423

lrtest(fit1.SVL.ols, fit.SVL.ols, fit1.SVL,fit.SVL) #fit1.SVL model of isometry w/ error preferred
lrtest(fit1.SVL.ols, fit1.SVL)
lrtest(fit1.SVL, fit.SVL)

#Mass
fit1.mass.ols <- pgls.Ives(tree, X = fem.morph[,1], y = mal.morph[,1], Vx = cxy, Vy = cxy, Cxy = cxy, fixed.b1 = 1)
fit1.mass.ols #converged; likelihood = -362.132

fit.mass.ols <- pgls.Ives(tree, X = fem.morph[,1], y = mal.morph[,1], Vx = cxy, Vy = cxy, Cxy = cxy)
fit.mass.ols #converged; slope = 0.9531646, likelihood = -361.456

fit1.mass <- pgls.Ives(tree, X = fem.morph[,1], y = mal.morph[,1], Vx = fem.morph[,5], Vy = mal.morph[,5], Cxy = cxy, fixed.b1 = 1)
fit1.mass #converged; likelihood = -346.781

fit.mass <- pgls.Ives(tree, X = fem.morph[,1], y = mal.morph[,1], Vx = fem.morph[,5], Vy = mal.morph[,5], Cxy = cxy)
fit.mass #converged; slope = 0.9753095, likelihood = -346.508

lrtest(fit1.mass.ols, fit.mass.ols, fit1.mass, fit.mass)
lrtest(fit1.mass.ols,fit1.mass) #fit1.mass model of isometry w/ error preferred
lrtest(fit1.mass,fit.mass)

#Relative Leg Length
fit1.L.ols <- pgls.Ives(tree, X = fem.morph[,3], y = mal.morph[,3], Vx = cxy, Vy = cxy, Cxy = cxy, fixed.b1 = 1)
fit1.L.ols #converged; int = 0.04384241,likelihood = 363.125

fit.L.ols <- pgls.Ives(tree, X = fem.morph[,3], y = mal.morph[,3], Vx = cxy, Vy = cxy, Cxy = cxy)
fit.L.ols #converged slope = 0.8512256; likelihood = 374.069105

fit1.L <- pgls.Ives(tree, X = fem.morph[,3], y = mal.morph[,3], Vx = fem.morph[,7], Vy = mal.morph[,7], Cxy = cxy, fixed.b1 = 1)
fit1.L # converged; int = 0.01426456, likelihood = 377.811433

fit.L <- pgls.Ives(tree, X = fem.morph[,3], y = mal.morph[,3], Vx = fem.morph[,7], Vy = mal.morph[,7], Cxy = cxy)
fit.L #converged; slope = 0.8922699; likelihood = 382.453996

lrtest(fit1.L.ols, fit.L.ols, fit1.L, fit.L) #Hypoallometry with error
lrtest(fit1.L,fit.L)

plot(fem.morph[,3], mal.morph[,3], pch = 16, xlab = expression(bold("Log Female Relative L"*''[com])), ylab = expression(bold("Log Male Relative L"*''[com])), main = "D", asp = 1)
abline(fit.L$beta, lwd = 2, col = "black")
abline(0,1, lwd = 2, lty= 2, col = "red")

#Relative Muscle Volume
fit1.muscle.ols <- pgls.Ives(tree, X = fem.morph[,4], y = mal.morph[,4], Vx = cxy, Vy = cxy, Cxy = cxy, fixed.b1 = 1)
fit1.muscle.ols #converged; int = -0.1222523, likelihood = -380.503508

fit.muscle.ols <- pgls.Ives(tree, X = fem.morph[,4], y = mal.morph[,4], Vx = cxy, Vy = cxy, Cxy = cxy)
fit.muscle.ols #converged; slope = 0.9282677, likelihood = -378.737565

fit1.muscle <- pgls.Ives(tree, X = fem.morph[,4], y = mal.morph[,4], Vx = fem.morph[,8], Vy = mal.morph[,8], Cxy = cxy, fixed.b1 = 1)
fit1.muscle #converged; int = -0.1279357, lik = -371.062725

fit.muscle <- pgls.Ives(tree, X = fem.morph[,4], y = mal.morph[,4], Vx = fem.morph[,8], Vy = mal.morph[,8], Cxy = cxy)
fit.muscle #converged; slope = 0.9319139, likelihood = -368.870322

lrtest(fit1.muscle.ols, fit.muscle.ols, fit1.muscle,fit.muscle)
lrtest(fit1.muscle.ols, fit1.muscle)
lrtest(fit1.muscle, fit.muscle) #Hypoallometry with error

plot(fem.morph[,4], mal.morph[,4], pch = 16, xlab = expression(bold("Log Female Relative Muscle Volume")), ylab = expression(bold("Log Male Relative Muscle Volume")), main = "", asp = 1)
abline(fit.muscle$beta, lwd = 2, col = "black")
abline(0,1, lwd = 2, lty= 2, col = "red")

#Energy
fit1.nrg.ols <- pgls.Ives(tree, X = fem.perf[,2], y = mal.perf[,2], Vx = cxy, Vy = cxy, Cxy = cxy, fixed.b1 = 1)
fit1.nrg.ols #converged; int = 0.127638, lik = 3.823156

fit.nrg.ols <- pgls.Ives(tree, X = fem.perf[,2], y = mal.perf[,2], Vx = cxy, Vy = cxy, Cxy = cxy)
fit.nrg.ols #converged; slope = 0.8613932, lik = 10.16419

fit1.nrg <- pgls.Ives(tree, X = fem.perf[,2], y = mal.perf[,2], Vx = fem.perf[,4], Vy = mal.perf[,4], Cxy = cxy, fixed.b1 = 1)
fit1.nrg #converged; int = , lik = 

fit.nrg <- pgls.Ives(tree, X = fem.perf[,2], y = mal.perf[,2], Vx = fem.perf[,4], Vy = mal.perf[,4], Cxy = cxy)
fit.nrg #converged; slope = , lik = 

lrtest(fit1.nrg.ols,fit.nrg.ols, fit1.nrg, fit.nrg) #sig
lrtest(fit.nrg.ols, fit.nrg)
lrtest(fit.nrg.ols, fit1.nrg)
#hypoallometry ols

lrtest(fit.nrg.ols, fit.nrg)

plot(fem.perf[,2], mal.perf[,2], pch = 16, xlab = expression(bold("Log Female Jumping Energy")), ylab = expression(bold("Log Male Jumping Energy")), main = "", asp = 1)
abline(fit.nrg$beta, lwd = 2, col = "black")
abline(0,1, lwd = 2, lty= 2, col = "red")

#Velocity

fit1.vel.ols <- pgls.Ives(tree, X = fem.perf[,1], y = mal.perf[,1], Vx = cxy, Vy = cxy, Cxy = cxy, fixed.b1 = 1)
fit1.vel.ols #converged; int = 0.05718727, likelihood = 213.675755

fit.vel.ols <- pgls.Ives(tree, X = fem.perf[,1], y = mal.perf[,1], Vx = cxy, Vy = cxy, Cxy = cxy)
fit.vel.ols #converged; slope = 0.866272, likelihood = 219.448684

fit1.vel <- pgls.Ives(tree, X = fem.perf[,1], y = mal.perf[,1], Vx = fem.perf[,3], Vy = mal.perf[,3], Cxy = cxy, fixed.b1 = 1)
fit1.vel #converged; int = 0.03427616, lik = 219.198361

fit.vel <- pgls.Ives(tree, X = fem.perf[,1], y = mal.perf[,1], Vx = fem.perf[,3], Vy = mal.perf[,3], Cxy = cxy)
fit.vel #converged; slope = 0.9216289, likelihood = 220.614765

lrtest(fit1.vel.ols, fit.vel.ols, fit1.vel, fit.vel)
lrtest(fit.vel.ols, fit.vel)
fit.vel.ols #hypoallometry no error

plot(fem.perf[,1], mal.perf[,1], pch = 16, xlab = expression(bold("Log Female Relative L"*''[com])), ylab = expression(bold("Log Male Relative L"*''[com])), main = "D", asp = 1)
abline(fit.vel$beta, lwd = 2, col = "black")
abline(0,1, lwd = 2, lty= 2, col = "red")
