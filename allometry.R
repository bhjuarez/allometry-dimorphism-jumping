library(ape)
library(geiger)
library(nlme)
library(geomorph)
library(phytools)
library(RColorBrewer)

source("/Users/Bryan/Dropbox/frogs/allometry/Code/imputation.R", echo = F, print.eval = F)

#Replacements used by Feng et al. 2017 PNAS (Chronogram_309sp.tre, found in their supplement)
dat[which(dat$species == "Incilius valliceps"), 1] <- "Incilius nebulifer"
dat[which(dat$species == "Leptodactylus ocellatus"), 1] <- "Leptodactylus latrans"
dat[which(dat$species == "Xenopus fraseri"), 1] <- "Xenopus kobeli"

dat[which(dat$species == "Colostethus jacobuspetersi"), 1] <- "Hyloxalus jacobuspetersi"
dat[which(dat$species == "Engystomops pustulosus"), 1] <- "Physalaemus pustulosus"
dat[which(dat$species == "Hylarana guentheri"), 1] <- "Sylvirana guentheri"
dat[which(dat$species == "Hylarana latouchii"), 1] <- "Papurana latouchii"
dat[which(dat$species == "Pedostibes hosii"), 1] <- "Rentapia hosii"
dat[which(dat$species == "Phrynoidis asper"), 1] <- "Phrynoidis aspera"
dat[which(dat$species == "Platymantis pelewensis"), 1] <- "Cornufer pelewensis"
dat[which(dat$species == "Bufo Bufo"), 1] <- "Bufo bufo" #capitalization
dat[which(dat$species == "Odorrana hosii "), 1] <- "Odorrana hosii" #extra space

# Calculate proxies -------------------------------------------------------
#for males and females
fem.dat <- dat[which(dat$sex == "F"),]
attach(fem.dat)
source("/Users/Bryan/Dropbox/frogs/Allometry/Code/calc_appxs.R")
FV <- V.proxy
FE <- E.proxy
FV.v <- V.proxy.v
FE.v <- E.proxy.v
F.perf <- as.matrix(cbind(FV, FE, FV.v, FE.v))
F.morph <- as.matrix(
  cbind(tapply(log(fem.dat$mass/1000), species, mean), 
        tapply(log(fem.dat$svl/1000), species, mean), 
        tapply(log((0.67 
                    * fem.dat$sac
                    + fem.dat$toe
                    + fem.dat$calc
                    + fem.dat$tibl 
                    + fem.dat$femur) / fem.dat$svl), species, mean), 
        tapply(log(
          (2 * pi * (fem.dat$diam1low / 2 / 1000) * (fem.dat$diam2low / 2 / 1000) * (fem.dat$tibl / 1000)) + (2 * pi * (fem.dat$diam1up / 2 / 1000) * (fem.dat$diam2up / 2 / 1000) * (fem.dat$femur / 1000)) / fem.dat$mass/1000
        ), species, mean),
        tapply(log(fem.dat$mass/1000), species, var), 
        tapply(log(fem.dat$svl/1000), species, var), 
        tapply(log((0.67 
                    * fem.dat$sac
                    + fem.dat$toe
                    + fem.dat$calc
                    + fem.dat$tibl 
                    + fem.dat$femur) / fem.dat$svl), species, var), 
        tapply(log(
          (2 * pi * (fem.dat$diam1low / 2 / 1000) * (fem.dat$diam2low / 2 / 1000) * (fem.dat$tibl / 1000)) + (2 * pi * (fem.dat$diam1up / 2 / 1000) * (fem.dat$diam2up / 2 / 1000) * (fem.dat$femur / 1000)) / fem.dat$mass/1000
        ), species, var)
  )
)
colnames(F.morph) <- c("mass", "SVL", "rel.L", "rel.muscle", "mass.v", "SVL.v", "rel.L.v", "rel.muscle.v")
detach(fem.dat)

mal.dat <- dat[which(dat$sex == "M"),]
attach(mal.dat)
source("/Users/Bryan/Dropbox/frogs/Allometry/Code/calc_appxs.R")
MV <- V.proxy
ME <- E.proxy
MV.v <- V.proxy.v
ME.v <- E.proxy.v
M.perf <- as.matrix(cbind(MV, ME, MV.v, ME.v))
M.morph <- as.matrix(
  cbind(tapply(log(mal.dat$mass/1000), species, mean), 
        tapply(log(mal.dat$svl/1000), species, mean), 
        tapply(log((0.67 
                    * mal.dat$sac
                    + mal.dat$toe
                    + mal.dat$calc
                    + mal.dat$tibl 
                    + mal.dat$femur) / mal.dat$svl), species, mean), 
        tapply(log(
          (2 * pi * (mal.dat$diam1low / 2 / 1000) * (mal.dat$diam2low / 2 / 1000) * (mal.dat$tibl / 1000)) + (2 * pi * (mal.dat$diam1up / 2 / 1000) * (mal.dat$diam2up / 2 / 1000) * (mal.dat$femur / 1000)) / mal.dat$mass/1000
        ), species, mean),
        tapply(log(mal.dat$mass/1000), species, var), 
        tapply(log(mal.dat$svl/1000), species, var), 
        tapply(log((0.67 
                    * mal.dat$sac
                    + mal.dat$toe
                    + mal.dat$calc
                    + mal.dat$tibl 
                    + mal.dat$femur) / mal.dat$svl), species, var), 
        tapply(log(
          (2 * pi * (mal.dat$diam1low / 2 / 1000) * (mal.dat$diam2low / 2 / 1000) * (mal.dat$tibl / 1000)) + (2 * pi * (mal.dat$diam1up / 2 / 1000) * (mal.dat$diam2up / 2 / 1000) * (mal.dat$femur / 1000)) / mal.dat$mass/1000
        ), species, var)
  )
)
colnames(M.morph) <- c("mass", "SVL", "rel.L", "rel.muscle", "mass.v", "SVL.v", "rel.L.v", "rel.muscle.v")
detach(mal.dat)

# Data matching ----------------------------------------------------------------

tree <- read.tree("/Users/Bryan/Dropbox/frogs/dimorphism/Data/Chronogram_309sp.tre")
tree$tip.label <- sub("_.*", "", sub("_", " ", tree$tip.label))

micro.dat <- read.csv("/Users/Bryan/Dropbox/frogs/dimorphism/Data/micro_FINAL.csv")
#updated taxonomy
micro.dat[which(micro.dat$species == "Incilius valliceps"), 1] <- "Incilius nebulifer"
micro.dat[which(micro.dat$species == "Leptodactylus ocellatus"), 1] <- "Leptodactylus latrans"
micro.dat[which(micro.dat$species == "Xenopus fraseri"), 1] <- "Xenopus kobeli"

micro <- micro.dat$micro.separate
names(micro) <- micro.dat$species

#Remove single species that does not jump (: Pipa pipa) from all data
del <- micro.dat$species[which(micro.dat$locomotion == "Remove")] 
micro <- micro[!(names(micro) %in% del)]

F.perf <- F.perf[!(rownames(F.perf) %in% del),]
M.perf <- M.perf[!(rownames(M.perf) %in% del),]
F.morph <- F.morph[!(rownames(F.morph) %in% del),]
M.morph <- M.morph[!(rownames(M.morph) %in% del),]

fem.perf <- treedata(tree, F.perf, sort = T)$data
mal.perf <- treedata(tree, M.perf, sort = T)$data

fem.morph <- treedata(tree, F.morph, sort = T)$data
mal.morph <- treedata(tree, M.morph, sort = T)$data

tree <- treedata(tree, mal.perf, sort = T)$phy

# Tree plot ---------------------------------------------------------------


pdf("Figure_1_tree.pdf", width = 8.5, height = 10)

layout(matrix(1:2, 1, 2), widths = c(0.7, 0.3))

brewcols <- brewer.pal(8, "Dark2")
cols <- c(rep(brewcols[1], 64), rep(brewcols[2], 61), rep(brewcols[3], 2), rep(brewcols[4], 1), rep(brewcols[5], 11), rep(brewcols[6], 3), rep(brewcols[7], 2), rep(brewcols[8], 2))
plot(tree, show.tip.label = T, cex = 0.01, tip.color = "white", x.lim = c(0,2.5))
legend(x = 0, y = 150, legend = c("Ranoidea", "Hyloidea", "Myobatrachidae", "Heleophrynidae", "Pelobatoidea", "Pipoidea", "Discoglossoidea", "Leiopelmatoidea"), fill = brewcols, text.font = 2)
cladelabels(tree, node=229,"Ranoidea")
cladelabels(tree, node=169,"Hyloidea")
cladelabels(tree, node=155,"Pelobatoidea", offset = 2.5, cex = 0.8)
cladelabels(tree, node=152,"Pipoidea", cex = 0.65)

plot.new()
plot.window(xlim=c(-0.01,0.01),
            ylim=get("last_plot.phylo",envir=.PlotPhyloEnv)$y.lim)
text(rep(0,length(tree$tip.label)),1:Ntip(tree),
     tree$tip.label,font=4, cex = 0.4, col = rev(cols))

dev.off()

# Histograms --------------------------------------------------------------

tiff("Figure_3_histograms.tiff", width = 8.5, height = 10, units = "in", res = 472)
par(mfrow = c(3,2), mar = c(5, 4, 4, 2) + 0.1)

#SVL
ave <- mean(mal.morph[,2] - fem.morph[,2])
se <- sqrt(var(mal.morph[,2] - fem.morph[,2])/146)
upp <- mean(mal.morph[,2] - fem.morph[,2]) + 1.96 * sqrt(var(mal.morph[,2] - fem.morph[,2])/146)
low <- mean(mal.morph[,2] - fem.morph[,2]) - 1.96 * sqrt(var(mal.morph[,2] - fem.morph[,2])/146)
ave; se; upp; low

border <- max(hist(mal.morph[,2] - fem.morph[,2], plot = F)$counts)
hist(mal.morph[,2] - fem.morph[,2], main="a", xlab="Log M/F SVL", font.lab=2)
abline(v = 0, lwd = 4)
abline(v = upp, lwd = 3, col = "red")
abline(v = low, lwd = 3, col = "red")
legend(x = "topleft", "95% CI", text.col = "red", cex = 0.8, , text.font = 2)

#Mass
ave <- mean(mal.morph[,1] - fem.morph[,1])
se <- sqrt(var(mal.morph[,1] - fem.morph[,1])/146)
upp <- mean(mal.morph[,1] - fem.morph[,1]) + 1.96 * sqrt(var(mal.morph[,1] - fem.morph[,1])/146)
low <- mean(mal.morph[,1] - fem.morph[,1]) - 1.96 * sqrt(var(mal.morph[,1] - fem.morph[,1])/146)
ave; se; upp; low

border <- max(hist(mal.morph[,1] - fem.morph[,1], plot = F)$counts)
hist(mal.morph[,1] - fem.morph[,1], main="b", xlab="Log M/F Mass", font.lab=2)
abline(v = 0, lwd = 4)
abline(v = upp, lwd = 3, col = "red")
abline(v = low, lwd = 3, col = "red")
legend(x = "topleft", "95% CI", text.col = "red", cex = 0.8, , text.font = 2)

#Rel Muscle Volume
ave <- mean(   mal.morph[,4] - fem.morph[,4])
se <- sqrt(var(mal.morph[,4] - fem.morph[,4])/146)
upp <- mean(   mal.morph[,4] - fem.morph[,4]) + 
  1.96 * sqrt(var(mal.morph[,4] - fem.morph[,4])/146)
low <- mean(   mal.morph[,4] - fem.morph[,4]) - 
  1.96 * sqrt(var(mal.morph[,4] - fem.morph[,4])/146)
ave; se; upp; low

border <- max(hist(mal.morph[,4] - fem.morph[,4], plot = F)$counts)
hist(mal.morph[,4] - fem.morph[,4], main="c", xlab="Log M/F Relative Muscle Volume", font.lab = 2)
abline(v = 0, lwd = 4)
abline(v = upp, lwd = 3, col = "red")
abline(v = low, lwd = 3, col = "red")
legend(x = "topleft", "95% CI", text.col = "red", cex = 0.8, , text.font = 2)

#Rel Lcom
ave <- mean(   mal.morph[,3] - fem.morph[,3])
se <- sqrt(var(mal.morph[,3] - fem.morph[,3])/146)
upp <- mean(   mal.morph[,3] - fem.morph[,3]) + 
  1.96 * sqrt(var(mal.morph[,3] - fem.morph[,3])/146)
low <- mean(   mal.morph[,3] - fem.morph[,3]) - 
  1.96 * sqrt(var(mal.morph[,3] - fem.morph[,3])/146)
ave; se; upp; low

border <- max(hist(mal.morph[,3] - fem.morph[,3], plot = F)$counts)
hist(mal.morph[,3] - fem.morph[,3], main="d", xlab= expression(bold("Log M/F Relative L"*''[com])), font.lab = 2)
abline(v = 0, lwd = 4)
abline(v = upp, lwd = 3, col = "red")
abline(v = low, lwd = 3, col = "red")
legend(x = "topleft", "95% CI", text.col = "red", cex = 0.8, , text.font = 2)

#Mass-specific energy
ave <- mean(   mal.perf[,2] - fem.perf[,2])
se <- sqrt(var(mal.perf[,2] - fem.perf[,2])/146)
upp <- mean(   mal.perf[,2] - fem.perf[,2]) + 
  1.96 * sqrt(var(mal.perf[,2] - fem.perf[,2])/146)
low <- mean(   mal.perf[,2] - fem.perf[,2]) - 
  1.96 * sqrt(var(mal.perf[,2] - fem.perf[,2])/146)
ave; se; upp; low

border <- max(hist(mal.perf[,2] - fem.perf[,2], plot = F)$counts)
hist(mal.perf[,2] - fem.perf[,2], main = "e", xlab="Log M/F Mass-specific Energy", font.lab = 2)
abline(v = 0, lwd = 4)
abline(v = upp, lwd = 3, col = "red")
abline(v = low, lwd = 3, col = "red")
legend(x = "topleft", "95% CI", text.col = "red", cex = 0.8, , text.font = 2)

#Velocity
ave <- mean(   mal.perf[,1] - fem.perf[,1])
se <- sqrt(var(mal.perf[,1] - fem.perf[,1])/146)
upp <- mean(   mal.perf[,1] - fem.perf[,1]) + 
  1.96 * sqrt(var(mal.perf[,1] - fem.perf[,1])/146)
low <- mean(   mal.perf[,1] - fem.perf[,1]) - 
  1.96 * sqrt(var(mal.perf[,1] - fem.perf[,1])/146)
ave; se; upp; low

border <- max(hist(mal.perf[,1] - fem.perf[,1], plot = F)$counts)
hist(mal.perf[,1] - fem.perf[,1], main="f", xlab="Log M/F Velocity", font.lab = 2)
abline(v = 0, lwd = 4)
abline(v = upp, lwd = 3, col = "red")
abline(v = low, lwd = 3, col = "red")
legend(x = "topleft", "95% CI", text.col = "red", cex = 0.8, , text.font = 2)

dev.off()

# PGLS ------------------------------------------------------------

#for CI of slope for models below
summary(fit)$tTable[2,1]
summary(fit)$tTable[2,1] + 1.96*summary(fit)$tTable[2,2]
summary(fit)$tTable[2,1] - 1.96*summary(fit)$tTable[2,2]

tiff("Figure_4_regressions.tiff", width = 8.5, height = 10, units = "in", res = 472)

#SVL
fit <- gls(M~F, correlation=corBrownian(phy=tree), data=data.frame(M=mal.morph[,2],F=fem.morph[,2]))
anova(fit)
summary(fit)  #0.96
T <- (coef(fit)[2]-1) / summary(fit)$tTable[2,2]
prob.T <- dt(T, df = (nrow(mal.perf) - 1))
T
prob.T # NOT sig diff from 1

#for r^2
gdf <-geomorph.data.frame(M=mal.morph[,2],F=fem.morph[,2])
#fit <-lm.rrpp(M~F, Cov = vcv.phylo(tree), data = gdf, iter=9999)
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
#fit <- lm.rrpp(M~F, Cov = vcv.phylo(tree), data = gdf, iter=9999)
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
#fit <- lm.rrpp(M~F, Cov = vcv.phylo(tree), data = gdf, iter=9999)
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
#fit <- lm.rrpp(M~F, Cov = vcv.phylo(tree), data = gdf, iter=9999)
anova(fit)

dev.off()

#energy
fit <- gls(SD~SSD, correlation=corBrownian(phy=tree), data=data.frame(SD=mal.perf[,2]-fem.perf[,2],SSD=mal.morph[,1]-fem.morph[,1]))
anova(fit)
summary(fit) #0.87

gdf <- geomorph.data.frame(SD=mal.perf[,2]-fem.perf[,2],SSD=mal.morph[,1]-fem.morph[,1])
fit <- lm.rrpp(SD~SSD, Cov = vcv.phylo(tree), data = gdf, iter=9999)
anova(fit)

#velocity
fit <- gls(SD~SSD, correlation=corBrownian(phy=tree), data=data.frame(SD=mal.perf[,1]-fem.perf[,1],SSD=mal.morph[,2]-fem.morph[,2]))
anova(fit)
summary(fit) #0.87

gdf <- geomorph.data.frame(SD=mal.perf[,1]-fem.perf[,1],SSD=mal.morph[,2]-fem.morph[,2])
fit <- lm.rrpp(SD~SSD, Cov = vcv.phylo(tree), data = gdf, iter=9999)
anova(fit)

#Mass
physignal(mal.morph[,1], tree,iter = 9999)
physignal(fem.morph[,1], tree, iter = 9999)

#SVL
physignal(mal.morph[,2], tree,iter = 9999)
physignal(fem.morph[,2], tree, iter = 9999)

#Relative Lcom
physignal(mal.morph[,3], tree,iter = 9999)
physignal(fem.morph[,3], tree, iter = 9999)

#Relative muscle volume
physignal(mal.morph[,4], tree,iter = 9999)
physignal(fem.morph[,4], tree, iter = 9999)

#Performance
physignal(mal.perf, tree,iter = 9999)
physignal(fem.perf, tree, iter = 9999)

#Lcom
compare.multi.evol.rates(as.matrix(cbind(mal.morph[,3], fem.morph[,3])), gp = c("M", "F"), phy = tree, iter = 9999)

#Performance
compare.multi.evol.rates(as.matrix(cbind(mal.perf[,1:2], fem.perf[,1:2])), gp = c("M", "M", "F", "F"), phy = tree, iter = 9999)

plot(tree, cex = 0.2)
nodelabels(cex = 0.3)

