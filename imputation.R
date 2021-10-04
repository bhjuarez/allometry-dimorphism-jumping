dat <- read.csv("/Users/Bryan/Dropbox/frogs/Data/morph/museums_morph.csv")

#remove tag mass
dat$mass <- dat$mass - 0.1 * dat$X.tags
dat$mass <- log(dat$mass)

# Imputation from Baken and Adams 2019 Ecology and Evolution
# Estimating missing leg muscle diameters using rest of dataset

      estimate.dist <- function(A, phy = NULL){
          if(any(is.na(A)) == FALSE){stop("No missing data.")}
          spec.NA <- which(complete.cases(A) == F)
          land.NA <- which(colSums(is.na(A)) > 0)  
          A2 <- A
          for (i in 1:length(spec.NA)){ 
            missing.dist <- which(is.na(A2[spec.NA[i], ]))  
            x <- A2[-spec.NA, -missing.dist] 
            y <- A2[-spec.NA, missing.dist]
            S12 <- cov(x, y)
            pls <- svd(S12) 
            U <- pls$u
            V <- pls$v
            XScores <- x%*%U
            YScores <- y%*%V
            beta <- coef(lm(YScores ~ XScores))
            miss.xsc <- c(1,A2[spec.NA[i], -missing.dist]%*%U)
            miss.ysc <- miss.xsc%*%beta
            pred.val <- miss.ysc%*%t(V)
            for (j in 1:length(V)){
              A2[spec.NA[i], missing.dist[j]] <- pred.val[j]    
            }
          }
          return(A2)
      }
      
dat2 <- as.matrix(dat[, 9:19])
dat2 <- estimate.dist(dat2)
dat2 <- as.data.frame(dat2)

dat <- cbind(dat[, 1:8], dat2, dat[, 20:24])# imputed ___% of the data
dat$mass <- exp(dat$mass)

