#Mean species performance
Mass <- mass / 1000 #convert to kilograms
Svl <- svl / 1000 #convert to m

#Mean species proxies

CSA <- ((2 * 2 * (1/3) * pi * (diam1low/1000 / 2) * (diam2low/1000 / 2) * (tibl/1000 / 2)) + (2 * 2 * (1/3) * pi * (diam1up/1000 / 2) * (diam2up/1000 / 2) * (femur/1000 / 2))) ^ (2 / 3)
L <- (0.67 * sac + femur + tibl + calc + toe) / 1000

V.proxy <- tapply(log(sqrt(2 * L * CSA / Mass)), species, mean)
E.proxy <- tapply(log(L * CSA / Mass), species, mean)

V.proxy.v <- tapply(log(sqrt(2 * L * CSA / Mass)), species, var)
E.proxy.v <- tapply(log(L * CSA / Mass), species, var)
Mass <- tapply(log(Mass), species, mean) #kg

V.proxy <- V.proxy[complete.cases(V.proxy)]
E.proxy <- E.proxy[complete.cases(E.proxy)]
V.proxy.v <- V.proxy.v[complete.cases(V.proxy.v)]
E.proxy.v <- E.proxy.v[complete.cases(E.proxy.v)]
Mass <- Mass[complete.cases(Mass)]

