## code to determine design values and informative priors
## this procedure involves processing data from the ENIGH 2018 survey

## BEGIN SETUP ##

## Process to obtain data file from ENIGH website
## Step 1: Go to this website: https://www.inegi.org.mx/programas/enigh/nc/2018/#Datos_abiertos
## Step 2: Click the "Download CSV" button to download a .zip file
## Step 3: Open .zip file
## Step 4: Open the "conjunto_de_datos_concentradohogar_enigh_2018_ns" subfolder
## Step 5: Open the "conjunto_de_datos" subfolder
## Step 6: Extract the lone .csv file: "conjunto_de_datos_concentradohogar_enigh_2018_ns.csv" and move it to the working directory

## read data set
enigh <- read.csv("conjunto_de_datos_concentradohogar_enigh_2018_ns.csv")

## create a numeric character for region
enigh$region <- as.numeric(substr(as.character(enigh$ubica_geo),1,nchar(as.character(enigh$ubica_geo))-3))

## keep only the households from Aguascalientes (Region = 1):
aguas <- enigh[enigh$region %in% c(1),]

## keep only columns of interest
aguas <- aguas[c("region","factor","est_socio","educa_jefe", "tot_integ", "alimentos", "sexo_jefe")]

## convert column titles to English
names(aguas)[4] <- "education"
names(aguas)[5] <- "total_people"
names(aguas)[6] <- "total_food"
names(aguas)[7] <- "sex"

## remove households with 0 quarterly expenditure on food
aguas_full <- aguas ## save all data for later
aguas <- aguas[aguas$total_food > 0,]

## keep only individuals with estimated socioeconomic class 2
aguas <- aguas[aguas$est_socio ==2 ,]

## create simplified weighting factor 
aguas$factor2 <- round(aguas$factor/75)

## repeat the observations according to new weighting factor
aguas_long <- aguas[1,]
for (i in 1:nrow(aguas)){
  if (i %% 100 == 0){
    print(i)
  }
  for (j in 1:aguas$factor2[i]){
    aguas_long <- rbind(aguas_long, aguas[i,])
  }
}
aguas_long <- aguas_long[-1,]

## calculate food expense per person in thousands of pesos (adjusted for 2 percent inflation)
aguas_long$food <- 0.001*aguas_long$total_food/aguas_long$total_people*(1.02)^2

## split based on sex of main household provider
aguas_F <- aguas_long[aguas_long$sex == 2,]
aguas_M <- aguas_long[aguas_long$sex ==1,]

## remove households with more than 20000 pesos per person per month;
## this is three households with female providers, and 2 with male providers
aguas_F <- subset(aguas_F, aguas_F$food <= 20)
aguas_M <- subset(aguas_M, aguas_M$food <= 20)

## save food expenditure data for both groups
write.csv(aguas_F[c("food")], "aguas2018_food_F.csv", row.names = FALSE)
write.csv(aguas_M[c("food")], "aguas2018_food_M.csv", row.names = FALSE)

## find the median food expenditure per person in upper income household from same region
aguas_soc4 <- aguas_full[aguas_full$est_socio == 4,]

## exclude households with no quarterly food expenditure (none for this example)
aguas_soc4 <- aguas_soc4[aguas_soc4$total_food > 0,]

## calculate food expense per person in thousands of pesos (adjusted for 2 percent inflation)
aguas_soc4$food <- 0.001*aguas_soc4$total_food/aguas_soc4$total_people*(1.02)^2

food4_rep <- NULL
for (i in 1:nrow(aguas_soc4)){
  food4_rep <- c(food4_rep, rep(aguas_soc4$food[i], aguas_soc4$factor[i]))
}

## confirmation that this expense is 4.29 thousand pesos (design value for threshold)
median(food4_rep)

## we now obtain design values for alpha1, beta1, alpha2, and beta2
## first load the .csv files just created
## save the food data in vector form
y1 <- read.csv("aguas2018_food_F.csv")$food 
y2 <- read.csv("aguas2018_food_M.csv")$food

set.seed(5)
burnin <- 1000
nchains <- 1
nthin <- 2
ndraws <- nthin*100000
sum_y1 <- sum(y1)
sum_logy1 <- sum(log(y1))
n1 <- length(y1)
mu0 <- 2; tau0 <- 0.25; kappa0 <- 2; lambda0 <- 0.25
model1.fit <- jags.model(file="JAGS_gamma.txt",
                         data=list(n=n1, sum_y = sum_y1, sum_logy = sum_logy1, 
                                   tau0 = tau0, mu0 = mu0, zero = 0,
                                   kappa0 = kappa0, lambda0 = lambda0), 
                         n.chains = nchains)

update(model1.fit, burnin)
model1.samples <- coda.samples(model1.fit, c("alpha", "beta"), n.iter=ndraws, thin=nthin)

alpha.1 <- unlist(model1.samples[,1])
beta.1 <- unlist(model1.samples[,2])

set.seed(4)
sum_y2 <- sum(y2)
sum_logy2 <- sum(log(y2))
n2 <- length(y2)
model2.fit <- jags.model(file="JAGS_gamma.txt",
                         data=list(n=n2, sum_y = sum_y2, sum_logy = sum_logy2, 
                                   tau0 = tau0, mu0 = mu0, zero = 0,
                                   kappa0 = kappa0, lambda0 = lambda0), 
                         n.chains = nchains)

update(model2.fit, burnin)
model2.samples <- coda.samples(model2.fit, c("alpha", "beta"), n.iter=ndraws, thin=nthin)

alpha.2 <- unlist(model2.samples[,1])
beta.2 <- unlist(model2.samples[,2])

## find posterior means for design values
mean(alpha.1) ## should be 2.11
mean(beta.1) ## should be 0.69
mean(alpha.2) ## should be 2.43
mean(beta.2) ## should be 0.79

## save posterior draws to .csv files for reference
write.csv(alpha.1, "alpha1s_2018.csv", row.names = FALSE)
write.csv(beta.1, "beta1s_2018.csv", row.names = FALSE)
write.csv(alpha.2, "alpha2s_2018.csv", row.names = FALSE)
write.csv(beta.2, "beta2s_2018.csv", row.names = FALSE)

## now we inflate the variances of these approximately gamma
## marginal posteriors to find informative priors for the numerical study

## define factor by which to inflate the variance
var_inf <- 10

## consider alpha1
## find the posterior mode and variance
a1dens <- density(alpha.1, n = 1024, bw = "SJ")
a1_mode <- a1dens$x[which.max(a1dens$y)]
a1_var <- var(alpha.1)

## solve system of equations to get desired gamma distribution
b1_star <- (-a1_mode - sqrt(a1_mode^2 + 4*var_inf*a1_var))/(-2*var_inf*a1_var); b1_star
a1_star <- 1 + a1_mode*b1_star; a1_star
informs <- c(a1_star, b1_star)

## consider beta1
## find the posterior mode and variance
a1dens <- density(beta.1, n = 1024, bw = "SJ")
a1_mode <- a1dens$x[which.max(a1dens$y)]
a1_var <- var(beta.1)

## solve system of equations to get desired gamma distribution
b1_star <- (-a1_mode - sqrt(a1_mode^2 + 4*var_inf*a1_var))/(-2*var_inf*a1_var); b1_star
a1_star <- 1 + a1_mode*b1_star; a1_star
informs <- rbind(informs, c(a1_star, b1_star))

## consider alpha2
## find the posterior mode and variance
a1dens <- density(alpha.2, n = 1024, bw = "SJ")
a1_mode <- a1dens$x[which.max(a1dens$y)]
a1_var <- var(alpha.2)

## solve system of equations to get desired gamma distribution
b1_star <- (-a1_mode - sqrt(a1_mode^2 + 4*var_inf*a1_var))/(-2*var_inf*a1_var); b1_star
a1_star <- 1 + a1_mode*b1_star; a1_star
informs <- rbind(informs, c(a1_star, b1_star))

## consider beta2
## find the posterior mode and variance
a1dens <- density(beta.2, n = 1024, bw = "SJ")
a1_mode <- a1dens$x[which.max(a1dens$y)]
a1_var <- var(beta.2)

## solve system of equations to get desired gamma distribution
b1_star <- (-a1_mode - sqrt(a1_mode^2 + 4*var_inf*a1_var))/(-2*var_inf*a1_var); b1_star
a1_star <- 1 + a1_mode*b1_star; a1_star
informs <- rbind(informs, c(a1_star, b1_star))

## output informative gamma distribution hyperparameters to .csv file
write.csv(informs, "informs_gamma.csv", row.names = FALSE)

## repeat this process for the Weibull distribution
set.seed(6)
burnin <- 1000
nchains <- 1
nthin <- 2
ndraws <- nthin*100000
n1 <- length(y1)
mu0 <- 2; tau0 <- 1; kappa0 <- 2; lambda0 <- 1
model1.fit <- jags.model(file="JAGS_weibull.txt",
                         data=list(n=n1, y = y1, 
                                   tau0 = tau0, mu0 = mu0,
                                   kappa0 = kappa0, lambda0 = lambda0), 
                         n.chains = nchains)

update(model1.fit, burnin)
model1.samples <- coda.samples(model1.fit, c("nu", "l"), n.iter=ndraws, thin=nthin)

nu.1 <- unlist(model1.samples[,2])
lambda.1 <- unlist(model1.samples[,1])

set.seed(7)
n2 <- length(y2)
model2.fit <- jags.model(file="JAGS_weibull.txt",
                         data=list(n=n2, y = y2, 
                                   tau0 = tau0, mu0 = mu0,
                                   kappa0 = kappa0, lambda0 = lambda0), 
                         n.chains = nchains)

update(model2.fit, burnin)
model2.samples <- coda.samples(model2.fit, c("nu", "l"), n.iter=ndraws, thin=nthin)

nu.2 <- unlist(model2.samples[,2])
lambda.2 <- unlist(model2.samples[,1])

## find posterior means for design values
mean(nu.1) ## should be 1.41
mean(lambda.1) ## should be 3.39
mean(nu.2) ## should be 1.49
mean(lambda.2) ## should be 3.42

## save posterior draws to .csv files for reference
write.csv(nu.1, "nu1s_2018.csv", row.names = FALSE)
write.csv(lambda.1, "lambda1s_2018.csv", row.names = FALSE)
write.csv(nu.2, "nu2s_2018.csv", row.names = FALSE)
write.csv(lambda.2, "lambda2s_2018.csv", row.names = FALSE)

## now we inflate the variances of these approximately gamma
## marginal posteriors to find informative priors for the numerical study

## define factor by which to inflate the variance
var_inf <- 100

## consider nu1
## find the posterior mode and variance
a1dens <- density(nu.1, n = 1024, bw = "SJ")
a1_mode <- a1dens$x[which.max(a1dens$y)]
a1_var <- var(nu.1)

## solve system of equations to get desired gamma distribution
b1_star <- (-a1_mode - sqrt(a1_mode^2 + 4*var_inf*a1_var))/(-2*var_inf*a1_var); b1_star
a1_star <- 1 + a1_mode*b1_star; a1_star
informs <- c(a1_star, b1_star)

## consider lambda1
## find the posterior mode and variance
a1dens <- density(lambda.1, n = 1024, bw = "SJ")
a1_mode <- a1dens$x[which.max(a1dens$y)]
a1_var <- var(lambda.1)

## solve system of equations to get desired gamma distribution
b1_star <- (-a1_mode - sqrt(a1_mode^2 + 4*var_inf*a1_var))/(-2*var_inf*a1_var); b1_star
a1_star <- 1 + a1_mode*b1_star; a1_star
informs <- rbind(informs, c(a1_star, b1_star))

## consider nu2
## find the posterior mode and variance
a1dens <- density(nu.2, n = 1024, bw = "SJ")
a1_mode <- a1dens$x[which.max(a1dens$y)]
a1_var <- var(nu.2)

## solve system of equations to get desired gamma distribution
b1_star <- (-a1_mode - sqrt(a1_mode^2 + 4*var_inf*a1_var))/(-2*var_inf*a1_var); b1_star
a1_star <- 1 + a1_mode*b1_star; a1_star
informs <- rbind(informs, c(a1_star, b1_star))

## consider lambda2
## find the posterior mode and variance
a1dens <- density(lambda.2, n = 1024, bw = "SJ")
a1_mode <- a1dens$x[which.max(a1dens$y)]
a1_var <- var(lambda.2)

## solve system of equations to get desired gamma distribution
b1_star <- (-a1_mode - sqrt(a1_mode^2 + 4*var_inf*a1_var))/(-2*var_inf*a1_var); b1_star
a1_star <- 1 + a1_mode*b1_star; a1_star
informs <- rbind(informs, c(a1_star, b1_star))

## output informative Weibull distribution hyperparameters to .csv file
write.csv(informs, "informs_weibull.csv", row.names = FALSE)