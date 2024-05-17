## code to reproduce Figure 1 in the main text
## this file should be run after "01-food-data-2020.R"

## BEGIN SETUP ##

## load necessary packages
require(ggplot2)
require(cowplot)
require(rjags)
require(coda)

## load the .csv files created using "01-food-data-2020.R" 
## save the food data in vector form
y1 <- read.csv("aguas_food_F.csv")$food 
y2 <- read.csv("aguas_food_M.csv")$food

## save the data as data frame (for plotting later)
foodF <- read.csv("aguas_food_F.csv")
foodM <- read.csv("aguas_food_M.csv")

## use JAGS to obtain posterior draws via MCMC
set.seed(1)
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

set.seed(2)
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

## obtain posterior of theta
threshold <- 4.82
theta.1 <- 1 - pgamma(threshold, alpha.1, beta.1)
theta.2 <- 1 - pgamma(threshold, alpha.2, beta.2)
theta <- theta.1/theta.2

## find good histogram breaks based on the data
breaks1 <- hist(y1, plot=FALSE)$breaks
breaks2 <- hist(y2, plot=FALSE)$breaks
binwidth2 <-0.64
bins2 <- seq(0, 19.2, binwidth2)

## create data frames for plotting
x_fit1 <- seq(floor(qgamma(0.001,mean(alpha.1),mean(beta.1))),
              ceiling(qgamma(0.999,mean(alpha.1),mean(beta.1))), by = 0.01)
fit1 <- data.frame(xfit = x_fit1, yfit = dgamma(x_fit1,mean(alpha.1),mean(beta.1)))

x_fit2 <- seq(floor(qgamma(0.001,mean(alpha.2),mean(beta.2))),
              ceiling(qgamma(0.999,mean(alpha.2),mean(beta.2))), by = 0.01)
fit2 <- data.frame(xfit = x_fit2, yfit = dgamma(x_fit2,mean(alpha.2),mean(beta.2)))

hist_lower2 <- 0
hist_upper2 <- 20

group1_name <- "Female Household Provider"
n1 <- length(y1)

use_pctl <- 1
metric_plot <- paste("Tail Probability", sep="")
metric_value1 <- round(mean(ifelse(foodF$food >= 4.82,1,0)),3)
percentile <- 0.9

## make the plot for the female provider group data (top left corner)
plot1a <-
  ggplot(data=foodF, aes(foodF$food)) + theme_bw() +
  geom_histogram(aes(y = (stat(count) / sum(count))/binwidth2), breaks = bins2,
                 col="#1B86B4", 
                 fill="#077DAA", 
                 alpha = 0, size = 1) + 
  coord_cartesian(xlim = c(hist_lower2, hist_upper2)) +
  labs(title=paste(group1_name, "\n")) +
  labs(x=bquote(atop(' ', atop(textstyle('Food Expenditure per Person (MXN $1000)'),
                               textstyle('n = '*.(n1)*", "*.(metric_plot)*' ('*hat(theta)[1]*  
                                           ') = '*.(metric_value1))))), y="Density\n") + 
  theme(plot.title = element_text(hjust = 0.5,size=16,face="bold")) + 
  theme(plot.subtitle = element_text(hjust = 0.5,size=14)) +
  geom_area(data = fit1, aes(x=xfit, y=yfit), fill="#94D3FF", col="#077DAA", alpha=0.45, size = 1) +
  geom_segment(aes(x = 4.82, y = 0, 
                   xend=4.82, yend = (Inf*(use_pctl))), color="grey16", size=1.5) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) + ylim(0, 0.3) +
  annotate(geom="text", x=8, y=0.28, label=expression(kappa*' = 4.82'),
           color="black", size = 6)

group2_name <- "Male Household Provider"
n2 <- length(y2)
metric_value2 <- round(mean(ifelse(foodM$food >= 4.82,1,0)),3)

## make the plot for the male provider group (bottom left corner)
plot1c <- ggplot(data=foodM, aes(foodM$food)) + theme_bw() +
  geom_histogram(aes(y = (stat(count) / sum(count))/binwidth2), breaks = bins2,
                 col="#9E8203", 
                 fill="#FED758", 
                 alpha = 0, size = 1) + 
  coord_cartesian(xlim = c(hist_lower2, hist_upper2)) +
  labs(title=paste(group2_name,"\n")) +
  labs(x=bquote(atop(' ', atop(textstyle('Food Expenditure per Person (MXN $1000)'),
                               textstyle('n = '*.(n2)*", "*.(metric_plot)*' ('*hat(theta)[2]*  
                                           ') = '*.(metric_value2))))), y="Density\n") + 
  theme(plot.title = element_text(hjust = 0.5,size=16,face="bold")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  geom_area(data = fit2, aes(x=xfit, y=yfit), fill="#FED758", col="#9E8203", alpha=0.45, size = 1) + 
  geom_segment(aes(x = 4.82, y = 0, 
                   xend=4.82, yend = (Inf*(use_pctl))), color="grey16", size=1.5) + ylim(0, 0.3) +
  annotate(geom="text", x=8, y=0.28, label=expression(kappa*' = 4.82'),
           color="black", size = 6)

fig1.row1 <- plot_grid(plot1a + theme(plot.margin=unit(c(0.25,0.25,0.45,0.25),"cm")), 
                       plot1c + theme(plot.margin=unit(c(0.45,0.25,0.25,0.25),"cm")),
                       rel_widths = c(1, 1))

## output as .pdf file for the article
pdf(file = "Figure1BAText.pdf",   # The directory you want to save the file in
    width = 9, # The width of the plot in inches (12.41)
    height = 4.85) # The height of the plot in inches (10.7)

fig1.row1

dev.off()