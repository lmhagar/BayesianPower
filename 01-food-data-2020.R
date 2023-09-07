## code to prepare ENIGH 2020 food data set

## BEGIN SETUP ##

## Process to obtain data file from ENIGH website
## Step 1: Go to this website: https://www.inegi.org.mx/programas/enigh/nc/2020/#Datos_abiertos
## Step 2: Click the "Download CSV" button to download a .zip file
## Step 3: Open .zip file
## Step 4: Open the "conjunto_de_datos_concentradohogar_enigh_2020_ns" subfolder
## Step 5: Open the "conjunto_de_datos" subfolder
## Step 6: Extract the lone .csv file: "conjunto_de_datos_concentradohogar_enigh_2020_ns.csv" and move it to the working directory

## read data set
enigh <- read.csv("conjunto_de_datos_concentradohogar_enigh_2020_ns.csv")

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

## calculate food expense per person in thousands of pesos
aguas_long$food <- 0.001*aguas_long$total_food/aguas_long$total_people

## split based on sex of main household provider
aguas_F <- aguas_long[aguas_long$sex == 2,]
aguas_M <- aguas_long[aguas_long$sex ==1,]

## remove households with more than 20000 pesos per person per month;
## this is three households with female providers, and 2 with male providers
aguas_F <- subset(aguas_F, aguas_F$food <= 20)
aguas_M <- subset(aguas_M, aguas_M$food <= 20)

## save food expenditure data for both groups

write.csv(aguas_F[c("food")], "aguas_food_F.csv", row.names = FALSE)
write.csv(aguas_M[c("food")], "aguas_food_M.csv", row.names = FALSE)

## find the median food expenditure per person in upper income household from same region
aguas_soc4 <- aguas_full[aguas_full$est_socio == 4,]

## exclude households with no quarterly food expenditure (none for this example)
aguas_soc4 <- aguas_soc4[aguas_soc4$total_food > 0,]

## calculate food expense per person in thousands of pesos
aguas_soc4$food <- 0.001*aguas_soc4$total_food/aguas_soc4$total_people

food4_rep <- NULL
for (i in 1:nrow(aguas_soc4)){
  food4_rep <- c(food4_rep, rep(aguas_soc4$food[i], aguas_soc4$factor[i]))
}

## confirmation that this expense is 4.82 thousand pesos
median(food4_rep)