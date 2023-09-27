# BayesianPower
Supplementary code for Fast Power Curve Approximation with Posterior Analyses manuscript 

The code files for this manuscript are divided into 3 groups.

Group 1: extract and process data from the ENIGH 2020 website and create Figure 1
- 01-food-data-2020: processes the ENIGH 2020 food data used in the main text
- 02-code-for-figure-1: reproduces Figure 1 in the main text

Group 2: extract and process data from the ENIGH 2018 website and use this data for design purposes
- 03-design-values-priors: processes the ENIGH 2018 food data and uses the relevant posteriors
                           to return design values for gamma and Weibull models along with informative priors

Group 3: code to reproduce Figure 2 to 4 in the main text
- 04-power-curve-figure-2: code to implement our power curve approximation procedure with the gamma tail
                           probability example and hypothesis tests facilitated posterior probabilities.
                           This code was also used to produce Figures D.1 and D.2 in the supplement.
- 05-power-curve-figure-3: code to implement our power curve approximation procedure with the gamma tail
                           probability example and hypothesis tests facilitated credible intervals.
- 06-power-curve-figure-4: code to implement our power curve approximation procedure with the Weibull tail
                           probability example and hypothesis tests facilitated posterior probabilities.

There are also 2 files "informs_gamma.csv" and "informs_weibull.csv" that contain the informative prior specifications.
The "JAGS_gamma.txt" and "JAGS_weibull.txt" files are used to approximate the relevant posteriors using JAGS
for confirmation analysis and prior fitting.
