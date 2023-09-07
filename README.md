# BayesianPower
Supplementary code for Practicable Power Curve Approximation for Bayesian Equivalence Tests manuscript 

The code files for this manuscript are divided into 3 groups.

Group 1: extract and process data from the ENIGH 2020 website and create Figure 1
- 01-food-data-2020: processes the ENIGH 2020 food data used in the main text
- 02-code-for-figure-1: reproduces Figure 1 in the main text

Group 2: extract and process data from the ENIGH 2018 website and use this data for design purposes
- 03-design-values-priors: processes the ENIGH 2018 food data and uses the relevant posteriors
                           to return design values for gamma and Weibull models along with informative priors

Group 3: code to reproduce Figure 3 to 5 in the main text
- 04-power-curve-figure-3: code to implement our power curve approximation procedure with the gamma tail
                           probability example and equivalence tests facilitated posterior probabilities.
                           This code was also used to produce Figures C.1 and C.2 in the supplement.
- 05-power-curve-figure-4: code to implement our power curve approximation procedure with the gamma tail
                           probability example and equivalence tests facilitated credible intervals.
- 06-power-curve-figure-5: code to implement our power curve approximation procedure with the Weibull tail
                           probability example and equivalence tests facilitated posterior probabilities.
