# A-method-to-estimate-actual-infrastructure-induced-mortality-by-integrating-sampling-biases

This repository contains data files and code for the publication *A method to estimate actual infrastructure-induced mortality by integrating sampling biases*.

This study estimates the total number of roadkills occurring on a road by accounting for the three main biases present in carcass surveys: carcass-location bias, carcass-persistence bias, and carcass-observation bias.

Installation & Reproducibility
Here we use **renv** to manage package dependencies. To replicate the exact computational environment:

1. Clone or download this repository.
2. Open the project in RStudio (`Supplementary_Material.Rproj` file).
3. Run the following command in the R console to install all necessary libraries:
```r
renv::restore()
```

ABSTRACT
1. Human infrastructures are among the most impactful wildlife threats. Although estimates of animal mortality by these structures exist over a given period, they typically do not account for several detection biases (i.e., difference between recorded and true mortality). Consequently, true mortality rates may be severely underestimated, as well as their impact on populations and species.
   
2. We present a hierarchical Bayesian latent-state modelling framework that sequentially accounts for three main processes that produce biases in estimating mortality abundance: the probability that a hit animal dies on the surveyed area (carcass location probability), the probability that the carcass remains on the surveyed area until the survey is conducted (carcass persistence probability), and the probability that the carcass is observed during the survey process (carcass observation probability). We employ a comprehensive simulation study where we test the effects of variability in species characteristics, sampling design, latent-state parameters, and prior information on the ability of our model to estimate mortality abundance on roads as total number of roadkills. We then apply our framework on a case study to estimate the total number of roadkills per km in Mediterranean ecosystems while evaluating the cross-efficiency of different sampling methods.
   
3.  Our framework accurately recovers the total number of roadkills from simulated census data for most simulation scenarios. We detected the highest disagreement between modelling outcomes and simulated data when variability in simulated carcass persistence probability was high and Bayesian priors were highly diffuse. In the case study, our results show notably high roadkill numbers (e.g., estimating 48.92 per km passerines based on 8.04 observed counts), along with substantial variation across different vertebrate groups. Furthermore, our case study confirms that walking and cycling surveys outperform driving surveys in carcass observation rate and provide complementary information between them, observing partially distinct sets of species and carcass sizes.
   
4. Our modelling framework offers an efficient approach to estimate mortality rates for a wide range of taxa. Optimizing application requires extensive fieldwork for bias estimation and integration. We provide a checklist to help managers to assess when infrastructure-related mortality can be assessed most robustly to prioritize conservation efforts.
