# IPLGP

In this package, we provide a tool to select parental lines for multiple traits in plant breeding. To confirm the potential of a set of parental lines, this package provides two useful tools. One is *D-score*, which is a new criteria we proposed in our study for determine the potential of a set of parental line. The other is the simulation test of the breeding process. There are three strategy in our simulation test, which is (i) GEBV-O considers only genomic estimated breeding values (GEBVs) of the candidate individuals; (ii) GD-O considers only genomic diversity (GD) of the candidate individuals; and (iii) GEBV-GD considers both GEBV and GD.   
  
## Installation
  
IPLGP can be installed form GitHub by the following command:  
```install_github
# library(devtools)  
install_github("py-chung/IPLGP", dependencies = TRUE, force = TRUE)
```
  
And IPLGP can be installed form CRAN by the following command:
```install.packages
install.packages("IPLGP")
```
  
## Main functions
  
+ `GA.Dscore()` Fonction for getting a set with highest D-score by genetic algorithm. 
+ `GBLUP.fit()` Fonction for getting the fitting values of a set of individuals by GBLUP.
+ `geno.d()` Fonction for getting the design matrix and kinship matrix of dominance effects.
+ `phe.sd()` Fonction for getting the standardize phenotypic values.
+ `simu.gamete()` Fonction for simulating the genotype of a gamete.
+ `simu.GDO()` Fonction for simulating the progeny with GD-O strategy.
+ `simu.GEBVO()` Fonction for simulating the progeny with GEBV-O strategy.
+ `simu.GEBVGD()` Fonction for simulating the progeny with GEBV-GD strategy.
+ `output.best()` Fonction for obtaining the GEBV average curves and the summary statistics form the output of simulation test.
+ `output.gain()` Fonction for obtaining the genetic gain average for each target trait.
  
More information can be seen in the following file:  
[Package ‘IPLGP’](https://cran.r-project.org/web/packages/IPLGP/IPLGP.pdf)
  
## Citing this package
  
For more imformation about our method, please check our published article:  
+ Chung PY, Liao CT. 2020. Identification of superior parental lines for biparental crossing via genomic prediction. PLoS ONE 15(12):e0243159. doi:[10.1371/journal.pone.0243159](https://journals.plos.org/plosone/article/authors?id=10.1371/journal.pone.0243159)
  
If you use IPLGP in your research, we would appreciate your citation of the study article.
