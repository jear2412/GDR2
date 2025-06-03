# Generalized Decomposition Priors on R2


_Abstract:_ The adoption of continuous shrinkage priors in high-dimensional linear models has gained widespread attention due to their practical and theoretical advantages. Among them, the R2D2 prior has gained popularity for its intuitive specification of the proportion of explained variance (R2) and its theoretically grounded properties. The R2D2 prior allocates variance among regression terms through a Dirichlet decomposition. However, this approach inherently limits the dependency structure among variance components to the negative dependence modeled by the Dirichlet distribution, which is fully determined by the mean. This limitation hinders the prior’s ability to capture more nuanced or positive dependency patterns that may arise in real-world data.
To address this, we propose the Generalized Decomposition R2 (GDR2) prior, which replaces the Dirichlet decomposition with the more flexible Logistic-Normal distribution and its variants. By allowing richer dependency structures, the GDR2 prior accommodates more realistic and adaptable competition among variance components, enhancing the expressiveness and applicability of R2
-based priors in practice. Through simulations and real-world benchmarks, we demonstrate that the GDR2 prior improves out-of-sample predictive performance and parameter recovery compared to the R2D2 prior. Our framework bridges the gap between flexibility in variance decomposition and practical implementation, advancing the utility of shrinkage priors in complex regression settings.

_Keywords:_ Bayesian inference, prior specification, shrinkage priors, variance decomposition, regularization


How to use this repo? 

The main point of this repo is that the reader gets familiarized with how to use the stan code for the Generalized Decomposition prior. The reader can take a look at example.html for this.

_Citation:_  Aguilar, J.E, Bürkner, P. C. (2025).  [Generalized Decomposition Priors on R2](https://projecteuclid.org/journals/bayesian-analysis/advance-publication/Generalized-Decomposition-Priors-on-R2/10.1214/25-BA1524.full). *Bayesian Analysis*. Available for free also on [arXiv](https://arxiv.org/abs/2401.10180). doi:10.1214/25-BA1524

_Funding Statement:_
Funded by Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany’s Excellence Strategy EXC 2075-390740016 and DFG Project 500663361. We acknowledge the support by the Stuttgart Center for Simulation Science (SimTech). We acknowledge the computing time provided on the Linux HPC cluster at Technical University Dortmund (LiDO3), partially funded in the course of the Large-Scale Equipment Initiative by DFG Project 271512359.

_Acknowledgments:_ We would like to thank Aki Vehtari and his research group at Aalto University for their valuable comments and insightful discussions on an earlier version of this work. In particular, we are grateful to David Kohns for his detailed and constructive feedback.

_Contact_ javier.aguilarr at icloud.com
