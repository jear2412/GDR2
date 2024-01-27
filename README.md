# GDR2
Generalized Decomposition Priors on R2

_Abstract:_ The adoption of continuous shrinkage priors in high-dimensional linear models has gained momentum, driven by their theoretical and practical advantages. One of these shrinkage priors is the R2D2 prior, which comes with intuitive hyperparameters and well understood theoretical properties.  
The core idea is to specify a prior on the percentage of explained variance $R^2$ and to conduct a Dirichlet decomposition to distribute the explained variance among all the regression terms of the model.  Due to the properties of the Dirichlet distribution, the competition among variance components tends to gravitate towards negative dependence structures, fully determined by the individual components' means. Yet, in reality, specific coefficients or groups may compete differently for the total variability than the Dirichlet would allow for. In this work we address this limitation by proposing a generalization of the R2D2 prior, which we term the \textit{Generalized Decomposition R2} (GDR2) prior. 

Our new prior provides great flexibility in expressing dependency structures as well as enhanced shrinkage properties. Specifically, we explore the capabilities of variance decomposition via logistic normal distributions.
Through extensive simulations and real-world case studies, we demonstrate that GDR2 priors yield strongly improved out-of-sample predictive performance and parameter recovery compared to R2D2 priors with similar hyper-parameter choices.

_Keywords:_ Bayesian inference, prior specification, shrinkage priors, variance decomposition, regularization

_Citation:_  Aguilar, J.E, Bürkner, P. C. (2024). Generalized Decomposition Priors on R2, [arXiv preprin](https://arxiv.org/abs/2401.10180)
