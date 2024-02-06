This repository contains all the code to reproduce the results in [Bayesian Causal Inference with Gaussian Process Networks](https://arxiv.org/abs/2402.00623).
- /A_Thaliana contains the data, code and results for the application to the A. Thaliana dataset.
- /Local_approx contains the implementation and simulation of the local approximation.
- /Simulation_plots contains all generated plots from the simulation study.

The .Stan files contain code for the Gaussian (Gauss.stan) and additive GP (Add.stan) models.
- BayesStanFns.R and GPscore.R contain functions to compute the score according to the GP and GP^2 models respectively.
- Fourier_fns.R contains functions to generate non-linear data for the simulations.
- LinGaussian_estimate.R generates the simulation results for the linear-Gaussian model for causal inference from observational data.
- MC_estimate.R generates the simulation results for causal inference with the GPN model.
- Fxsampling_fns.R contains the functions to implement Bayesian causal inference with the GPN model.

Reference
---------

```
@misc{giudice2024bayesian,
      title={Bayesian Causal Inference with Gaussian Process Networks}, 
      author={Enrico Giudice and Jack Kuipers and Giusi Moffa},
      year={2024},
      eprint={2402.00623},
      archivePrefix={arXiv},
      primaryClass={stat.ML}
}
```
