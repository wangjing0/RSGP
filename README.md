##  Reward Sensitive Gaussian Process (RSGP)
<p align="center">
  <img src="RSGP.png" height="400" >
</p>
### Algorithm for RSGP
<p align="center">
  <img src="algo1.png" height="400" >
</p>
### Algorithm for MCMC
<p align="center">
  <img src="algo2.png" height="400" >
</p>

### Algorithm for Gradient search
<p align="center">
  <img src="algo3.png" height="400" >
</p>


  1) simulate data from the model 
  2) feed the data back to recover model parameters 

  	RSGPsim.m


## Auxiliary Functions

Reward sensitive covariance matrix

  	covSEisoRew.m
    covScaleRew.m

## References
Gaussian process for machine learning (GPML) toolkit

http://www.gaussianprocess.org/gpml

ShadedErrorBar

https://github.com/raacampbell/shadedErrorBar
