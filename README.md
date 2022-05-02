##  Reward Sensitive Gaussian Process (RSGP)

[Article](https://elifesciences.org/articles/55872)
<p align="center">
  <img src="RSGP.png" height="500" >
</p>


<image src="https://github.com/wangjing0/RSGP/blob/master/Algo.pdf" width="1000"  />

### Algorithm for RSGP
<p align="center">
  <img src="algo1.png" width="1000" >
</p>

### Algorithm for Gradient search
<p align="center">
  <img src="algo2.png" width="1000" >
</p>

### Algorithm for MCMC
<p align="center">
  <img src="algo3.png" width="1000" >
</p>

## Code and Auxiliary Functions
  1) simulate data from the model 
  2) feed the data back to recover model parameters 

  	RSGPsim.m


Reward sensitive covariance matrix

  	covSEisoRew.m
    covScaleRew.m

## References
Gaussian process for machine learning (GPML) toolkit

http://www.gaussianprocess.org/gpml

ShadedErrorBar

https://github.com/raacampbell/shadedErrorBar
