# Data Generation for Causal Inference Simulations
Generates synthetic data to apply simulations for causal inference.

The basic model used in this function is a partially linear regression model extensions: 

![img](http://latex.codecogs.com/svg.latex?Y%3D%5Ctheta_%7B0%7DD%2Bg_%7B0%7D%28X%29%2BU%2C%5C%5C%0D%0AD%3Dm_%7B0%7D%28X%29%2BV%2C%5C%5C%0D%0A%5Ctheta_%7B0%7D%3Dt_%7B0%7D%28Z%29%2BW%0D%0A)

The data generating process creates data of the following form.
Note that all variables are randomly generated which is why the distribution might slightly change every time.

![Data Distribution](https://github.com/QuantLet/Data_Generation/blob/master/DataGen_Distribution_Plot_different_theta.png)
