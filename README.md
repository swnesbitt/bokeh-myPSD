# myPSD
## A t-matrix based interactive tool for understanding weather radar polarimetric variables in rain

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/swnesbitt/bokeh-myPSD/master?urlpath=/proxy/5006/bokeh-app)

This repository runs the myPSD (bokeh)[https://docs.bokeh.org/en/latest/] app, which enables interactive exploration of simulation of radar polarimetric variables at weather radar wavelengths varying the parameters of a 3-parameter normalized Gamma rain particle size distributions (e.g., Testud et al. 2001) and the standard deviation of the droplet canting angle.  Under the hood, the tool uses `pytmatrix` to forward simulate the radar variables. Thanks to Jussi Leinonen of MeteoSwiss for developing [pytmatrix](https://github.com/jleinonen/pytmatrix)!
   
Click on the Binder button above to run the `binder` app!

## Contributions

Developed by [Steve Nesbitt](https://swnesbitt.com), Professor at the University of Illinois Urbana-Champaign.

If you'd like to make this tool better, feel free to create a pull request!
