# Running myPSDexplorer on binder

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/swnesbitt/bokeh-myPSD/master?urlpath=/proxy/5006/bokeh-app)

This repository runs the myPSD bokeh app, which enables interactive exploration of simulation of radar polarimetric variables at weather radar wavelengths varying the parameters of a 3-parameter normalized Gamma rain particle size distributions (e.g., Testud et al. 2001) and the standard deviation of the droplet canting angle.  Under the hood, the tool uses `pytmatrix` to forward simulate the radar variables. Thanks to Jussi Leinonen of MeteoSwiss for developing [pytmatrix](https://github.com/jleinonen/pytmatrix)!
   
When people click on the Binder button above, they should be directed to the running Bokeh app.

## Contributions

Developed by [Steve Nesbitt](https://swnesbitt.com), Professor at the University of Illinois Urbana-Champaign.

If you'd like to make this tool better, feel free to create a pull request!
