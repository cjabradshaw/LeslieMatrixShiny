# LeslieMatrixShiny

## Preamble
This app projects a user-defined Leslie (age-classified) matrix to examine population changes through time. This fully customisable app includes a density-feedback function on survival relative to desired initial population size and carrying capacity, stochastic projections with user-defined variance for survival, a generationally scaled catastrophe function, a single 'pulse' disturbance function, and a 'press' disturbance function with a user-defined time application.

A detailed instructions tab (tab H) is included for guidance, but a brief sequence description is included below. User- defined settings in each tab are carried over to subsequent tabs.

<ol type="A">
  <li>SET-UP: set matrix dimensions (longevity), age x -specific survival (<em>s</em><sub>x</sub>) and fertility (<em>f</em><sub>x</sub>) probabilities, offspring sex ratio, % variance around survival/fertility probabilties, and whether lifespan is abrupt or diffuse.</li>

<li>MATRIX PROPERTIES : shows Leslie matrix according to settings in tab A, as well as the dominant eigen value <em>λ</em> instaneous rate of population change <em>r</em> generation length <em>G</em>, and reproduction number <em>R</em><sub>0</sub> (number of ♀ offspring/adult ♀).</li>

<li>DENSITY FEEDBACK : set initial population size and carrying capacity <em>K</em>, as well as the three coefficients (<em>a</em>, <em>b</em>, <em>c</em>) from a logistic power function to define the relationship between a survival modifier <em>S</em><sub>mod</sub> and population size.</li>

<li>PROJECT : deterministic projection of the population, setting the number of years (or generations) to project the population, initial population size, and whether to invoke the density-feedback function set in the previous tab.</li>

<li>STOCHASTIC : stochastic projection of the population based on previous settings (including the % variances set in the first tab); the user can set the number of iterations to repeat the stochastic resampling, the quasi-extinction threshold (population size below which it is considered functionally extinct), and whether to invoke a generationally scaled catatastrophic mortality probability (the magnitude and variance of which can be set by the user).</li>

<li>SINGLE PULSE : a single 'pulse' disturbance, where the user can set the disturbance to be either a proportion of the total population that is removed, or a fixed number of individuals removed, at the time (year) the user wishes to invoke the pulse.</li>

<li>PRESS : a press disturbance, where the user can set the disturbance to be either a proportion of the total population that is removed, or a fixed number of individuals removed, during the interval over which the user wishes to invoke the press.</li></ol>

This  Github repository provides all the 'under-the-bonnet' R code for the app.
