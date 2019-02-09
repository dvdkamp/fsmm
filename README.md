# fsmm: Fuel Stick Moisture Model
Simulates moisture levels and temperature for standardized fuel sticks

This repository provides code and example forcing data for the model described by: 

[van der Kamp, D.W., Moore, R.D., and McKendry, I.G. (2017) A model for simulating the moisture
content of standardized fuel sticks of various sizes. Agriculture and Forest Meteorology. 236: 123-
134.](https://doi.org/10.1016/j.agrformet.2017.01.013)

For those without access you can download the [final manuscript](https://drive.google.com/file/d/1ILxx_vfGJNLtrObvQFuGoauspsppLOb9/view?usp=sharing).

It also provides example fuel moisture observations that were used to optimize and evaluate the model



## Dependencies
The model requires the package [`deSolve`](https://cran.r-project.org/web/packages/deSolve/index.html) to be installed:

```R
install.packages("deSolve")
require(deSolve)
```

## Fortran is required

This core of the model is written in Fortran and called by R using the `deSolve` package. The Shared Object/DLL is generated with `R CMD SHLIB` and is loaded with `dyn.load`

As such you will be required to have Fortran installed. 

## Shortwave partitioning and downwelling longwave radiation not yet implemented. 

I have yet to implement these two features for this version of the model. The model will be updated soon with these additional features. See section A.2 of van der kamp *et al* (2017) for details on these two features. 

## Running the model

Running `run.fsmm()` will run the model with default arguments and produce modelled fuel moisture and temperature. 

#### Required arguments for `run.fsmm`

  * `fsmm.model.input`: a data.frame containing the required forcing variables. The default is to load the example forcing data.frame: `data/fsmm.example.forcing.vars.csv`. This dataframe has to have the following column names:
    * `Alt_s`: Sun Altitude (rad)
    * `Azi_s`: Sun Azimuth. (0 rad at noon, Westerly directions: 0-pi rad Easterly directions: pi to 2pi rad)
      * note: For now, the user will need to calculate these sun position variables themselves. the packages `oce` and `suncalc` can help here.
    * `k.down.dir`: Downwelling direct shortwave radiation (Wm<sup>-2</sup>) 
    * `k.down.diff`: Downwelling diffuse shortwave radiation (Wm<sup>-2</sup>) 
    * `temp`: Air temperature (<sup>o</sup>C)
    * `L_dn`: Downwelling longwave radiation (Wm<sup>-2</sup>) 
    * `rh`: Relative humidity (%)
    * `wind.speed`: Wind speed (m/s)
    * `precip`: Precipitation (mm)
  * `fsmm.model.pars`: a vector storing the model parameter values. The default is the optimal parameter values for site "BC1" listed in Table 1 of van der Kamp *et al.* (2017). The parameters need to be provided **in the following order:**
    * `A` : constant for EMC equation
    * `B` : constant for EMC equation
    * `K_s`: Diffusivity (m^2 s^-1)
    * `M_MAX`: Maximum allowable moisture content (fraction)
    * `DR_O`: volume fraction of the stick taken up by the outer layer
  * `RHO_S`: stick density (kg m-3). default value is 400 kg m-3.
  * `r`: radius of stick (m). Default value is 0.0065 m.
  * `L`: length of stick (m). Deafult value is 0.41 m.
  * `m_i`: initial moisture content of stick (weight fraction). Default value is 0.1.

#### Output

A matrix with seven columns:
* `time` Model time (sec)       
* `fuel.temp.mod.outer` Temperature of the outer layer (<sup>o</sup>C)
* `fuel.mois.mod.outer`: Moisture of the outer layer (weight fraction)
* `fuel.mois.mod.core`: Moisture of the stick core (weight fraction)
* `fuel.temp.mod.core`: Temperature of the stick core (<sup>o</sup>C)
* `fuel.mois.mod`: Average moisture of the entire stick (weight fraction)     
* `fuel.temp.mod`: Average temperature of the entire stick (<sup>o</sup>C)

## Example fuel moisture observations

The csv file: `data/fsmm.example.obs.fuel.moisture.csv` contains hourly observations of fuel moisture for two sites: sites "BC1" and "BC2" mentioned in van der Kamp *et al* (2017)
