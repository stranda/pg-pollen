# Spatiotemporal model for fossil pollen counts
The R code in this repo runs a spatiotemporal multinomial model to estimate tree pollen relative abundances across eastern North America over the last 21,000 years (i.e., since the Last Glacial Maximum, or LGM). Additionally, we use prediction output from this model to calculate biotic velocity (i.e., migration rate) of trees since the LGM. Here is the order in which the files should be run, along with descriptions.

### 1. pollen_migrated_ST.R
* Read in fossil pollen data, which is already downloaded from neotomadb.org
* Data cleaning - remove undesired taxa, locations, ages
* Construct time bins, allowing pollen data to be binned into discrete time groupings
* Restructure data into arrays

### 2. estimate_pars_ST_*.R
* Read data from pollen_migrated_ST.R
* Specify model parameters, priors, and initial conditions before running model

### 3. model_assess_ST.R
* Read model output from estimate_pars_ST.R
* Create trace plots to check parameter estimation for tau^2, theta1, theta2, beta, and rho
* Create maps of estimated relative proportions (i.e., pis, calculated from eta parameter)

### 4. build_grid.R
* Create regular grid across spatial extent of study, on which predictions will be made in the next step

### 5. predict_ST_*.R
* Read model output from estimate_pars_ST_*.R and grid from build_grid.R
* Run prediction function

### 6. plot_predictions.R
* Read output from predict_ST.R
* Create maps of relative abundance predictions across time and taxa

### 7. biotic_velocity.R
* Read output from predict_ST.R
* Resample and mask prediction output to match spatial extent and resolution of interest
* Run and plot output of biotic velocity function

### 8. prep_pollen_for_ABC.R
* Read output from predict_ST.R
* Interpolate pollen predictions to 30-year time resolution
* Mask out ice and water from pollen prediction rasters
* Weight relative abundance predictions by proportion of cell covered by land (for cells with partial glacial coverage)
* Convert to rasterstacks
