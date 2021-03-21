# ORCA025.L75-GJM189
This repository hold the code, as well as namelist files for running the NEMO global configuration named ORCA025.L75-GJM189
## REFERENCE CODE
This configuration is based on the NEMO code, version 3.5beta than can be downloaded using

```
svn co -r 4521 http://forge.ipsl.jussieu.fr/nemo/svn/trunk NEMO_3.5beta
``` 
With respect to this revision of NEMO, some specific developments were added. The directory [NEMOGCM](./NEMOGCM) holds the
exact version of the NEMO code and specific developments used for this configuration.

## BRIEF DESCRIPTION
### Overview
  * This simulation uses the ORCA025.L75 NEMO model configuration : a global 1/4 degree configuration with 75 vertical levels. 
    * The grid-cell size varied from about 28 km at the equator down to 12 km in the Arctic and even down to less than 7 km in the Southern most part of the model.
    * Vertical grid spacing is 1m near the surface, increasing to 60m at 500m and up to 200m below 4000m.
    * time step for model integration is set to 1080 second, although it may have been reduced for short period, in order to fix unstabilities.
  * This simulation cover the period 1958-2017.
  * It uses the DFS5.2 atmospheric forcing set.

### Parameterizations:
 1. use filtered free surface (key_dynspg_flt)
 2. use vector form advection scheme for dynamics, with corrections.
 3. use TVD advection scheme for tracers.
 4. use biharmonic viscosity scaled with the cube of the mesh size.
 5. use laplacian isopycnal diffusivity for tracers.
 6. use TKE vertical mixing parameterization with enhanced vertical diffusion for deep convection. use tidal mixing parameterization.
 7. use LIM2 ice model
 8. use BBL (bottom boundary layer) parameterization.
 9. use free slip lateral condition except in the Idonesian through-flow area, Mediterannean Sea, West coast of Greenland ( near cape desolation) and in the northern part of Nares Strait.

### Forcing:
  1. Atmospheric forcing is DFS5.2, with ```CORE bulk formulae``` and relative winds (ocean surface velocity components taken into account for the wind stress computation). 
  2. SSS restoring toward Levitus98 with a piston velocity of 167 mm/day ( 60 days/10 meters).
  3. Run-off from Dai-Trenberth including climatological iceberg contribution from Da Silva
 

### Run time files:
   Most of the run time files are indicated in the namelist files, except for:
   
   * bathymetry : ```ORCA025_bathy_etopo1_gebco1_smoothed_coast_corrected_bering_may11.nc```
   * coordinates : ```coordinates_ORCA_R025_lombok+ombai_v2.nc```
   * bottom friction : ```orca025_bfr_coef_GRD100.nc```
   * Ice initialisation : ```Init_Ice_GLORYS1V1_NSIDC_BOOTSTRAP_y1989m01_new.nc```
   * AABW damping mask : ```ORCA025.L75_dmp_mask.nc```

Data files as well as forcing file are available at [MEOM-IGE OpenDap Serveur](https://ige-meom-opendap.univ-grenoble-alpes.fr/thredds/catalog/meomopendap/extract/ORCA025.L75/ORCA025.L75-GJM189/catalog.html).


