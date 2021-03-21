# Historic of ORCA025.L75-GJM189 run

## Two phases in the run :
### Initial work in 2014
 This configuration was set and simulation was performed in 2014, covering the period 1958-2012, using DFS5.2 forcing set.
This initial integration was using `dimg` output, and used the namelist parameters for controling the output.

### Completion of the run in 2019:
 Years 2013-2017 were produced after DFS5.2 was upgraded for recent years. For this last leg in the integration, although we
keep exactly the same physical   parameters and code, we use the XIOS server for the model output, and no more the `dimg`output.  
This is the reason why, for this configuration both namelistio (for `dimg`) and iodef.xml (for XIOS) files are provided.
In all cases the results were archived in the same netcdf format.

