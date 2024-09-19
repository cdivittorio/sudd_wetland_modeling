# sudd_wetland_modeling
This repository contains input data and code that can be used to model the Sudd Wetland hydrology (South Sudan, Nile River Basin), from 2000 to 2015. 

For information on the satellite-derived dynamic flood extents, please refer to this paper: https://doi.org/10.1016/j.rse.2017.11.001<br /> 
For information on the hydrologic model and data inputs, please refer to this paper: https://doi.org/10.1016/j.ejrh.2021.100922

For access to the actual flood maps in geotiff or netcdf form, please contact Courtney Di Vittorio at divittoc@wfu.edu

## Description of files:

### *SuddHydrologicData.xlsx* - contains all data inputs to run the lumped hydrologic balance model, including
	Flooded Area
		1. Connected Flooded Area (km(<sup>2</sup>)) - The monthly net connected flooded area of the Sudd derived from MODIS satellite imagery. The connectivity algorithm identifies flooded pixels that are physically connected to the main wetland water body. This is the flooded area used in the publication
		2. Full Flooded Area (km(<sup>2</sup>)) - The monthly net flooded area of the Sudd derived from MODIS satellite imagery, before applying the connectivity algorithm.
	Averaged Monthly Precipitation (mm), averaged over the Sudd connected flooded area, from different sources described in the hydrologic model publication, including:
		1. CRU
		2. TRMM 3B43
		3. PERSIANN
		4. ARCv2
	Averaged Monthly Potential Evapotranspiration (PET) (mm), averaged over the Sudd connected flooded area, from different sources described in the hydrologic model publication, including:
		1. Sutcliffe and Parks Climatology
		2. CRU (Penmen Monteith)
		3. Hargreaves (using CRU data)
  	Estimated Inflows (bcm/month)
   	Estimated Outflows (bcm/month)
        
      
### *hydrologicModel.m* - MATLAB code that runs the calibrated hydrologic model. Instructions are provided within code.
