# Forcings
#CO2 data is from ftp://data.iac.ethz.ch/CMIP6/input4MIPs/UoM/GHGConc/CMIP/yr/atmos/UoM-CMIP-1-1-0/GHGConc/gr3-GMNHSH/v20160701/mole_fraction_of_carbon_dioxide_in_air_input4MIPs_GHGConcentrations_CMIP_UoM-CMIP-1-1-0_gr3-GMNHSH_0000-2014.csv
#The data is saved as [File 4- Global and hemispheric means with annual resolution.csv] in "./forcings/data"

#Volcanic data is from https://www.dropbox.com/sh/fm5588buvfbf0x8/AACoMfme2uhygDdvO5o0cSy8a?dl=0
#and ftp://iacftp.ethz.ch/pub_read/luo/CMIP6/ECHAM6/
#EVA(eVolv2k) reconstructed zonal mean AOD (mid-visible, i.e., 550 nm), covering the 500 BCE to 1900 Ce time period:
#For 1900 (or 1850) to present, CMIP_ECHAM6_T63_radiation_AOD.nc is used to fill in the forcing table
#volcanic data at varying locations is weighted by cosine of their corresponding latitudes. 
#detailed computation of the volcanic forcing can be found in Global_mean_volcanic.R


#The solar forcing data in the table is computed from from SATIRE-H (Holocene) at http://www2.mps.mpg.de/projects/sun-climate/data.html
#Irradiance from SATIRE-H (Holocene) is recorded on a decadal basis from 9495BC - 1939AD and daily from 1940AD onwards. Hence, the data prior to 1939AD is interpolated in a splined manner and after 1940AD, a annual mean is computed. For consistence of resolution of the data, annual data from 1640AD to 2000AD is also spline interpolated. 
#solar forcing is saved as TSI.rds and TSI.csv. A moving average for TSI is also computed and saved as TSI_ma.csv and TSI_ma.rds
