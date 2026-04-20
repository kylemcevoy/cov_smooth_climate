"""
data pre-processing of SST
sst.mnmean_1854_2024.nc file contains the NOAA Extended Reconstructed SST V5
obtained from https://psl.noaa.gov/data/gridded/data.noaa.ersst.v5.html
downloaded Nov 19, 2024 covering dates 1854-01-01 to 2024-10-01
see paper for full citation.

this file loads the SST dataset, subsets to July, and removes the grid level 
monthly mean to get anomalies, subsets to observations with latitudes in 
[-60, 60], and writes out to a .csv file for each month.
"""

import xarray as xr
import calendar

def detrend_ssts(sst, month):
    sst_month = sst.sel(time=(sst['time.month'] == month))
    
    # remove a linear trend, including a mean term, from the ssts at each location
    sst_poly_coef = sst_month.polyfit(dim='time', deg=1)['polyfit_coefficients']
    sst_trend = xr.polyval(sst_month.time, sst_poly_coef)
    sst_anoms = sst_month - sst_trend
    
    return sst_anoms

# July SST detrended anoms.
month = 7
months = [month.lower() for month in calendar.month_abbr]

# path to the ERSSTv5 data
sst_path = '/home/data/ERSSTv5/monthly/sst/sst.mnmean_1854_2024.nc'
save_path = '/home/data/projects/clim_smooth/'

sst = xr.open_mfdataset(sst_path)
sst = sst['sst']
sst = sst.astype('float64')
sst = sst.sel(time=slice('1959', '2008'), lat=slice(60, -60))

sst_anoms = detrend_ssts(sst, month)
    
(sst_anoms.to_dataframe(name='sst_anoms').
    dropna().
    reset_index().
    to_csv(save_path + f'sst_detrend_anom_{months[month]}.csv'))
