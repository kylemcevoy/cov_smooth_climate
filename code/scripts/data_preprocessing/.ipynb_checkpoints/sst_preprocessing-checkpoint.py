"""

@author: Kyle McEvoy

data pre-processing of SST
sst.mnmean_1854_2024.nc file contains the NOAA Extended Reconstructed SST V5
obtained from https://psl.noaa.gov/data/gridded/data.noaa.ersst.v5.html
downloaded Nov 19, 2024 covering dates 1854-01-01 to 2024-10-01
see paper for full citation.

this file loads the SST dataset, removes the grid level monthly mean to get 
anomalies, subsets to observations with latitudes in [-60, 60] and months january
and july, and writes out to a .csv file for each month.
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

months = [month.lower() for month in calendar.month_abbr]

# path to the ERSSTv5 data
sst_path = '/home/data/ERSSTv5/monthly/sst/sst.mnmean_1854_2024.nc'
save_path = '/home/data/projects/clim_smooth/'

sst = xr.open_mfdataset(sst_path)
sst = sst['sst']
sst = sst.astype('float64')
sst = sst.sel(time=slice('1959', '2008'), lat=slice(60, -60))

sst_anom_list = []
for month in range(1, 13):
    sst_anoms = detrend_ssts(sst, month)
    sst_anom_list.append(sst_anoms)
    
    (sst_anoms.to_dataframe(name='sst_anoms').
    dropna().
    reset_index().
    to_csv(save_path + f'sst_detrend_anom_{months[month]}.csv'))

sst_anoms = xr.concat(sst_anom_list, dim='time').sortby('time')
sst_anoms_seas = sst_anoms.resample(time='QS-DEC').mean('time')
sst_anoms_jja = sst_anoms_seas.sel(time=(sst_anoms_seas['time.season'] == 'JJA'))

(sst_anoms_jja.to_dataframe(name='sst_anoms').
 dropna().
 reset_index().
 to_csv(save_path + 'sst_detrend_anom_jja.csv'))
