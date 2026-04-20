import numpy as np
import pandas as pd
import xarray as xr

data_path = '/home/data/ERA5/month/t2m/t2m.nc'
save_path = '/home/data/projects/clim_smooth/'

# time bounds
year_start = '1959'
year_end = '2008'

month = 7
# latitude bounds for the box containing the Great Plains region
north, south = 45, 32.5
# longitude bounds for the box (in 0 to 360)
west, east = 255, 275

t2m = xr.open_mfdataset(data_path)
t2m = t2m['t2m']

t2m_sub = t2m.sel(latitude=slice(north, south),
                  longitude=slice(west, east),
                  time=slice(year_start, year_end))

t2m_sub_jul = t2m_sub.sel(time=t2m_sub['time.month'] == month)

weights = np.cos(np.deg2rad(t2m_sub_jul.latitude))
t2m_sub_jul_weighted = t2m_sub_jul.weighted(weights)
t2m_GP_mean = t2m_sub_jul_weighted.mean(dim=('latitude', 'longitude'))

t2m_GP_mean_pd = t2m_GP_mean.to_pandas()
t2m_GP_mean_pd.to_csv(save_path + 't2m_GP_mean_jul.csv')

t2m_sub_JJA = t2m_sub.resample(time='QS-DEC').mean('time')
t2m_sub_JJA = t2m_sub_JJA.sel(time=t2m_sub_JJA['time.season'] == 'JJA')

weights = np.cos(np.deg2rad(t2m_sub_JJA.latitude))
t2m_sub_JJA_weighted = t2m_sub_JJA.weighted(weights)
t2m_GP_mean_JJA = t2m_sub_JJA_weighted.mean(dim=('latitude', 'longitude'))

t2m_GP_mean_JJA_pd = t2m_GP_mean_JJA.to_pandas()
t2m_GP_mean_JJA_pd.to_csv(save_path + 't2m_GP_mean_jja.csv')