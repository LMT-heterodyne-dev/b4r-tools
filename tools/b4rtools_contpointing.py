import xarray as xr
import numpy as np
from astropy.stats import sigma_clipped_stats
import fmflow as fm
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.stats import sigma_clip

def findSkyLevel(data_sci):
    Psky = fm.array(data_sci['array'].copy().values/(data_sci['integtime'].copy().values)[:,None]).mean('ch')
    date = data_sci['date'].copy().values
    t_sci = np.array([s[:-4] for s in date], 'datetime64[us]')
    flag1 = t_sci < t_sci[-1]
    #flag2 = data_sci['bufpos'] == 'ON'

    flag = flag1 #& flag2

    model_sky = models.Polynomial1D(2)
    model_sky.c0 = 0.0
    model_sky.c1 = -1.0
    model_sky.c2 = 1.0
    pfit = fitting.LinearLSQFitter()
    opfit = fitting.FittingWithOutlierRemoval(pfit, sigma_clip,niter=15, sigma=2.0)
    #fitted_sky = pfit(model_sky, t_sci[flag], Psky[flag])
    fitted_sky,filtered_data = opfit(model_sky, t_sci[flag]-t_sci[0], Psky[flag])
    #opfit = fitting.FittingWithOutlierRemoval(pfit, sigma_clip,niter=3, sigma=3.0)
    #fitted_sky,filtered_data = opfit(model_sky,t_sci[flag][~filtered_data],Psky[flag][~filtered_data])
    #print(filtered_data)
    #print(fitted_sky.c1,fitted_sky.c0)

    return Psky, fitted_sky, flag, t_sci, filtered_data


def main(path_sci,path_cal,path_ant):
    data_sci = xr.open_dataset(path_sci)
    data_cal = xr.open_dataset(path_cal)
    data_ant = xr.open_dataset(path_ant)

    Psky,fitted_sky,flag, t_sci, filtered_data = findSkyLevel(data_sci)


    #Psky = findSkyLevel(data_sci)



    plt.close()
    plt.plot(t_sci[flag]-t_sci[flag][0],Psky[flag],'o')
    plt.plot(t_sci[flag][filtered_data]-t_sci[0],Psky[flag][filtered_data],'o')
    plt.plot(t_sci[flag]-t_sci[0],fitted_sky(t_sci[flag]-t_sci[0]))
    plt.show()

    # TSYS

    # R-sky













##
