#%%
import numpy as np
import astropy.units as u
from astropy import constants as const
from astropy.cosmology import z_at_value
from astropy.cosmology import Planck15 as cosmo

#%%

#interval=100*u.Myr

class AstroObj(object):
    """
    A initial class that can evolved in time

    Args:
    mass (float): mass in Msun
    rate (float): rate in Msun/yr 
    z    (float): redshift of the source
    """
    def __init__(self, mass, rate, z):
        self.mass = mass*u.Msun
        self.rate = rate*u.Msun/u.yr
        self.z    = z

    def next_step_in_time(self, t_i):
        interval = (t_i - cosmo.age(self.z)).to(u.Myr)
        newmass = self.mass + self.rate * interval
        self.mass = newmass
        self.rate = self.rate
        newz = cosmo.age(self.z) + interval
        if cosmo.age(0)==newz:
            self.z=0
        else:
            self.z=z_at_value(cosmo.age, newz).value
        return


# %%
blackhole=AstroObj(mass=1e8, rate=1, z=0.5)
# %%
interval=100*u.Myr
age_obs,age_z0=cosmo.age(blackhole.z).to(u.Myr).value,cosmo.age(0).to(u.Myr).value
time=np.arange(age_obs,age_z0,interval.value)
time[-1]=cosmo.age(0).to(u.Myr).value
# %%
print(np.round(age_obs,2),np.round(blackhole.z,1),np.round(np.log10(blackhole.mass.value),1))
for binning in time:
    blackhole.next_step_in_time(binning*u.Myr)
    print(np.round(binning,2),np.round(blackhole.z,1),np.round(np.log10(blackhole.mass.value),1))


# %%
