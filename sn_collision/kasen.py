from scipy import optimize, interpolate, integrate

import astropy.constants as const
import astropy.units as u
import numpy as np

from dust_extinction.dust import extinction_cal

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

class SNCompanionCollision(object):
    
    def __init__(self, a=1e13, M=1.4*const.M_sun.cgs.value, v=1e9, kappa=0.2):
        """
        parameters:
            a: binary separation in cm (default value: 1e13)
            M: ejecta mass in gram (default value: Chandrasekar mass)
            v: expansion velocity in cm/s (default value: 1e9)
            kappa: opacity in cm^2/gram (default value: 0.2)
        """
        self.a = a / 1e13
        self.M = M / (1.4*const.M_sun.cgs.value)
        self.v = v / 1e9
        self.kappa = kappa / 0.2
    
    def _AngularDependence(self, theta):
        """
        This function uses the parameterization of eq. (3) in Brown et al. 2012, ApJ, 749, 18
        parameters:
            theta: viewing angle in degrees
        return paramters:
            angular dependence factor
        """
        theta_rad = theta * 0.01745
        return (0.5 * np.cos(theta_rad) + 0.5) *\
            (0.14 * theta_rad * theta_rad - 0.4 * theta_rad + 1)
        
    def _DerivedEquation(self, t):
        """
        Eqs. (22) and (25) in Kasen 2010, ApJ, 708, 1025
        parameter:
            t: in units of days after explosion
        return parameters:
            luminosity: luminosity in ergs/s
            temperature: temperature in kelvin
        """
        luminosity = 1e43 * self.a * self.M**(0.25) * self.v**(1.75) * self.kappa**(-0.75) * t**(-0.5)
        temperature = 2.5e4 * self.a**(0.25) * self.kappa**(-35. / 36) * t**(-37. / 72)
        return luminosity, temperature
    
    def _PlanckFunction(self, wavelength, temperature):
        """
        parameters:
            wavelength: wavelength in Angstrom
            temperature: thermal temperature in Kelvin
        return paramters:
            flux density in ergs/s/cm^2/ster/Angstrom
        """
        factor = const.h.cgs.value * const.c.cgs.value / (wavelength * 1e-8 * const.k_B.cgs.value * temperature)
        flux = 2. * const.h.cgs.value * const.c.cgs.value * const.c.cgs.value / (wavelength * 1e-8)**5 /\
            (np.exp(factor) - 1) * 1e-8
        return flux
    
    def Spectrum(self, t, theta):
        """
        parameters:
            t: in units of days after explosion
            theta: viewing angle in degrees
        return parameters:
            a function to calculate spectrum at a given wavelength in units of ergs/s/ster/Angstrom
        """
        luminosity, temperature = self._DerivedEquation(t)
        f = self._AngularDependence(theta)
        func = lambda wavelength: f * luminosity * np.pi * self._PlanckFunction(wavelength, temperature) /\
                                  (const.sigma_sb.cgs.value * temperature**4)
        return func

class Filter(object):
    
    def __init__(self, name, wavelength, transmission):
        """
        parameters:
            name: filter name
            wavelength: wavelength in Angstrom
            tranmission: in the units of per Angstrom per photon
        """
        self.name = name
        self.wavelength = wavelength
        self.transmission = transmission
        norm = integrate.quad(self._TransmissionInterpolation(),
                              self.wavelength.min(),
                              self.wavelength.max())[0]
        self.transmission /= norm
        self.wavelength_eff = 1 /\
            integrate.quad(interpolate.interp1d(self.wavelength,
                                                self.transmission,
                                                bounds_error=False,
                                                fill_value=0),
                           self.wavelength.min(),
                           self.wavelength.max())[0]
        
    def _TransmissionInterpolation(self):
        return interpolate.interp1d(self.wavelength, 
                                    self.transmission * self.wavelength, 
                                    bounds_error=False, 
                                    fill_value=0)
    
    def TransmissionCurve(self, wv):
        return self._TransmissionInterpolation()(wv)

class Spectrum(object):
    
    def __init__(self, wavelength, flux):
        self.wavelength = wavelength
        self.flux = flux
        
    def extinction(self, EBV, RV=3.1):
        A = extinction_cal.calALambda(self.wavelength, RV, EBV)
        self.flux *= 10**(-A/2.5)
        
    def redshift(self, z):
        self.wavelength *= (1 + z)
        self.flux /= (1 + z)
        
    def SyntheticPhotometry(self, filter_curve):
        func = lambda wv: interpolate.interp1d(self.wavelength,
                                               self.flux * filter_curve.TransmissionCurve(self.wavelength),
                                               bounds_error=False,
                                               fill_value=0)(wv)
        return integrate.quad(func, self.wavelength.min(), self.wavelength.max())[0]

def ObservedMagnitude(model, theta, filter_curve, 
                      t = 1.0, redshift=0, mu=0,
                      local_EBV=0, galactic_EBV=0,
                      ):
    wv = filter_curve.wavelength / (1 + redshift)
    spec_func = model.Spectrum(t, theta)
    spec = Spectrum(wv, spec_func(wv))
    spec.extinction(local_EBV)
    spec.redshift(redshift)
    dist = 10**(mu/5 + 1) * u.Quantity(1, unit=u.pc).to(u.cm).value
    spec.flux /= (4. * np.pi * dist * dist)
    spec.extinction(galactic_EBV)
    obs_flux = spec.SyntheticPhotometry(filter_curve)  # F_lambda
    obs_flux *= filter_curve.wavelength_eff * (filter_curve.wavelength_eff * 1e-8 / const.c.cgs.value) * 1e23  # F_nu in Jy
    return -2.5 * np.log10(obs_flux / 3631)  # AB mag