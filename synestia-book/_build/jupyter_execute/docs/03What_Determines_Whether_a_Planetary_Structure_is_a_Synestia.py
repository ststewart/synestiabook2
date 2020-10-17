#!/usr/bin/env python
# coding: utf-8

# Before reading this Jupyter Notebook, it may be helpful to review [angular velocity](./Angular_Velocity.ipynb).

# # What determines whether a planetary object is a synestia?
# Imagine that you and a friend are sitting on the edge of a playground merry-go-round. You are sitting closer to the center and are strapped in place, while your friend is at the very edge, free of restraints. If you join hands and hold on to each other while a third person spins the merry-go-round faster and faster, you'll reach a point where your hands cannot hold on. Your friend will be subject to their own speed of rotation about the merry-go-round as they fly off. This is not physically applicable, but you can use it as an analogy to help wrap your head around what happens when a planetary object transforms into a synestia.
# 
# ## Corotation Limit: The Boundary Between Planet and Synestia
# A synestia is a planetary body (planet, for example) that is rotating so rapidly that its outer portion cannot rotate at the same rate and moves at an angular velocity dictated by the gravity and gas (pressure) fields. The body exceeds a threshold called the <i>corotation limit</i>, the threshold above which planets turn into synestias.
# 
# Beyond the corotation limit, the body can no longer rotate together at the same angular velocity (<i>corotate</i>). Part of the body (located roughly near the body's equator) assumes a flared shape and will rotate slower than the interior, which will continue to corotate rapidly. This rapid rotation causes oblateness, or a bulge along the equator, in the planet-like region. The corotation limit is dependent on a number of factors, namely, the mass of the body, the angular momentum of the body (which determines the body's rate of rotation), and the structure of the body's interior [dependent on temperature and composition of the body (e.g. rock versus iron)]. Because it is dependent on so many parameters that can be specific to a single body, the corotation limit is hard to define via an equation. However, we can develop some intuition for how these different parameters affect the corotation limit. Let's add one more piece to your understanding of the corotation limit: the threshold above which planets turn into synestias.
# 
# ## Imagining the Corotation Limit as a Geosynchronous Orbit
# Consider a planet which has a moon stably orbiting it (the orbit of the moon is Keplerian). If this moon were also in a <i>geosynchronous</i> orbit, then the time it takes the moon to orbit the planet is the same as the time it takes the planet to complete a rotation. When a planet's equatorial radius is equal to the orbital radius of its geosynchronous moon, the part of the planet at the location of the geosynchronous moon can no longer orbit at the planet's rotational rate. This is similar to when a planet has exceeded the corotation limit, and the outer portion of the planet (anything beyond the orbital radius of the geosynchronous satellite) can no longer corotate with the interior of the planet.
# 
# ### The Dependence of the Corotation Limit on Angular Momentum
# If the planet rotated at a faster rate, the equator would bulge, as if the planet had been squashed down at the poles. The faster the spin rate of a planet, the more oblate it is.
# 
# For a planet with a shorter length of day, the shorter amount of time a geosynchronous moon has to complete an orbit. If we are comparing planet-moon systems with the same angular momentum (one with a longer day, one with a shorter day), then for a moon to maintain geosynchronous orbit around a planet with a shorter day, the moon must travel a shorter distance during the planet's day. This moon orbits at a smaller radius, closer to the center of the planet.
# 
# If a planet spins so fast that its oblate equator reaches where its geosynchronous moon tightly orbits, then what? The moon cannot keep up with the rotation of the planet. If the planet were to spin any faster, then even the equator of the planet could not corotate with the rest of the planet! That is the essence of the boundary between planet and synestia.
# 
# So, think of it like this: if a planet's equatorial radius exceeds the orbital radius of a geosynchronous satellite, then the planet morphs into a synestia. Please see the interactive below. Pay attention to how the planet's equator (at z = 0 km on orange body) and moon (red dot) move in relation to one another as the planet's rotation speeds up (length of day decreases) and reaches the corotation limit (body turns cyan).

# ```{margin} Running a code cell
# Access interactive features by 'Launch CoLab' or 'Launch Binder' from the rocket logo at the top of the page. When the interactive environment is ready, place your cursor in the code cell and press shift-return to execute the code. CoLab launches more quickly.
# ```
# Click the + symbol to see the code that generates the next interactive feature.

# In[1]:


import math
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import warnings

G = 6.67408e-11 #mks

def Maclaurin(ecc):
    return 2.*(3. - 2.*ecc**2)*np.arcsin(ecc)*((1. - ecc**2)**0.5)/(ecc**3) - 6.*(1. - ecc**2)/(ecc**2) #- (omega**2)/(math.pi*G*rho)

def semi_axes(e, rho, M_Earth):
    # 0 <= e < 1
    # rho >= 0
    radius_avg_cubed = 3.*M_Earth/(4.*math.pi*rho) #m^3
    c = ((1. - e**2)*radius_avg_cubed)**(1./3.)
    a = c/((1. - e**2)**0.5)
    return c,a

def synestia_shape(T,rho,M_Earth):
    warnings.catch_warnings()
    warnings.simplefilter("ignore")
    omega = 2.*math.pi/(T*3600.) #s^-1 #convert period to seconds
    r_geosynch_km = (G*M_Earth/omega**2)**(1./3.)/1e3
    radius_avg_cubed = 3.*M_Earth/(4.*math.pi*rho) #m^3
    omega_kep_avg = (G*M_Earth/radius_avg_cubed)**0.5
    Macla_ratio = (omega**2)/(math.pi*G*rho)

    ecc_interp = np.linspace(0.,0.952887,100)
    omg_ratio = Maclaurin(ecc_interp)
    f = interp1d(omg_ratio, ecc_interp)
    ecc_planet = f(Macla_ratio)
    if Macla_ratio > 0.449331:
        ecc_planet = 0.952887 #max. eccentricity possible for stable planet

    ecc_test = np.linspace(0.,1.,100)
    y_test = Maclaurin(ecc_test)

    n = 50
    theta = np.linspace(0.,math.pi,n)
    phi = np.arange(0.,math.pi,math.pi/n)
    c,a = semi_axes(ecc_planet, rho, M_Earth)
    v_a_km = omega*a/1e3
    omega_kep_a = (G*M_Earth/a**3)**0.5
    v_kep_a_km = omega_kep_a*a/1e3

    THETA, PHI = np.meshgrid(theta, phi)
    x = a*np.sin(THETA)*np.cos(PHI)
    z = c*np.cos(THETA)

    ind_x_R = np.where(np.amax(x))
    edge_ellip_x = x[ind_x_R]
    edge_ellip_z = z[ind_x_R]
    circ_x = (radius_avg_cubed**(1./3.))*np.linspace(-1,1,100)
    circ_y = (((radius_avg_cubed**(1./3.))**2 - circ_x**2)**0.5)
    min_x = -1.6*(radius_avg_cubed**(1./3.))/1e3
    max_x = 2.*(radius_avg_cubed**(1./3.))/1e3

    plt.figure(figsize=(16,4))
    if (omega > omega_kep_a):
        plt.fill_between(edge_ellip_x[0]/1e3,-edge_ellip_z[0]/1e3,edge_ellip_z[0]/1e3,color='cyan')
        plt.fill_between(-edge_ellip_x[0]/1e3,-edge_ellip_z[0]/1e3,edge_ellip_z[0]/1e3,color='cyan')
        plt.annotate('Planet has become\na synestia!',(0,0))
    else:
        plt.annotate('Still\na planet.',(0,0))
        plt.fill_between(edge_ellip_x[0]/1e3,-edge_ellip_z[0]/1e3,edge_ellip_z[0]/1e3,color='orange')
        plt.fill_between(-edge_ellip_x[0]/1e3,-edge_ellip_z[0]/1e3,edge_ellip_z[0]/1e3,color='orange')
    plt.plot(circ_x/1e3,circ_y/1e3,'--',color='grey')
    plt.plot(circ_x/1e3,-circ_y/1e3,'--',color='grey')
    plt.plot(r_geosynch_km,0,'ro')
    plt.annotate('  Planet\n  equatorial\n  linear\n  velocity\n  {0:.1f} km/s\n'.format(v_a_km),(a/1e3,0))
    plt.xlabel('x (km)')
    plt.ylabel('z (km)')
    plt.axis('scaled')
    if (r_geosynch_km < max_x):
        plt.xlim([min_x,max_x])
    else:
        plt.xlim(xmin=min_x)
    plt.show()
    
    print('For reference, {0:.1f} km/s would be the Keplerian linear velocity at the equator of the planet.'.format(v_kep_a_km))   
    
from ipywidgets import interactive,FloatLogSlider,IntSlider
style = {'description_width': 'initial'}
layout = {'width': '400px'}
interactive_plot = interactive(synestia_shape,
         T=IntSlider(value=24, min=1, max=24, step=1, description='Length of Day for Planet (hr)',
                        continuous_update=True, readout=True, readout_format='.0f', style=style, layout=layout),
         rho=IntSlider(value=3300., min=1000., max=3700., step=100., description='Density of Planet (kg/m$^3$)',
                      continuous_update=True, readout=True, readout_format='.0f', style=style, layout=layout),
         M_Earth=FloatLogSlider(value=5.97219e24, base=10, min=22, max=27, step=1, description='Mass of Planet (kg)',
                      continuous_update=False, readout=True, readout_format='.1e', style=style, layout=layout)
        )
output = interactive_plot.children[-1]
interactive_plot


# <i>Caption</i>. Earth-mass, uniformly-dense silicate body. Many of the finer details about a body's interior structure are complicated, so this calculation is simplified and assumes we have a "homogeneous" spherical body [same material (no separation of a mantle and core), density, and temperature throughout]. For reference, Earth's present-day mass is 6.0 x 10$^{24}$ kg, and the average density of a rock (silicate) under standard conditions is 3300 kg/m$^3$.

# ## Weighing Dependencies of the Corotation Limit on Planet Density and Mass
# In addition to increasing the rotation rate of a body, decreasing the density of a material also helps a body exceed the corotation limit. Decreasing the density of the material can indicate that the material is hotter. It is easier to exceed the corotation limit when the body is warmer.
# 
# However, it is easiest to cross the threshold between planet and synestia when decreasing the length of day for the planet (the time it takes for a planet to complete one rotation) as opposed to changing the planet's mass or density. The corotation limit is most sensitive to the rate of rotation of a body. The faster a planet spins, the more likely it will transition into a synestia. This is why it is easier for planets with greater angular momentum, particularly those involved in giant impacts, to become synestias.
