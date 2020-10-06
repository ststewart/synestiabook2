#!/usr/bin/env python
# coding: utf-8

# In[3]:


#move this notebook to folder above syndef to run
from syndef import synfits #import synestia snapshot (impact database)
import numpy as np
import matplotlib.pyplot as plt

test_rxy=np.linspace(7e6,60e6,100) #m
test_z=np.linspace(0.001e6,30e6,50) #m
rxy=np.log10(test_rxy/1e6) #Mm log10
z=np.log10(test_z/1e6) #Mm log10
TESTRXY,TESTZ=np.meshgrid(test_rxy,test_z) #2-D grid of rxy, z for color plot
#y=np.zeros(np.shape(rxy)) #array of zeros for residual fit

#rho1=synfits.resfuncspl(synfits.SNAP_Canup.rhomidfit[1],rxy,y)
#rho2=synfits.resfuncspl(synfits.SNAP_CukStewart.rhomidfit[1],rxy,y)
#rho3=synfits.resfuncspl(synfits.SNAP_Quintana.rhomidfit[1],rxy,y)

snaprho1=np.log10(synfits.SNAP_Canup.rho[synfits.SNAP_Canup.ind_outer])
snaprho2=np.log10(synfits.SNAP_CukStewart.rho[synfits.SNAP_CukStewart.ind_outer])
snaprho3=np.log10(synfits.SNAP_Quintana.rho[synfits.SNAP_Quintana.ind_outer])

snaprhomid1=np.log10(synfits.SNAP_Canup.rho[synfits.SNAP_Canup.ind_outer_mid])
snaprhomid2=np.log10(synfits.SNAP_CukStewart.rho[synfits.SNAP_CukStewart.ind_outer_mid])
snaprhomid3=np.log10(synfits.SNAP_Quintana.rho[synfits.SNAP_Quintana.ind_outer_mid])

snaprxy1=np.log10(synfits.SNAP_Canup.rxy[synfits.SNAP_Canup.ind_outer]/1e6)
snaprxy2=np.log10(synfits.SNAP_CukStewart.rxy[synfits.SNAP_CukStewart.ind_outer]/1e6)
snaprxy3=np.log10(synfits.SNAP_Quintana.rxy[synfits.SNAP_Quintana.ind_outer]/1e6)

snaprxymid1=np.log10(synfits.SNAP_Canup.rxy[synfits.SNAP_Canup.ind_outer_mid]/1e6)
snaprxymid2=np.log10(synfits.SNAP_CukStewart.rxy[synfits.SNAP_CukStewart.ind_outer_mid]/1e6)
snaprxymid3=np.log10(synfits.SNAP_Quintana.rxy[synfits.SNAP_Quintana.ind_outer_mid]/1e6)

snapz1=np.log10(synfits.SNAP_Canup.z[synfits.SNAP_Canup.ind_outer]/1e6)
snapz2=np.log10(synfits.SNAP_CukStewart.z[synfits.SNAP_CukStewart.ind_outer]/1e6)
snapz3=np.log10(synfits.SNAP_Quintana.z[synfits.SNAP_Quintana.ind_outer]/1e6)

const1=10.5#10 to 11; 10.55 (fiducial)
const2=0.86#0.85 to 0.9; 0.86 (fiducial)
const3=1e38 #0.9e35 (fiducial) / 1.5e33 (underestimate) / 1.1e37 (cross) / 1e38 (overestimate)
const4=-5.1 #-4.7 (fiducial) / -4.5 (underestimate) / -5 (cross) / -5.1 (overestimate)
test_z_s=const1*np.power(TESTRXY,const2) #scale height fit in m
test_rho_g=const3*np.power(TESTRXY,const4)*np.exp(-np.power(TESTZ/test_z_s,2))
test_rho_gmid=const3*np.power(test_rxy,const4)

plt.figure(figsize=(16,5))
plt.subplot(131)
#plt.plot(rxy,rho1,'b')
plt.plot(snaprxymid1,snaprhomid1,'r.')
plt.plot(np.log10(test_rxy/1e6),np.log10(test_rho_gmid),'k')
plt.xlabel('log r$_{xy}$ (Mm)')
plt.ylabel('log midplane density (kg/m$^3$)')
plt.title('Canup')
plt.xlim([.8,2])
plt.ylim([-2,3])

plt.subplot(132)
#plt.plot(rxy,rho2,'b')
plt.plot(snaprxymid2,snaprhomid2,'r.')
plt.plot(np.log10(test_rxy/1e6),np.log10(test_rho_gmid),'k')
plt.xlabel('log r$_{xy}$ (Mm)')
plt.ylabel('log midplane density (kg/m$^3$)')
plt.title('Cuk and Stewart')
plt.xlim([.8,2])
plt.ylim([-2,3])

plt.subplot(133)
#plt.plot(rxy,rho3,'b')
plt.plot(snaprxymid3,snaprhomid3,'r.')
plt.plot(np.log10(test_rxy/1e6),np.log10(test_rho_gmid),'k')
plt.xlabel('log r$_{xy}$ (Mm)')
plt.ylabel('log midplane density (kg/m$^3$)')
plt.title('Quintana')
plt.xlim([.8,2])
plt.ylim([-2,3])
plt.show()
plt.close()

plt.figure(figsize=(16,5))
plt.subplot(131)
plt.pcolor(np.log10(TESTRXY/1e6),np.log10(TESTZ/1e6),np.log10(test_rho_g))
plt.scatter(snaprxy1,snapz1,c=snaprho1)
plt.xlabel('log r$_{xy}$ (Mm)')
plt.ylabel('log z (Mm)')
plt.colorbar(label='log density (kg/m$^3$)')
plt.xlim([.8,2])

plt.subplot(132)
plt.pcolor(np.log10(TESTRXY/1e6),np.log10(TESTZ/1e6),np.log10(test_rho_g))
plt.scatter(snaprxy2,snapz2,c=snaprho2)
plt.xlabel('log r$_{xy}$ (Mm)')
plt.ylabel('log z (Mm)')
plt.colorbar(label='log density (kg/m$^3$)')
plt.xlim([.8,2])

plt.subplot(133)
plt.pcolor(np.log10(TESTRXY/1e6),np.log10(TESTZ/1e6),np.log10(test_rho_g))
plt.scatter(snaprxy3,snapz3,c=snaprho3)
plt.xlabel('log r$_{xy}$ (Mm)')
plt.ylabel('log z (Mm)')
plt.colorbar(label='log density (kg/m$^3$)')
plt.xlim([.8,2])
plt.show()
plt.close()


# In[ ]:




