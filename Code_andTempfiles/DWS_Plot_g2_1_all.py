import matplotlib.pyplot as plt
import numpy as np
import DWS_modules.DWS_DataTreatment_1 as DWS
import pylab as pl
import matplotlib.pyplot as plt
    
        
#--------------------------- Main program ------------------------------------   
#path='C:\Users\sdh\Documents\Facilities\DWS\Extract Parameters from data Files'
mypath = 'C:\\Users\\sdh\\Documents\\Projects\\2 DWS\\J_J'
#path='C:\Users\sdh\Documents\Projects\2 µRheoSANS\P2 wlm\DWS\170926c'
#path='C:\Users\Mathias Reufer\Desktop\US\NIST\NetDrive\Production\Final Tests before Shipping\Final Tests\Micelle\Run1'


fnameBase='\\020-24A-Tscan_'
minNr=0
maxNr=3
maxCRtolerated=500

plt.clf()
plt.xscale('log')
plt.yscale('linear')
plt.axis([1e-8, 100, 0, 1.1])

for i in range(minNr,maxNr):
    fname=mypath+fnameBase+str(i)+'.dat'
    MTDuration,EchoDuration,wavelength,T,L,labs,n,R,lstar,CRmean,g2,CRhistory=DWS.loadDWSdata(fname) 
    LagTime=g2[:,0]
    ICF=g2[:,1]
    plt.plot(LagTime, ICF, 'ro')
    print(i)
    
plt.show()