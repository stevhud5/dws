import matplotlib.pyplot as plt
import numpy as np
import DWS_modules.DWS_DataTreatment_1 as DWS
import pylab as pl
import matplotlib.pyplot as plt
    
        
#--------------------------- Main program ------------------------------------   
#path='C:\Users\sdh\Documents\Facilities\DWS\Extract Parameters from data Files'
mypath = 'C:\\Users\\sdh\\Documents\\Projects\\2 DWS\\J_J\\Code_andTempfiles'
#path='C:\Users\sdh\Documents\Projects\2 µRheoSANS\P2 wlm\DWS\170926c'
#path='C:\Users\Mathias Reufer\Desktop\US\NIST\NetDrive\Production\Final Tests before Shipping\Final Tests\Micelle\Run1'


fnameBase='\\019-24A-Tscan_'
minNr=0
maxNr=2
maxCRtolerated=500
fig_dpi = 120
write_mean = True

plt.clf()
plt.figure(dpi=fig_dpi)
plt.xscale('log')
plt.yscale('linear')
plt.ylabel('ICF')
plt.axis([1e-8, 100, 0, 1.1])

ICF_array = []
for j in range(minNr,maxNr+1):
    fname=mypath+fnameBase+str(j)+'.dat'
    MTDuration,EchoDuration,wavelength,T,L,labs,n,R,lstar,CRmean,g2,CRhistory=DWS.loadDWSdata(fname) 
    LagTime=g2[:,0]
    ICF=g2[:,1]
    ICF_array.append(ICF)
    plt.plot(LagTime, ICF, 'ro')
    print(j,', ',T)

ICF_mean = np.mean(ICF_array, axis=0)
ICF_sd = np.std(ICF_array, axis=0)
ICF_grad = np.gradient(ICF_array, axis=0)

plt.plot(LagTime, ICF_mean, 'b+')
plt.show()


plt.clf()
plt.figure(dpi=fig_dpi)
plt.xscale('log')
plt.yscale('linear')
plt.ylabel('ICF_sd')
plt.axis([1e-8, 100, 0, np.amax(ICF_sd)])
#ICF_slope =  
plt.plot(LagTime, ICF_sd, 'b+')
plt.show()

plt.clf()
plt.figure(dpi=fig_dpi)
plt.xscale('log')
plt.yscale('linear')
plt.ylabel('ICF_grad')
plt.axis([1e-8, 100, np.amin(ICF_grad[1]), np.amax(ICF_grad[1])])
#ICF_slope =  
plt.plot(LagTime, ICF_grad[1], 'b+')
plt.show()

if write_mean:
    #get header to first file
    fname=mypath+fnameBase+str(minNr)+'.dat'
    with open(fname, 'r') as f_in:
        header = ''
        for index in range (21):
            header = header+f_in.readline()
        
    fname = mypath+fnameBase+str(minNr)+'_'+str(maxNr)+'mean.dat'
    with open(fname, 'w') as f_out:
        f_out.write(header)
        for i in range(len(LagTime)):
            #f_out.write(LagTime[i]+'\t'+ICF_mean[i]+'\n')
            f_out.write('%8.6G\t%8.6G\n' % (LagTime[i],ICF_mean[i]))
            
        f_out.write('\n\nCount Rate  History (kHz) \nElapsed Time (s)\tCount Rate (kHz)\n\n')
        for i in range(len(CRhistory)):
            f_out.write('%10.6f\t%10.6f\n' % (CRhistory[i,0],CRhistory[i,1]))
            pass
        
        f_out.write('\nMeasurement notes:\n\n')
        
