import scipy
import scipy.optimize   
import scipy.special
import itertools,numpy as np
import pandas as pd 
from scipy.signal import savgol_filter
    
  
def loadDWSdata(path):
    with open(path) as t_in:
        header_=np.genfromtxt(itertools.islice(t_in, 17),skip_header=6,delimiter='\t',usecols=(1))
    with open(path) as t_in:
        data_=np.genfromtxt(itertools.islice(t_in, 3000),skip_header=21,delimiter='\t',comments="Count Rate  History (kHz)", skip_footer=4)   
    nanPos_=pd.isnull(data_).any(1).nonzero()[0]     # get all positions of NaN 
    g2_=data_[0:min(nanPos_)]               # save valid part of g2 (ignore NaN)
    CRhistory_=data_[max(nanPos_)+1:]       # save valid part of CR history (ignore NaN)
    return (header_[0],header_[1],header_[2],header_[3],header_[4],header_[5],header_[7],header_[8],header_[9],header_[10],g2_,CRhistory_)
# MTDuration,EchoDuration,wavelength,T,L,labs,n,R,lstar,CR,g2,CRhistory=loadDWSdata(fname)
    

def loadDWSdataMR(path):
    with open(path) as t_in:
        header_=np.genfromtxt(itertools.islice(t_in, 16),skip_header=5,delimiter='\t',usecols=(1))
    data_=np.loadtxt(path,skiprows=21,delimiter='\t',)  
    return (header_[0],header_[1],header_[2],header_[3],header_[4],header_[5],header_[7],header_[8],header_[10],header_[9],data_)
# Transmission:     MTDuration,EchoDuration,wavelength,T,L,labs,n,R,lstar,CR,dataMRA=loadDWSdataMR(fnameMR)
# Backscattering:   MTDuration,EchoDuration,wavelength,T,Gamma,labs,n,R,lstar,Cr,dataMRA=loadDWSdataMR(fnameMR)


def writeDWSdata(path,header_,data_):
    pass
