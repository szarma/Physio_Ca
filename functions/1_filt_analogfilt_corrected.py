import numpy as np

import matplotlib.pyplot as plt

from scipy.signal import butter, sosfiltfilt, ricker, cwt

Tzac=1000  # zacetek vizualizacije v izhodnih slikah, v sekundah
Tkon=1600  # konec vizualizacije v izhodnih slikah, v sekundah



data=np.loadtxt("data.txt")
time   = data[:-1,0]  #za time uporabi celoten 0-ti stolpec minus zadnji zapis
signal = data[:-1,1:len(data[0])] #za signal pa uporabi vse podatke (-zadnji zapis), razen podatke v prvem stolpcu

# sampling=np.loadtxt("sampling.txt")
sampling = 1./np.diff(time).mean()
Tzac1=int(Tzac*sampling)
Tkon1=int(Tkon*sampling)

class Filter(object):

    """docstring for Filter."""
    
    sampling=np.loadtxt("sampling.txt")

    #def __init__(self, file, cell, fs=sampling):
    def __init__(self, signal, cell, fs=sampling):

        #self.data = np.loadtxt(file).transpose()[cell]
        self.data = signal[:,cell]

        self.nyq = 0.5*fs

        self.time = np.arange(0, (np.size(self.data-1))/fs, 1/fs)


    def bandpass(self, lowcut, highcut, order=5):

        low = lowcut / self.nyq

        high = highcut / self.nyq

        sos = butter(order, [low, high], analog=False, btype='band', output='sos')

        y = sosfiltfilt(sos, self.data)

        return y



    def lowpass(self, cutoff, order=5):

        normal_cutoff = cutoff / self.nyq

        sos = butter(order, normal_cutoff, btype='low', analog=False, output='sos')

        y = sosfiltfilt(sos, self.data)

        return y



    def highpass(self, cutoff, order=1):

        normal_cutoff = cutoff / self.nyq

        sos = butter(order, normal_cutoff, btype='high', analog=False, output='sos')

        y = sosfiltfilt(sos, self.data)

        return y



    def plot(self, *args):

        plt.clf()

        plt.plot(self.time, self.data)

        for arg in args:

            plt.plot(self.time, arg)

        plt.show()



#filter = Filter("data1.txt", 0)
filt_sig=np.zeros((len(signal[0]),len(signal)))
for rep in range(len(signal[0])):
    filter = Filter(signal,rep)
    
    #fast_component = filter.bandpass(0.03, 0.49)
    fast_component = filter.highpass(0.045)
    
    filt_sig[rep]=fast_component

    #slow_component = filter.lowpass(0.01)
    #print len(fast_component)
    
    #filter.plot(slow_component)
    print rep
    filter.plot(fast_component)

    plt.subplot(311)
    plt.plot(time,signal[:,rep],linewidth=0.5,color='dimgrey')
    plt.title(r'trace $ %ld $'%(float(rep)/1.0),fontsize=8)
    plt.rc('font', size=6)
    plt.subplot(312)
    miny=np.min(signal[(Tzac1):(Tkon1),rep])
    maxy=np.max(signal[(Tzac1):(Tkon1),rep])
    plt.plot(time,signal[:,rep],linewidth=0.67,color='dimgrey')
    plt.xlim([Tzac,Tkon])
    plt.ylim([miny-0.05,maxy+0.05])
    plt.subplot(313)
    miny=np.min(fast_component[Tzac1:Tkon1])
    maxy=np.max(fast_component[Tzac1:Tkon1])
    plt.plot(time,fast_component,linewidth=0.67,color='dimgrey')
    plt.xlim([Tzac,Tkon])
    plt.ylim([miny*1.05-0.05,maxy*1.05+0.05])
    plt.xlabel('time (s)', fontsize=6)
    plt.savefig("filt_traces\\filt_trace_%d_new.jpg"%(rep+0),dpi=200,bbox_inches = 'tight')
    plt.close()  

print len(filt_sig),len(filt_sig[0])
filt_sig1=filt_sig.transpose()
datoteka=open("DATAfilt_new.txt","w+")
for i in range(len(filt_sig1)):
    for j in range(len(filt_sig1[0]-1)):
        #print>>datoteka,filt_sig1[j][i],  
        print>>datoteka,filt_sig1[i][j],
    #print>>datoteka,filt_sig1[len(filt_sig1[0])-1][i]   
    print>>datoteka,filt_sig1[i][len(filt_sig1[0])-1]  
datoteka.close()

