import numpy as np
import LORAparameters as LORA

def find_event_time(info):
    
    earliest=9e9
    earliest_sec=9e9
    earliest_nsec=9e9
    first_det=0
    
    for i in np.arange(LORA.nLORA):
        if  info[i]['gps']>2:
            if (info[i]['gps'])<=earliest_sec:
                earliest_sec=info[i]['gps']


    for i in np.arange(LORA.nLORA):
        if  info[i]['gps']>2:
            if (info[i]['gps'])==earliest_sec:
                if info[i]['nsec']<=earliest_nsec:
                    earliest_nsec=info[i]['nsec']
                    first_det=i+1

    return first_det,earliest_sec,earliest_nsec



def find_peaks(detectors):
    peaks=np.zeros([LORA.nLORA])
    peak_args=np.zeros([LORA.nLORA])

    for i in np.arange(LORA.nLORA):
        peaks[i]=np.max(detectors[i].counts)
        peak_args[i]=np.argmax(detectors[i].counts)

    return peaks, peak_args
