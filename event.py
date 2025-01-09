import numpy as np
import LORAparameters as LORA
import detector as det
import math
from matplotlib import pyplot as plt
import scipy.interpolate as ip
from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
from scipy.interpolate import Rbf
import scipy.stats as stats
import process_func as prf
from scipy.optimize import curve_fit
from scipy.special import gamma
from scipy.optimize import minimize




atm_flag=0
atm_l=11
f_int=np.zeros([atm_l])
f_logN_Ref=np.zeros([atm_l])
err1=np.zeros([atm_l])
f_X0=np.zeros([atm_l])
err2=np.zeros([atm_l])
f_lamb=np.zeros([atm_l])
err3=np.zeros([atm_l])

class Event:
    def __init__(self, name,type):
        self.name = name  # LOFAR event ID
        self.type = type  # V1 or V2
    
    theta=0
    elevation=0
    phi=0
    fit_theta=0
    fit_elevation=0
    fit_phi=0
    fit_theta_err=0
    fit_phi_err=0
    x_core=0
    y_core=0
    x_core_err=0
    y_core_err=0
    z_core=0
    UTC_min=0
    nsec_min=0
    energy=0
    energy_err=0
    Rm=0
    fit_elevation_err=0
    fit_phi_err=0
    Ne=0
    Ne_err=0
    CorCoef_xy=0
    Ne_RefA=0
    NeErr_RefA=0
    Energy_RefA=0
    EnergyErr_RefA=0
    direction_flag=0
    event_flag=0


def read_attenuation():
    try:
        data=np.genfromtxt(LORA.atm_file)
        return data
    except:
        print('problem reading atm file')
        print('test23')
        return 0
 
def find_counts(detector):

    if detector.number<=20:
       
    
    
        # find background
        background=detector.counts[0:int(LORA.BG_No_Bin)]
        background_mean=np.average(background)
        background_rms=np.std(background)
        detector.trace_rms=background_rms
        detector.trace_mean=background_mean
        detector.peak=np.max(detector.counts)
    
        if(background_mean<300 and background_mean>0):
            detector.corrected_threshold=detector.threshold-background_mean
        
            #print 'finding threshold from real background: {0}  {1}  {2}'.format(detector.corrected_threshold,detector.threshold,background_mean)
        
        else:
            detector.corrected_threshold=detector.threshold-detector.sec_mean
            #print 'finding threshold from second background: {0}  {1}  {2}'.format(detector.corrected_threshold,detector.threshold,detector.sec_mean)
        
        if detector.corrected_threshold<0:
            detector.corrected_threshold=detector.threshold-detector.sec_mean
        '''
        if detector.corrected_threshold<0:
            #print('~*~***~*~*~*~**~')
            #print('what the heck is going on with this?????')
            #print(detector.threshold,detector.sec_mean,background_mean)
        '''
        if background_rms<10.0:
            corrected=detector.counts-background_mean
            peak=np.max(corrected)
            max_bin=np.argmax(corrected)
            if peak<LORA.Max_ADC_Count:
                BIN_S=int(max_bin-detector.B_min) # start integration
                BIN_E=int(max_bin+(int(LORA.Sig_Time_Window/2.5))) # end integration
                total_count=np.sum(corrected[BIN_S:BIN_E])-np.sum(corrected[0:(int(LORA.Sig_Time_Window/2.5))])
            else:
                total_count=0
                BIN_S=int(max_bin-detector.B_min) # start integration
                BIN_E=int(max_bin+(int(LORA.Sig_Time_Window/2.5))) # end integration
                total_count=np.sum(corrected[BIN_S:BIN_E])-np.sum(corrected[0:(int(LORA.Sig_Time_Window/2.5))])

            if total_count>0:
                detector.trace_int_counts=total_count

            detector.corrected_peak=peak

    else:
        #print('running counts V2')
        counts_hold=detector.counts*-1.0
        background=counts_hold[0:50]
        background_mean=np.average(background)
        background_rms=np.std(background)
        detector.trace_rms=background_rms
        detector.trace_mean=background_mean
        detector.peak=np.max(counts_hold)
        
        #print(detector.number,background_mean,background_rms)
        '''
        if(background_mean<300):
            detector.corrected_threshold=detector.threshold-background_mean
           
            print('finding threshold from real background: {0}  {1}  {2}'.format(detector.corrected_threshold,detector.threshold,background_mean))
           
        else:
            detector.corrected_threshold=detector.threshold-detector.sec_mean
            print('finding threshold from second background: {0}  {1}  {2}'.format(detector.corrected_threshold,detector.threshold,detector.sec_mean))
        
        
        if detector.corrected_threshold<0:
            detector.corrected_threshold=detector.threshold-detector.sec_mean
           
        if detector.corrected_threshold<0:
            print('~*~***~*~*~*~**~')
            print('what the heck is going on with this?????')
            print(detector.threshold,detector.sec_mean,background_mean)
        '''
        
        detector.corrected_threshold=detector.threshold

        if background_rms<10.0:
            corrected=counts_hold-background_mean
            corrected_mean=np.average(corrected[0:50])
            peak=np.max(corrected)
            #print(peak)
            max_bin=np.argmax(corrected)
            total_count=0
            BIN_S=int(max_bin-20) # start integration
            BIN_E=int(max_bin+(int(LORA.Sig_Time_Window_V2/5.))) # end integration
            total_count=np.sum(corrected[BIN_S:BIN_E])-corrected_mean#-np.sum(corrected[0:(int(LORA.Sig_Time_Window_V2/5.0))])
            #print('total count {0}'.format(total_count))

            '''
            if peak<LORA.Max_ADC_Count_V2:
                BIN_S=int(max_bin-20) # start integration
                BIN_E=int(max_bin+(int(LORA.Sig_Time_Window_V2/5.0))) # end integration
                total_count=np.sum(corrected[BIN_S:BIN_E])-np.sum(corrected[0:(int(LORA.Sig_Time_Window_V2/5.0))])
            else:
                total_count=0
                BIN_S=int(max_bin-20) # start integration
                BIN_E=int(max_bin+(int(LORA.Sig_Time_Window_V2/5.))) # end integration
                total_count=np.sum(corrected[BIN_S:BIN_E])-np.sum(corrected[0:(int(LORA.Sig_Time_Window_V2/5.0))])
            '''
            if total_count>0:
                #print('good count {0}'.format(total_count))
                detector.trace_int_counts=total_count

            detector.corrected_peak=peak
        
        
        
        
        
def find_counts_backgroundV2(detector):

    if detector.number<=20:
   
        temp=remove_noise_via_fft(detector.counts)
        new_trace = remove_simple_baseline(temp[0])
        charge, peak = return_charge_and_peak(new_trace)
        print(charge,peak)
        detector.total_counts_backgroundV2=charge
        '''
        # find background
        background=detector.counts[0:int(LORA.BG_No_Bin)]
        background_mean=np.average(background)
        background_rms=np.std(background)
        detector.trace_rms=background_rms
        detector.trace_mean=background_mean
        detector.peak=np.max(detector.counts)

        if(background_mean<300 and background_mean>0):
            detector.corrected_threshold=detector.threshold-background_mean
    
            #print 'finding threshold from real background: {0}  {1}  {2}'.format(detector.corrected_threshold,detector.threshold,background_mean)
    
        else:
            detector.corrected_threshold=detector.threshold-detector.sec_mean
            #print 'finding threshold from second background: {0}  {1}  {2}'.format(detector.corrected_threshold,detector.threshold,detector.sec_mean)
    
        if detector.corrected_threshold<0:
            detector.corrected_threshold=detector.threshold-detector.sec_mean
    
        if detector.corrected_threshold<0:
            print('~*~***~*~*~*~**~')
            print('what the heck is going on with this?????')
            print(detector.threshold,detector.sec_mean,background_mean)
    
        if background_rms<10.0:
            corrected=detector.counts-background_mean
            peak=np.max(corrected)
            max_bin=np.argmax(corrected)
            if peak<LORA.Max_ADC_Count:
                BIN_S=int(max_bin-detector.B_min) # start integration
                BIN_E=int(max_bin+(int(LORA.Sig_Time_Window/2.5))) # end integration
                total_count=np.sum(corrected[BIN_S:BIN_E])-np.sum(corrected[0:(int(LORA.Sig_Time_Window/2.5))])
            else:
                total_count=0
                BIN_S=int(max_bin-detector.B_min) # start integration
                BIN_E=int(max_bin+(int(LORA.Sig_Time_Window/2.5))) # end integration
                total_count=np.sum(corrected[BIN_S:BIN_E])-np.sum(corrected[0:(int(LORA.Sig_Time_Window/2.5))])

            if total_count>0:
                detector.trace_int_counts=total_count

            detector.corrected_peak=peak

        '''
       
        
        


def remove_noise_via_fft(trace, wpre_index=30, wpost_index=500):
    peak_index=np.argmax(trace)

    if peak_index>3000:
        peak_index = 3000

    tw_start = peak_index - wpre_index
    tw_stop = peak_index + wpost_index

    if tw_start<0:
        tw_stop += -1.*tw_start
        tw_start=0
    elif tw_stop>len(trace):
        tw_start -= (tw_stop - len(trace))
        tw_stop = len(trace)

    tw_start= int(tw_start)
    tw_stop= int(tw_stop)

    assert( (tw_stop-tw_start) == (wpre_index+wpost_index) )

    ontw_trace=trace[tw_start: tw_stop]
    #assert(len(ontw_trace)<1023)
#     zero_pad=np.zeros(1024-len(ontw_trace))
#     ontw_trace=np.append(ontw_trace,zero_pad)

    offtw_trace=trace[tw_stop:]
    #assert(len(offtw_trace)<4095)
    #zero_pad=np.zeros(4096-len(offtw_trace))
    #offtw_trace=np.append(offtw_trace,zero_pad)

    spec_offtw=np.fft.rfft(offtw_trace)
    spec_ontw=np.fft.rfft(ontw_trace)
    spec_full=np.fft.rfft(trace)

    # freq_offtw=np.fft.rfftfreq(len(offtw_trace),d=2.5e-9)#2.5ns sample spacing
    # freq_ontw=np.fft.rfftfreq(len(ontw_trace),d=2.5e-9)#2.5ns sample spacing
    freq_full=np.fft.rfftfreq(len(trace),d=2.5e-9)#2.5ns sample spacing

    # ratio_with_offtw=np.array([])
    # store_loc=0
    # for kk in range(len(freq_ontw)):
    #     for ll in range(store_loc,len(freq_offtw)):
    #         if freq_ontw[kk]==freq_offtw[ll]:
    #             store_loc=ll
    #             break
    #     r=np.log10(np.absolute(spec_offtw[store_loc])/np.absolute(spec_ontw[kk]))
    #     ratio_with_offtw=np.append(ratio_with_offtw,r)
    #
    # ratio_with_full=np.array([])
    # store_loc=0
    # for kk in range(len(freq_ontw)):
    #     for ll in range(store_loc,len(freq_full)):
    #         if freq_ontw[kk]==freq_full[ll]:
    #             store_loc=ll
    #             break
    #     r=np.log10(np.absolute(spec_full[store_loc])/np.absolute(spec_ontw[kk]))
    #     ratio_with_full=np.append(ratio_with_full,r)

    spec_full[freq_full>75e6]=0+0j
    rev_trace=np.fft.irfft(spec_full)
    rev_trace=remove_simple_baseline(rev_trace)

    return rev_trace, spec_ontw, spec_offtw
    #, spec_full
    # , ratio_with_offtw, ratio_with_full


def remove_simple_baseline(trace, wpre_index=30, wpost_index=320, offtw_len=100):
        peak_index=np.argmax(trace)
        if peak_index>3000:
            peak_index = 3000

        tw_start = peak_index - wpre_index
        tw_stop = peak_index + wpost_index

        if tw_start<0:
            tw_stop += -1.*tw_start
            tw_start=0
        elif tw_stop>len(trace):
            tw_start -= (tw_stop - len(trace))
            tw_stop = len(trace)

        tw_start= int(tw_start)
        tw_stop= int(tw_stop)

        assert( (tw_stop-tw_start) == (wpre_index+wpost_index) )

        #find offwtrace
        # offwtrace = trace[tw_start:0:-1]
        #if len(offwtrace)<offtw_len:
            #print ("offwtracelen:",len(offwtrace))
        #    offwtrace=np.concatenate([offwtrace, trace[tw_stop:]])
        # offwtrace=offwtrace[:offtw_len]


        offwtrace1 = trace[tw_start:0:-1][:offtw_len]
        offwtrace2 = trace[tw_stop:][:offtw_len]
        offwtrace=np.concatenate([offwtrace1, offwtrace2])

        offwmean = np.mean(offwtrace)

        final_trace = trace - offwmean

        return final_trace

def return_charge_and_peak(trace, wpre_index=30, wpost_index=320):
        peak_index=np.argmax(trace)
        if peak_index>3000:
            peak_index = 3000

        tw_start = peak_index - wpre_index
        tw_stop = peak_index + wpost_index

        if tw_start<0:
            tw_stop += -1.*tw_start
            tw_start=0
        elif tw_stop>len(trace):
            tw_start -= (tw_stop - len(trace))
            tw_stop = len(trace)

        tw_start= int(tw_start)
        tw_stop= int(tw_stop)

        assert( (tw_stop-tw_start) == (wpre_index+wpost_index) )

        peak= trace[peak_index]
        charge = np.sum( trace[tw_start : tw_stop] )

        return charge, peak


def retrive_sat_signal(detector):
    test=0


def get_arrival_time(detector):

    if detector.number<=20:
        cut=LORA.Det_Thres*detector.trace_rms+detector.trace_mean
        flag=0
        for i in np.arange(LORA.nTrace):
            if detector.counts[i]>cut and flag==0:
                if i<400:
                    continue
                else:
                    detector.threshold_time=i*2.5*10  # unit of 0.1 ns
                    flag=1
                    
    else:
        
        cut=1*detector.threshold*.5
        
        
        
        flag=0
        print('number {0},counts {1}'.format(detector.number,detector.trace_int_counts))
        print(np.max(np.abs(detector.counts)),cut)
        for i in np.arange(LORA.nTrace):
            #print(detector.threshold_time_no_trig)
            #if -1*detector.counts[i]>(45) and flag==0:
            #    if i>50:
            #        detector.threshold_time_no_trig=i*5.0*10
            #print('no trigger crossing ',i,cut,-1* detector.counts[i])


            if -1*detector.counts[i]>cut and flag==0:
                if i<50:
                
                    continue
                else:
                    print('found crossing ',i,cut,-1* detector.counts[i])
                    detector.threshold_time=i*5.0*10  # unit of 0.1 ns
                    flag=1



def get_event_timestamp(detector,lasa):
    print('_______event timestamp______')
    print(lasa.CTP)
    print(lasa.sec_flag)
    if lasa.sec_flag!=1 and lasa.CTP[1]>0:
        detector.event_time_stamp=10*(lasa.sync[0]+lasa.quant[1]+(1.0*detector.ctd/lasa.CTP[1])*(1000000000.0-lasa.quant[1]+lasa.quant[2]))
        print('doing real time stamp ')
        print(detector.event_time_stamp)
    else:
        print('doing est. time stamp ')

        detector.event_time_stamp=10*detector.nsec#10*(1.0*detector.ctd/200000000.0 )*(1000000000.0)
        print(detector.event_time_stamp)

def get_event_timestamp_V2(detector,lasa):
    
    #print detector.number, detector.number%2
    if detector.number<=20:
        if detector.number%2==1:
            if lasa.sec_flag!=1 and lasa.CTP_M[1]>0:
                detector.event_time_stamp=10*(lasa.sync_M[0]+lasa.quant_M[1]+(1.0*detector.ctd/lasa.CTP_M[1])*(1000000000.0-lasa.quant_M[1]+lasa.quant_M[2]))
                #print detector.event_time_stamp
            else:
                #print 'flagged event'
                detector.event_time_stamp=10*detector.nsec
                #print detector.event_time_stamp

        
        if detector.number%2==0 and lasa.CTP_S[1]>0:
            if lasa.sec_flag!=1:
                detector.event_time_stamp=10*(lasa.sync_S[0]+lasa.quant_S[1]+(1.0*detector.ctd/lasa.CTP_S[1])*(1000000000.0-lasa.quant_S[1]+lasa.quant_S[2]))
                #print detector.event_time_stamp
            else:
                #print 'flagged event'
                detector.event_time_stamp=10*detector.nsec
                #print detector.event_time_stamp
                
    else:
        if lasa.sec_flag==0:  # hack for the moment because something is off here
            detector.event_time_stamp=10*(lasa.sync[0]+lasa.quant[1]+(1.0*detector.ctd/lasa.CTP[1])*(1000000000.0-lasa.quant[1]+lasa.quant[2]))
            #detector.event_time_stamp=10*((1.0*detector.ctd/lasa.CTP[1]))*1e9
            #print(lasa.sync[0],lasa.quant[1],detector.ctd,detector.event_time_stamp,lasa.CTP[1],lasa.quant[2])
        else:
            print('flagged event')
            detector.event_time_stamp=10*detector.nsec
            #print detector.event_time_stamp



def cal_event_timestamp(detectors,lasa):
    #print '_____________________________________'
    #print 'lasa number: {0}'.format(lasa.number)
    lasa_ind=int(lasa.number-1)
    #print('lasa index: {0}'.format(lasa_ind))
    #print('trigger condition: {0}'.format(detectors[4*lasa_ind].trigg_condition))
    trigg_cond=detectors[4*lasa_ind].trigg_condition
    thresh_times=np.asarray([detectors[4*int(lasa.number-1)].threshold_time,detectors[4*int(lasa.number-1)+1].threshold_time,detectors[4*int(lasa.number-1)+2].threshold_time,detectors[4*int(lasa.number-1)+3].threshold_time])
    #print(thresh_times)
    thresh_use=1.0*thresh_times[thresh_times!=0]
    #print('thresh_use: ',thresh_use)
    args=np.argsort(thresh_use)
    #print(int(args[1]-1))

    if len(thresh_use)>1:
        print('trigger condition: ',trigg_cond)
        print('sorted: {0}'.format(np.sort(thresh_use)))
        #print thresh_use
        #print args
        #print 'use index: {0}'.format(args[int(trigg_cond)-1])
    
        if lasa.number<=5:
            try:
                trigg_time=thresh_use[args[int(trigg_cond)-1]]
            except: #this is in here becasue in a few cases the length of thresh_use was less than trigg. condition
                trigg_time=np.sort(thresh_use)[len(thresh_use)-1]
            print('V1: ',thresh_times,trigg_time)
 
    

            for i in np.arange(4):
                if detectors[4*int(lasa_ind)+i].threshold_time>0:
                    detectors[4*int(lasa_ind)+i].cal_time=detectors[4*int(lasa_ind)+i].threshold_time-trigg_time+detectors[lasa_ind*4].event_time_stamp+det.cable_delay[4*int(lasa_ind)+i]
                    detectors[4*int(lasa_ind)+i].final_event_time=detectors[4*int(lasa_ind)+i].cal_time
                    print(detectors[4*int(lasa_ind)+i].number,  detectors[4*int(lasa_ind)+i].final_event_time)

        else:
            try:
                trigg_time=thresh_use[args[0]]
                print('found trigger time')
            except: #this is in here becasue in a few cases the length of thresh_use was less than trigg. condition
                trigg_time=np.sort(thresh_use)[len(thresh_use)-1]
            print('V2: ',thresh_times,trigg_time)
            for i in np.arange(4):
                print(4*int(lasa_ind)+i,detectors[4*int(lasa_ind)+i].threshold_time)
                if detectors[4*int(lasa_ind)+i].threshold_time>0:
                    detectors[4*int(lasa_ind)+i].cal_time=detectors[4*int(lasa_ind)+i].threshold_time-trigg_time+detectors[lasa_ind*4].event_time_stamp+det.cable_delay[4*int(lasa_ind)+i]
                    detectors[4*int(lasa_ind)+i].final_event_time=detectors[4*int(lasa_ind)+i].cal_time
                    print('cal time: ',int(detectors[4*int(lasa_ind)+i].cal_time))
                    print('event time: ',int(detectors[4*int(lasa_ind)+i].final_event_time))
                    

'''
def cal_event_timestamp(detectors,lasa):
    #print '_____________________________________'
    #print 'lasa number: {0}'.format(lasa.number)
    lasa_ind=int(lasa.number-1)
    print('lasa index: {0}'.format(lasa_ind))
    print('trigger condition: {0}'.format(detectors[4*lasa_ind].trigg_condition))
    trigg_cond=detectors[4*lasa_ind].trigg_condition
    thresh_times=np.asarray([detectors[4*int(lasa.number-1)].threshold_time,detectors[4*int(lasa.number-1)+1].threshold_time,detectors[4*int(lasa.number-1)+2].threshold_time,detectors[4*int(lasa.number-1)+3].threshold_time])
    thresh_use=1.0*thresh_times[thresh_times!=0]
    print('thresh_use: ',thresh_use)
    args=np.argsort(thresh_use)
    #print int(args[detectors[0].trigg_condition]-1)
    if len(thresh_use)>1:
        print('trigger condition: ',trigg_cond)
        print('sorted: {0}'.format(np.sort(thresh_use)))
        #print thresh_use
        #print args
        #print 'use index: {0}'.format(args[int(trigg_cond)-1])
        
        if lasa.number<=5:
            try:
                trigg_time=thresh_use[args[int(trigg_cond)-1]]
            except: #this is in here becasue in a few cases the length of thresh_use was less than trigg. condition
                trigg_time=np.sort(thresh_use)[len(thresh_use)-1]
            print('V1: ',thresh_times,trigg_time)
        else:
            try:
                trigg_time=thresh_use[args[0]]
            except: #this is in here becasue in a few cases the length of thresh_use was less than trigg. condition
                trigg_time=np.sort(thresh_use)[len(thresh_use)-1]
            print('V2: ',thresh_times,trigg_time)
        
        

        for i in np.arange(4):
            if detectors[4*int(lasa_ind)+i].threshold_time>0:
            
                detectors[4*int(lasa_ind)+i].cal_time=detectors[4*int(lasa_ind)+i].threshold_time-trigg_time+detectors[lasa_ind*4].event_time_stamp+det.cable_delay[4*int(lasa_ind)+i]
                detectors[4*int(lasa_ind)+i].final_event_time=detectors[4*int(lasa_ind)+i].cal_time
                print(detectors[4*int(lasa_ind)+i].number,  detectors[4*int(lasa_ind)+i].final_event_time)

                # maybe this has to be corrected for wrap-around seconds

            #print '{0}   {1}  {2}  {3}  {4}  {5}'.format(4*int(lasa_ind)+i,int(detectors[4*int(lasa_ind)+i].cal_time),detectors[4*int(lasa_ind)+i].threshold_time, trigg_time,detectors[4*int(lasa_ind)+i].event_time_stamp,det.cable_delay[4*int(lasa_ind)+i])


'''

def do_arrival_time_diff(detectors):
    event_times=np.zeros([len(detectors)])
    event_weight=np.zeros([len(detectors)])

    ind_0=1000
    time_min=1e11
    for i in np.arange(len(detectors)):
        #print(detectors[i].trace_int_counts)
        if detectors[i].trace_int_counts>5:
            event_times[i]=detectors[i].final_event_time
            event_weight[i]=detectors[i].trace_int_counts/detectors[i].gain
            #print(i, time_min,event_times[i],event_weight[i])
            if event_times[i]<time_min and event_times[i]!=0:
                time_min=event_times[i]
                ind_0=i
    #print ind_0, event_times[ind_0]

    for i in np.arange(len(detectors)):
        if(event_times[i]>0 and event_weight[i]>0):
            detectors[i].cdt=(event_times[i]-event_times[ind_0])*0.1*(1.e-9*LORA.vel_light)

        #print(i,detectors[i].cdt)
def do_arrival_direction(detectors,event):

    P=0
    Q=0
    R=0
    S=0
    W=0
    T1=0
    S1=0
    S2=0
    S3=0
    S4=0
    S5=0
    S6=0
    
    counter=0
    
    for i in np.arange(len(detectors)):
        if detectors[i].cdt>=0:
            counter=counter+1
            S=S+detectors[i].x_cord**2
            W=W+detectors[i].y_cord**2
            
            T1=T1+detectors[i].x_cord
            S1=S1+detectors[i].x_cord*detectors[i].cdt
            S2=S2+detectors[i].x_cord*detectors[i].y_cord
            S3=S3+detectors[i].y_cord*detectors[i].cdt
            S4=S4+detectors[i].y_cord
            S5=S5+detectors[i].cdt
            S6=S6+1
    if counter>2:
        P=(S1*S2)/S
        Q=(T1*S2)/S
        R=-((S2*S2)/S)+W

        t0=(T1*S1*R-R*S*S5+T1*P*S2-T1*S2*S3-S*S4*P+S*S4*S3)/(T1*Q*S2-T1*S2*S4+R*T1*T1-S*S4*Q+S*S4*S4-R*S*S6)
        m=(P-t0*Q-S3+t0*S4)/R
        l=(-S1/S)-((P*S2)/(R*S))+((t0*Q*S2)/(R*S))+((S2*S3)/(R*S))-((t0*S2*S4)/(R*S))+((t0*T1)/S)
        #print(l,m)
        n=np.sqrt(1.0-(l*l+m*m))
    
    
    
        #print('direction params:  {0} {1}  {2}  {3}'.format(t0,m,l,n))

        if l*l+m*m<1:

            theta=(np.arcsin(np.sqrt(l*l+m*m)))*(180.0/np.pi)    #Zenith in degrees (from vertical direction +Z)
            phi=(np.arccos(m/np.sqrt(l*l+m*m)))*(180.0/np.pi)    #in degrees (Eastward from North)
        else:
            theta=0.0
            phi=(np.arccos(m/np.sqrt(l*l+m*m)))*(180.0/np.pi)    #in degrees (Eastward from North)

        
        if l<0:
            phi=360.0-phi



        event.theta=theta
        event.phi=phi
        event.elevation=90.0-theta
        print('theta: {0:.2f}   phi: {1:.2f}    el: {2:.2f}'.format(event.theta,event.phi,event.elevation))
    else:
        print('not enough stations to reconstruct direction')
        event.direction_flag=1

def do_COM_core(detectors,event):

    #print 'COM core'
    x=0
    Event_Size=0
    SumX=0
    SumY=0
    
    for i in np.arange(len(detectors)):
    
        x=detectors[i].trace_int_counts/detectors[i].gain/(LORA.Det_Area*np.cos(event.theta*np.pi/180.0))
        
        if(detectors[i].final_event_time>0 and x>LORA.Density_Cut):
        
            Event_Size=Event_Size+detectors[i].trace_int_counts/detectors[i].gain
            SumX=SumX+detectors[i].trace_int_counts/detectors[i].gain*detectors[i].x_cord
            SumY=SumY+detectors[i].trace_int_counts/detectors[i].gain*detectors[i].y_cord
        

    if event.direction_flag==0:
        print('Event size: {0}'.format(Event_Size))
        event.x_core=SumX/Event_Size
        event.y_core=SumY/Event_Size
        event.z_core=0     #we assume that all the LORA detectors are at z=0

        print('core: ({0:.2f}, {1:.2f}, {2:.2f})'.format(event.x_core,event.y_core,event.z_core))


def find_density(detectors,event):
     

    for i in np.arange(len(detectors)):
        detectors[i].density=detectors[i].trace_int_counts/detectors[i].gain/(LORA.Det_Area*np.cos(event.theta*np.pi/180.0))
        detectors[i].err_density=np.power(detectors[i].density,0.5) ;    #Assuming possionian error
        if detectors[i].density<=LORA.Density_Cut or detectors[i].threshold_time==0:
            detectors[i].density=0
            detectors[i].err_density=0



def func_plane(x,par):    #//Plane function to fit the shower front: (Parameters: t0,l,m)
    return par[0]-par[1]*x[0]-par[2]*x[1]

def theta_phi(theta,phi,psi,x0,y0,z0):

    #1st ROTATION: counterclockwise from X-axis(N) for angle 'phi' about Z-axis.
    x1= x0*np.cos(phi)+y0*np.sin(phi)
    y1=-x0*np.sin(phi)+y0*np.cos(phi)
    z1= z0
    #-------------xxx----------------
    #2nd ROTATION: clockwise from Z-axis for angle 'theta' about Y-axis(W).
    x2= x1*np.cos(theta)-z1*np.sin(theta) ;
    y2= y1
    z2= x1*np.sin(theta)+z1*np.cos(theta) ;
    #-------------xxx----------------
    #3rd ROTATION: counterclockwise from X-axis(N) for angle 'psi' about Z-axis.
    x3= x2*np.cos(psi)+y2*np.sin(psi)
    y3=-x2*np.sin(psi)+y2*np.cos(psi)
    z3= z2
    #-------------xxx----------------
    return x3,y3,z3
    
#NKG
def NKG(r, r_M, N_ch):
    s = 1.7
    C = gamma(4.5 - s) / (2 * np.pi * r_M**2 * gamma(s) * gamma(4.5 - 2 * s))
    return N_ch * C * (r / r_M)**(s - 2) * (1 + r / r_M)**(s - 4.5)
#NKG in Log
def NKG_log(r, r_M, N_ch):
    s = 1.7
    C = gamma(4.5 - s) / (2 * np.pi * r_M**2 * gamma(s) * gamma(4.5 - 2 * s))
    return np.log(N_ch * C * (r / r_M)**(s - 2) * (1 + r / r_M)**(s - 4.5))

# Log-Likelihood function
def log_likelihood(params, distance, detect_dens_filtered, sigma):
    r_M, N_ch = params
    model = NKG(distance, r_M, N_ch)
    return -0.5 * np.sum(((detect_dens_filtered - model) / sigma) ** 2)

# Log-Prior function
def log_prior(params):
    r_M, N_ch = params
    if 1 < r_M < 250 and 1e3 < N_ch < 1e9:
       
        return 0.0  
    return -np.inf  

# Log-Posterior function
def log_posterior(params, distance, detect_dens_filtered, sigma):
    lp = log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(params, distance, detect_dens_filtered, sigma)

def NKG_core_log(coords, xc, yc, r_M, N_ch):
    x, y = coords
    s = 1.7  
    r = np.sqrt((xc - x)**2 + (yc - y)**2)
    C = gamma(4.5 - s) / (2 * np.pi * r_M**2 * gamma(s) * gamma(4.5 - 2 * s))
    return np.log(N_ch * C * (r / r_M)**(s - 2) * (1 + r / r_M)**(s - 4.5))

def NKG_core(coords, xc, yc, r_M, N_ch):
    x, y = coords
    s = 1.7  
    r = np.sqrt((xc - x)**2 + (yc - y)**2)
    C = gamma(4.5 - s) / (2 * np.pi * r_M**2 * gamma(s) * gamma(4.5 - 2 * s))
    return N_ch * C * (r / r_M)**(s - 2) * (1 + r / r_M)**(s - 4.5)

# Log-likelihood core function
def likelihood_core(params, coords, detect_dens_filtered, sigma):
    xc, yc, r_M, N_ch = params
    model = NKG_core(coords, xc, yc, r_M, N_ch)
    return -0.5 * np.sum(((detect_dens_filtered - model) / sigma) ** 2)

# Log-prior core function
def prior_core(params):
    xc, yc, r_M, N_ch = params
    
    if -500 < xc < 500 and -500 < yc < 500 and 1 < r_M < 250 and 1e3 < N_ch < 1e9:
        return 0.0  
    return -np.inf  

# Log-posterior core function
def posterior_core(params, coords, detect_dens_filtered, sigma):
    lp = prior_core(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + likelihood_core(params, coords, detect_dens_filtered, sigma)

# Residual function (sum of squared residuals)
def residu(params, r, detect_dens_filtered):
    r_M, N_ch = params
    model = NKG(r, r_M, N_ch)
    return np.sum((model - detect_dens_filtered)**2)

#  function for minimize in log space
def residu_log(params, r, log_detect_dens, log_errors):
    r_M, N_ch = params
    residuals = (NKG_log(r, r_M, N_ch) - log_detect_dens)/log_errors
    return np.sum((residuals)**2)

#  function for minimize in log space
def residu_core(params, r, detect_dens):
    xc, yc, r_M, N_ch = params
    return np.sum((NKG_core(r, xc, yc, r_M, N_ch) - detect_dens)**2)

#  function for minimize in log space
def residu_core_log(params, r, log_detect_dens):
    xc, yc, r_M, N_ch = params
    return np.sum((NKG_core_log(r, xc, yc, r_M, N_ch) - log_detect_dens)**2)

# RMSE (Root Mean Square Error) function to measure fit quality
def rmse(params, r, detect_dens_filtered):
    r_M, N_ch = params
    model = NKG(r, r_M, N_ch)
    return np.sqrt(np.mean((model - detect_dens_filtered)**2))

def model_function(xy, t0, l, m):
    x, y = xy
    return t0 + l * x + m * y

def reduced_chi_squared_with_errors(observed, r, model, params, param_errors, errors=None):
    s = 1.7
    residuals = observed - model
    n = len(observed)  
    p = len(params)    

    # Assuming Poissonian error if no errors are provided
    if errors is None:
        errors = np.sqrt(observed)
    
    # Prevent division by zero in error propagation
    errors[errors == 0] = 1e-10  

    # --- Error Propagation ---
    
    r_M, N_ch = params
    r_M_err, N_ch_err = param_errors

    # Partial derivatives of the model wrt each parameter
    p1 = -(gamma(4.5 - s) * N_ch * r * (s-4.5) * (r / r_M)**(s - 2) * (1 + r / r_M)**(s - 5.5) / (2 * np.pi * r_M**4 * gamma(s) * gamma(4.5 - 2 * s)))
    p2 = -(gamma(4.5 - s) * N_ch * r * (s-2) * (r / r_M)**(s - 3) * (1 + r / r_M)**(s - 4.5) / (2 * np.pi * r_M**4 * gamma(s) * gamma(4.5 - 2 * s)))
    p3 = -(gamma(4.5 - s) * N_ch * (r / r_M)**(s - 2) * (1 + r / r_M)**(s - 4.5) / (np.pi * r_M**3 * gamma(s) * gamma(4.5 - 2 * s)))
    partial_r_M = p1 + p2 + p3       # ∂NKG/∂r_M
    partial_N_ch = model / N_ch      # ∂NKG/∂N_ch

    # Propagate the parameter errors to get uncertainty in the model predictions
    model_errors = np.sqrt(
        (partial_r_M * r_M_err) ** 2 +
        (partial_N_ch * N_ch_err) ** 2
    )

    # Calculate chi-squared using both observational and model errors
    total_errors = np.sqrt(errors**2 + model_errors**2)

    # Calculate chi-squared
    chi_squared = np.sum((residuals / total_errors) ** 2)

    # Degrees of freedom (dof = number of data points - number of fitted parameters)
    dof = n - p

    # Calculate reduced chi-squared
    reduced_chi_sq = chi_squared / dof

    return reduced_chi_sq   

def core_reduced_chi_squared_with_errors(observed, coords, model, params, param_errors, errors=None):
    s = 1.7
    x, y = coords
    residuals = observed - model
    n = len(observed)  
    p = len(params)    

    # Assuming Poissonian error if no errors are provided
    if errors is None:
        errors = np.sqrt(observed)
    
    # Prevent division by zero in error propagation
    errors[errors == 0] = 1e-10  

    # --- Error Propagation ---
    
    xc, yc, r_M, N_ch = params
    xc_err, yc_err, r_M_err, N_ch_err = param_errors

    r = np.sqrt((xc - x)**2 + (yc - y)**2)

    # Partial derivatives of the model wrt each parameter
    p1 = -(gamma(4.5 - s) * N_ch * r * (s-4.5) * (r / r_M)**(s - 2) * (1 + r / r_M)**(s - 5.5) / (2 * np.pi * r_M**4 * gamma(s) * gamma(4.5 - 2 * s)))
    p2 = -(gamma(4.5 - s) * N_ch * r * (s-2) * (r / r_M)**(s - 3) * (1 + r / r_M)**(s - 4.5) / (2 * np.pi * r_M**4 * gamma(s) * gamma(4.5 - 2 * s)))
    p3 = -(gamma(4.5 - s) * N_ch * (r / r_M)**(s - 2) * (1 + r / r_M)**(s - 4.5) / (np.pi * r_M**3 * gamma(s) * gamma(4.5 - 2 * s)))
    p4 = N_ch * gamma(4.5 - s) / (2 * np.pi * r_M**2 * gamma(s) * gamma(4.5 - 2 * s))*(s-4.5)*(xc-x)*(np.sqrt((xc - x)**2 + (yc - y)**2) / r_M)**(s - 2)*(1 + np.sqrt((xc - x)**2 + (yc - y)**2) / r_M)**(s - 5.5)/r_M*np.sqrt((xc - x)**2 + (yc - y)**2)
    p5 = N_ch * gamma(4.5 - s) / (2 * np.pi * r_M**2 * gamma(s) * gamma(4.5 - 2 * s))*(s-2)*(xc-x)*(np.sqrt((xc - x)**2 + (yc - y)**2) / r_M)**(s-3)*(1 + np.sqrt((xc - x)**2 + (yc - y)**2) / r_M)**(s - 4.5)/r_M*np.sqrt((xc - x)**2 + (yc - y)**2)
    p6 = N_ch * gamma(4.5 - s) / (2 * np.pi * r_M**2 * gamma(s) * gamma(4.5 - 2 * s))*(s-4.5)*(yc-x)*(np.sqrt((xc - x)**2 + (yc - y)**2) / r_M)**(s - 2)*(1 + np.sqrt((xc - x)**2 + (yc - y)**2) / r_M)**(s - 5.5)/r_M*np.sqrt((xc - x)**2 + (yc - y)**2)
    p7 = N_ch * gamma(4.5 - s) / (2 * np.pi * r_M**2 * gamma(s) * gamma(4.5 - 2 * s))*(s-2)*(yc-x)*(np.sqrt((xc - x)**2 + (yc - y)**2) / r_M)**(s-3)*(1 + np.sqrt((xc - x)**2 + (yc - y)**2) / r_M)**(s - 4.5)/r_M*np.sqrt((xc - x)**2 + (yc - y)**2)

    partial_xc  = p4 + p5            # ∂NKG/∂xc
    partial_yc  = p6 + p7            # ∂NKG/∂yc 
    partial_r_M = p1 + p2 + p3       # ∂NKG/∂r_M
    partial_N_ch = model / N_ch      # ∂NKG/∂N_ch

    # Propagate the parameter errors to get uncertainty in the model predictions
    model_errors = np.sqrt(
        (partial_xc * xc_err) ** 2 +
        (partial_yc * yc_err) ** 2 +
        (partial_r_M * r_M_err) ** 2 +
        (partial_N_ch * N_ch_err) ** 2
    )

    # Calculate chi-squared using both observational and model errors
    total_errors = np.sqrt(errors**2 + model_errors**2)

    # Calculate chi-squared
    chi_squared = np.sum((residuals / total_errors) ** 2)

    # Degrees of freedom (dof = number of data points - number of fitted parameters)
    dof = n - p

    # Calculate reduced chi-squared
    reduced_chi_sq = chi_squared / dof

    return reduced_chi_sq 


def minimize_core_fit(detect_coords, detect_dens, cx0, cy0, theta, phi, rM0, Nch0, n_iters, detect_cord_x, detect_cord_y):
    initial_guess = [cx0, cy0, rM0, Nch0]
    tolerance = 1e-6
    param_tolerance = 1e-3
    prev_rmse = np.inf
    prev_params = initial_guess

    # Filtered data (positive values)
    detectUVW = prf.GetUVW(detect_coords, cx0, cy0, 0, theta, phi, Binc=1.1837)
    distance = np.linalg.norm(detectUVW, axis=1)
    positive_mask = detect_dens > 0
    distance_filtered = distance[positive_mask]
    detect_dens_filtered = detect_dens[positive_mask]
    log_detect_density_filtered = np.log(detect_dens_filtered)
    x_filtered = detect_cord_x[0:40][positive_mask]
    y_filtered = detect_cord_y[0:40][positive_mask]

    coords = np.vstack((x_filtered, y_filtered))

    boundscf = ([-500, -500, 1, 1e3], [500, 500, 500, 1e10]) 
    bounds = ([-500, 500], [-500, 500], [1, 300], [1e3, 1e10])

    for i in range(n_iters):
        print(f"\nIteration {i + 1}")
        
        # Perform minimize
        result =  minimize(residu_core_log, initial_guess, args=(coords, log_detect_density_filtered),
                  method='L-BFGS-B', bounds=bounds)
        
        # Extract fitted parameters
        fitted_xc, fitted_yc, fitted_r_M, fitted_N_ch = result.x
        fitted_params = [fitted_xc, fitted_yc, fitted_r_M, fitted_N_ch]

        # Calculate RMSE for this iteration
        current_rmse = np.sqrt(np.mean((NKG_core_log(coords, fitted_xc, fitted_yc, fitted_r_M, fitted_N_ch) - log_detect_density_filtered)**2))
        
        print(f"Fitted xc: {fitted_xc:.6f}")
        print(f"Fitted yc: {fitted_yc:.6f}")
        print(f"Fitted r_M: {fitted_r_M:.6f}")
        print(f"Fitted N_ch: {fitted_N_ch:.6f}")
        print(f"RMSE: {current_rmse:.6f}")
        
        # Check convergence
        if abs(prev_rmse - current_rmse) < tolerance:
            print("Convergence reached: RMSE improvement below tolerance.")
            break
        
        param_change = np.linalg.norm(np.array(prev_params) - np.array(fitted_params))
        if param_change < param_tolerance:
            print("Convergence reached: Parameter change below tolerance.")
            break
        
        prev_rmse = current_rmse
        prev_params = fitted_params
        initial_guess = fitted_params

  
    
      # Perform least squares fitting using curve_fit to obtain 1-sigma errors
    try:
        popt, pcov = curve_fit(NKG_core_log, coords, log_detect_density_filtered, p0=fitted_params, bounds=boundscf, maxfev=20000)
        param_errors = np.sqrt(np.diag(pcov))
    except RuntimeError as e:
        print("Curve fitting did not converge:", e)
        popt = fitted_params
        param_errors = [np.inf] * len(popt)

    print(f"Final Fit xc: {popt[0]:.6f} ± {param_errors[0]:.6f}")
    print(f"Final Fit yc: {popt[1]:.6f} ± {param_errors[1]:.6f}")
    print(f"Final Fit r_M: {popt[2]:.6f} ± {param_errors[2]:.6f}")
    print(f"Final Fit N_ch: {popt[3]:.6f} ± {param_errors[3]:.6f}")

    return fitted_xc, fitted_yc, fitted_r_M, fitted_N_ch, param_errors[0], param_errors[1], param_errors[2], param_errors[3]


def Minimize_log_fit(detect_coords, detect_dens, core_x, core_y, theta, phi, r_M0, N_ch0, n_iterations):
    
    tolerance = 1e-6
    param_tolerance = 1e-3
    initial_guess = [r_M0, N_ch0]
    prev_rmse = np.inf
    prev_params = initial_guess

    detectUVW=prf.GetUVW(detect_coords, core_x, core_y, 0, theta, phi, Binc=1.1837)
    distance = np.linalg.norm(detectUVW, axis=1)

    # Filtered data (positive values)
    positive_mask = detect_dens > 0
    distance_filtered = distance[positive_mask]
    detect_dens_filtered = detect_dens[positive_mask]
    log_detect_density_filtered = np.log(detect_dens_filtered)

    #erros
    log_errors = 1 / (np.sqrt(detect_dens_filtered) * np.log(10))

    boundscf = ([1, 1e3], [500, 1e10])  # Bounds for r_M and N_ch
    bounds   = ([1, 500], [1e3, 1e10])

    for i in range(n_iterations):
        print(f"\nIteration {i + 1}")
        
        # Perform minimize
        result = minimize(residu_log, initial_guess, args=(distance_filtered, log_detect_density_filtered, log_errors), 
                        method='L-BFGS-B', bounds=bounds)
        
        # Extract fitted parameters
        fitted_r_M, fitted_N_ch = result.x
        fitted_params = [fitted_r_M, fitted_N_ch]

        # Calculate RMSE for this iteration
        current_rmse = np.sqrt(np.mean((NKG_log(distance_filtered, fitted_r_M, fitted_N_ch) - log_detect_density_filtered)**2))
        
        print(f"Fitted r_M: {fitted_r_M:.6f}")
        print(f"Fitted N_ch: {fitted_N_ch:.6f}")
        print(f"RMSE: {current_rmse:.6f}")
        
        # Check convergence
        if abs(prev_rmse - current_rmse) < tolerance:
            print("Convergence reached: RMSE improvement below tolerance.")
            break
        
        param_change = np.linalg.norm(np.array(prev_params) - np.array(fitted_params))
        if param_change < param_tolerance:
            print("Convergence reached: Parameter change below tolerance.")
            break
        
        prev_rmse = current_rmse
        prev_params = fitted_params
        initial_guess = fitted_params



    # Use the results from minimize as initial guesses for curve_fit
    initial_guesses = [fitted_r_M, fitted_N_ch]
    popt, pcov = curve_fit(NKG_log, distance_filtered, log_detect_density_filtered, p0=initial_guesses, bounds=boundscf, maxfev=5000)

    r_M_fit, N_ch_fit = popt
    param_errors = np.sqrt(np.diag(pcov))

    fitted_r_M_err, fitted_N_ch_err = param_errors

    print(f"Final Fit r_M: {r_M_fit:.6f} ± {fitted_r_M_err:.6f}")
    print(f"Final Fit N_ch: {N_ch_fit:.6f} ± {fitted_N_ch_err:.6f}")


    param_errors = [fitted_r_M_err, fitted_N_ch_err]  # Replace with actual values obtained from the covariance matrix or numerical methods

    # Calculate the model predictions for the fitted parameters
    model_fitted = NKG_log(distance_filtered, r_M_fit, N_ch_fit)

    # Calculate reduced chi-squared with propagated parameter errors
    reduced_chi_sq_with_errors = reduced_chi_squared_with_errors(detect_dens_filtered, distance_filtered, model_fitted, [r_M_fit, N_ch_fit], param_errors)

    # Print the reduced chi-squared result
    print(f"Reduced Chi-Squared (with parameter errors): {reduced_chi_sq_with_errors}")

    return reduced_chi_sq_with_errors, r_M_fit, N_ch_fit, fitted_r_M_err, fitted_N_ch_err



def NKG_fit_trial(detectors, event):
    
    theta=event.fit_theta
    phi=event.fit_phi
    psi=2*np.pi-phi
    
    x_shower=np.zeros([len(detectors)])
    y_shower=np.zeros([len(detectors)])
    x_pos=np.zeros([len(detectors)])
    y_pos=np.zeros([len(detectors)])
    z_shower=np.zeros([len(detectors)])

    temp_det_array=[]
    temp_den_array=[]
    temp_x_pos=[]
    temp_y_pos=[]
    detect_dens=[]

    # find core position in shower plane for first guess for fit. 4 densest detectors are used
    
    for i in np.arange(len(detectors)):
        x_shower[i],y_shower[i],z_shower[i]=theta_phi(theta,phi,psi,detectors[i].x_cord,detectors[i].y_cord,detectors[i].z_cord)
        x_pos[i]=detectors[i].x_cord
        y_pos[i]=detectors[i].y_cord
        if detectors[i].density>=LORA.Density_Cut and detectors[i].density<=LORA.Density_Cut_High:
            temp_den_array.append(detectors[i].density)
            temp_det_array.append(i)
            temp_x_pos.append(x_shower[i])
            temp_y_pos.append(x_shower[i])

    ind=np.argsort(np.asarray(temp_den_array))[::-1]
    temp_den_array=np.asarray(temp_den_array)
    temp_det_array=np.asarray(temp_det_array)
    temp_x_array=np.asarray(temp_x_pos)
    temp_y_array=np.asarray(temp_y_pos)

    temp_total_den=temp_den_array[ind[0]]+temp_den_array[ind[1]]+temp_den_array[ind[2]]+temp_den_array[ind[3]]
    x_core=(x_shower[temp_det_array[ind[0]]]*temp_den_array[ind[0]]+x_shower[temp_det_array[ind[1]]]*temp_den_array[ind[1]]+x_shower[temp_det_array[ind[2]]]*temp_den_array[ind[2]]+x_shower[temp_det_array[ind[3]]]*temp_den_array[ind[3]])/temp_total_den
    y_core=(y_shower[temp_det_array[ind[0]]]*temp_den_array[ind[0]]+y_shower[temp_det_array[ind[1]]]*temp_den_array[ind[1]]+y_shower[temp_det_array[ind[2]]]*temp_den_array[ind[2]]+y_shower[temp_det_array[ind[3]]]*temp_den_array[ind[3]])/temp_total_den

    # load in data
    detect_cord_x = np.loadtxt('all_coords.txt',usecols = 0, max_rows = len(detectors))
    detect_cord_y = np.loadtxt('all_coords.txt',usecols = 1, max_rows = len(detectors))
    detect_coords = np.loadtxt('all_coords.txt', max_rows = len(detectors))
    for i in np.arange(len(detectors)):
        detect_dens.append(detectors[i].density)
    detect_dens=np.array(detect_dens)
    n_iters = 5

    # fit NKG with core position
    xc_fit, yc_fit,  r_M_fit, N_ch_fit, fit_xc_err, fit_yc_err, fit_r_M_err, fit_N_ch_err = minimize_core_fit(detect_coords, detect_dens, x_core, y_core, theta, phi, LORA.rM, LORA.Nch, n_iters, detect_cord_x, detect_cord_y)
    chi, final_r_M_fit, final_N_ch_fit, final_fit_r_M_err, final_fit_N_ch_err = Minimize_log_fit(detect_coords, detect_dens, xc_fit, yc_fit, theta, phi, r_M_fit, N_ch_fit, n_iters)
    
    # do atm correction 

    atm_data=read_attenuation()
    size_theta=np.zeros([30])
    if len(atm_data)!=11:
        print('no atm corrections')
    else:
        f_no=len(atm_data)
        f_int=atm_data.T[0]
        f_logN_Ref=atm_data.T[1]
        err1=atm_data.T[2]
        f_X0=atm_data.T[3]
        err2=atm_data.T[4]
        f_lamb=atm_data.T[5]
        err3=atm_data.T[6]
        Lambda0=0
        err_Lambda0=0
        for k in np.arange(len(f_int)):
            size_theta[k]=f_logN_Ref[k]-(f_X0[k]/f_lamb[k])*(1/np.cos(theta)-1/np.cos(np.pi*LORA.Ref_angle/180.0))*0.4342944819 #log10(size) for attenuation curve kk at zenith angle 'theta'
            #print size_theta[k]
            
            if np.log10(final_N_ch_fit) >= size_theta[1]: ## typo?   1->k?
                Lambda0=f_lamb[k]    #  //Extropolation
                err_Lambda0=err3[k]
                break
            
            elif np.log10(final_N_ch_fit) >= size_theta[k]:
                Lambda0=(f_lamb[k]*(np.log10(final_N_ch_fit)-size_theta[k-1])+f_lamb[k-1]*(size_theta[k]-np.log10(final_N_ch_fit)))/(size_theta[k]-size_theta[k-1]) #  //Interpolation
                err_Lambda0=(err3[k]*(np.log10(final_N_ch_fit)-size_theta[k-1])+err3[k-1]*(size_theta[k]-np.log10(final_N_ch_fit)))/(size_theta[k]-size_theta[k-1]) #  //Interpolation
         
                break
    ####

    size_theta[f_no-1]=f_logN_Ref[f_no-1]-(f_X0[f_no-1]/f_lamb[f_no-1])*(1/np.cos(theta)-1/np.cos(np.pi*LORA.Ref_angle/180))*0.4342944819

    if np.log10(final_N_ch_fit) >= size_theta[f_no-1]:
        Lambda0=f_lamb[f_no-1] #; //Extropolation
        err_Lambda0=err3[f_no-1]
        print('Extrapolation: k={0}  s_theta={1}  s2={2}  lamb={3}'.format(k,size_theta[k],np.log10(final_N_ch_fit),Lambda0))

    log_size_Ref=np.log10(final_N_ch_fit)+(LORA.X0/Lambda0)*(1/np.cos(theta)-1/np.cos(np.pi*LORA.Ref_angle/180.0))*0.4342944819# ; //log10()
    size_Ref=np.power(10,log_size_Ref)
    err_size_Ref=np.sqrt(np.power(final_fit_N_ch_err/final_N_ch_fit,2)+np.power(np.log(10)*LORA.X0*(1/np.cos(theta)-1/np.cos(np.pi*LORA.Ref_angle/180.0))*0.4342944819,2)*np.power(err_Lambda0/np.power(Lambda0,2),2))*size_Ref
    
   
    #energy_Ref=pow(size_Ref,par_b)*pow(10,par_a)*pow(10,-6) ; //Energy(PeV) at Ref_angle: Formula from KASCADE simulation (2008)
    #err_energy_Ref=sqrt(pow(log(10)*err_a,2)+pow(log(size_Ref)*err_b,2)+pow(par_b*err_size_Ref/size_Ref,2))*energy_Ref ;    //error on energy at Ref_angle
    energy_Ref=np.power(size_Ref,LORA.par_b)*np.power(10.0,LORA.par_a)*np.power(10.0,-6.0) #; #//Energy(PeV) at Ref_angle: Formula from KASCADE simulation (2008)
    err_energy_Ref=np.sqrt(np.power(np.log(10)*LORA.err_a,2)+np.power(np.log10(size_Ref)*LORA.err_b,2)+np.power(LORA.par_b*err_size_Ref/size_Ref,2))*energy_Ref ##;    //error on energy at
    

    #calculate energy
    energy=np.power(final_N_ch_fit,LORA.par_b)*np.power(10.0,LORA.par_a)*np.power(10.0,-6.0)# ; //Energy (PeV): Formula from KASCADE simulation (2008)
    err_energy=np.sqrt(np.power(np.log(10)*LORA.err_a,2)+np.power(np.log10(final_N_ch_fit)*LORA.err_b,2)+np.power(LORA.par_b*final_fit_N_ch_err/final_N_ch_fit,2))*energy #;    //error on energy

    #find first time stamp
    min_utc=1e10
    min_nsec=1e10
    for i in np.arange(len(detectors)):
        if detectors[i].gps>1 and detectors[i].cal_time>1:
            if detectors[i].gps<min_utc:
                min_utc=detectors[i].gps
            if detectors[i].cal_time<min_nsec:
                min_nsec=detectors[i].cal_time

    #temproarry fix
    corr_coef_xy = 0

    # assign event values
    event.x_core=xc_fit
    event.y_core=yc_fit
    event.x_core_err=fit_xc_err
    event.y_core_err=fit_yc_err
    event.z_core=0
    event.UTC_min=min_utc
    event.nsec_min=min_nsec
    event.energy= energy
    event.energy_err=err_energy
    event.Rm=final_r_M_fit
    event.rM_err=final_fit_r_M_err
    event.Ne=final_N_ch_fit
    event.Ne_err=final_fit_N_ch_err
    event.CorCoef_xy=corr_coef_xy
    event.Ne_RefA=size_Ref
    event.NeErr_RefA=err_size_Ref
    event.Energy_RefA=energy_Ref
    event.EnergyErr_RefA=err_energy_Ref


def fit_arrival_direction_trial(detectors, event):
    
    cdts = []
    for i in np.arange(len(detectors)):
        cdts.append(detectors[i].cdt)
    cdts = np.array(cdts)
    print(cdts)
    detector_positions = np.loadtxt('all_coords.txt', max_rows = len(detectors))
    positive_val = cdts > 0
    x_positions = detector_positions[positive_val, 0]
    y_positions = detector_positions[positive_val, 1]
    cdts_positive = cdts[positive_val]
    
    xy_data = np.vstack((x_positions, y_positions)).T

    initial_guess = [0, 0, 0]

    params, covariance = curve_fit(model_function, xy_data.T, cdts_positive, p0=initial_guess)
    t0, l, m = params
    print(f"Fitted parameters: t0={t0}, l={l}, m={m}")

    # Calculate derived quantities
    n = np.sqrt(1.0 - (l**2 + m**2)) if (l**2 + m**2) < 1 else 0
    theta = np.arcsin(np.sqrt(l**2 + m**2)) * (180.0 / np.pi)  
    phi = np.arccos(m/np.sqrt(l**2 + m**2)) * (180.0 / np.pi)  

    if l < 0:
        phi = 360.0 - phi

    # Calculate uncertainties from the covariance matrix
    err_l = np.sqrt(covariance[1, 1])
    err_m = np.sqrt(covariance[2, 2])

    err_theta = np.sqrt((l * err_l / n)**2 + (m * err_m / n)**2) / np.tan(theta * np.pi / 180.0)
    err_theta *= (180.0 / np.pi)
    
    df_f2 = ((l**2 / (m * (l**2 + m**2)))**2 * err_m**2 +
                (l / (l**2 + m**2))**2 * err_l**2)
    err_phi = np.sqrt(df_f2 / np.tan(phi * np.pi / 180.0)**2)
    err_phi *= (180.0 / np.pi)

    if math.isnan(theta)==True or math.isnan(phi)==True:
        event.fit_theta=event.theta
        event.fit_phi=event.phi
        event.fit_elevation=90.0-theta
        event.fit_theta_err=err_theta
        event.fit_elevation_err=err_theta
        event.fit_phi_err=err_phi
        print('fitting doesn\'t work')

    else:
        event.fit_theta=theta
        event.fit_phi=phi
        event.fit_elevation=90.0-theta
        event.fit_theta_err=err_theta
        event.fit_elevation_err=err_theta
        event.fit_phi_err=err_phi
    
    print('fit    theta: {0:.2f}  phi: {1:.2f} '.format(theta,phi,err_theta,err_phi))