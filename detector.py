import numpy as np
import LORAparameters as LORA

cable_delay=np.asarray([0,27,-49,-296,(-206-73+18),(-726-73+18),(-503-73+18),(-732-73+18),(-765),(-611+13),(-850+4),(106),(-116+5),(-314+0),(-528+34),(-499+20),(-1601+26),(-986+15),(-1198+20),(-1013+6)])#,(-4000),(-4000),(-1500),(-1500),(-4000),(-4000),(-1500),(-1500),(-4000),(-4000),(-1500),(-1500),(-4000),(-4000),(-1500),(-1500),(-4000),(-4000),(-1500),(-1500)])




class Detector:
    def __init__(self, name):
        self.name = name

    x_cord=0
    y_cord=0
    z_cord=0

    B_min=0
    B_max=0
    fArea12=0
    err_fArea12=0
    fArea23=0
    err_fArea23=0
    gain=0
    number=0
    
    #information stored from root file
    trigg_pattern=0
    nsec=0
    total_counts=0
    total_counts_backgroundV2=0

    counts=np.zeros([LORA.nTrace])
    trigg_condition=0
    ctd=0
    pulse_width=0
    gps=0
    pulse_height=0
    trig=0
    # calculated times
    cal_time=0
    final_event_time=0
    cdt=-1
    sec_mean=0
    sec_sigma=0
     
    
    #information calculated here
    trace_background=0
    trace_rms=0
    trace_int_counts=0
    threshold_time=0   # unit of ns
    event_time_stamp=0  #cal timestamp
    threshold=0
    corrected_peak=0
    function=0
    peak=0
    correted_threshold=0
    
    density=0
    err_density=0

class Lasa:
    def __init__(self, name):
        self.name = name
    # 3 for timestamp of event, 1 second after, 2 seconds after
    GPS_time_stamp=np.zeros([3])
    number=0
    CTP=np.zeros([3])
    sync=np.zeros([3])
    quant=np.zeros([3])
    YMD=np.zeros([3])
    
    sec_flag=0
    
    GPS_time_stamp_M=np.zeros([3])
    number_M=0
    CTP_M=np.zeros([3])
    sync_M=np.zeros([3])
    quant_M=np.zeros([3])
    YMD_M=np.zeros([3])


    GPS_time_stamp_S=np.zeros([3])
    number_S=0
    CTP_S=np.zeros([3])
    sync_S=np.zeros([3])
    quant_S=np.zeros([3])
    YMD_S=np.zeros([3])


def load_event_information(info,detectors):
    for i in np.arange(len(detectors)):
        detectors[i].number=i+1
        detectors[i].trigg_pattern=info[i]['trigg_pattern']
        detectors[i].nsec=info[i]['nsec']
        detectors[i].counts=info[i]['counts']
        detectors[i].trigg_condition=info[i]['trigg_condition']
        detectors[i].ctd=info[i]['ctd']
        detectors[i].pulse_width=info[i]['pulse_width']
        detectors[i].gps=info[i]['gps']
        detectors[i].total_counts=info[i]['total_counts']

        detectors[i].pulse_height=info[i]['pulse_height']

    for i in np.arange(5):
        m1=int(detectors[i*4].trigg_pattern)>>0&1
        m2=int(detectors[i*4].trigg_pattern)>>2&1
        s1=int(detectors[i*4].trigg_pattern)>>4&1
        s2=int(detectors[i*4].trigg_pattern)>>6&1

        detectors[i*4+0].trig=m1
        detectors[i*4+1].trig=m2
        detectors[i*4+2].trig=s1
        detectors[i*4+3].trig=s2
    
    if len(detectors)>20:
        for i in np.arange(5,10):
            detectors[i*4+0].trig = int(detectors[i*4].trigg_pattern)>>8&1
            detectors[i*4+1].trig = int(detectors[i*4].trigg_pattern)>>9&1
            detectors[i*4+2].trig = int(detectors[i*4].trigg_pattern)>>10&1
            detectors[i*4+3].trig = int(detectors[i*4].trigg_pattern)>>11&1
        


def load_sec_information(info0,info1,info2,lasas,V):
    if V=='V1':
        for i in np.arange(len(lasas)):
            lasas[i].number=i+1

            if (info2[i]['GPS_time_stamp']!=info1[i]['GPS_time_stamp']+1) and (info1[i]['GPS_time_stamp']!=info0[i]['GPS_time_stamp']+1):
                lasas[i].sec_flag=1
      
            lasas[i].CTP=np.asarray([info0[i]['CTP'],info1[i]['CTP'],info2[i]['CTP']])
            lasas[i].GPS_time_stamp=np.asarray([info0[i]['GPS_time_stamp'],info1[i]['GPS_time_stamp'],info2[i]['GPS_time_stamp']])
            lasas[i].sync=np.asarray([info0[i]['sync'],info1[i]['sync'],info2[i]['sync']])
            lasas[i].quant=np.asarray([info0[i]['quant'],info1[i]['quant'],info2[i]['quant']])
            lasas[i].YMD=np.asarray([info0[i]['YMD'],info1[i]['YMD'],info2[i]['YMD']])
           
    if V=='V2':
     
          for i in np.arange(5):
            
                lasas[i].number=i+1

                if info0[i]['CTP_S']==1 or info1[i]['CTP_S']==1 or info2[i]['CTP_S']==1:
                    lasas[i].sec_flag=1
                    
                    
        
                if (info2[i]['GPS_time_stamp_M']!=info1[i]['GPS_time_stamp_M']+1) and (info1[i]['GPS_time_stamp_M']!=info0[i]['GPS_time_stamp_M']+1):
                    lasas[i].sec_flag=1
          
                lasas[i].CTP_M=np.asarray([info0[i]['CTP_M'],info1[i]['CTP_M'],info2[i]['CTP_M']])
                lasas[i].CTP_S=np.asarray([info0[i]['CTP_S'],info1[i]['CTP_S'],info2[i]['CTP_S']])

                lasas[i].GPS_time_stamp_M=np.asarray([info0[i]['GPS_time_stamp_M'],info1[i]['GPS_time_stamp_M'],info2[i]['GPS_time_stamp_M']])
                lasas[i].GPS_time_stamp_S=np.asarray([info0[i]['GPS_time_stamp_S'],info1[i]['GPS_time_stamp_S'],info2[i]['GPS_time_stamp_S']])
                
                lasas[i].sync_M=np.asarray([info0[i]['sync_M'],info1[i]['sync_M'],info2[i]['sync_M']])
                lasas[i].sync_S=np.asarray([info0[i]['sync_S'],info1[i]['sync_S'],info2[i]['sync_S']])

                lasas[i].quant_M=np.asarray([info0[i]['quant_M'],info1[i]['quant_M'],info2[i]['quant_M']])
                lasas[i].quant_S=np.asarray([info0[i]['quant_S'],info1[i]['quant_S'],info2[i]['quant_S']])

                lasas[i].YMD_M=np.asarray([info0[i]['YMD_M'],info1[i]['YMD_M'],info2[i]['YMD_M']])
                lasas[i].YMD_S=np.asarray([info0[i]['YMD_S'],info1[i]['YMD_S'],info2[i]['YMD_S']])
          '''
          for i in np.arange(5,10):
               
           
            
                lasas[i].number=i+1

                if info0[i]['CTP']==1 or info1[i]['CTP']==1 or info2[i]['CTP']==1:
                    lasas[i].sec_flag=1
                
        
                if (info2[i]['GPS_time_stamp']!=info1[i]['GPS_time_stamp']+1) and (info1[i]['GPS_time_stamp']!=info0[i]['GPS_time_stamp']+1):
                    lasas[i].sec_flag=1
                lasas[i].CTP=np.asarray([info0[i]['CTP'],info1[i]['CTP'],info2[i]['CTP']])
                lasas[i].GPS_time_stamp=np.asarray([info0[i]['GPS_time_stamp'],info1[i]['GPS_time_stamp'],info2[i]['GPS_time_stamp']])
                lasas[i].sync_M=np.asarray([info0[i]['sync'],info1[i]['sync'],info2[i]['sync']])
                lasas[i].quant_M=np.asarray([info0[i]['quant'],info1[i]['quant'],info2[i]['quant']])
                lasas[i].YMD=np.asarray([info0[i]['YMD'],info1[i]['YMD'],info2[i]['YMD']])

          '''

def load_sec_informationV2(info0,info1,info2,lasas,V):
    if V=='V1':
        for i in np.arange(len(lasas)):
            lasas[i].number=i+1

            if (info2[i]['GPS_time_stamp']!=info1[i]['GPS_time_stamp']+1) and (info1[i]['GPS_time_stamp']!=info0[i]['GPS_time_stamp']+1):
                lasas[i].sec_flag=1
      
            lasas[i].CTP=np.asarray([info0[i]['CTP'],info1[i]['CTP'],info2[i]['CTP']])
            lasas[i].GPS_time_stamp=np.asarray([info0[i]['GPS_time_stamp'],info1[i]['GPS_time_stamp'],info2[i]['GPS_time_stamp']])
            lasas[i].sync=np.asarray([info0[i]['sync'],info1[i]['sync'],info2[i]['sync']])
            lasas[i].quant=np.asarray([info0[i]['quant'],info1[i]['quant'],info2[i]['quant']])
            lasas[i].YMD=np.asarray([info0[i]['YMD'],info1[i]['YMD'],info2[i]['YMD']])
           
    if V=='V2':
    

          for i in np.arange(10):
               
           
            
                lasas[i].number=i+1

                if info0[i]['CTP']==1 or info1[i]['CTP']==1 or info2[i]['CTP']==1:
                    lasas[i].sec_flag=1
                
        
                if (info2[i]['GPS_time_stamp']!=info1[i]['GPS_time_stamp']+1) and (info1[i]['GPS_time_stamp']!=info0[i]['GPS_time_stamp']+1):
                    lasas[i].sec_flag=1
                lasas[i].CTP=np.asarray([info0[i]['CTP'],info1[i]['CTP'],info2[i]['CTP']])
                lasas[i].GPS_time_stamp=np.asarray([info0[i]['GPS_time_stamp'],info1[i]['GPS_time_stamp'],info2[i]['GPS_time_stamp']])
                lasas[i].sync_M=np.asarray([info0[i]['sync'],info1[i]['sync'],info2[i]['sync']])
                lasas[i].quant_M=np.asarray([info0[i]['quant'],info1[i]['quant'],info2[i]['quant']])
                #lasas[i].YMD=np.asarray([info0[i]['YMD'],info1[i]['YMD'],info2[i]['YMD']])

        





def load_positions(detectors):
    #file=open(LORA.det_cord_file,'r')
    file=LORA.det_cord_file


    cordinates=np.genfromtxt(file,usecols=(0, 1, 2))
    #file.close()
    for i in np.arange(len(detectors)):
        detectors[i].x_cord=cordinates[i][0]
        detectors[i].y_cord=cordinates[i][1]
        detectors[i].z_cord=cordinates[i][2]



def load_signal(detectors):
    
    #file=open(LORA.signal_retrieve_file,'r')
    file=LORA.signal_retrieve_file

    info=np.genfromtxt(file,skip_header=8,usecols=(1,2,3,4,5,6))
    #file.close()
    for i in np.arange(len(detectors)):
        detectors[i].fArea12=info[i][0]
        detectors[i].err_fArea12=info[i][1]
        detectors[i].fArea23=info[i][2]
        detectors[i].err_fArea23=info[i][3]
        detectors[i].B_min=info[i][4]
        detectors[i].B_ax=info[i][5]

def load_gain(detectors):
    
    #file=open(LORA.gain_cal_file,'r')
    file=LORA.gain_cal_file

    info=np.genfromtxt(file,usecols=(1))
    #file.close()
    for i in np.arange(len(detectors)):
        detectors[i].gain=info[i]



def load_log_information(info,detectors):
    #print info[0].keys()
    for i in np.arange(20):
        detectors[i].threshold=info[i]['threshold']/0.48
    if len(detectors)>20:
        for i in np.arange(20,40):
            detectors[i].threshold=info[i]['threshold']
            
def load_noise_information(info,detectors):
    #print info[0].keys()
    for i in np.arange(len(detectors)):
        detectors[i].sec_mean=info[i]['mean']
        detectors[i].sec_sigma=info[i]['sigma']


