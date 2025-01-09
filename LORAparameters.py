nDet=4
dtV1=2.5
nLORA=20
nLASA=5
nLasaB=10
nLasaA=5

nTrace=4000

nDetA=20
nDetB=40

Min_trigg_conf=8    #Minimum no. of detectors to accept an event
Age=1.7            #Staring value of shower age parameter
rM=30            #Starting value of Moliere radius
Nch=10**6        #Starting value number of charged aprticles
Rho_cut=10000    #Upper density cut in particles/m^2: considering only detectors with density<=Rho_cut
No_Bin_R=1600    #No. of bin for lateral density
min_R= 0      #Min. radius for ,,    ,,
max_R=350    #Max    ,,    ,,    ,,

X0=1024  #Vertical atmos. thickness in g/cm^2
Ref_angle=21     #Reference zenith angle (deg.) for calculating atmos. attenuation

par_a=1.23    #Eneregy reconstruction paramter (From Kickelbick 2008, Kascade thesis)
par_b=0.95    #  ,,  ,,  ,,
err_a=0.14    #Error on "par_a"
err_b=0.02    #Error on "par_b"

#det_cord_file='/Users/kmulrey/LOFAR/LORA/LORAprocessing/newLORA_processing/LORA_software_V2/data/Detector_Cord.dat'
#gain_cal_file='/Users/kmulrey/LOFAR/LORA/LORAprocessing/newLORA_processing/LORA_software_V2/data/gain_calib.dat'
##signal_retrieve_file='/Users/kmulrey/LOFAR/LORA/LORAprocessing/newLORA_processing/LORA_software_V2/data/signal_retrive.dat'
#atm_file='/Users/kmulrey/LOFAR/LORA/LORAprocessing/newLORA_processing/LORA_software_V2/data/atmos_attenuation.dat'

det_cord_file='/vol/astro7/lofar/lvdongen/LORA_V2_final_version/all_coords.txt'
gain_cal_file='/vol/astro7/lofar/lvdongen/LORA_V2_final_version/data/gain2.dat'
signal_retrieve_file='/vol/astro7/lofar/lvdongen/LORA_V2_final_version/data/signal_retrive.dat'
atm_file='/vol/astro7/lofar/lvdongen/LORA_V2_final_version/data/atmos_attenuation.dat'

event_id_offset=1262304000

# === Background traces calculation ============================================

BG_Min=0            #t_min (nsec) for background calculation
BG_Max=400          #t_max (nsec) for background calculation
BG_No_Bin=(BG_Max-BG_Min)/2.5    #No. of bins for background calculation
Window_Open=70    #START of the signal time window [: T_peak-Window_open] nsecs
Sig_Time_Window=1000    #END of the signal time window [: T_peak+Sig_Time_Window] nsecs



BG_Min_V2=0            #t_min (nsec) for background calculation
BG_Max_V2=100          #t_max (nsec) for background calculation
BG_No_Bin_V2=(BG_Max-BG_Min)/5    #No. of bins for background calculation
Window_Open_V2=70    #START of the signal time window [: T_peak-Window_open] nsecs
Sig_Time_Window_V2=1000    #END of the signal time window [: T_peak+Sig_Time_Window] nsecs



Max_ADC_Count=4075    #//Maximum ADC count (-background) cosidered as saturated
Max_ADC_Count_V2=16000    #//Maximum ADC count (-background) cosidered as saturated

Det_Thres=4

vel_light=2.99792458e8    #velocity of light in m/sec.
Det_Area=0.9        #Detector collection area in m^2
Density_Cut=1.0     #Considering only those detectors with density greater than this.
Density_Cut_High=10000 #Upper density cut in particles/m^2: considering only detectors with density<=Rho_cut





# fit values

Age=1.7
rM=30
No_Bin_R=1600    # No. of bin for lateral density
min_R=0        # Min. radius for ,,    ,,
max_R=350        # Max    ,,    ,,    ,,

#define    X0            1024    // Vertical atmos. thickness in g/cm^2
#define    Ref_angle    21        // Reference zenith angle (deg.) for calculating atmos. attenuation

par_a=1.23    # Eneregy reconstruction paramter (From Kickelbick 2008, Kascade thesis)
par_b=0.95    #  ,,  ,,  ,,
err_a=0.14    # Error on "par_a"
err_b=0.02    # Error on "par_b"

Bin_Size_X=1  # Bin size along X-axis for core distribution in meters
Bin_Size_Y=1  # Bin size along Y-axis for ''    ''    ''    ''


Ref_angle=21.0
X0=1024    #Vertical atmos. thickness in g/cm^2


#
