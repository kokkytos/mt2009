import pandas as pd
import os, math
import sys


pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 10000)
pd.set_option('display.precision', 10) 
#os.chdir("/media/leonidas/Hitachi/daily_viirs/2017_packed/mt2009/VIIRS_DNB") #todo:set relative path of current file

#date_time = 201001010000


def mt2009(date_obj):
    print(date_obj)
    date_time=int(date_obj)
    mean_earthsun_dist=149598022.6071
    mean_earthmoon_dist=384400.0
    radius_earth=6378.140

    print("------------------------------------------------------------------------")
    print("Datetime = {}".format(date_time))
    
    test_year = int(date_time/1.0e+08)
    
    
    #def lun_irrad_scl(date_time):
    
    
    
    if (test_year < 2010) or (test_year > 2030):
       print('Illegal Year (Valid for 2010-2030): ',test_year)
       exit()
       

       
    #Shorten date/time to match table and separate minutes out
    egroup = int(date_time/100)
    mins = date_time - egroup*100 
    
    path = 'DIST_2010-2030.dat'
    
    col_specification = [(1, 12), (21, 29), (33, 43), (50, 57)]
    df = pd.read_fwf(path, colspecs=col_specification, header=None)
    df.columns = ['group1','phase_angle1','earthsun_dist1','earthmoon_dist1']
    #print df
        
    if mins == 0:
        exact = True
        mydata = df.loc[df['group1'] == egroup]
    
        group1=mydata['group1']
        phase_angle1=mydata['phase_angle1'].values[0]
        earthsun_dist1=mydata['earthsun_dist1'].values[0]
        earthmoon_dist1=mydata['earthmoon_dist1'].values[0]
    else:
        exact = False
        mydata = df.loc[df['group1'] == egroup]
        idx=df.loc[df['group1'] == egroup].index.values.astype(int)[0]
        group1=mydata['group1']
        phase_angle1=mydata['phase_angle1'].values[0]
        earthsun_dist1=mydata['earthsun_dist1'].values[0]
        earthmoon_dist1=mydata['earthmoon_dist1'].values[0]
        
        #nextegroup=egroup+1
        mydata2 = df.iloc[[idx+1]]
        group2=mydata2['group1']
        phase_angle2=mydata2['phase_angle1'].values[0]
        earthsun_dist2=mydata2['earthsun_dist1'].values[0]
        earthmoon_dist2=mydata2['earthmoon_dist1'].values[0]
    
    
    
    if exact:
      group = group1
      current_phase_angle = phase_angle1
      current_earthsun_dist = earthsun_dist1
      current_earthmoon_dist = earthmoon_dist1
      
    else:
      frac = mins/60.0
      delta = phase_angle1 - phase_angle2
      current_phase_angle = phase_angle1 - (delta * frac)
      delta = earthsun_dist1 - earthsun_dist2
      current_earthsun_dist = earthsun_dist1 - (delta * frac)
      delta = earthmoon_dist1 - earthmoon_dist2
      current_earthmoon_dist = earthmoon_dist1 - (delta * frac)
      
      
    
    #Open lunar irradiance file and read header lines
    path = 'lunar_irrad_1AU_MeanME_CONVOLVED_DNB.dat'
    
    col_specification = [(6, 13), (17, 27)]
    df2 = pd.read_fwf(path, colspecs=col_specification, header=None, skiprows=3)
    df2.columns = ['phase','lunar_irrad']
    #with pd.option_context('display.precision', 10):
        #print df2
        
    
    #Determine convolved irradiance from table
    
    
    search = True
    
    for index, row in df2.iterrows():
        if not search:
            break
        phase = row['phase']
        lunar_irrad = row['lunar_irrad']
        if  current_phase_angle > phase:
            phase_prev = phase
            lunar_irrad_prev = lunar_irrad
        else:      
          frac = (phase - current_phase_angle) / (phase - phase_prev)
          delta = lunar_irrad - lunar_irrad_prev
          lunar_irrad_dnb_interp = lunar_irrad - (delta * frac)
          search = False
    
    #Compute scaling factor
    cos_phase_angle = math.cos(math.radians(current_phase_angle))
    T1 = (mean_earthsun_dist**2.0) + (mean_earthmoon_dist**2.0) + (2.0*mean_earthmoon_dist*mean_earthsun_dist*cos_phase_angle)
    T2 = current_earthsun_dist**2.0 + current_earthmoon_dist**2.0 + 2.0*current_earthmoon_dist*current_earthsun_dist*cos_phase_angle
    T3 = ((mean_earthmoon_dist-radius_earth) / (current_earthmoon_dist - radius_earth))**2.0
            
    SCALE_FACTOR = (T1/T2)*T3
    
    #Compute lunar irradiance
    
    lun_irrad_scl = lunar_irrad_dnb_interp * SCALE_FACTOR
    
    print ("Current Earth/Moon Distance (km) and percent of mean = %s %s" %(current_earthmoon_dist, 100.0*current_earthmoon_dist/mean_earthmoon_dist))
    print ("Current Earth/Sun Distance (km) and percent of mean = %s %s" %(current_earthsun_dist,100.0*current_earthsun_dist/mean_earthsun_dist))
    print ("Current lunar phase angle (degrees) = %s"  %current_phase_angle)
    print ("SCALE_FACTOR = %s" %SCALE_FACTOR)
    
    print ("Lunar Irradiance (mW/m^2-micron):%s" %lun_irrad_scl)
    
    return {'current_phase_angle':current_phase_angle, 'lun_irrad_scl':lun_irrad_scl}

if __name__ == '__main__':
    date_obj  = str(sys.argv[1])
    mt2009(date_obj)
