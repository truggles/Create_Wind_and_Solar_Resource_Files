# update info:

# 2019-02-27: 
# (1) Added the "variable total_sunlight_tolerance". This is used to scale the maximum allowable surface radiation; 
# (2) Modify the function "cal_incidence_angles". Panel tilt is now returned to calculet the adjusted diffuse sunlight;
# (3) Add funtions to determine if year_input is a leap year or not;

# 2019-04-25:
# (1) Change the calculation for leap year

# 2019-08-27:
# (1) Change the partitions to direct and diffuse radiation; 
# (2) Change the adjustments of in-panel diffuse radiation;
# (3) Change the calculations of solar declination angle;
# (4) Change the use of Beer's law to limit direct only;
# (5) Read 2m-temperature from wind files as well.

import cdms2 as cdms,MV2 as MV, numpy as np,cdutil
import os,sys,calendar,datetime
from math import degrees as deg, radians as rad

# helpers in a different directory
sys.path.append("/home/truggles/Create_Wind_and_Solar_Resource_Files/get_regional_CF_time_series_on_linux-mac") 

from helpers import get_prefix_name

if len(sys.argv) == 1:
    year = 2016
else:
    year = sys.argv[1]

###################### control variables ###################### 

case_name = get_prefix_name(int(year), True)+str(year)
isleap = calendar.isleap(int(year))
if isleap == True:
    leap_year = 1
else:
    leap_year = 0
data_path="/lustre/scratch/leiduan/MERRA2_data/Solar/"
data_path2="/lustre/scratch/leiduan/MERRA2_data/Wind/"
data_path3="/lustre/scratch/truggles/MERRA2_data/Temp/"
lat_num = 361
lon_num = 576
max_clearness_index = 1.0
total_sunlight_tolerance = 1.0
tilt_pv_input = 0.  # in units of degree
azim_pv_input = 0.
###############################################################

if leap_year == 0:
    hour_in_years=8760
    month_days={1:31,2:28,3:31,4:30,5:31,6:30,7:31,8:31,9:30,10:31,11:30,12:31}
    month_days_array = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
elif leap_year == 1:
    hour_in_years=8784
    month_days={1:31,2:29,3:31,4:30,5:31,6:30,7:31,8:31,9:30,10:31,11:30,12:31}
    month_days_array = np.array([31,29,31,30,31,30,31,31,30,31,30,31])

fd_cam=os.listdir(data_path)
temp = MV.array(np.zeros([hour_in_years,lat_num,lon_num])); temp.id='temp'; temp.units='1'; temp.missing_value = 1e20
R_TAMB = 20.          # Reference ambient temperature (degC)
R_TMOD = 25.          # Reference module temperature (degC)
R_IRRADIANCE = 1000.  # Reference irradiance (W/m2)
HF_free = 0.035       # Free-standing module, assuming no wind
HF_inte = 0.05        # Building-integrated module
# Assume CSI solar panel here. Numbers are derived directly from Ninja's code
k_1, k_2, k_3, k_4, k_5, k_6 = -0.017162, -0.040289, -0.004681, 0.000148, 0.000169, 0.000005
degree_to_radius = np.pi / 180.
radius_to_degree = 180. / np.pi




# main function starts here
tilt_pv = np.zeros([lat_num,lon_num])+tilt_pv_input*degree_to_radius
azim_pv = np.zeros([lat_num,lon_num])+azim_pv_input*degree_to_radius
df = np.zeros([lat_num,lon_num])
count_num=1
for file in fd_cam:
        if file[:-8] == case_name and file[-4:]=='.nc4':
            f=cdms.open(data_path+file)
            month = int(file[-8:-6])
            days  = int(file[-6:-4])
            if month == 1:
                days_ord = 0 + days
                position = (0 + (days-1)) *24
            else:
                days_ord = (np.sum(month_days_array[:month-1]) + days)
                position = (np.sum(month_days_array[:month-1]) + (days-1)) *24
            print (month, days, days_ord, position)
            
            # get data
            lat = f.getAxis('lat')
            lon = f.getAxis('lon')
            swgdn_tmp = MV.filled(f('SWGDN',squeeze=1),0.)
            swtdn_tmp = MV.filled(f('SWTDN',squeeze=1),0.)
            file2 = file[:20]+'slv'+file[-16:]
            fwind = cdms.open(data_path2+file2)
            t = fwind('T2M',squeeze=1)-273.15
            fwind.close()

            if len(lat) != lat_num or len(lon) != lon_num:
                print ('lat/lon number error')
                sys.exit
            
            #kt = np.zeros([24,lat_num,lon_num])
            #kt[swtdn_tmp != 0.] = (swgdn_tmp[swtdn_tmp != 0.]/swtdn_tmp[swtdn_tmp != 0.])
            #kt[kt<0.] = 0.
            for hr_idx in range(24):
                #zenith, solar_azi, ha = cal_solar_angles(lat, lon, int(year), month, days, hr_idx, days_ord)
                #incidence_rad, panel_tilt_rad = cal_incidence_angles(zenith, solar_azi, tilt_pv, azim_pv, 'h')
                #mask1 = MV.filled(MV.masked_equal(swtdn_tmp[hr_idx]   ,0)*0+1,0)
                #mask2 = MV.filled(MV.masked_greater_equal(zenith,np.pi/2)*0+1,0)
                #cosine_zenith = np.cos(zenith)*mask1*mask2
                ##cosine_zenith[(cosine_zenith>0)&(cosine_zenith<0.087)] = 0.087 ### do we want it here?
                #cosine_incide = np.cos(incidence_rad) * mask1*mask2
                #adjust_factor_dni = replace_inf(replace_nan(cosine_incide / cosine_zenith))
                #adjust_factor_dni[(adjust_factor_dni<1.)]=1.
                #
                #maximum_index = np.argwhere( swgdn_tmp[hr_idx]== np.max(swgdn_tmp[hr_idx]) )
                #base = swgdn_tmp[hr_idx][maximum_index[0,0]][maximum_index[0,1]]/swtdn_tmp[hr_idx][maximum_index[0,0]][maximum_index[0,1]]
                #cons = swtdn_tmp[hr_idx][maximum_index[0,0]][maximum_index[0,1]]
                #potential_max_solar = np.zeros([lat_num,lon_num])
                #potential_max_solar = np.where(cosine_zenith!=0., cons*base**(1/cosine_zenith), potential_max_solar) 
                #
                #kt2 = kt[hr_idx]
                #df = np.zeros([lat_num,lon_num])
                #df = np.where( kt2 <= 0.30, 1.020 - 0.254 * kt2 + 0.0123 * cosine_zenith, df)
                #df = np.where( (kt2 > 0.30) & (kt2 < 0.78), 1.4 - 1.749 * kt2 + 0.0177 * cosine_zenith, df)
                #df = np.where( kt2 >= 0.78, 0.486 * kt2 - 0.182 * cosine_zenith, df)
                #df = np.where( kt2> 1.0, 1.0, df)
                #dhi = df * swgdn_tmp[hr_idx]
                #dni = (swgdn_tmp[hr_idx] - dhi)
                #
                #z = zenith  #in radius
                #kappa   = 1.041
                #airmass = replace_inf(replace_nan(1./cosine_zenith))
                #delta   = dhi * airmass / cons
                #eps     = replace_inf(replace_nan(swgdn_tmp[hr_idx] / dhi + kappa * (z ** 3)) / (1 + kappa * (z ** 3)))
                #F1 = np.zeros([lat_num,lon_num])
                #F1[(1.000<eps)&(eps<=1.065)] = F11[1] + F12[1] * delta[(1.000<eps)&(eps<=1.065)] + F13[1]*z[(1.000<eps)&(eps<=1.065)]
                #F1[(1.065<eps)&(eps<=1.230)] = F11[2] + F12[2] * delta[(1.065<eps)&(eps<=1.230)] + F13[2]*z[(1.065<eps)&(eps<=1.230)]
                #F1[(1.230<eps)&(eps<=1.500)] = F11[3] + F12[3] * delta[(1.230<eps)&(eps<=1.500)] + F13[3]*z[(1.230<eps)&(eps<=1.500)]
                #F1[(1.500<eps)&(eps<=1.950)] = F11[4] + F12[4] * delta[(1.500<eps)&(eps<=1.950)] + F13[4]*z[(1.500<eps)&(eps<=1.950)]
                #F1[(1.950<eps)&(eps<=2.800)] = F11[5] + F12[5] * delta[(1.950<eps)&(eps<=2.800)] + F13[5]*z[(1.950<eps)&(eps<=2.800)]
                #F1[(2.800<eps)&(eps<=4.500)] = F11[6] + F12[6] * delta[(2.800<eps)&(eps<=4.500)] + F13[6]*z[(2.800<eps)&(eps<=4.500)]
                #F1[(4.500<eps)&(eps<=6.200)] = F11[7] + F12[7] * delta[(4.500<eps)&(eps<=6.200)] + F13[7]*z[(4.500<eps)&(eps<=6.200)]
                #F1[(6.200<eps)]              = F11[8] + F12[8] * delta[(6.200<eps)]              + F13[8]*z[(6.200<eps)] 
                #F2 = np.zeros([lat_num,lon_num])
                #F2[(1.000<eps)&(eps<=1.065)] = F21[1] + F22[1] * delta[(1.000<eps)&(eps<=1.065)] + F23[1]*z[(1.000<eps)&(eps<=1.065)]
                #F2[(1.065<eps)&(eps<=1.230)] = F21[2] + F22[2] * delta[(1.065<eps)&(eps<=1.230)] + F23[2]*z[(1.065<eps)&(eps<=1.230)]
                #F2[(1.230<eps)&(eps<=1.500)] = F21[3] + F22[3] * delta[(1.230<eps)&(eps<=1.500)] + F23[3]*z[(1.230<eps)&(eps<=1.500)]
                #F2[(1.500<eps)&(eps<=1.950)] = F21[4] + F22[4] * delta[(1.500<eps)&(eps<=1.950)] + F23[4]*z[(1.500<eps)&(eps<=1.950)]
                #F2[(1.950<eps)&(eps<=2.800)] = F21[5] + F22[5] * delta[(1.950<eps)&(eps<=2.800)] + F23[5]*z[(1.950<eps)&(eps<=2.800)]
                #F2[(2.800<eps)&(eps<=4.500)] = F21[6] + F22[6] * delta[(2.800<eps)&(eps<=4.500)] + F23[6]*z[(2.800<eps)&(eps<=4.500)]
                #F2[(4.500<eps)&(eps<=6.200)] = F21[7] + F22[7] * delta[(4.500<eps)&(eps<=6.200)] + F23[7]*z[(4.500<eps)&(eps<=6.200)]
                #F2[(6.200<eps)]              = F21[8] + F22[8] * delta[(6.200<eps)]              + F23[8]*z[(6.200<eps)] 
                #
                #dni_adjust  = dni * adjust_factor_dni  
                #dhi_adjust  = dhi * ( (1-F1)*(1+np.cos(panel_tilt_rad))/2 + F2*np.sin(panel_tilt_rad) )
                #dhi_adjust2 = dhi * F1 * adjust_factor_dni
                #dni_adjust2 = np.minimum(dni_adjust+dhi_adjust2, potential_max_solar * total_sunlight_tolerance)
                #rad_adjust  = replace_nan(dni_adjust2 + dhi_adjust)

                # Now calculate solar capacity factor
                # Huld, T. et al., 2010. Mapping the performance of PV modules, effects of module type and data averaging. Solar Energy, 84(2), p.324-338. DOI: 10.1016/j.solener.2009.12.002
                #T_ = np.array((1*t[hr_idx] + HF_free * rad_adjust) - R_TMOD )
                #G_ = np.array(rad_adjust / R_IRRADIANCE)
                T_ = np.array(1*t[hr_idx])
                index = int(position+hr_idx)
                #temp[index,G_==0.] = 0.
                #temp[index,G_ >0.] = G_[G_>0.] * (1+k_1*np.log(G_[G_>0.])+k_2*(np.log(G_[G_>0.]))**2 + T_[G_>0.]*(k_3+k_4*(np.log(G_[G_>0.]))+k_5*(np.log(G_[G_>0.]))**2) + k_6*(T_[G_>0.]**2))
                temp[index] = T_

#temp[temp<0]=0
#temp[temp>1]=1

fout=cdms.open(data_path3+case_name+'_temp.nc','w')
fout.write(temp)
fout.close()
temp_annual = cdutil.averager(temp,axis=0,weights='equal')
temp_annual.id='temp_annual'
gout=cdms.open(case_name+'_temp_annual.nc','w')
gout.write(temp_annual)
gout.close()
