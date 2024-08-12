
import numpy as np
from def_function_final_version import fltr_day, NPZ_model_stratified_euphoticlight_oneND_concentration
import matplotlib.dates as mdates
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator



#load MLD
df1      =   pd.read_csv('interpolated_mld_daily.csv')
MLD      =   df1['Averaged MLD'].values
MLZ      =   MLD[351:11674]
date_new =   df1['date'].values#MLD date 
date_mld =   pd.to_datetime(date_new[351:11674],format='%d/%m/%Y') #have a check 
lf_MLZ   =   fltr_day(MLZ, 30)


#input the parameters
df2      = pd.read_csv('parameter_input.csv')
parame   = df2['default']

year=np.arange(1991,2022,1)
Lat        = 31    #Latitude (location)

#Define initial conditions and parameters 
# Solar Constant (W/m^2), atmospheric attenuation, PAR fraction (Brock 1981)

SolarK     = parame[0]
AtmAtt     = parame[1]
#AtmAtt     = parame[1]*0.9
ParFrac    = parame[2]
#ParFrac    = parame[2]*0.9


#surface paramters 
# Light (water extinction, phyto extinction), from Fasham et al. (1990);
Kdw        = parame[3]    # Water extinction
Kdp        = parame[4]   # Phytoplankton extinction (concentration specific)
#Kdp        = parame[4]*0.9
Ks         = parame[5]    # 1/2 saturation constant (I think)
#Ks         = parame[5]*0.9
Nd         = parame[6]   # Deep nutrient concentrations
#Nd         = parame[6]*0.9
Nut        = parame[7] # Intital nutrient concentartion
#Nut       = parame[7]*0.9
Mprev      = parame[8]     # Mixing depth
#Mprev      = parame[8]*0.9


# Phytoplankton growth 
alpha      = parame[9] # alpha [h^-1/W m^-2] #classic value
#alpha      = parame[9]*0.9
Vmax       = parame[10] # Vmax [d^-1] #averaged value from Schartau &   Oschli 2003
#Vmax       = parame[10]*0.9
P          = parame[11]# Phytoplankton starting values
#P         = parame[11]*0.9
m          = parame[12]   # mortality [d^-1]#
#m          = parame[12]*0.9

# Zooplankton parameters

Gamma      = parame[16]   # growth efficiency 
#Gamma      = parame[16]*0.9
Z          = parame[17] # statrting value#65*C_to_N*mg_to_g*gN_to_mmolN
#Z          = parame[17]*0.9
CN         = parame[18] # death rate
#CN         = parame[18]*0.9


#subsurface paramters 

a         = parame[19]#1#1.6#2from PASQUERO ET AL 2005
#a         = parame[19]*0.9
k_e       = parame[20]#3#1.6#0.7#1
#k_e       = parame[20]*0.9
frac      = parame[21]
#frac      = parame[21]*0.9
#theta     = parame[22]#0.2
theta     = parame[22]
beta      = parame[23]#
#beta      = parame[23]*0.9
#k_cn      = parame[24]*0.9
k_cn      = parame[24]#0.2



alpha_sub = parame[25] #0.256#alpha*7
#alpha_sub = parame[25]*0.9
Vmax_sub  = parame[26]#0.27#Vmax/2 #0.6 #to make growth rate closer to 0.2 adjust alpha_sub and m_sub
#Vmax_sub  = parame[26]*0.9
a_sub     = parame[27]#1.575#a#a/4
#a_sub     = parame[27]*0.9
k_e_sub   = parame[28]#1.6#k_e#k_e/4
#k_e_sub   = parame[28]*0.9
Gamma_sub = parame[29]#Gamma0.1
#Gamma_sub = parame[29]*0.9
m_sub     = parame[30]#0.053
#m_sub     = parame[30]*0.9
CN_sub    = parame[31]#0.013
#CN_sub    = parame[31]*0.9
Kdp_sub = parame[4]
#Kdp_sub = parame[4]*0.9
P_sub   = parame[32]
#P_sub   = parame[32]*0.9
Z_sub   = parame[33]
#Z_sub   = parame[33]*0.9
lf_MLZ=fltr_day(MLZ, 30)
#lf_MLZ=0.9*fltr_day(MLZ, 30)
M_EU=parame[34]
#M_EU=parame[34]*0.9

para_sub=(alpha_sub, Kdp_sub, Vmax_sub, a_sub, k_e_sub, Gamma_sub, m_sub, CN_sub)
initial_conditions2 = (P, Nut, Z, Mprev, lf_MLZ, P_sub, Z_sub, M_EU)
parameters_surface=(SolarK, AtmAtt, ParFrac, Kdw, Kdp, Vmax, alpha, Ks, Nd, m, Gamma, CN, frac, theta, beta,k_cn, a,k_e)
Pstr, Nutstr, Zstr, Io_noon, Nut_Gro, Lig_Gro, Daystr, I, I_sub, Nd_str, Pstr_sub,  Zstr_sub, Depth_EU,Lig_Gro_sub,graz_sub_str=NPZ_model_stratified_euphoticlight_oneND_concentration(Lat, year, initial_conditions2, parameters_surface,para_sub)


t=pd.to_datetime(date_new[11674],format='%d/%m/%Y') #have a check 

darker_grey = '#696969'

fig = plt.figure(figsize=[25,20])
plt.clf()
widths = [10]
heights = [2,2,2,2,2]
spec5 = fig.add_gridspec(ncols=1, nrows=5, width_ratios=widths,
                        height_ratios=heights, wspace=0.08, hspace=0.12) 

ax_current = fig.add_subplot(spec5[0,0])
#ax_current.plot(date_mld,Pstr+Nutstr+Zstr+Nd_str+Pstr_sub+Zstr_sub, color=darker_grey)
ax_current.plot(date_mld,I, '-' , color=darker_grey)
ax_current.plot(date_mld,I_sub,'--',color=darker_grey)
#ax1.set_title('Model output at BATS', fontsize=20, color='k')
ax_current.set_ylabel('Light [W m$^{-2}$]', fontsize=15, color='k')
ax_current.yaxis.set_tick_params(labelsize=15)
ax_current.xaxis.set_tick_params(labelsize=15)
ax_current.set_xlim(date_mld.min(),date_mld.max())
ax_current.set_ylim(0,260)
major_interval = 50
ax_current.tick_params(which='major', direction='in', length=7)
ax_current.yaxis.set_major_locator(MultipleLocator(major_interval))
ax_current.yaxis.set_minor_locator(AutoMinorLocator(5))
ax_current.tick_params(which='minor', direction='in', length=4)
ax_current.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax_current.xaxis.set_major_locator(mdates.YearLocator())
plt.setp(ax_current,xticklabels=[])
ax_current.set_title('(a)', fontsize=16)

ax_current = fig.add_subplot(spec5[1,0])
ax_current.plot(date_mld,lf_MLZ, '-', color=darker_grey)
ax_current.plot(date_mld,Depth_EU, '--',color=darker_grey)
ax_current.set_ylabel('Depth [m]', fontsize=15, color='k')
ax_current.invert_yaxis() 
#ax2.plot(Daystr,fl_mld, label='MLD 3 day filtered', color='k')
ax_current.yaxis.set_tick_params(labelsize=15)
ax_current.set_xlim(date_mld.min(),date_mld.max())
ax_current.legend()
major_interval = 100
ax_current.tick_params(which='major', direction='in', length=7)
ax_current.yaxis.set_major_locator(MultipleLocator(major_interval))
ax_current.yaxis.set_minor_locator(AutoMinorLocator(10))
ax_current.tick_params(which='minor', direction='in', length=4)
ax_current.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax_current.xaxis.set_major_locator(mdates.YearLocator())
plt.setp(ax_current,xticklabels=[])
ax_current.set_title('(b)', fontsize=16)
ax_current.set_ylim(320,0)


ax_current = fig.add_subplot(spec5[2,0])
ax_current.plot(date_mld, Pstr,'-', color=darker_grey)
ax_current.plot(date_mld, Pstr_sub,'--',color=darker_grey)
ax_current.set_ylabel('Phyto [mmol m$^{-3}$]', fontsize=15, color='k')
ax_current.yaxis.set_tick_params(labelsize=15)
ax_current.set_xlim(date_mld.min(),date_mld.max())
#ax_current.set_ylim(0,200)
ax_current.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax_current.xaxis.set_major_locator(mdates.YearLocator())
plt.setp(ax_current,xticklabels=[])
ax_current.set_ylim(0,1.25)
major_interval = 0.25
ax_current.tick_params(which='major', direction='in', length=7)
ax_current.yaxis.set_major_locator(MultipleLocator(major_interval))
ax_current.yaxis.set_minor_locator(AutoMinorLocator(5))
ax_current.tick_params(which='minor', direction='in', length=4)
ax_current.set_title('(c)', fontsize=16)

ax_current = fig.add_subplot(spec5[3,0])
ax_current.plot(date_mld,Zstr,'-', color=darker_grey)
ax_current.plot(date_mld,Zstr_sub,'--',color=darker_grey)
ax_current.set_ylabel('Zoo [mmol m$^{-3}$]', fontsize=15, color='k')
ax_current.yaxis.set_tick_params(labelsize=15)
ax_current.set_xlim(date_mld.min(),date_mld.max())
ax_current.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax_current.xaxis.set_major_locator(mdates.YearLocator())
plt.setp(ax_current,xticklabels=[])
ax_current.set_ylim(0,1.25)
major_interval = 0.25
ax_current.tick_params(which='major', direction='in', length=7)
ax_current.yaxis.set_major_locator(MultipleLocator(major_interval))
ax_current.yaxis.set_minor_locator(AutoMinorLocator(5))
ax_current.tick_params(which='minor', direction='in', length=4)
ax_current.set_title('(d)', fontsize=16)

ax_current = fig.add_subplot(spec5[4,0])
ax_current.plot(date_mld,Nutstr,'-', color=darker_grey)
ax_current.plot(date_mld,Nd_str,'--', color=darker_grey)
ax_current.set_ylabel('Nutrients [mmol m$^{-3}$]', fontsize=15, color='k')
ax_current.yaxis.set_tick_params(labelsize=15)
ax_current.set_xlim(date_mld.min(),t)
ax_current.set_xlabel('Year', fontsize=16)
ax_current.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax_current.xaxis.set_major_locator(mdates.YearLocator())
ax_current.set_title('(e)', fontsize=16)
ax_current.xaxis.set_tick_params(labelsize=15)
ax_current.set_ylim(0,11)
major_interval = 2
ax_current.tick_params(which='major', direction='in', length=7)
ax_current.yaxis.set_major_locator(MultipleLocator(major_interval))
ax_current.yaxis.set_minor_locator(AutoMinorLocator(2))
ax_current.tick_params(which='minor', direction='in', length=4)
#fig.savefig('./paper/fig02.png',bbox_inches='tight', dpi=400)
 
