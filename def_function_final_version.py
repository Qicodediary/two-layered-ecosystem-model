#final version of one-layer, fixed two-layer and final two layer models

import numpy as np
import math as mt

import pandas as pd
#===============================================================================================
#
#                    Define a function to calculate the number of days in each year
#
#===============================================================================================


#This function allows you to input a year to gain the number of days in that year 
#for example: year_days(2020) should return 366 because 2020 is a leap year which has 366 days

def year_days(year):
    #to test if it's a leap year, leap years have 366 days, other years have 365 days
    # A leap year is divisible by 4
    if year % 4 == 0:
        # If year is divisible by 100, it must also be divisible by 400 to be a leap year
        if year % 100 == 0:
            Dyear = 366
        # If year is not divisible by 100, it is a leap year
        else:
            Dyear=366
    # If year is not divisible by 4, it is not a leap year
    else:
        Dyear=365
    return Dyear

def fltr_day(ts, ll):
    NN=len(ts)
    ts_out=np.full(NN, np.nan)
    LL=mt.floor(ll/2)# make average from the middle
    for n in range(NN):
        if n<LL:
            ts_out[n]=np.nanmean(ts[n:LL*2+1]) # return itself
        elif n>NN-LL:
            ts_out[n]=np.nanmean(ts[n-LL*2:n+1]) # return itself
        else:
            ts_out[n]=np.nanmean(ts[n-LL:n+LL+1])
    return ts_out

#how to use:ts_out=fltr_5day(ts,5) means filter the ts every five days.

def NPZ_model_2_dynamic_Nd_fixedbry_holling(Lat, year, initial_conditions, parameters):
    P, Nut, Z, Mprev, MLZ= initial_conditions

    # Unpack parameters
    SolarK, AtmAtt, ParFrac, Kdw, Kdp, Vmax, alpha, Ks, Nd, m, c, kappa, Gamma, CN, frac, theta, beta,k_cn,bry,a,k_e = parameters
    
    #check how many years
    NN=len(year)
    
    #define initial zero value to calculate the total days (time step) of this simulation
    LL=0
    for n in range(NN):
        LL += year_days(year[n])
    
    # Initialize arrays
    Pstr, Nutstr, Zstr, Io_noon, Nut_Gro, Lig_Gro, Daystr, I, Nd_str, cn_term, Depth_EU  = [np.full(LL, np.nan) for _ in range(11)]



# Loop through each year of year

    Ltom3=1000 #1L=0.001m3
    Day=-1
    for Yr in range(NN):
        y_day=year_days(year[Yr-1])
    # Loop through the days
        for i in range(y_day):
            Day         += 1
            Daystr[Day] = Day
        #COMPUTE THE LIGHT ENVIRONMENT
            #COMPUTE EXTINCTION COEFFICIENT Kd
            #NOTE Error in the Matlab scrip p83 of Miller and Wheeler
            #P need to be converted to Chl
            Kd         = Kdw + Kdp * P #
   
            # Compute surface PAR at latitude for times of day 
            # and dawn to noon and production rate integration down to MLD
            # Decination = angle of sun above equator on each day
            D1          = 23.45*mt.sin(mt.radians((360.*(284.+i+1)/365.)))  
            #=================================================================================
            # Angle(deg) between south (i.e. noon) and setting sun
            W1          = mt.degrees(mt.acos(-1.*(mt.tan(mt.radians(Lat)) * mt.tan(mt.radians(D1)))))
            
            #=================================================================================
            # One-half daylength, hours to noon
            #=================================================================================
            # Earth rotates 15 = degrees / hour =360/24
            L1          = W1/15.
            #=================================================================================
            # Distance of Earth from sun relative to average, a minor effect
            Rx          = 1./mt.sqrt(1.+0.033*np.cos(np.deg2rad((360.*(i+1)/365.))))
            VtofD       = L1/40.
            TofD        = 12.01-L1-VtofD
            SGr         = 0.
            I_d2n=0
            # Loop over dawn to noon in 40 steps
            for j in range(40):
                TofD   = TofD + VtofD   # for the following see Brock (1981)
                W2     = (TofD-12.)*15.
                CosZen = np.sin(np.deg2rad(D1))*np.sin(np.deg2rad(Lat)) \
                         +np.cos(np.deg2rad(D1))*np.cos(np.deg2rad(Lat)) \
                         *np.cos(np.deg2rad(W2))
                Isurf  = SolarK*CosZen/(Rx*Rx)
                Io     = Isurf*AtmAtt*ParFrac
                # Save Noon PAR for plotting
                #if j == 19:
                #    Io_noon[Day] = Io
                sm=0
                for k in range(1,int(MLZ[Day])+1):
                    Iz = Io*mt.exp(-1.*Kd*k) #LORENZEN 1972
                    sm=sm+Iz #accumulated from surface to MLZ
                    #COMPUTE Phytoplankton Growth Rate for each sum from 0 to MLZ and dawn to non             
                                # Denman-Pena function for phyto growth, scales (0-1)
                    Gr     = (1. - mt.exp(-alpha*Iz/Vmax))
                    SGr    += Gr
                I_d2n=I_d2n+sm/MLZ[Day] # calculated MLZ-averaged irradiance from dawn to noon
            I[Day]=2*I_d2n/40 # to get dawn to dusk

                    
                
        #COMPUTE Phytoplankton Growth Rate over the day (LIGHT EFFECT)
            AveGr = 2.*Vmax*SGr/(MLZ[Day]*40.)  # (2* to get dawn to dusk)
            Lig_Gro[Day] = AveGr
            
        #COMPUTE Phytoplankton Growth rate NUTRIENT EFFECT)
            NL    = Vmax*Nut/(Ks+Nut) # (NUTRIENT EFFECT)
            Nut_Gro[Day] = NL
            
            # Take minimium growth rate
            G = min(NL, AveGr)
            
        #COMPUTE mixing due to mixed layer deepening
            zeta = MLZ[Day] - Mprev if MLZ[Day] > Mprev else 0.
            
        #COMPUTE Grazing rate accounting for grazing threshold
        
            #graz = c * (P - P0) / (kappa + P - P0) if P > P0 else 0.
            graz = a*k_e*P*P/(a+k_e*P*P)
        
        #COMPUTE Phytoplankton state variable

        
            DP      =    G * P - m*P - graz*Z 
            DN      =  - G*P + theta*m*P + (1-Gamma-beta)*graz*Z + k_cn*CN*Z + (frac*MLZ[Day]+zeta)*( Nd- Nut)/MLZ[Day]
            DZ      =    graz * Gamma * Z - CN*Z

            
            Nd = Nd*(bry-Mprev)/(bry-MLZ[Day])- (frac*MLZ[Day]+zeta)*(Nd- Nut)/(bry-MLZ[Day])+ (1-k_cn)*CN*Z*MLZ[Day]/(bry-MLZ[Day]) +beta*graz*Z*MLZ[Day]/(bry-MLZ[Day]) + (1-theta)*m*P*MLZ[Day]/(bry-MLZ[Day])

            
            P        = P*Mprev/MLZ[Day] + DP
            Nut      = Nut*Mprev/MLZ[Day] +DN
            Z        = Z*Mprev/MLZ[Day]+DZ


            Mprev = MLZ[Day]
            
            Pstr[Day], Nutstr[Day], Zstr[Day], Nd_str[Day],  = P*MLZ[Day] , Nut*MLZ[Day], Z*MLZ[Day], Nd*(bry-MLZ[Day])
            
    return Pstr, Nutstr, Zstr, Io_noon, Nut_Gro, Lig_Gro, Daystr, I, Nd_str





def NPZ_model_stratified_euphoticlight_oneND(Lat, year, initial_conditions, parameters, parameters_sub):
    #this function is the final verison of permenant stratifed ocean model with two-layer (with uphotic zone always deeper than MLD)
    P, Nut, Z, Mprev, MLZ, P_sub, Z_sub, M_EU= initial_conditions

    # Unpack parameters
    SolarK, AtmAtt, ParFrac, Kdw, Kdp, Vmax, alpha, Ks, Nd, m, Gamma, CN, frac, theta, beta,k_cn, a,k_e = parameters
    alpha_sub, Kdp_sub, Vmax_sub, a_sub, k_e_sub, Gamma_sub, m_sub, CN_sub= parameters_sub
    
    #check how many years
    NN=len(year)
    
    #define initial zero value to calculate the total days (time step) of this simulation
    LL=0
    for n in range(NN):
        LL += year_days(year[n])
    
    # Initialize arrays
    Pstr, Nutstr, Zstr, Io_noon, Nut_Gro, Lig_Gro, Daystr, I, Nd_str, cn_term, Depth_EU, I_sub, Lig_Gro_sub,Pstr_sub, Nut_sub_str, Zstr_sub,Vmax_sub_str,graz_sub_str= [np.full(LL, np.nan) for _ in range(18)]



# Loop through each year of year

    Day=-1
    for Yr in range(NN):
        y_day=year_days(year[Yr-1])
    # Loop through the days
        for i in range(y_day):
            Day         += 1
            Daystr[Day] = Day
            #=================================================================================
            #                    COMPUTE THE LIGHT ENVIRONMENT
            #COMPUTE EXTINCTION COEFFICIENT Kd
            #=================================================================================
            Kd         = Kdw + Kdp * P # first layer
            #Depth_EU[Day] = 300  ###13.8
            #Depth_EU[Day] = 9.2/Kd+50
            #if Depth_EU[Day]< MLZ[Day]:
            #    Depth_EU[Day]=MLZ[Day]+100

            Depth_EU[Day] = 13.8/Kd


           
            #Kd_sub = Kdw
            Kd_sub     = Kdw + Kdp_sub*P_sub #second layer

            # Compute surface PAR at latitude for times of day 
            # and dawn to noon and production rate integration down to MLD
            # Decination = angle of sun above equator on each day
            D1          = 23.45*mt.sin(mt.radians((360.*(284.+i+1)/365.)))  
            #=================================================================================
            # Angle(deg) between south (i.e. noon) and setting sun
            W1          = mt.degrees(mt.acos(-1.*(mt.tan(mt.radians(Lat)) * mt.tan(mt.radians(D1)))))
            
            #=================================================================================
            # One-half daylength, hours to noon
            #=================================================================================
            # Earth rotates 15 = degrees / hour =360/24
            L1          = W1/15.
            #=================================================================================
            # Distance of Earth from sun relative to average, a minor effect
            Rx          = 1./mt.sqrt(1.+0.033*np.cos(np.deg2rad((360.*(i+1)/365.))))
            VtofD       = L1/40.
            TofD        = 12.01-L1-VtofD
          
            
            SGr         = 0.
            SGr2        = 0
            SGr_sub     = 0
            cn_sgr2     = 0 
            cn_sgr_sub     = 0 
            
            I_d2n=0 #layer 1 initial value
            I_d2n_sub=0 #laer 2 initial value 
            I_mlz=0 
            # Loop over dawn to noon in 40 steps
            for j in range(40):
                TofD   = TofD + VtofD   # for the following see Brock (1981)
                W2     = (TofD-12.)*15.
                CosZen = np.sin(np.deg2rad(D1))*np.sin(np.deg2rad(Lat)) \
                         +np.cos(np.deg2rad(D1))*np.cos(np.deg2rad(Lat)) \
                         *np.cos(np.deg2rad(W2))
                Isurf  = SolarK*CosZen/(Rx*Rx)
                Io     = Isurf*AtmAtt*ParFrac
                #================================================================================
                #
                #                    first layer 
                #
                #=================================================================================
                sm=0 #first layer initial


                for k in range(1,int(MLZ[Day])+1):
                    
                    Iz = Io*mt.exp(-1.*Kd*k) #LORENZEN 1972
                    sm=sm+Iz #accumulated from surface to MLZ
                    #COMPUTE Phytoplankton Growth Rate for each sum from 0 to MLZ and dawn to non             
                                # Denman-Pena function for phyto growth, scales (0-1)
                    Gr     = (1. - mt.exp(-alpha*Iz/Vmax))
                    SGr    += Gr
                I_d2n=I_d2n+sm/int(MLZ[Day]) # calculated MLZ-averaged irradiance from dawn to noon
                #=================================================================================
                #                    second layer: 2cases: bd_layer <mld and when bd_layer >mld
                #
                #=================================================================================

                sm_sub=0
                I_m=Io*mt.exp(-1.*Kd*int(MLZ[Day]))  #subsurface light caluclation 
                for k2 in range(int(MLZ[Day])+1, int(Depth_EU[Day])+1):
                    Iz2 = I_m*np.exp(-1.*Kd_sub*k2) #LORENZEN 1972
                    sm_sub=sm_sub+Iz2 #accumulated from surface to MLZ
                    #COMPUTE Phytoplankton Growth Rate for each sum from 0 to MLZ and dawn to non             
                                # Denman-Pena function for phyto growth, scales (0-1)
                    Gr_sub     = (1. - mt.exp(-alpha_sub*Iz2/Vmax_sub))
                    SGr_sub    += Gr_sub

                I_d2n_sub=I_d2n_sub+sm_sub/(Depth_EU[Day]-MLZ[Day]) # calculated MLZ-averaged irradiance from dawn to noon
                cn_sgr_sub +=1
            #out of loop 40: 
            
            #out of loop 40: 
            
            #================================================================================
            #
            #                    first layer : averaged light
            #
            #=================================================================================   
        
            I[Day]=2*I_d2n/40 # to get dawn to dusk for surface layer light 

            #=================================================================================
            #                    second layer : averaged light
            #
            #=================================================================================

            I_sub[Day]=2*I_d2n_sub/cn_sgr_sub #layer 2: case2 



       
            #================================================================================
            #
            #                    first layer 
            #
            #=================================================================================

            #COMPUTE Phytoplankton Growth Rate over the day (LIGHT EFFECT)
            AveGr = 2.*Vmax*SGr/(MLZ[Day]*40)  # (2* to get dawn to dusk)
            Lig_Gro[Day] = AveGr

            #COMPUTE Phytoplankton Growth rate NUTRIENT EFFECT)
            NL    = Vmax*Nut/(Ks+Nut) # (NUTRIENT EFFECT)
            Nut_Gro[Day] = NL
        
            # Take minimium growth rate
            G = min(NL, AveGr)

          

            #================================================================================
            #
            #                    second layer : must be limited by light so only calcualte light growth rate
            #
            #=================================================================================  

            AveGr_sub= 2.*Vmax_sub*SGr_sub/((Depth_EU[Day]-MLZ[Day]) *cn_sgr_sub)  # (2* to get dawn to dusk)
                
            Lig_Gro_sub[Day] = AveGr_sub
            G_sub = AveGr_sub

                                 
            

                                    
            
        #COMPUTE mixing due to mixed layer deepening
            zeta = MLZ[Day] - Mprev if MLZ[Day] > Mprev else 0.
            
        #COMPUTE Grazing rate accounting for grazing threshold
        

            graz = a*k_e*P*P/(a+k_e*P*P)
            #eta=np.sqrt(a/k_e)
            #graz = a*P*P/(eta*eta+P*P)
            
            #eta_sub=np.sqrt(a_sub/k_e_sub)
            #graz_sub=a_sub*P_sub*P_sub/(eta_sub*eta_sub+P_sub*P_sub)
            graz_sub = a_sub*k_e_sub*P_sub*P_sub/(a_sub+k_e_sub*P_sub*P_sub)
            graz_sub_str[Day]=graz_sub
        #COMPUTE Phytoplankton state variable
            #================================================================================
            #
            #                    first layer core function
            #
            #=================================================================================   
        
            DP      =    G * P - m*P - graz*Z 
            DN      =  - G*P + theta*m*P + (1-Gamma-beta)*graz*Z + k_cn*CN*Z + (frac*MLZ[Day]+zeta)*(Nd- Nut)/MLZ[Day]
            DZ      =    graz * Gamma * Z - CN*Z
           
           
            
            #================================================================================
            #
            #                    second layer core function
            #
            #================================================================================= 
            ly2_diff=Depth_EU[Day]-MLZ[Day] 

            DP_sub = G_sub * P_sub - m_sub*P_sub - graz_sub*Z_sub
            DN_sub= - G_sub*P_sub +  m_sub*P_sub + (1-Gamma_sub)*graz_sub*Z_sub +CN_sub*Z_sub 
            DZ_sub = graz_sub * Gamma_sub * Z_sub - CN_sub*Z_sub 
            
            #Nut_sub = Nut_sub*(bry-Mprev)/(bry-MLZ[Day]) + DN_sub
            P_sub= P_sub*(M_EU-Mprev)/(Depth_EU[Day]-MLZ[Day]) + DP_sub 
            Z_sub = Z_sub*(M_EU-Mprev)/(Depth_EU[Day]-MLZ[Day]) + DZ_sub

                
            Nd = Nd*(M_EU-Mprev)/(Depth_EU[Day]-MLZ[Day])- (frac*MLZ[Day]+zeta)*(Nd- Nut)/(Depth_EU[Day]-MLZ[Day])+ (1-k_cn)*CN*Z*MLZ[Day]/(Depth_EU[Day]-MLZ[Day]) +beta*graz*Z*MLZ[Day]/(Depth_EU[Day]-MLZ[Day]) + (1-theta)*m*P*MLZ[Day]/(Depth_EU[Day]-MLZ[Day]) + DN_sub
            P        = P*Mprev/MLZ[Day] + DP
            Nut      = Nut*Mprev/MLZ[Day] +DN
            Z        = Z*Mprev/MLZ[Day]+DZ
           
            Pstr[Day], Nutstr[Day], Zstr[Day]  = P*MLZ[Day] , Nut*MLZ[Day], Z*MLZ[Day]
            Pstr_sub[Day], Zstr_sub[Day], Nd_str[Day] = P_sub*ly2_diff , Z_sub*ly2_diff, Nd*ly2_diff
            Mprev = MLZ[Day]
            M_EU = Depth_EU[Day]
            


    return Pstr, Nutstr, Zstr, Io_noon, Nut_Gro, Lig_Gro, Daystr, I, I_sub, Nd_str, Pstr_sub,  Zstr_sub, Depth_EU,Lig_Gro_sub,graz_sub_str



def NPZ_model_stratified_euphoticlight_oneND_concentration(Lat, year, initial_conditions, parameters, parameters_sub):
    P, Nut, Z, Mprev, MLZ, P_sub, Z_sub, M_EU= initial_conditions

    # Unpack parameters
    SolarK, AtmAtt, ParFrac, Kdw, Kdp, Vmax, alpha, Ks, Nd, m, Gamma, CN, frac, theta, beta,k_cn, a,k_e = parameters
    alpha_sub, Kdp_sub, Vmax_sub, a_sub, k_e_sub, Gamma_sub, m_sub, CN_sub= parameters_sub
    
    #check how many years
    NN=len(year)
    
    #define initial zero value to calculate the total days (time step) of this simulation
    LL=0
    for n in range(NN):
        LL += year_days(year[n])
    
    # Initialize arrays
    Pstr, Nutstr, Zstr, Io_noon, Nut_Gro, Lig_Gro, Daystr, I, Nd_str, cn_term, Depth_EU, I_sub, Lig_Gro_sub,Pstr_sub, Nut_sub_str, Zstr_sub,Vmax_sub_str,graz_sub_str= [np.full(LL, np.nan) for _ in range(18)]



# Loop through each year of year

    Day=-1
    for Yr in range(NN):
        y_day=year_days(year[Yr-1])
    # Loop through the days
        for i in range(y_day):
            Day         += 1
            Daystr[Day] = Day
            #=================================================================================
            #                    COMPUTE THE LIGHT ENVIRONMENT
            #COMPUTE EXTINCTION COEFFICIENT Kd
            #=================================================================================
            Kd         = Kdw + Kdp * P # first layer
            #Depth_EU[Day] = 300  ###13.8
            #Depth_EU[Day] = 9.2/Kd

            Depth_EU[Day] = 13.8/Kd


           
            
            Kd_sub     = Kdw + Kdp_sub*P_sub #second layer
            #print(Kd_sub)

            # Compute surface PAR at latitude for times of day 
            # and dawn to noon and production rate integration down to MLD
            # Decination = angle of sun above equator on each day
            D1          = 23.45*mt.sin(mt.radians((360.*(284.+i+1)/365.)))  
            #=================================================================================
            # Angle(deg) between south (i.e. noon) and setting sun
            W1          = mt.degrees(mt.acos(-1.*(mt.tan(mt.radians(Lat)) * mt.tan(mt.radians(D1)))))
            
            #=================================================================================
            # One-half daylength, hours to noon
            #=================================================================================
            # Earth rotates 15 = degrees / hour =360/24
            L1          = W1/15.
            #=================================================================================
            # Distance of Earth from sun relative to average, a minor effect
            Rx          = 1./mt.sqrt(1.+0.033*np.cos(np.deg2rad((360.*(i+1)/365.))))
            VtofD       = L1/40.
            TofD        = 12.01-L1-VtofD
          
            
            SGr         = 0.
            SGr2        = 0
            SGr_sub     = 0
            cn_sgr2     = 0 
            cn_sgr_sub     = 0 
            
            I_d2n=0 #layer 1 initial value
            I_d2n_sub=0 #laer 2 initial value 
            I_mlz=0 
            # Loop over dawn to noon in 40 steps
            for j in range(40):
                TofD   = TofD + VtofD   # for the following see Brock (1981)
                W2     = (TofD-12.)*15.
                CosZen = np.sin(np.deg2rad(D1))*np.sin(np.deg2rad(Lat)) \
                         +np.cos(np.deg2rad(D1))*np.cos(np.deg2rad(Lat)) \
                         *np.cos(np.deg2rad(W2))
                Isurf  = SolarK*CosZen/(Rx*Rx)
                Io     = Isurf*AtmAtt*ParFrac
                #================================================================================
                #
                #                    first layer 
                #
                #=================================================================================
                sm=0 #first layer initial


                for k in range(1,int(MLZ[Day])+1):
                    
                    Iz = Io*np.exp(-1.*Kd*k) #LORENZEN 1972
                    sm=sm+Iz #accumulated from surface to MLZ
                    #COMPUTE Phytoplankton Growth Rate for each sum from 0 to MLZ and dawn to non             
                                # Denman-Pena function for phyto growth, scales (0-1)
                    Gr     = (1. - mt.exp(-alpha*Iz/Vmax))
                    SGr    += Gr
                I_d2n=I_d2n+sm/int(MLZ[Day]) # calculated MLZ-averaged irradiance from dawn to noon
                #=================================================================================
                #                    second layer: 2cases: bd_layer <mld and when bd_layer >mld
                #
                #=================================================================================

                sm_sub=0
                I_m=Io*np.exp(-1.*Kd*int(MLZ[Day]))  #subsurface light caluclation 
                for k2 in range(int(MLZ[Day])+1, int(Depth_EU[Day])+1):
                    #Iz2 = Io*np.exp(-1.*Kd_sub*k2) #LORENZEN 1972
                    Iz2 = I_m*np.exp(-1.*Kd_sub*k2) #LORENZEN 1972
                    sm_sub=sm_sub+Iz2 #accumulated 
                    Gr_sub     = (1. - mt.exp(-alpha_sub*Iz2/Vmax_sub))
                    SGr_sub    += Gr_sub

                I_d2n_sub=I_d2n_sub+sm_sub/(Depth_EU[Day]-MLZ[Day]) # calculated MLZ-averaged irradiance from dawn to noon
                cn_sgr_sub +=1
            #out of loop 40: 
            
            #out of loop 40: 
            
            #================================================================================
            #
            #                    first layer : averaged light
            #
            #=================================================================================   
        
            I[Day]=2*I_d2n/40 # to get dawn to dusk for surface layer light 

            #=================================================================================
            #                    second layer : averaged light
            #
            #=================================================================================

            I_sub[Day]=2*I_d2n_sub/cn_sgr_sub #layer 2: case2 



       
            #================================================================================
            #
            #                    first layer 
            #
            #=================================================================================

            #COMPUTE Phytoplankton Growth Rate over the day (LIGHT EFFECT)
            AveGr = 2.*Vmax*SGr/(MLZ[Day]*40)  # (2* to get dawn to dusk)
            Lig_Gro[Day] = AveGr

            #COMPUTE Phytoplankton Growth rate NUTRIENT EFFECT)
            NL    = Vmax*Nut/(Ks+Nut) # (NUTRIENT EFFECT)
            Nut_Gro[Day] = NL
        
            # Take minimium growth rate
            G = min(NL, AveGr)

          

            #================================================================================
            #
            #                    second layer : must be limited by light so only calcualte light growth rate
            #
            #=================================================================================  

            AveGr_sub= 2.*Vmax_sub*SGr_sub/((Depth_EU[Day]-MLZ[Day]) *cn_sgr_sub)  # (2* to get dawn to dusk)
                
            Lig_Gro_sub[Day] = AveGr_sub
            G_sub = AveGr_sub

                                 
            

                                    
            
        #COMPUTE mixing due to mixed layer deepening
            zeta = MLZ[Day] - Mprev if MLZ[Day] > Mprev else 0.
            
        #COMPUTE Grazing rate accounting for grazing threshold
        

            graz = a*k_e*P*P/(a+k_e*P*P)
            #eta=np.sqrt(a/k_e)
            #graz = a*P*P/(eta*eta+P*P)
            
            #eta_sub=np.sqrt(a_sub/k_e_sub)
            #graz_sub=a_sub*P_sub*P_sub/(eta_sub*eta_sub+P_sub*P_sub)
            graz_sub = a_sub*k_e_sub*P_sub*P_sub/(a_sub+k_e_sub*P_sub*P_sub)
            graz_sub_str[Day]=graz_sub
        #COMPUTE Phytoplankton state variable
            #================================================================================
            #
            #                    first layer core function
            #
            #=================================================================================   
        
            DP      =    G * P - m*P - graz*Z - (MLZ[Day] - Mprev)*P/MLZ[Day]
            DN      =  - G*P + theta*m*P + (1-Gamma-beta)*graz*Z + k_cn*CN*Z*Z + (frac*MLZ[Day]+zeta)*(Nd- Nut)/MLZ[Day]- (MLZ[Day] - Mprev)*Nut/MLZ[Day]
            DZ      =    graz * Gamma * Z - CN*Z*Z- (MLZ[Day] - Mprev)*Z/MLZ[Day]
           
           
            
            #================================================================================
            #
            #                    second layer core function
            #
            #================================================================================= 


            DP_sub = G_sub * P_sub - m_sub*P_sub - graz_sub*Z_sub -((Depth_EU[Day]-MLZ[Day])-(M_EU-Mprev))*P_sub/(Depth_EU[Day]-MLZ[Day])
            DN_sub= - G_sub*P_sub +  m_sub*P_sub + (1-Gamma_sub)*graz_sub*Z_sub +CN_sub*Z_sub*Z_sub -((Depth_EU[Day]-MLZ[Day])-(M_EU-Mprev))*Nd/(Depth_EU[Day]-MLZ[Day])
            DZ_sub = graz_sub * Gamma_sub * Z_sub - CN_sub*Z_sub*Z_sub-((Depth_EU[Day]-MLZ[Day])-(M_EU-Mprev))*Z_sub/(Depth_EU[Day]-MLZ[Day])
            
            #Nut_sub = Nut_sub*(bry-Mprev)/(bry-MLZ[Day]) + DN_sub
            P_sub= P_sub + DP_sub 
            Z_sub = Z_sub + DZ_sub 
#*(M_EU-Mprev)/(Depth_EU[Day]-MLZ[Day])
                
            Nd = Nd - (frac*MLZ[Day]+zeta)*(Nd- Nut)/(Depth_EU[Day]-MLZ[Day])+ (1-k_cn)*CN*Z*Z*MLZ[Day]/(Depth_EU[Day]-MLZ[Day]) +beta*graz*Z*MLZ[Day]/(Depth_EU[Day]-MLZ[Day]) + (1-theta)*m*P*MLZ[Day]/(Depth_EU[Day]-MLZ[Day]) + DN_sub
            
            P        = P + DP
            Nut      = Nut +DN
            Z        = Z + DZ
           
            Pstr[Day], Nutstr[Day], Zstr[Day]  = P , Nut, Z
            Pstr_sub[Day], Zstr_sub[Day], Nd_str[Day] = P_sub , Z_sub, Nd
            Mprev = MLZ[Day]
            M_EU = Depth_EU[Day]
            


    return Pstr, Nutstr, Zstr, Io_noon, Nut_Gro, Lig_Gro, Daystr, I, I_sub, Nd_str, Pstr_sub,  Zstr_sub, Depth_EU,Lig_Gro_sub,graz_sub_str





def NPZ_model_stratified_euphoticlight_oneND_qua(Lat, year, initial_conditions, parameters, parameters_sub):
    #this function is the final verison of permenant stratifed ocean model with two-layer (with uphotic zone always deeper than MLD)
    P, Nut, Z, Mprev, MLZ, P_sub, Z_sub, M_EU= initial_conditions

    # Unpack parameters
    SolarK, AtmAtt, ParFrac, Kdw, Kdp, Vmax, alpha, Ks, Nd, m, Gamma, CN, frac, theta, beta,k_cn, a,k_e = parameters
    alpha_sub, Kdp_sub, Vmax_sub, a_sub, k_e_sub, Gamma_sub, m_sub, CN_sub= parameters_sub
    
    #check how many years
    NN=len(year)
    
    #define initial zero value to calculate the total days (time step) of this simulation
    LL=0
    for n in range(NN):
        LL += year_days(year[n])
    
    # Initialize arrays
    Pstr, Nutstr, Zstr, Io_noon, Nut_Gro, Lig_Gro, Daystr, I, Nd_str, cn_term, Depth_EU, I_sub, Lig_Gro_sub,Pstr_sub, Nut_sub_str, Zstr_sub,Vmax_sub_str,graz_sub_str= [np.full(LL, np.nan) for _ in range(18)]



# Loop through each year of year

    Day=-1
    for Yr in range(NN):
        y_day=year_days(year[Yr-1])
    # Loop through the days
        for i in range(y_day):
            Day         += 1
            Daystr[Day] = Day
            #=================================================================================
            #                    COMPUTE THE LIGHT ENVIRONMENT
            #COMPUTE EXTINCTION COEFFICIENT Kd
            #=================================================================================
            Kd         = Kdw + Kdp * P # first layer
            #Depth_EU[Day] = 300  ###13.8
            #Depth_EU[Day] = 9.2/Kd

            Depth_EU[Day] = 13.8/Kd


           
            
            Kd_sub     = Kdw + Kdp_sub*P_sub #second layer
            #print(Kd_sub)

            # Compute surface PAR at latitude for times of day 
            # and dawn to noon and production rate integration down to MLD
            # Decination = angle of sun above equator on each day
            D1          = 23.45*mt.sin(mt.radians((360.*(284.+i+1)/365.)))  
            #=================================================================================
            # Angle(deg) between south (i.e. noon) and setting sun
            W1          = mt.degrees(mt.acos(-1.*(mt.tan(mt.radians(Lat)) * mt.tan(mt.radians(D1)))))
            
            #=================================================================================
            # One-half daylength, hours to noon
            #=================================================================================
            # Earth rotates 15 = degrees / hour =360/24
            L1          = W1/15.
            #=================================================================================
            # Distance of Earth from sun relative to average, a minor effect
            Rx          = 1./mt.sqrt(1.+0.033*np.cos(np.deg2rad((360.*(i+1)/365.))))
            VtofD       = L1/40.
            TofD        = 12.01-L1-VtofD
          
            
            SGr         = 0.
            SGr2        = 0
            SGr_sub     = 0
            cn_sgr2     = 0 
            cn_sgr_sub     = 0 
            
            I_d2n=0 #layer 1 initial value
            I_d2n_sub=0 #laer 2 initial value 
            I_mlz=0 
            # Loop over dawn to noon in 40 steps
            for j in range(40):
                TofD   = TofD + VtofD   # for the following see Brock (1981)
                W2     = (TofD-12.)*15.
                CosZen = np.sin(np.deg2rad(D1))*np.sin(np.deg2rad(Lat)) \
                         +np.cos(np.deg2rad(D1))*np.cos(np.deg2rad(Lat)) \
                         *np.cos(np.deg2rad(W2))
                Isurf  = SolarK*CosZen/(Rx*Rx)
                Io     = Isurf*AtmAtt*ParFrac
                #================================================================================
                #
                #                    first layer 
                #
                #=================================================================================
                sm=0 #first layer initial


                for k in range(1,int(MLZ[Day])+1):
                    
                    Iz = Io*np.exp(-1.*Kd*k) #LORENZEN 1972
                    sm=sm+Iz #accumulated from surface to MLZ
                    #COMPUTE Phytoplankton Growth Rate for each sum from 0 to MLZ and dawn to non             
                                # Denman-Pena function for phyto growth, scales (0-1)
                    Gr     = (1. - mt.exp(-alpha*Iz/Vmax))
                    SGr    += Gr
                I_d2n=I_d2n+sm/int(MLZ[Day]) # calculated MLZ-averaged irradiance from dawn to noon
                #=================================================================================
                #                    second layer: 2cases: bd_layer <mld and when bd_layer >mld
                #
                #=================================================================================

                sm_sub=0
                I_m=Io*np.exp(-1.*Kd*int(MLZ[Day]))  #subsurface light caluclation 
                for k2 in range(int(MLZ[Day])+1, int(Depth_EU[Day])+1):
                    #Iz2 = Io*np.exp(-1.*Kd_sub*k2) #LORENZEN 1972
                    Iz2 = I_m*np.exp(-1.*Kd_sub*k2) #LORENZEN 1972
                    sm_sub=sm_sub+Iz2 #accumulated 
                    Gr_sub     = (1. - mt.exp(-alpha_sub*Iz2/Vmax_sub))
                    SGr_sub    += Gr_sub

                I_d2n_sub=I_d2n_sub+sm_sub/(Depth_EU[Day]-MLZ[Day]) # calculated MLZ-averaged irradiance from dawn to noon
                cn_sgr_sub +=1
            #out of loop 40: 
            
            #out of loop 40: 
            
            #================================================================================
            #
            #                    first layer : averaged light
            #
            #=================================================================================   
        
            I[Day]=2*I_d2n/40 # to get dawn to dusk for surface layer light 

            #=================================================================================
            #                    second layer : averaged light
            #
            #=================================================================================

            I_sub[Day]=2*I_d2n_sub/cn_sgr_sub #layer 2: case2 



       
            #================================================================================
            #
            #                    first layer 
            #
            #=================================================================================

            #COMPUTE Phytoplankton Growth Rate over the day (LIGHT EFFECT)
            AveGr = 2.*Vmax*SGr/(MLZ[Day]*40)  # (2* to get dawn to dusk)
            Lig_Gro[Day] = AveGr

            #COMPUTE Phytoplankton Growth rate NUTRIENT EFFECT)
            NL    = Vmax*Nut/(Ks+Nut) # (NUTRIENT EFFECT)
            Nut_Gro[Day] = NL
        
            # Take minimium growth rate
            G = min(NL, AveGr)

          

            #================================================================================
            #
            #                    second layer : must be limited by light so only calcualte light growth rate
            #
            #=================================================================================  

            AveGr_sub= 2.*Vmax_sub*SGr_sub/((Depth_EU[Day]-MLZ[Day]) *cn_sgr_sub)  # (2* to get dawn to dusk)
                
            Lig_Gro_sub[Day] = AveGr_sub
            G_sub = AveGr_sub

                                 
            

                                    
            
        #COMPUTE mixing due to mixed layer deepening
            zeta = MLZ[Day] - Mprev if MLZ[Day] > Mprev else 0.
            
        #COMPUTE Grazing rate accounting for grazing threshold
        

            graz = a*k_e*P*P/(a+k_e*P*P)
            #eta=np.sqrt(a/k_e)
            #graz = a*P*P/(eta*eta+P*P)
            
            #eta_sub=np.sqrt(a_sub/k_e_sub)
            #graz_sub=a_sub*P_sub*P_sub/(eta_sub*eta_sub+P_sub*P_sub)
            graz_sub = a_sub*k_e_sub*P_sub*P_sub/(a_sub+k_e_sub*P_sub*P_sub)
            graz_sub_str[Day]=graz_sub
        #COMPUTE Phytoplankton state variable
            #================================================================================
            #
            #                    first layer core function
            #
            #=================================================================================   
        
            DP      =    G * P - m*P - graz*Z 
            DN      =  - G*P + theta*m*P + (1-Gamma-beta)*graz*Z + k_cn*CN*Z*Z + (frac*MLZ[Day]+zeta)*(Nd- Nut)/MLZ[Day]
            DZ      =    graz * Gamma * Z - CN*Z*Z
           
           
            
            #================================================================================
            #
            #                    second layer core function
            #
            #================================================================================= 
            ly2_diff=Depth_EU[Day]-MLZ[Day] 

            DP_sub = G_sub * P_sub - m_sub*P_sub - graz_sub*Z_sub
            DN_sub= - G_sub*P_sub +  m_sub*P_sub + (1-Gamma_sub)*graz_sub*Z_sub +CN_sub*Z_sub*Z_sub
            DZ_sub = graz_sub * Gamma_sub * Z_sub - CN_sub*Z_sub*Z_sub
            
            #Nut_sub = Nut_sub*(bry-Mprev)/(bry-MLZ[Day]) + DN_sub
            P_sub= P_sub*(M_EU-Mprev)/(Depth_EU[Day]-MLZ[Day]) + DP_sub 
            Z_sub = Z_sub*(M_EU-Mprev)/(Depth_EU[Day]-MLZ[Day]) + DZ_sub

                
            Nd = Nd*(M_EU-Mprev)/(Depth_EU[Day]-MLZ[Day])- (frac*MLZ[Day]+zeta)*(Nd- Nut)/(Depth_EU[Day]-MLZ[Day])+ (1-k_cn)*CN*Z*Z*MLZ[Day]/(Depth_EU[Day]-MLZ[Day]) +beta*graz*Z*MLZ[Day]/(Depth_EU[Day]-MLZ[Day]) + (1-theta)*m*P*MLZ[Day]/(Depth_EU[Day]-MLZ[Day]) + DN_sub
            P        = P*Mprev/MLZ[Day] + DP
            Nut      = Nut*Mprev/MLZ[Day] +DN
            Z        = Z*Mprev/MLZ[Day]+DZ
           
            Pstr[Day], Nutstr[Day], Zstr[Day]  = P*MLZ[Day] , Nut*MLZ[Day], Z*MLZ[Day]
            Pstr_sub[Day], Zstr_sub[Day], Nd_str[Day] = P_sub*ly2_diff , Z_sub*ly2_diff, Nd*ly2_diff
            Mprev = MLZ[Day]
            M_EU = Depth_EU[Day]
            


    return Pstr, Nutstr, Zstr, Io_noon, Nut_Gro, Lig_Gro, Daystr, I, I_sub, Nd_str, Pstr_sub,  Zstr_sub, Depth_EU,Lig_Gro_sub,graz_sub_str


def NPZ_model_stratified_euphoticlight_oneND_qua_paper(Lat, year, initial_conditions, parameters, parameters_sub):
    #this function is the final verison of permenant stratifed ocean model with two-layer (with uphotic zone always deeper than MLD)
    P, Nut, Z, Mprev, MLZ, P_sub, Z_sub, M_EU= initial_conditions

    # Unpack parameters
    SolarK, AtmAtt, ParFrac, Kdw, Kdp, Vmax, alpha, Ks, Nd, m, Gamma, CN, frac, theta, beta,k_cn, a,k_e = parameters
    alpha_sub, Kdp_sub, Vmax_sub, a_sub, k_e_sub, Gamma_sub, m_sub, CN_sub= parameters_sub
    
    #check how many years
    NN=len(year)
    
    #define initial zero value to calculate the total days (time step) of this simulation
    LL=0
    for n in range(NN):
        LL += year_days(year[n])
    
    # Initialize arrays
    Pstr, Nutstr, Zstr, Io_noon, Nut_Gro, Lig_Gro, Daystr, I, Nd_str, cn_term, Depth_EU, I_sub, Lig_Gro_sub,Pstr_sub, Nut_sub_str, Zstr_sub,Vmax_sub_str,graz_sub_str= [np.full(LL, np.nan) for _ in range(18)]



# Loop through each year of year

    Day=-1
    for Yr in range(NN):
        y_day=year_days(year[Yr-1])
    # Loop through the days
        for i in range(y_day):
            Day         += 1
            Daystr[Day] = Day
            #=================================================================================
            #                    COMPUTE THE LIGHT ENVIRONMENT
            #COMPUTE EXTINCTION COEFFICIENT Kd
            #=================================================================================
            Kd         = Kdw + Kdp * P # first layer
            #Depth_EU[Day] = 300  ###13.8
            #Depth_EU[Day] = 9.2/Kd

            Depth_EU[Day] = 13.8/Kd


           
            
            Kd_sub     = Kdw + Kdp_sub*P_sub #second layer
            #print(Kd_sub)

            # Compute surface PAR at latitude for times of day 
            # and dawn to noon and production rate integration down to MLD
            # Decination = angle of sun above equator on each day
            D1          = 23.45*mt.sin(mt.radians((360.*(284.+i+1)/365.)))  
            #=================================================================================
            # Angle(deg) between south (i.e. noon) and setting sun
            W1          = mt.degrees(mt.acos(-1.*(mt.tan(mt.radians(Lat)) * mt.tan(mt.radians(D1)))))
            
            #=================================================================================
            # One-half daylength, hours to noon
            #=================================================================================
            # Earth rotates 15 = degrees / hour =360/24
            L1          = W1/15.
            #=================================================================================
            # Distance of Earth from sun relative to average, a minor effect
            Rx          = 1./mt.sqrt(1.+0.033*np.cos(np.deg2rad((360.*(i+1)/365.))))
            VtofD       = L1/40.
            TofD        = 12.01-L1-VtofD
          
            
            SGr         = 0.
            SGr2        = 0
            SGr_sub     = 0
            cn_sgr2     = 0 
            cn_sgr_sub     = 0 
            
            I_d2n=0 #layer 1 initial value
            I_d2n_sub=0 #laer 2 initial value 
            I_mlz=0 
            # Loop over dawn to noon in 40 steps
            for j in range(40):
                TofD   = TofD + VtofD   # for the following see Brock (1981)
                W2     = (TofD-12.)*15.
                CosZen = np.sin(np.deg2rad(D1))*np.sin(np.deg2rad(Lat)) \
                         +np.cos(np.deg2rad(D1))*np.cos(np.deg2rad(Lat)) \
                         *np.cos(np.deg2rad(W2))
                Isurf  = SolarK*CosZen/(Rx*Rx)
                Io     = Isurf*AtmAtt*ParFrac
                #================================================================================
                #
                #                    first layer 
                #
                #=================================================================================
                sm=0 #first layer initial


                for k in range(1,int(MLZ[Day])+1):
                    
                    Iz = Io*np.exp(-1.*Kd*k) #LORENZEN 1972
                    sm=sm+Iz #accumulated from surface to MLZ
                    #COMPUTE Phytoplankton Growth Rate for each sum from 0 to MLZ and dawn to non             
                                # Denman-Pena function for phyto growth, scales (0-1)
                    Gr     = (1. - mt.exp(-alpha*Iz/Vmax))
                    SGr    += Gr
                I_d2n=I_d2n+sm/int(MLZ[Day]) # calculated MLZ-averaged irradiance from dawn to noon
                #=================================================================================
                #                    second layer: 2cases: bd_layer <mld and when bd_layer >mld
                #
                #=================================================================================

                sm_sub=0
                I_m=Io*np.exp(-1.*Kd*int(MLZ[Day]))  #subsurface light caluclation 
                for k2 in range(int(MLZ[Day])+1, int(Depth_EU[Day])+1):
                    #Iz2 = Io*np.exp(-1.*Kd_sub*k2) #LORENZEN 1972
                    Iz2 = I_m*np.exp(-1.*Kd_sub*k2) #LORENZEN 1972
                    sm_sub=sm_sub+Iz2 #accumulated 
                    Gr_sub     = (1. - mt.exp(-alpha_sub*Iz2/Vmax_sub))
                    SGr_sub    += Gr_sub

                I_d2n_sub=I_d2n_sub+sm_sub/(Depth_EU[Day]-MLZ[Day]) # calculated MLZ-averaged irradiance from dawn to noon
                cn_sgr_sub +=1
            #out of loop 40: 
            
            #out of loop 40: 
            
            #================================================================================
            #
            #                    first layer : averaged light
            #
            #=================================================================================   
        
            I[Day]=2*I_d2n/40 # to get dawn to dusk for surface layer light 

            #=================================================================================
            #                    second layer : averaged light
            #
            #=================================================================================

            I_sub[Day]=2*I_d2n_sub/cn_sgr_sub #layer 2: case2 



       
            #================================================================================
            #
            #                    first layer 
            #
            #=================================================================================

            #COMPUTE Phytoplankton Growth Rate over the day (LIGHT EFFECT)
            AveGr = 2.*Vmax*SGr/(MLZ[Day]*40)  # (2* to get dawn to dusk)
            Lig_Gro[Day] = AveGr

            #COMPUTE Phytoplankton Growth rate NUTRIENT EFFECT)
            NL    = Vmax*Nut/(Ks+Nut) # (NUTRIENT EFFECT)
            Nut_Gro[Day] = NL
        
            # Take minimium growth rate
            G = min(NL, AveGr)

          

            #================================================================================
            #
            #                    second layer : must be limited by light so only calcualte light growth rate
            #
            #=================================================================================  

            AveGr_sub= 2.*Vmax_sub*SGr_sub/((Depth_EU[Day]-MLZ[Day]) *cn_sgr_sub)  # (2* to get dawn to dusk)
                
            Lig_Gro_sub[Day] = AveGr_sub
            G_sub = AveGr_sub

                                 
            

                                    
            
        #COMPUTE mixing due to mixed layer deepening
            zeta = MLZ[Day] - Mprev if MLZ[Day] > Mprev else 0.
            
        #COMPUTE Grazing rate accounting for grazing threshold
        

            graz = a*k_e*P*P/(a+k_e*P*P)
            #eta=np.sqrt(a/k_e)
            #graz = a*P*P/(eta*eta+P*P)
            
            #eta_sub=np.sqrt(a_sub/k_e_sub)
            #graz_sub=a_sub*P_sub*P_sub/(eta_sub*eta_sub+P_sub*P_sub)
            graz_sub = a_sub*k_e_sub*P_sub*P_sub/(a_sub+k_e_sub*P_sub*P_sub)
            graz_sub_str[Day]=graz_sub
        #COMPUTE Phytoplankton state variable
            #================================================================================
            #
            #                    first layer core function
            #
            #=================================================================================   
        
            DP      =    G * P - m*P - graz*Z - (MLZ[Day] - Mprev)*P/MLZ[Day]
            DN      =  - G*P + theta*m*P + (1-Gamma-beta)*graz*Z + k_cn*CN*Z*Z + (frac*MLZ[Day]+zeta)*(Nd- Nut)/MLZ[Day]- (MLZ[Day] - Mprev)*Nut/MLZ[Day]
            DZ      =    graz * Gamma * Z - CN*Z*Z- (MLZ[Day] - Mprev)*Z/MLZ[Day]
           
           
            
            #================================================================================
            #
            #                    second layer core function
            #
            #================================================================================= 
            ly2_diff=Depth_EU[Day]-MLZ[Day] 

            DP_sub = G_sub * P_sub - m_sub*P_sub - graz_sub*Z_sub -((Depth_EU[Day]-MLZ[Day])-(M_EU-Mprev))*P_sub/(Depth_EU[Day]-MLZ[Day])
            DN_sub= - G_sub*P_sub +  m_sub*P_sub + (1-Gamma_sub)*graz_sub*Z_sub +CN_sub*Z_sub*Z_sub -((Depth_EU[Day]-MLZ[Day])-(M_EU-Mprev))*Nd/(Depth_EU[Day]-MLZ[Day])
            DZ_sub = graz_sub * Gamma_sub * Z_sub - CN_sub*Z_sub*Z_sub-((Depth_EU[Day]-MLZ[Day])-(M_EU-Mprev))*Z_sub/(Depth_EU[Day]-MLZ[Day])
            
            #Nut_sub = Nut_sub*(bry-Mprev)/(bry-MLZ[Day]) + DN_sub
            P_sub= P_sub + DP_sub 
            Z_sub = Z_sub + DZ_sub 
#*(M_EU-Mprev)/(Depth_EU[Day]-MLZ[Day])
                
            Nd = Nd - (frac*MLZ[Day]+zeta)*(Nd- Nut)/(Depth_EU[Day]-MLZ[Day])+ (1-k_cn)*CN*Z*Z*MLZ[Day]/(Depth_EU[Day]-MLZ[Day]) +beta*graz*Z*MLZ[Day]/(Depth_EU[Day]-MLZ[Day]) + (1-theta)*m*P*MLZ[Day]/(Depth_EU[Day]-MLZ[Day]) + DN_sub
            
            P        = P + DP
            Nut      = Nut +DN
            Z        = Z + DZ
           
            Pstr[Day], Nutstr[Day], Zstr[Day]  = P*MLZ[Day] , Nut*MLZ[Day], Z*MLZ[Day]
            Pstr_sub[Day], Zstr_sub[Day], Nd_str[Day] = P_sub*ly2_diff , Z_sub*ly2_diff, Nd*ly2_diff
            Mprev = MLZ[Day]
            M_EU = Depth_EU[Day]
            


    return Pstr, Nutstr, Zstr, Io_noon, Nut_Gro, Lig_Gro, Daystr, I, I_sub, Nd_str, Pstr_sub,  Zstr_sub, Depth_EU,Lig_Gro_sub,graz_sub_str

