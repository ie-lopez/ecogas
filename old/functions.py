def sfrdelayed(tm,t,tb,ab,f_burst,sfrti,ti):
    #Scaling parameter
    scale=sfrti/(ti*np.exp(-ti/tm)/tm**2)
    # SFR for each component
    sfr = t * np.exp(-t / tm) / tm**2
    z0=cosmo.age(0).to(u.Myr).value
    if t>z0-ab:
        sfr_burst = np.exp(-t / tb)
        sfr_burst *= (f_burst / (1.-f_burst) * np.sum(sfr) / np.sum(sfr_burst))        
        sfr += sfr_burst
    return sfr*scale#*u.Msun/u.yr

def mst_to_radius(sf,mass):
    #From Ichikawa+12
    a=0.103
    logr=0.874
    if sf==0:
        #star forming
        m=np.log10(mass)-10
    elif sf==1:
        #quiescent
        m=np.log10(mass)-11
    radius=(10**(a*m+logr))*u.kpc
    return radius

def is_sf(ssfr,z):
    epoch_norm_ssfr=ssfr-0.3*z+9.62
    if epoch_norm_ssfr<-1.2:
        #is quiescent
        tipo=1
    elif epoch_norm_ssfr>-1.2:
        #is star forming
        tipo=0
    return tipo

def Mste_to_Mhalo(z,Mste):
    #From Girelli+20 - Table2
    if z<0.2:
        A,Ma,B,G=0.0494,11.81,0.94,0.726
    elif z>=0.2 and z<0.5:
        A,Ma,B,G=0.0429,11.87,0.99,0.669   
    elif z>=0.5 and z<0.8:
        A,Ma,B,G=0.0348,12.07,0.86,0.622
    elif z>=0.8 and z<1.1:
        A,Ma,B,G=0.0429,12.03,1.04,0.657
    elif z>=1.1 and z<1.5:
        A,Ma,B,G=0.0325,12.11,0.87,0.659
    elif z>=1.5 and z<2:
        A,Ma,B,G=0.0285,12.21,0.94,0.624
    elif z>=2 and z<2.5:
        A,Ma,B,G=0.0297,12.23,1.31,0.604
    elif z>=2.5 and z<3:
        A,Ma,B,G=0.0294,12.33,1.13,0.583
    elif z>=3 and z<=4.1:
        A,Ma,B,G=0.0330,12.55,1.05,0.626
    elif z<0:
        print('Z lower than 0')
        return
    else:
        print('Invalid Z')
        return
    Ma=10**Ma
    range_mhalos=10**np.arange(10,15,0.1)
    Mhalo=range_mhalos
    range_mstars=2*A*Mhalo/((Mhalo/Ma)**(-B)+(Mhalo/Ma)**G)
    index=np.argmin(abs(range_mstars-Mste.value))
    Mhalo=range_mhalos[index]*u.Msun
    return Mhalo

def mhalo_at_r(Mhalo_total,r,z):
    #Use Jaffe profile to obtain Mhalo at given radius
    ##To obtain profile's parameters 
    #Critic density of the Universe at z
    p200=(3*cosmo.H(z)**2/(8*np.pi*const.G)).to(u.kg/u.m**3)
    #r200=r_vir
    r200=((3*Mhalo_total/(4*np.pi*p200))**(1/3)).to(u.kpc)
    #Concentration from c-Mh relation Duffy+08
    Mpivot=2e12*u.Msun
    alpha,beta,gamma=5.71,-0.084,-0.47
    c200=alpha*(Mhalo_total/Mpivot)**beta*(1+z)**gamma
    ##Profile
    rs=r200/c200
    ps=Mhalo_total/(4*np.pi*rs**3*(r200/(rs+r200)))
    Mhr=4*np.pi*ps*rs**3*(r/(r+rs))
    return Mhr

# def mhalo_at_r(Mhalo_total,r,z):
#     #Use Isothermical profile to obtain Mhalo at given radius
#     ##To obtain profile's parameters 
#     #Critic density of the Universe at z
#     p200=(3*cosmo.H(z)**2/(8*np.pi*const.G)).to(u.kg/u.m**3)
#     #r200=r_vir
#     r200=((3*Mhalo_total/(4*np.pi*p200))**(1/3)).to(u.kpc)
#     #Concentration from c-Mh relation Duffy+08
#     Mpivot=2e12*u.Msun
#     alpha,beta,gamma=5.71,-0.084,-0.47
#     c200=alpha*(Mhalo_total/Mpivot)**beta*(1+z)**gamma
#     ##Profile
#     rs=r200/c200
#     ps=Mhalo_total/(4*np.pi*rs**3*(np.log(1+c200)))
#     Mhr=4*np.pi*ps*rs**3*(np.log(1+r/rs))
#     return Mhr

def mass_to_z0(cosmo,z,mbh,mbh_dot,mstellar,mstellar_dot,mgas,mode,tm,tb,ab,f_burst):
    #Correct the stellar mass and black hole mass from observed Z to Z=0
    #3 different modes: 1 for constant rates, 2 for variables rates with limited Gas Mass, 3 for variable rates with limited and Ebh<Ecool
    #Variables rates:
        #For stellar mass, recover the SFR(t) using the SFH delayed model from Cigale
        #For black hole mass use the BHAR(t) using the minimun between Bondi accretion rate and Eddington accretion rate
    #Rest of variables:
        #mbh,mbh_p,mstellar,mstellar_p: Black Hole Mass, rate, Stellar Mass, SFR - at observed z
        #z: Observed redshift
        #am,tb,ab,f_burst: Parameters from the SFH model of Cigale.
    
    #Create array with bin in time
    age_obs,age_z0=cosmo.age(z).to(u.Myr).value,cosmo.age(0).to(u.Myr).value
    interval=100 #Bin=10 Myr
    time=np.arange(age_obs,age_z0,interval)
    mst_t=[]*u.Msun
    mbh_t=[]*u.Msun
    SFH=[]*u.Msun/u.yr
    BHH=[]*u.Msun/u.yr
    mst=mstellar.copy()
    mbh=mbh.copy()
    e_bh=0*u.J
    epg_last=0*u.J
    rate_ST=mstellar_dot.to(u.Msun/u.Myr)
    rate_BH=mbh_dot.to(u.Msun/u.Myr)
    if mode==1:
        for binning in time:
            mst_t=np.append(mst_t,mst+(rate_ST*interval*u.Myr).to(u.Msun))
            mbh_t=np.append(mbh_t,mbh+(rate_BH*interval*u.Myr).to(u.Msun))
            SFH=np.append(SFH,rate_ST.to(u.Msun/u.yr))
            BHH=np.append(BHH,rate_BH.to(u.Msun/u.yr))
            mst=mst_t[-1] #Last element of the array
            mbh=mbh_t[-1]
    elif mode==2:
        for binning in time:
            rate_ST=sfrdelayed(tm,binning,tb,ab,f_burst,mstellar_dot,age_obs).to(u.Msun/u.Myr)
            rate_BH=sfrdelayed(tm,binning,tb,ab,f_burst,mbh_dot,age_obs).to(u.Msun/u.Myr)
            mst_t=np.append(mst_t,mst+(rate_ST*interval*u.Myr).to(u.Msun))
            mbh_t=np.append(mbh_t,mbh+(rate_BH*interval*u.Myr).to(u.Msun))
            SFH=np.append(SFH,rate_ST.to(u.Msun/u.yr))
            BHH=np.append(BHH,rate_BH.to(u.Msun/u.yr))
            mst=mst_t[-1] #Last element of the array
            mbh=mbh_t[-1]
    elif mode==3:
        for binning in time:
            test=False
            if test==True:
                rate_ST=0*u.Msun/u.Myr
                rate_BH=0*u.Msun/u.Myr
            else:
                sf=is_sf(np.log10(rate_ST.to(u.Msun/u.yr).value/mst.value),z)
                r_gax=mst_to_radius(sf,mst.value)
                z_at_binning=z_at_value(cosmo.age,binning*u.Myr)
                Mhalo_total=Mste_to_Mhalo(z_at_binning,mst)
                mhalo=mhalo_at_r(Mhalo_total,r_gax,z_at_binning)
                epg=(const.G*(mhalo+mst+mgas+mbh)**2/r_gax).to(u.J)
                #print(binning.round(0),np.log10(mbh.value).round(1),np.log10(mst.value).round(1),np.log10(mhalo.value).round(1),np.log10(mgas.value).round(1),np.log10(epg.value).round(1),np.log10(e_bh.value).round(1))
                if e_bh.value<=epg.value:
                    rate_BH=sfrdelayed(tm,binning,tb,ab,f_burst,mbh_dot,age_obs).to(u.Msun/u.Myr)
                else:
                    #diff=abs(epg-epg_last) #epg - previous epg
                    #rate_BH=(diff/(interval*u.Myr*0.1*const.c**2)).to(u.Msun/u.Myr)     
                    rate_BH=0*u.Msun/u.Myr
                rate_ST=sfrdelayed(tm,binning,tb,ab,f_burst,mstellar_dot,age_obs).to(u.Msun/u.Myr)
            mst_t=np.append(mst_t,mst+(rate_ST*interval*u.Myr).to(u.Msun))
            mbh_t=np.append(mbh_t,mbh+(rate_BH*interval*u.Myr).to(u.Msun))
            e_bh+=(0.1*(rate_BH*interval*u.Myr).to(u.Msun)*const.c**2).to(u.J)
            epg_last=epg.copy()
            mgas=mgas-(rate_ST*interval*u.Myr).to(u.Msun)-(rate_BH*interval*u.Myr).to(u.Msun)
            if mgas<0:
                mgas=0*u.Msun
            SFH=np.append(SFH,rate_ST.to(u.Msun/u.yr))
            BHH=np.append(BHH,rate_BH.to(u.Msun/u.yr))
            mst=mst_t[-1] #Last element of the array
            mbh=mbh_t[-1]
    return time, mst_t, mbh_t, SFH, BHH