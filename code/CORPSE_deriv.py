expected_params={	'vmaxref': 'Relative maximum enzymatic decomp rates (length 3)',
        	'Ea':	'Activation energy (length 3)',
        	'kC':	'Michaelis-Menton parameter (length 3)',
        	'gas_diffusion_exp': 'Determines suppression of decomp at high soil moisture',
            'substrate_diffusion_exp':'Determines suppression of decomp at low soil moisture',
        	'minMicrobeC':	   'Minimum microbial biomass (fraction of total C)',
        	'Tmic': 'Microbial lifetime at 20C (years)',
        	'et':		  'Fraction of turnover not converted to CO2',
        	'eup': 'Carbon uptake efficiency (length 3)',
            'nup': 'Nitrogen uptake efficiency (length 3)',
        	'tProtected':	'Protected C turnover time (years)',
        	'protection_rate':'Protected carbon formation rate (year-1) (length 3)',
            'CN_microbe': 'C:N ratio of microbial biomass',
            'frac_N_turnover_min': 'Fraction of microbial biomass N turnover that is mineralized',
            'frac_turnover_slow': 'Fraction of microbial biomass N turnover that goes to slow pool',
            'max_immobilization_rate': 'Maximum N immobilization rate (fraction per day)',
            'new_resp_units':'If true, vmaxref has units of 1/years and assumes optimal soil moisture has a relative rate of 1.0'
            }

chem_types = ['Fast','Slow','Necro']

# This sets up pools like "uFastC" for unprotected fast C and "pNecroN" for protected necromass N
expected_pools = ['u'+t+'C' for t in chem_types]+\
                 ['p'+t+'C' for t in chem_types]+\
                 ['u'+t+'N' for t in chem_types]+\
                 ['p'+t+'N' for t in chem_types]+\
                 ['livingMicrobeC','livingMicrobeN','CO2','inorganicN']


#All soils: slope=0.4833,intercept=2.3282
#Alfisols: slope=0.5945, intercept=2.2788
def prot_clay(claypercent,slope=0.4833,intercept=2.3282,BD=1.15,porosity=0.4):
    ''' Calculate protection rate as a function of clay content, based on sorption isotherms from Mayes et al (2012) Table 3
    Calculates Qmax in mgC/kg soil from Mayes et al 2012, converted to g/m3 using bulk density
    Typically used as relative value for calculating protection_rate parameter.
    claypercent: Soil % clay (out of 100)
    slope: Either for all soils or a soil order, from Mayes paper
    intercept: Either for all soils or a soil order, from Mayes paper
    BD: Soil bulk density in g/cm3
    '''
    from numpy import log10,where,atleast_1d
    prot=where(atleast_1d(claypercent)!=0.0,1.0*(10**(slope*log10(claypercent)+intercept)*BD*1e-6),0.0)
    return prot


def check_params(params):
    '''params: dictionary containing parameter values. Must contain all the fields in expected_params'''

    from numpy import iterable,array
    unused_params=expected_params.copy()
    for k in params.keys():
        if k not in expected_params:
            raise ValueError('Parameter set contains unexpected parameter %s'%k)
        unused_params.pop(k)
        if iterable(params[k]):
            params[k]=array(params[k])
    if len(unused_params)>0:
        for k in unused_params.keys():
            print ('Missing parameter: %s [%s]'%(k,unused_params[k]))
        raise ValueError('Missing parameters: %s'%unused_params.keys())


from numpy import zeros,size,where,atleast_1d,zeros_like
def CORPSE_deriv(SOM,T,theta,params,claymod=1.0):
    '''Calculate rates of change for all CORPSE pools
       T: Temperature (K)
       theta: Soil water content (fraction of saturation)

       Returns same data structure as SOM'''

    # if any(T<0.0):
    #     raise ValueError('T must be >=0')

    theta=atleast_1d(theta)
    T=atleast_1d(T)

    theta[theta<0]=0.0
    theta[theta>1]=1.0

    et=params['et']
    eup=params['eup']

    # Calculate maximum potential C decomposition rate
    decomp=decompRate(SOM,T,theta,params)

    # Microbial turnover
    microbeTurnover=atleast_1d((SOM['livingMicrobeC']-params['minMicrobeC']*(sumCtypes(SOM,'u')))/params['Tmic']);   # kg/m2/yr
    if isinstance(microbeTurnover,float):
        microbeTurnover=max(0.0,microbeTurnover)
    else:
        microbeTurnover[microbeTurnover<0.0]=0.0

    maintenance_resp=microbeTurnover*(1.0-et)
    overflow_resp=zeros_like(maintenance_resp)

    deadmic_C_production=microbeTurnover*et   # actual fraction of microbial turnover
    deadmic_N_production=microbeTurnover*et/params['CN_microbe']

    # C and N available for microbial growth
    carbon_supply=zeros_like(microbeTurnover)
    nitrogen_supply=zeros_like(microbeTurnover)
    for t in chem_types:
        carbon_supply=carbon_supply+decomp[t+'C']*params['eup'][t]
        nitrogen_supply=nitrogen_supply+decomp[t+'N']*params['nup'][t]

    IMM_N_max=atleast_1d(params['max_immobilization_rate']*365*SOM['inorganicN']/(SOM['inorganicN']+params['max_immobilization_rate']))

    dmicrobeC=zeros_like(microbeTurnover)
    dmicrobeN=zeros_like(microbeTurnover)
    CN_imbalance_term=zeros_like(IMM_N_max)


    # Growth is nitrogen limited, with not enough mineral N to support it with max immobilization
    # loc_Nlim is a vector of True/False that tells the code where this condition applies
    loc_Nlim=(carbon_supply - maintenance_resp)>((nitrogen_supply+IMM_N_max)*params['CN_microbe'])
    CN_imbalance_term[loc_Nlim] = -IMM_N_max[loc_Nlim]
    dmicrobeC[loc_Nlim] =  ((nitrogen_supply[loc_Nlim]+IMM_N_max[loc_Nlim])*params['CN_microbe'] - microbeTurnover[loc_Nlim]*et)
    dmicrobeN[loc_Nlim] = (nitrogen_supply[loc_Nlim]+IMM_N_max[loc_Nlim] - microbeTurnover[loc_Nlim]*et/params['CN_microbe'])
    overflow_resp[loc_Nlim]=carbon_supply[loc_Nlim]-maintenance_resp[loc_Nlim] - (nitrogen_supply[loc_Nlim]+IMM_N_max[loc_Nlim])*params['CN_microbe'];

    # Growth must be supported by immobilization of some mineral nitrogen, but is ultimately carbon limited
    loc_immob=(carbon_supply - maintenance_resp >= nitrogen_supply*params['CN_microbe']) & (carbon_supply - maintenance_resp < (nitrogen_supply+IMM_N_max)*params['CN_microbe'])
    CN_imbalance_term[loc_immob] = -((carbon_supply[loc_immob]-maintenance_resp[loc_immob])/params['CN_microbe'] - nitrogen_supply[loc_immob])
    dmicrobeC[loc_immob] = (carbon_supply[loc_immob] - microbeTurnover[loc_immob])
    dmicrobeN[loc_immob] = ((carbon_supply[loc_immob]-maintenance_resp[loc_immob])/params['CN_microbe'] - microbeTurnover[loc_immob]*et/params['CN_microbe'])

    # Growth is carbon limited and extra N is mineralized
    loc_Clim=~(loc_Nlim | loc_immob)
    dmicrobeC[loc_Clim]=(carbon_supply[loc_Clim] - microbeTurnover[loc_Clim]) ;
    dmicrobeN[loc_Clim]=((carbon_supply[loc_Clim]-maintenance_resp[loc_Clim])/params['CN_microbe'] - microbeTurnover[loc_Clim]*params['et']/params['CN_microbe'])
    CN_imbalance_term[loc_Clim] = nitrogen_supply[loc_Clim] - (carbon_supply[loc_Clim]-maintenance_resp[loc_Clim])/params['CN_microbe'];

    # CO2 production and cumulative CO2 produced by cohort
    CO2prod=maintenance_resp+overflow_resp
    for t in chem_types:
        CO2prod=CO2prod+decomp[t+'C']*(1.0-eup[t])


    # Update protected carbon
    protectedCturnover = dict([(t,SOM['p'+t+'C']/params['tProtected']) for t in chem_types])
    protectedCprod =     dict([(t,SOM['u'+t+'C']*params['protection_rate'][t]*claymod) for t in chem_types])
    protectedNturnover = dict([(t,SOM['p'+t+'N']/params['tProtected']) for t in chem_types])
    protectedNprod =     dict([(t,SOM['u'+t+'N']*params['protection_rate'][t]*claymod) for t in chem_types])

    derivs=SOM.copy()
    for k in derivs.keys():
        derivs[k]=0.0
    derivs['livingMicrobeC']=dmicrobeC
    derivs['livingMicrobeN']=dmicrobeN
    derivs['CO2']=CO2prod
    derivs['inorganicN']=CN_imbalance_term
    for t in chem_types:
        derivs['inorganicN'] += decomp[t+'N']*(1-params['nup'][t])

    for t in chem_types:
        derivs['u'+t+'C']=-decomp[t+'C']+protectedCturnover[t]-protectedCprod[t]
        derivs['p'+t+'C']=protectedCprod[t]-protectedCturnover[t]
        derivs['u'+t+'N']=-decomp[t+'N']+protectedNturnover[t]-protectedNprod[t]
        derivs['p'+t+'N']=protectedNprod[t]-protectedNturnover[t]

    derivs['uNecroC']=derivs['uNecroC']+deadmic_C_production*(1.0-params['frac_turnover_slow'])
    derivs['uSlowC']=derivs['uSlowC']+deadmic_C_production*params['frac_turnover_slow']
    turnover_N_min=deadmic_N_production*params['frac_N_turnover_min'];
    turnover_N_slow=deadmic_N_production*params['frac_turnover_slow'];
    derivs['uNecroN']=derivs['uNecroN']+deadmic_N_production-turnover_N_min-turnover_N_slow
    derivs['uSlowN']+=turnover_N_slow
    derivs['inorganicN']+=turnover_N_min
    return derivs


# Decomposition rate
def decompRate(SOM,T,theta,params):

    # This only really needs to be calculated once
    if params['new_resp_units']:
        theta_resp_max=params['substrate_diffusion_exp']/(params['gas_diffusion_exp']*(1.0+params['substrate_diffusion_exp']/params['gas_diffusion_exp']))
        aerobic_max=theta_resp_max**params['substrate_diffusion_exp']*(1.0-theta_resp_max)**params['gas_diffusion_exp']

    else:
        aerobic_max=1.0

    vmax=Vmax(T,params)

    decompRate={}
    dodecomp=atleast_1d((sumCtypes(SOM,'u')!=0.0)&(theta!=0.0)&(SOM['livingMicrobeC']!=0.0))
    for t in chem_types:
        if dodecomp.any():
            drate=where(dodecomp,vmax[t]*theta**params['substrate_diffusion_exp']*(SOM['u'+t+'C'])*SOM['livingMicrobeC']/(sumCtypes(SOM,'u')*params['kC'][t]+SOM['livingMicrobeC'])*(1.0-theta)**params['gas_diffusion_exp']/aerobic_max,0.0)
        decompRate[t+'C']=drate
        decompRate[t+'N']=where(SOM['u'+t+'C']>0,drate*SOM['u'+t+'N']/SOM['u'+t+'C'],0.0)

    return decompRate

def Vmax(T,params):
    '''Vmax function, normalized to Tref=293.15
    T is in K'''

    Tref=293.15;
    Rugas=8.314472;

    from numpy import exp

    Vmax=dict([(t,params['vmaxref'][t]*exp(-params['Ea'][t]*(1.0/(Rugas*T)-1.0/(Rugas*Tref)))) for t in chem_types]);
    return Vmax

def sumCtypes(SOM,prefix,suffix='C'):
    out=SOM[prefix+chem_types[0]+suffix]
    if len(chem_types)>1:
        for t in chem_types[1:]:
            out=out+SOM[prefix+t+suffix]

    return out
