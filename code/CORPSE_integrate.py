
import CORPSE_deriv
import pandas

fields=CORPSE_deriv.expected_pools

# This is a function that translates the CORPSE model pools to/from the format that the equation solver expects
# The solver will call it multiple times and passes it a list of parameters that needs to be converted to a named "dictionary" that CORPSE expects
def fsolve_wrapper(SOM_list,T,theta,inputs,clay,params):
    from numpy import asarray,concatenate

    # Make an empty dictionary and fill it with the right values
    # Values are in the correct order because we use the same "fields" list of pools
    SOM_dict={}
    if len(SOM_list)==len(fields)*2:
        # Running in isotope labeled mode
        for n in range(len(fields)):
            SOM_dict[fields[n]]=asarray([SOM_list[n],SOM_list[n+len(fields)]])
    else:
        for n in range(len(fields)):
            SOM_dict[fields[n]]=asarray(SOM_list[n])

    # Call the CORPSE model function that returns the derivative (with time) of each pool
    deriv=CORPSE_deriv.CORPSE_deriv(SOM_dict,T,theta,params,claymod=CORPSE_deriv.prot_clay(clay)/CORPSE_deriv.prot_clay(20))

    for pool in inputs.keys():
        deriv[pool]+=inputs[pool]

    # Put the CORPSE pools back into a list that the equation solver can deal with
    if len(SOM_list)==len(fields)*2:
        # With isotope label, concatenate one after the other, put back together later
        vals=[deriv[f][0] for f in fields]+[deriv[f][1] for f in fields]
    else:
        vals=concatenate([deriv[f] for f in fields])

    return vals

# The ordinary differential equation (ODE) integrating function also wants to send the current time to the function it's integrating
# Our model doesn't have an explicit dependence on time, but we need a separate function that can deal with the extra argument.
# We just ignore the time argument and pass the rest to the same function we used for the numerical solver
def ode_wrapper(SOM_list,time,T,theta,inputs,clay,params):
    return fsolve_wrapper(SOM_list,T,theta,inputs,clay,params)

def arrayify_dict(d):
    from numpy import atleast_1d
    return dict(((v,atleast_1d(d[v])) for v in d))

def run_CORPSE_ODE(T,theta,inputs,clay,initvals,params,times):
    # Use ODE integrator to actually integrate the model. Currently set up for constant temperature, moisture, and inputs
    # import time
    # t0=time.time()
    from scipy.integrate import odeint

    if not isinstance(initvals['livingMicrobeC'],float) and len(initvals['livingMicrobeC'])==2:
        # With isotope label, concatenate one after the other, put back together later
        ivals=[initvals[f][0] for f in fields]+[initvals[f][1] for f in fields]
    else:
        ivals=[initvals[f] for f in fields]

    # Runs the ODE integrator
    result=odeint(ode_wrapper,ivals,times,
        args=(T+273.15,theta,inputs,clay,params))

    # Store the output in a pandas DataFrame (similar to R's dataframes)
    if not isinstance(initvals['livingMicrobeC'],float) and len(initvals['livingMicrobeC'])==2:
        result_unlabeled=pandas.DataFrame(result[:,:len(fields)],columns=fields)
        result_labeled=pandas.DataFrame(result[:,len(fields):],columns=fields)
        return result_unlabeled,result_labeled
    else:
        result_df=pandas.DataFrame(result[:,:len(fields)],columns=fields)
        return result_df

    # print('Time elapsed: %1.1f s'%(time.time()-t0))

def run_CORPSE_iterator(T,theta,inputs,clay,initvals,params,times):
    # Run model with explicit timestepping. This allows arbitrary time series of T, theta, and inputs but may be slower/less accurate depending on time step length
    # Allows running a vector of points together which can be more efficient
    # To represent multiple soil compartments such as rhizosphere, litter layer, multiple soil layers, etc: Run with a vector of cells representing the different compartments
    #  This requires C and N inputs to be properly divided across cells, and may also require an added step to transfer C and N across compartments if there is mixing.
    from numpy import zeros,atleast_1d,floor

    # Check number of cells being integrated
    nsteps=len(times)
    if len(T.shape)>1: # Assume temperature is being specified for each cell individually
        npoints=T.shape[1]
    else:
        npoints=1
    nrecords=nsteps

    # Set up pools
    SOM={}
    SOM_out={}
    for field in initvals.keys():
        SOM_out[field]=zeros((npoints,nrecords))
        if len(atleast_1d(initvals['uFastC'])) == 1: 
            SOM[field]=zeros(npoints)+initvals[field]
        else:
            SOM[field]=initvals[field]
    
    SOM['livingMicrobeN']=SOM['livingMicrobeC']/params['CN_microbe']
    SOM_out['livingMicrobeN']=zeros((npoints,nrecords))

    # Iterate through simulations
    for step in range(nsteps):
        if step==nsteps-1:
            dt=times[step]-times[step-1]
        else:
            dt=times[step+1]-times[step]
        if len(T.shape)>1:
            T_step=T[step,:]
        else:
            T_step=T[step]
        if len(theta.shape)>1:
            theta_step=theta[step,:]
        else:
            theta_step=theta[step]
        # In this case, T, theta, clay, and all the pools in SOM are vectors containing one value per geographical location
        deriv=CORPSE_deriv.CORPSE_deriv(SOM,T_step+273.15,theta_step,params,claymod=CORPSE_deriv.prot_clay(clay)/CORPSE_deriv.prot_clay(20))

        # Inputs and N uptake calculations below can be improved in the future to include plant C allocation to mycorrhizae (via FUN), plant N uptake and demand, etc

        # Inorganic N is lost at some fixed rate that accounts for things like leaching, denitrification, and plant uptake
        deriv['inorganicN']-=SOM['inorganicN']*params['iN_loss_rate']

        # Since we have carbon/nitrogen inputs, these also need to be added to those rates of change with time
        for pool in inputs.keys():
            if len(atleast_1d(inputs[pool]).shape)>1:
                deriv[pool]+=inputs[pool][step,:]
            else:
                deriv[pool]+=inputs[pool]

        # Here we update the pools, using a simple explicit time step calculation
        for field in deriv.keys():
            SOM[field]=SOM[field]+deriv[field]*dt

        # Print some output to track progress
        if floor((step*dt*365))%365==0:
            print('Time = %d'%(step*dt))
        for field in SOM.keys():
            SOM_out[field][:,step]=SOM[field]

    return SOM_out

if __name__ == '__main__':
    # Test simulation
    import numpy
        
    SOM_init={'uFastC':0.1,
        'uSlowC':10.0,
        'uNecroC':0.1,
        'pFastC':0.1,
        'pSlowC':0.1,
        'pNecroC':10.0,
        'livingMicrobeC':0.01,
        'uFastN':0.1e-1,
        'uSlowN':20.0e-1,
        'uNecroN':0.1e-1,
        'pFastN':10.0e-1,
        'pSlowN':0.1e-1,
        'pNecroN':10.0e-1,
        # 'livingMicrobeN':0.01/8.0,
        'inorganicN':0.1,
        'CO2':0.0}


    # Set model parameters
    # Note that carbon types have names, in contrast to previous version
    params={
        'vmaxref':{'Fast':9.0,'Slow':.5,'Necro':4.5}, #Relative maximum enzymatic decomp rates (year-1)
        'Ea':{'Fast':5e3,'Slow':30e3,'Necro':3e3},    # Activation energy (Controls temperature sensitivity via Arrhenius relationship)
        'kC':{'Fast':0.01,'Slow':0.01,'Necro':0.01},    # Michaelis-Menton parameter
        'gas_diffusion_exp':0.6,  # Determines suppression of decomp at high soil moisture
        'substrate_diffusion_exp':1.5, # Dimensionless exponent. Determines suppression of decomp at low soil moisture
        'minMicrobeC':1e-3,       #Minimum microbial biomass (fraction of total C)
        'Tmic':0.25,       # Microbial lifetime (years)
        'et':0.6,          # Fraction of turnover not converted to CO2 (dimensionless)
        'eup':{'Fast':0.6,'Slow':0.05,'Necro':0.6}, # Carbon uptake efficiency (dimensionless fraction)
        'tProtected':75.0,    # Protected C turnover time (years)
        'protection_rate':{'Fast':0.1,'Slow':0.0001,'Necro':1.5}, # Protected carbon formation rate (year-1)
        'new_resp_units':True,
        'frac_N_turnover_min':0.2,
        'frac_turnover_slow':0.2,
        'nup':{'Fast':0.9,'Slow':0.6,'Necro':0.9},
        'CN_microbe':8.0,
        'max_immobilization_rate':3.65,
        'substrate_diffusion_exp':1.5,
        'new_resp_units':True,
        'iN_loss_rate':10.0, # Loss rate from inorganic N pool (year-1). >1 since it takes much less than a year for it to be removed
        'Ohorizon_transfer_rates':{'uFastC':0.1,'uSlowC':0.1,'uNecroC':0.1,'uFastN':0.1,'uSlowN':0.1,'uNecroN':0.1}
    }
    SOM_init['livingMicrobeN']=SOM_init['livingMicrobeC']/params['CN_microbe']

    times=numpy.linspace(0,10,365*10) # Units of years. 10 years at daily time step
    T=numpy.zeros(len(times))+20
    theta=numpy.zeros(len(times))+0.6
    clay=numpy.atleast_1d(30)

    total_inputs = 1.0 # kgC/m2
    fastfrac_AM  = 0.4
    fastfrac_ECM = 0.1
    ECM_pct=50
    fastfrac_site = fastfrac_ECM*ECM_pct/100 + fastfrac_AM*(1-ECM_pct/100)
    litter_CN_AM = 30
    litter_CN_ECM= 50
    litter_CN_site      = litter_CN_ECM*ECM_pct/100 + litter_CN_AM*(1-ECM_pct/100)
    inputs={'uFastC':total_inputs*fastfrac_site,
        'uSlowC':total_inputs*(1-fastfrac_site),
        'uFastN':total_inputs*fastfrac_site/litter_CN_site,
        'uSlowN':total_inputs*(1-fastfrac_site)/litter_CN_site} # Can contain any model pools. Pools not included have zero inputs

    results_iterator = run_CORPSE_iterator(T=T,theta=theta,inputs=inputs,clay=clay,initvals=SOM_init,params=params,times=times)
    results_ODE      = run_CORPSE_ODE     (T=T[0],theta=theta[0],inputs=inputs,clay=clay,initvals=SOM_init,params=params,times=times)

    import matplotlib.pyplot as plt
    f,a=plt.subplots(2,2,clear=True,num='Results')
    # The "squeeze" removes an extra empty dimension from the output of the iterator function that the plot command doesn't like
    a[0,0].plot(times,results_iterator['uSlowC'].squeeze(),'b-',label='Slow unprotected')
    a[0,0].plot(times,CORPSE_deriv.sumCtypes(results_iterator,'p').squeeze(),'g-',label='Protected')
    # Uncomment these lines to compare the two integration methods
    # a[0,0].plot(times,results_ODE['uSlowC'],'r--')
    # a[0,0].plot(times,CORPSE_deriv.sumCtypes(results_ODE,'p'),'r--')
    a[0,0].set(title='Slower C pools',xlabel='Time (years)',ylabel='C stock (kg m$^{-2}$)')
    a[0,0].legend()

    a[1,0].plot(times,results_iterator['uSlowN'].squeeze(),'b-',label='Slow unprotected')
    a[1,0].plot(times,CORPSE_deriv.sumCtypes(results_iterator,'p','N').squeeze(),'g-',label='Protected')
    # a[1,0].plot(times,results_ODE['uSlowN'],'r--')
    # a[1,0].plot(times,CORPSE_deriv.sumCtypes(results_ODE,'p','N'),'r--')
    a[1,0].set(title='Slower N pools',xlabel='Time (years)',ylabel='N stock (kg m$^{-2}$)')

    a[0,1].plot(times,results_iterator['uFastC'].squeeze(),'b-',label='Fast')
    a[0,1].plot(times,results_iterator['livingMicrobeC'].squeeze(),'g-',label='Live microbe')
    a[0,1].plot(times,results_iterator['uNecroC'].squeeze(),'m-',label='Necromass')
    a[0,1].set(title='Faster C pools',xlabel='Time (years)',ylabel='C stock (kg m$^{-2}$)')
    a[0,1].legend()

    a[1,1].plot(times,results_iterator['uFastN'].squeeze(),'b-',label='Fast')
    a[1,1].plot(times,results_iterator['livingMicrobeN'].squeeze(),'g-',label='Live microbe')
    a[1,1].plot(times,results_iterator['uNecroN'].squeeze(),'m-',label='Necromass')
    a[1,1].set(title='Faster N pools',xlabel='Time (years)',ylabel='N stock (kg m$^{-2}$)')


    plt.show()