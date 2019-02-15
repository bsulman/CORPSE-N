
import CORPSE_deriv
import pandas

fields=CORPSE_deriv.expected_pools

# This is a function that translates the CORPSE model pools to/from the format that the equation solver expects
# The solver will call it multiple times and passes it a list of parameters that needs to be converted to a named "dictionary" that CORPSE expects
def fsolve_wrapper(SOM_list,T,theta,inputs,clay,params):
    from numpy import asarray

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
        vals=[deriv[f] for f in fields]

    return vals

# The ordinary differential equation (ODE) integrating function also wants to send the current time to the function it's integrating
# Our model doesn't have an explicit dependence on time, but we need a separate function that can deal with the extra argument.
# We just ignore the time argument and pass the rest to the same function we used for the numerical solver
def ode_wrapper(SOM_list,time,T,theta,inputs,clay,params):
    return fsolve_wrapper(SOM_list,T,theta,inputs,clay,params)


def run_CORPSE_ODE(T,theta,inputs,clay,initvals,params,times):
    # Use ODE integrator to actually integrate the model
    # import time
    # t0=time.time()
    from scipy.integrate import odeint

    if len(initvals['livingMicrobeC'])==2:
        # With isotope label, concatenate one after the other, put back together later
        ivals=[initvals[f][0] for f in fields]+[initvals[f][1] for f in fields]
    else:
        ivals=[initvals[f] for f in fields]

    # Runs the ODE integrator
    result=odeint(ode_wrapper,ivals,times,
        args=(T+273.15,theta,inputs,clay,params))

    # Store the output in a pandas DataFrame (similar to R's dataframes)
    if len(initvals['livingMicrobeC'])==2:
        result_unlabeled=pandas.DataFrame(result[:,:len(fields)],columns=fields)
        result_labeled=pandas.DataFrame(result[:,len(fields):],columns=fields)
        return result_unlabeled,result_labeled
    else:
        result_df=pandas.DataFrame(result[:,:len(fields)],columns=fields)
        return result_df

    # print('Time elapsed: %1.1f s'%(time.time()-t0))
