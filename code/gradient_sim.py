import pandas,numpy

# Code for integrating the model using python ODE solver
import CORPSE_integrate

# Initial values for all pools in the model
# The {} brackets create a "dictionary" where everything is indexed by a name rather than a number
# Pools starting with "u" are unprotected and pools starting with "p" are protected
# These are assumed to be ugC/g mineral as in the incubation. But the units are not that important as long as they are
# consistent throughout

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

# ECM gradient plots
nplots=20
nclays=2
nclimates=3
# Environmental conditions
# Gradient of mycorrhizal association
ECM_pct  = numpy.linspace(0,100,nplots)          # Percent ECM basal area
MAT      = numpy.linspace(5,20,nclimates)      # degrees C
clay     = numpy.linspace(10,70,nclays)        # percent clay

fastfrac_AM  = 0.4
fastfrac_ECM = 0.1
fastfrac_site = fastfrac_ECM*ECM_pct/100 + fastfrac_AM*(1-ECM_pct/100)
litter_CN_AM = 30
litter_CN_ECM= 50
litter_CN_site      = litter_CN_ECM*ECM_pct/100 + litter_CN_AM*(1-ECM_pct/100)
total_inputs = 0.5 # kgC/m2
inputs={'uFastC':total_inputs*fastfrac_site,
        'uSlowC':total_inputs*(1-fastfrac_site),
        'uFastN':total_inputs*fastfrac_site/litter_CN_site,
        'uSlowN':total_inputs*(1-fastfrac_site)/litter_CN_site} # gC/year. Can contain any model pools.


theta=0.5   # fraction of saturation

times=numpy.arange(0,1000,10)  # Time steps to evaluate. Defining these in day units but need to convert to years for actual model simulations.
                      # The ODE solver uses an adaptive timestep but will return these time points



# Run the simulations
protC=numpy.zeros((nplots,nclays,nclimates))
protN=numpy.zeros((nplots,nclays,nclimates))
unprotC=numpy.zeros((nplots,nclays,nclimates))
unprotN=numpy.zeros((nplots,nclays,nclimates))
n=0
from CORPSE_deriv import sumCtypes
for plotnum in range(nplots):
    for claynum in range(nclays):
        for climnum in range(nclimates):
            print('Sim {simnum:d} of {totsims:d}. %ECM = {ecmpct:1.1f}, %clay = {claypct:1.1f}, MAT = {mat:1.1f}'.format(
                    simnum=n,totsims=nplots*nclays*nclimates,ecmpct=ECM_pct[plotnum],claypct=clay[claynum],mat=MAT[climnum]))
            result = CORPSE_integrate.run_CORPSE_ODE(T=MAT[climnum],theta=theta,inputs=dict([(k,inputs[k][plotnum]) for k in inputs]),clay=clay[claynum],initvals=SOM_init,params=params,times=times)
            protC[plotnum,claynum,climnum]=sumCtypes(result.iloc[-1],'p')
            protN[plotnum,claynum,climnum]=sumCtypes(result.iloc[-1],'p','N')
            unprotC[plotnum,claynum,climnum]=sumCtypes(result.iloc[-1],'u')
            unprotN[plotnum,claynum,climnum]=sumCtypes(result.iloc[-1],'u','N')
            n+=1

# Plot the results
import matplotlib.pyplot as plt

def totalCarbon(SOM):
    from CORPSE_deriv import sumCtypes
    return sumCtypes(SOM,'u')+sumCtypes(SOM,'p')+SOM['livingMicrobeC']


def totalNitrogen(SOM):
    from CORPSE_deriv import sumCtypes
    return sumCtypes(SOM,'u','N')+sumCtypes(SOM,'p','N')+SOM['livingMicrobeN']


plt.figure('C and N for one sim',figsize=(4,5.3));plt.clf()

plt.subplot(211)
plt.plot(times,totalCarbon(result),c='k',label='Total C')
plt.plot(times,result['pNecroC'],label='pNecroC')
plt.plot(times,result['uSlowC'],label='uSlowC')

plt.xlabel('Time (days)')
plt.ylabel('Total carbon')
plt.title('Total C stock')
plt.legend(fontsize='small')

plt.subplot(212)
plt.plot(times,totalNitrogen(result),c='k',label='Total N')
plt.plot(times,result['pNecroN'],label='pNecroN')
plt.plot(times,result['uSlowN'],label='uSlowN')

plt.xlabel('Time (days)')
plt.ylabel('Total nitrogen')
plt.title('Total N')
plt.legend(fontsize='small')



protCfrac=protC/(protC+unprotC)
protNfrac=protN/(protN+unprotN)

norm=plt.Normalize(5,20)
cmap=plt.get_cmap('cool')
markers=['o','s']

plt.figure('Protected fraction of C and N',figsize=(6,8));plt.clf()
plt.subplot(211)
for claynum in [0,1]:
    for climnum in range(len(MAT)):
        plt.plot(ECM_pct,protCfrac[:,claynum,climnum],marker=markers[claynum],c=cmap(norm(MAT[climnum])),label='Clay={claypct:1.1f}%, MAT={mat:1.1f}C'.format(claypct=clay[claynum],mat=MAT[climnum]))

plt.xlabel('ECM percent (%)')
plt.ylabel('Protected C fraction')
plt.legend(fontsize='small')
plt.title('Protected SOM C fraction')

plt.subplot(212)
for claynum in [0,1]:
    for climnum in range(len(MAT)):
        plt.plot(ECM_pct,protNfrac[:,claynum,climnum],marker=markers[claynum],c=cmap(norm(MAT[climnum])),label='Clay={claypct:1.1f}%, MAT={mat:1.1f}C'.format(claypct=clay[claynum],mat=MAT[climnum]))

plt.xlabel('ECM percent (%)')
plt.ylabel('Protected N fraction')
# plt.legend()
plt.title('Protected SOM N fraction')





plt.figure('N stock and C:N',figsize=(6,8));plt.clf()
plt.subplot(311)
for claynum in [0,1]:
    for climnum in range(len(MAT)):
        plt.plot(ECM_pct,(protC+unprotC)[:,claynum,climnum],ms=4,marker=markers[claynum],c=cmap(norm(MAT[climnum])),label='Clay={claypct:1.1f}%, MAT={mat:1.1f}C'.format(claypct=clay[claynum],mat=MAT[climnum]))
        plt.plot(ECM_pct,(protC)[:,claynum,climnum],ms=4,marker=markers[claynum],mfc='w',c=cmap(norm(MAT[climnum])))


plt.xlabel('ECM percent (%)')
plt.ylabel('Total C stock')
# plt.legend()
plt.title('Total C stock')

plt.subplot(312)
for claynum in [0,1]:
    for climnum in range(len(MAT)):
        plt.plot(ECM_pct,(protN+unprotN)[:,claynum,climnum],ms=4,marker=markers[claynum],c=cmap(norm(MAT[climnum])),label='Clay={claypct:1.1f}%, MAT={mat:1.1f}C'.format(claypct=clay[claynum],mat=MAT[climnum]))

plt.xlabel('ECM percent (%)')
plt.ylabel('Total N stock')
# plt.legend()
plt.title('Total N stock')

plt.subplot(313)
for claynum in [0,1]:
    for climnum in range(len(MAT)):
        plt.plot(ECM_pct,((protC+unprotC)/(protN+unprotN))[:,claynum,climnum],ms=4,marker=markers[claynum],c=cmap(norm(MAT[climnum])),label='Clay={claypct:1.1f}%, MAT={mat:1.1f}C'.format(claypct=clay[claynum],mat=MAT[climnum]))

plt.xlabel('ECM percent (%)')
plt.ylabel('C:N ratio')
# plt.legend()
plt.title('C:N ratio')
plt.legend(fontsize='small')


plt.figure('Myco effect vs decomp rate');plt.clf()
protCfracdiff = protCfrac[-1,:,:]-protCfrac[0,:,:]
# plt.scatter(unprotC/total_inputs,protCfrac)
for claynum in [0,1]:
    for climnum in range(len(MAT)):
        plt.plot(unprotC[:,claynum,climnum]/total_inputs,protCfrac[:,claynum,climnum],marker=markers[claynum],c=cmap(norm(MAT[climnum])),label='Clay={claypct:1.1f}%, MAT={mat:1.1f}C'.format(claypct=clay[claynum],mat=MAT[climnum]))

plt.xlabel('Unprotected C turnover time (years)')
plt.ylabel('Protected C fraction')



plt.show()
