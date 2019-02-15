import pandas,numpy

# Code for integrating the model using python ODE solver
import CORPSE_integrate

# Initial values for all pools in the model
# The {} brackets create a "dictionary" where everything is indexed by a name rather than a number
# Pools starting with "u" are unprotected and pools starting with "p" are protected
# These are assumed to be ugC/g mineral as in the incubation. But the units are not that important as long as they are
# consistent throughout

SOM_init={'uFastC':[0.0,0.0],
     'uSlowC':[0.0,0.0],
     'uNecroC':[0.0,440.0], # From the paper, 440 ugC/g mineral of 13C enriched residues
     'pFastC':[0.0,0.0],
     'pSlowC':[0.0,0.0],
     'pNecroC':[0.0,0.0],
     'livingMicrobeC':[100.0,0.0], # From the paper, 100 ug E. coli/g mineral
     'CO2':[0.0,0.0]}

# Set model parameters
# Note that carbon types have names, in contrast to previous version
params={
  'vmaxref':{'Fast':9.0,'Slow':.25,'Necro':10.5}, #Relative maximum enzymatic decomp rates (year-1)
  'Ea':{'Fast':5e3,'Slow':30e3,'Necro':3e3},    # Activation energy (Controls temperature sensitivity via Arrhenius relationship)
  'kC':{'Fast':0.01,'Slow':0.01,'Necro':0.01},    # Michaelis-Menton parameter
  'gas_diffusion_exp':0.6,  # Determines suppression of decomp at high soil moisture
  'substrate_diffusion_exp':1.5, # Dimensionless exponent. Determines suppression of decomp at low soil moisture
  'minMicrobeC':1e-3,       #Minimum microbial biomass (fraction of total C)
  'Tmic':0.25,       # Microbial lifetime (years)
  'et':0.6,          # Fraction of turnover not converted to CO2 (dimensionless)
  'eup':{'Fast':0.6,'Slow':0.05,'Necro':0.6}, # Carbon uptake efficiency (dimensionless fraction)
  'tProtected':75.0,    # Protected C turnover time (years)
  'protection_rate':{'Fast':0.3,'Slow':0.001,'Necro':20.5}, # Protected carbon formation rate (year-1)
  'new_resp_units':True
}

# Environmental conditions
T=20        # degrees C
theta=0.5   # fraction of saturation
inputs={'uFastC':[0.0,0.0],'uSlowC':[0.0,0.0]} # gC/year. Can contain any model pools. Setting to zero for now.
times=numpy.arange(0,30,1)  # Time steps to evaluate. Defining these in day units but need to convert to years for actual model simulations.
                      # The ODE solver uses an adaptive timestep but will return these time points

# We should set these parameters in a more logical way, but for now using different %clay values to differentiate protection capacity of the minerals
clay_feldspar=1.0
clay_AlOH=90.0

# Run the simulations
# Living microbes, lower protection capacity
result_mic_feld = CORPSE_integrate.run_CORPSE_ODE(T=T,theta=theta,inputs=inputs,clay=clay_feldspar,initvals=SOM_init,params=params,times=times/365)
# No living microbes
SOM_init_nomicrobes=SOM_init.copy()
SOM_init_nomicrobes['livingMicrobeC']=[0.0,0.0]
result_nomic_feld = CORPSE_integrate.run_CORPSE_ODE(T=T,theta=theta,inputs=inputs,clay=clay_feldspar,initvals=SOM_init_nomicrobes,params=params,times=times/365)

# Living microbes, higher protection capacity
result_mic_AlOH = CORPSE_integrate.run_CORPSE_ODE(T=T,theta=theta,inputs=inputs,clay=clay_AlOH,initvals=SOM_init,params=params,times=times/365)
# No living microbes, higher protection capacity
result_nomic_AlOH = CORPSE_integrate.run_CORPSE_ODE(T=T,theta=theta,inputs=inputs,clay=clay_AlOH,initvals=SOM_init_nomicrobes,params=params,times=times/365)

SOM_init_nores=SOM_init.copy()
SOM_init_nores['uNecroC']=[0.0,0.0]
result_mic_nores_feld = CORPSE_integrate.run_CORPSE_ODE(T=T,theta=theta,inputs=inputs,clay=clay_feldspar,initvals=SOM_init_nores,params=params,times=times/365)
result_mic_nores_AlOH = CORPSE_integrate.run_CORPSE_ODE(T=T,theta=theta,inputs=inputs,clay=clay_AlOH,initvals=SOM_init_nores,params=params,times=times/365)

# Plot the results
import matplotlib.pyplot as plt
from CORPSE_deriv import totalCarbon

# This figure emulates Fig. 1 in the manuscript main text
plt.figure('Microbial biomass and respiration',figsize=(4,5.3));plt.clf()

plt.subplot(211)
plt.plot(times,result_mic_feld[1]['CO2']/totalCarbon(result_mic_feld[1].iloc[0])*100,c='b',label='Feldspar, E. coli')
# plt.plot(times,result_nomic_feld[1]['CO2']/totalCarbon(result_nomic_feld[1].iloc[0])*100,c='C1',ls='--',label='Feldspar, no E. coli')
plt.plot(times,result_mic_AlOH[1]['CO2']/totalCarbon(result_mic_AlOH[1].iloc[0])*100,c='r',label='am-Al(OH)$_3$, E. coli')
# plt.plot(times,result_nomic_AlOH[1]['CO2']/totalCarbon(result_nomic_AlOH[1].iloc[0])*100,c='C2',ls='--',label='AlOH, no E. coli')
plt.xlabel('Time (days)')
plt.ylabel('Cumulative $^{13}$C-CO$_2$ (%)')
plt.title('Respired C')
plt.legend()

# Currently not directly comparable to manuscript figure, which shows isotope labeled % of microbial biomass, not total microbial biomass
plt.subplot(212)
plt.plot(times,result_mic_feld[1]['livingMicrobeC']/(result_mic_feld[0]['livingMicrobeC']+result_mic_feld[1]['livingMicrobeC'])*100,c='b')
# plt.plot(times,result_nomic_feld[1]['livingMicrobeC']/(result_nomic_feld[0]['livingMicrobeC']+result_nomic_feld[1]['livingMicrobeC'])*100,c='C1',ls='--')
plt.plot(times,result_mic_AlOH[1]['livingMicrobeC']/(result_mic_AlOH[0]['livingMicrobeC']+result_mic_AlOH[1]['livingMicrobeC'])*100,c='r')
# plt.plot(times,result_nomic_AlOH[1]['livingMicrobeC']/(result_nomic_AlOH[0]['livingMicrobeC']+result_nomic_AlOH[1]['livingMicrobeC'])*100,c='C2',ls='--')
plt.xlabel('Time (days)')
plt.ylabel('Microbial biomass $^{13}$C (%)')
plt.title('Microbial biomass C')

plt.tight_layout()


# This figure emulates Fig. S3 in the manuscript
plt.figure('E coli C respired');plt.clf()
plt.plot(times,result_mic_feld[0]['CO2']/totalCarbon(result_mic_feld[0])*100,c='b',label='Residue and E. coli with feldspar')
plt.plot(times,result_mic_nores_feld[0]['CO2']/totalCarbon(result_mic_nores_feld[0])*100,c='c',label='E. coli control with feldspar')
plt.plot(times,result_mic_AlOH[0]['CO2']/totalCarbon(result_mic_AlOH[0])*100,c='r',label='Residue and E. coli with am-Al(OH)$_3$')
plt.plot(times,result_mic_nores_AlOH[0]['CO2']/totalCarbon(result_mic_nores_AlOH[0])*100,c='orange',label='E. coli control with am-Al(OH)$_3$')

plt.xlabel('Time (days)')
plt.ylabel('E. coli C-CO$_2$ respired (%)')
plt.legend(fontsize='small')


# Figure showing all carbon pools: microbial, protected, necromass
plt.figure('Carbon pools',figsize=(4.5,9));plt.clf()
plt.subplot(311)
plt.plot(times,result_mic_feld[0]['livingMicrobeC'],'b-',label='Feldspar $^{12}$C')
plt.plot(times,result_mic_feld[1]['livingMicrobeC'],'b--',label='Feldspar $^{13}$C')
plt.plot(times,result_nomic_feld[1]['livingMicrobeC'],'b:',label='Feldspar $^{13}$C (no microbes)')
plt.plot(times,result_mic_AlOH[0]['livingMicrobeC'],'r-',label='am-Al(OH)$_3$ $^{12}$C')
plt.plot(times,result_mic_AlOH[1]['livingMicrobeC'],'r--',label='am-Al(OH)$_3$ $^{13}$C')
plt.plot(times,result_nomic_AlOH[1]['livingMicrobeC'],'r:',label='am-Al(OH)$_3$ $^{13}$C (no microbes)')
plt.title('Live microbe C')
plt.ylabel('Microbe C ($\mu$g C g mineral$^{-1}$)')
plt.ylim(-5,460)
plt.legend(fontsize='small')

plt.subplot(312)
plt.plot(times,result_mic_feld[0]['uNecroC'],'b-')
plt.plot(times,result_mic_feld[1]['uNecroC'],'b--')
plt.plot(times,result_nomic_feld[1]['uNecroC'],'b:')
plt.plot(times,result_mic_AlOH[0]['uNecroC'],'r-')
plt.plot(times,result_mic_AlOH[1]['uNecroC'],'r--')
plt.plot(times,result_nomic_AlOH[1]['uNecroC'],'r:')
plt.title('Necromass C')
plt.ylabel('Necromass C ($\mu$g C g mineral$^{-1}$)')
plt.ylim(-5,460)

plt.subplot(313)
plt.plot(times,result_mic_feld[0]['pNecroC'],'b-')
plt.plot(times,result_mic_feld[1]['pNecroC'],'b--')
plt.plot(times,result_nomic_feld[1]['pNecroC'],'b:')
plt.plot(times,result_mic_AlOH[0]['pNecroC'],'r-')
plt.plot(times,result_mic_AlOH[1]['pNecroC'],'r--')
plt.plot(times,result_nomic_AlOH[1]['pNecroC'],'r:')
plt.title('Mineral-associated C')
plt.ylabel('Mineral-associated C ($\mu$g C g mineral$^{-1}$)')
plt.ylim(-5,460)
plt.xlabel('Time (days)')

plt.tight_layout()



plt.show()
