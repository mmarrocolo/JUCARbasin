import numpy       as np
from Weigfun import Weigfun

def hillslope_func(timestep, Par, forcing, Fluxes, States):
	# HBVpareto Calculates values of 3 objective functions for HBV model

	I_max = Par[0]
	Ce = Par[1]
	Su_max = Par[2]
	beta = Par[3]
	D = Par[4]
	Kf = Par[5]
	Qo = forcing['Qo']
	Prec = forcing['Prec']
	Etp = forcing['Etp']


	t_max = len(Prec)
	Si = States[:,0]
	Su = States[:,1]
	Sf = States[:,2]

	Ei_dt = Fluxes[:,0]
	Ea_dt = Fluxes[:,1]
	Qf_dt = Fluxes[:,2]
	Qus_dt = Fluxes[:,3]

	dt = 1
	i = timestep
	P_dt = Prec[i] * dt
	Ep_dt = Etp[i] * dt
    
	'''UPDATE'''
	# Interception Reservoir
	if P_dt > 0:
		Si[i] = Si[i] + P_dt*dt
		Pe_dt = np.maximum(0, (Si[i] - I_max) / dt)
		Si[i] = Si[i] - Pe_dt*dt
		Ei_dt[i] = 0
              
	else:
		Pe_dt = np.maximum(0, (Si[i] - I_max)/dt)
		Ei_dt[i] = np.minimum(Ep_dt, Si[i] / dt)
		Si[i] = Si[i] - Pe_dt*dt - Ei_dt[i]*dt      
	if i < t_max-1:
		Si[i+1] = Si[i]		  

	# Unsaturated Reservoir
	
	
	
	if Pe_dt > 0:
		Cr = (Su[i] / Su_max) ** beta
		Qiu_dt = (1 - Cr) * Pe_dt # flux from Ir to Ur
		Su[i] = Su[i] + Qiu_dt*dt
		Quf_dt = Cr * Pe_dt
	else:
		Quf_dt = 0

	# Transpiration
	Ep_dt = max(0, Ep_dt - Ei_dt[i])
	Ea_dt[i] = Ep_dt * (Su[i] / (Su_max * Ce))
	Ea_dt[i] = min(Su[i] / dt, Ea_dt[i])
	Su[i] = Su[i] - Ea_dt[i]*dt

	# Preferential Percolation
	Qus_dt = D * Quf_dt

	# Fast Reservoir
	Sf[i] = Sf[i] + (1-D) * Quf_dt
	Qf_dt[i] = dt * Kf * Sf[i]
	Sf[i] = Sf[i] - Qf_dt[i]*dt

	# Save output
	States[:,0] = Si
	States[:,1] = Su
	States[:,2] = Sf

	Fluxes[:,0] = Ei_dt
	Fluxes[:,1] = Ea_dt
	Fluxes[:,2] = Qf_dt
	Fluxes[:,3] = Qus_dt

	return(Fluxes, States)