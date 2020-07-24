import sympy as sp
import numpy as np

## Define Symbols for Input Paraeters
t , t_0 , Om_0 , Omdot_0 , Om_QNM , h_0 , tau = sp.symbols("t t_0 Om_0 Om_0dot Om_Q h_0 tau")

## Computed Quantities:

## Matching Quantities
tanh_arg = -sp.Rational(1,2)*sp.log( (Om_QNM**4 - Om_0**4)/(2*tau*(Om_0**3)*Omdot_0) - 1 )
tanh_term = sp.tanh(tanh_arg)

##Frequency Calculation
Om4p = (Om_0**4 - (Om_QNM**4)*tanh_term)/(1 - tanh_term)
Om4m = (Om_QNM**4 - Om_0**4)/(1 - tanh_term)
kappa_p = (Om_0**4 + Om4m*(1 - tanh_term))**(sp.Rational(1,4))
kappa_m = (Om_0**4 - Om4m*(1 + tanh_term))**(sp.Rational(1,4))
t_p = t_0 - tau*tanh_arg
Omega = ( Om4p + Om4m*sp.tanh( (t - t_p) / tau ) )**sp.Rational(1,4)

## Phase Calculation
## used here, the alternative definition of arctanh
## arctanh(x) = 0.5*ln( (1+x)/(1-x) )

arctanhp = sp.Rational(1,2)*kappa_p*tau*sp.log( (1 + Omega/kappa_p)*(1 - Om_0/kappa_p) / ( (1 - Omega/kappa_p)*(1 + Om_0/kappa_p) ) )
arctanhm = sp.Rational(1,2)*kappa_m*tau*sp.log( (1 + Omega/kappa_m)*(1 - Om_0/kappa_m) / ( (1 - Omega/kappa_m)*(1 + Om_0/kappa_m) ) )
arctanp = kappa_p*tau*( sp.atan(Omega/kappa_p) - sp.atan(Om_0/kappa_p) )
arctanm = kappa_m*tau*( sp.atan(Omega/kappa_m) - sp.atan(Om_0/kappa_m) )

##Amplitude Calculation
A_0 = 4*h_0*(Om_0**2)
A_p = A_0*sp.cosh(tanh_arg)
psi4_amp = A_p/sp.cosh( (t - t_p)/tau )
strain_amp = psi4_amp/(4*(Omega**2))

## Define frequency function
frequency = sp.lambdify( [t,t_0,Om_0,Omdot_0,Om_QNM,tau], Omega )

## Define phase function
phase = sp.lambdify([t,t_0,Om_0,Omdot_0,Om_QNM,tau], arctanhp+arctanp-arctanhm-arctanm)

## Define amplitude function
strain_amplitude = sp.lambdify([t,t_0,Om_0,Omdot_0,Om_QNM,h_0,tau], strain_amp)