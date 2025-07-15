import camb
import numpy as np
import matplotlib.pyplot as plt

print(camb.__file__)

outdir = '/home/simran/particle_burst_dr6'

plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18

A = 0.01
omega = 10
phi = 0

ks = np.logspace(-4.5, 0, 500)

def Pk_oscillations(A_s, n_s, k_0, A, omega, phi, k, spacing):
    #eq 27 - 30 from 2309.17287
    orig = A_s*(k/k_0)**(n_s - 1)

    if spacing == 'linear':

        return orig*(1 + A*np.sin(omega*(k/k_0) + 2.*np.pi*phi))

    if spacing == 'log':

        return orig*(1 + A*np.sin(omega*(np.log(k/k_0)) + 2.*np.pi*phi))


'''
Linear
'''


cosmo_params_camb_lin = camb.CAMBparams()
cosmo_params_camb_lin.set_cosmology(H0 = 67.16, tau = 0.05357, ombh2 = 0.02237*0.6716**2, omch2 = 0.1206*0.6716**2)

cosmo_params_camb_lin.InitPower.set_params(r=0, A_osc = A, omega_osc = omega, phi_osc = phi, spacing_osc = 1,
                                            As = np.exp(3.044)/(1e10), ns = 0.965)
scalar_pk_camb_lin = cosmo_params_camb_lin.scalar_power(ks)

'''
Log
'''


cosmo_params_camb_log = camb.CAMBparams()
cosmo_params_camb_log.set_cosmology(H0 = 67.16, tau = 0.05357, ombh2 = 0.02237*0.6716**2, omch2 = 0.1206*0.6716**2)

cosmo_params_camb_log.InitPower.set_params(r=0, A_osc = A, omega_osc = omega, phi_osc = phi, spacing_osc = 2,
                                            As = np.exp(3.044)/(1e10), ns = 0.965)
scalar_pk_camb_log = cosmo_params_camb_log.scalar_power(ks)



plt.figure(figsize = (10, 8))

plt.plot(ks, Pk_oscillations(np.exp(3.044)/(1e10), 0.965, 0.05, A, omega, phi, ks, 'linear'), color = 'coral', label = 'Python')
plt.plot(ks, scalar_pk_camb_lin, color = 'purple', ls = '--', label = 'CAMB')

plt.xlabel("$k$ [Mpc$^{-1}]$", fontsize = 22)
plt.ylabel("$P_{\\mathcal{R}}(k)$", fontsize = 22)
plt.loglog()
plt.legend(loc = 'best', fontsize = 18)
plt.tight_layout()
plt.savefig('{}/linear_comparing_CAMB_to_pythonfunc.png'.format(outdir))

plt.figure(figsize = (10, 8))

plt.plot(ks, Pk_oscillations(np.exp(3.044)/(1e10), 0.965, 0.05, A, omega, phi, ks, 'log'), color = 'coral', label = 'Python')
plt.plot(ks, scalar_pk_camb_log, color = 'purple', ls = '--', label = 'CAMB')

plt.xlabel("$k$ [Mpc$^{-1}]$", fontsize = 22)
plt.ylabel("$P_{\\mathcal{R}}(k)$", fontsize = 22)
plt.loglog()
plt.legend(loc = 'best', fontsize = 18)
plt.tight_layout()
plt.savefig('{}/log_comparing_CAMB_to_pythonfunc.png'.format(outdir))