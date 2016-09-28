# Copyright (C) 2016  Heather Fong, Kipp Cannon
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

# copied from calculations_updated.py, restructured to input different rho for select t_k
# sample mass distribution directly from template bank
# new P_tj calculation: using Voronoi diagram
# use Chebyshev nodes to compute values of rho
# Edited with Kipp, July 25

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import h5py
from scipy.misc import logsumexp
from time import time
import random
from scipy.optimize import curve_fit
from scipy.interpolate import lagrange, interp1d
from scipy.integrate import quad
from glue.text_progress_bar import ProgressBar
import sys

def open_hdf5(filename):
    f = h5py.File(filename, 'r')
    return f['rho'].value, f['ln_P_jk'].value, f['templates'].value

def find_rho_interp(start, end, rho, num=100):
    new_rho = np.linspace(start, end, num)
    new_rho = np.append(new_rho, rho)
    new_rho.sort()
    return new_rho[np.where(new_rho==start)[0][0]:]

def P_rho(log_const=0):
    return lambda x: np.exp(log_const + np.log(192./x**4))

def set_const(lnP, rho, rho_interp, logp_jk):
    idx = min(np.where(lnP==lnP[-1])[0])
    rho_idx = rho[idx]
    idx_interp = np.where(rho_interp==rho_idx)[0]
    logP_jk = np.full(len(rho_interp), lnP[-1])
    logP_jk[:idx_interp] = logp_jk(rho_interp)[:idx_interp]
    return logP_jk
    
start_time = time()

pop = sys.argv[1]

x = np.array([10,1100,5245, 6150,10100,50100,102000])

low_SNR, high_SNR = 4, 60
#rho, lnP, templates = open_hdf5("logP_vs_rho_20160917_"+pop+".hdf5")
rho, lnP, templates = open_hdf5("logP_vs_rho_"+pop+"_160923.hdf5")
rho_interp = find_rho_interp(low_SNR, max(rho), rho)
P_rho1 = P_rho()

fits = np.zeros(len(lnP[0]))

lagrangepoly = []

#try:
#    np.loadtxt("P_signal_recovered_by_tj_numericalintonly.txt",unpack=True)
#    ans = raw_input("File already exists, do you want to overwrite? ")
#    if ans != "yes":
#        save_data = False
#    if ans == "yes":
#        f = open("P_signal_recovered_by_tj.txt","w")
#        save_data = True

f = open("P_signal_recovered_by_tj_morerho_"+pop+"50rho.txt","w")
save_data = False
    
progress = ProgressBar(max=len(lnP[0]))
for k in range(len(lnP[0])):
    #for k in x:
    progress.increment(text=str(k)+"/"+str(len(lnP[0])))
    #polyfits = np.array(lagrange(rho, lnP[:,k]))
    #logp_jk = lagrange(rho, lnP[:,k])
    logp_jk = interp1d(rho, lnP[:,k], kind='cubic')
    logP_jk = logp_jk(rho_interp)
    if k in x:
        plt.plot(rho_interp, np.exp(logP_jk), '-')
        plt.plot(rho, np.exp(lnP[:,k]), 'o')
    #logP_jk = set_const(lnP[:,k], rho, rho_interp, logp_jk)
    #lagrangepoly.append(polyfits)
    P_rho2 = P_rho(log_const=logP_jk[-1])
    fits[k] = np.trapz(np.exp(logP_jk)*P_rho1(rho_interp), rho_interp) + quad(P_rho2, high_SNR, np.inf)[0]
    #fits[:,k] = polyfits(rho_interp)
    #print templates[k]
    if save_data:
        f.write(str(fits[k])+"\n")

if save_data:
    f.close()

print "Elapsed time:", time()-start_time

if sum(fits) > 1.05 or sum(fits) < 0.95:
    print "fits not normalized:", sum(fits)

#plt.yscale('log')
plt.legend(loc="lower left")
plt.xlabel(r'$\rho$')
plt.ylabel(r'$\text{log}P$')
plt.title(pop)
plt.savefig('lnP_vs_rho_'+pop+'.pdf')
    
#lagrangepoly = np.array(lagrangepoly)

f2 = h5py.File("Pjk_"+pop+".hdf5","w")
f2.create_dataset('P_jk',data=fits)
f2.create_dataset('templates',data=templates)
f2.close()
    
#for i in x:
#    plt.plot(rho_interp, fits[:,i],label=templates[i])
#    plt.plot(rho, lnP[:,i], 'o')

#plt.legend(loc="lower left", fontsize=6)
#plt.grid()
#plt.xlabel(r'$\rho$')
#plt.ylabel(r'$P_{jk}$, probability that signal is recovered by the given template')
#plt.savefig("lagrangian_fit_"+pop+".pdf")

#print "Elapsed time:", time()-start_time
