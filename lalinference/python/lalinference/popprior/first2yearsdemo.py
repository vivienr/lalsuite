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

import numpy as np
import h5py
import sqlite3
import lalinference.popprior.uberbank_database as ud
from glue.text_progress_bar import ProgressBar
from lalinference import popprior

start_time = time.time()

x = ud.Bank(sqlite3.connect("uberbank_database.sqlite")
plots=False
datetime=time.strftime('%Y%m%d')
bank = x.get_templates()
save_data = True

f, popt = popprior.source_population("Kiziltan2012_DNS.txt") # find probability density of masses for BNS source population
mass = load_bank(bank) # load mass
#mass = load_bank(bank) # load mass and chi
overlap = load_overlaps(bank)
print "Overlap data loaded. Total time elapsed:", time.time()-start_time

############################################
############################################

# CALCULATE PROBABILITIES

n = 10
k = np.arange(1,n+1)
rho = 30*np.cos((2*k-1)*np.pi/(2*n)) + 30 # min(rho)~0, max(rho)~80
rho.sort()
rho = np.insert(rho,0,0)

ln_p_j, p_j_templates, p_j_indices = find_ln_p_j_voronoi(mass, f, popt) #p_j_indices = indices of the templates for m1, m2 arrays (eg. p_j_indices[0] = 1, which means p_j_template[0] corresponds to the 1th template in m1, m2)
t_k = np.array(range(len(mass)))
ln_p_jk = np.log(np.zeros((len(rho), len(t_k)))) # P(signal t_j is recovered by template t_k)

#p_j_indices.tolist().sort(key=ln_p_j.__getitem__) # doing loop in order of p_j (smallest to largest) for numerical accuracy
order = np.argsort(ln_p_j)
ln_p_j, p_j_templates, p_j_indices = ln_p_j[order], p_j_templates[order], p_j_indices[order]

progress = ProgressBar(max=len(p_j_indices))
for i in range(len(p_j_indices)): # loop over all signal population
        progress.increment(text=str(p_j_templates[order][i]))
        ovrlp = overlap[p_j_indices[i]]
        for r in range(len(rho)): # loop over all rho
                ln_p_jk[r,:] = np.logaddexp(ln_p_jk[r,:], ln_p_j[order][i]+ln_p_k(ovrlp, rho[r], t_k))
        
print "ln(P_jk) computed for all templates. Time elapsed:", time.time()-start_time
                
# Save data to hdf5 file
if save_data:
        directory = ""
        filename = "logP_vs_rho_"+datetime
        f = h5py.File(directory+filename+".hdf5","w")
        f.create_dataset('rho', data=rho)
        f.create_dataset('ln_P_jk', data=ln_p_jk)
        f.create_dataset('masses', data=mass)
        f.close()
        print "Data saved. Time elapsed:", time.time()-start_time
