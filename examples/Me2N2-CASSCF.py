import sys
from pyscf import gto, scf, symm, future
from pyscf import mcscf
import numpy as np
import Me2N2_struct

#############
#   Input   #
#############              
localization = 'meta_lowdin'                # 'iao' or 'meta_lowdin' or 'boys'


one_bath_orb_per_bond = False        # Sun & Chan, JCTC 10, 3784 (2014) [ http://dx.doi.org/10.1021/ct500512f ]
casci_energy_formula = False         # CASCI or DMET energy formula

basis = '6-31g' 

for distance in np.arange(1.2, 1.25, 0.1): # Ni---H2O distance
	print ('-----------------------------------------------')
	print ('----  Me-N=N-Me at',distance,'angstroms  --------')
	print ('-----------------------------------------------')	
	mol = Me2N2_struct.structure( distance, basis)
	mf = scf.RHF( mol )
	mf.verbose = 3
	mf.scf()
	mc = mcscf.CASSCF(mf,8,8)
	cas_list = [12,14,15,16,17,18,19,20] 	
	mo = mc.sort_mo(cas_list)
	the_energy = mc.kernel(mo)[0]
	print ("----Energy at ", distance," angstrom:", the_energy)
		

    
