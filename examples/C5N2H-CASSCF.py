import sys
sys.path.append('../../QC-DMET/src')
import localintegrals, dmet, qcdmet_paths
from pyscf import gto, scf, symm, future
from pyscf import mcscf
import numpy as np
import C5N2H_struct
from pyscf.tools import mocube

#############
#   Input   #
#############              
localization = 'meta_lowdin'                # 'iao' or 'meta_lowdin' or 'boys'


casci_energy_formula = False         # CASCI or DMET energy formula

basis = '6-31g' 
DMguess = None
for distance in np.arange(2.9, 2.95, 0.1): 
	print ('-----------------------------------------------')
	print ('----  C5H11-N=N-H at',distance,'angstroms  --------')
	print ('-----------------------------------------------')	
	mol = C5N2H_struct.structure( distance, basis)
	xyz = np.asarray(mol.atom)
#	for atom in xyz:
#		print(atom[0],atom[1],atom[2],atom[3])	
	mf = scf.RHF( mol )
	mf.verbose = 3
	mf.verbose = 3
	mf.max_cycle = 500
	mf.scf(DMguess)
	DMguess = mf.make_rdm1()
	MO = mf.mo_coeff
	#print(mf.mo_occ.sum()/2)
	mc = mcscf.CASSCF(mf,2,2)
	cas_list = [29,30]
	mo = mc.sort_mo(cas_list)
	the_energy = mc.kernel(mo)[0]
	#MOo = mo
	#MOnat = mc.cas_natorb()[0]
	
	print ("----Energy at ", distance," angstrom:", the_energy)
	
	#PRINT DMET NATURAL ORBITALS
	'''for mo in range(17,38):
		mo_coeff = MO[:,mo].reshape(mol.nao_nr())
		name = 'MO_' + str(mo) + '.cube'
		id = mo + 1
		mocube.mo(mol, name, id, mo_coeff, nx=60, ny=60, nz=60)'''
		

    
