import sys
sys.path.append('../../QC-DMET/src')
import localintegrals, dmet, qcdmet_paths
from pyscf import gto, scf, symm, future
from pyscf import mcscf
import numpy as np

localization = 'meta_lowdin'                # 'iao' or 'meta_lowdin' or 'boys'

bondlengths = np.arange(0.6, 2.75, 0.1)
energies = []

for bondlength in bondlengths:

	print ('-----------------------------------------------')
	print ('----  H10 at',bondlength,'angstroms seperation --------')
	print ('-----------------------------------------------')	
	
	nat = 10
	mol = gto.Mole()
	mol.atom = []
	r = 0.5 * bondlength / np.sin(np.pi/nat)
	for i in range(nat):
		theta = i * (2*np.pi/nat)
		mol.atom.append(('H', (r*np.cos(theta), r*np.sin(theta), 0)))

	mol.basis = '6-31g'
	mol.build(verbose=0)

	mf = scf.RHF(mol)
	mf.max_cycle = 1000
	mf.scf()
	mc = mcscf.CASSCF(mf,16,10)
	#cas_list = [12,14,15,16,17,18,19,20] 	
	#mo = mc.sort_mo(cas_list)
	the_energy = mc.kernel()[0]
	OccNum = mc.cas_natorb()[2]
	
	#PRINT DMET ORBITALS
	'''for mo in range(0,20):
		mo_coeff = MO[:,mo].reshape(mol.nao_nr())
		name = 'MO_cas44_' + str(mo) + '.cube'
		mocube.mo(mol, name, mo_coeff, nx=60, ny=60, nz=60)	

	#PRINT DMET NATURAL ORBITALS
	for mo in range(0,20):
		mo_coeff = MOnat[:,mo].reshape(mol.nao_nr())
		name = 'MO_nat_cas44_' + str(mo) + '.cube'
		mocube.mo(mol, name, mo_coeff, nx=60, ny=60, nz=60)'''		
	print ("----Energy at ", bondlength," angstrom:", the_energy)
	print ("Occupation number:",OccNum)	  

