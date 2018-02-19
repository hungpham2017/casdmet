import sys
sys.path.append('../../QC-DMET/src')
import localintegrals, dmet, qcdmet_paths
from pyscf import gto, scf, symm, future
from pyscf import mcscf
import numpy as np

localization = 'meta_lowdin'                # 'iao' or 'meta_lowdin' or 'boys'

bondlengths = np.arange(2.7, 2.75, 0.1)
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
        
	if ( True ):
		myInts = localintegrals.localintegrals( mf, range( mol.nao_nr() ), localization )
		myInts.molden( 'hydrogen-loc.molden' )
		myInts.TI_OK = True # Only s functions

		atoms_per_imp = 2 # Impurity size = 1 atom
		assert ( nat % atoms_per_imp == 0 )
		orbs_per_imp = int(myInts.Norbs * atoms_per_imp / nat)

		impurityClusters = []
		for cluster in range( nat // atoms_per_imp ):
			impurities = np.zeros( [ myInts.Norbs ], dtype=int )
			for orb in range( orbs_per_imp ):
				impurities[ orbs_per_imp*cluster + orb ] = 1
			impurityClusters.append( impurities )
		isTranslationInvariant = True
		method = 'CASSCF'
		SCmethod = 'NONE'
		theDMET = dmet.dmet( myInts, impurityClusters, isTranslationInvariant, method, SCmethod, doDET=True )
		theDMET.impCAS = (6,6)
		theDMET.CASlist = [2,3,4,6,7,8]		
		the_energy = theDMET.doselfconsistent()
		energies.append(the_energy)		
		print ("----Energy at ", bondlength," angstrom:", the_energy)

print ("Bondlengths =", bondlengths)
print ("Energies =", energies)

