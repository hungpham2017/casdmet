import sys
sys.path.append('../../QC-DMET/src')
import localintegrals, dmet, qcdmet_paths
from pyscf import gto, scf, symm, future
from pyscf import mcscf
import numpy as np
import C5N2H_struct
from pyscf.tools import mocube
from functools import reduce

#############
#   Input   #
#############              
localization = 'meta_lowdin'                # 'iao' or 'meta_lowdin' or 'boys'


casci_energy_formula = True

basis = '6-31g' 
DMguess = None
for distance in np.arange(1.2, 1.25, 0.1): 
	print ('-----------------------------------------------')
	print ('----  C5H11-N=N-H at',distance,'angstroms  --------')
	print ('-----------------------------------------------')	
	mol = C5N2H_struct.structure( distance, basis)
	mf = scf.RHF( mol )
	mf.verbose = 3
	mf.max_cycle = 500
	mf.scf(DMguess)
	DMguess = mf.make_rdm1()
    
	if ( True ):
		myInts = localintegrals.localintegrals( mf, range( mol.nao_nr() ), localization )
		unit_sizes = None
		if ( basis == '6-31g' ):
			unit_sizes = np.array([ 20, 13, 13, 13, 13, 15]) # HCC, 4(CH2), CH3 
		assert( np.sum( unit_sizes ) == mol.nao_nr() )

		impurityClusters = []
		if ( casci_energy_formula ): # Do only 1 impurity at the edge
			num_orb_in_imp = np.sum( unit_sizes[ 0 : 1 ] )
			impurity_orbitals = np.zeros( [ mol.nao_nr() ], dtype=int )
			impurity_orbitals[ 0 : num_orb_in_imp ] = 1
			impurityClusters.append( impurity_orbitals )
		else: # Partition
			jump = 0
			for fragment in range(unit_sizes.size):	
				impurity_orbitals = np.zeros( [ mol.nao_nr() ], dtype=int )	
				num_orb_in_imp = unit_sizes[fragment]
				if (fragment > 0): 
					impurity_orbitals[ jump : jump + num_orb_in_imp ] = -1
				else:
					impurity_orbitals[ jump : jump + num_orb_in_imp ] = 1
				impurityClusters.append( impurity_orbitals )
				jump += num_orb_in_imp
			
				
		isTranslationInvariant = False
		method = 'CASSCF'
		SCmethod = 'NONE'
		theDMET = dmet.dmet( myInts, impurityClusters, isTranslationInvariant, method, SCmethod, doDET=True )
		theDMET.impCAS = (4,4)
		theDMET.CASlist = [18,19,21,22]
		theDMET.CC_E_TYPE  = 'CASCI'				
		the_energy = theDMET.doselfconsistent()		
		
		'''X1 = np.linalg.inv(theDMET.ao2loc)  # since C~ = X.C (X is the inverse of transformation matrix)
		X2 = np.linalg.inv(theDMET.loc2dmet)
		#PRINT RHF in embedding space
		for mo in range(10,40):
		    mo_coeff = reduce(np.dot,(X1,X2[:,:40], theDMET.MOmf))[:,mo].reshape(mol.nao_nr())
		    name = 'MO_nat_casdmet44_' + str(mo) + '.cube'
		    id = mo + 1
		    mocube.mo(mol, name, id, mo_coeff, nx=60, ny=60, nz=60)'''			
		
		#PRINT LOCALIZED ORBITALS
		'''for mo in range(0,48):
		    mo_coeff = X1[:,mo].reshape(mol.nao_nr())
		    name = 'MO_loc_' + str(mo) + '.cube'
		    mocube.mo(mol, name, mo_coeff, nx=60, ny=60, nz=60)	

			
		#PRINT DMET ORBITALS
		for mo in range(0,48):
		    mo_coeff = np.dot(X1,X2)[:,mo].reshape(mol.nao_nr())
		    name = 'MO_dmet_' + str(mo) + '.cube'
		    mocube.mo(mol, name, mo_coeff, nx=60, ny=60, nz=60)

		#PRINT DMET ORBITALS
		for mo in range(0,34):
		    mo_coeff = reduce(np.dot,(X1,X2[:,:34], theDMET.MO))[:,mo].reshape(mol.nao_nr())
		    name = 'MO_casdmet44_' + str(mo) + '.cube'
		    mocube.mo(mol, name, mo_coeff, nx=60, ny=60, nz=60)	

		#PRINT DMET NATURAL ORBITALS
		for mo in range(0,34):
		    mo_coeff = reduce(np.dot,(X1,X2[:,:34], theDMET.MOnat))[:,mo].reshape(mol.nao_nr())
		    name = 'MO_nat_casdmet44_' + str(mo) + '.cube'
		    mocube.mo(mol, name, mo_coeff, nx=60, ny=60, nz=60)'''	
			
		print ("----Energy at ", distance," angstrom:", the_energy)
		print ("Occupation number:",theDMET.OccNum)
		

    
