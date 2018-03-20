import sys
sys.path.append('../../QC-DMET/src')
import localintegrals, dmet, qcdmet_paths
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

for distance in np.arange(1.8, 1.85, 0.1): # Ni---H2O distance
	print ('-----------------------------------------------')
	print ('----  Me-N=N-Me at',distance,'angstroms  --------')
	print ('-----------------------------------------------')	
	mol = Me2N2_struct.structure( distance, basis)
	mf = scf.RHF( mol )
	mf.verbose = 3
	mf.scf()
    
	if ( True ):
		myInts = localintegrals.localintegrals( mf, range( mol.nao_nr() ), localization )
		myInts.molden( 'Me2N2.molden' )
		
		unit_sizes = None
		if ( basis == '6-31g' ):
			unit_sizes = np.array([ 18, 15, 15]) # N2, CH3, CH3
		assert( np.sum( unit_sizes ) == mol.nao_nr() )

		impurityClusters = []
		if ( casci_energy_formula ): # Do only 1 impurity at the edge
			num_orb_in_imp = np.sum( unit_sizes[ 0 : 1 ] )
			impurity_orbitals = np.zeros( [ mol.nao_nr() ], dtype=int )
			impurity_orbitals[ 0 : num_orb_in_imp ] = 1
			impurityClusters.append( impurity_orbitals )
		else: # Partition
			jump = 0
			for fragment in range(3):	
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
		theDMET.impCAS = (8,8)
		theDMET.CASlist = [12,13,14,15,17,18,19,20]	
		the_energy = theDMET.doselfconsistent()		
		
		X1 = np.linalg.inv(theDMET.ao2loc)  # since C~ = X.C (X is the inverse of transformation matrix)
		#PRINT LOCALIZED ORBITALS
		'''for mo in range(0,48):
		    mo_coeff = X1[:,mo].reshape(mol.nao_nr())
		    name = 'MO_loc_' + str(mo) + '.cube'
		    mocube.mo(mol, name, mo_coeff, nx=60, ny=60, nz=60)	

		X2 = np.linalg.inv(theDMET.loc2dmet)			
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
		

    
