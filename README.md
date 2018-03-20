# casdmet
The modified QC-DMET to use with CASSCF solver (CASDMET)

### OVERVIEW:
- Install PySCF, QC-DMET
- Copy casdmet.py amd pyscf_casscf.py to the QC-DMET home

### NOTE:
- The data in the paper (in the SI) was computed with PySCF-1.3. The authors noticed that there will be a small
difference if using PySCF-1.4 or newer. This may due to some changes in the mcscf module of PySCF.
