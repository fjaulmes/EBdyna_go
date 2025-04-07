Hi, here Selwyn.

Nothing much has been documented about the FINESSE equilibrium etc.
Here you find files that have been used to solve the GS-equation.
The best documentation of the code is found in the .pdf (my report), appendix C.
Note that I've not modelled TAE's, so you'll need to look at the thesis of Fabien Jaulmes for this subject:
F. Jaulmes, \Kinetic behaviour of ions in tokamak inductive scenarios," Ph.D. dissertation,
Eindhoven University of Technology, 2016. [Online]: http://repository.tue.nl/841352

Note, the NBI-distributions have been stored in:
\\RIJNH\Departments\Fusiefysica\CPP\CODES\EBdyna_go\NBI-distributions

The NBI-distribution that I've been using was called: initial_NBI_60keV_transp_D_distribution_all.
I believe Fabien has used this one as well (this is the one the once recommended me).
Here below the info about the other NBI-distributions in the email I've previously sent to Egbert Westerhof.

List of other distributions:
	MB_D
	Deuterium
	0.8×〖10〗^6 deeltjes
	[200 – 30000] eV               (thermal + small 30 keV bump)
	Allen binnen ψ ̅=0.6 (0 op magnetische as, 1 op LCFS)
	With rotation profile
	With `current’ for each particle?
	MB_W
	Wolfraam
	1×〖10〗^6 deeltjes
	[200 – 30000] eV (thermal + small 30 keV bump)
	Allen binnen ψ ̅=0.7 (0 op magnetische as, 1 op LCFS)
	With rotation profile
	With `current’ for each particle?
	Initial_flat_D_40keV
	Deuterium
	1.2×〖10〗^6 deeltjes
	[39.75 40.25] keV                             flat
	Full tokamak
	Initial_flat_D40keV_pre
	Same as Initial_flat_D_40keV, but without those on prompt loss orbits
	MB_D_f
	Deuterium
	2.47×〖10〗^6 deeltjes
	[0 20] keV,          flat
	Allen binnen ψ ̅=0.7 (0 op magnetische as, 1 op LCFS)
	With `current’ for each particle?
	MB_D_pre
	Deuterium
	7.4×〖10〗^6 deeltjes
	[0 20] keV,          flat
	Allen binnen ψ≈0.9 (0 op magnetische as, 1 op LCFS)
	NBI60keV_transp_D
	Deuterium
	2×〖10〗^6 deeltjes
	[0 60] keV,          NBI distribution
	Allen binnen ψ ̅≈0.9 (0 op magnetische as, 1 op LCFS)
	NB: Lijken geclusterd op theta=constant lijnen te liggen
	NBI60keV_pre
	Deuterium
	1.99×〖10〗^6 deeltjes
	[0 60] keV,          NBI distribution
	Allen binnen ψ ̅≈0.9 (0 op magnetische as, 1 op LCFS)
	NBI60keV_R_pre
	Deuterium
	1.99×〖10〗^6 deeltjes
	[0 60] keV,          NBI distribution
	Allen binnen ψ ̅≈0.9 (0 op magnetische as, 1 op LCFS)
	NBI60keV_transpM_D
	Deuterium
	2×〖10〗^6 deeltjes (exactly the same as NBI60keV_transp_D (1999964) )
	[0 60] keV,          NBI distribution
	Allen binnen ψ ̅≈0.9 (0 op magnetische as, 1 op LCFS)
	NBI60keV_transpR_D
	Deuterium
	1×〖10〗^6 deeltjes
	[0 60] keV,          NBI distribution
	Allen binnen ψ ̅≈0.9 (0 op magnetische as, 1 op LCFS)
	W35_pre
	Wolfraam
	1×〖10〗^6 deeltjes
	[1 15] keV            rapidly decreasing 
	Allen binnen ψ ̅≈0.6 (0 op magnetische as, 1 op LCFS)

The distribution ` NBI60keV_transpM_D’ contains the same particles (energy and position) as  `NBI60keV_transp_D’.
Every *_R and *_pre-file appears to be an output of a simulation.

(NB: Names given by Fabien are not always clear to me.)
