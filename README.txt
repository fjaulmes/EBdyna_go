Fabien Jaulmes, IPP Prague, 17th January 2018



Dear user and contributor,

This is the (GIT) repository for the EBdyna_go code.
This code is intended as an orbit-following code for tokamak full 3D toroidal configurations.
It aims at simulating trajectories of any confined particles (impurities, thermal plasma, fast particles) in a tokamak plasma.
It also offers you the possibility to represent dynamically evolving perturbations such as sawteeth or TAE.
The time scale that are reasonbale to simulate are typically collisionless ones, let us say roughly of order 1ms in JET.

This code was developped during my first years of contact with Nuclear Fusion research 
(during my PhD in The Netherlands, at the Dutch Institute For Fundamental Energy Research).
As such, it has a somewhat heavy "legacy" content that still needs to be cleaned and upgraded in terms of clarity and coding standards
(on this topic, you can read the "README" contribution of Selwyn below).

However, since I still find EBdyna_go is highly flexible and of high precision, scientific quality and interest, I would
like to continue its development. In doing so, I intend to use more and more the highly parallelized structure introduced 
by Selwyn during its Master of Science project (2016-2017 at DIFFER). Other significant improvements Selwyn has 
implemented include a more stable orbit-following algorithm, an analytical model for the sawtooth reconnection 
(developped by Hugo de Blank) and the Biot Savart (BS) solver for RMPs and TFR.

Unfortunately, these improvements came along a few other designs changes, such as shortening of variable names and heavy use
of structures and functions.
Whilst these are considered good coding practices, it removes a lot of the readability and "interactive rewriting of code"
that MATLAB is good at and that I enjoyed just for my own use of EBdyna.

But modernity shall prevail.
So, I am now undertaking the (long term, maybe 2 or 3 years) project of wrapping all the qualities of the new, powerful
EBdyna_go and make them more flexible and clearer in order to have them embrace all the possibilities 
(stabilitiy calculation, TAEs growth rate...) of the "PhD legacy" part of the code.

I create this repository with the hope that I will not be alone during the course of this work and that other researchers
or students can bring in significant contributions.

I have written this code with the constant goal in mind for it to be a tool to bring productive discussions and insights on plasma
phenomenons that are very complex to describe analyticaly and where the orbit topology matters in a crucial way.

I sincerely hope that this goal is not only my own but that of the international effort dedicated to the study of Nuclear
Fusion as a clean and lond-term energy source for mankind.

Long live this GIT repository and please contact me should you have any desire to work with me on this.
My current working email is : jaulmes@ipp.cas.cz


Fabien


 

------------------------------------------------------------------------------------
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
