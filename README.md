# crowding-kMC

Simple fortran code to model crowding effect in Hsp70 chaperone binding to amyloid fibrils. 
The rejection-free kinetic Monte Carlo algorithm is used to calculate binding isotherms with a minimal discrete model taking into account 
NC chaperones (NC =[Hsc70]/ (0.01µM) in the original work) and fibrils represented as a linear array of NS binding sites (NS=[α-syn]/ (0.01µM) in the original work). 
Each free chaperone can bind to a free binding site with a rate kon and bound chaperones can dissociate with rate koff from fibril binding sites.
In this framework, the effects of crowding can be easily introduced by assuming that an occupied site disfavor further occupation of 
neighboring sites by increasing their unbinding rate and hence their dissociation constant. 
For the sake of simplicity, we assume here that this energetic penalty V is constant and has a range of n sites, i.e. it is applied to all other n proximal sites. 
The unbinding rate of a chaperone bound to site j is therefore expressed as koff(j)=koff_0*exp(q_j*V) where q_j is the number of occupied sites in a range of width n around site j.
The unperturbed rates kon and koff_0 are chosen in order to reproduce the experimental dissociation constants observed at low binder concentrations. 
The equilibrium distribution of the model was obtained by calculating time averages over simulation lengths significantly exceeding the microscopic binding/unbinding times (τ >> 1/ kb , 1/ ku).
