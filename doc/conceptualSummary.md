\page summary_page Short summary of workflow

## Summary of Millepede workflow 

This page is intended to provide a quick, first-glance overview of how Millepede functions.
High-energy physics detector alignment is used as an example.

A more in-depth derivation of the mathematics and detailed functionality is available in the \subpage draftman_page. 

### Short summary of the principle

Millepede expresses a least-square fit problem as a (huge) system of linear equations, populated using data from a large number of individual measurements. The system is comprised of 
- *global* parameters which affect all of the measurements and are of interest to the user, 
- and *local* parameters which affect the cost function (e.g. a fit $\chi^2$) for small, self-contained sets of measurements, but are not of immediate interest to the user. 

In a typical HEP detector alignment scenario, this corresponds to 
- (thousands of) alignment parameters describing the deviations of the detector geometry from the reference assumption
- (billions of) track parameters describing the trajectories of the individual particles traversing the detector (tracks). 

In this example, the alignment parameters take the role of global parameters that the user wishes to determine, while the track parameters are local parameters.

"Out of the box", the dimensionality of such a system would rank in the billions, rendering it inaccessible to direct solution. 

Because the local parameters only affect individual measurements, they can typically be solved efficiently as a (large) set of low-dimensionality problems **if** the global parameters are assumed to be fixed. In HEP, this corresponds to performing individual track fits. 

Millepede exploits this to reduce the dimensionality of the system.
This is achieved by obtaining a *special* solution for the local parameters, for which the global parameters are treated as constant. 
This special solution is then used to reformulate and partition the overall system.
This leads to a sub-system of dimensionality of the number of global parameters that can be solved to determine the same. 

This solution is *exact* and fully accounts for the correlations between global and local parameters. The effect of the dimensional reduction is that the final values of the local parameters (which are not of interest to the user) are not explicitly determined. 

The reduced system is solved using a choice of several techniques from direct matrix inversion over diagonalisation to analyse the problem structure to fast, approximate approaches. 


### The workflow

\htmlonly <style>div.image img[src="fig_1.svg"]{width:600px;}</style> \endhtmlonly 
\image html fig_1.svg 

Millepede performs the work described above in two steps. 

First, the **input information** describing the problem is written into a custom binary format. This step is called *Mille*, and frequently implemented within the software frameworks of the experiments. 

This is combined with one or several **steering files** in plain text format, which allow the user to configure the program behaviour, define constraints on (combinations of) parameters, and configure starting values. 

This information enables the second step, *pede*, to run. This standalone program constructs the minimisation problem from the binary input files guided by the steering files, and then performs the reformulation and solution strategy described above. 

#### Using Mille to dump input files

The following information is required to construct the minimisation problem, for each measurement entering the system, 
- the residual of the given measurement. In HEP detector alignment, this is a local track fit residual - the difference between the predicted and measured location of an individual measurement on a track. 
- The uncertainty of the residual. In HEP alignment, this would be the hit position uncertainty.
- The derivatives of the residual by the global parameters affecting it. In HEP tracker alignment, this would be the impact of the alignment parameters on the locally predicted hit location. 
- The derivatives of the residual by the *local* parameters affecting it. In HEP tracker alignment, this would be the impact of the track parameters on the locally predicted hit location. 

Parameters are identified using **labels**, which are simply integer numbers. The labels do not have to be contiguous, and should be chosen in a way that allows them to be easily mapped back to their meaning in the problem at hand. 
In HEP, this can involve making the label of alignment parameters encode the particular detector component and the particular alignment degree of freedom it represents.   

A reference *Mille* implementation is provided in this package, which can store its output in custom C / ROOT / CSV format (pede can only read the first two). However any implementation producing the specified format (see Mille.f90) can be used. 

The [GeneralBrokenLines](https://gitlab.desy.de/claus.kleinwort/general-broken-lines) track fit is frequently used in HEP detector alignment to supply a track model well-suited to fast solution in the context of Millepede, and capable of correctly describing the impact of multiple scattering along the trajectory as required for an accurate alignment. It has built-in support for writing *Mille* binaries. 

#### Using pede to solve the problem 

After dumping the binary input, the user should prepare steering files to describe the desired solution strategy, implement appropriate constraints to avoid degeneracies in the problem, and fine-tune the solution to their liking. See \subpage option_page for a range of supported arguments. The \subpage draftman_page documents the steering in more detail. The steering file is also responsible for telling `pede` which input binaries to process. 

Then, `pede` can be invoked using the steering file as argument. 

Depending on the size of the problem, the run time and memory requirements can vary. A comparably small problem with a few hundred degrees of freedom can be solved in seconds, while a full-scale CMS detector alignment with 200.0000 alignment parameters takes around a day of run-time and on the order of 100GB of memory. 

#### Extracting the solution

After a successful run, `pede` writes the resulting global parameter values into a `millepede.res` file. For each parameter, the label and value, as well as additional information depending on the run mode (e.g. uncertainties, number of measurements affecting the parameter) are written out. This can be parsed by the client for further use.

Additional output files provide detail on the program execution, and can be used for troubleshooting and diagnostics. While `pede` as a fortran program does not return a classical return code, a `millepede.end` file is generated specifying an exit code (see \subpage exit_code_page) that can be used for a first-glance check of the outcome. 
