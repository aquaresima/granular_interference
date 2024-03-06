# Granular Interference

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> GrInt.jl

The code is meant to reproduce the results of the journal article:

___Shear profile in a dense packing of large grains___

Alessio Quaresima, Andrea Plati, Andrea Gnoli, Alberto Petri

The manuscript is under review and it is accessible on arXiv at the [link](https://arxiv.org/abs/2401.10062)


## Reproducibility of the analysis (DrWatson)
To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```
   This will install all necessary packages for you to be able to run the scripts and everything should work out of the box, including correctly finding local paths.

2. Download the processed data from : [analysis_data](https://owncloud.gwdg.de/index.php/s/ionezChYGMKuASM) and placed them in the folder 'data'.
3. Verify the path of the data corresponds to the one set for "analysis_path" in 'conf.yaml'

4. Run files '4_paper_figures.jl' and '5_supplementary_figures.jl'


# Raw data

Raw data are available upon request and the packages provides the script to process them. 
If you want to process your own data with the package set the 'full_matrices' path in the çonfl.yaml' file and run scripts: 
   '1_compute_variance.jl' 
   '2_compute_correlations.jl 
   '3_fit_correlations.jl'
   
 

