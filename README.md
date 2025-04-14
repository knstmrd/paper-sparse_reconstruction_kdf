# Reproducibility repository for "Sparse Reconstruction of Multi-Dimensional Kinetic Distributions"

This repository contains the code to reproduce the results in "Reproducibility repository for "Sparse Reconstruction of Multi-Dimensional Kinetic Distributions".

## Usage
Clone the repository; once cloned, navigate to the directory, run `julia --project=.`, and in the Julia interpreter run using `Pkg; Pkg.resolve(); Pkg.instantiate()` to install the required packages.

To run a simulation that reconstructs a 2-dimensional VDF, run `julia --project=. simulations/sparsity_promoting.jl`, this
will output results to the `out` directory.

The `sparsity_promoting.jl` file has two main routines: `main` and `constraint_only_main`. The second one only finds a solution
that satisfies the constraints without any further optimization.
`main` takes in as parameters
* `n_v_arr`: array with a list of number of Diracs to be used
* `λ1_arr`: array with a list of regularization strengths to be used
* `n_subs_arr`: number of subintervals in each dimension into which to subdivide the domain of the function to compute
    entropy using histograms
* `max_mom_c`: maximum conserved moment
* `max_next_mom`: maximum predicted moment
* `vdf_name`: the name of the VDF being reconstructed, used to name subdirectories for output
* `vdf_func`: a function which takes in a vector of moment indices and returns the computed moment (see examples in `src/basics.jl`)

Several CSV files are created: one storing the moments of the original (given) VDF and the ones obtained via the reconstruction,
one storing the computed solution points, and one storing additional information (number of distinct Diracs, pre-optimization
and post-optimization values of entropy).

## Citing
The pre-print can be cited as:
```bibtex
@online{oblapenko2024entropyconservative,
  title={Sparse Reconstruction of Multi-Dimensional Kinetic Distributions},
  author={Oblapenko, Georgii and Torrilhon, Manuel and Herty, Michael},
  year={2025},
  eprint={2504.06014},
  archivePrefix={arXiv},
  primaryClass={physics.comph-ph},
  eprintclass={physics.flu-dyn},
  journal={arXiv preprint 2504.06014},
  doi={10.48550/arXiv.2504.06014}
}
```

This repository can be cited as:
```bibtex
@misc{oblapenko2024entropyconservativeRepro,
  title={Reproducibility repository for "Sparse Reconstruction of Multi-Dimensional Kinetic Distributions"},
  author={Oblapenko, Georgii and Torrilhon, Manuel and Herty, Michael},
  year={2025},
  howpublished={\url{https://github.com/knstmrd/paper-sparse_reconstruction_kdf}},
  doi={10.5281/zenodo.15114064}
}
```