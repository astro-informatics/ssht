# SSHT Python Documentation

This guide is intended to explanin the python interface of SSHT. For a description of the workings of SSHT see [here](http://astro-informatics.github.io/ssht/ "SSHT documentation")

## ssht_forward

~~~~{.python}
flm = ssht_forward(f, int L, Spin=0, Method='MW', Reality=False)
~~~~

Performs the forward spherical harmonic transform.

###Inputs

* `f` the signal on the sphere, `numpy.ndarray` type `complex` or `real`. NB different for `'MW_pole'` sampling.
* `L` the band limit of the signal, non-zero positive integer
* `Spin` the spin of the signal, non-negative integer (default = 0)
* `Method` the sampling scheme used, string:
    1. `'MW'`         [McEwen & Wiaux sampling (default)]
    2. `'MW_pole'`    [McEwen & Wiaux sampling with the south pole as a seperate double.]
    3. `'MWSS'`       [McEwen & Wiaux symmetric sampling]
    4. `'DH'`         [Driscoll & Healy sampling]
    5. `'GL'`         [Gauss-Legendre sampling]
* `Reality`  determines if the signal is real or complex, Boolian (default = False)

### Output

`flm` the spherical harmonic transform of `f`, 1D `numpy.ndarray` type `complex`

#### Note on `'MW_pole'` sampling

This is the same as the `'MW'` sampling however the south pole is expessed as a double only not a vector. Therefore the size of the array is one smaller on the \(\theta\) direction. `f` is now a tuple containing either:
* `(f_array, f_sp, phi_sp)` if complex transform
* `(f_array, f_sp)` if real transform


## ssht_inverse

~~~~
f = ssht_inverse(flm, L, Spin=0, Method='MW', Reality=False)
~~~~

Performs the inverse spherical harmonic transform.

### Inputs

* `flm` the spherical harmonic transform of `f`, `numpy.ndarray` type `complex`
* `L` the band limit of the signal, non-zero positive integer
* `Spin` the spin of the signal, non-negative integer (default = 0)
* `Method` the sampling scheme used, string:
    1. `'MW'`         [McEwen & Wiaux sampling (default)]
    2. `'MW_pole'`    [McEwen & Wiaux sampling with the south pole as a seperate double.]
    3. `'MWSS'`       [McEwen & Wiaux symmetric sampling]
    4. `'DH'`         [Driscoll & Healy sampling]
    5. `'GL'`         [Gauss-Legendre sampling]
* `Reality`  determines if the signal is real or complex, Boolian (default = False)

### Output

`f` the signal on the sphere, 2D `numpy.ndarray` type `complex` or `real`. NB different for `'MW_pole'` sampling.

#### Note on `'MW_pole'` sampling

This is the same as the `'MW'` sampling however the south pole is expessed as a double only not a vector. Therefore the size of the array is one smaller on the \(\theta\) direction. `f` is now a tuple containing either:
* `(f_array, f_sp, phi_sp)` if complex transform
* `(f_array, f_sp)` if real transform

## ssht_ind2elm

~~~~
(n_theta, n_phi) = ssht_sample_shape(int L, Method='MW')
~~~~

Computes the index in the `flm` array of a particular harmonic coeeficient \(\ell \) and \(m\).

### Inputs

* `L` the band limit of the signal, non-zero positive integer
* `Method` the sampling scheme used, string:
    1. `'MW'`         [McEwen & Wiaux sampling (default)]
    2. `'MW_pole'`    [McEwen & Wiaux sampling with the south pole as a seperate double.]
    3. `'MWSS'`       [McEwen & Wiaux symmetric sampling]
    4. `'DH'`         [Driscoll & Healy sampling]
    5. `'GL'`         [Gauss-Legendre sampling]

### Output

Tuple containing `(n_theta, n_phi)`
* `n_theta` the number of samples in the \(\theta\) direction, integer
* `n_phi` the number of samples in the \(\phi\) direction, integer

### Output

* Tuple containing `n_theta, n_phi`, integers describin the number of samples in the \(\theta\) and \(\phi\) directions.

## ssht_sample_positions

~~~~
(thetas, phis) = ssht_sample_positions(int L, Method = 'MW', Grid=False)
~~~~

Computes the positions on the sphere of the samples.

### Inputs

* `L` the band limit of the signal, non-zero positive integer
* `Method` the sampling scheme used, string:
    1. `'MW'`         [McEwen & Wiaux sampling (default)]
    2. `'MW_pole'`    [McEwen & Wiaux sampling with the south pole as a seperate double.]
    3. `'MWSS'`       [McEwen & Wiaux symmetric sampling]
    4. `'DH'`         [Driscoll & Healy sampling]
    5. `'GL'`         [Gauss-Legendre sampling]
* `Grid` describes if the output is a vector of the smaple positions or a 2D array the same shape as the signal on the sphere.

### Outputs

Tuple containing `(thetas, phis)`
* `thetas` positions of the samples in the \(\theta\) direction
* `phis` positions of the samples in the \(\theta\) direction
