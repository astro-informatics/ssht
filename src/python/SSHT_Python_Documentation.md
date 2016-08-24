# SSHT Python Documentation

This guide is intended to explanin the python interface of SSHT. For a description of the workings of SSHT see [here](http://astro-informatics.github.io/ssht/ "SSHT documentation")

## pyssht.forward

~~~~{.python}
flm = pyssht.forward(f, int L, Spin=0, Method='MW', Reality=False)
~~~~

Performs the forward spherical harmonic transform.

#### Inputs

* `f` the signal on the sphere, `numpy.ndarray` type `complex` or `real`, ndim 2. NB different for `'MW_pole'` sampling.
* `L` the band limit of the signal, non-zero positive integer
* `Spin` the spin of the signal, non-negative integer (default = 0)
* `Method` the sampling scheme used, string:
    1. `'MW'`         [McEwen & Wiaux sampling (default)]
    2. `'MW_pole'`    [McEwen & Wiaux sampling with the south pole as a seperate double.]
    3. `'MWSS'`       [McEwen & Wiaux symmetric sampling]
    4. `'DH'`         [Driscoll & Healy sampling]
    5. `'GL'`         [Gauss-Legendre sampling]
* `Reality`  determines if the signal is real or complex, Boolian (default = False)

#### Output

`flm` the spherical harmonic transform of `f`, 1D `numpy.ndarray` type `complex`

#### Note on `'MW_pole'` sampling

This is the same as the `'MW'` sampling however the south pole is expessed as a double only not a vector. Therefore the size of the array is one smaller on the \(\theta\) direction. `f` is now a tuple containing either:
* `(f_array, f_sp, phi_sp)` if complex transform
* `(f_array, f_sp)` if real transform



## pyssht.inverse

~~~~
f = pyssht.inverse(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, L, Spin=0, Method='MW', Reality=False)
~~~~

Performs the inverse spherical harmonic transform.

#### Inputs

* `flm` the spherical harmonic transform of `f`, `numpy.ndarray` type `complex`, ndim 1
* `L` the band limit of the signal, non-zero positive integer
* `Spin` the spin of the signal, non-negative integer (default = 0)
* `Method` the sampling scheme used, string:
    1. `'MW'`         [McEwen & Wiaux sampling (default)]
    2. `'MW_pole'`    [McEwen & Wiaux sampling with the south pole as a seperate double.]
    3. `'MWSS'`       [McEwen & Wiaux symmetric sampling]
    4. `'DH'`         [Driscoll & Healy sampling]
    5. `'GL'`         [Gauss-Legendre sampling]
* `Reality`  determines if the signal is real or complex, Boolian (default = False)

#### Output

`f` the signal on the sphere, 2D `numpy.ndarray` type `complex` or `real`. NB different for `'MW_pole'` sampling.

#### Note on `'MW_pole'` sampling

This is the same as the `'MW'` sampling however the south pole is expessed as a double only not a vector. Therefore the size of the array is one smaller on the \(\theta\) direction. `f` is now a tuple containing either:
* `(f_array, f_sp, phi_sp)` if complex transform
* `(f_array, f_sp)` if real transform



## pyssht.forward_adjoint

~~~~
f = pyssht.forward_adjoint(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, L, Spin=0, Method='MW', Reality=False)
~~~~

Performs the adjoint of the forward spherical harmonic transform.

#### Inputs

* `flm` the spherical harmonic transform of `f`, `numpy.ndarray` type `complex`, ndim 1
* `L` the band limit of the signal, non-zero positive integer
* `Spin` the spin of the signal, non-negative integer (default = 0)
* `Method` the sampling scheme used, string:
    1. `'MW'`         [McEwen & Wiaux sampling (default)]
    3. `'MWSS'`       [McEwen & Wiaux symmetric sampling]
* `Reality`  determines if the signal is real or complex, Boolian (default = False)

#### Output

`f` the signal on the sphere, 2D `numpy.ndarray` type `complex` or `real`.


## pyssht.inverse_adjoint

~~~~{.python}
flm = pyssht.inverse_adjoint(f, int L, Spin=0, Method='MW', Reality=False)
~~~~

Performs the adjoint of the inverse spherical harmonic transform.

#### Inputs

* `f` the signal on the sphere, `numpy.ndarray` type `complex` or `real`, ndim 2.
* `L` the band limit of the signal, non-zero positive integer
* `Spin` the spin of the signal, non-negative integer (default = 0)
* `Method` the sampling scheme used, string:
    1. `'MW'`         [McEwen & Wiaux sampling (default)]
    3. `'MWSS'`       [McEwen & Wiaux symmetric sampling]
* `Reality`  determines if the signal is real or complex, Boolian (default = False)

#### Output

`flm` the spherical harmonic transform of `f`, 1D `numpy.ndarray` type `complex`


## pyssht.elm2ind

~~~~
index = pyssht.elm2ind( int el, int m)
~~~~

Computes the index in the `flm` array of a particular harmonic coeeficient \(\ell \) and \(m\).

#### Inputs

* `el` the scale parameter of the spherical harmonic coefficients, integer from \(0\) to \(L-1\), where \(L\) is the band limit.
* `em` the azimuthal parameter, integer from -el to el.

#### Output

Index of the coficiant in `flm` array, integer

## pyssht.ind2elm

~~~~
(el, em) = pyssht.ind2elm(int ind)
~~~~

Computes harmonic coeeficient \(\ell \) and \(m\) from the index in the `flm` array.

#### Inputs

* `ind` index of the `flm` array

#### Output

Tuple containing `(el, em)`
* `el` the scale parameter of the spherical harmonic coefficients, integer from \(0\) to \(L-1\), where \(L\) is the band limit.
* `em` the azimuthal parameter, integer from -el to el.


## pyssht.sample_shape

~~~~
(n_theta, n_phi) = pyssht.sample_shape(int L, Method='MW')
~~~~

Outputs a tuple with the shape of the array used for storing the data on the sphere for different sampling schemes.

#### Inputs

* `L` the band limit of the signal, non-zero positive integer
* `Method` the sampling scheme used, string:
    1. `'MW'`         [McEwen & Wiaux sampling (default)]
    2. `'MW_pole'`    [McEwen & Wiaux sampling with the south pole as a seperate double.]
    3. `'MWSS'`       [McEwen & Wiaux symmetric sampling]
    4. `'DH'`         [Driscoll & Healy sampling]
    5. `'GL'`         [Gauss-Legendre sampling]

#### Output

Tuple containing `(n_theta, n_phi)`
* `n_theta` the number of samples in the \(\theta\) direction, integer
* `n_phi` the number of samples in the \(\phi\) direction, integer


## pyssht.sample_positions

~~~~
(thetas, phis) = pyssht.sample_positions(int L, Method = 'MW', Grid=False)
~~~~

Computes the positions on the sphere of the samples.

#### Inputs

* `L` the band limit of the signal, non-zero positive integer
* `Method` the sampling scheme used, string:
    1. `'MW'`         [McEwen & Wiaux sampling (default)]
    2. `'MW_pole'`    [McEwen & Wiaux sampling with the south pole as a seperate double.]
    3. `'MWSS'`       [McEwen & Wiaux symmetric sampling]
    4. `'DH'`         [Driscoll & Healy sampling]
    5. `'GL'`         [Gauss-Legendre sampling]
* `Grid` describes if the output is a vector of the smaple positions or a 2D array the same shape as the signal on the sphere, default `False`

#### Outputs

Tuple containing `(thetas, phis)`
* `thetas` positions of the samples in the \(\theta\) direction
* `phis` positions of the samples in the \(\theta\) direction

## pyssht.s2_to_cart

~~~~
(x, y, z) = pyssht.s2_to_cart(theta, phi)
~~~~

Computes the \(x\), \(y\), and \(z\) coordinates from \(\theta\) and \(\phi\) on the sphere.

#### Inputs

* `theta` \(\theta\) values, type `numpy.ndarray`
* `phi` \(\phi\) values, type `numpy.ndarray`

#### Output

Tuple containing `(x, y, z)`
* `x` the \(x\) coordinate of each point, type `numpy.ndarray`
* `y` the \(y\) coordinate of each point, type `numpy.ndarray`
* `z` the \(z\) coordinate of each point, type `numpy.ndarray`

## pyssht.spherical_to_cart

~~~~
(x, y, z) = pyssht.spherical_to_cart(r, theta, phi)
~~~~

Computes the \(x\), \(y\), and \(z\) coordinates from the sphrical coordinates \(r\), \(\theta\) and \(\phi\).

#### Inputs

* `r` \(r\) values, type `numpy.ndarray`
* `theta` \(\theta\) values, type `numpy.ndarray`
* `phi` \(\phi\) values, type `numpy.ndarray`

#### Output

Tuple containing `(x, y, z)`
* `x` the \(x\) coordinate of each point, type `numpy.ndarray`
* `y` the \(y\) coordinate of each point, type `numpy.ndarray`
* `z` the \(z\) coordinate of each point, type `numpy.ndarray`

## pyssht.theta_phi_to_ra_dec

~~~~
(dec, ra) = pyssht.theta_phi_to_ra_dec(theta, phi, Degrees=False)
~~~~

Computes the Right Assension and declination from an array of \(\theta\) and \(\phi\) values.

#### Inputs

* `theta` \(\theta\) values, type `numpy.ndarray`
* `phi` \(\phi\) values, type `numpy.ndarray`
* `Degrees` defines if the output is in degrees or radians

#### Output

Tuple contaning `(dec, ra)`
* `dec` the declination angle, `numpy.ndarray`
* `ra` the Right Assension, `numpy.ndarray`


## pyssht.plot_sphere

~~~~
pyssht.plot_sphere(f, L, Method='MW', Close=True, Parametric=False,\
                    Parametric_Saling=[0.0,0.5], Output_File=None,\
                    Show=True, Color_Bar=True, Units=None, Color_Range=None,\ Axis=True)
~~~~

Plots data on to a sphere. It is really slow and not very good!

#### Inputs

* `f` the signal on the sphere, `numpy.ndarray` type `complex` or `real`, ndim 2. NB different for `'MW_pole'` sampling.
* `L` the band limit of the signal, non-zero positive integer
* `Method` the sampling scheme used, string:
    1. `'MW'`         [McEwen & Wiaux sampling (default)]
    2. `'MW_pole'`    [McEwen & Wiaux sampling with the south pole as a seperate double.]
    3. `'MWSS'`       [McEwen & Wiaux symmetric sampling]
    4. `'DH'`         [Driscoll & Healy sampling]
    5. `'GL'`         [Gauss-Legendre sampling]
* `Close` if true the full sphere is plotted (without a gap after the last \(\phi\) position), default `True`
* `Parametric` the radius of the object at a certain point is defined by the function (not just the color), default `False`
* `Parametric_Saling` used if `Parametric=True`, defines the radius of the shape at a particular angle `r = norm(f)*Parametric_Saling[1] + Parametric_Saling[0]`, default `[0.0,0.5]`
* `Output_File` if set saves the plot to a file of that name
* `Show` if `True` shows you the plot, default `False`
* `Color_Bar` if `True` shows the color bar, default `True` 
* `Units` is set puts a label on the color bar
* `Color_Range` if set saturates the color bar in that range, else the function min and max is used
* `Axis` if `True` shows the 3d axis, default `True`

#### Outputs

None

## pyssht.plot_mollweide

~~~~
m = pyssht.plot_mollweide(f, L, Method="MW", Close=True)
~~~~

Plots the data in a mollweide projection. Better but not super fast.

#### Inputs

* `f` the signal on the sphere, `numpy.ndarray` type `complex` or `real`, ndim 2. NB different for `'MW_pole'` sampling.
* `L` the band limit of the signal, non-zero positive integer
* `Method` the sampling scheme used, string:
    1. `'MW'`         [McEwen & Wiaux sampling (default)]
    2. `'MW_pole'`    [McEwen & Wiaux sampling with the south pole as a seperate double.]
    3. `'MWSS'`       [McEwen & Wiaux symmetric sampling]
    4. `'DH'`         [Driscoll & Healy sampling]
    5. `'GL'`         [Gauss-Legendre sampling]
* `Close` if true the full sphere is plotted (without a gap after the last \(\phi\) position), default `True`


#### Outputs

Python axes object. This means the user can call any of the other python ploting functions to costermise the plot.

## pyssht.dl_beta_recurse

~~~~
dl = pyssht.dl_beta_recurse(np.ndarray[ double, ndim=2, mode="c"] dl not None,\
            double beta, int L, int el, \
            np.ndarray[ double, ndim=1, mode="c"] sqrt_tbl not None,\
            np.ndarray[ double, ndim=1, mode="c"] signs not None)
~~~~

Function to recursively calculate the small Wigner D matricies.

Documentation to be added.

## pyssht.generate_dl

~~~~
dl_array = pyssht.generate_dl(double beta, int L)
~~~~

Generates the small Wigner D matricies up to a given band limit for a given \(\beta\)

#### Inputs

* `beta` angle to calculate Wigner D matrix at, type `double`
* `L` the band limit of the signal, non-zero positive integer

#### Outputs

Numpy ndarray `dl_array`, type `float_` of the small Wigner D matricies

## pyssht.rotate_flms

~~~~
flm_rotated = pyssht.rotate_flms(
                np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None,\
                double alpha, double beta, double gamma, int L, dl_array=None,\
                M=None, Axisymmetric=False, Keep_dl=False)
~~~~

Function to rotate a set of spherical harmonic coefficients by the set of Eular angles \(\alpha, \beta, \gamma \) using the \(z,y,z\) convention.

#### Inputs

* `flm` the spherical harmonic transform of `f`, `numpy.ndarray` type `complex`, ndim 1
* `alpha` rotation angle \(\alpha\), type `double`
* `beta` rotation angle \(\beta\), type `double`
* `gamma` rotation angle \(\gamma\), type `double`
* `L` the band limit of the signal, non-zero positive integer
* `dl_array` if set should be the precomuted small Wigner D matirix for angle \(\beta\) and harmonic band limit `L`. If not set this is calculated in the function.
* `M` if set is the azimuthal band limit of the function to be rotated, default `M=L`.
* `Axisymmetric` set if the function is axisymetic and axisytric harmonic coefficients are parsed.
* `Keep_dl` if set the output is changed to allow one to keep the computed `dl_array`

#### Output

If `Keep_dl` is not set the output is the rotated set of spherical harmonic coeficiants. If it is the output is a tuple `(flm_rotated, dl_array)`, ie the rotated harmonic coefficients and the small Wigner D matirix computed for that band limit and \(\alpha\) value.
