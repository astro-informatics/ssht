# SSHT Python Documentation

This guide is intended to explain the python interface of SSHT. For a
description of the workings of SSHT see
[here](http://astro-informatics.github.io/ssht/ "SSHT documentation"). The
python package also offers an interface to some of the functionality from
[ducc0](https://pypi.org/project/ducc0/), including forward and inverse
transforms.

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
    2. `'MW_pole'`    [McEwen & Wiaux sampling with the south pole as a separate double.]
    3. `'MWSS'`       [McEwen & Wiaux symmetric sampling]
    4. `'DH'`         [Driscoll & Healy sampling]
    5. `'GL'`         [Gauss-Legendre sampling]
* `Reality`  determines if the signal is real or complex, Boolean (default = False)
* `backend` the backend that runs the transforms:
    1. `'SSHT'` this package
    2. `'ducc'` interface to [ducc0](https://pypi.org/project/ducc0/). "MW_pole"
       is not available in this backend.
* `nthreads`: number of threads when calling into the `'ducc'` backend. Ignored otherwise.

#### Output

`flm` the spherical harmonic transform of `f`, 1D `numpy.ndarray` type `complex`

#### Note on `'MW_pole'` sampling

This is the same as the `'MW'` sampling however the south pole is expressed as a double only not a vector. Therefore the size of the array is one smaller on the \(\theta\) direction. `f` is now a tuple containing either:

* `(f_array, f_sp, phi_sp)` if complex transform
* `(f_array, f_sp)` if real transform

## pyssht.inverse

~~~~{.python}
f = pyssht.inverse(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, L, Spin=0, Method='MW', Reality=False)
~~~~

Performs the inverse spherical harmonic transform.

#### Inputs

* `flm` the spherical harmonic transform of `f`, `numpy.ndarray` type `complex`, ndim 1
* `L` the band limit of the signal, non-zero positive integer
* `Spin` the spin of the signal, non-negative integer (default = 0)
* `Method` the sampling scheme used, string:
    1. `'MW'`         [McEwen & Wiaux sampling (default)]
    2. `'MW_pole'`    [McEwen & Wiaux sampling with the south pole as a separate double.]
    3. `'MWSS'`       [McEwen & Wiaux symmetric sampling]
    4. `'DH'`         [Driscoll & Healy sampling]
    5. `'GL'`         [Gauss-Legendre sampling]
* `Reality`  determines if the signal is real or complex, Boolean (default = False)
* `backend` the backend that runs the transforms:
    1. `'SSHT'` this package
    2. `'ducc'` interface to [ducc0](https://pypi.org/project/ducc0/). "MW_pole"
       is not available in this backend.
* `nthreads`: number of threads when calling into the `'ducc'` backend. Ignored otherwise.

#### Output

`f` the signal on the sphere, 2D `numpy.ndarray` type `complex` or `real`. NB different for `'MW_pole'` sampling.

#### Note on `'MW_pole'` sampling

This is the same as the `'MW'` sampling however the south pole is expressed as a double only not a vector. Therefore the size of the array is one smaller on the \(\theta\) direction. `f` is now a tuple containing either:

* `(f_array, f_sp, phi_sp)` if complex transform
* `(f_array, f_sp)` if real transform

## pyssht.forward_adjoint

~~~~{.python}
f = pyssht.forward_adjoint(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, L, Spin=0, Method='MW', Reality=False)
~~~~

Performs the adjoint of the forward spherical harmonic transform.

#### Inputs

* `flm` the spherical harmonic transform of `f`, `numpy.ndarray` type `complex`, ndim 1
* `L` the band limit of the signal, non-zero positive integer
* `Spin` the spin of the signal, non-negative integer (default = 0)
* `Method` the sampling scheme used, string:
    1. `'MW'`         [McEwen & Wiaux sampling (default)]
    2. `'MWSS'`       [McEwen & Wiaux symmetric sampling]
* `Reality`  determines if the signal is real or complex, Boolean (default = False)
* `backend` the backend that runs the transforms:
    1. `'SSHT'` this package
    2. `'ducc'` interface to [ducc0](https://pypi.org/project/ducc0/).
* `nthreads`: number of threads when calling into the `'ducc'` backend. Ignored otherwise.

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
    2. `'MWSS'`       [McEwen & Wiaux symmetric sampling]
* `Reality`  determines if the signal is real or complex, Boolean (default = False)
* `backend` the backend that runs the transforms:
    1. `'SSHT'` this package
    2. `'ducc'` interface to [ducc0](https://pypi.org/project/ducc0/).
* `nthreads`: number of threads when calling into the `'ducc'` backend. Ignored otherwise.

#### Output

`flm` the spherical harmonic transform of `f`, 1D `numpy.ndarray` type `complex`

## pyssht.elm2ind

~~~~{.python}
index = pyssht.elm2ind( int el, int m)
~~~~

Computes the index in the `flm` array of a particular harmonic coefficient \(\ell \) and \(m\).

#### Inputs

* `el` the scale parameter of the spherical harmonic coefficients, integer from \(0\) to \(L-1\), where \(L\) is the band limit.
* `em` the azimuthal parameter, integer from -el to el.

#### Output

Index of the coefficient in `flm` array, integer

## pyssht.ind2elm

~~~~{.python}
(el, em) = pyssht.ind2elm(int ind)
~~~~

Computes harmonic coefficient \(\ell \) and \(m\) from the index in the `flm` array.

#### Inputs

* `ind` index of the `flm` array

#### Output

Tuple containing `(el, em)`

* `el` the scale parameter of the spherical harmonic coefficients, integer from \(0\) to \(L-1\), where \(L\) is the band limit.
* `em` the azimuthal parameter, integer from -el to el.

## pyssht.theta_to_index

~~~~{.python}
p = pyssht.theta_to_index(double theta, int L, str Method="MW")
~~~~

Outputs the \(\theta\) index (the first) in the 2 dimensional array used to store spherical images. The index returned is that of the closest \(\theta\) sample smaller then the angle given on input.

#### Inputs

* `theta` the angle \(\theta\)
* `L` the band limit of the signal, non-zero positive integer
* `Method` the sampling scheme used, string:
    1. `'MW'`         [McEwen & Wiaux sampling (default)]
    2. `'MW_pole'`    [McEwen & Wiaux sampling with the south pole as a separate double.]
    3. `'MWSS'`       [McEwen & Wiaux symmetric sampling]
    4. `'DH'`         [Driscoll & Healy sampling]
    5. `'GL'`         [Gauss-Legendre sampling]

#### Output

Int `p` of corresponding to the angle \(\theta\).

## pyssht.phi_to_index

~~~~{.python}
q = pyssht.phi_to_index(double phi, int L, str Method="MW")
~~~~

Outputs the \(\phi\) index (the second) in the 2 dimensional array used to store spherical images. The index returned is that of the closest \(\phi\) sample smaller then the angle given on input.

#### Inputs

* `phi` the angle \(\phi\)
* `L` the band limit of the signal, non-zero positive integer
* `Method` the sampling scheme used, string:
    1. `'MW'`         [McEwen & Wiaux sampling (default)]
    2. `'MW_pole'`    [McEwen & Wiaux sampling with the south pole as a separate double.]
    3. `'MWSS'`       [McEwen & Wiaux symmetric sampling]
    4. `'DH'`         [Driscoll & Healy sampling]
    5. `'GL'`         [Gauss-Legendre sampling]

#### Output

Int `q` of corresponding to the angle \(\phi\).

## pyssht.sample_length

~~~~{.python}
n = pyssht.sample_length(int L, Method = 'MW')
~~~~

Outputs a size of the array used for storing the data on the sphere for different sampling schemes.

#### Inputs

* `L` the band limit of the signal, non-zero positive integer
* `Method` the sampling scheme used, string:
    1. `'MW'`         [McEwen & Wiaux sampling (default)]
    2. `'MW_pole'`    [McEwen & Wiaux sampling with the south pole as a separate double.]
    3. `'MWSS'`       [McEwen & Wiaux symmetric sampling]
    4. `'DH'`         [Driscoll & Healy sampling]
    5. `'GL'`         [Gauss-Legendre sampling]

#### Output

Int `n` equal to the collapsed 1 dimensional size of the 2 dimensional array used to store data on the sphere.

## pyssht.sample_shape

~~~~{.python}
(n_theta, n_phi) = pyssht.sample_shape(int L, Method='MW')
~~~~

Outputs a tuple with the shape of the array used for storing the data on the sphere for different sampling schemes.

#### Inputs

* `L` the band limit of the signal, non-zero positive integer
* `Method` the sampling scheme used, string:
    1. `'MW'`         [McEwen & Wiaux sampling (default)]
    2. `'MW_pole'`    [McEwen & Wiaux sampling with the south pole as a separate double.]
    3. `'MWSS'`       [McEwen & Wiaux symmetric sampling]
    4. `'DH'`         [Driscoll & Healy sampling]
    5. `'GL'`         [Gauss-Legendre sampling]

#### Output

Tuple containing `(n_theta, n_phi)`

* `n_theta` the number of samples in the \(\theta\) direction, integer
* `n_phi` the number of samples in the \(\phi\) direction, integer

## pyssht.sample_positions

~~~~{.python}
(thetas, phis) = pyssht.sample_positions(int L, Method = 'MW', Grid=False)
~~~~

Computes the positions on the sphere of the samples.

#### Inputs

* `L` the band limit of the signal, non-zero positive integer
* `Method` the sampling scheme used, string:
    1. `'MW'`         [McEwen & Wiaux sampling (default)]
    2. `'MW_pole'`    [McEwen & Wiaux sampling with the south pole as a separate double.]
    3. `'MWSS'`       [McEwen & Wiaux symmetric sampling]
    4. `'DH'`         [Driscoll & Healy sampling]
    5. `'GL'`         [Gauss-Legendre sampling]
* `Grid` describes if the output is a vector of the sample positions or a 2D array the same shape as the signal on the sphere, default `False`

#### Outputs

Tuple containing `(thetas, phis)`

* `thetas` positions of the samples in the \(\theta\) direction
* `phis` positions of the samples in the \(\theta\) direction

## pyssht.s2_to_cart

~~~~{.python}
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

## pyssht.cart_to_s2

~~~~{.python}
thetas, phis = cart_to_s2(x, y, z)
~~~~

Computes the \(\theta\) and \(\phi\) on the coordinates on the sphere from \(x\), \(y\), and \(z\) coordinates.

#### Inputs

* `x` \(x\) values, type `numpy.ndarray`
* `y` \(y\) values, type `numpy.ndarray`
* `z` \(z\) values, type `numpy.ndarray`

#### Output

Tuple containing `(theta, phi)`

* `theta` the \(\theta\) coordinate of each point, type `numpy.ndarray`
* `phi` the \(\phi\) coordinate of each point, type `numpy.ndarray`

## pyssht.spherical_to_cart

~~~~{.python}
(x, y, z) = pyssht.spherical_to_cart(r, theta, phi)
~~~~

Computes the \(x\), \(y\), and \(z\) coordinates from the spherical coordinates \(r\), \(\theta\) and \(\phi\).

#### Inputs

* `r` \(r\) values, type `numpy.ndarray`
* `theta` \(\theta\) values, type `numpy.ndarray`
* `phi` \(\phi\) values, type `numpy.ndarray`

#### Output

Tuple containing `(x, y, z)`

* `x` the \(x\) coordinate of each point, type `numpy.ndarray`
* `y` the \(y\) coordinate of each point, type `numpy.ndarray`
* `z` the \(z\) coordinate of each point, type `numpy.ndarray`

## pyssht.cart_to_spherical

~~~~{.python}
r, theta, phi = pyssht.cart_to_spherical(x, y, z)
~~~~

Computes the \(\r\), \(\theta\) and \(\phi\) on the spherical coordinates from \(x\), \(y\), and \(z\) coordinates.

#### Inputs

* `x` \(x\) values, type `numpy.ndarray`
* `y` \(y\) values, type `numpy.ndarray`
* `z` \(z\) values, type `numpy.ndarray`

#### Output

Tuple containing `(r, theta, phi)`

* `r` the \(r\) coordinate of each point, type `numpy.ndarray`
* `theta` the \(\theta\) coordinate of each point, type `numpy.ndarray`
* `phi` the \(\phi\) coordinate of each point, type `numpy.ndarray`

## pyssht.theta_phi_to_ra_dec

~~~~{.python}
(dec, ra) = pyssht.theta_phi_to_ra_dec(theta, phi, Degrees=False)
~~~~

Computes the Right Assension and declination from an array of \(\theta\) and \(\phi\) values.

#### Inputs

* `theta` \(\theta\) values, type `numpy.ndarray`
* `phi` \(\phi\) values, type `numpy.ndarray`
* `Degrees` defines if the output is in degrees or radians

#### Output

Tuple containing `(dec, ra)`

* `dec` the declination angle, `numpy.ndarray`
* `ra` the Right Assension, `numpy.ndarray`

## pyssht.ra_dec_to_theta_phi

~~~~{.python}
(theta, phi) = pyssht.ra_dec_to_theta_phi(ra, dec, Degrees=False)
~~~~

Computes the \(\theta\) and \(\phi\) values from an array of Right Assension and declination values.

#### Inputs

* `dec` the declination angle, `numpy.ndarray`
* `ra` the Right Assension, `numpy.ndarray`
* `Degrees` defines if the input is in degrees or radians, if degrees they are converted

#### Output

Tuple containing `(theta, phi])`

* `theta` \(\theta\) values, type `numpy.ndarray`
* `phi` \(\phi\) values, type `numpy.ndarray`

## pyssht.make_rotation_matrix

~~~~{.python}
rot_much = pyssht.make_rotation_matrix(list rot)
~~~~

Computes the 3 by 3 rotation matrix from the Euler angles given on input

#### Inputs

* `rot` List of length 3. Each element are the Euler angles `[alpha, beta, gamma]`

#### Output

3 by 3 `rot_matrix` the rotation matrix type `ndarray` dtype `float`

## pyssht.rot_cart

~~~~{.python}
x_p, y_p, z_p = pyssht.rot_cart(x, y, z, list rot)
~~~~

Computes the rotations of the cartesian coordinates given a set of Euler angles. The inputs can be any shape `ndarray`s. For speed if the arrays are 1 or 2 dimensional it is recommended to use `pyssht.rot_cart_1D` or `pyssht.rot_cart_2D`.

#### Inputs

* `x` \(x\) values, type `numpy.ndarray`
* `y` \(y\) values, type `numpy.ndarray`
* `z` \(z\) values, type `numpy.ndarray`
* `rot` List of length 3. Each element are the Euler angles `[alpha, beta, gamma]`

#### Output

Tuple containing `(x_p, y_p, z_p)` the rotated coordinates the same shape and type as the inputs.

## pyssht.rot_cart_1d and pyssht.rot_cart_2d

~~~~{.python}
(x_p, y_p, z_p) = pyssht.rot_cart_1d(np.ndarray[np.float_t, ndim=1] x, np.ndarray[np.float_t, ndim=1] y, np.ndarray[np.float_t, ndim=1] z, list rot)
~~~~

~~~~{.python}
(x_p, y_p, z_p) = pyssht.rot_cart_2d(np.ndarray[np.float_t, ndim=2] x, np.ndarray[np.float_t, ndim=2] y, np.ndarray[np.float_t, ndim=2] z, list rot)
~~~~

Computes the rotations of the cartesian coordinates given a set of Euler angles. The inputs can be any shape `ndarray`s. Same as `pyssht.rot_cart` except optimised for arrays that are 1 or 2 dimensional.

#### Inputs

* `x` \(x\) values, type `numpy.ndarray`, dtype `float`, ndim 1 or 2
* `y` \(y\) values, type `numpy.ndarray`, dtype `float`, ndim 1 or 2
* `z` \(z\) values, type `numpy.ndarray`, dtype `float`, ndim 1 or 2
* `rot` List of length 3. Each element are the Euler angles `[alpha, beta, gamma]`

#### Output

Tuple containing `(x_p, y_p, z_p)` the rotated coordinates the same shape and type as the inputs.

## pyssht.plot_sphere

~~~~{.python}
pyssht.plot_sphere(
    f, L, Method='MW', Close=True, Parametric=False,
    Parametric_Scaling=[0.0,0.5], Output_File=None,
    Show=True, Color_Bar=True, Units=None, Color_Range=None,
    Axis=True
)
~~~~

Plots data on to a sphere. It is really slow and not very good!

#### Inputs

* `f` the signal on the sphere, `numpy.ndarray` type `complex` or `real`, ndim 2. NB different for `'MW_pole'` sampling.
* `L` the band limit of the signal, non-zero positive integer
* `Method` the sampling scheme used, string:
    1. `'MW'`         [McEwen & Wiaux sampling (default)]
    2. `'MW_pole'`    [McEwen & Wiaux sampling with the south pole as a separate double.]
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

## pyssht.mollweide_projection

~~~~{.python}
f_plot, mask = pyssht.mollweide_projection(
    f, L, resolution=500, rot=None,
    zoom_region=[np.sqrt(2.0)*2,np.sqrt(2.0)],
    Method="MW"
)
~~~~

Creates an `ndarray` of the mollweide projection of a spherical image and a mask array. This is useful for plotting results, not to be used for analysis on the plane. Elements in the signal `f` that are `NaN`s are marked in the mask. This allows one to plot these regions the color of their choice.

Here is an example of using the function to plot real spherical data.

~~~~{.python}
f_plot, mask = pyssht.mollweide_projection(f, L, Method="MW") # make projection
plt.figure() # start figure
imgplot = plt.imshow(f_real_plot,interpolation='nearest') # plot the projected image
plt.colorbar(imgplot,fraction=0.025, pad=0.04) # plot color bar (these extra keywords make the bar a reasonable size)
plt.imshow(mask_real, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.) # plot the NaN regions in grey
plt.gca().set_aspect("equal") # ensures the region is the correct proportions
plt.axis('off') # removes axis (looks better)
plt.show()
~~~~

#### Inputs

* `f` the signal on the sphere, `numpy.ndarray` type `complex` or `real`, ndim 2.
* `L` the band limit of the signal, non-zero positive integer
* `resolution` size of the projected image, default 500
* `rot` If the image should be rotated before projecting, None. `rot` should be a list of length 1 or 3. If 1 then the image is rotated around the \(z\) axis by that amount. If 3 then the image is rotated by the Euler angles given in the list.
* `zoom_region` the region of the sphere to be plotted, default `[np.sqrt(2.0)*2,np.sqrt(2.0)]` is the full sphere.
* `Method` the sampling scheme used, string:
    1. `'MW'`         [McEwen & Wiaux sampling (default)]
    2. `'MWSS'`       [McEwen & Wiaux symmetric sampling]
    3. `'DH'`         [Driscoll & Healy sampling]
    4. `'GL'`         [Gauss-Legendre sampling]

#### Outputs

If the input is real:

Tuple containing:

* `f_plot` the projection of the image as a 2 dimensional `ndarray` of type `float`. Masked regions and regions not in the sphere are `NaN`s to make them clear when plotted
* `mask` the projection of the masked regions (`NaN`s in input `f`) as a 2 dimensional `ndarray` of type `float`. Masked regions have a value `0.0` and regions not in the sphere are `NaN`s to make them clear when plotted.

If the input is complex:

Tuple containing:

* `f_plot_real` the projection of the real part of the image as a 2 dimensional `ndarray` of type `float`. Masked regions and regions not in the sphere are `NaN`s to make them clear when plotted
* `mask_real` the projection of the masked regions (`NaN`s in input `f.real`) as a 2 dimensional `ndarray` of type `float`. Masked regions have a value `0.0` and regions not in the sphere are `NaN`s to make them clear when plotted.
* `f_plot_imag` the projection of the imaginary part of the image as a 2 dimensional `ndarray` of type `float`. Masked regions and regions not in the sphere are `NaN`s to make them clear when plotted
* `mask_imag` the projection of the masked regions (`NaN`s in input `f.imag`) as a 2 dimensional `ndarray` of type `float`. Masked regions have a value `0.0` and regions not in the sphere are `NaN`s to make them clear when plotted.

## pyssht.equatorial_projection

~~~~{.python}
f_proj_real, mask_real, (f_proj_imag, mask_imag)\
            = pyssht.equatorial_projection(f, int L, int resolution=500,\
             rot=None, list zoom_region=[-1,-1], str Method="MW", \
             str Projection="MERCATOR", int Spin=0)
~~~~

Creates `ndarray`s of the projections of a spherical image and a mask array. This is useful for plotting results and performing analysis on the plane. All the spherical samples that fall in one planar pixel is averaged, if no samples fall in a pixel then the pixel is assigned the value of the closest spherical sample. Elements in the signal `f` that are `NaN`s are marked in the mask. This allows one to plot these regions the color of their choice.

There are two projections supported.

1. The Mercator projection often used in maps.
2. The Sinusoidal projection a simple equal area projection.

Here is an example of using the function to plot real spherical data using the Mercator projection.

~~~~{.python}
f_proj, mask \
        = pyssht.equatorial_projection(f, L, resolution=500, Method="MW", \
        Projection="MERCATOR")
plt.figure() # start figure
imgplot = plt.imshow(f_proj,interpolation='nearest')# plot the projected image (north part)
plt.colorbar(imgplot) # plot color bar
plt.imshow(mask, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.) # plot the NaN regions in grey
plt.axis('off') # removes axis

plt.show()
~~~~

#### Inputs

* `f` the signal on the sphere, `numpy.ndarray` type `complex` or `real`, ndim 2.
* `L` the band limit of the signal, non-zero positive integer
* `resolution` size of the projected image, default 500
* `rot` If the image should be rotated before projecting, None. `rot` should be a list of length 1 or 3. If 1 then the image is rotated around the \(z\) axis by that amount. If 3 then the image is rotated by the Euler angles given in the list.
* `zoom_region` the region of the sphere to be plotted in radians. The first element is the angle left and right of the centre, default is `np.pi` for both projections. The second element is up and down of the equator, default is `np.pi/2` for the Sinusoidal projection and `7*np.pi/16` for the Mercator projection.
* `Method` the sampling scheme used, string:
    1. `'MW'`         [McEwen & Wiaux sampling (default)]
    2. `'MWSS'`       [McEwen & Wiaux symmetric sampling]
    3. `'DH'`         [Driscoll & Healy sampling]
    4. `'GL'`         [Gauss-Legendre sampling]
* `Projection` string describing which of the projections to use. Use `"MERCATOR"` for the Mercator projection and `"SINE"` for the Sinusoidal projection, default is `"MERCATOR"`
* `Spin` the spin of the signal. If the signal has non-zero spin then on projection the signal must be rotated to account for the changing direction of the definition of the signal. By setting this to a non-zero integer will ensure this rotation is performed.

#### Outputs

If the input is real:

Tuple containing:

* `f_plot` the projection of the image as a 2 dimensional `ndarray` of type `float`. Masked regions and regions not in the sphere are `NaN`s to make them clear when plotted
* `mask` the projection of the masked regions (`NaN`s in input `f`) as a 2 dimensional `ndarray` of type `float`. Masked regions have a value `0.0` and regions not in the sphere are `NaN`s to make them clear when plotted.

If the input is complex:

Tuple containing:

* `f_plot_real` the projection of the real part of the image as a 2 dimensional `ndarray` of type `float`. Masked regions and regions not in the sphere are `NaN`s to make them clear when plotted
* `mask_real` the projection of the real part of the masked regions (`NaN`s in input `f`) as a 2 dimensional `ndarray` of type `float`. Masked regions have a value `0.0` and regions not in the sphere are `NaN`s to make them clear when plotted.
* `f_plot_imag` the projection of the imaginary part of the image as a 2 dimensional `ndarray` of type `float`. Masked regions and regions not in the sphere are `NaN`s to make them clear when plotted
* `mask_imag` the projection of the imaginary part of the masked regions (`NaN`s in input `f`) as a 2 dimensional `ndarray` of type `float`. Masked regions have a value `0.0` and regions not in the sphere are `NaN`s to make them clear when plotted.

## pyssht.polar_projection

~~~~{.python}
f_proj_north_real, mask_north_real, f_proj_south_real, mask_south_real,\
(f_proj_north_imag, mask_north_imag, f_proj_south_imag, mask_south_imag\ )
        = pyssht.polar_projection(f, int L, int resolution=500, rot=None,\
        float zoom_region=-1, str Method="MW", str Projection="OP",int Spin=0):
~~~~

Creates an two `ndarray`s of the polar projection of a spherical image and a mask array. This is useful for plotting results and performing analysis on the plane. All the spherical samples that fall in one planar pixel is averaged, if no samples fall in a pixel then the pixel is assigned the value of the closest spherical sample. Elements in the signal `f` that are `NaN`s are marked in the mask. This allows one to plot these regions the color of their choice.

All the projections are centred around a pole. There are three projections supported.

1. The Gnomic projection, defined by drawing a line from the centre of the circle trough the sphere on to the plane.
2. The Stereographic projection, defined by drawing a line starting at the opposite pole through the sphere to the plane.
3. The Orthographic, defined by a vertical projection on the plane.

Here is an example of using the function to plot real spherical data.

~~~~{.python}
f_proj_north, mask_north, f_proj_south, mask_south \
        = pyssht.polar_projection(f, L, resolution=500, Method="MW", \
        Projection="OP")
plt.figure() # start figure
imgplot = plt.imshow(f_proj_north,interpolation='nearest')# plot the projected image (north part)
plt.colorbar(imgplot) # plot color bar
plt.imshow(mask_north, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.) # plot the NaN regions in grey
plt.axis('off') # removes axis (looks better)

plt.figure()
imgplot = plt.imshow(f_proj_south,interpolation='nearest')# plot the projected image (south part)
plt.colorbar(imgplot)
plt.imshow(mask_south, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
plt.title("orthographic projection south")
plt.axis('off')
plt.show()
~~~~

#### Inputs

* `f` the signal on the sphere, `numpy.ndarray` type `complex` or `real`, ndim 2.
* `L` the band limit of the signal, non-zero positive integer
* `resolution` size of the projected image, default 500
* `rot` If the image should be rotated before projecting, default None. `rot` should be a list of length 1 or 3. If 1 then the image is rotated around the \(z\) axis by that amount. If 3 then the image is rotated by the Euler angles given in the list.
* `zoom_region` the region of the sphere to be plotted in radians, default `np.pi/2` is the full half sphere for the orthographic and stereographic projections and `np.pi/4` for the gnomic projection as the equator is at infinity in this projection.
* `Method` the sampling scheme used, string:
    1. `'MW'`         [McEwen & Wiaux sampling (default)]
    2. `'MWSS'`       [McEwen & Wiaux symmetric sampling]
    3. `'DH'`         [Driscoll & Healy sampling]
    4. `'GL'`         [Gauss-Legendre sampling]
* `Projection` string describing which of the projections to use. Use `"GP"` for the Gnomic projection, `"SP"` for the Stereographic projection and `"OP"` for the Orthographic projection, default is `"OP"`
* `Spin` the spin of the signal. If the signal has non-zero spin then on projection the signal must be rotated to account for the changing direction of the definition of the signal. By setting this to a non-zero integer will ensure this rotation is performed.

#### Outputs

If the input is real:

Tuple containing:

* `f_plot_north` the projection of the north part of the image as a 2 dimensional `ndarray` of type `float`. Masked regions and regions not in the sphere are `NaN`s to make them clear when plotted
* `mask_north` the projection of the north part of the masked regions (`NaN`s in input `f`) as a 2 dimensional `ndarray` of type `float`. Masked regions have a value `0.0` and regions not in the sphere are `NaN`s to make them clear when plotted.
* `f_plot_south` the projection of the south part of the image as a 2 dimensional `ndarray` of type `float`. Masked regions and regions not in the sphere are `NaN`s to make them clear when plotted
* `mask_south` the projection of the south part of the masked regions (`NaN`s in input `f`) as a 2 dimensional `ndarray` of type `float`. Masked regions have a value `0.0` and regions not in the sphere are `NaN`s to make them clear when plotted.

If the input is complex:

Tuple containing:

* `f_plot_north_real` the projection of the north part of the real part of the image as a 2 dimensional `ndarray` of type `float`. Masked regions and regions not in the sphere are `NaN`s to make them clear when plotted
* `mask_north_real` the projection of the north part of the real part of the masked regions (`NaN`s in input `f`) as a 2 dimensional `ndarray` of type `float`. Masked regions have a value `0.0` and regions not in the sphere are `NaN`s to make them clear when plotted.
* `f_plot_south_real` the projection of the south part of the real part of the image as a 2 dimensional `ndarray` of type `float`. Masked regions and regions not in the sphere are `NaN`s to make them clear when plotted
* `mask_south_real` the projection of the south part of the real part of the masked regions (`NaN`s in input `f`) as a 2 dimensional `ndarray` of type `float`. Masked regions have a value `0.0` and regions not in the sphere are `NaN`s to make them clear when plotted.
* `f_plot_north_imag` the projection of the north part of the imaginary part of the image as a 2 dimensional `ndarray` of type `float`. Masked regions and regions not in the sphere are `NaN`s to make them clear when plotted
* `mask_north_imag` the projection of the north part of the imaginary part of the masked regions (`NaN`s in input `f`) as a 2 dimensional `ndarray` of type `float`. Masked regions have a value `0.0` and regions not in the sphere are `NaN`s to make them clear when plotted.
* `f_plot_south_imag` the projection of the south part of the imaginary part of the image as a 2 dimensional `ndarray` of type `float`. Masked regions and regions not in the sphere are `NaN`s to make them clear when plotted
* `mask_south_imag` the projection of the south part of the imaginary part of the masked regions (`NaN`s in input `f`) as a 2 dimensional `ndarray` of type `float`. Masked regions have a value `0.0` and regions not in the sphere are `NaN`s to make them clear when plotted.

## pyssht.dl_beta_recurse

~~~~{.python}
dl = pyssht.dl_beta_recurse(np.ndarray[ double, ndim=2, mode="c"] dl not None,\
            double beta, int L, int el, \
            np.ndarray[ double, ndim=1, mode="c"] sqrt_tbl not None,\
            np.ndarray[ double, ndim=1, mode="c"] signs not None)
~~~~

Compute the el-th plane of the Wigner small-d functions (from the (el-1)-th plane) using Risbo's method.

#### Inputs

* `dl` the Wigner plane for all m and n, indexed dl[m][n] of size (2*L-1)*(2*L-1)
* `beta` angle to calculate Wigner D matrix at, type `double`
* `L` the band limit of the signal, non-zero positive integer
* `el` is the current harmonic degree (i.e. dl input should already be computed for el-1, and dl output will be computed for el)
* `sqrt_tbl` precomputed square-roots from \(0\) to \(2*(L-1)+1\)
* `signs` precomputed \((-1)^m\) signs from \(m=0\) to \(L\)

#### Outputs

Numpy ndarray `dl`, type `float_` of the Wigner plane for all m and n

## pyssht.dln_beta_recurse

~~~~{.python}
dl = pyssht.dln_beta_recurse(np.ndarray[ double, ndim=1, mode="c"] dl not None,\
            np.ndarray[ double, ndim=1, mode="c"] dlm1 not None, double beta,\
            int L, int el, int n, np.ndarray[ double, ndim=1, mode="c"] sqrt_tbl not None,\
            np.ndarray[ double, ndim=1, mode="c"] signs not None)
~~~~

Compute the el-th line of the Wigner small-d functions for given n (from the (el-1)-th and (el-2)-th lines) using 3-term recursion of Kostelec.

#### Inputs

* `dl` the Wigner line for el for non-negative m and given n of size L
* `dlm1` is the line for el-1 and dlp1 is the line computed for el+1
* `beta` angle to calculate Wigner D matrix at, type `double`
* `L` the band limit of the signal, non-zero positive integer
* `el` el is the current harmonic degree
* `n` the third index in Wigner D matrices
* `sqrt_tbl` precomputed square-roots from \(0\) to \(2*(L-1)+1\)
* `signs` precomputed \((-1)^m\) signs from \(m=0\) to \(L\)

#### Outputs

Numpy ndarray `dl`, type `float_` the Wigner line for el for non-negative m and given n of size

## pyssht.generate_dl

~~~~{.python}
dl_array = pyssht.generate_dl(double beta, int L)
~~~~

Generates the small Wigner D matrices up to a given band limit for a given \(\beta\)

#### Inputs

* `beta` angle to calculate Wigner D matrix at, type `double`
* `L` the band limit of the signal, non-zero positive integer

#### Outputs

Numpy ndarray `dl_array`, type `float_` of the small Wigner D matrices

## pyssht.rotate_flms

~~~~{.python}
flm_rotated = pyssht.rotate_flms(
                np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None,\
                double alpha, double beta, double gamma, int L, dl_array=None,\
                M=None, Axisymmetric=False, Keep_dl=False)
~~~~

Function to rotate a set of spherical harmonic coefficients by the set of Euler angles \(\alpha, \beta, \gamma \) using the \(z,y,z\) convention.

#### Inputs

* `flm` the spherical harmonic transform of `f`, `numpy.ndarray` type `complex`, ndim 1
* `alpha` rotation angle \(\alpha\), type `double`
* `beta` rotation angle \(\beta\), type `double`
* `gamma` rotation angle \(\gamma\), type `double`
* `L` the band limit of the signal, non-zero positive integer
* `dl_array` if set should be the precomputed small Wigner D matrix for angle \(\beta\) and harmonic band limit `L`. If not set this is calculated in the function. (This parameter is ignored when using the `ducc` backend.)
* `M` if set is the azimuthal band limit of the function to be rotated, default `M=L`.
* `Axisymmetric` set if the function is axisymmetric and axisymmetric harmonic coefficients are parsed.
* `Keep_dl` if set the output is changed to allow one to keep the computed `dl_array`. (This parameter is ignored when using the `ducc` backend.)
* `backend` the backend that runs the transforms:
    1. `'SSHT'` this package
    2. `'ducc'` interface to [ducc0](https://pypi.org/project/ducc0/)

#### Output

If `Keep_dl` is not set the output is the rotated set of spherical harmonic coefficients. If it is the output is a tuple `(flm_rotated, dl_array)`, ie the rotated harmonic coefficients and the small Wigner D matrix computed for that band limit and \(\alpha\) value.

## pyssht.guassian_smoothing

~~~~{.python}
fs_lm = pyssht.guassian_smoothing(np.ndarray[ double complex, ndim=1, mode="c"] f_lm not None, int L, sigma_in=None, bl_in = None)
~~~~

Smooths a set of harmonic coefficients either with a precomputed smoothing kernel `bl` or with a Gaussian given on input.

#### Inputs

* `f_lm` the spherical harmonic transform of `f`, `numpy.ndarray` type `complex`, ndim 1
* `L` the band limit of the signal, non-zero positive integer
* `sigma_in` the input sigma of the Gaussian to smooth the signal with, default `None`
* `bl_in` the smoothing kernel to smooth the signal with, default `None`

#### Output

`fs_lm` the smoothed harmonic coefficients.

## pyssht.create_ylm

~~~~{.python}
ylm = pyssht.create_ylm(thetas, phis, int L, int Spin=0, str recursion='Kostelec')
~~~~

Computes spherical harmonic functions for all el and all 0<=|m|<= el using various recursions.

#### Inputs

* `thetas` positions of the samples in the \(\theta\) direction
* `phis` positions of the samples in the \(\phi\) direction
* `L` the band limit of the signal, non-zero positive integer
* `Spin` the spin of the signal, non-negative integer (default = 0)
* `recursion` the recursion scheme used, string:
    1. `'Kostelec'` [3-term recursion, e.g. Kostelec (default)]
    2. `'Risbo'` [Risbo recursion]
    3. `'NumericalRecipes'` [Numerical Recipes]

#### Output

`ylm` the spherical harmonics indexed `ylm[ind][theta][phi]`, where `ind = pyssht.elm2ind(el, m)`
