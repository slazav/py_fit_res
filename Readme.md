## py_fit_res -- python library for fitting linear and non-linear resonances

<img align=right height="250"
src="https://raw.githubusercontent.com/slazav/py_fit_res/master/examples/title.png">

Features:

* Linear oscillator, Duffing oscillator, oscillator in ballistic B-phase
* Fitting of X and Y data components vs frequency
* Fitting of coordinate or velocity response
* Optional fit parameters for constant and linear backgrounds
* A single fit of multiple frequency sweeps and multiple drives is possible
* Internal data scaling (to have all parameters of the order of 1 during fitting)
* Non-linear resonance proccessing according with
  ![PDF](https://raw.githubusercontent.com/slazav/py_fit_res/master/fit_res_note.pdf)

Comments and suggestions are very welcome!

V.Zavjalov (vl.zavjalov at gmail dot com), 04.05.2024

----
### Linear oscillator

Basic usage:

Do the fit and return result object:
```
fit = fit_res.fit_lin(FF,XX,YY, <parameters>)
```

Calculate model function:
```
vv = fit.func(ff)
```

Get resonance frequency, width, and amplitude:
```
f0  = fit.get_f0(ff)
df  = fit.get_df(ff)
amp = fit.get_amp(ff)
```

#### Parameters of `fit_res.fit_lin` function

*   `FF`    -- Frequency array, numpy array
*   `XX,YY` -- X and Y data components, numpy arrays
*   `DD=1`  -- Drive array or number. Use to fit multiple drives simultaneously
*   `coord=1`     -- switch between coordinate/velocity response (`1i*f` factor)
*   `cbg0=1` -- use background idependent on both drive and freqeuncy (2 extra parameters)
*   `cbg=1`  -- use background proportional to drive and independent on frequency (2 extra parameters)
*   `lbg=1`  -- use background proportional to drive and linear in frequency (2 extra parameters)
*   `do_fit=1`    -- do actual fitting or return initial conditions (useful for tests)
*   `fit_displ=None, fit_maxiter=10000` -- parameters passed to `scipy.optimize.minimize`

No not use `cbg0` and `cbg` together unless you have multiple drives in your data

#### Fitting model and free parameters

* Coordinate response:
```
XX+iYY = amp*DD * df*f0 / (f0^2 - FF^2 + i*FF*df) [+cbg0] [+cbg*DD] [+lbg*DD*(FF-f0)]
```

* Velocity response:
```
XX+iYY = amp*DD * i*df*FF / (f0^2 - FF^2 + i*FF*df) [+cbg0] [+cbg*DD] [+lbg*DD*(FF-f0)]
```
where

*   `amp` -- complex amplitude per unit drive (same units as `XX/DD` and `YY/DD`)
*   `f0, df` -- resonance frequency and width (same units as `FF`)
*   `cbg0` -- complex constant drive-independent background per unit drive (only if `cbg0=1`)
*   `cbg`  -- complex constant background per unit drive (only if `cbg=1`)
*   `lbg`  -- complex linear background per unit drive (only if `lbg=1`)

For each parameter two functions are defined:

* `get_<name>(p=None)` -- extract value of parameter `<name>` from parameter array `p`
  (if not None) or from fit result.

* `get_<name>_e(e=None)` -- extract uncertainty of parameter `<name>` from uncertainty
  array `e` (if not None) or from the fit result.

#### Functions

* `fit.func(FF,DD=1,p=None)` -- calculate the model function for parameter array `p`
  (if not None) or for the fit result.

* `fit.func_bg(FF,DD=1,p=None)` -- calculate background part of the model function.

----
### Duffing oscillator

```
fit = fit_res.fit_duff(FF,XX,YY, <parameters>)
```

All parameters and usage is same as for the linear oscillator. The only difference is
the extra fitting parameter `a` which is a factor in the Duffing force `F = -a*x^3`.
Additional functions: `get_a`, `get_a_e`.

Some notes:
* If frequency units are switched from rad/s to Hz then units of `A` changes as well.
  Also note amplitude and drive units for understanding quantitative value of `A`.

* The response with hysteresis should be fitted properly (see examples below),
  order of the frequency array is important for calculations.

* In real life you can rarely observe pure Duffing systems, in many cases this fit
  will be useful only for small non-linearities.

----
### Oscillator in ballistic B phase

```
fit = fit_res.fit_bphase(FF,XX,YY, <parameters>)
```

All parameters and usage is same as for the linear oscillator. The only difference is
the extra fitting parameter `v0` which is a characteristic velocity (in XX and YY units).
Additional functions: `get_v0`, `get_v0_e`.

For the non-linear damping function the following approximation is used:
```
df_n(|v|) = df / (1 + 0.447*(|v|/v0)**1.16)
```

----
### Arbitrary non-linear oscillator
```
fit = fit_res.fit_nonlin(FF,XX,YY, ffunc=None, dfunc=None, <parameters>)
```

Use arbitrary non-linear functions `f0n(|x|)` and `dfn(|v|)`.

The functions can be calculated as integral transforms of
the non-linear force acting on the oscillator using following interface:

```
transform(func, x, npts=100, dfunc0=None)
itransform(func, dfunc, x, npts=100)
```

See the PDF note for the theory

----
#### Examples

See `examples` folder.

As experimental data I use measurements of vibrating wires done during
my work in Lancaster University. I have chosen data which can illustrate
linear/duffing/B-phase behaviour, but one can not expect exact match
with theoretical models.

Example 1: Linear oscillator, multiple frequency sweeps at a single drive

![example 1](https://raw.githubusercontent.com/slazav/py_fit_res/master/examples/example1.png)

Example 2: Linear oscillator, multiple frequency sweeps at multiple drives

![example 2](https://raw.githubusercontent.com/slazav/py_fit_res/master/examples/example2.png)

Example 3: Duffing oscillator, multiple frequency sweeps at multiple drives.
Note that experimental system is not exactly a Duffing oscillator, perfect match is not expected.

![example 3](https://raw.githubusercontent.com/slazav/py_fit_res/master/examples/example3.png)

Example 4: Oscillator in ballistic B-phase, multiple frequency sweeps at multiple drives.

![example 4](https://raw.githubusercontent.com/slazav/py_fit_res/master/examples/example4.png)

Example 5: Same data as in the example 4. 

Examples 6 and 7: more realistic examples for the note (see below)

#### See also:


* `https://github.com/slazav/fit_res` -- my c++ command line filter for fitting linear resonances.
* fit_res modules in `https://github.com/slazav/exp_py` -- older version of the same code.
* `https://raw.githubusercontent.com/slazav/py_fit_res/master/fit_res_note.pdf` -- the note with some theory and examples :
