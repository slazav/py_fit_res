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

Comments and suggestions are very welcome!

V.Zavjalov (vl.zavjalov at gmail dot com), 04.05.2024

----
### Linear oscillator

Do the fit and return result object:
```
fit = fit_res.fit_lin(FF,XX,YY, <parameters>)
```

Calculate model function:
```
vv = fit.func(ff)
```

#### Parameters

*   `FF`    -- Frequency array, numpy array
*   `XX,YY` -- X and Y data components, numpy arrays
*   `DD=1`  -- Drive array or number. Use to fit multiple drives simultaneously
*   `coord=1`     -- switch between coordinate/velocity response (`1i*f` factor)
*   `const_bg=1`  -- use constant background (2 extra fitting parameters)
*   `linear_bg=1` -- use linear background (2 extra fitting parameters)
*   `do_fit=1`    -- do actual fitting or return initial conditions (useful for tests)
*   `fit_displ=None, fit_maxiter=10000` -- parameters passed to `scipy.optimize.minimize`


#### Fitting model and free parameters

* Coordinate response:
```
XX+iYY = Amp * DD * df*f0 / (f0^2 - FF^2 + i*FF*df) [ + cbg*DD] [ + lbg*DD*(FF-f0)]
```

* Velocity response:
```
XX+iYY = Amp*DD * i*df*FF / (f0^2 - FF^2 + i*FF*df) [ + Cbg*DD] [ + lbg*DD*(FF-f0)]
```
where

*   `Amp` -- complex amplitude per unit drive (same units as `XX/DD` and `YY/DD`)
*   `f0, df` -- resonance frequency and width (same units as `FF`)
*   `cbg` -- complex constant background per unit drive (only if `const_bg=1`)
*   `lbg` -- complex linear background per unit drive (only if `linear_bg=1`)


#### Parameter array (4, 6, or 8 values)

  `fit.pars` -- `Re(Amp) Im(Amp) F0 dF [Re(cbg) Im(cbg)] [Re(lbg) Im(lbg)]`


#### Uncertainties array

  `fit.errs` -- same order as for `fit.pars`


#### Functions

* `fit.func(FF,DD=1)` -- calculate the model function
* `fit.func_bg(FF,DD=1)` -- background part of the model function
* `fit.fitfunc(par, FF, DD=1)` -- calculate the model function using custom parameters


----
### Duffing oscillator

```
fit = fit_res.fit_duff(FF,XX,YY, <parameters>)
```

All parameters and usage is same as for the linear oscillator. The only difference is
the extra fitting parameter A which is a factor in the Duffing force `F = -A*x^3`.

Some notes:
* If frequency units are switched from rad/s to Hz then units of `A` changes as well.
  Also note amplitude and drive units for understanding quantitative value of `A`.

* The response with hysteresis should be fitted properly (see examples below),
  order of the frequency array is important for calculations.

* In real life you can rarely observe pure Duffing systems, in many cases this fit
  will be useful only for small non-linearities.

Parameter array:

`fit.pars` -- `Re(Amp) Im(Amp) F0 dF A [Re(cbg) Im(cbg)] [Re(lbg) Im(lbg)]`

----
### Oscillator in ballistic B phase

```
fit = fit_res.fit_bphase(FF,XX,YY, <parameters>)
```

All parameters and usage is same as for the linear oscillator. The only difference is
the extra fitting parameter v0 which is a characteristic velocity (in XX and YY units).

Parameter array:

`fit.pars` -- `Re(Amp) Im(Amp) F0 dF V0 [Re(cbg) Im(cbg)] [Re(lbg) Im(lbg)]`

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

#### See also:

* fit_res modules in `https://github.com/slazav/exp_py`: older version of the same code.
* `https://github.com/slazav/fit_res`: my c++ command line filter for fitting linear resonances.
