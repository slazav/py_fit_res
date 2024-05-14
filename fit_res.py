# fit_res -- a python library fir fitting resonance data
# V.Zavjalov, 2023-2024

import numpy
import math
import scipy.optimize


###############################################################
# Initial conditions for a linear resonance with linear background (8 pars)
#    V(f) = f0*df*(C+iD)/(f0^2-f^2 + i*f*df) + (A+iB) + (E+iF)*(f-f0)
def init8(FF,XS,YS,coord):

  # points with min/max freq (note that frequency could be non-monotonic)
  ifmin = numpy.argmin(FF)
  ifmax = numpy.argmax(FF)

  xmin=XS[ifmin]; xmax=XS[ifmax]
  ymin=YS[ifmin]; ymax=YS[ifmax]
  fmin=FF[ifmin]; fmax=FF[ifmax]

  # A,B - in the middle between these points:
  A = (xmin+xmax)/2
  B = (ymin+ymax)/2

  # E,F - slope of the line connecting first and line points
  E = (xmax-xmin)/(fmax-fmin);
  F = (ymax-ymin)/(fmax-fmin);

  # furthest point from line connecting min and max point
  # it should be near resonance:
  dist = numpy.hypot(XS - xmin - (FF-fmin)*E, YS - ymin - (FF-fmin)*F)
  ires = numpy.argmax(dist)
  F0 = FF[ires]

  # min/max freq where distance > dmax/sqrt(2),
  # this is resonance width:
  ii = dist > dist[ires]/math.sqrt(2)
  dF = numpy.max(FF[ii]) - numpy.min(FF[ii])
  if dF == 0:
    dF = abs(FF[max(0,ires-1)] - FF[min(ires+1,FF.size-1)])

  # amplitude
  C = -(YS[ires]-B);
  D =   XS[ires]-A;

  return (C,D,F0,dF,A,B,E,F)

###############################################################
# Linear oscillator - a base class for fitting and returning results

class fit_lin:
  bg_offset = 4; # offset of background parameters in pars array

  # constructor; do the fit
  def __init__(self, FF, XX, YY, DD=1,
        coord=0, cbg0=0, cbg=1, lbg=0,
        do_fit=1, fit_displ=None, fit_maxiter=10000):

    self.pars=[]; # list of free parameters
    self.errs=[]; # list of parameter uncertainties
    self.coord=coord; # use coordinates/velocities
    self.cbg0=cbg0; # use background idependent on both drive and freqeuncy (2 extra parameters)
    self.cbg=cbg;   # use background proportional to drive and independent on frequency (2 extra parameters)
    self.lbg=lbg;  # use background proportional to drive and linear in frequency (2 extra parameters)

    # scaling factors
    self.fsc = numpy.max(FF)
    self.asc = numpy.max(numpy.hypot(XX, YY))
    self.dsc = numpy.max(DD)
    #self.fsc = 1
    #self.asc = 1
    #self.dsc = 1

    self.pars = self.par_init(FF/self.fsc, XX/self.asc, YY/self.asc, DD/self.dsc)
    self.errs = [0]*len(self.pars)
    self.sc = self.par_scales(self.asc, self.dsc, self.fsc)
    if len(self.pars) != len(self.sc):   raise Exception('len(pars)!=len(sc)')
    e = 0

    if do_fit:
      res = scipy.optimize.minimize(self.minfunc, self.pars, (FF/self.fsc,XX/self.asc,YY/self.asc,DD/self.dsc),
        options={'disp': fit_displ, 'maxiter': fit_maxiter})

      # Parameter uncertainty which corresponds to res.fun
      # which is relative RMS difference between function and the model.
      # df = d2f/dx2 dx2 -> dx = dqrt(0.5*df*H^-1)
      self.errs = numpy.sqrt(0.5*res.fun*numpy.diag(res.hess_inv)).tolist()
      self.pars = res.x.tolist()
      e = res.fun/FF.size

    # convert parameters to original scale
    for i in range(len(self.pars)):
      self.pars[i] *= self.sc[i]
      self.errs[i] *= self.sc[i]

    self.fsc = 1
    self.dsc = 1
    self.asc = 1

  #### Get parameters

  # complex amplitude per unit drive
  def get_amp(self, p = None):
    if p is None: p = self.pars
    return p[0] + 1j*p[1]

  # resonance frequency
  def get_f0(self, p = None):
    if p is None: p = self.pars
    return p[2]

  # resonance width
  def get_df(self, p = None):
    if p is None: p = self.pars
    return p[3]

  # complex background independent on drive and frequency
  def get_cbg0(self, p = None):
    if p is None: p = self.pars
    n=self.bg_offset
    if self.cbg0: return p[n] + 1j*p[n+1];
    else: return 0

  # complex constant background per unit drive
  def get_cbg(self, p = None):
    if p is None: p = self.pars
    n=self.bg_offset
    if self.cbg0: n+=2
    if self.cbg:  return p[n] + 1j*p[n+1];
    else: return 0

  # complex linear background per unit drive
  def get_lbg(self, p = None):
    if p is None: p = self.pars
    n=self.bg_offset
    if self.cbg0: n+=2
    if self.cbg:  n+=2
    if self.lbg: return p[n] + 1j*p[n+1];
    else: return 0

  # Uncertainties

  # complex amplitude per unit drive
  def get_amp_e(self, e = None):
    if e is None: e = self.pars
    return e[0] + 1j*e[1]

  # resonance frequency
  def get_f0_e(self, e = None):
    if e is None: e = self.pars
    return e[2]

  # resonance width
  def get_df_e(self, e = None):
    if e is None: e = self.pars
    return e[3]

  # complex background independent on drive and frequency
  def get_cbg0_e(self, e = None):
    if e is None: e = self.pars
    n=self.bg_offset
    if self.cbg0: return e[n] + 1j*e[n+1];
    else: return 0

  # complex constant background per unit drive
  def get_cbg_e(self, e = None):
    if e is None: e = self.pars
    n=self.bg_offset
    if self.cbg0: n+=2
    if self.cbg:  return e[n] + 1j*e[n+1];
    else: return 0

  # complex linear background per unit drive
  def get_lbg_e(self, e = None):
    if e is None: e = self.pars
    n=self.bg_offset
    if self.cbg0: n+=2
    if self.cbg:  n+=2
    if self.lbg: return e[n] + 1j*e[n+1];
    else: return 0

  #### Functions

  # Background part of the function
  def func_bg(self, FF, DD=1, p=None):
    VV = numpy.zeros_like(FF, dtype=complex)
    if self.cbg0: VV += self.get_cbg0(p)
    if self.cbg: VV += DD*self.get_cbg(p)
    if self.lbg: VV += DD*self.get_lbg(p)*(FF-self.get_f0(p))
    return VV

  # Function for fitting:
  def func(self, FF, DD=1, p=None):
    AM = self.get_amp(p)
    F0 = self.get_f0(p)
    dF = self.get_df(p)

    VV = dF*F0*AM*DD/(F0**2 - FF**2 + 1j*FF*dF)
    if not self.coord: VV *= 1j*FF/F0
    VV += self.func_bg(FF, DD, p)
    return VV

  # function for minimization
  def minfunc(self, par, FF, XX, YY, DD):
    VV = self.func(FF,DD, par)
    return numpy.linalg.norm(XX + 1j*YY - VV)

  ####

  # find initial conditions for scaled data, fill pars list
  def par_init(self, FF,XX,YY,DD):
    (C,D,F0,dF,A,B,E,F) = init8(FF, XX/DD, YY/DD, self.coord)
    p = [C,D,F0,dF]
    if not self.coord:
      p[0] =  D
      p[1] = -C
    if self.cbg0: p.extend((A*numpy.min(DD), B*numpy.min(DD)))
    if self.cbg: p.extend((A,B))
    if self.lbg: p.extend((E,F))
    return p

  # parameter scaling
  def par_scales(self, asc,dsc,fsc):
    ret = [asc/dsc, asc/dsc, fsc, fsc]
    if self.cbg0: ret.extend([self.asc]*2)
    if self.cbg:  ret.extend([self.asc/self.dsc]*2)
    if self.lbg:  ret.extend([self.asc/self.dsc/self.fsc]*2)
    return ret

###############################################################
# Duffing oscillator - child of Linear oscillator class
class fit_duff(fit_lin):
  bg_offset = 5; # offset of background parameters in pars array

  def __init__(*args, **kargs): fit_lin.__init__(*args, **kargs)

  # Duffing parameter
  def get_a(self, p=None):
    if p is None: p = self.pars
    return p[4]

  def get_a_e(self, e=None):
    if e is None: e = self.errs
    return e[4]

  # Function for fitting:
  def func(self, FF, DD, p=None):
    AM = self.get_amp(p);
    F0 = self.get_f0(p)
    dF = self.get_df(p)
    a  = self.get_a(p)

    VV = numpy.zeros_like(FF, dtype=complex)
    if AM!=0:
      for i in range(FF.size):
        if isinstance(DD, (list, tuple, numpy.ndarray)): d=DD[i]
        else: d = DD
        d *= dF*F0*AM
        pp = [9/16.0*a**2,  3/2.0*(F0**2-FF[i]**2)*a, (FF[i]**2-F0**2)**2 + (FF[i]*dF)**2, -abs(d)**2]
        # find only real roots (1 or 3)
        V = numpy.sqrt(numpy.roots(pp))
        V = numpy.real(V[numpy.imag(V)==0])
        # add phase
        V = d/(F0**2 - FF[i]**2 + 1j*FF[i]*dF + 3/4.0*a*V**2)
        if i>0 and len(V)>1:
          ii=numpy.argmin(abs(V-VV[i-1]))
          VV[i] = V[ii]
        else:
          ii=numpy.argmin(abs(V))
          VV[i] = V[ii]

      if not self.coord: VV *= 1j*FF/F0

    VV += self.func_bg(FF, DD, p)
    return VV

  # parameter scaling
  def par_scales(self, asc,dsc,fsc):
    sc = fit_lin.par_scales(self,asc,dsc,fsc)
    sc.insert(4, fsc**2/asc**2)
    return sc

  # find initial conditions for scaled data, fill pars list
  def par_init(self, FF,XX,YY,DD):
    p = fit_lin.par_init(self,FF,XX,YY,DD)
    # some reasonable Duffing parameter:  a*V^2 ~ df*f0
    a = p[2]*p[3] / numpy.max(numpy.hypot(XX, YY))**2
    p.insert(4, -a)
    return p


###############################################################
# Oscillator in ballisic B-phase - child of Linear oscillator class
class fit_bphase(fit_lin):
  bg_offset = 5; # offset of background parameters in pars array

  def __init__(*args, **kargs): fit_lin.__init__(*args, **kargs)

  def get_v0(self, p=None):
    if p is None: p = self.pars
    return p[4]

  def get_v0_e(self, e=None):
    if e is None: e = self.errs
    return e[4]

  # Model function:
  def func(self, FF, DD, p=None):
    AM = self.get_amp(p);
    F0 = self.get_f0(p)
    dF = self.get_df(p)
    v0 = self.get_v0(p)

    def vfunc(x,d,f):
      dFx = dF / (1 + 0.477*(x/v0)**1.16)
      return d / (F0**2 - FF[i]**2 + 1j*FF[i]*dFx)

    VV = numpy.zeros_like(FF, dtype=complex)
    if AM!=0 and v0>0:
      x=0
      for i in range(FF.size):
        if isinstance(DD, (list, tuple, numpy.ndarray)): d=DD[i]
        else: d = DD
        d *= 1j*dF*FF[i]*AM

        # find |V|
        def zfunc(x): return abs(vfunc(x,d,FF[i])) - x
        res = scipy.optimize.root_scalar(zfunc, x0=x)
        #if not res.converged: continue
        x = res.root

        VV[i] = vfunc(x,d,FF[i])

      if self.coord: VV /= 1j*FF/F0

    VV += self.func_bg(FF, DD, p)
    return VV

  # Better function for minimization.
  # (function from the base class works too, but slower and less stable)
  def minfunc(self, par, FF, XX, YY, DD):
    VV = XX + 1j*YY - self.func_bg(FF,DD,par)
    AM = self.get_amp(par);
    F0 = self.get_f0(par)
    dF = self.get_df(par)
    v0 = self.get_v0(par)
    if self.coord:  VV*= 1j*FF/F0; #  -> vel

    dFx = dF / (1 + 0.477*(numpy.abs(VV)/v0)**1.16)

    VVc = DD*AM*dF*1j*FF / (F0**2 - FF**2 + 1j*FF*dFx)
    return numpy.linalg.norm(VV - VVc)

  # parameter scaling
  def par_scales(self, asc,dsc,fsc):
    sc = fit_lin.par_scales(self,asc,dsc,fsc)
    sc.insert(4, asc)
    return sc

  # find initial conditions for scaled data, fill pars list
  def par_init(self, FF,XX,YY,DD):
    p = fit_lin.par_init(self,FF,XX,YY,DD)
    # some reasonable Duffing parameter:  a*V^2 ~ df*f0
    v0 = numpy.max(numpy.hypot(XX, YY))
    p.insert(4, v0)
    return p


###############################################################
# Oscillator with arbitrary non-linear functions f0n(|x|) and dFn(|v|)
class fit_nonlin(fit_lin):

  def __init__(self, *args, ffunc=None, dfunc=None, **kargs):
    self.ffunc=ffunc
    self.dfunc=dfunc
    fit_lin.__init__(self, *args, **kargs)

  # Function for fitting:
  def func(self, FF, DD, p=None):
    AM = self.get_amp(p);
    F0 = self.get_f0(p)
    dF = self.get_df(p)

    # coordinate!
    def cfunc(x, d, f):
      if not self.ffunc is None: F0x = F0*self.ffunc(self.asc*x)
      else: F0x = F0
      if not self.dfunc is None: dFx = dF*self.dfunc(self.asc*x*f/F0)
      else: dFx = dF
      return d / (F0x**2 - f**2 + 1j*f*dFx)

    VV = numpy.zeros_like(FF, dtype=complex)
    if AM!=0:

      x=0
      for i in range(FF.size):
        if isinstance(DD, (list, tuple, numpy.ndarray)): d=DD[i]
        else: d = DD
        d *= dF*F0*AM

        # for finding x=|X|
        def zfunc(x): return x - abs(cfunc(x,d,FF[i]))
        res = scipy.optimize.root_scalar(zfunc, x0=x)
        if not res.converged: print('not converged')
        x = res.root
        VV[i] = cfunc(x,d,FF[i])

      if not self.coord: VV *= 1j*FF/F0

    VV += self.func_bg(FF, DD, p)
    return VV


  # function for minimization
  def minfunc(self, par, FF, XX, YY, DD):
    CC = XX + 1j*YY - self.func_bg(FF,DD, par)
    AM = self.get_amp(par);
    F0 = self.get_f0(par)
    dF = self.get_df(par)
    # -> coord
    if not self.coord:  CC*= F0/(1j*FF)
    if not self.ffunc is None: F0x = F0*self.ffunc(self.asc*numpy.abs(CC))
    else: F0x = F0
    if not self.dfunc is None: dFx = dF*self.dfunc(self.asc*numpy.abs(CC)*FF/F0)
    else: dFx = dF
    CCc = DD*AM*dF*F0 / (F0x**2 - FF**2 + 1j*FF*dFx)
    return numpy.linalg.norm(CC - CCc)


###############################################################
### f(x) -> w^2(|x|), f(v) -> delta(|v|) transformation
### x==0: F(x) = -f'(0)
### else: F(x) = - 2/x int_0^2pi func(x*cos(t))*cos(t) dt/2pi

def transform(func, x, npts=100, dfunc0=None):
  tt=numpy.linspace(0, 2*math.pi, npts)
  ct=numpy.cos(tt)

  # transform one point
  def t1(x):
    if x==0:
      if not dfunc0 is None: return -dfunc0
      else: return None
    ii = func(abs(x)*ct)*ct
    return - 2/abs(x) * numpy.trapz(ii,tt)/2/math.pi

  # transform an numpy array point by point
  if isinstance(x, (numpy.ndarray, tuple, list)):
    r = numpy.zeros_like(x, dtype=float)
    for i in range(r.size): r[i] = t1(x[i])
  else:
    r = t1(x)
  return r

def itransform(func, dfunc, x, npts=100):

  tt=numpy.linspace(0, math.pi/2, npts)
  ct=numpy.cos(tt)

  # transform one point
  def t1(x):
    cx = abs(x)*ct
    ii = cx * func(cx) + cx**2/2 * dfunc(cx)
    f = -numpy.trapz(ii,tt)
    if x>0: return f
    else: return -f

  # transform an numpy array point by point
  if isinstance(x, (numpy.ndarray, tuple, list)):
    r = numpy.zeros_like(x, dtype=float)
    for i in range(r.size): r[i] = t1(x[i])
  else:
    r = t1(x)
  return r
