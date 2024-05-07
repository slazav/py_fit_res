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
        coord=0, const_bg=1, linear_bg=0,
        do_fit=1, fit_displ=None, fit_maxiter=10000):

    self.pars=[]; # list of free parameters
    self.errs=[]; # list of parameter uncertainties
    self.coord=coord; # use coordinates/velocities
    self.const_bg=const_bg;   # use constant background (2 extra parameters)
    self.linear_bg=linear_bg; # use linear background (2 extra parameters)

    # scaling factors
    self.fsc = numpy.max(FF)
    self.asc = numpy.max(numpy.hypot(XX, YY))
    self.dsc = numpy.max(DD)
    #self.fsc = 1
    #self.asc = 1
    #self.dsc = 1


    self.do_init(FF/self.fsc, XX/self.asc, YY/self.asc, DD/self.dsc)

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
    else:
      self.errs = [0]*len(self.pars)
      e = 0

    self.do_unscale()

  # Background part of the function
  def fitfunc_bg(self, par, FF, DD=1):
    F0 = par[2]
    VV = numpy.zeros_like(FF, dtype=complex)
    n=self.bg_offset
    if self.const_bg:
      VV += par[n] + 1j*par[n+1]; n+=2
    if self.linear_bg:
      VV += (par[n] + 1j*par[n+1])*(FF-F0); n+=2
    return DD*VV

  # Function for fitting:
  def fitfunc(self, par, FF, DD=1):
    AM = par[0] + 1j*par[1];
    F0 = par[2];
    dF = par[3];

    VV = dF*F0*AM*DD/(F0**2 - FF**2 + 1j*FF*dF)
    if not self.coord: VV *= 1j*FF/F0
    VV += self.fitfunc_bg(par, FF, DD)
    return VV

  # function for minimization
  def minfunc(self, par, FF, XX, YY, DD):
    VV = self.fitfunc(par, FF,DD)
    return numpy.linalg.norm(XX + 1j*YY - VV)

  # function representing the fit result
  def func(self, FF,DD=1):
    return self.fitfunc(self.pars, FF, DD)

  # only background
  def func_bg(self, FF,DD=1):
    return self.fitfunc_bg(self.pars, FF, DD)

  # find initial conditions for scaled data, fill pars list
  def do_init(self, FF,XX,YY,DD):
    (C,D,F0,dF,A,B,E,F) = init8(FF, XX/DD, YY/DD, self.coord)
    self.pars = [C,D,F0,dF]
    if not self.coord:
      self.pars[0] =  D
      self.pars[1] = -C
    if self.const_bg:  self.pars.extend((A,B))
    if self.linear_bg: self.pars.extend((E,F))

  # convert parameters to original scale
  def do_unscale(self):
    for n in (0,1):
      self.pars[n]*=self.asc/self.dsc
      self.errs[n]*=self.asc/self.dsc
    for n in (2,3):
      self.pars[n]*=self.fsc
      self.errs[n]*=self.fsc
    n=self.bg_offset
    if self.const_bg:
      self.pars[n]*=self.asc/self.dsc
      self.errs[n]*=self.asc/self.dsc
      self.pars[n+1]*=self.asc/self.dsc
      self.errs[n+1]*=self.asc/self.dsc
      n+=2
    if self.linear_bg:
      self.pars[n]*=self.asc/self.dsc/self.fsc
      self.errs[n]*=self.asc/self.dsc/self.fsc
      self.pars[n+1]*=self.asc/self.dsc/self.fsc
      self.errs[n+1]*=self.asc/self.dsc/self.fsc
      n+=2

###############################################################
# Duffing oscillator - child of Linear oscillator class
class fit_duff(fit_lin):
  bg_offset = 5; # offset of background parameters in pars array

  def __init__(*args, **kargs): fit_lin.__init__(*args, **kargs)

  # Function for fitting:
  def fitfunc(self, par, FF, DD):
    AM = par[0] + 1j*par[1];
    F0 = par[2];
    dF = par[3];
    a  = par[4];

    VV = numpy.zeros_like(FF, dtype=complex)
    if AM!=0:
      for i in range(FF.size):
        if isinstance(DD, (list, tuple, numpy.ndarray)): d=DD[i]
        else: d = DD
        d *= dF*F0*AM
        p = [9/16.0*a**2,  3/2.0*(F0**2-FF[i]**2)*a, (FF[i]**2-F0**2)**2 + (FF[i]*dF)**2, -abs(d)**2]
        # find only real roots (1 or 3)
        V = numpy.sqrt(numpy.roots(p))
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

    VV += self.fitfunc_bg(par, FF, DD)
    return VV


  # find initial conditions for scaled data, fill pars list
  def do_init(self, FF,XX,YY,DD):
    fit_lin.do_init(self,FF,XX,YY,DD)
    # some reasonable Duffing parameter:  a*V^2 ~ df*f0
    a = self.pars[2]*self.pars[3] / numpy.max(numpy.hypot(XX, YY))**2
    self.pars.insert(4, -a)

  # convert parameters to original scale
  def do_unscale(self):
    fit_lin.do_unscale(self)
    self.pars[4]*=self.fsc**2/self.asc**2
    self.errs[4]*=self.fsc**2/self.asc**2


###############################################################
# Oscillator in ballisic B-phase - child of Linear oscillator class
class fit_bphase(fit_lin):
  bg_offset = 5; # offset of background parameters in pars array

  def __init__(*args, **kargs): fit_lin.__init__(*args, **kargs)

  # Function for fitting:
  def fitfunc(self, par, FF, DD):
    AM = par[0] + 1j*par[1];
    F0 = par[2];
    dF = par[3];
    v0  = par[4];

    VV = numpy.zeros_like(FF, dtype=complex)
    if AM!=0 and v0>0:
      VV0 = 0
      E0 = 2
      while 1:
        # magic function
        dFx = dF / (1 + 0.447*(numpy.abs(VV0)/v0)**1.16)

        # velocity response
        VV = 1j*AM*DD*dF*FF / (F0**2 - FF**2 + 1j*FF*dFx)

        VV[numpy.isnan(VV)] = 0
        VV[numpy.isinf(VV)] = 0
        E = numpy.max(abs(VV-VV0)/abs(VV))
        if E < 1e-6: break
        if E > E0: break # avoid infinite growth of V
        VV0 = VV
        E0 = E

      if self.coord: VV /= 1j*FF/F0

    VV += self.fitfunc_bg(par, FF, DD)
    return VV


  # find initial conditions for scaled data, fill pars list
  def do_init(self, FF,XX,YY,DD):
    fit_lin.do_init(self,FF,XX,YY,DD)
    # some reasonable Duffing parameter:  a*V^2 ~ df*f0
    v0 = numpy.max(numpy.hypot(XX, YY))
    self.pars.insert(4, v0)

  # convert parameters to original scale
  def do_unscale(self):
    fit_lin.do_unscale(self)
    self.pars[4]*=self.asc
    self.errs[4]*=self.asc

