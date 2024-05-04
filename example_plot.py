# A function for making example plots

import numpy
import matplotlib.pyplot as plt

def plot_fsweeps(fname, fit, FF,XX,YY,DD=1, npts=200, sh=0, plot_bg=1, plot_lin=0, both_dir=0):

  (fig, (ax,ay,aa)) = plt.subplots(1,3)

  n=0
  for d in numpy.unique(DD):
    # plot data (each sweep separately with different colors)
    if isinstance(DD, (list, tuple, numpy.ndarray)): ii=DD==d
    else: ii = numpy.ones_like(FF, dtype=bool)

    fmt = 'C%d.-'%(n%10)
    ax.plot(FF[ii], sh*n+XX[ii], fmt, label="data")
    ay.plot(FF[ii], sh*n+YY[ii], fmt, label="data")
    aa.plot(FF[ii], numpy.hypot(XX[ii],YY[ii]), fmt, label="data")

    # frequency greed for calculation
    ff=numpy.linspace(min(FF), max(FF), npts)
    if both_dir: ff = numpy.append(ff, numpy.flip(ff))

    # fit result
    if fit != None:
      vv = fit.func(ff, d)
      ax.plot(ff, sh*n+numpy.real(vv), 'k-', linewidth=1, label="fit")
      ay.plot(ff, sh*n+numpy.imag(vv), 'k-', linewidth=1, label="fit")
      aa.plot(ff, sh*n+numpy.abs(vv),  'k-', linewidth=1, label="fit")

    # plot also linear Lorentzian if available
    if plot_lin and fit.func_lin:
      vv = fit.func_lin(ff, d)
      ax.plot(ff, sh*n+numpy.real(vv), 'k--', linewidth=0.7, label="lin")
      ay.plot(ff, sh*n+numpy.imag(vv), 'k--', linewidth=0.7, label="lin")
      aa.plot(ff, sh*n+numpy.abs(vv),  'k--', linewidth=0.7, label="lin")

    # plot background
    if plot_bg:
      vv = fit.func_bg(ff, d)
      ax.plot(ff, sh*n+numpy.real(vv), 'k--', linewidth=0.7, label="bg")
      ay.plot(ff, sh*n+numpy.imag(vv), 'k--', linewidth=0.7, label="bg")
      aa.plot(ff, sh*n+numpy.abs(vv),  'k--', linewidth=0.7, label="bg")

    if n==0: aa.legend()

    n+=1

  ax.set_xlabel("F")
  ay.set_xlabel("F")
  aa.set_xlabel("F")

  ax.set_title("X")
  ay.set_title("Y")
  aa.set_title("hypot(X,Y)")


  plt.gcf().set_size_inches(12, 6)
  plt.savefig(fname, dpi=100)

