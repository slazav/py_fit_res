%define oname fit_res

Name:         python3-module-%oname
Version:      1.1
Release:      alt1

Summary:      python library for fitting linear and non-linear resonances
Group:        Sciences/Physics
License:      GPL
Url:          https://github.com/slazav/py_fit_res
Packager:     Vladislav Zavjalov <slazav@altlinux.org>
BuildArch:    noarch

Source:       %name-%version.tar

BuildRequires(pre): rpm-build-python3


%description
Python library for fitting linear and non-linear resonances

%prep
%setup

%install
install -d %buildroot%python3_sitelibdir
cp -fR %oname.py %buildroot%python3_sitelibdir/

%files
%python3_sitelibdir/*

%changelog
* Fri May 10 2024 Vladislav Zavjalov <slazav@altlinux.org> 1.1-alt1
- changes in the interface:
  - rename background parameters: const_bg -> cbg, linear_bg -> lbg
  - add cbg0 parameter (constant, drive-independent background)
  - get_* methods for accessing fit parameters
- internal modifications:
  - rewrite parameter scaling, parameter initializing, background calculation
  - fit_nonlin: remove repeated calculations
  - fit_bphase: use scipy.optimize.root_scalar for calculating response,
    use a separate minimization function for fitting
- add spec-file for building and installing altlinux rpm package
- examples:
  - always use local library; simplify code
  - add examples 5,6,7 for fit_nonlin fitting
    (examples 6 and 7 are used for the note to be published)

* Thu May 09 2024 Vladislav Zavjalov <slazav@altlinux.org> 1.0-alt1
- v 1.0. Fitting response of a linear oscillator, Duffing oscillator,
and oscillator in ballistic B-phase. Basic documentation and examples.


