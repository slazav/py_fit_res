%define oname fit_res

Name:         python3-module-%oname
Version:      1.0
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
* Thu May 09 2024 Vladislav Zavjalov <slazav@altlinux.org> 1.0-alt1
- v 1.0. Fitting response of a linear oscillator, Duffing oscillator,
and oscillator in ballistic B-phase. Basic documentation and examples.


