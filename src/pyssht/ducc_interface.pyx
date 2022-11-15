import numpy as np
cimport numpy as np

cdef import_ducc0():
    import ducc0
    major, minor, patch = ducc0.__version__.split('.')
    if int(major) < 1 and int(minor) < 15:
        raise RuntimeError("pyssht requires ducc0>=0.16")
    return ducc0

cdef Py_ssize_t _nalm(Py_ssize_t lmax, Py_ssize_t mmax):
    return ((mmax + 1) * (mmax + 2)) // 2 + (mmax + 1) * (lmax - mmax)

cdef _get_theta(Py_ssize_t L, str Method):
    cdef double pi=3.1415926535897932384626433832795028841971
    if Method == 'MW' or Method == 'MW_pole':
        return pi*(2.*np.arange(L)+1) / ( 2.0 * float(L) - 1.0 )
             
    if Method == 'MWSS':
        return pi*np.arange(L+1)/float(L)

    if Method == 'DH':
        return pi*(2*np.arange(2*L)+1.) / ( 4.0 * float(L) )
           
    if Method == 'GL':
        return import_ducc0().misc.GL_thetas(L)

cdef np.ndarray _get_lidx(Py_ssize_t L):
    res = np.arange(L)
    return res*(res+1)


cdef _extract_real_alm(flm, Py_ssize_t L):
    res = np.empty((_nalm(L-1, L-1),), dtype=np.complex128)
    cdef complex[:] myres = res
    cdef complex[:] myflm = flm
    cdef Py_ssize_t ofs=0, m, i
    cdef Py_ssize_t[:] mylidx = _get_lidx(L)
    for m in range(L):
        for i in range(m,L):
           myres[ofs-m+i] = myflm[mylidx[i]+m]
        ofs += L-m
    return res

def _build_real_flm(alm, Py_ssize_t L):
    res = np.empty((L*L), dtype=np.complex128)
    cdef Py_ssize_t ofs=0, m, i
    cdef complex[:] myres=res
    cdef complex[:] myalm=alm
    cdef Py_ssize_t[:] lidx = _get_lidx(L)
    cdef double mfac
    for m in range(L):
        mfac = (-1)**m
        for i in range(m,L):
            myres[lidx[i]+m] = myalm[ofs-m+i]
            myres[lidx[i]-m] = mfac*(myalm[ofs-m+i].real - 1j*myalm[ofs-m+i].imag)
        ofs += L-m
    return res

cdef _extract_complex_alm(flm, Py_ssize_t L, Py_ssize_t Spin):
    res = np.empty((2, _nalm(L-1, L-1),), dtype=np.complex128)
    cdef Py_ssize_t ofs=0, m, i
    cdef double mfac, sfac=(-1)**abs(Spin)
    cdef complex[:,:] myres=res
    cdef complex[:] myflm=flm
    cdef Py_ssize_t[:] lidx = _get_lidx(L)
    cdef complex fp, fm
    if Spin >= 0:
        for m in range(L):
            mfac = (-1)**m
            for i in range(m,L):
                fp = myflm[lidx[i]+m]
                fm = mfac * (myflm[lidx[i]-m].real - 1j*myflm[lidx[i]-m].imag)
                myres[0, ofs-m+i] = 0.5*(fp+fm)
                myres[1, ofs-m+i] = -0.5j*(fp-fm)
            ofs += L-m
    else:
        for m in range(L):
            mfac = (-1)**m
            for i in range(m,L):
                fp = mfac*sfac*(myflm[lidx[i]-m].real - 1j*myflm[lidx[i]-m].imag)
                fm = sfac*myflm[lidx[i]+m]
                myres[0, ofs-m+i] = 0.5*(fp+fm)
                myres[1, ofs-m+i] = -0.5j*(fp-fm)
            ofs += L-m
    return res

cdef _build_complex_flm(alm, Py_ssize_t L, Py_ssize_t Spin):
    res = np.empty((L*L), dtype=np.complex128)
    cdef Py_ssize_t ofs=0, m, i
    cdef complex fp, fm
    cdef complex[:] myres=res
    cdef complex[:,:] myalm=alm
    cdef Py_ssize_t[:] lidx = _get_lidx(L)
    cdef double mfac, sfac=(-1)**abs(Spin)
    if Spin >= 0:
        for m in range(L):
            mfac = (-1)**m
            for i in range(m,L):
                fp = myalm[0, ofs-m+i] + 1j*myalm[1, ofs-m+i]
                fm = myalm[0, ofs-m+i] - 1j*myalm[1, ofs-m+i]
                myres[lidx[i]+m] = fp
                myres[lidx[i]-m] = mfac*(fm.real - 1j*fm.imag)
            ofs += L-m
    else:
        for m in range(L):
            mfac = (-1)**m
            for i in range(m,L):
                fp = myalm[0, ofs-m+i] + 1j*myalm[1, ofs-m+i]
                fm = myalm[0, ofs-m+i] - 1j*myalm[1, ofs-m+i]
                myres[lidx[i]+m] = sfac*fm
                myres[lidx[i]-m] = sfac*mfac*(fp.real -1j *fp.imag)
            ofs += L-m
    return res


def rotate_flms(flm, alpha, beta, gamma, L, int nthreads = 1):
    ducc0 = import_ducc0()
    alm = _extract_complex_alm(flm, L, 0)
    for i in range(2):
        alm[i] = ducc0.sht.rotate_alm(
            alm[i], L-1, gamma, beta, alpha, nthreads=nthreads)
    return _build_complex_flm(alm, L, 0)


def inverse(np.ndarray flm, Py_ssize_t L, Py_ssize_t Spin, str Method, bint Reality, int nthreads = 1):
    ducc0 = import_ducc0()
    gdict = {"DH":"F1", "MW":"MW", "MWSS":"CC", "GL":"GL"}
    theta = _get_theta(L, Method)
    ntheta = theta.shape[0]
    nphi = 2*L-1
    if Method == 'MWSS':
        nphi += 1
    if Reality:
        return ducc0.sht.experimental.synthesis_2d(
            alm=_extract_real_alm(flm, L).reshape((1,-1)),
            ntheta=ntheta,
            nphi=nphi,
            lmax=L-1,
            nthreads=nthreads,
            spin=0,
            geometry=gdict[Method])[0]
    elif Spin == 0:
        alm = _extract_complex_alm(flm, L,0)
        flmr = _build_real_flm(alm[0], L)
        flmi = _build_real_flm(alm[1], L)
        return inverse(flmr, L, 0, Method, True) + 1j*inverse(flmi, L, 0, Method, True)
    else:
        tmp=ducc0.sht.experimental.synthesis_2d(
            alm=_extract_complex_alm(flm, L, Spin),
            ntheta=ntheta,
            nphi=nphi,
            lmax=L-1,
            nthreads=nthreads,
            spin=abs(Spin),
            geometry=gdict[Method])
        res = -1j*tmp[1] if Spin >=0 else 1j*tmp[1]
        res -= tmp[0]
        return res


def inverse_adjoint(f, L, Spin, Method, Reality, int nthreads = 1):
    ducc0 = import_ducc0()
    gdict = {"DH":"F1", "MW":"MW", "MWSS":"CC", "GL":"GL"}
    theta = _get_theta(L, Method)
    ntheta = theta.shape[0]
    if ntheta != f.shape[0]:
        raise RuntimeError("ntheta mismatch")
    nphi = f.shape[1]
    if Reality:
        return _build_real_flm(ducc0.sht.experimental.adjoint_synthesis_2d(
            map=f.reshape((-1,f.shape[0],f.shape[1])),
            lmax=L-1,
            nthreads=nthreads,
            spin=0,
            geometry=gdict[Method])[0], L)
    elif Spin == 0:
        flmr = inverse_adjoint(f.real, L, Spin, Method, True)
        flmi = inverse_adjoint(f.imag, L, Spin, Method, True)
        alm = np.empty((2,_nalm(L-1, L-1)), dtype=np.complex128)
        alm[0] = _extract_real_alm(flmr, L)
        alm[1] = _extract_real_alm(flmi, L)
        return _build_complex_flm(alm, L, 0)
    else:
        map = f.astype(np.complex128).view(dtype=np.float64).reshape((f.shape[0],f.shape[1],2)).transpose((2,0,1))
        if  Spin < 0:
            map[1]*=-1
        res = _build_complex_flm(ducc0.sht.experimental.adjoint_synthesis_2d(
            map=map,
            lmax=L-1,
            nthreads=nthreads,
            spin=abs(Spin),
            geometry=gdict[Method]), L, Spin)
        res *= -1
        return res


def forward(f, L, Spin, Method, Reality, int nthreads = 1):
    ducc0 = import_ducc0()
    gdict = {"DH":"F1", "MW":"MW", "MWSS":"CC", "GL":"GL"}
    theta = _get_theta(L, Method)
    ntheta = theta.shape[0]
    if ntheta != f.shape[0]:
        raise RuntimeError("ntheta mismatch")
    nphi = f.shape[1]
    if Reality:
        return _build_real_flm(ducc0.sht.experimental.analysis_2d(
            map=f.reshape((-1,f.shape[0],f.shape[1])),
            lmax=L-1,
            nthreads=nthreads,
            spin=0,
            geometry=gdict[Method])[0], L)
    elif Spin == 0:
        flmr = forward(f.real, L, Spin, Method, True)
        flmi = forward(f.imag, L, Spin, Method, True)
        alm = np.empty((2,_nalm(L-1, L-1)), dtype=np.complex128)
        alm[0] = _extract_real_alm(flmr, L)
        alm[1] = _extract_real_alm(flmi, L)
        return _build_complex_flm(alm, L, 0)
    else:
        map = f.astype(np.complex128).view(dtype=np.float64).reshape((f.shape[0],f.shape[1],2)).transpose((2,0,1))
        if Spin < 0:
            map[1]*=-1
        res = _build_complex_flm(ducc0.sht.experimental.analysis_2d(
            map=map,
            lmax=L-1,
            nthreads=nthreads,
            spin=abs(Spin),
            geometry=gdict[Method]), L, Spin)
        res *= -1
        return res


def forward_adjoint(np.ndarray flm, Py_ssize_t L, Py_ssize_t Spin, str Method, bint Reality, int nthreads = 1):
    ducc0 = import_ducc0()
    gdict = {"DH":"F1", "MW":"MW", "MWSS":"CC", "GL":"GL"}
    theta = _get_theta(L, Method)
    ntheta = theta.shape[0]
    nphi = 2*L-1
    if Method == 'MWSS':
        nphi += 1
    if Reality:
        return ducc0.sht.experimental.adjoint_analysis_2d(
            alm=_extract_real_alm(flm, L).reshape((1,-1)),
            ntheta=ntheta,
            nphi=nphi,
            lmax=L-1,
            nthreads=nthreads,
            spin=0,
            geometry=gdict[Method])[0]
    elif Spin == 0:
        alm = _extract_complex_alm(flm, L, 0)
        flmr = _build_real_flm(alm[0], L)
        flmi = _build_real_flm(alm[1], L)
        return forward_adjoint(flmr, L, 0, Method, True) + 1j*forward_adjoint(flmi, L, 0, Method, True)
    else:
        tmp=ducc0.sht.experimental.adjoint_analysis_2d(
            alm=_extract_complex_alm(flm, L, Spin),
            ntheta=ntheta,
            nphi=nphi,
            lmax=L-1,
            nthreads=nthreads,
            spin=abs(Spin),
            geometry=gdict[Method])
        res = -1j*tmp[1] if Spin >=0 else 1j*tmp[1]
        res -= tmp[0]
        return res
