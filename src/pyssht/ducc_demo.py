import numpy as np
import pyssht as ssht
from time import time


def _l2error(a, b):
    return np.sqrt(np.sum(np.abs(a - b) ** 2) / np.sum(np.abs(a) ** 2))


def test_SHT(L, Method, Reality, Spin, nthreads=1):
    if Reality:
        # Generate random flms (of real signal).
        flm = np.zeros((L * L), dtype=complex)

        # Impose reality on flms.
        for el in range(L):
            m = 0
            ind = ssht.elm2ind(el, m)
            flm[ind] = np.random.randn()
            for m in range(1, el + 1):
                ind_pm = ssht.elm2ind(el, m)
                ind_nm = ssht.elm2ind(el, -m)
                flm[ind_pm] = np.random.randn() + 1j * np.random.randn()
                flm[ind_nm] = (-1) ** m * np.conj(flm[ind_pm])
    else:
        flm = np.random.randn(L * L) + 1j * np.random.randn(L * L)

    t0 = time()
    f = ssht.inverse(
        flm,
        L,
        Reality=Reality,
        Method=Method,
        Spin=Spin,
        backend="ducc",
        nthreads=nthreads,
    )
    flm_syn = ssht.forward(
        f,
        L,
        Reality=Reality,
        Method=Method,
        Spin=Spin,
        backend="ducc",
        nthreads=nthreads,
    )
    tducc = time() - t0
    t0 = time()
    f2 = ssht.inverse(flm, L, Reality=Reality, Method=Method, Spin=Spin)
    flm_syn2 = ssht.forward(f2, L, Reality=Reality, Method=Method, Spin=Spin)
    tssht = time() - t0
    return (_l2error(f, f2) + _l2error(flm_syn, flm_syn2), tssht / tducc)


def test_rot(L, nthreads=1):
    gamma = np.pi / 5
    beta = np.pi / 7
    alpha = -np.pi / 3

    flm = np.random.randn(L * L) + 1j * np.random.randn(L * L)

    t0 = time()
    flm_rot = ssht.rotate_flms(
        flm, alpha, beta, gamma, L, backend="ducc", nthreads=nthreads
    )
    tducc = time() - t0
    t0 = time()
    flm_rot2 = ssht.rotate_flms(flm, alpha, beta, gamma, L)
    tssht = time() - t0
    return (_l2error(flm_rot, flm_rot2), tssht / tducc)


L_list = [32, 64, 128, 256, 512, 1024]
nthreads = 1

print("flm rotation tests:")
for L in L_list:
    res = test_rot(L, nthreads=nthreads)
    print("L={:4}:  L2 error={:e}, speedup factor={:f}".format(L, res[0], res[1]))
L_list = [32, 64, 128, 256, 512, 1024, 2048]
print("SHT tests:")
for Method in ["MW", "MWSS", "GL", "DH"]:
    for L in L_list:
        for Reality in [False, True]:
            for Spin in [0] if Reality else [0, 1]:
                res = test_SHT(L, Method, Reality, Spin, nthreads=nthreads)
                print(
                    "{:4}, L={:4}, Reality={:5}, Spin={:1}:  L2 error={:e}, speedup factor={:f}".format(
                        Method, L, str(Reality), Spin, res[0], res[1]
                    )
                )
