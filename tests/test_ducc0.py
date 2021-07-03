import numpy as np
import pyssht as ssht
from pytest import approx, fixture, mark, importorskip

skip_ducc0 = importorskip("ducc0", minversion="0.16")


@fixture
def rng():
    return np.random.default_rng()


@mark.parametrize("L", [5, 10, 15])
@mark.parametrize("Method", ["MW", "MWSS", "GL", "DH"])
def test_real_ssht_vs_ducc0(rng: np.random.Generator, L, Method, nthreads=1):
    flm = rng.random((L * L, 2), dtype="float64") @ np.array([1, 1j])
    pos_index = [ssht.elm2ind(el, m) for el in range(L) for m in range(1, el + 1)]
    neg_index = [ssht.elm2ind(el, -m) for el in range(L) for m in range(1, el + 1)]
    flm[neg_index] = np.conjugate(flm[pos_index])

    ssht.select_ducc_backend(nthreads=nthreads)
    ssht_image = ssht.inverse(flm, L, Reality=True, Method=Method, Spin=0)
    ssht_coeffs = ssht.forward(ssht_image, L, Reality=True, Method=Method, Spin=0)
    ssht.select_ssht_backend()
    ducc0_image = ssht.inverse(flm, L, Reality=True, Method=Method, Spin=0)
    ducc0_coeffs = ssht.forward(ducc0_image, L, Reality=True, Method=Method, Spin=0)
    assert ssht_image == approx(ducc0_image)
    assert ssht_coeffs == approx(ducc0_coeffs)


@mark.parametrize("L", [5, 10, 15])
@mark.parametrize("Method", ["MW", "MWSS", "GL", "DH"])
@mark.parametrize("Spin", [0, 1, 2])
def test_complex_ssht_vs_ducc0(rng: np.random.Generator, L, Method, Spin, nthreads=1):
    flm = rng.random((L * L, 2), dtype="float64") @ np.array([1, 1j])

    ssht.select_ducc_backend(nthreads=nthreads)
    ssht_image = ssht.inverse(flm, L, Reality=False, Method=Method, Spin=Spin)
    ssht_coeffs = ssht.forward(ssht_image, L, Reality=False, Method=Method, Spin=Spin)
    ssht.select_ssht_backend()
    ducc0_image = ssht.inverse(flm, L, Reality=False, Method=Method, Spin=Spin)
    ducc0_coeffs = ssht.forward(ducc0_image, L, Reality=False, Method=Method, Spin=Spin)
    assert ssht_image == approx(ducc0_image)
    assert ssht_coeffs == approx(ducc0_coeffs)


@mark.parametrize("L", [5, 10, 15])
def test_rot(rng, L, nthreads=1):
    flm = rng.random((L * L, 2), dtype="float64") @ np.array([1, 1j])

    gamma = np.pi / 5
    beta = np.pi / 7
    alpha = -np.pi / 3

    ssht.select_ducc_backend(nthreads=nthreads)
    ssht_rot = ssht.rotate_flms(flm, alpha, beta, gamma, L)
    ssht.select_ssht_backend()
    ducc0_rot = ssht.rotate_flms(flm, alpha, beta, gamma, L)
    assert ssht_rot == approx(ducc0_rot)
