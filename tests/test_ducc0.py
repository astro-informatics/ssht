import numpy as np
from pytest import approx, fixture, mark, importorskip

import pyssht as ssht

skip_ducc0 = importorskip("ducc0", minversion="0.18")


@fixture
def rng():
    return np.random.default_rng()


@fixture(params=[5, 10, 15])
def order(request):
    return request.param


@fixture(params=["MW", "MWSS", "DH", "GL"])
def method(request):
    return request.param


@fixture(params=[-2, -1, 0, 1, 2])
def spin(request):
    return request.param


@fixture
def complex_coeffs(order, rng):
    coeffs = rng.random((order * order, 2), dtype="float64")
    return coeffs @ np.array([1, 1j])


@fixture
def real_coeffs(order, complex_coeffs):
    pos_index = [ssht.elm2ind(el, m) for el in range(order) for m in range(1, el + 1)]
    neg_index = [ssht.elm2ind(el, -m) for el in range(order) for m in range(1, el + 1)]
    complex_coeffs[neg_index] = np.conjugate(complex_coeffs[pos_index])
    return complex_coeffs


@fixture
def real_image(order, method, real_coeffs):
    return ssht.inverse(real_coeffs, order, Reality=True, Method=method, Spin=0)


@fixture
def complex_image(complex_coeffs, order, method, spin):
    return ssht.inverse(complex_coeffs, order, Reality=False, Method=method, Spin=spin)


def test_real_inverse_ssht_vs_ducc0(real_coeffs, order, method, nthreads=1):
    ssht_image = ssht.inverse(real_coeffs, order, Reality=True, Method=method, Spin=0)
    ducc0_image = ssht.inverse(
        real_coeffs,
        order,
        Reality=True,
        Method=method,
        Spin=0,
        backend="ducc",
        nthreads=nthreads,
    )
    assert ssht_image == approx(ducc0_image)


def test_real_forward_ssht_vs_ducc0(real_image, order, method, nthreads=1):
    ssht_coeffs = ssht.forward(real_image, order, Reality=True, Method=method, Spin=0)
    ducc0_coeffs = ssht.forward(
        real_image,
        order,
        Reality=True,
        Method=method,
        Spin=0,
        backend="ducc",
        nthreads=nthreads,
    )
    assert ssht_coeffs == approx(ducc0_coeffs)


def test_real_inverse_adjoint_ssht_vs_ducc0(real_image, order, method, nthreads=1):
    from pyssht.exceptions import ssht_input_error

    try:
        ssht_adj_coeffs = ssht.inverse_adjoint(
            real_image, order, Reality=True, Method=method, Spin=0
        )
    except ssht_input_error:
        assert method not in ("MW", "MWSS")
        return
    ducc0_adj_coeffs = ssht.inverse_adjoint(
        real_image,
        order,
        Reality=True,
        Method=method,
        Spin=0,
        backend="ducc",
        nthreads=nthreads,
    )
    assert ssht_adj_coeffs == approx(ducc0_adj_coeffs)


def test_complex_inverse_ssht_vs_ducc0(complex_coeffs, order, method, spin, nthreads=1):
    ssht_image = ssht.inverse(
        complex_coeffs, order, Reality=False, Method=method, Spin=spin
    )
    ducc0_image = ssht.inverse(
        complex_coeffs,
        order,
        Reality=False,
        Method=method,
        Spin=spin,
        backend="ducc",
        nthreads=nthreads,
    )
    assert ssht_image == approx(ducc0_image)


def test_complex_forward_ssht_vs_ducc0(complex_image, order, method, spin, nthreads=1):
    ssht_coeffs = ssht.forward(
        complex_image, order, Reality=False, Method=method, Spin=spin
    )
    ducc0_coeffs = ssht.forward(
        complex_image,
        order,
        Reality=False,
        Method=method,
        Spin=spin,
        backend="ducc",
        nthreads=nthreads,
    )
    assert ssht_coeffs == approx(ducc0_coeffs)


def test_complex_inverse_adjoint_ssht_vs_ducc0(
    complex_image, order, method, spin, nthreads=1
):
    from pyssht.exceptions import ssht_input_error

    try:
        ssht_adj_coeffs = ssht.inverse_adjoint(
            complex_image, order, Reality=False, Method=method, Spin=spin
        )
    except ssht_input_error:
        assert method not in ("MW", "MWSS")
        return

    ducc0_adj_coeffs = ssht.inverse_adjoint(
        complex_image,
        order,
        Reality=False,
        Method=method,
        Spin=spin,
        backend="ducc",
        nthreads=nthreads,
    )
    assert ssht_adj_coeffs == approx(ducc0_adj_coeffs)


@mark.parametrize("method", ["MW", "MWSS"])
def test_real_forward_adjoint(rng: np.random.Generator, method, order):
    shape = ssht.sample_shape(order, Method=method)
    f = rng.standard_normal(shape, dtype="float64")
    flm = ssht.forward(f, order, Reality=True, Method=method)
    f = ssht.inverse(flm, order, Reality=True, Method=method)

    f_prime = rng.standard_normal(shape, dtype="float64")
    flm_prime = ssht.forward(f_prime, order, Reality=True, Method=method)
    f_prime = ssht.forward_adjoint(
        flm_prime, order, Reality=True, Method=method, backend="ducc"
    )

    assert flm_prime.conj() @ flm == approx(f_prime.flatten().conj() @ f.flatten())


@mark.parametrize("method", ["MW", "MWSS"])
def test_forward_adjoint(rng: np.random.Generator, spin, method, order):
    flm = rng.standard_normal((order * order, 2), dtype="float64") @ np.array([1, 1j])
    flm[0 : spin * spin] = 0.0
    f = ssht.inverse(flm, order, Spin=spin, Method=method, backend="ducc")

    flm_prime = rng.standard_normal((order * order, 2), dtype="float64") @ np.array(
        [1, 1j]
    )
    flm_prime[0 : spin * spin] = 0.0
    f_prime = ssht.forward_adjoint(
        flm_prime, order, Spin=spin, Method=method, backend="ducc"
    )

    assert flm_prime.conj() @ flm == approx(f_prime.flatten().conj() @ f.flatten())


def test_rot(complex_coeffs, order):
    gamma = np.pi / 5
    beta = np.pi / 7
    alpha = -np.pi / 3

    ssht_rot = ssht.rotate_flms(complex_coeffs, alpha, beta, gamma, order)
    ducc0_rot = ssht.rotate_flms(
        complex_coeffs, alpha, beta, gamma, order, backend="ducc"
    )
    assert ssht_rot == approx(ducc0_rot)
