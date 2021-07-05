import numpy as np
from pytest import approx, fixture, mark, raises

import pyssht as ssht


@fixture
def rng(request) -> np.random.Generator:
    return np.random.default_rng(getattr(request.config.option, "randomly_seed", None))


@mark.parametrize("method", ["MW", "MWSS", "DH", "GL"])
def test_real_back_and_forth(rng: np.random.Generator, method, L=128):
    shape = ssht.sample_shape(L, Method=method)
    f = rng.standard_normal(shape, dtype="float64")
    flm = ssht.forward(f, L, Reality=True, Method=method)
    f = ssht.inverse(flm, L, Reality=True, Method=method)
    flm_rec = ssht.forward(f, L, Reality=True, Method=method)
    assert flm_rec == approx(flm)


def test_real_back_and_forth_mw_pole(rng: np.random.Generator, method="MW_pole", L=128):
    shape = ssht.sample_shape(L, Method=method)
    f = (
        rng.standard_normal(shape, dtype="float64"),
        rng.standard_normal(dtype="float64"),
    )
    flm = ssht.forward(f, L, Reality=True, Method=method)
    f = ssht.inverse(flm, L, Reality=True, Method=method)
    flm_rec = ssht.forward(f, L, Reality=True, Method=method)
    assert flm_rec == approx(flm)


@mark.parametrize("method", ["MW", "MWSS"])
def test_real_inverse_adjoint(rng: np.random.Generator, method, L=128):
    shape = ssht.sample_shape(L, Method=method)
    f = rng.standard_normal(shape, dtype="float64")
    flm = ssht.forward(f, L, Reality=True, Method=method)
    f = ssht.inverse(flm, L, Reality=True, Method=method)

    f_prime = rng.standard_normal(shape, dtype="float64")
    flm_prime = ssht.inverse_adjoint(f_prime, L, Reality=True, Method=method)

    assert flm_prime.conj() @ flm == approx(f_prime.flatten().conj() @ f.flatten())


@mark.parametrize("method", ["MW", "MWSS"])
def test_real_forward_adjoint(rng: np.random.Generator, method, L=128):
    shape = ssht.sample_shape(L, Method=method)
    f = rng.standard_normal(shape, dtype="float64")
    flm = ssht.forward(f, L, Reality=True, Method=method)
    f = ssht.inverse(flm, L, Reality=True, Method=method)

    f_prime = rng.standard_normal(shape, dtype="float64")
    flm_prime = ssht.forward(f_prime, L, Reality=True, Method=method)
    f_prime = ssht.forward_adjoint(flm_prime, L, Reality=True, Method=method)

    assert flm_prime.conj() @ flm == approx(f_prime.flatten().conj() @ f.flatten())


@mark.parametrize("method", ["MW", "MWSS", "DH", "GL", "MW_pole"])
@mark.parametrize("spin", range(-2, 3))
def test_back_and_forth(rng: np.random.Generator, spin, method, L=128):
    flm = rng.standard_normal((L * L, 2), dtype="float64") @ [1, 1j]
    flm[0 : spin * spin] = 0.0
    f = ssht.inverse(flm, L, Spin=spin, Method=method)
    flm_rec = ssht.forward(f, L, Spin=spin, Method=method)
    assert flm_rec == approx(flm)


@mark.parametrize("method", ["MW", "MWSS"])
@mark.parametrize("spin", range(-2, 3))
def test_inverse_adjoint(rng: np.random.Generator, spin, method, L=128):
    flm = rng.standard_normal((L * L, 2), dtype="float64") @ [1, 1j]
    flm[0 : spin * spin] = 0.0
    f = ssht.inverse(flm, L, Spin=spin, Method=method)

    shape = ssht.sample_shape(L, Method=method)
    f_prime = rng.standard_normal((*shape, 2), dtype="float64") @ [1, 1j]
    flm_prime = ssht.inverse_adjoint(f_prime, L, Spin=spin, Method=method)

    assert flm_prime.conj() @ flm == approx(f_prime.flatten().conj() @ f.flatten())


@mark.parametrize("method", ["MW", "MWSS"])
@mark.parametrize("spin", range(-2, 3))
def test_forward_adjoint(rng: np.random.Generator, spin, method, L=128):
    flm = rng.standard_normal((L * L, 2), dtype="float64") @ [1, 1j]
    flm[0 : spin * spin] = 0.0
    f = ssht.inverse(flm, L, Spin=spin, Method=method)

    flm_prime = rng.standard_normal((L * L, 2), dtype="float64") @ [1, 1j]
    flm_prime[0 : spin * spin] = 0.0
    f_prime = ssht.forward_adjoint(flm_prime, L, Spin=spin, Method=method)

    assert flm_prime.conj() @ flm == approx(f_prime.flatten().conj() @ f.flatten())


def test_wrong_array_dimension_in_forward_transform(L=64):
    f = np.zeros(np.multiply.reduce(ssht.sample_shape(L)), dtype="float64")
    with raises(ssht.ssht_input_error):
        ssht.forward(f, L, Reality=True)


def test_wrong_array_dimension_in_inverse_transform(L=64):
    flm = np.zeros((L * L, 3, 2), dtype="float64") @ [1, 1j]
    with raises(ssht.ssht_input_error):
        ssht.inverse(flm, L)


def test_forward_method_type(L=64):
    f = np.zeros(np.multiply.reduce(ssht.sample_shape(L)), dtype="float64")
    with raises(ssht.ssht_input_error):
        ssht.forward(f, L, Method="DJ", Reality=True)


def test_forward_real_spin_nonzero(L=64):
    f = np.zeros(np.multiply.reduce(ssht.sample_shape(L)), dtype="float64")
    with raises(ssht.ssht_spin_error):
        ssht.forward(f, L, Reality=True, Spin=2)


def test_inverse_method_type(L=64):
    flm = np.zeros((L * L, 3, 2), dtype="float64") @ [1, 1j]
    with raises(ssht.ssht_input_error):
        ssht.inverse(flm, L, Method="DJ")


def test_inverse_real_spin_nonzero(L=64):
    flm = np.zeros((L * L, 3, 2), dtype="float64") @ [1, 1j]
    with raises(ssht.ssht_input_error):
        ssht.inverse(flm, L, Reality=True, Spin=2)


def test_rotation(rng: np.random.Generator, method="MW", L=256):
    shape = ssht.sample_shape(L, Method=method)
    f = rng.standard_normal(shape, dtype="float64")
    f[:] = 0
    f[:, 20:60] = 1
    f[20:60] = 10
    flm = ssht.forward(f, L, Reality=True, Method=method)
    f = ssht.inverse(flm, L, Reality=True, Method=method)
    rotations = np.pi * rng.uniform(0, 1, size=3)
    rotated_flm = ssht.rotate_flms(flm, *rotations, L)
    rotated_image = ssht.rotate_image(f, rotations)
    from_rotated_flm = ssht.inverse(rotated_flm, L, Reality=True, Method=method)

    def level(a: np.ndarray, where: float, range: float = 0.2) -> np.ndarray:
        return np.logical_and(a > where - range * where, a < where + range * where)

    def ratio(a: np.ndarray, b: np.ndarray) -> float:
        return 2 * (a != b).sum() / (a.sum() + b.sum())

    assert ratio(level(rotated_image, 1), level(from_rotated_flm, 1)) < 0.1
    assert ratio(level(rotated_image, 10), level(from_rotated_flm, 10)) < 0.1
