from typing import Text
from dataclasses import dataclass


class MethodParameter:
    """Base class for method parameters."""

    def asDict(self):
        """Dict where keys are Titled."""
        from dataclasses import fields

        return {field.name.title(): getattr(self, field.name) for field in fields(self)}


@dataclass
class SSHT(MethodParameter):
    method: Text = "MW"
    spin: int = 0
    reality: bool = False

    @classmethod
    def factory(
        cls, method: Text = "MW", spin: int = 0, reality: bool = False
    ) -> MethodParameter:
        from pyssht.exceptions import ssht_input_error, ssht_spin_error

        if method.upper() not in ("MW", "MW_POLE", "MWSS", "DH", "GL"):
            raise ssht_input_error(
                "Method is not recognised, Methods are: MW, MW_pole, MWSS, DH and GL"
            )
        if spin != 0 and reality:
            raise ssht_spin_error(
                "Reality set to True and Spin is not 0. However, spin signals must be complex."
            )
        return cls(method.upper(), spin=int(spin), reality=bool(reality))


@dataclass
class Ducc(MethodParameter):
    method: Text = "MW"
    spin: int = 0
    reality: bool = False
    nthreads: int = 1

    @classmethod
    def factory(
        cls,
        method: Text = "MW",
        spin: int = 0,
        reality: bool = False,
        nthreads: int = 1,
    ) -> MethodParameter:
        from pyssht.exceptions import ssht_input_error, ssht_spin_error

        if method.upper() not in ("MW", "MWSS", "DH", "GL"):
            raise ssht_input_error(
                "Method is not recognised, Methods are: MW, MW_pole, MWSS, DH and GL"
            )
        if spin != 0 and reality:
            raise ssht_spin_error(
                "Reality set to True and Spin is not 0. However, spin signals must be complex."
            )
        return cls(
            method.upper(),
            spin=int(spin),
            reality=bool(reality),
            nthreads=int(nthreads),
        )

    def asDict(self):
        """Dict where keys are Titled."""
        result = super().asDict()
        result["nthreads"] = result.pop("Nthreads", 1)
        return result


def method(name: Text = "MW", backend: Text = "SSHT", **kwargs) -> MethodParameter:
    Class = Ducc if backend.lower().startswith("ducc") else SSHT
    return Class.factory(method=name.upper(), **kwargs)
