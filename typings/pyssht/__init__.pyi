import numpy as np
from typing import Optional, Tuple

def forward(
    f: np.ndarray,
    L: int,
    Spin: int = 0,
    Method: str = "MW",
    Reality: bool = False,
    backend: str = "SSHT",
    **kwargs,
) -> np.ndarray:
    pass

def inverse(
    flm: np.ndarray,
    L: int,
    Spin: int = 0,
    Method: str = "MW",
    Reality: bool = False,
    backend: str = "SSHT",
    **kwargs,
) -> np.ndarray:
    pass

def elm2ind(el: int, m: int) -> int:
    pass

def ind2elm(ind: int) -> Tuple[int, int]:
    pass

def rotate_image(image, rot_list) -> np.ndarray:
    pass

def rotate_flms(
    f_lm: np.ndarray,
    alpha: float,
    beta: float,
    gamma: float,
    L: int,
    dl_array_in: np.ndarray = None,
    M_in: Optional[int] = None,
    Axisymmetric: bool = False,
    Keep_dl: bool = False,
    backend: str = "SSHT",
    **kwargs,
) -> np.ndarray:
    pass
