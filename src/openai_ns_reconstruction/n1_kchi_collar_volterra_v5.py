"""Directed Kchi collar propagation on ordered physical-x cells.

The caller owns all input functions.  ``chi`` and ``chi_Z`` are supplied as
cellwise interval enclosures, as are the local row-5 sources ``h5`` and
``h5_Z``.  The F/Fx pair is solved from

    F = H(Fc,Fxc) + integral K0*h5*dx - 2 integral K0*chi*F*dx,

with ``K0(x,t)=t*(1-t**2/x**2)/2``.  Ordered cell prefixes are expanded with
directed interval arithmetic; an explicit factorial majorant for omitted
Volterra words is kept in a separate symmetric error branch.  Ordinary Z1
uses ``chi_Z*F + chi*F_Z`` and ``h5_Z``.  No midpoint or reset is used.

U, K, P and Ux remain the no-potential K0 transport of the incoming state.
For physical radius ``X`` the transport distance is exactly
``d=2*log(X/a)``.  The Kchi pair always uses physical x and physical ``dx``.
``incoming_Z`` is mandatory; a Z-independent inlet must pass six explicit
zero intervals.  When using direct ``Fc``/``Fxc`` aliases, all six state
components and all six Z1 components are mandatory.  This module is standalone and
does not import or read an ancestor producer or receipt.

v5 repair: physical cell endpoints are required to be exact singleton
intervals, and Z0/Z1 state aliases are kept in separate tables.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
FIELDS = ("F", "U", "K", "P", "Fx", "Ux")
LEGACY_FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
KCHI_FIELDS = ("F", "Fx")
K0_FIELDS = ("U", "K", "P", "Ux")
PRODUCT_RULE_Z0 = ((0, 0, 1, "chi * F"),)
PRODUCT_RULE_Z1 = ((1, 0, 1, "chi_Z * F"), (0, 1, 1, "chi * F_Z"))
PRODUCT_RULE = {"Z0": PRODUCT_RULE_Z0, "Z1": PRODUCT_RULE_Z1}

__all__ = [
    "SOURCE_POWER_BASIS", "FIELDS", "LEGACY_FIELDS", "KCHI_FIELDS",
    "K0_FIELDS", "PRODUCT_RULE_Z0", "PRODUCT_RULE_Z1", "PRODUCT_RULE",
    "k0_kernel", "k0_x_kernel", "propagate", "build", "enclose",
]


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_") and hasattr(value, "a") and hasattr(value, "b")


def _raw(value: Any, upper: bool) -> Any:
    raw = getattr(value, "_mpi_", None)
    return raw[1 if upper else 0] if raw is not None else value


def _iv(c: Any, value: Any) -> Any:
    if _is_interval(value):
        return value
    if isinstance(value, Mapping):
        if "interval" in value:
            return _iv(c, value["interval"])
        for lo_key, hi_key in (("lower", "upper"), ("lo", "hi"), ("left", "right")):
            if lo_key in value and hi_key in value:
                lo, hi = _iv(c, value[lo_key]), _iv(c, value[hi_key])
                return c.make_mpf((_raw(lo, False), _raw(hi, True)))
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)) and len(value) == 2:
        lo, hi = _iv(c, value[0]), _iv(c, value[1])
        return c.make_mpf((_raw(lo, False), _raw(hi, True)))
    return c.mpf(value)


def _bounds(c: Any, value: Any) -> tuple[Any, Any]:
    value = _iv(c, value)
    lo, hi = value._mpi_
    return c.make_mpf((lo, lo)), c.make_mpf((hi, hi))


def _lo(c: Any, value: Any) -> Any:
    return _bounds(c, value)[0]


def _hi(c: Any, value: Any) -> Any:
    return _bounds(c, value)[1]


def _zero(c: Any) -> Any:
    return c.mpf(0)


def _one(c: Any) -> Any:
    return c.mpf(1)


def _symmetric(c: Any, radius: Any) -> Any:
    r = _hi(c, radius)
    if r < 0:
        raise ValueError("symmetric radius must be nonnegative")
    return c.make_mpf(((-r)._mpi_[0], r._mpi_[1]))


def _nonnegative(c: Any, value: Any, label: str) -> Any:
    value = _iv(c, value)
    lo, hi = value._mpi_
    if c.make_mpf((hi, hi)) < 0:
        raise ValueError(f"{label} has negative upper endpoint")
    if c.make_mpf((lo, lo)) < 0:
        lo = _zero(c)._mpi_[0]
    return c.make_mpf((lo, hi))


def _abs_upper(c: Any, value: Any) -> Any:
    lo, hi = _bounds(c, value)
    if lo >= 0:
        return hi
    if hi <= 0:
        return -lo
    return max(-lo, hi)


def _max_upper(c: Any, values: Sequence[Any]) -> Any:
    result = _zero(c)
    for value in values:
        upper = _hi(c, value)
        if upper > result:
            result = upper
    return result


def _sum(c: Any, values: Sequence[Any]) -> Any:
    result = _zero(c)
    for value in values:
        result = result + value
    return result


def _as_ordered(value: Any, label: str) -> list[Any]:
    if isinstance(value, Mapping):
        parsed = []
        for key, item in value.items():
            try:
                index = int(key)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{label} keys must be integer cell indices") from exc
            if index < 0:
                raise ValueError(f"{label} cell index must be nonnegative")
            parsed.append((index, item))
        parsed.sort(key=lambda pair: pair[0])
        if [index for index, _ in parsed] != list(range(len(parsed))):
            raise ValueError(f"{label} cells must be consecutive from zero")
        return [item for _, item in parsed]
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return list(value)
    raise TypeError(f"{label} must be an ordered cell sequence")


def _cell_endpoints(c: Any, raw: Any, index: int) -> tuple[Any, Any]:
    if isinstance(raw, Mapping):
        left = next((raw[key] for key in ("left", "x_left", "a", "lo", "x0", "start") if key in raw), None)
        right = next((raw[key] for key in ("right", "x_right", "b", "hi", "x1", "end") if key in raw), None)
        if left is None or right is None:
            raise ValueError(f"x_cells[{index}] needs left/right endpoints")
    elif isinstance(raw, Sequence) and not isinstance(raw, (str, bytes, bytearray)) and len(raw) == 2:
        left, right = raw
    else:
        raise TypeError(f"x_cells[{index}] must be an endpoint pair")
    left, right = _iv(c, left), _iv(c, right)
    # The propagation formulas use fixed physical endpoints and a fixed
    # causal seam.  Do not let endpoint uncertainty create a signed dx or
    # an interval-valued target whose endpoint is not a physical point.
    if _lo(c, left) != _hi(c, left) or _lo(c, right) != _hi(c, right):
        raise ValueError("physical x cell endpoints must be exact singleton intervals")
    if _lo(c, left) < 0:
        raise ValueError("physical x cells must be nonnegative and ordered")
    if _lo(c, right) <= _lo(c, left):
        raise ValueError("physical x cells must have positive width")
    return left, right


def _normalise_cells(c: Any, x_cells: Any) -> list[tuple[Any, Any]]:
    raw = _as_ordered(x_cells, "x_cells")
    if not raw:
        raise ValueError("at least one physical x cell is required")
    cells = [_cell_endpoints(c, item, i) for i, item in enumerate(raw)]
    for index in range(1, len(cells)):
        if _lo(c, cells[index - 1][1]) != _lo(c, cells[index][0]):
            raise ValueError("physical x cells must be contiguous and ordered")
    return cells


def _call_cell_provider(c: Any, provider: Any, cell: tuple[Any, Any], index: int) -> Any:
    if not callable(provider):
        return provider
    try:
        return provider(cell[0], cell[1], index)
    except TypeError:
        try:
            return provider(cell, index)
        except TypeError:
            try:
                return provider(cell[0], cell[1])
            except TypeError:
                whole_cell = c.make_mpf((_raw(cell[0], False), _raw(cell[1], True)))
                return provider(whole_cell)


def _normalise_cell_data(c: Any, raw: Any, cells: Sequence[tuple[Any, Any]], label: str) -> list[Any]:
    if raw is None:
        raise TypeError(f"{label} is required")
    if callable(raw):
        return [_iv(c, _call_cell_provider(c, raw, cell, i)) for i, cell in enumerate(cells)]
    if not isinstance(raw, (Mapping, Sequence)) or _is_interval(raw):
        return [_iv(c, raw) for _ in cells]
    values = _as_ordered(raw, label)
    if len(values) != len(cells):
        raise ValueError(f"{label} must have one interval for each physical x cell")
    return [_iv(c, value) for value in values]


def _lookup(mapping: Mapping[Any, Any], names: Sequence[str], label: str) -> Any:
    for name in names:
        if name in mapping:
            return mapping[name]
    lowered = {str(key).lower(): value for key, value in mapping.items()}
    for name in names:
        if name.lower() in lowered:
            return lowered[name.lower()]
    raise KeyError(label)


Z0_STATE_ALIASES = {
    "F": ("F", "F1", "Fc", "Fc1"),
    "U": ("U", "Uz1", "Uc"),
    "K": ("K", "K1", "Kc"),
    "P": ("P", "P1", "Pc"),
    "Fx": ("Fx", "F_x", "d_x_F", "d_x_F1", "Fxc"),
    "Ux": ("Ux", "U_x", "d_x_U", "d_x_Uz1", "Uxc"),
}
Z1_STATE_ALIASES = {
    # Canonical F/U/K/P/Fx/Ux names are row-local for an explicit
    # incoming_Z mapping.  The suffixed and legacy names are Z1-specific;
    # no Z1 component is filled from an omitted Z0 component.
    "F": ("F", "F_Z", "Fc_Z", "incoming_Fc_Z", "F1", "F1_Z"),
    "U": ("U", "U_Z", "Uc_Z", "incoming_Uc_Z", "Uz1", "Uz1_Z"),
    "K": ("K", "K_Z", "Kc_Z", "incoming_Kc_Z", "K1", "K1_Z"),
    "P": ("P", "P_Z", "Pc_Z", "incoming_Pc_Z", "P1", "P1_Z"),
    "Fx": ("Fx", "Fx_Z", "F_x_Z", "d_x_F_Z", "d_x_F1_Z",
           "Fxc_Z", "incoming_Fxc_Z", "d_x_F1"),
    "Ux": ("Ux", "Ux_Z", "U_x_Z", "d_x_U_Z", "d_x_Uz1_Z",
           "Uxc_Z", "incoming_Uxc_Z", "d_x_Uz1"),
}


def _state_values(c: Any, raw: Any, label: str, *, z1: bool = False) -> dict[str, Any]:
    if isinstance(raw, Mapping):
        aliases = Z1_STATE_ALIASES if z1 else Z0_STATE_ALIASES
        result = {}
        for field, names in aliases.items():
            try:
                result[field] = _iv(c, _lookup(raw, names, f"{label}.{field}"))
            except KeyError:
                raise ValueError(f"{label} missing {field}") from None
        return result
    if isinstance(raw, Sequence) and not isinstance(raw, (str, bytes, bytearray)) and len(raw) == 6:
        return {field: _iv(c, value) for field, value in zip(FIELDS, raw)}
    raise TypeError(f"{label} must be a six-state mapping or sequence")


def _incoming_pair(c: Any, incoming: Any, incoming_Z: Any) -> tuple[dict[str, Any], dict[str, Any]]:
    values = _state_values(c, incoming, "incoming", z1=False)
    if incoming_Z is None:
        raise TypeError("incoming_Z is required; pass explicit zero intervals for a Z-independent inlet")
    return values, _state_values(c, incoming_Z, "incoming_Z", z1=True)


def k0_kernel(c: Any, x: Any, t: Any) -> Any:
    """Physical kernel ``K0(x,t)=t*(1-t**2/x**2)/2``."""
    x, t = _iv(c, x), _iv(c, t)
    if _lo(c, x) <= 0:
        raise ValueError("K0 target x must be strictly positive")
    return t * (1 - (t / x) ** 2) / 2


def k0_x_kernel(c: Any, x: Any, t: Any) -> Any:
    """Physical derivative ``d_x K0=t**3/x**3``."""
    x, t = _iv(c, x), _iv(c, t)
    if _lo(c, x) <= 0:
        raise ValueError("d_x K0 target x must be strictly positive")
    return t ** 3 / x ** 3


def _kernel_bound(c: Any, target: tuple[Any, Any], source: tuple[Any, Any], derivative: bool) -> Any:
    x = c.make_mpf((_raw(target[0], False), _raw(target[1], True)))
    t = c.make_mpf((_raw(source[0], False), _raw(source[1], True)))
    value = k0_x_kernel(c, x, t) if derivative else k0_kernel(c, x, t)
    upper = _abs_upper(c, value)
    return c.make_mpf((_zero(c)._mpi_[0], upper._mpi_[1]))


def _mass(c: Any, target_x: Any, source: tuple[Any, Any], derivative: bool) -> Any:
    x = _iv(c, target_x)
    a, b = source
    if _lo(c, x) <= 0:
        raise ValueError("positive target x required for physical Kchi mass")
    value = ((b ** 4 - a ** 4) / (4 * x ** 3) if derivative else
             (b ** 2 - a ** 2) / 4 - (b ** 4 - a ** 4) / (8 * x ** 2))
    return _nonnegative(c, value, "Kchi cell mass")


def _homogeneous(c: Any, seam: Any, target: Any, value: Any, derivative: Any,
                 z_value: Any, z_derivative: Any) -> tuple[Any, Any, Any, Any]:
    a, x = _iv(c, seam), _iv(c, target)
    if _lo(c, a) <= 0 or _lo(c, x) <= 0:
        raise ValueError("the Kchi collar seam and target must be strictly positive")
    ratio = a / x
    h = value + derivative * a * (1 - ratio ** 2) / 2
    hx = derivative * ratio ** 3
    hz = z_value + z_derivative * a * (1 - ratio ** 2) / 2
    hzx = z_derivative * ratio ** 3
    return h, hx, hz, hzx


def _prefix_base(c: Any, cells: Sequence[tuple[Any, Any]], h5: Sequence[Any],
                 h5_z: Sequence[Any], incoming: Mapping[str, Any],
                 incoming_z: Mapping[str, Any], branch: str) -> dict[str, Any]:
    """Build profile and endpoint boxes for the zeroth Volterra word."""

    seam = cells[0][0]
    profile = {field: [] for field in KCHI_FIELDS}
    profile_z = {field: [] for field in KCHI_FIELDS}
    endpoint = {field: [] for field in KCHI_FIELDS}
    endpoint_z = {field: [] for field in KCHI_FIELDS}
    for i, cell in enumerate(cells):
        target_profile = c.make_mpf((_raw(cell[0], False), _raw(cell[1], True)))
        target_endpoint = cell[1]
        h, hx, hz, hzx = _homogeneous(
            c, seam, target_profile, incoming["F"], incoming["Fx"],
            incoming_z["F"], incoming_z["Fx"],
        )
        he, hxe, hze, hzee = _homogeneous(
            c, seam, target_endpoint, incoming["F"], incoming["Fx"],
            incoming_z["F"], incoming_z["Fx"],
        )
        if branch == "homogeneous":
            profile["F"].append(h)
            profile["Fx"].append(hx)
            profile_z["F"].append(hz)
            profile_z["Fx"].append(hzx)
            endpoint["F"].append(he)
            endpoint["Fx"].append(hxe)
            endpoint_z["F"].append(hze)
            endpoint_z["Fx"].append(hzee)
            continue

        # The forcing branch starts with zero inlet and receives h5 once.
        # ``width`` below is physical dx, not a normalized cell measure.
        source_f = _zero(c)
        source_fx = _zero(c)
        source_fz = _zero(c)
        source_fxz = _zero(c)
        for j, source_cell in enumerate(cells[: i + 1]):
            width = source_cell[1] - source_cell[0]
            k = _kernel_bound(c, cell, source_cell, False)
            kx = _kernel_bound(c, cell, source_cell, True)
            source_f = source_f + k * h5[j] * width
            source_fx = source_fx + kx * h5[j] * width
            source_fz = source_fz + k * h5_z[j] * width
            source_fxz = source_fxz + kx * h5_z[j] * width
        endpoint_f = _zero(c)
        endpoint_fx = _zero(c)
        endpoint_fz = _zero(c)
        endpoint_fxz = _zero(c)
        for j, source_cell in enumerate(cells[: i + 1]):
            endpoint_f = endpoint_f + _mass(c, target_endpoint, source_cell, False) * h5[j]
            endpoint_fx = endpoint_fx + _mass(c, target_endpoint, source_cell, True) * h5[j]
            endpoint_fz = endpoint_fz + _mass(c, target_endpoint, source_cell, False) * h5_z[j]
            endpoint_fxz = endpoint_fxz + _mass(c, target_endpoint, source_cell, True) * h5_z[j]
        profile["F"].append(source_f)
        profile["Fx"].append(source_fx)
        profile_z["F"].append(source_fz)
        profile_z["Fx"].append(source_fxz)
        endpoint["F"].append(endpoint_f)
        endpoint["Fx"].append(endpoint_fx)
        endpoint_z["F"].append(endpoint_fz)
        endpoint_z["Fx"].append(endpoint_fxz)
    return {"profile": profile, "profile_Z": profile_z,
            "endpoint": endpoint, "endpoint_Z": endpoint_z}


def _word(c: Any, cells: Sequence[tuple[Any, Any]], chi: Sequence[Any],
          chi_z: Sequence[Any], previous: Mapping[str, Sequence[Any]],
          previous_z: Mapping[str, Sequence[Any]]) -> dict[str, Any]:
    """Apply one ordered Volterra word to interval profiles and endpoints."""

    out_profile = {field: [] for field in KCHI_FIELDS}
    out_profile_z = {field: [] for field in KCHI_FIELDS}
    out_endpoint = {field: [] for field in KCHI_FIELDS}
    out_endpoint_z = {field: [] for field in KCHI_FIELDS}
    for i, target_cell in enumerate(cells):
        profile_f = _zero(c)
        profile_fx = _zero(c)
        profile_fz = _zero(c)
        profile_fxz = _zero(c)
        endpoint_f = _zero(c)
        endpoint_fx = _zero(c)
        endpoint_fz = _zero(c)
        endpoint_fxz = _zero(c)
        for j, source_cell in enumerate(cells[: i + 1]):
            width = source_cell[1] - source_cell[0]
            k = _kernel_bound(c, target_cell, source_cell, False)
            kx = _kernel_bound(c, target_cell, source_cell, True)
            value_z = chi[j] * previous_z["F"][j] + chi_z[j] * previous["F"][j]
            profile_f = profile_f + (-2) * k * chi[j] * previous["F"][j] * width
            # d_x acts on the outer kernel.  It therefore consumes the
            # previous F word, not the previous Fx word.
            profile_fx = profile_fx + (-2) * kx * chi[j] * previous["F"][j] * width
            profile_fz = profile_fz + (-2) * k * value_z * width
            profile_fxz = profile_fxz + (-2) * kx * value_z * width
            mass = _mass(c, target_cell[1], source_cell, False)
            mass_x = _mass(c, target_cell[1], source_cell, True)
            endpoint_f = endpoint_f + (-2) * mass * chi[j] * previous["F"][j]
            endpoint_fx = endpoint_fx + (-2) * mass_x * chi[j] * previous["F"][j]
            endpoint_fz = endpoint_fz + (-2) * mass * value_z
            endpoint_fxz = endpoint_fxz + (-2) * mass_x * value_z
        out_profile["F"].append(profile_f)
        out_profile["Fx"].append(profile_fx)
        out_profile_z["F"].append(profile_fz)
        out_profile_z["Fx"].append(profile_fxz)
        out_endpoint["F"].append(endpoint_f)
        out_endpoint["Fx"].append(endpoint_fx)
        out_endpoint_z["F"].append(endpoint_fz)
        out_endpoint_z["Fx"].append(endpoint_fxz)
    return {"profile": out_profile, "profile_Z": out_profile_z,
            "endpoint": out_endpoint, "endpoint_Z": out_endpoint_z}


def _sum_words(c: Any, words: Sequence[Mapping[str, Any]], key: str,
               field: str, index: int) -> Any:
    return _sum(c, [word[key][field][index] for word in words])


def _factorial(c: Any, n: int) -> Any:
    result = _one(c)
    for index in range(2, n + 1):
        result = result * index
    return result


def _tail_bound(c: Any, cells: Sequence[tuple[Any, Any]], chi: Sequence[Any],
                chi_z: Sequence[Any], words: Sequence[Mapping[str, Any]],
                max_order: int, branch: str) -> dict[str, Any]:
    """Explicit factorial remainder for one branch and both derivative rows."""

    total_width = _sum(c, [cell[1] - cell[0] for cell in cells])
    # The Volterra operator only sees ordered source cells j <= i.  Excluding
    # noncausal t>x pairs is essential: their independent boxes can create a
    # huge artificial negative K0 interval and an invalid tail majorant.
    k_values = [_kernel_bound(c, target, source, False)
                for i, target in enumerate(cells)
                for source in cells[: i + 1]]
    kx_values = [_kernel_bound(c, target, source, True)
                 for i, target in enumerate(cells)
                 for source in cells[: i + 1]]
    kmax = _max_upper(c, k_values)
    kxmax = _max_upper(c, kx_values)
    chimax = _max_upper(c, [_abs_upper(c, value) for value in chi])
    chizmax = _max_upper(c, [_abs_upper(c, value) for value in chi_z])
    alpha = 2 * chimax * kmax * total_width
    alpha_x = 2 * chimax * kxmax * total_width
    alpha_z = 2 * chizmax * kmax * total_width
    alpha_xz = 2 * chizmax * kxmax * total_width
    base_f = _max_upper(c, [_abs_upper(c, value) for value in words[0]["profile"]["F"]])
    base_fz = _max_upper(c, [_abs_upper(c, value) for value in words[0]["profile_Z"]["F"]])

    def zero_rate(rate: Any) -> bool:
        return _lo(c, rate) == 0 and _hi(c, rate) == 0

    def value_tail(base: Any, rate: Any) -> Any:
        """Tail of repeated K0 words with rate ``rate``."""

        if zero_rate(rate):
            return _zero(c)
        return (base * c.exp(rate) * rate ** (max_order + 1)
                / _factorial(c, max_order + 1))

    def z_insertion_tail(base: Any, rate: Any, insertion: Any) -> Any:
        """Tail with one chi_Z insertion among repeated K0 words."""

        if zero_rate(rate):
            # For N=0 the first omitted word is exactly A_Z B.  Returning
            # zero here would silently lose chi_Z when chi itself is zero.
            return base * insertion if max_order == 0 else _zero(c)
        return (base * insertion * c.exp(rate) * rate ** max_order
                / _factorial(c, max_order))

    def outer_x_tail(base: Any, inner_rate: Any, outer_rate: Any) -> Any:
        """Tail with one outer d_x K0 and K0 words inside it."""

        if zero_rate(outer_rate):
            return _zero(c)
        # The outer derivative is one factor alpha_x; all inner words use
        # alpha.  alpha_x is not used as a repeated rate.
        return (base * outer_rate * c.exp(inner_rate) * inner_rate ** max_order
                / _factorial(c, max_order + 1))

    def inner_z_x_tail(base: Any, inner_rate: Any, outer_rate: Any,
                       insertion: Any, outer_insertion: Any) -> Any:
        """A conservative Fx_Z tail, split by Z insertion location."""

        # Z insertion in the outer differentiated kernel: one alpha_x_Z,
        # with the remaining K0 words controlled by alpha.
        outer = z_insertion_tail(base, inner_rate, outer_insertion)
        if zero_rate(outer_rate):
            inner = _zero(c)
        elif max_order == 0:
            # This larger bound covers all n>=1 inner-insertion terms; the
            # n=1 term actually has no inner insertion.
            inner = base * outer_rate * insertion * c.exp(inner_rate)
        else:
            # (n-1)/n! <= 1/((N+1)(n-2)!), n>=N+1.
            inner = (base * outer_rate * insertion * c.exp(inner_rate)
                     * inner_rate ** (max_order - 1)
                     / ((max_order + 1) * _factorial(c, max_order - 1)))
        return outer + inner

    f_tail = value_tail(base_f, alpha)
    fz_base = value_tail(base_fz, alpha)
    fz_extra = z_insertion_tail(base_f, alpha, alpha_z)
    # d_x differentiates only the outer kernel, so alpha controls every
    # inner word even when d_x K0 has a larger or smaller supremum.
    fx_tail = outer_x_tail(base_f, alpha, alpha_x)
    fxz_base = outer_x_tail(base_fz, alpha, alpha_x)
    fxz_extra = inner_z_x_tail(base_f, alpha, alpha_x, alpha_z, alpha_xz)
    fz_tail = fz_base + fz_extra
    fxz_tail = fxz_base + fxz_extra
    return {
        "F": _symmetric(c, f_tail), "Fx": _symmetric(c, fx_tail),
        "F_Z": _symmetric(c, fz_tail), "Fx_Z": _symmetric(c, fxz_tail),
        "radius": {"F": f_tail, "Fx": fx_tail,
                   "F_Z": fz_tail, "Fx_Z": fxz_tail},
        "factorial": {
            "max_order": max_order, "word_tail_starts_at": max_order + 1,
            "formula": "B*exp(alpha)*alpha^(N+1)/(N+1)!",
            "Fx_formula": "B*alpha_x*exp(alpha)*alpha^N/(N+1)!",
            "Z_formula": "B_Z*exp(alpha)*alpha^(N+1)/(N+1)! + B*alpha_Z*exp(alpha)*alpha^N/N!",
            "Fx_Z_formula": "B_Z*alpha_x*exp(alpha)*alpha^N/(N+1)! + outer_alpha_x_Z + inner_alpha_x_alpha_Z",
            "alpha": alpha, "alpha_x": alpha_x,
            "alpha_Z": alpha_z, "alpha_x_Z": alpha_xz,
            "chi_sup": chimax, "chi_Z_sup": chizmax,
            "K0_sup": kmax, "d_x_K0_sup": kxmax,
            "physical_width": total_width,
        },
        "branch": branch, "symmetric_error": True, "error_only": True,
    }


def _branch(c: Any, cells: Sequence[tuple[Any, Any]], chi: Sequence[Any],
            chi_z: Sequence[Any], h5: Sequence[Any], h5_z: Sequence[Any],
            incoming: Mapping[str, Any], incoming_z: Mapping[str, Any],
            *, branch: str, max_order: int) -> dict[str, Any]:
    base = _prefix_base(c, cells, h5, h5_z, incoming, incoming_z, branch)
    words: list[dict[str, Any]] = [base]
    for _ in range(max_order):
        previous = words[-1]
        words.append(_word(c, cells, chi, chi_z,
                           previous["profile"], previous["profile_Z"]))
    endpoints = {}
    for index, cell in enumerate(cells):
        endpoints[str(index)] = {
            "x_left": cell[0], "x_right": cell[1],
            "F": _sum_words(c, words, "endpoint", "F", index),
            "Fx": _sum_words(c, words, "endpoint", "Fx", index),
            "F_Z": _sum_words(c, words, "endpoint_Z", "F", index),
            "Fx_Z": _sum_words(c, words, "endpoint_Z", "Fx", index),
        }
    tail = _tail_bound(c, cells, chi, chi_z, words, max_order, branch)
    final = endpoints[str(len(cells) - 1)]
    enclosure = {field: final[field] + tail[field]
                 for field in ("F", "Fx", "F_Z", "Fx_Z")}
    return {
        "branch": branch, "words": words, "word_count": len(words),
        "max_order": max_order, "prefix_endpoint_enclosures": endpoints,
        "finite_prefix": {field: final[field] for field in ("F", "Fx", "F_Z", "Fx_Z")},
        "tail": tail, "enclosure": enclosure,
        "symmetric_tail_separate": True, "homogeneous_source_separate": True,
    }


def _transport_k0(c: Any, cells: Sequence[tuple[Any, Any]], incoming: Mapping[str, Any],
                 incoming_z: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    """Transport U/K/P/Ux by the existing K0 no-reset formulas."""

    # X=a*exp(d/2) is the physical collar map; do not replace this logarithmic
    # distance by hb*(X-a), which belongs to a normalized phase coordinate.
    d = 2 * c.log(cells[-1][1] / cells[0][0])
    e = c.exp(-d)
    eh = c.exp(-d / 2)

    def row(values: Mapping[str, Any]) -> dict[str, Any]:
        return {
            "U": values["U"] + d * values["Ux"],
            "K": e * values["K"] - (1 - e) * values["Ux"],
            "P": values["P"],
            "Ux": eh * values["Ux"],
        }

    return row(incoming), row(incoming_z)


def _parse_args(args: tuple[Any, ...], kwargs: dict[str, Any]) -> dict[str, Any]:
    if args:
        names = ("x_cells", "chi", "chi_Z", "incoming", "incoming_Z", "h5", "h5_Z")
        if len(args) > len(names):
            raise TypeError("at most seven positional propagation arguments are accepted")
        for name, value in zip(names, args):
            if name in kwargs:
                raise TypeError(f"{name} supplied both positionally and by keyword")
            kwargs[name] = value
    aliases = {
        "cells": "x_cells", "physical_x_cells": "x_cells",
        "chi_cells": "chi", "chi_z": "chi_Z", "chiZ": "chi_Z",
        "incoming_state": "incoming", "incoming_state_Z": "incoming_Z",
        "h5_cells": "h5", "h5_z": "h5_Z", "h5Z": "h5_Z",
        "N": "max_order", "picard_order": "max_order",
        "transport_hb": "hb", "phase_scale": "hb",
    }
    for old, new in aliases.items():
        if old in kwargs:
            if new in kwargs:
                raise TypeError(f"{old} and {new} supplied together")
            kwargs[new] = kwargs.pop(old)
    return kwargs


def propagate(c: Any, *args: Any, **kwargs: Any) -> dict[str, Any]:
    """Return a directed enclosure for one physical Kchi collar prefix.

    Required data are ``x_cells``, ``chi``, ``chi_Z``, ``incoming``,
    ``incoming_Z``, ``h5`` and ``h5_Z``.  A Z-independent inlet must still
    pass six explicit zero intervals in ``incoming_Z``.  Cell data are whole-cell intervals or
    callbacks evaluated on interval endpoints; no midpoint is constructed.
    ``max_order`` is the finite Neumann prefix order (default 8).  The
    omitted words are returned only through each branch's symmetric tail.
    """

    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf"):
        raise TypeError("c must be the caller's MPIntervalContext")
    options = _parse_args(args, dict(kwargs))
    x_cells = options.pop("x_cells", None)
    if x_cells is None:
        raise TypeError("x_cells are required")
    cells = _normalise_cells(c, x_cells)
    chi_raw = options.pop("chi", None)
    chi_z_raw = options.pop("chi_Z", None)
    h5_raw = options.pop("h5", None)
    h5_z_raw = options.pop("h5_Z", None)
    if chi_raw is None or chi_z_raw is None or h5_raw is None or h5_z_raw is None:
        raise TypeError("chi, chi_Z, h5 and h5_Z are required caller-owned cell data")
    chi = _normalise_cell_data(c, chi_raw, cells, "chi")
    chi_z = _normalise_cell_data(c, chi_z_raw, cells, "chi_Z")
    h5 = _normalise_cell_data(c, h5_raw, cells, "h5")
    h5_z = _normalise_cell_data(c, h5_z_raw, cells, "h5_Z")
    incoming_raw = options.pop("incoming", None)
    incoming_z_raw = options.pop("incoming_Z", None)
    # Also accept the natural caller spelling Fc/Fxc and their ordinary Z1
    # rows.  Direct aliases are all-or-nothing: silently zeroing an omitted
    # U/K/P/Ux history would erase caller-owned data.
    if incoming_raw is None:
        direct_names = {
            "F": ("Fc", "incoming_Fc", "F", "F1"),
            "U": ("Uc", "incoming_Uc", "U", "Uz1"),
            "K": ("Kc", "incoming_Kc", "K", "K1"),
            "P": ("Pc", "incoming_Pc", "P", "P1"),
            "Fx": ("Fxc", "incoming_Fxc", "Fx", "d_x_F", "d_x_F1"),
            "Ux": ("Uxc", "incoming_Uxc", "Ux", "d_x_U", "d_x_Uz1"),
        }
        direct = {}
        for field, names in direct_names.items():
            for name in names:
                if name in options:
                    direct[field] = options.pop(name)
                    break
        if direct:
            missing = [field for field in FIELDS if field not in direct]
            if missing:
                raise ValueError("direct incoming spelling requires all six states; missing "
                                 + ", ".join(missing))
            incoming_raw = direct
    if incoming_z_raw is None:
        direct_z_names = {
            "F": ("Fc_Z", "incoming_Fc_Z", "F_Z"),
            "U": ("Uc_Z", "incoming_Uc_Z", "U_Z"),
            "K": ("Kc_Z", "incoming_Kc_Z", "K_Z"),
            "P": ("Pc_Z", "incoming_Pc_Z", "P_Z"),
            "Fx": ("Fxc_Z", "incoming_Fxc_Z", "Fx_Z"),
            "Ux": ("Uxc_Z", "incoming_Uxc_Z", "Ux_Z"),
        }
        direct_z = {}
        for field, names in direct_z_names.items():
            for name in names:
                if name in options:
                    direct_z[field] = options.pop(name)
                    break
        if direct_z:
            missing_z = [field for field in FIELDS if field not in direct_z]
            if missing_z:
                raise ValueError("direct incoming_Z spelling requires all six states; missing "
                                 + ", ".join(missing_z))
            incoming_z_raw = direct_z
    if incoming_raw is None:
        raise TypeError("incoming six-state history is required")
    incoming, incoming_z = _incoming_pair(c, incoming_raw, incoming_z_raw)
    max_order = options.pop("max_order", 8)
    if isinstance(max_order, bool) or not isinstance(max_order, int) or max_order < 0:
        raise ValueError("max_order must be a nonnegative integer")
    if "hb" in options:
        raise ValueError("hb is a normalized-phase scale; physical K0 transport uses d=2*log(X/a)")
    if options:
        raise TypeError("unexpected propagation options: " + ", ".join(sorted(options)))
    if _lo(c, cells[0][0]) <= 0:
        raise ValueError("Kchi collar seam x must be strictly positive")
    if _lo(c, cells[0][0]) != 2 or _hi(c, cells[0][0]) != 2:
        raise ValueError("physical K0 collar transport is certified only at seam x=2")

    zero_state = {field: _zero(c) for field in FIELDS}
    zero_z_state = {field: _zero(c) for field in FIELDS}
    homogeneous = _branch(
        c, cells, chi, chi_z,
        [_zero(c) for _ in cells], [_zero(c) for _ in cells],
        incoming, incoming_z, branch="homogeneous", max_order=max_order,
    )
    forcing = _branch(
        c, cells, chi, chi_z, h5, h5_z,
        zero_state, zero_z_state, branch="forcing", max_order=max_order,
    )
    k0_end, k0_end_z = _transport_k0(c, cells, incoming, incoming_z)
    last_h = homogeneous["enclosure"]
    last_f = forcing["enclosure"]
    state = {
        "F": last_h["F"] + last_f["F"], "U": k0_end["U"],
        "K": k0_end["K"], "P": k0_end["P"],
        "Fx": last_h["Fx"] + last_f["Fx"], "Ux": k0_end["Ux"],
    }
    state_z = {
        "F": last_h["F_Z"] + last_f["F_Z"], "U": k0_end_z["U"],
        "K": k0_end_z["K"], "P": k0_end_z["P"],
        "Fx": last_h["Fx_Z"] + last_f["Fx_Z"], "Ux": k0_end_z["Ux"],
    }
    legacy = {"F1": state["F"], "Uz1": state["U"], "K1": state["K"],
              "P1": state["P"], "d_x_F1": state["Fx"], "d_x_Uz1": state["Ux"]}
    legacy_z = {"F1": state_z["F"], "Uz1": state_z["U"], "K1": state_z["K"],
                "P1": state_z["P"], "d_x_F1": state_z["Fx"], "d_x_Uz1": state_z["Ux"]}
    return {
        "schema": "n1_kchi_collar_volterra_v5",
        "cells": [{"x_left": cell[0], "x_right": cell[1], "dx": cell[1] - cell[0]}
                  for cell in cells],
        "x_cells": len(cells), "seam_x": cells[0][0], "exit_x": cells[-1][1],
        "chi_cells": chi, "chi_Z_cells": chi_z,
        "h5_cells": h5, "h5_Z_cells": h5_z,
        "incoming": incoming, "incoming_Z": incoming_z,
        "homogeneous": homogeneous, "forcing": forcing,
        "homogeneous_branch": homogeneous, "forcing_branch": forcing,
        "k0_transport": {"state": k0_end, "state_Z": k0_end_z,
                          "physical_distance": 2 * c.log(cells[-1][1] / cells[0][0]),
                          "reset": False, "F_and_Fx_excluded": True},
        "state": state, "state_Z": state_z,
        "exit_state": state, "exit_state_Z": state_z,
        "legacy_state": legacy, "legacy_state_Z": legacy_z,
        "ordinary_Z_orders": [0, 1], "ordinary_product_rule": PRODUCT_RULE,
        "source_measure": "physical dx = x_right - x_left exactly once",
        "kernel": "Kchi = sum_n (-2)^n K0 (chi K0)^n",
        "equation": "L0 F + 2 chi(x,Z) F = h5(x,Z), L0=d_xx+3*d_x/x",
        "ordered_volterra_prefix": True, "factorial_tail": True,
        "directed_interval_arithmetic": True, "midpoint_values_used": False,
        "homogeneous_and_forcing_separate": True, "core_history_reset": False,
        "chi_variable_across_cells": True, "chi_Z_included": True,
        "h5_Z_included": True, "allword_tail_added": False,
        "conditions": {
            "physical_x_cells_ordered": True, "physical_x_cells_contiguous": True,
            "seam_x_positive": True, "caller_owned_chi_and_chi_Z": True,
            "caller_owned_h5_and_h5_Z": True,
            "incoming_Fc_Fxc_nonzero_allowed": True,
            "incoming_Z_required": True,
            "direct_aliases_require_all_six_components": True,
            "transport_physical_log_ratio": True,
        },
        "limits": [
            "Cellwise chi and h5 boxes must enclose the corresponding functions on the whole physical cell.",
            "The factorial tail is a separate symmetric enclosure; it is not a signed field value.",
            "No complete all-word feedback, gauge compensation, or terminal identity is claimed.",
            "K0 transport of U,K,P,Ux uses the certified physical relation d=2*log(X/a).",
        ],
    }


build = propagate
enclose = propagate
