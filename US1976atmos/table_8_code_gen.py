import numpy as np
from scipy.interpolate import interp1d, PchipInterpolator
import json
from typing import List
import pickle
from importlib.resources import files


input_file = "table_8.json"
output_file = "table_8.h"

r0 = 6.356766e3  # km
g0 = 9.80665  # m/s^2
Rstr = 8.31432e3  # N*m/(kmol*K)
M0 = 28.9644  # kg/kmol page 9
Na = 6.022169e26
T7 = 186.8673  # K
Lmb = [-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0]  # Table 4 K/km
# Calculated form table 4
Tbn = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 215.65, 186.946]

_data_file = "_table_8.pickle"
try:
    with open(_data_file, "rb") as handle:
        _table_eight = pickle.load(handle)
except FileNotFoundError:
    with open(input_file) as infile:
        _table_eight = json.load(infile)

    molar_mass = [28.0134, 15.999, 15.999*2, 39.948, 4.002602, 1.008]
    # Some sums we fill in
    z = np.asarray(_table_eight['altitude Z'])
    _table_eight["sum_ni"] = np.zeros_like(z)
    _table_eight["sum_mass_density_g/m^3"] = np.zeros_like(z)
    for i, key in enumerate(['N2', 'O', 'O2', 'Ar', 'He', 'H']):
        nd = 1000.0 * np.asarray(_table_eight[key]) / Na  # moles / m^3
        _table_eight["sum_ni"] += nd
        rho = molar_mass[i] * nd
        _table_eight[f'{key}_mass_density_g/m^3'] = rho
        _table_eight["sum_mass_density_g/m^3"] += rho
        y = np.asarray(np.log10(_table_eight[key]))
        idx = np.isfinite(y)
        try:
            _table_eight[f'{key}_zero_under'] = z[idx][0]
        except IndexError:
            _table_eight[f'{key}_zero_under'] = -5000.0
        _table_eight[f'{key}_log10_number_density'] = PchipInterpolator(
            z[idx], y[idx], extrapolate=True
        )
        _table_eight[f'{key}_log10_mass_density'] = PchipInterpolator(
            z[idx], y[idx]*molar_mass[i], extrapolate=True
        )

    _table_eight['log10_sum_ni_z'] = PchipInterpolator(
        z, np.log10(_table_eight["sum_ni"]), extrapolate=True
    )
    _table_eight['log10_sum_rho_z'] = PchipInterpolator(
        z, np.log10(_table_eight["sum_mass_density_g/m^3"]), 
        extrapolate=True
    )

    _data_file = "_table_8.pickle"

    with open(_data_file, "wb") as handle:
        pickle.dump(_table_eight, handle, protocol=pickle.HIGHEST_PROTOCOL)

    del z


def variable_string(v_name: str, data: List[float]) -> str:
    x_str = [str(x) for x in data]
    out = f"double {v_name}[] = " + "{" + ",".join(x_str) + "};\n"
    return out


def H_to_Z(H: np.ndarray) -> np.ndarray:
    "Take in geopotential height and spits out geometric height. "
    return (r0*H)/(r0-H)


def Z_to_H(Z: np.ndarray) -> np.ndarray:
    "Take in geometric height and spits out geopotentiol height. "
    return (r0*Z)/(r0+Z)

H8 = [79, 79.5, 80, 80.5, 81, 81.5, 82, 82.5, 83, 83.5, 84, 84.5]

mr8 = [1.0,
       0.999996,
       0.999988,
       0.999969,
       0.999938,
       0.999904,
       0.999864,
       0.999822,
       0.999778,
       0.999731,
       0.999681,
       0.999679]
MtoM0_H = PchipInterpolator(H8, mr8)

# Generate data for interpolation in the lower atmosphere
H_b = np.asarray([0, 11, 20, 32, 47, 51, 71, 84.8520])
T_b = [288.15]
for i, lapse_rate in enumerate([-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0]):
    T_b.append(T_b[i] + lapse_rate*(H_b[i+1] - H_b[i]))

t_H_0_79 = interp1d(H_b, T_b, kind='linear', bounds_error=False,
                    fill_value='extrapolate')

T1000 = 1000


def t_K_79_h86(H: np.ndarray) -> np.ndarray:
    return t_H_0_79(H) * MtoM0_H(H)


def t_Z_86_91(Z: np.ndarray) -> np.ndarray:
    return np.ones_like(Z) * 186.8643


def t_Z_91_110(Z: np.ndarray) -> np.ndarray:
    Tc = 263.1905
    A = -76.3232
    a = -19.9429
    return Tc + A*np.sqrt(1 - ((Z - 91.0)/a)**2.0)


def t_Z_110_120(Z: np.ndarray) -> np.ndarray:
    T9 = 240.0
    Z9 = 110.0
    L9 = 12.0
    return T9 + L9 * (Z - Z9)


def xi_Z(Z: np.ndarray) -> np.ndarray:
    return (Z - 120) * (r0 + 120) / (r0 + Z)


def t_Z_120_1000(Z: np.ndarray) -> np.ndarray:
    T10 = 360.0
    lmbda = 0.01875
    return T1000 - (T1000 - T10) * np.exp(-lmbda * xi_Z(Z))


def _to_t_H(Z: np.ndarray) -> np.ndarray:
    " I'm the only function that should make this call. "
    return _t_H_2(Z_to_H(Z))


def t_Z(Z: np.ndarray) -> np.ndarray:
    return np.piecewise(Z, [Z < 86.0,
                            np.logical_and(86.0 <= Z, Z < 91.0,),
                            np.logical_and(91.0 <= Z, Z < 110.0,),
                            np.logical_and(110.0 <= Z, Z < 120.0,),
                            np.logical_and(120.0 <= Z, Z < 1000.0,),
                            Z >= 1000.0],
                            [_to_t_H,
                             t_Z_86_91,
                             t_Z_91_110,
                             t_Z_110_120,
                             t_Z_120_1000,
                             T1000])


def _to_t_Z(H: np.ndarray) -> np.ndarray:
    return t_Z(H_to_Z(H))


def t_H(H: np.ndarray) -> np.ndarray:
    h86 = Z_to_H(86)
    hm5 = Z_to_H(-5)
    return np.piecewise(H, [H < hm5,
                            np.logical_and(hm5 <= H, H <=79),
                            np.logical_and(79 < H, H <= h86),
                            H > h86],
                            [0.0,
                             t_H_0_79,
                             t_K_79_h86,
                             _to_t_Z])


def _t_H_2(H: np.ndarray) -> np.ndarray:
    "You probably shouldn't be calling this."
    hm5 = Z_to_H(-5)
    return np.piecewise(H, [H < hm5, hm5 <= H],
                        [0.0, t_H_0_79])


def equation33(Pb, Lb, Hb, H, Tb):
    # See Equations 33a and 33b, page 12 of the US 1976 Standard Atmosphere
    if Lb == 0.0:
        return equation33b(Pb, Hb, Tb, H)
    else:  # Equation 33a
        exponent = g0 * M0 / (Rstr * Lb)*1000
        return Pb * (Tb / (Tb + Lb * (H-Hb)))**exponent


def equation33b(Pb, Hb, Tb, H):
    # See Equation 33d, page 12 of the US 1976 Standard Atmosphere report
    numerator = -g0 * M0 * (H - Hb) * 1000
    denomenator = Rstr * Tb
    return Pb * np.exp(numerator / denomenator)


def sum_ni(Z):
    return 10.0**(_table_eight['log10_sum_ni_z'](Z))


def _p_H0(H: np.ndarray) -> np.ndarray:
    Hb = 0
    Lb = Lmb[0]
    Pb = 1013.25
    Tb = Tbn[0]
    return equation33(Pb, Lb, Hb, H, Tb)


def _p_H1(H: np.ndarray) -> np.ndarray:
    Hb = 11
    Lb = 0.0
    Pb = (_p_H0(Hb))
    Tb = Tbn[1]
    return equation33(Pb, Lb, Hb, H, Tb)


def _p_H2(H: np.ndarray) -> np.ndarray:
    Hb = 20
    Lb = 1.0
    Pb = (_p_H1(Hb))
    Tb = Tbn[2]
    return equation33(Pb, Lb, Hb, H, Tb)


def _p_H3(H: np.ndarray) -> np.ndarray:
    Hb = 32
    Lb = 2.8
    Pb = (_p_H2(Hb))
    Tb = Tbn[3]
    return equation33(Pb, Lb, Hb, H, Tb)


def _p_H4(H: np.ndarray) -> np.ndarray:
    Hb = 47
    Lb = 0.0
    Pb = (_p_H3(Hb))
    Tb = Tbn[4]
    return equation33(Pb, Lb, Hb, H, Tb)


def _p_H5(H: np.ndarray) -> np.ndarray:
    Hb = 51
    Lb = -2.8
    Pb = (_p_H4(Hb))
    Tb = Tbn[5]
    return equation33(Pb, Lb, Hb, H, Tb)


def _p_H6(H: np.ndarray) -> np.ndarray:
    Hb = 71
    Lb = -2.0
    Pb = (_p_H5(Hb))
    Tb = Tbn[6]
    return equation33(Pb, Lb, Hb, H, Tb)


def _p_Z_86(Z: np.ndarray) -> np.ndarray:
    # See page 12.
    T = t_Z(Z)
    return sum_ni(Z*1000) * Rstr * T / 100 / 1000


def _to_P_Z(H: np.ndarray) -> np.ndarray:
    " You shouldn't be callnig this function. "
    return _p_Z_86(H_to_Z(H))


def p_H(H: np.ndarray) -> np.ndarray:
    h86 = Z_to_H(86)
    hm5 = Z_to_H(-5)
    return np.piecewise(H, [H < hm5,                            # 1
                            np.logical_and(hm5 <= H, H < 11),   # 2
                            np.logical_and(11 <= H, H < 20),    # 3
                            np.logical_and(20 <= H, H < 32),    # 4
                            np.logical_and(32 <= H, H < 47),    # 5
                            np.logical_and(47 <= H, H < 51),    # 6
                            np.logical_and(51 <= H, H < 71),    # 7
                            np.logical_and(71 <= H, H < h86),   # 8
                            H >= h86,                           # 9    
                            ],
                            [0,       # 1
                             _p_H0,   # 2
                             _p_H1,   # 3
                             _p_H2,   # 4
                             _p_H3,   # 5
                             _p_H4,   # 6
                             _p_H5,   # 7
                             _p_H6,   # 8
                             _to_P_Z  # 9
                            ])
    # I'm the only function thta should call _to_P_Z


def _to_P_H(Z: np.ndarray) -> np.ndarray:
    return p_H(Z_to_H(Z))


def p_Z(Z: np.ndarray) -> np.ndarray:
    return np.piecewise(Z, [Z < 86, Z >= 86],
                        [_to_P_H, _p_Z_86])


if __name__ == "__main__":
    with open(input_file) as infile:
        data = json.load(infile)

    Z = np.asarray(data['altitude Z'])

    pressures = p_Z(Z/1000.0)

    with open(output_file, "w") as outfile:
        outfile.write(variable_string('pressures', pressures))

        for i, key in enumerate(['N2','O','O2','Ar','He','H']):
            x = [x * 1000 / Na for x in data[key]]
            outfile.write(variable_string(key, x))

    print("Task Complete")