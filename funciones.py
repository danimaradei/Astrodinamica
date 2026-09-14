import numpy as np

# ---------------------------------------------------------------
# Constantes del problema (Guia 1 y enunciado del punto 7)
# ---------------------------------------------------------------
Rt = 6371.0          # radio terrestre [km] (asumido en la Guia 1)
mu = 398600.0        # parametro gravitacional terrestre [km^3/s^2]

h_p = 500.0          # altura de perigeo [km]
h_a = 1000.0         # altura de apogeo [km]


# ---------------------------------------------------------------
# Punto 1: ecuacion de Kepler (Algoritmo 2, Vallado)
# ---------------------------------------------------------------
class NonEllipticOrbitException(Exception):
    pass


def KepEqtnE(M, e, abstol=1e-6, maxiter=200):
    """Devuelve la anomalia excentrica E dada la anomalia media M y e.
    Newton-Raphson, valor inicial segun Vallado Alg. 2."""
    if -np.pi < M < 0 or M > np.pi:
        E = M - e
    else:
        E = M + e

    for _ in range(maxiter):
        dE = (M - E + e * np.sin(E)) / (1 - e * np.cos(E))
        E += dE
        if abs(dE) < abstol:
            break
    return E


# ---------------------------------------------------------------
# Punto 2: conversion entre anomalia verdadera y excentrica
# (Algoritmos 5 y 6, Vallado). SOLO caso eliptico/circular:
# para este ejercicio no hace falta parabolico/hiperbolico (punto 3).
# ---------------------------------------------------------------
def TrAnomaly2EccAnomaly(e, nu):
    if 0 <= e < 1:
        sinE = np.sin(nu) * np.sqrt(1 - e**2) / (1 + e * np.cos(nu))
        cosE = (e + np.cos(nu)) / (1 + e * np.cos(nu))
        return np.arctan2(sinE, cosE)
    else:
        raise NonEllipticOrbitException


def EccAnomaly2TrAnomaly(e, E):
    if 0 <= e < 1:
        sinNu = np.sin(E) * np.sqrt(1 - e**2) / (1 - e * np.cos(E))
        cosNu = (np.cos(E) - e) / (1 - e * np.cos(E))
        return np.arctan2(sinNu, cosNu)
    else:
        raise NonEllipticOrbitException


# ---------------------------------------------------------------
# Punto 4: RV2COE (Algoritmo 9, Vallado).
# Version rustica: SOLO caso general (e != 0, orbita no ecuatorial).
# No hace falta mas para este ejercicio (nuestra orbita es polar,
# i=90 => vector nodal bien definido, y e != 0).
# ---------------------------------------------------------------
def RV2COE(r_vec, v_vec, mu):
    r_vec = np.array(r_vec, dtype=float)
    v_vec = np.array(v_vec, dtype=float)

    r = np.linalg.norm(r_vec)
    v = np.linalg.norm(v_vec)

    h_vec = np.cross(r_vec, v_vec)
    h = np.linalg.norm(h_vec)

    n_vec = np.cross([0, 0, 1], h_vec)
    n = np.linalg.norm(n_vec)

    e_vec = ((v**2 - mu / r) * r_vec - np.dot(r_vec, v_vec) * v_vec) / mu
    e = np.linalg.norm(e_vec)

    energy = v**2 / 2 - mu / r
    a = -mu / (2 * energy)

    i = np.arccos(np.clip(h_vec[2] / h, -1.0, 1.0))

    Omega = np.arccos(np.clip(n_vec[0] / n, -1.0, 1.0))
    if n_vec[1] < 0:
        Omega = 2 * np.pi - Omega

    omega = np.arccos(np.clip(np.dot(n_vec, e_vec) / (n * e), -1.0, 1.0))
    if e_vec[2] < -1e-9:   # en vez de < 0
        omega = 2 * np.pi - omega

    nu = np.arccos(np.clip(np.dot(e_vec, r_vec) / (e * r), -1.0, 1.0))
    if np.dot(r_vec, v_vec) < 0:
        nu = 2 * np.pi - nu

    return a, e, i, Omega, omega, nu


# ---------------------------------------------------------------
# Punto 5: COE2RV (Algoritmo 10, Vallado)
# ---------------------------------------------------------------
def COE2RV(coe, mu):
    a, e, i, Omega, omega, nu = coe

    p = a * (1 - e**2) # semi-latus rectum

    # Posicion y velocidad en el marco PERIFOCAL (PQW): P apunta al perigeo, Q esta 90 grados adelante en el plano orbital, W = P x Q (normal)
    # !!!!!!!! BUSCAR EN LAS DIAPOS FOTO DE LAS FORMULAS DE R,V EN EL MARCO PERIFOCAL
    r_pqw = np.array([p * np.cos(nu) / (1 + e * np.cos(nu)), p * np.sin(nu) / (1 + e * np.cos(nu)), 0.0]) # r_pqw = (r_p,0,0) porque nu=0
    v_pqw = np.sqrt(mu / p) * np.array([-np.sin(nu), e + np.cos(nu), 0.0])

    # Matriz de rotacion PQW -> IJK (inercial): secuencia 3-1-3 
    cO, sO = np.cos(Omega), np.sin(Omega)
    ci, si = np.cos(i), np.sin(i)
    cw, sw = np.cos(omega), np.sin(omega)

    R = np.array([
        [cO*cw - sO*sw*ci, -cO*sw - sO*cw*ci, sO*si],
        [sO*cw + cO*sw*ci, -sO*sw + cO*cw*ci, -cO*si],
        [sw*si, cw*si, ci  ]
    ])

    r_ijk = R @ r_pqw
    v_ijk = R @ v_pqw
    return r_ijk, v_ijk


# ---------------------------------------------------------------
# Punto 6: simulador de orbitas elipticas/circulares (2 cuerpos)
# ---------------------------------------------------------------
def simular(r_vec0, v_vec0, t, mu):
    """t: vector de tiempos [s]. Devuelve X (6xN) con [r;v] en cada t."""
    N = len(t)
    X = np.zeros((6, N))
    X[:, 0] = np.concatenate([r_vec0, v_vec0])

    a, e, i, Omega, omega, nu0 = RV2COE(r_vec0, v_vec0, mu)
    n_medio = np.sqrt(mu / a**3)

    nu = np.zeros(N)
    nu[0] = nu0

    for k in range(1, N):
        dt = t[k] - t[k - 1]

        E_prev = TrAnomaly2EccAnomaly(e, nu[k - 1])
        M_prev = E_prev - e * np.sin(E_prev)
        M_new = (M_prev + n_medio * dt) % (2 * np.pi)
        E_new = KepEqtnE(M_new, e)
        nu[k] = EccAnomaly2TrAnomaly(e, E_new)

        r_k, v_k = COE2RV((a, e, i, Omega, omega, nu[k]), mu)
        X[:, k] = np.concatenate([r_k, v_k])

    return X


# ---------------------------------------------------------------
# Punto 7: caso particular (orbita polar, perigeo en el ecuador)
# ---------------------------------------------------------------
if __name__ == "__main__":
    # --- elementos orbitales calculados "a mano" (ver discusion previa) ---
    rp = Rt + h_p
    ra = Rt + h_a
    a0 = (rp + ra) / 2
    e0 = (ra - rp) / (ra + rp)
    i0 = np.deg2rad(90)        # orbita polar
    Omega0 = np.deg2rad(90)    # GST(50) + long_nodo(40)
    omega0 = np.deg2rad(0)     # perigeo coincide con el nodo ascendente
    nu0 = np.deg2rad(0)        # t=0 en el perigeo

    print(f"a = {a0:.3f} km   e = {e0:.5f}")

    # --- posicion y velocidad iniciales (COE2RV, punto 5) ---
    r0, v0 = COE2RV((a0, e0, i0, Omega0, omega0, nu0), mu)
    print(f"r0 = {r0} km")
    print(f"v0 = {v0} km/s")

    # --- propagacion con el simulador del punto 6 ---
    T = 2 * np.pi * np.sqrt(a0**3 / mu)     # periodo orbital [s]
    t = np.linspace(0, T, 200)
    X = simular(r0, v0, t, mu)

    # --- verificacion: recalcular COE en cada instante con RV2COE ---
    print("\n t[s]      a[km]      e        i[deg]   Omega[deg]  omega[deg]  nu[deg]")
    for k in [0, 50, 100, 150, 199]:
        a_k, e_k, i_k, O_k, w_k, nu_k = RV2COE(X[0:3, k], X[3:6, k], mu)
        print(f"{t[k]:7.1f}  {a_k:9.3f}  {e_k:7.5f}  {np.rad2deg(i_k):7.2f}  "
              f"{np.rad2deg(O_k):9.2f}  {np.rad2deg(w_k):9.2f}  {np.rad2deg(nu_k):7.2f}")