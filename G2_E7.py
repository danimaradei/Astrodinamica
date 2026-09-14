import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(suppress=True, precision=4)

# Datos
Rt = 6371 # Radio de la Tierra en km
hp = 500 # altura perigeo en km
ha = 1000 # altura apogeo en km
mu = 398600.4418  # parametro gravitacional terrestre en km^3/s^2

# Calculos iniciales 
#!!!!!! Clase 3 Diapo 15
rp = Rt + hp # radio perigeo
ra = Rt + ha # radio apogeo
a = (rp + ra) / 2 # semieje mayor
e = (ra - rp) / (ra + rp) # excentricidad

print("Radio perigeo (km):", rp)
print("Radio apogeo (km):", ra)
print("Semieje mayor (km):", a)
print("Excentricidad:", e)

# Razonamiento
i = np.deg2rad(90) # Inclinación: angulo entre el plano ecuatorial y el orbital. Como la órbita es polar, i=90°
Omega = np.deg2rad(90) # RAAN = theta_GST + lambda_nodo = 50 + 40
omega = np.deg2rad(0) # Argumento del perigeo: angulo en el plano orbital, desde linea de nodos al perigeo
nu = np.deg2rad(0) # Anomalia verdadera: angulo en el plano orbital, desde el perigeo a la posición del satelite. En t=0, es igual a 0 por estar en el perigeo

# COE 2 RV
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


from funciones import simular, RV2COE, mu

r0, v0 = COE2RV((a, e, i, Omega, omega, nu), mu)

print("r0:", r0)
print("v0:", v0)

# Comparacion
T = 2 * np.pi * np.sqrt(a**3 / mu) # periodo orbital en segundos
t = np.linspace(0, T, 200) # 200 instantes de tiempo entre 0 y T

X = simular(r0, v0, t, mu)

for k in [0, 50, 100, 150, 199]:
    a_k, e_k, i_k, O_k, w_k, nu_k = RV2COE(X[0:3, k], X[3:6, k], mu)
    print(f"t={t[k]:7.1f}s  a={a_k:.3f}  e={e_k:.5f}  i={np.rad2deg(i_k):.2f}  "
          f"Omega={np.rad2deg(O_k):.2f}  omega={np.rad2deg(w_k):.2f}  nu={np.rad2deg(nu_k):.2f}")

# "¿Se obtienen los mismos elementos orbitales?"
# Acá la comparación es entre dos fuentes distintas de los elementos orbitales:
# los que calculamos a mano, geométricamente, al principio
# y los que salen de aplicarle RV2COE que devolvió el simulador, en distintos instantes de la propagación
# La respuesta es sí: en los 5 instantes que imprimiste coinciden exactamente con los valores iniciales
# Lo único que cambia es nu, que barre 0 a 360 a lo largo de la vuelta — como tiene que ser
# porque nu es la posición del satélite dentro de la órbita en cada instante

# Grafico: lo hice con MATLAB
