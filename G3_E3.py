# Daniela Maradei 62020 - Actividad 3: Guía 3A - Ejercicio 3

from geomodels import GeoidModel

# Coordenadas Centro Espacial
lat_gd = -31.524075 # deg
lon    = -64.463522 # deg
H_msl  = 730 # m

# A: Ondulación del geoide N
geoid = GeoidModel("egm96-5", path=r"C:\geographiclib\geoids")
N = geoid(lat_gd, lon) # altura del geoide sobre el elipsoide en m
print(f"A: N = {N:.3f} m")

# B: Altura elipsoidal h_ellp
h_ellp = H_msl + N
print(f"B: h_ellp = {h_ellp:.3f} m")

# C: Posición ECEF de la estación
import numpy as np

# Parametros de referencia WGS-84 (guía)
a  = 6378137 # Semi-eje mayor de esferoide en m
f  = 1/298.2572235630 # Achatamiento de esferoide
e2 = 2*f - f**2 # Excentricidad al cuadrado

phi = np.radians(lat_gd)
lam = np.radians(lon)

den = np.sqrt(1 - e2*np.sin(phi)**2)
C = a/den # Vallado Ec. 3-7
S = a*(1 - e2)/den #Vallado Ec. 3-7

# Coordenadas cartesianas de la estación en ECEF
x = (C + h_ellp) * np.cos(phi) * np.cos(lam)
y = (C + h_ellp) * np.cos(phi) * np.sin(lam)
z = (S + h_ellp) * np.sin(phi)

r_ecef = np.array([x, y, z])
print(f"C: r_ECEF = [{r_ecef[0]:.3f}, {r_ecef[1]:.3f}, {r_ecef[2]:.3f}] m")
print(f"   |r| = {np.linalg.norm(r_ecef):.3f} m")

# D: Longitud y latitud geocéntricas desde el vector ECEF
x, y, z = r_ecef
lon_gc = np.degrees(np.arctan2(y, x))
lat_gc = np.degrees(np.arctan2(z, np.hypot(x, y)))   # hypot = sqrt(x² + y²)
print(f"D: lon_gc = {lon_gc:.6f} deg")
print(f"   lat_gc = {lat_gc:.6f} deg")

# E: Diferencia entre latitud geodésica y geocéntrica
dphi = lat_gd - lat_gc
print(f"E: dphi = {dphi:.6f} deg")

# F: Interpretación del signo de N
print("F: Interpretacion del signo de N")
print("   La ondulacion (N = h_ellp - H_msl) mide la diferencia de altura entre el geoide y el elipsoide de referencia. Al obtener un valor de +27.187m, se comprueba que en esta ubicacion el geoide EGM96 se encuentra fisicamente por encima del elipsoide WGS-84. Como consecuencia, la altura elipsoidal de la estacion resulta mayor que su altura sobre el nivel medio del mar (757.187m en vez de 730m). Si se hubiera utilizado directamente la elevacion ortometrica H_msl como altura elipsoidal para las transformaciones, se habria subestimado la posicion de la estacion introduciendo un error de mas de 27m sobre la normal al elipsoide.")
