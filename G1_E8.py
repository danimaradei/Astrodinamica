import matplotlib.pyplot as plt
import numpy as np

# Datos
r_s = 6.957*10**8 # Radio del sol en metros
r_p = 8*r_s # Radio perihelio en metros
r_a = 108.2*10**9 # Radio afelio en metros
pg_sol = 1.327*10**20 # Parámetro gravitatorio del sol en m^3/s^2

# Punto A
a = (r_p + r_a)/2 # Semieje mayor en metros
c = a - r_p # Distancia focal en metros
e = c/a # Excentricidad
p = a*(1 - e**2) # Semi-parámetro de la órbita final en metros
b = np.sqrt(a*p) # Semieje menor en metros

print("Punto A")
print(f"Semieje mayor (a): {a:.3e} m")
print(f"Semieje menor (b): {b:.3e} m")
print(f"Excentricidad (e): {e:.3f}")
print(f"Semi-parametro (p): {p:.3e} m")

# Punto B
em = -pg_sol/(2*a) # Energía mecánica específica en m^2/s^2
h = np.sqrt(pg_sol*p) # Momento angular específico en m^2/s

print("\nPunto B")
print(f"Energia mecanica especifica (em): {em:.3f} m^2/s^2")
print(f"Momento angular especifico (h): {h:.3e} m^2/s")

# Punto C
# Sonda
vs_p = (np.sqrt(pg_sol*(2/r_p - 1/a))) * 3.6 # Velocidad en perihelio en km/h
vs_a = (np.sqrt(pg_sol*(2/r_a - 1/a))) * 3.6 # Velocidad en afelio en km/h

print("\nPunto C")
print(f"Velocidad en perihelio (vs_p): {vs_p:.3f} km/h")
print(f"Velocidad en afelio (vs_a): {vs_a:.3f} km/h")

# Punto C - Comparacion con planetas
# Mercurio
mercurio_p = 212328
mercurio_a = 139896
# Venus
venus_p = 126936
venus_a = 125208
# Tierra
tierra_p = 109044
tierra_a = 105444

print("\nComparacion con planetas (km/h)")
print(f"Sonda: perihelio {vs_p:,.0f} | afelio {vs_a:,.0f}")
print(f"Mercurio: perihelio {mercurio_p:,.0f} | afelio {mercurio_a:,.0f}")
print(f"Venus: perihelio {venus_p:,.0f} | afelio {venus_a:,.0f}")
print(f"Tierra: perihelio {tierra_p:,.0f} | afelio {tierra_a:,.0f}")

# Grafico
nombres = ["Sonda", "Mercurio", "Venus", "Tierra"]
v_p = [vs_p, mercurio_p, venus_p, tierra_p]
v_a = [vs_a, mercurio_a, venus_a, tierra_a]
colores = ["crimson", "gray", "goldenrod", "royalblue"]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11,5), sharey=True)

ax1.bar(nombres, v_p, color=colores)
ax1.set_title("Velocidad en Perihelio")
ax1.set_ylabel("Velocidad (km/h)")

ax2.bar(nombres, v_a, color=colores)
ax2.set_title("Velocidad en Afelio")

for ax, valores in [(ax1, v_p), (ax2, v_a)]:
    for i, v in enumerate(valores):
        ax.text(i, v, f"{v:,.0f}", ha="center", va="bottom", fontsize=8)

plt.tight_layout()
plt.show()