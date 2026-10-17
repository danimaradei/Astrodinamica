# Actividad 3B: Guía 3 - Parte 2 - Ejercicio 3
# Maradei Daniela - 62020

from datetime import datetime, timedelta
import numpy as np
import erfa
import matplotlib.pyplot as plt

# DATOS
# Sirio (J2000.0), direccion fija en GCRF
ALPHA = np.radians(15.0 * (6 + 45 / 60 + 8.9173 / 3600)) # 06h 45m 8.9173s
DELTA = -np.radians(16 + 42 / 60 + 58.02 / 3600) # -16° 42' 58.02"

PHI = np.radians(-34.60) # phi gd: latitud geodesica - Buenos Aires
LAM = np.radians(-58.38) # lambda: longitud - Buenos Aires

D_AT = 37.0 # TAI - UTC (s)
TT_MENOS_TAI = 32.184 # (s)
ARCSEC = np.pi / (180 * 3600)

# EOP a las 00:00 UTC del 18, 19, 20 y 21 de octubre de 2026
# Segundos desde el 2026-10-18 00:00:00 UTC
EOP_T = np.array([0.0, 86400.0, 172800.0, 259200.0])
EOP_DUT1 = np.array([-0.03716, -0.03804, -0.03911, -0.04106])
EOP_XP = np.array([0.1543, 0.1531, 0.1519, 0.1508])
EOP_YP = np.array([0.3205, 0.3204, 0.3204, 0.3204])

JD1 = 367 * 2026 - int(7 * (2026 + int((10 + 9) / 12)) / 4) + int(275 * 10 / 9) + 18 + 1721013.5  # Alg. 14
# Ventana: 18/10 00:00 local (UTC-3) a 20/10 00:00 local
T_INI = 3 * 3600.0 # 03:00 UTC del 18/10
T_FIN = 2 * 86400.0 + 3 * 3600.0 # 03:00 UTC del 20/10
PASO = 60.0 # (s)


# MATRICES DE ROTACION
def rot1(a):
    a = np.atleast_1d(a)
    c, s = np.cos(a), np.sin(a)
    m = np.zeros(a.shape + (3, 3))
    m[..., 0, 0] = 1
    m[..., 1, 1], m[..., 1, 2] = c, s
    m[..., 2, 1], m[..., 2, 2] = -s, c
    return m

def rot2(a):
    a = np.atleast_1d(a)
    c, s = np.cos(a), np.sin(a)
    m = np.zeros(a.shape + (3, 3))
    m[..., 0, 0], m[..., 0, 2] = c, -s
    m[..., 1, 1] = 1
    m[..., 2, 0], m[..., 2, 2] = s, c
    return m

def rot3(a):
    a = np.atleast_1d(a)
    c, s = np.cos(a), np.sin(a)
    m = np.zeros(a.shape + (3, 3))
    m[..., 0, 0], m[..., 0, 1] = c, s
    m[..., 1, 0], m[..., 1, 1] = -s, c
    m[..., 2, 2] = 1
    return m


# ERA
def era(jd1, jd2):
    f = np.mod(jd1, 1.0) + np.mod(jd2, 1.0)
    t = (jd1 - 2451545.0) + jd2
    return 2 * np.pi * np.mod(f + 0.7790572732640 + 0.00273781191135448 * t, 1.0)


# TRANSFORMACIONES
# Ec 19
U_GCRF = np.array([np.cos(DELTA) * np.cos(ALPHA),
                   np.cos(DELTA) * np.sin(ALPHA),
                   np.sin(DELTA)])    
M_LOCAL = (rot2(np.pi / 2 - PHI) @ rot3(LAM))[0] # Ec 20

def elevacion(t):
    """Elevacion [rad] de Sirio para t = segundos UTC desde 2026-10-18 00:00."""
    t = np.atleast_1d(np.asarray(t, dtype=float))
    dut1 = np.interp(t, EOP_T, EOP_DUT1)
    xp = np.interp(t, EOP_T, EOP_XP) * ARCSEC
    yp = np.interp(t, EOP_T, EOP_YP) * ARCSEC

    jd2_ut1 = t / 86400.0 + dut1 / 86400.0 # JD_UT1 = JD1 + jd2
    tt1 = 2451545.0
    tt2 = (JD1 - 2451545.0) + t / 86400.0 + (D_AT + TT_MENOS_TAI) / 86400.0 # JD_TT = tt1 + tt2

    pn_t = erfa.c2i06a(tt1, tt2) # = [PN]^T (GCRS -> CIRS)
    r_t = rot3(era(JD1, jd2_ut1)) # [R]^T = ROT3(ERA)
    w_t = rot1(-yp) @ rot2(-xp) # [W]^T

    u = np.einsum('nij,j->ni', pn_t, U_GCRF) # CIRS
    u = np.einsum('nij,nj->ni', r_t, u) # TIRS
    u = np.einsum('nij,nj->ni', w_t, u) # ITRF (Ec. 15)
    s, e, z = (M_LOCAL @ u.T) # SEZ (Ec. 20)
    return np.arctan2(z, np.hypot(s, e)) # El (Ec. 21)


# CRUCES DEL HORIZONTE Y CULMINACIONES
def biseccion(f, a, b, it=60):
    fa = f(a)
    for _ in range(it):
        m = 0.5 * (a + b)
        fm = f(m)
        if fa * fm <= 0:
            b = m
        else:
            a, fa = m, fm
    return 0.5 * (a + b)

def derivada(t, h=1.0):
    return float((elevacion(t + h) - elevacion(t - h))[0] / (2 * h))

t = np.arange(T_INI, T_FIN + 0.5 * PASO, PASO)
el = elevacion(t)
sobre = el > 0

# Cruces del horizonte
cruces = []
for i in range(len(t) - 1):
    if sobre[i] != sobre[i + 1]:
        tc = biseccion(lambda x: float(elevacion(x)[0]), t[i], t[i + 1])
        cruces.append((tc, 'sube' if sobre[i + 1] else 'baja'))

# Intervalos con El > 0
intervalos = []
inicio = t[0] if sobre[0] else None
for tc, tipo in cruces:
    if tipo == 'sube':
        inicio = tc
    else:
        intervalos.append((inicio, tc))
        inicio = None
if inicio is not None:
    intervalos.append((inicio, t[-1]))

# Culminaciones (maximos locales de El)
culminaciones = []
for i in range(1, len(t) - 1):
    if el[i] > 0 and el[i] > el[i - 1] and el[i] >= el[i + 1]:
        tm = biseccion(derivada, t[i - 1], t[i + 1])
        culminaciones.append((tm, float(np.degrees(elevacion(tm)[0]))))


# RESULTADOS (hora local UTC-3)
def local(ts):
    return datetime(2026, 10, 18) + timedelta(seconds=float(ts) - 3 * 3600)

def fmt(ts):
    return local(ts).strftime('%d/%m/%Y %H:%M:%S')

print('Intervalos con El > 0 (hora local, UTC-3):')
for a, b in intervalos:
    nota_a = ' (inicio de la ventana)' if a == t[0] else ''
    nota_b = ' (fin de la ventana)' if b == t[-1] else ''
    print(f'  {fmt(a)}{nota_a}  ->  {fmt(b)}{nota_b}')
print('Culminaciones:')
for tm, em in culminaciones:
    print(f'  {fmt(tm)}   El = {em:.3f} deg')


# GRAFICO
x = (t - 3 * 3600.0) / 3600.0 # horas locales desde el 18/10 00:00
el_deg = np.degrees(el)

fig, ax = plt.subplots(figsize=(11, 4.5))
ax.plot(x, el_deg, color='tab:blue', lw=1.6, label='Elevación de Sirio')
ax.axhline(0, color='k', lw=0.8)
ax.fill_between(x, 0, el_deg, where=el_deg > 0, color='tab:blue', alpha=0.15)
if culminaciones:
    xc = [(tm - 3 * 3600.0) / 3600.0 for tm, _ in culminaciones]
    ec = [em for _, em in culminaciones]
    ax.plot(xc, ec, 'o', color='tab:red', label='Culminaciones')
    for xi, ei, (tm, _) in zip(xc, ec, culminaciones):
        ax.annotate(f'{local(tm).strftime("%d/%m %H:%M")}\n{ei:.2f}°', (xi, ei),
                    textcoords='offset points', xytext=(0, 8), ha='center', fontsize=8)
ax.set_xlim(0, 48)
ax.set_ylim(-80, 90)
ticks = np.arange(0, 49, 6)
ax.set_xticks(ticks)
ax.set_xticklabels([(datetime(2026, 10, 18) + timedelta(hours=int(h))).strftime('%d/%m\n%H:%M') for h in ticks])
ax.set_ylabel('Elevación [°]')
ax.set_xlabel('Hora local (UTC−3)')
ax.set_title('Elevación de Sirio desde Buenos Aires, 18 al 20 de octubre de 2026')
ax.grid(alpha=0.3)
ax.legend(loc='lower right')
fig.tight_layout()
fig.savefig('G3P2_E3_elevacion_sirio.png', dpi=150)
plt.show()
