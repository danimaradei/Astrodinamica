import math

def KepEqtnE(M, e, abstol=1e-6):
    """
    Resuelve la Ecuación de Kepler para órbitas elípticas (Algoritmo 2 de Vallado).
    
    Parámetros:
        M (float): Anomalía media en radianes.
        e (float): Excentricidad (0 <= e < 1).
        abstol (float): Tolerancia de convergencia.
        
    Retorna:
        tuple: (E, iterations) -> Anomalía excéntrica en rad y cantidad de iteraciones.
    """
    # 1. Normalizar M entre -pi y pi para seleccionar la estimación inicial de Vallado
    M_norm = math.atan2(math.sin(M), math.cos(M))
    
    # 2. Estimación inicial (E0) según Algoritmo 2 de Vallado
    if (-math.pi < M_norm < 0) or (M_norm > math.pi):
        E = M - e
    else:
        E = M + e
        
    iterations = 0
    
    # 3. Bucle de iteración por Newton-Raphson
    while True:
        iterations += 1
        
        # Cálculo del paso delta = (M - E + e*sin(E)) / (1 - e*cos(E))
        delta = (M - E + e * math.sin(E)) / (1.0 - e * math.cos(E))
        E = E + delta
        
        # Criterio de parada
        if abs(delta) < abstol:
            break
            
    return E, iterations


def GetTOF(E1, E2, e, a, mu=398600.4418):
    """
    Calcula el tiempo de vuelo (TOF) entre dos anomalías excéntricas.
    
    Parámetros:
        E1 (float): Anomalía excéntrica inicial en radianes.
        E2 (float): Anomalía excéntrica final en radianes.
        e (float): Excentricidad.
        a (float): Semieje mayor en km.
        mu (float): Parámetro gravitacional estándar (km^3/s^2). Por defecto: Tierra.
        
    Retorna:
        float: Tiempo de vuelo en segundos.
    """
    # Movimiento medio
    n = math.sqrt(mu / (a ** 3))
    
    # Convertir E1 y E2 a anomalías medias M1 y M2
    M1 = E1 - e * math.sin(E1)
    M2 = E2 - e * math.sin(E2)
    
    delta_M = M2 - M1
    
    # Tiempo de vuelo: TOF = delta_M / n
    tof = delta_M / n
    return tof


# ==============================================================================
# EJECUCIÓN Y TESTS DE CONVERGENCIA (PUNTO 3)
# ==============================================================================
if __name__ == "__main__":
    print("=== 1. TEST DE FUNCIONAMIENTO BÁSICO ===")
    a_test = 7000.0        # Semieje mayor de 7000 km
    e_test = 0.1           # Excentricidad de 0.1
    M_test = math.pi / 4   # M = 45 grados (en rad)
    
    # Prueba de la Ecuación de Kepler
    E_calculado, iters = KepEqtnE(M_test, e_test, abstol=1e-6)
    print(f"Para M = {M_test:.4f} rad y e = {e_test}:")
    print(f" -> E calculado = {E_calculado:.6f} rad ({iters} iteraciones)")
    
    # Prueba del Tiempo de Vuelo entre E1 = 0 rad y E2 = E_calculado
    tof = GetTOF(0.0, E_calculado, e_test, a_test)
    print(f" -> Tiempo de vuelo (E1=0 a E2={E_calculado:.4f}): {tof:.2f} segundos\n")
    
    print("=== 2. ANÁLISIS DE CONVERGENCIA SEGÚN LA EXCENTRICIDAD ===")
    print("Tolerancia fija: abstol = 1e-6 | M = 0.5 rad\n")
    
    M_eval = 0.5
    excentricidades = [0.1, 0.5, 0.8, 0.9, 0.99, 0.999]
    
    print(f"{'Excentricidad (e)':<20} | {'Iteraciones':<12} | {'E Resultado (rad)':<20}")
    print("-" * 58)
    for e in excentricidades:
        E_res, num_it = KepEqtnE(M_eval, e, abstol=1e-6)
        print(f"{e:<20} | {num_it:<12} | {E_res:<20.6f}")