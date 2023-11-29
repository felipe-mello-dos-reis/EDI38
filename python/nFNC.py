import numpy as np
import matplotlib.pyplot as plt
from math import pi

from nFNC_functions import *

# Specifications for concrete
def main():
    fck = 30  # MPa
    gamma_c = 1.4
    sigma_cd = 0.85 * fck / gamma_c  # MPa
    gamma_conc = 2500  # kg/m3
    
    # Specifications for steel
    f_yk = 500  # MPa
    gamma_s = 1.15
    E_s = 210  # GPa
    f_yd = f_yk / gamma_s  # MPa
    epsilon_yd = f_yd / E_s  # per mil
    gamma_aco = 7850  # kg/m3
    
    # Design specifications
    c = 0.02  # m
    b = 0.20  # m
    h = 0.32  # m
    d = 0.05  # m
    nc = 3
    nb = [2, 2, 2]
    phi = 0.010  # m
    y_s = np.linspace(-(h / 2 - d), (h / 2 - d), nc)
    
    # Applied forces
    N_d = 0.8160  # MN
    M_d = 0.0373 # MN*m
    
    # Initializing variables
    epsilon_c2 = 0
    epsilon_cu = 0
    x_lim = 0
    n = 0
    
    # Tolerances
    tol_J = 1e-9
    tol_k = 1e-9
    tol_f = 1e-9
    i = 0
    it_max = 1e6
    
    # Problem specifications
    if fck <= 50:
        x_lim = 0.45 * d
        n = 2
        epsilon_c2 = 2
        epsilon_cu = 3.5
    elif 50 < fck < 90:
        x_lim = 0.35 * d
        n = 1.4 + 23.4 * ((90 - fck) / 100) ** 4
        epsilon_c2 = 2 + 0.085 * (fck - 50) ** 0.53
        epsilon_cu = 2.6 + 35 * ((90 - fck) / 100) ** 4
    
    # Plotting resistant forces
    y_t = h / 2
    y_b = -h / 2
    
    
    # Constants
    epsilon_0 = 0
    k = 0
    epsilon_0_it = [epsilon_0]
    k_it = [k]
    
    # Calculate N, M, epsilon_t, and epsilon_b
    N_s = Ns(y_s, nb, phi, f_yd, epsilon_yd, epsilon_0, k)
    M_s = Ms(y_s, nb, phi, f_yd, epsilon_yd, epsilon_0, k)
    
    N_c = Nc(epsilon_0, k, epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k)
    M_c = Mc(epsilon_0, k, epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k)
    
    N_r = N_c + N_s
    M_r = M_c + M_s
    
    epsilon_t = epsilon_0 + y_t * k
    epsilon_b = epsilon_0 + y_b * k
    
    # Newton-Raphson method
    i = 0
    f_ad = np.sqrt(((N_d - N_r) / (sigma_cd * b * h))**2 + ((M_d - M_r) / (sigma_cd * b * h**2))**2)
    
    while abs(f_ad) > tol_f and i < it_max:
        i += 1
        epsilon_t = epsilon_0 + k * y_t
        epsilon_b = epsilon_0 + k * y_b
    
        N_s = Ns(y_s, nb, phi, f_yd, epsilon_yd, epsilon_0, k)
        M_s = Ms(y_s, nb, phi, f_yd, epsilon_yd, epsilon_0, k)
    
        N_c = Nc(epsilon_0, k, epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k)
        M_c = Mc(epsilon_0, k, epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k)
    
        N_r = N_c + N_s
        M_r = M_c + M_s
    
        J_c = Jc(epsilon_t, epsilon_b, epsilon_0, epsilon_c2, sigma_cd, n, b, k, h, tol_k)
        J_s = Js(E_s, epsilon_yd, phi, y_s, nb, epsilon_0, k)
        J = J_c + J_s
    
        J_ad = J[0, 0] * J[1, 1] / (sigma_cd * b * h * sigma_cd * b * h * h**2) - J[0, 1] * J[1, 0] / ((sigma_cd * b * h**2)**2)
    
        if abs(J_ad) > tol_J:
            J_inv = 1 / (J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]) * np.array([[J[1, 1], -J[1, 0]], [-J[0, 1], J[0, 0]]])
            f = np.array([[N_d - N_r], [M_d - M_r]])
            new_it = np.array([[epsilon_0], [k]]) - np.dot(J_inv, f)
            epsilon_0 = new_it[0, 0]
            epsilon_0_it.append(epsilon_0)
            k = new_it[1, 0]
            k_it.append(k)
            f_ad = np.sqrt(((N_d - N_r) / (sigma_cd * b * h))**2 + ((M_d - M_r) / (sigma_cd * b * h**2))**2)
        else:
            print("No solution exists!")
            break
        
    # Check ELU
    ELU(epsilon_0, k, y_b, y_t, y_s, epsilon_c2, epsilon_cu)
    
    # Plot epsilon_0, k
    plot_epsilon_0_k(epsilon_c2, epsilon_cu, h, y_t, y_b, y_s, epsilon_0_it, k_it)
    plot_N_r_M_r(epsilon_c2, epsilon_cu, sigma_cd, n, b, f_yd, epsilon_yd, h, y_t, y_b, y_s, phi, nb, tol_k, epsilon_0_it, k_it)
    
    # Plot N_r, M_r
    # plt.plot(N_r, M_r, marker='o', linestyle='-', color='b')
    # plt.xlabel('N_r')
    # plt.ylabel('M_r')
    # plt.title('N_r - M_r Diagram')
    # plt.grid(True)
    # plt.show()

if __name__ == "__main__":
    main()