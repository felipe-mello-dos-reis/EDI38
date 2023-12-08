import numpy as np

def Entrada_de_dados_2021():
    # Specifications for concrete
    fck = 65  # MPa
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
    h = 0.50  # m
    d = 0.03  # m
    nc = 4
    nb = [2, 2, 2, 2]
    # A_s = 0.4*b*h*sigma_cd/f_yd
    # A_s = 0.4*b*h*sigma_cd/f_yd
    # phi = np.sqrt(A_s*4/(np.pi*sum(nb)))  # m
    phi = 0.010 # m
    y_s = np.linspace(-(h / 2 - d), (h / 2 - d), nc)

    # Applied forces
    N_d = 500*1e-3  # MN
    M_d = 5000*1e-5  # MN*m

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

    # Constants
    y_t = h / 2
    y_b = -h / 2
    epsilon_0 = 0
    k = 0
    epsilon_0_it = [epsilon_0]
    k_it = [k]

    epsilon_t = epsilon_0 + y_t * k
    epsilon_b = epsilon_0 + y_b * k

    # Return all variables
    return fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b