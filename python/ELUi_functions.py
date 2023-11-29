import numpy as np
import matplotlib.pyplot as plt
from math import pi
import nFNC_functions
from Entrada_de_dados import Entrada_de_dados

def Verificacao_DF(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e):
    e = M_d/N_d
    F = N_d*b*h
    L = l_e/2
    dL = L/m
    tol_DF = 1e-10
    # Arbitra-se uma flecha tentativa inicial f = 0
    f = 0
    i=0
    N = F/(b*h)
    M = N*(e+f)
    Rompeu, epsilon_0, k, epsilon_0_it, k_it, f_ad = nFNC_functions.Verificacao(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N, M, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b)
    if Rompeu:
        print('Failure due to lack of resistant capacity')
        return True, e, f
    while (True):
        M_d = N_d*(e+f)
        Rompeu, epsilon_0, k, epsilon_0_it, k_it, f_ad = nFNC_functions.Verificacao(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b)
        N = N_d*np.ones(m+1)
        y = np.zeros(m+1)
        M = np.zeros(m+1)
        _r = np.zeros(m+1)
        
        for i in range(m+1):
            if (i == 0):
                y[i] = 0
                _r[i] = k/1000
                M[i] = N[i]*(e+f)
            elif (i == 1):
                y[i] = _r[0]*dL**2/2
                M[i] = N[i]*(e+f-y[i])
                Rompeu, epsilon_0, k, epsilon_0_it, k_it, f_ad = nFNC_functions.Verificacao(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N[i], M[i], epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b)
                _r[i] = k/1000
            else:
                y[i] = _r[i-1]*dL**2+2*y[i-1]-y[i-2]
                M[i] = N[i]*(e+f-y[i])
                Rompeu, epsilon_0, k, epsilon_0_it, k_it, f_ad = nFNC_functions.Verificacao(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N[i], M[i], epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b)
                _r[i] = k/1000
            if Rompeu:
                print('Falha por falta de capacidade resistente da secao')
                return True, e, f
        
        if abs(y[-1] - f) <= tol_DF:
            print('ELUi ok!')
            return False, e, f
        else:
            f = y[-1]
            
            
def Curva_de_projeto_ELUi(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e):
    N_c = nFNC_functions.Nc(epsilon_c2, 0, epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k)
    # M_c = nFNC_functions.Mc(epsilon_c2, 0, epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k)
    N_s = nFNC_functions.Ns(y_s, nb, phi, f_yd, epsilon_yd, epsilon_c2, 0)
    # M_s = nFNC_functions.Ms(y_s, nb, phi, f_yd, epsilon_yd, epsilon_c2, 0)
    e_values = np.linspace(0,1,101)
    N_cr_values = np.zeros(len(e_values))
    M_cr_values = np.zeros(len(e_values))
    f_values = np.zeros(len(e_values))
    N = 0
    j = 0
    while (j < len(e_values)):
        dN = (N_c + N_s)/1000
        N = -dN
        Rompeu = False
        while (Rompeu == False) and (N <= N_c + N_s):
            N = N + dN
            M = e_values[j]*N
            Rompeu, _, f_values[j] = Verificacao_DF(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N, M, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e)
        N_cr_values[j] = N
        M_cr_values[j] = N_cr_values[j]*(e_values[j]+f_values[j])
        j = j + 1
    plt.figure(figsize=(8, 6))
    plt.plot(e_values/h, N_cr_values*1000, '-r', linewidth=2)
    plt.xlabel('$e/h$', fontsize=12)
    plt.ylabel('$N_cr \; (kN)$', fontsize=12)
    plt.xlim(min(e_values) * 1.1, max(e_values) * 1.1)
    plt.ylim(min(N_cr_values) * 1.1, max(N_cr_values) * 1.1)
    plt.title('Curve e/h versus N_cr', fontsize=12)
    # plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True)
    plt.show()        
    return 

def Normal_critica(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e, e):
    N_c = nFNC_functions.Nc(epsilon_c2, 0, epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k)
    # M_c = nFNC_functions.Mc(epsilon_c2, 0, epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k)
    N_s = nFNC_functions.Ns(y_s, nb, phi, f_yd, epsilon_yd, epsilon_c2, 0)
    # M_s = nFNC_functions.Ms(y_s, nb, phi, f_yd, epsilon_yd, epsilon_c2, 0)
    N_d = 0
    f = 0
    dN = sigma_cd*b*h/10
    passo = dN/(sigma_cd*b*h)
    M_d = N_d*(e+f)
    N = [N_d]
    M = [M_d]
    # Rompeu, _, _, _, _, _  = nFNC_functions.Verificacao(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b)
    while dN/(sigma_cd*b*h) > 1e-5:
        N_d = N_d + dN
        M_d = N_d*(e)
        # Rompeu, _, _, _, _, _  = nFNC_functions.Verificacao(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b)
        Rompeu, _, f = Verificacao_DF(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e)
        F = 100*f/h
        if Rompeu == False:
            N.append(N_d)
            M.append(N_d*(e+f))
            # M_d = N_d*(e+f)
        else:
            N_d = N_d - dN
            dN = dN/10
            passo = dN/(sigma_cd*b*h)
    nu = np.array(N)/(sigma_cd*b*h)
    mu = np.array(M)/(sigma_cd*b*h**2)
    # plt.figure(figsize=(8, 6))
    # plt.plot(nu, mu, '-r', linewidth=2)
    # plt.xlabel('$ nu $', fontsize=12)
    # plt.ylabel('$ mu $', fontsize=12)
    # plt.xlim(min(nu) * 1.1, max(nu) * 1.1)
    # plt.ylim(min(mu) * 1.1, max(mu) * 1.1)
    # plt.title('Trajetoria de equilibrio', fontsize=12)
    # # plt.gca().set_aspect('equal', adjustable='box')
    # plt.grid(True)
    # plt.show()
    N_cr = nu[-1]*sigma_cd*b*h
    M_cr = mu[-1]*sigma_cd*b*h**2
    return nu[-1]#, M_cr

fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b = Entrada_de_dados()
m = 5
l_e = 20*h
# Rompeu, e, f = Verificacao_DF(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e)
# print(Rompeu, e, f)
i = 1
e = np.zeros(15)
N_cr = np.zeros(15)
while i <= 14:
    e[i] = i*h/10
    N_cr[i] = Normal_critica(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e, e[i])
    i = i + 1

for i in range(len(e)):
    if i != 0:
        print(e[i]/h,' : ', N_cr[i])

# Curva_de_projeto_ELUi(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e)

    