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
    # An initial tentative deflection is assumed, f = 0
    f = 0
    i=0
    N = F/(b*h)
    M = N*(e+f)
    Rompeu, epsilon_0, k, epsilon_0_it, k_it, f_ad = nFNC_functions.Verificacao(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N, M, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b)
    if Rompeu:
        # print('Failure due to lack of resistant capacity')
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
                # print('Failure due to lack of resistant capacity')
                return True, e, f
        
        if abs(y[-1] - f) <= tol_DF:
            # print('ELUi ok!')
            return False, e, f
        else:
            f = y[-1]
            

def Curva_viga_DF(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e):
    e = M_d/N_d
    F = N_d*b*h
    L = l_e/2
    dL = L/m
    tol_DF = 1e-10
    # An initial tentative deflection is assumed, f = 0
    f = 0
    i=0
    N = F/(b*h)
    M = N*(e+f)
    Rompeu, epsilon_0, k, epsilon_0_it, k_it, f_ad = nFNC_functions.Verificacao(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N, M, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b)
    if Rompeu:
        # print('Failure due to lack of resistant capacity')
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
                # print('Failure due to lack of resistant capacity')
                return
        
        if abs(y[-1] - f) <= tol_DF:
            # print('ELUi ok!')
            z = np.linspace(0,L,m+1)
            plt.figure(figsize=(8, 6))
            plt.plot(100*y, z, '-b', linewidth=2)
            plt.xlabel('$ y (cm) $', fontsize=12)
            plt.ylabel('$ z (m) $', fontsize=12)
            plt.xlim(max(100*y) * -1.1, max(100*y) * 1.1)
            plt.ylim(0, max(z) * 1.1)
            plt.title('Discretizacao do pilar em $m$ secoes', fontsize=12)
            plt.gca().set_aspect('equal', adjustable='box')
            plt.grid(True)
            plt.show()
            return
        else:
            f = y[-1]
    
    
def Compressao_uniforme(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e):
    N_r_epsilon_c2 = nFNC_functions.Nc(epsilon_c2, k, epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k) + nFNC_functions.Ns(y_s, nb, phi, f_yd, epsilon_yd, epsilon_c2, k)
    epsilon_0_l = 0
    epsilon_0_r = epsilon_c2
    EI_epsilon_0_i = 0
    i = 0
    while i < it_max and abs(epsilon_0_r-epsilon_0_l)/2 > 1e-9:
        EI_epislon_0_l = -(nFNC_functions.Jc(epsilon_t, epsilon_b, epsilon_0_l, epsilon_c2, sigma_cd, n, b, k, h, tol_k)[1][1] + nFNC_functions.Js(E_s, epsilon_yd, phi, y_s, nb, epsilon_0_l, k)[1][1])
        N_r_l = nFNC_functions.Nc(epsilon_0_l, k, epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k) + nFNC_functions.Ns(y_s, nb, phi, f_yd, epsilon_yd, epsilon_0_l, k)
        N_cr_l = (np.pi/l_e)**2*EI_epislon_0_l*1e3
        f_l = N_cr_l - N_r_l
        
        epsilon_0_i = (epsilon_0_l+epsilon_0_r)/2
        EI_epsilon_0_i = -(nFNC_functions.Jc(epsilon_t, epsilon_b, epsilon_0_i, epsilon_c2, sigma_cd, n, b, k, h, tol_k)[1][1] + nFNC_functions.Js(E_s, epsilon_yd, phi, y_s, nb, epsilon_0_i, k)[1][1])
        N_r_i = nFNC_functions.Nc(epsilon_0_i, k, epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k) + nFNC_functions.Ns(y_s, nb, phi, f_yd, epsilon_yd, epsilon_0_i, k)
        N_cr_i = (np.pi/l_e)**2*EI_epsilon_0_i*1e3
        f_i = N_cr_i - N_r_i

        i += 1
        if (f_i*f_l) >= 0:
            epsilon_0_l = epsilon_0_i
        else:
            epsilon_0_r = epsilon_0_i
        
    
    # if N_d >= N_r_epsilon_c2:
    if i == it_max:
        N_cr = N_r_epsilon_c2
    else:
        N_cr = (np.pi/l_e)**2*EI_epsilon_0_i*1e3
    return N_cr
        
# def Curva_de_projeto_ELUi(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e):
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
    plt.ylabel('$N_cr  (kN)$', fontsize=12)
    plt.xlim(min(e_values) * 1.1, max(e_values) * 1.1)
    plt.ylim(min(N_cr_values) * 1.1, max(N_cr_values) * 1.1)
    plt.title('Curve e/h versus N_cr', fontsize=12)
    # plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True)
    plt.show()        
    return 

def Curva_de_projeto_ELUi(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e, e):
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
    while dN/(sigma_cd*b*h) > 1e-5:
        N_d = N_d + dN
        M_d = N_d*(e)
        Rompeu, _, f = Verificacao_DF(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e)
        F = 100*f/h
        if Rompeu == False:
            N.append(N_d)
            M.append(N_d*(e+f))
        else:
            N_d = N_d - dN
            dN = dN/10
    nu = np.array(N)/(sigma_cd*b*h)
    mu = np.array(M)/(sigma_cd*b*h**2)
    plt.figure(figsize=(8, 6))
    plt.plot(nu, mu, '-b', linewidth=2)
    plt.xlabel('$ \\nu $', fontsize=12)
    plt.ylabel('$ \\mu $', fontsize=12)
    plt.xlim(min(nu) * 1.1, max(nu) * 1.1)
    plt.ylim(min(mu) * 1.1, max(mu) * 1.1)
    plt.title('Trajetoria de equilibrio', fontsize=12)
    # plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True)
    plt.show()
    N = np.array(N)
    M = np.array(M)
    # N_cr = N[-1]
    # M_cr = M[-1]
    return N, M

# def Pilar_padrao_M_i(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e, e, N_d):
    M_d = 0
    f = 0
    k_i = 0
    xi = (l_e/np.pi)**2
    dM = sigma_cd*b*h**2/10
    N = [N_d]
    M_i = [M_d]
    k = [k_i]
    _r = [k_i/1000]
    M_e = [N_d*(e + xi*k_i/1000)]
    while dM/(sigma_cd*b*h**2) > 1e-8:
        M_d = M_d + dM
        Rompeu, epsilon_0_i, k_i, _, _, _ = nFNC_functions.Verificacao(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k_i, epsilon_0_it, k_it, epsilon_t, epsilon_b)
        if Rompeu == False:
            N.append(N_d)
            M_i.append(M_d)
            k.append(k_i)
            _r.append(k_i/1000)
            M_e.append(N_d*(e + xi*k_i/1000))
        else:
            M_d = M_d - dM
            dM = dM/10
    N = np.array(N)
    M_i = np.array(M_i)
    M_e = np.array(M_e)
    _r = np.array(_r)
    k = np.array(k)
    plt.figure(figsize=(8, 6))
    plt.plot(_r, M_i, '-b', label='M_int', linewidth=2)
    plt.plot(_r, M_e, '-r', label='M_ext', linewidth=2)
    plt.xlabel('$ \\frac{1}{r} (m^{-1})$', fontsize=12)
    plt.ylabel('$ M (MN \\cdot m)$', fontsize=12)
    # plt.xlim(min(nu) * 1.1, max(nu) * 1.1)
    # plt.ylim(min(mu) * 1.1, max(mu) * 1.1)
    plt.title('Metodo do Pilar Padrao', fontsize=12)
    plt.legend()
    # plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True)
    plt.show()
    # N_cr = N[-1]
    # M_cr = M[-1]
    return M_i, M_e, _r

def Pilar_padrao_M_i(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e, N_d, e):
    M_d = 0
    i = 0
    k_i = 0
    xi = (l_e/np.pi)**2
    dM = sigma_cd*b*h**2/100
    N = [N_d]
    M_i = []
    k = []
    _r = []
    M_e = []
    while dM > sigma_cd*b*h**2/100*1e-9 and i < it_max:
        Rompeu, epsilon_0_i, k_i, _, _, _ = nFNC_functions.Verificacao(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d + dM, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k_i, epsilon_0_it, k_it, epsilon_t, epsilon_b)
        if Rompeu == False:
            M_d += dM
            N.append(N_d)
            M_i.append(M_d)
            k.append(k_i)
            _r.append(k_i/1000)
            M_e.append(N_d*(e + xi*k_i/1000))
        else:
            dM = dM/2
    N = np.array(N)
    M_i = np.array(M_i)
    M_e = np.array(M_e)
    _r = np.array(_r)
    k = np.array(k)
    plt.figure(figsize=(8, 6))
    plt.plot(_r, M_i, '-b', label='M_int', linewidth=2)
    plt.plot(_r, M_e, '-r', label='M_ext', linewidth=2)
    plt.xlabel('$ \\frac{1}{r} (m^{-1})$', fontsize=12)
    plt.ylabel('$ M (MN \\cdot m)$', fontsize=12)
    # plt.xlim(min(nu) * 1.1, max(nu) * 1.1)
    # plt.ylim(min(mu) * 1.1, max(mu) * 1.1)
    plt.title('Metodo do Pilar Padrao', fontsize=12)
    plt.legend()
    # plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True)
    plt.show()
    # N_cr = N[-1]
    # M_cr = M[-1]
    return M_i, M_e, _r

def Pilar_padrao_sol(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e, N_d, e):
    M_d = 0
    i = 0
    k_i = 0
    xi = (l_e/np.pi)**2
    dM = sigma_cd*b*h**2/1e4
    N = [N_d]
    M_i = []
    k = []
    _r = []
    M_e = []
    while dM > sigma_cd*b*h**2/100*1e-9 and i < it_max:
        Rompeu, epsilon_0_i, k_i, _, _, _ = nFNC_functions.Verificacao(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d + dM, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k_i, epsilon_0_it, k_it, epsilon_t, epsilon_b)
        if Rompeu == False:
            M_d += dM
            N.append(N_d)
            M_i.append(M_d)
            k.append(k_i)
            _r.append(k_i/1000)
            M_e.append(N_d*(e + xi*k_i/1000))
        else:
            dM = dM/2
    N = np.array(N)
    M_i = np.array(M_i)
    M_e = np.array(M_e)
    _r = np.array(_r)
    k = np.array(k)
    sol = []
    for i in range(len(M_i)):
        if abs(M_i[i]-M_e[i])/(sigma_cd*b*h**2) < 1e-4:
            sol.append([M_i[i],k[i],_r[i]])
    if sol == []:
        return None
    else:
        return sol



def Normal_critica(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e, e):
    N_c = nFNC_functions.Nc(epsilon_c2, 0, epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k)
    # M_c = nFNC_functions.Mc(epsilon_c2, 0, epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k)
    N_s = nFNC_functions.Ns(y_s, nb, phi, f_yd, epsilon_yd, epsilon_c2, 0)
    # M_s = nFNC_functions.Ms(y_s, nb, phi, f_yd, epsilon_yd, epsilon_c2, 0)
    N_d = 0
    f = 0
    dN = sigma_cd*b*h/10
    M_d = N_d*(e+f)
    N = [N_d]
    M = [M_d]
    while dN/(sigma_cd*b*h) > 1e-5:
        N_d = N_d + dN
        M_d = N_d*(e)
        Rompeu, _, f = Verificacao_DF(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e)
        F = 100*f/h
        if Rompeu == False:
            N.append(N_d)
            M.append(N_d*e)
            # M.append(N_d*(e+f))
        else:
            N_d = N_d - dN
            dN = dN/10
    nu = np.array(N)/(sigma_cd*b*h)
    mu = np.array(M)/(sigma_cd*b*h**2)
    N_cr = nu[-1]*sigma_cd*b*h
    M_cr = mu[-1]*sigma_cd*b*h**2
    return N_cr, M_cr

# fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b = Entrada_de_dados()
# m = 5
# l_e = 20*h
# e = h/3
# M_i,M_e,_r = Pilar_padrao_M_i(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e, N_d, e)
# theta = h*1e3*_r

# plt.figure(figsize=(8, 6))
# plt.plot(theta, M_i/(sigma_cd*b*h**2), '-b', label='mu_int', linewidth=2)
# plt.plot(theta, M_e/(sigma_cd*b*h**2), '-r', label='mu_ext', linewidth=2)
# plt.xlabel('$ \\theta$', fontsize=12)
# plt.ylabel('$ \\mu $', fontsize=12)
# # plt.xlim(min(nu) * 1.1, max(nu) * 1.1)
# # plt.ylim(min(mu) * 1.1, max(mu) * 1.1)
# plt.title('Metodo do Pilar Padrao', fontsize=12)
# plt.legend()
# # plt.gca().set_aspect('equal', adjustable='box')
# plt.grid(True)
# plt.show()

# sol = Pilar_padrao_sol(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e, N_d, e)
# if sol != None:
#     print('M (MNm) | k (m) | _r (m)')
#     print('------------------------')
#     for i in range(len(sol)):
#         print('{:.4f} | {:.3f} | {:.4f}'.format(sol[i][0],sol[i][1],sol[i][2]))

# print(Compressao_uniforme(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e))
# # Rompeu, e, f = Verificacao_DF(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e)
# # print(Rompeu, e, f)
# i = 1
# e = np.zeros(15)
# N_cr = np.zeros(15)
# while i <= 14:
#     e[i] = i*h/10
#     N_cr[i] = Normal_critica(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e, e[i])
#     i = i + 1

# for i in range(len(e)):
#     if i != 0:
#         print(e[i]/h,' : ', N_cr[i])

# Curva_de_projeto_ELUi(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e)

    