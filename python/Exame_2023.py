import numpy as np
import matplotlib.pyplot as plt
from math import pi
from nFNC_functions import *
from Entrada_de_dados_2023 import *
from ELUi_functions import *

fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b = Entrada_de_dados_2023()
m = 10

l_e = 18*2
N_d = 6500*1e-3
M_d = 1807000*1e-5
e = M_d/N_d
nc = 9
nb = 2*np.ones(nc)
nb[0] = 40
nb[-1] = 40
phi = 0.040 # m
y_s = np.linspace(-(h / 2 - d), (h / 2 - d), nc)

Rompeu, epsilon_0_DF, k_DF, f_DF = Verificacao_DF(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e)
print('Diferencas Finitas')
print('------------------')
print('epsilon_0: {:.5f}'.format(epsilon_0_DF))
print('k: {:.5f}'.format(k_DF))
print('f: {:.5f}'.format(f_DF))

_,_,_,_,_,_ = Pilar_padrao_M_i(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e, N_d, e)

sol = Pilar_padrao_sol(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e, N_d, e)


if sol != None:
    print('M (kNm)    | _r (m) | f (cm)  | epsilon_0  | k (m)')
    print('--------------------------------------------------')
    for i in range(len(sol)):
        print('{:.4f} | {:.4f} | {:.4f} | {:.5f}   | {:.3f}'.format(sol[i][0]*1e3,sol[i][1],sol[i][2]*1e2,sol[i][3],sol[i][4]))
else:
    print('Nao foi possivel equilibrar o pilar')
    
l_e = 25*2
N_d = 7200*1e-3
M_d = 2268000*1e-5
e = M_d/N_d
nc = 9
nb = 2*np.ones(nc)
nb[0] = 40
nb[-1] = 40
phi = 0.040 # m
y_s = np.linspace(-(h / 2 - d), (h / 2 - d), nc)


Rompeu, epsilon_0_DF, k_DF, f_DF = Verificacao_DF(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e)
print('Diferencas Finitas')
print('------------------')
print('epsilon_0: {:.5f}'.format(epsilon_0_DF))
print('k: {:.5f}'.format(k_DF))
print('f: {:.5f}'.format(f_DF))

_,_,_,_,_,_ = Pilar_padrao_M_i(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e, N_d, e)

sol = Pilar_padrao_sol(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e, N_d, e)

if sol != None:
    print('M (kNm)    | _r (m) | f (cm)  | epsilon_0  | k (m)')
    print('--------------------------------------------------')
    for i in range(len(sol)):
        print('{:.4f} | {:.4f} | {:.4f} | {:.5f}   | {:.3f}'.format(sol[i][0]*1e3,sol[i][1],sol[i][2]*1e2,sol[i][3],sol[i][4]))
else:
    print('Nao foi possivel equilibrar o pilar')
    
l_e = 25*2
N_d = 7200*1e-3
M_d = 2268000*1e-5
e = M_d/N_d
nc = 9
nb = 2*np.ones(nc)
nb[0] = 40
nb[1] = 30
nb[-1] = 40
phi = 0.040 # m
y_s = np.linspace(-(h / 2 - d), (h / 2 - d), nc)


Rompeu, epsilon_0_DF, k_DF, f_DF = Verificacao_DF(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e)
print('Diferencas Finitas')
print('------------------')
print('epsilon_0: {:.5f}'.format(epsilon_0_DF))
print('k: {:.5f}'.format(k_DF))
print('f: {:.5f}'.format(f_DF))


_,_,_,_,_,_ = Pilar_padrao_M_i(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e, N_d, e)

sol = Pilar_padrao_sol(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b, m, l_e, N_d, e)

if sol != None:
    print('M (kNm)    | _r (m) | f (cm)  | epsilon_0  | k (m)')
    print('--------------------------------------------------')
    for i in range(len(sol)):
        print('{:.4f} | {:.4f} | {:.4f} | {:.5f}   | {:.3f}'.format(sol[i][0]*1e3,sol[i][1],sol[i][2]*1e2,sol[i][3],sol[i][4]))
else:
    print('Nao foi possivel equilibrar o pilar')
    
