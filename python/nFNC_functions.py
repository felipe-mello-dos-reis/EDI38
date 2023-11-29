import matplotlib.pyplot as plt
import numpy as np


def plot_epsilon_0_k(epsilon_c2, epsilon_cu, h, y_t, y_b, y_s, epsilon_0_it, k_it):
    c = epsilon_c2
    u = epsilon_cu
    A = [0, c]
    B = [u / h, u * y_t / h]
    C = [(u + 10) / (y_t - min(y_s)), u - y_t * ((u + 10) / (y_t - min(y_s)))]
    D = [0, -10]
    E = [(u + 10) / (y_b - max(y_s)), u - y_b * ((u + 10) / (y_b - max(y_s)))]
    F = [-u / h, -u * y_b / h]
    points = np.array([A, B, C, D, E, F, A]).T
    plt.figure(figsize=(8, 6))
    plt.plot(points[0], points[1], '-b', linewidth=2)
    plt.xlabel('$\kappa \; (1/m)$', fontsize=12)
    plt.ylabel('$\epsilon_0$', fontsize=12)
    plt.xlim(1.1 * np.min(points[0]), 1.1 * np.max(points[0]))
    plt.ylim(1.1 * np.min(points[1]), 1.1 * np.max(points[1]))
    plt.title('Feasible region for design', fontsize=12)
    # plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True)
    # plt.hold(True)
    plt.plot(k_it, epsilon_0_it, '-g')
    plt.plot(k_it[-1], epsilon_0_it[-1], '-gx')
    # plt.hold(False)
    plt.show()

def plot_N_r_M_r(epsilon_c2, epsilon_cu, sigma_cd, n, b, f_yd, epsilon_yd, h, y_t, y_b, y_s, phi, nb, tol_k, epsilon_0_it, kappa_it):
    c = epsilon_c2
    u = epsilon_cu
    A = [0, c]
    B = [u / h, u * y_t / h]
    C = [(u + 10) / (y_t - min(y_s)), u - y_t * ((u + 10) / (y_t - min(y_s)))]
    D = [0, -10]
    E = [(u + 10) / (y_b - max(y_s)), u - y_b * ((u + 10) / (y_b - max(y_s)))]
    F = [-u / h, -u * y_b / h]
    points = np.array([A, B, C, D, E, F, A]).T

    num_points = 100
    kappa_values = []
    epsilon_0_values = []
    for i in range(6):
        kappa_values = np.concatenate([kappa_values, np.linspace(points[0, i], points[0, i + 1], num_points)])
        epsilon_0_values = np.concatenate([epsilon_0_values, np.linspace(points[1, i], points[1, i + 1], num_points)])

    N_r_values = np.zeros(len(kappa_values))
    M_r_values = np.zeros(len(epsilon_0_values))

    for j in range(len(kappa_values)):
        N_c = Nc(epsilon_0_values[j], kappa_values[j], epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k)
        M_c = Mc(epsilon_0_values[j], kappa_values[j], epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k)
        N_s = Ns(y_s, nb, phi, f_yd, epsilon_yd, epsilon_0_values[j], kappa_values[j])
        M_s = Ms(y_s, nb, phi, f_yd, epsilon_yd, epsilon_0_values[j], kappa_values[j])
        N_r_values[j] = (N_c + N_s) * 1e3
        M_r_values[j] = (M_c + M_s) * 1e3


    N_r_values_it = np.zeros(len(kappa_it))
    M_r_values_it = np.zeros(len(epsilon_0_it))

    for j in range(len(epsilon_0_it)):
        N_c = Nc(epsilon_0_it[j], kappa_it[j], epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k)
        M_c = Mc(epsilon_0_it[j], kappa_it[j], epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k)
        N_s = Ns(y_s, nb, phi, f_yd, epsilon_yd, epsilon_0_it[j], kappa_it[j])
        M_s = Ms(y_s, nb, phi, f_yd, epsilon_yd, epsilon_0_it[j], kappa_it[j])
        N_r_values_it[j] = (N_c + N_s) * 1e3
        M_r_values_it[j] = (M_c + M_s) * 1e3
        
    plt.figure(figsize=(8, 6))
    plt.plot(N_r_values, M_r_values, '-r', linewidth=2)
    plt.plot(N_r_values_it, M_r_values_it, '-g', linewidth=2)
    plt.plot(N_r_values_it[-1], M_r_values_it[-1], '-gx', linewidth=2)
    plt.xlabel('$N_r \; (kN)$', fontsize=12)
    plt.ylabel('$M_r \; (kNm)$', fontsize=12)
    plt.xlim(min(N_r_values) * 1.1, max(N_r_values) * 1.1)
    plt.ylim(min(M_r_values) * 1.1, max(M_r_values) * 1.1)
    plt.title('Resistant forces', fontsize=12)
    # plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True)
    plt.show()

def N_r_M_r(epsilon_c2, epsilon_cu, sigma_cd, n, b, f_yd, epsilon_yd, h, y_t, y_b, y_s, phi, nb, tol_k):
    c = epsilon_c2
    u = epsilon_cu
    A = [0, c]
    B = [u / h, u * y_t / h]
    C = [(u + 10) / (y_t - min(y_s)), u - y_t * ((u + 10) / (y_t - min(y_s)))]
    D = [0, -10]
    E = [(u + 10) / (y_b - max(y_s)), u - y_b * ((u + 10) / (y_b - max(y_s)))]
    F = [-u / h, -u * y_b / h]
    points = np.array([A, B, C, D, E, F, A]).T

    num_points = 100
    kappa_values = []
    epsilon_0_values = []
    for i in range(6):
        kappa_values = np.concatenate([kappa_values, np.linspace(points[0, i], points[0, i + 1], num_points)])
        epsilon_0_values = np.concatenate([epsilon_0_values, np.linspace(points[1, i], points[1, i + 1], num_points)])

    N_r_values = np.zeros(len(kappa_values))
    M_r_values = np.zeros(len(epsilon_0_values))

    for j in range(len(kappa_values)):
        N_c = Nc(epsilon_0_values[j], kappa_values[j], epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k)
        M_c = Mc(epsilon_0_values[j], kappa_values[j], epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k)
        N_s = Ns(y_s, nb, phi, f_yd, epsilon_yd, epsilon_0_values[j], kappa_values[j])
        M_s = Ms(y_s, nb, phi, f_yd, epsilon_yd, epsilon_0_values[j], kappa_values[j])
        N_r_values[j] = (N_c + N_s)
        M_r_values[j] = (M_c + M_s)


    # N_r_values_it = np.zeros(len(kappa_it))
    # M_r_values_it = np.zeros(len(epsilon_0_values))

    # for j in range(len(epsilon_0_values)):
    #     N_c = Nc(epsilon_0_values[j], kappa_values[j], epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k)
    #     M_c = Mc(epsilon_0_values[j], kappa_values[j], epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k)
    #     N_s = Ns(y_s, nb, phi, f_yd, epsilon_yd, epsilon_0_values[j], kappa_values[j])
    #     M_s = Ms(y_s, nb, phi, f_yd, epsilon_yd, epsilon_0_values[j], kappa_values[j])
    #     N_r_values_it[j] = (N_c + N_s) * 1e3
    #     M_r_values_it[j] = (M_c + M_s) * 1e3
    return N_r_values, M_r_values
     

def I0(epsilon, epsilon_c2, sigma_cd, n):
    if epsilon < 0:
        I_0 = 0
    elif 0 <= epsilon <= epsilon_c2:
        I_0 = sigma_cd * (epsilon - epsilon_c2 * (1 - (1 - epsilon / epsilon_c2) ** (n + 1)) / (n + 1))
    else:
        I_0 = sigma_cd * (epsilon - epsilon_c2 / (n + 1))
    return I_0

def I1(epsilon, epsilon_c2, sigma_cd, n):
    if epsilon < 0:
        I_1 = 0
    elif 0 <= epsilon <= epsilon_c2:
        I_1 = sigma_cd * (epsilon**2 / 2 + epsilon_c2**2 * (
                -1 * (1 - (1 - epsilon / epsilon_c2) ** (n + 1)) / (n + 1) + 1 * (1 - (1 - epsilon / epsilon_c2) ** (n + 2)) / (n + 2)))
    else:
        I_1 = sigma_cd * (epsilon**2 / 2 - epsilon_c2**2 / ((n + 1) * (n + 2)))
    return I_1

def DI0(epsilon_t, epsilon_b, epsilon_c2, sigma_cd, n):
    return I0(epsilon_t, epsilon_c2, sigma_cd, n) - I0(epsilon_b, epsilon_c2, sigma_cd, n)

def DI1(epsilon_t, epsilon_b, epsilon_c2, sigma_cd, n):
    return I1(epsilon_t, epsilon_c2, sigma_cd, n) - I1(epsilon_b, epsilon_c2, sigma_cd, n)

def DJm(epsilon_t, epsilon_b, epsilon_0, epsilon_c2, sigma_cd, n, m):
    return Jm(epsilon_t, epsilon_0, epsilon_c2, sigma_cd, n, m) - Jm(epsilon_b, epsilon_0, epsilon_c2, sigma_cd, n, m)

def Jm(epsilon, epsilon_0, epsilon_c2, sigma_cd, n, m):
    return sigmac(epsilon, epsilon_c2, sigma_cd, n) * (epsilon - epsilon_0)**m

def Dc(epsilon, epsilon_c2, sigma_cd, n):
    if epsilon < 0:
        D_c = 0
    elif 0 <= epsilon <= epsilon_c2:
        D_c = sigma_cd * n * (1 - epsilon / epsilon_c2)**(n - 1) / epsilon_c2
    else:
        D_c = 0
    return D_c

def Nc(epsilon_0, k, epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k):
    if abs(k * h) < tol_k:
        N_c = sigmac(epsilon_0, epsilon_c2, sigma_cd, n) * b * h
    else:
        epsilon_b = epsilon_0 + y_b * k
        epsilon_t = epsilon_0 + y_t * k
        N_c = b * DI0(epsilon_t, epsilon_b, epsilon_c2, sigma_cd, n) / k
    return N_c

def Mc(epsilon_0, k, epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k):
    if abs(k * h) < tol_k:
        M_c = sigmac(epsilon_0, epsilon_c2, sigma_cd, n) * S(b, y_b, y_t)
    else:
        epsilon_b = epsilon_0 + y_b * k
        epsilon_t = epsilon_0 + y_t * k
        M_c = b * (DI1(epsilon_t, epsilon_b, epsilon_c2, sigma_cd, n) - epsilon_0 * DI0(epsilon_t, epsilon_b, epsilon_c2, sigma_cd, n)) / k**2
    return M_c

def sigmac(epsilon, epsilon_c2, sigma_cd, n):
    if epsilon < 0:
        sigma_c = 0
    elif 0 <= epsilon <= epsilon_c2:
        sigma_c = sigma_cd * (1 - (1 - epsilon / epsilon_c2)**n)
    else:
        sigma_c = sigma_cd
    return sigma_c

def Ns(y_s, nb, phi, f_yd, epsilon_yd, epsilon_0, k):
    A_si = [np.pi * phi**2 * nbi / 4 for nbi in nb]
    epsilon_s = epsilon_0 + y_s * k
    N_si = np.zeros_like(y_s)
    for i in range(len(y_s)):
        N_si[i] = A_si[i] * sigmas(epsilon_s[i], f_yd, epsilon_yd)
    return np.sum(N_si)

def Ms(y_s, nb, phi, f_yd, epsilon_yd, epsilon_0, k):
    A_si = [np.pi * phi**2 * nbi / 4 for nbi in nb]
    epsilon_s = epsilon_0 + y_s * k
    M_si = np.zeros_like(y_s)
    for i in range(len(y_s)):
        M_si[i] = A_si[i] * sigmas(epsilon_s[i], f_yd, epsilon_yd) * y_s[i]
    return np.sum(M_si)

def sigmas(epsilon_s, f_yd, epsilon_yd):
    if epsilon_s < -epsilon_yd:
        sigma_s = -f_yd
    elif -epsilon_yd <= epsilon_s <= epsilon_yd:
        sigma_s = f_yd / epsilon_yd * epsilon_s
    else:
        sigma_s = f_yd
    return sigma_s

def Js(E_s, epsilon_yd, phi, y_s, nb, epsilon_0, k):
    A_si = [np.pi * phi**2 * nbi / 4 for nbi in nb]
    epsilon_s = epsilon_0 + y_s * k
    J_s = np.zeros((2, 2))
    for i in range(len(A_si)):
        if epsilon_s[i] < -epsilon_yd:
            D_si = 0
        elif -epsilon_yd <= epsilon_s[i] <= epsilon_yd:
            D_si = E_s
        else:
            D_si = 0
        J_s -= D_si * A_si[i] * np.array([[1, y_s[i]], [y_s[i], y_s[i]**2]])
    return J_s

def Jc(epsilon_t, epsilon_b, epsilon_0, epsilon_c2, sigma_cd, n, b, k, h, tol_k):
    J_c = np.zeros((2, 2))
    if abs(k * h) < tol_k:
        J_c[0, 0] = -Dc(epsilon_0, epsilon_c2, sigma_cd, n) * b * h
        J_c[0, 1] = 0
        J_c[1, 0] = 0
        J_c[1, 1] = -Dc(epsilon_0, epsilon_c2, sigma_cd, n) * b * h**3 / 12
    else:
        J_c[0, 0] = -b * DJm(epsilon_t, epsilon_b, epsilon_0, epsilon_c2, sigma_cd, n, 0) / k
        J_c[0, 1] = -b * (DJm(epsilon_t, epsilon_b, epsilon_0, epsilon_c2, sigma_cd, n, 1) - DI0(epsilon_t, epsilon_b, epsilon_c2, sigma_cd, n)) / k**2
        J_c[1, 0] = -b * (DJm(epsilon_t, epsilon_b, epsilon_0, epsilon_c2, sigma_cd, n, 1) - DI0(epsilon_t, epsilon_b, epsilon_c2, sigma_cd, n)) / k**2
        J_c[1, 1] = -b * (DJm(epsilon_t, epsilon_b, epsilon_0, epsilon_c2, sigma_cd, n, 2) - 2 * (DI1(epsilon_t, epsilon_b, epsilon_c2, sigma_cd, n) - epsilon_0 * DI0(epsilon_t, epsilon_b, epsilon_c2, sigma_cd, n))) / k**3
    return J_c

def S(b, y_b, y_t):
    return b * (y_t**2 - y_b**2) / 2

def ELU(epsilon_0, k, y_b, y_t, y_s, epsilon_c2, epsilon_cu):
    epsilon_b = epsilon_0 + y_b * k
    epsilon_t = epsilon_0 + y_t * k
    epsilon_s = epsilon_0 + y_s * k
    epsilon_max = max(epsilon_b, epsilon_t)
    epsilon_min = min(epsilon_b, epsilon_t)
    epsilon_s_min = min(epsilon_s)
    if epsilon_max > epsilon_cu:
        # print("ELU ultrapassado: polo epsilon_cu!")
        return True    
    if epsilon_s_min < -10:
        # print("ELU ultrapassado: polo epsilon_s!")
        return True
    if epsilon_max - (epsilon_cu - epsilon_c2) * (epsilon_max - epsilon_min) / epsilon_cu > epsilon_c2:
        # print("ELU ultrapassado: polo epsilon_c2!")
        return True
    if (epsilon_max <= epsilon_cu) and (epsilon_s_min >= -10) and (epsilon_max - (epsilon_cu - epsilon_c2) * (epsilon_max - epsilon_min) / epsilon_cu <= epsilon_c2):
        # print("ELU ok!")
        return False
		
def Verificacao(fck, gamma_c, sigma_cd, gamma_conc, f_yk, gamma_s, E_s, f_yd, epsilon_yd, gamma_aco, c, b, h, d, nc, nb, phi, y_s, N_d, M_d, epsilon_c2, epsilon_cu, x_lim, n, tol_J, tol_k, tol_f, i, it_max, y_t, y_b, epsilon_0, k, epsilon_0_it, k_it, epsilon_t, epsilon_b):
	## Newton-Raphson method
	# Calculate N, M for first iteration

	N_s = Ns(y_s, nb, phi, f_yd, epsilon_yd, epsilon_0, k)
	M_s = Ms(y_s, nb, phi, f_yd, epsilon_yd, epsilon_0, k)

	N_c = Nc(epsilon_0, k, epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k)
	M_c = Mc(epsilon_0, k, epsilon_c2, sigma_cd, n, y_b, y_t, b, h, tol_k)

	N_r = N_c + N_s
	M_r = M_c + M_s
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
			# print("No solution exists!")
			return True, epsilon_0, k, epsilon_0_it, k_it, f_ad

	return ELU(epsilon_0, k, y_b, y_t, y_s, epsilon_c2, epsilon_cu), epsilon_0, k, epsilon_0_it, k_it, f_ad


