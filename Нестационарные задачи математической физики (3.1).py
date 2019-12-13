import math

#третья смешанная задача в полуполосе 
# \partial_t u - \partial_{x^2} u = f(x, t), 0 < x < l, 0 < t < T
# u(x, 0) = I(x)
# \partial_x u(0, t) + s1 u(0, t) = mu_1(t)
# \partial_x u(l, t) + s2 u(l, t) = mu_2(t)
# с шагом h (по X) и tau (по T)
def heat_equation(l, T, f, I, mu_1, mu_2, s1, s2, h, tau):
    solution = [] #сетка решения
    
    #число слоев по X
    nx = int(l / h) + 1
    #число слоев по T
    nt = int(T / tau) + 1
    
    solution.append([I(i * h) for i in range(nx)]) #заполняем начальное условие
    
    #цикл по слоям
    for i in range(1, nt):
        line = [0.0] * nx
        
        #серединные точки по разностной схеме:
        for j in range(1, nx - 1):
            S = tau / (h * h)
            line[j] = S * solution[-1][j + 1] + (1 - 2 * S) * solution[-1][j] + S * solution[-1][j - 1] + tau * f(j * h, i * tau)
        
        #крайние типа тоже:
        # \partial_x u(0, t) + s1 u(0, t) = mu_1(t)
        # (u[1] - u[0]) / h + s1 * u[0] = mu_1(t)
        # u[1] = (mu1*h - u1)/(-1 + s1*h)
        
        # \partial_x u(1, t) + s2 u(0, t) = mu_2(t)
        # (u[-1] - u[-2]) / h + s2 * u[-1] = mu_2(t)
        # u[-1] = (mu2*h + u[-2])/(1 + s2*h)
        
        line[0] = (mu_1(tau * i) * h - line[1]) / (-1 + s1 * h)
        line[-1] = (mu_2(tau * i) * h + line[-2]) / (1 + s2 * h)
        solution.append(line)
        
    return solution
    
#точное решение
def exact_solution(x, t):
    return math.sin(t * (1 + 2 * x - 3 * x * x))

#правая часть (f)   
def op_exact(x, t):
    return (1 + 6 * t + (2 - 3 * x) * x) * math.cos(t + t * (2 - 3 * x) * x) + 4 * t * t * (1 - 3 * x) * (1 - 3 * x) * math.sin(t + t * (2 - 3 * x) * x)

#частная производная по x от точного решения
def dx_exact_solution(x, t):
    return t * (2 - 6 * x) * math.cos(t + t * (2 - 3 * x) * x)
    
#считает ошибку
def calculate_error(mat1, mat2):
    err = 0.0
    for i in range(len(mat1)):
        for j in range(len(mat1[0])):
            err = max(err, abs(mat1[i][j] - mat2[i][j]))
    return err
    
#заполняет матрицу точного решения
#0 < x < l, 0 < t < T
# с шагом h и tau
# ее следовало бы в L2 считать, но посчитаем в C
def fill_exact_matrix(l, T, h, tau):
    nx = int(l / h) + 1
    nt = int(T / tau) + 1
    
    return [[exact_solution(j * h, i * tau) for j in range(nx)] for i in range(nt)]

#схема устойчива при tau / (h*h) < 1/2 
#так, если h = 0.1, то tau <= 0.005
#тестирование: h = 0.1, tau = 0.005
#              h = 0.1, tau = 0.0025
#              h = 0.05, tau = 0.00125

#устанавливаем параметры:
sigma1 = sigma2 = 0
#задаем граничные условия третьего рода:
m1 = lambda t: dx_exact_solution(0, t) + sigma1 * exact_solution(0, t)
m2 = lambda t: dx_exact_solution(1, t) + sigma2 * exact_solution(1, t)

print("sigma1 = sigma2 =", sigma1)
#решаем численно h = 0.1, tau = 0.005
num = heat_equation(1.0, 0.1, op_exact, lambda x: exact_solution(x, 0), m1, m2, sigma1, sigma2, 0.1, 0.005)
exa = fill_exact_matrix(1.0, 0.1, 0.1, 0.005)
print("h = 0.1, tau = 0.005, error:", calculate_error(num, exa))
#решаем численно h = 0.1, tau = 0.0025
num = heat_equation(1.0, 0.1, op_exact, lambda x: exact_solution(x, 0), m1, m2, sigma1, sigma2, 0.1, 0.0025)
exa = fill_exact_matrix(1.0, 0.1, 0.1, 0.0025)
print("h = 0.1, tau = 0.0025, error:", calculate_error(num, exa))
#решаем численно h = 0.05, tau = 0.00125
num = heat_equation(1.0, 0.1, op_exact, lambda x: exact_solution(x, 0), m1, m2, sigma1, sigma2, 0.05, 0.00125)
exa = fill_exact_matrix(1.0, 0.1, 0.05, 0.00125)
print("h = 0.05, tau = 0.00125, error:", calculate_error(num, exa))

#устанавливаем параметры:
sigma1 = sigma2 = 10
#задаем граничные условия третьего рода:
m1 = lambda t: dx_exact_solution(0, t) + sigma1 * exact_solution(0, t)
m2 = lambda t: dx_exact_solution(1, t) + sigma2 * exact_solution(1, t)

print("sigma1 = sigma2 =", sigma1)
#решаем численно h = 0.2, tau = 0.005
num = heat_equation(1.0, 0.1, op_exact, lambda x: exact_solution(x, 0), m1, m2, sigma1, sigma2, 0.2, 0.005)
exa = fill_exact_matrix(1.0, 0.1, 0.2, 0.005)
print("h = 0.2, tau = 0.005, error:", calculate_error(num, exa))
#решаем численно h = 0.2, tau = 0.0025
num = heat_equation(1.0, 0.1, op_exact, lambda x: exact_solution(x, 0), m1, m2, sigma1, sigma2, 0.2, 0.0025)
exa = fill_exact_matrix(1.0, 0.1, 0.2, 0.0025)
print("h = 0.2, tau = 0.0025, error:", calculate_error(num, exa))


#решаем численно h = 0.01, tau = 0.00125
#num = heat_equation(1.0, 0.1, op_exact, lambda x: exact_solution(x, 0), m1, m2, sigma1, sigma2, 0.05, 0.00125)
#exa = fill_exact_matrix(1.0, 0.1, 0.05, 0.00125)
#print("h = 0.05, tau = 0.00125, error:", calculate_error(num, exa))
#ерунда с этими параметрами почему то...



#устанавливаем параметры:
sigma1 = sigma2 = 100
#задаем граничные условия третьего рода:
m1 = lambda t: dx_exact_solution(0, t) + sigma1 * exact_solution(0, t)
m2 = lambda t: dx_exact_solution(1, t) + sigma2 * exact_solution(1, t)

print("sigma1 = sigma2 =", sigma1)
#решаем численно h = 0.1, tau = 0.005
num = heat_equation(1.0, 0.1, op_exact, lambda x: exact_solution(x, 0), m1, m2, sigma1, sigma2, 0.1, 0.005)
exa = fill_exact_matrix(1.0, 0.1, 0.1, 0.005)
print("h = 0.1, tau = 0.005, error:", calculate_error(num, exa))
#решаем численно h = 0.1, tau = 0.0025
num = heat_equation(1.0, 0.1, op_exact, lambda x: exact_solution(x, 0), m1, m2, sigma1, sigma2, 0.1, 0.0025)
exa = fill_exact_matrix(1.0, 0.1, 0.1, 0.0025)
print("h = 0.1, tau = 0.0025, error:", calculate_error(num, exa))
#решаем численно h = 0.05, tau = 0.00125
num = heat_equation(1.0, 0.1, op_exact, lambda x: exact_solution(x, 0), m1, m2, sigma1, sigma2, 0.05, 0.00125)
exa = fill_exact_matrix(1.0, 0.1, 0.05, 0.00125)
print("h = 0.05, tau = 0.00125, error:", calculate_error(num, exa))