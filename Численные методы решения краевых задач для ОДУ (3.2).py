import math

#решаем задачу: u'' = f(x, u), u(0) = mu1, u(l) = mu2 с шагом h
def solve_diff_equ(f, l, mu1, mu2, h):
    N = int(l / h + 0.5) #число точек + 1
    
    #разностная схема:
    # (u[i-1] - 2 u[i] + u[i]) / h^2 == f(x[i], u[i]), u[0] = mu1, u[N] = mu2
    # нелинейное уравнение; решаем итерационными методами
    # (u[i-1][s] - 2 u[i][s] + u[i][s]) == h^2 f(x[i], u[i][s-1]), u[0] = mu1, u[N] = mu2
    # про разностную схему: http://scask.ru/q_book_dig_m.php?id=124
    
    # задаем начальное приближение
    u = [0] * (N + 1)
    u[0] = mu1
    u[N] = mu2
    
    # счетчик итераций
    it = 0
    while True:
        #каждый раз решаем систему методом прогонки
        xi = dict()
        th = dict()
        xi[1] = 0
        th[1] = mu1
        
        # (u[i-1][s] - 2 u[i][s] + u[i][s]) == h^2 f(x[i], u[i][s-1]), u[0] = mu1, u[N] = mu2
        # -u[i-1][s] + 2 u[i][s] - u[i][s]) == -h^2 f(x[i], u[i][s-1])
        # alpha = beta = 1, gamma = 2, phi[i] = -h^2 f(x[i], u[i][s-1])
        #рассчитываем xi и theta по формулам из методички
        for i in range(1, N):
            xi[i + 1] = 1 / (2 - xi[i])
            th[i + 1] = (-h * h * f(i * h, u[i]) + th[i]) / (2 - xi[i])
        
        #новое приближение  
        new_u = [0] * (N + 1)
        
        new_u[N] = mu2
        for i in reversed(range(0, N)):
            new_u[i] = xi[i + 1] * new_u[i + 1] + th[i + 1]
        
        #если новое решение слабо отличается от старого (L2-норма не превышает 10 в -9 степени), задача решена
        if sum((new_u[i] - u[i]) * (new_u[i] - u[i]) for i in range(N + 1)) < 1e-9: 
            return new_u
            
        #если сходится слабо, то и так сойдет...
        if it > 100: 
            return new_u
            
        it += 1
        #новое приближение становится старым
        u = new_u.copy()

#решим уравнение
sol = solve_diff_equ(lambda xx, uu: -math.exp(uu), 1, 0, 0, 0.1)

#выведем ответ
print(sol)

#выведем ответ в нормальном виде:
for i in range(11):
    print("u(" + str(0.1 * i) + ") =", sol[i])