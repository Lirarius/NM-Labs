#Метод Якоби

import numpy as np
import math
from timeit import default_timer as timer

#решение системы A x = b методом Якоби c точностью eps
def jacobi(A, b, eps):
    #начальное приближение
    x = [1] * len(b)
    
    #число итераций
    it = 0
    while True:
        #новое приближение
        xn = [0] * len(b)
        #проводим расчет
        for i in range(len(x)):
            s1 = 0
            s2 = 0
            
            for j in range(i):
                s1 += A[i][j] * x[j]
            
            for j in range(i + 1, len(b)):
                s2 += A[i][j] * x[j]
            
            xn[i] = (b[i] - s1 - s2) / A[i][i]
        
        #если массивы мало отличаются (точность достигнута), возвращаем результат
        if (math.sqrt(sum([(xn[i] - x[i]) * (xn[i] - x[i]) for i in range(len(b))])) < eps):
            return xn
            
        it += 1
        #много итераций - сходимость плохая.
        if (it > 100):
            return xn
            
        x = xn

#решение системы A x = b методом Якоби c точностью eps и векторизацией
def jacobi_vec(A, b, eps):
    x = np.array([1] * len(b))
    U = np.triu(A, k = 1) #возьмем верхедиагональную матрицу у A
    L = np.tril(A, k = -1) #возьмем нижнедиагональную матрицу у A
    D = np.diag(A) #возьмем диагональ матрицы A
    D = np.array([1 / d for d in D]) #обратим числа
    D = np.diagflat(D) #возьмем матрицу обратную к диагональной матрице из A
    
    it = 0
    
    while True:
        xn = np.dot(D, b - np.dot(L + U, x)) #векторизуем
        #если массивы мало отличаются (точность достигнута), возвращаем результат
        if (np.linalg.norm(xn - x) < eps):
            return xn
        
        it += 1
        #много итераций - сходимость плохая.
        if (it > 100):
            return xn
            
        x = xn
        

#выберем матрицу, (осторожно: не сходится для любой матрицы)
A = [[1, 0.2, 0.3], [0.3, 1, 0.06], [0.07, 0.08, 1]]
#выберем столбец правых членов
b = [1, 2, 3]
#решим систему методом Якоби с точностью 1e-6
sol = jacobi(A, b, 1e-6)
print("solution", sol)
#проверим, (должно выть меньше 1e-6)
print(np.linalg.norm(np.dot(np.matrix(A), np.array(sol)) - b), ("bad!", "ok!")[np.linalg.norm(np.dot(np.matrix(A), np.array(sol)) - b) < 1e-6])

#решим систему методом Якоби с точностью 1e-6 и векторизацией
sol = list(jacobi_vec(np.matrix(A), np.array(b), 1e-6))
print("solution", sol)
#проверим (должно выть меньше 1e-6)
print(np.linalg.norm(np.dot(np.matrix(A), np.array(sol)) - b), ("bad!", "ok!")[np.linalg.norm(np.dot(np.matrix(A), np.array(sol)) - b) < 1e-6])

#сравним время работы (10000 раз решим систему)
start = timer()
for i in range(10000):
    sol = jacobi(A, b, 1e-6)
end = timer()

loops = end - start

start = timer()
#преобразуем в numpy-типы
AA = np.matrix(A)
bb = np.array(b)
for i in range(10000):
    sol = jacobi_vec(AA, bb, 1e-6)
end = timer()

vectorization = end - start

print("time of loops", loops)
print("time of vectorization", vectorization)