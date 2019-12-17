#Метод Зейделя

import numpy as np

def seidel(A, b, eps):
    
    x = np.array([0] * len(b))
    U = np.triu(A, k = 1)
    L = np.tril(A, k = -1)
    D = np.diagflat(np.diag(A))
    
    i = 0
    while True:
        y = np.dot(np.linalg.inv(L + D), np.dot(-U, x) + b)
        if (np.linalg.norm(x - y) < eps):
            return y
        x = y
        
        i += 1
        if (i > 200):
            return x
            
A = np.array([[1, 0.3, 0.1], [0.2, 1, 0.1], [0.05, 0.05, 1]])
b = np.array([1, 1, 1])
x = seidel(A, b, 1e-9)

print("solution:", x)
print("check: ||A x - b|| =", np.linalg.norm(np.dot(A, x) - b))