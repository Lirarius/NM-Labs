#Решение полной задачи на собственные значения

import numpy as np

def qr(mat):

    proj = lambda w, v: np.dot(w, v) / np.dot(w, w) * w
    c = np.transpose(mat) 
    u = []
    for k in range(len(mat)):
        u.append(np.copy(c[k]))
        for j in range(k):
            u[-1] -= proj(u[j], c[k])
            
    e = [t / np.linalg.norm(t) for t in u]
    Q = np.transpose(np.matrix(e))
    R = np.zeros([len(mat), len(mat)])
    
    for i in range(len(mat)):
        for j in range(i, len(mat)):
            R[i][j] = np.dot(e[i], c[j])
            
    return Q, R

def eigenvalues_qr(mat):

    a = np.copy(mat)
    qq = np.identity(len(a))
    
    while True:

        q, r = qr(a)
        ak = np.squeeze(np.asarray(np.dot(r, q)))
        qq = np.dot(qq, q)
        
        if np.linalg.norm(np.diag(a) - np.diag(ak)) < 1e-15:
            return np.diag(ak) 
            
        a = ak
        
def generate_matrix(n, alpha):

    a = np.zeros([n + 1, n + 1])
    a[0][0] = 2
    a[0][1] = -1 - alpha
    a[n - 1][n] = -1 + alpha
    a[n][n] = 2
    for i in range(1, n):
        a[i][i] = 2
        a[i][i + 1] = -1 - alpha
        a[i][i - 1] = -1 + alpha
    return a

for al in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:

    test_matrix = generate_matrix(10, al)
    test_eigenvalues = sorted(eigenvalues_qr(test_matrix))
    np_eigenvalues, _ = np.linalg.eig(test_matrix)
    np_eigenvalues = sorted(np_eigenvalues)
    print("alpha:", al)
    print("eigenvalues:", test_eigenvalues)
    print("eigenvalues (numpy):", np_eigenvalues)
    print("difference:", np.linalg.norm(np.array(np_eigenvalues) - np.array(test_eigenvalues)))