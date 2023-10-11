# GEPP for matrix operations
# AX = B , where A is a n x n square matrix, X is an n x p matrix , and B is an n x p matrix


import numpy as np

def pivot(A, k, P):
    n = A.shape[0]
    q = max(range(k, n), key=lambda i: abs(A[i, k]))
    if q == k:
        return

    for j in range(k, n):
        A[k, j], A[q, j] = A[q, j], A[k, j]
        P[k, j], P[q, j] = P[q, j], P[k, j]

    ## print pivot and what was swapped
    print(f"R_{k+1} <-> R_{q+1}")

def eliminate(A, k, i, L):
    n = A.shape[1]
    m_ik = A[i, k] / A[k, k]
    if A[k, k] < 0:
        m_ik = m_ik * -1
    L[i,k] = m_ik
    for j in range(k, n):
        A[i, j] = A[i, j] - m_ik * A[k, j]
        ## print row elimination


def complete_forward_elimination(U, k, L):
    n = U.shape[0]
    for i in range(k + 1, n):
        eliminate(U, k, i, L)


def LUP_factorization(A):
    n = A.shape[0]
    # initialize L to identity matrix
    # initialize P to identity matrix
    # copy A to U
    L = np.identity(n)
    P = np.identity(n)
    U = A.copy()
    # use gaussian elimination to find L and U
    # L is a lower triangular matrix with 1s on the diagonal
    # the lower part is the multipliers used in gaussian elimination
    # U is an upper triangular matrix with the values after gaussian elimination
    for k in range(n - 1):
        pivot(U, k, P)
        complete_forward_elimination(U, k, L)
    return L, U, P

def backward_substitution(A, b):
    n = A.shape[0]
    x = np.zeros((n, 1))
    for i in range(n - 1, -1, -1):
        x[i, 0] = (b[i, 0] - sum(A[i, j] * x[j, 0] for j in range(i + 1, n))) / A[i, i] 
    return x

def forward_substitution(A, b):
    n = A.shape[0]
    x = np.zeros((n, 1))
    for i in range(n):
        x[i, 0] = (b[i, 0] - sum(A[i, j] * x[j, 0] for j in range(i))) / A[i, i] 
    return x

    

# general gaussian elimination with partial pivoting
def ggepp(A,B):
    n = A.shape[0]
    # initialize L to identity matrix
    # initialize P to identity matrix
    # copy A to U
    L = np.identity(n)
    P = np.identity(n)
    U = A.copy()

    p = B.shape[1]
    # use gaussian elimination to find L and U
    # L is a lower triangular matrix with 1s on the diagonal
    # the lower part is the multipliers used in gaussian elimination
    # U is an upper triangular matrix with the values after gaussian elimination
    for k in range(n - 1):
        pivot(U, k, P)
        complete_forward_elimination(U, k, L)
    # solve for y in Ly = Pb for all columns in B
    Y = np.zeros((n, p))
    X = np.zeros((n, p))
    for i in range(p):
        y = forward_substitution(L, np.matmul(P, B[:,i:i+1]))
        # solve for x in Ux = y
        x = backward_substitution(U, y)
        Y[:,i:i+1] = y
        X[:,i:i+1] = x
    return X
        

def pretty_print(matrix):
    if not isinstance(matrix, np.ndarray):
        raise ValueError("Input should be a NumPy array.")
    rows, cols = matrix.shape
    v_borders = ("[", "]")
    for i in range(rows):
        print(f"{v_borders[0]:>2}", end=" ")
        for j in range(cols):
            print(f"{matrix[i, j]:2.5f}", end="  ")
        print(f"{v_borders[1]:>2}")
    print()


def hilbert_matrix(n):
    H = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            H[i, j] = 1 / (i + j + 1)
    return H

def random_matrix(n, p):
    return np.random.randn(n, p)

def main():
    matrix_A = np.array([[-3, 0, 24, 6],
                        [3, 9, -9, 24],
                        [6, 6, -6, 24],
                        [2, 5, 1, 29]
                        ])
    

    # matrix_B = np.array([[-6, 0, 24, 6],
    #                 [3, 9, -9, 24],
    #                 [6, 3, -5, 24],
    #                 [2, 3, 1, 4]
    #                 ])
    
    matrix_B = np.array([[-15, 1],
                    [39, 0],
                    [30, 0],
                    [31, 0]
                    ])

    print("A:")
    pretty_print(matrix_A)

    L, U, P = LUP_factorization(matrix_A)
    print("L:")
    pretty_print(L)
    print("U:")
    pretty_print(U)
    print("P:")
    pretty_print(P)

    # solve
    # Ax = b
    X = ggepp(matrix_A, matrix_B)


    H = hilbert_matrix(10)

    print("H:")
    pretty_print(H)

    X_t = random_matrix(10, 4)

    print("X_t:")
    pretty_print(X_t)

    B = np.matmul(H, X_t)
    
    print("B:")
    pretty_print(B)

    X_c = ggepp(H, B)

    print("X_c:")
    pretty_print(X_c)



main()

