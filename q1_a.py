
# GEPP
# Gaussian Elimination With Partial Pivoting

import numpy as np
import math

def pivot(A, b, k):
    n = A.shape[0]
    q = max(range(k, n), key=lambda i: abs(A[i, k]))
    if q == k:
        return

    for j in range(k, n):
        A[k, j], A[q, j] = A[q, j], A[k, j]
    temp = b[k, 0]
    b[k, 0] = b[q, 0]
    b[q, 0] = temp

    ## print pivot and what was swapped
    print(f"R_{k+1} <-> R_{q+1}")
    display_augmented_matrix(A, b)

def eliminate(A, b, k, i):
    n = A.shape[1]
    m_ik = A[i, k] / A[k, k]
    if A[k, k] < 0:
        m_ik = m_ik * -1
    for j in range(k, n):
        A[i, j] = A[i, j] - m_ik * A[k, j]
        ## print row elimination
    b[i, 0] = b[i, 0] - m_ik * b[k, 0]
    print(f"R_{i+1} - {m_ik} * R_{k+1} -> R_{i+1}")


def complete_forward_elimination(A, b, k):
    n = A.shape[0]
    for i in range(k + 1, n):
        eliminate(A, b, k, i)
        ## print row elimination
    display_augmented_matrix(A, b)

def backward_substitution(A, b):
    n = A.shape[0]
    x = np.zeros((n, 1))
    for i in range(n - 1, -1, -1):
        x[i, 0] = (b[i, 0] - sum(A[i, j] * x[j, 0] for j in range(i + 1, n))) / A[i, i]   

    # print unsolved equatoins
    equations = []
    for i in range(n):
        equation = ''
        for j in range(n):
            if A[i, j] != 0 and j != n - 1:
                equation += f"{A[i, j]}(x_{j}) + "
            if A[i, j] != 0 and j == n - 1:
                equation += f"{A[i, j]}(x_{j})"
        equation += f" = {b[i, 0]}"
        print(equation)
        equations.append(equation)
    
    # start from the last, and solve for x
    # then move backwards and solve for x using the new value
    # 1. show equation
    # plug in values if there are values to be plugged in, except for the one we are solving for
    # 2. solve for x
    # 3. show solution
    print()
    for i in range(n - 1, -1, -1):
        

        print(f"{equations[i]}")
        for j in range(n):
            
            if j != i:
                equations[i] = equations[i].replace(f"(x_{j})", f"({x[j, 0]})")
        print(f"{equations[i]}")
        print(f"x_{i} = {x[i, 0]}")
        print()
        

    return x


def gepp(A, b):
  """
  A : should be an n x n matrix
  b : should be a 1 x n column vecotr
  """

  # loop over every column
  # determine the row with the largest value in the column
  # swap that row with the row that corresponds with the column
  # for each row below, we get the factor
  # the factor is the current row leading col, divided by the largest val leading col row
  # we subtract each element in the row by the element in the same col in the largest val row times the factor
  # we do the same for b


  # Steps

  # Find pivot

  # pivot

  # Get factor

  # eliminate

  # move on

  if A.shape[0] != A.shape[1]:
    raise ValueError("Matrix A must be a square matrix.")
  
  if A.shape[0] != b.shape[0]:
    raise ValueError("Matrix A and vector b must have the same number of rows.")
  
  if b.shape[1] != 1:
    raise ValueError("Vector b must be a column vector.")
  
  n = A.shape[0]
  
  
  for k in range(n - 1):
    pivot(A, b, k)
    complete_forward_elimination(A, b, k)
  x = backward_substitution(A, b)
  print("Ax = b")



  return x


      




  


def pretty_print(matrix):
    if not isinstance(matrix, np.ndarray):
        raise ValueError("Input should be a NumPy array.")
    rows, cols = matrix.shape
    v_borders = ("[", "]")
    for i in range(rows):
        print(f"{v_borders[0]:>2}", end=" ")
        for j in range(cols):
            print(f"{matrix[i, j]:8.1f}", end="  ")
        print(f"{v_borders[1]:>2}")
    print()

import numpy as np

# def display_augmented_matrix(matrix, vector):
#     if vector.ndim == 1: 
#         vector = vector[:, np.newaxis]
#     for i in range(matrix.shape[0]):
#         print("\n[", end="   ")
#         for j in range(matrix.shape[1]):
#             print(f"{matrix[i, j]:10.2f}", end="   ")
#         print(f"| {vector[i, 0]:10.2f}  ]", end="\n")

def display_augmented_matrix(matrix, vector):
    if vector.ndim == 1: 
        vector = vector[:, np.newaxis]
    rows, cols = matrix.shape
    
    max_width_matrix = max([len(f"{item:.2f}") for row in matrix for item in row])
    max_width_vector = max([len(f"{item:.2f}") for item in vector.flatten()])
    max_width = max(max_width_matrix, max_width_vector)
    print()
    for i in range(rows):
        if i == 0:
            print("⎡", end="")
        elif i == rows - 1:
            print("⎣", end="")
        else:
            print("⎢", end="")

        for j in range(cols):
            print(f"{matrix[i, j]: {max_width}.0f}", end="")
        
        print(f" |{vector[i, 0]:{max_width}.0f}", end="")

        if i == 0:
            print(" ⎤")
        elif i == rows - 1:
            print(" ⎦")
        else:
            print(" ⎥")


def display_equation(A, x, b):
    rows_A, cols_A = A.shape
    
    max_width_A = max([len(f"{item:.2f}") for row in A for item in row])
    max_width_x = max([len(f"{item:.2f}") for item in x.flatten()])
    max_width_b = max([len(f"{item:.2f}") for item in b.flatten()])
    max_width = max(max_width_A, max_width_x, max_width_b)
    
    brackets = np.array([["⎡", "⎤"], ["⎢", "⎥"], ["⎣", "⎦"]])
    
    for i in range(rows_A):
        bracket_type = 1 if 0 < i < rows_A-1 else 0 if i == 0 else 2
        
        a_row = " ".join([f"{val: {max_width}.2f}" for val in A[i, :]])
        x_val = f"{x[i, 0]: {max_width}.2f}"
        b_val = f"{b[i, 0]: {max_width}.2f}"
        
        equal_sign = " = " if i == rows_A // 2 else "   "
        
        print(f"{brackets[bracket_type, 0]}{a_row}{brackets[bracket_type, 1]} {brackets[bracket_type, 0]}{x_val}{brackets[bracket_type, 1]}{equal_sign}{brackets[bracket_type, 0]}{b_val}{brackets[bracket_type, 1]}")


def main():
    matrix_A = np.array([[-3, 0, 24, 6],
                        [3, 9, -9, 24],
                        [6, 6, -6, 24],
                        [2, 5, 1, 29]
                        ])

    vector_b = np.array([[-15],
                        [39],
                        [30],
                        [31]
                        ])

    pretty_print(matrix_A)
    pretty_print(vector_b)

    display_augmented_matrix(matrix_A, vector_b)

    display_equation(matrix_A, gepp(matrix_A, vector_b), vector_b)


main()