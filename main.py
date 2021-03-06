import numpy as np
from numpy import sqrt, arctan, sin, cos
from numpy.linalg import inv
from matplotlib import pyplot as plt

from euler import euler
from functions import y, f1, g1, f2, g2, f3, g3, derivative
from visual import draw_func

def gen_c_matrix(f1, g1, f2, g2, y0, z0, a, b, h):
    n = int((b - a) / h) + 1
    m = 2
    C = [[], []]

    for z1 in euler(f1, g1, y0, z0, a, b, h):
        C[0].append(z1)
#    for z2 in euler(f2, g2, y0, z0, a, b, h):
#        C[1].append(z2)
    for z3 in euler(f3, g3, y0, z0, a, b, h):
        C[1].append(z3)

    return np.array(C)

def calc_zhat(C, n):
    theta0 = np.array([0.3, 10])
    D = 0.01
    M = 0
    error = np.random.normal(M, sqrt(D), size=(n,))

    z = C.dot(theta0) + error
    return z

def calc_theta_ls(zhat, C):
    theta = inv(C.T.dot(C)).dot(C.T).dot(zhat)
    return theta

def calc_ellipse(C):

    eig = np.linalg.eig(inv(C.T.dot(C)))
    a1 = eig[0][0] ** 0.5
    a2 = eig[0][1] ** 0.5
    phi1 = eig[1][0][:2]
    phi2 = eig[1][1][:2]

    alpha = abs(arctan(phi1[1]/phi1[0]))
    t = np.linspace(0., 2 * np.pi, 100)                         # изменение параметра t от 0 до 2pi
    ellipse_matrix = np.array([a1 * np.cos(t), a2 * np.sin(t)]) # параметрическое ур-е эллипса

    rotate_matrix = np.array([[cos(alpha), sin(alpha)],  # матрица поворота против часовой стрелке
                              [-sin(alpha), cos(alpha)]])

    inv_matrix = np.array([[-1, 0],                              # матрица инверсии по оси y x
                          [0, -1]])
    S = inv_matrix.dot(rotate_matrix)                           # итоговоя матрица перехода к базису
    ellipse_rotated = np.zeros((2, ellipse_matrix.shape[1]))

    for i in range(ellipse_rotated.shape[1]):                   # поворачиваем эллипс
        ellipse_rotated[:, i] = np.dot(S, ellipse_matrix[:, i]) 

    return ellipse_rotated

def main():
    a = 0
    b = 1000
    h = 1
    n = int((b - a) / h) + 1
    y0 = z0 = 0
    C = gen_c_matrix(f1, g1, f3, g3, y0, z0, a, b, h)
    C = C.T

    x = np.linspace(a, b, n)
    
    draw_func(x, y(x), 'y(t)')
    draw_func(x, derivative(x), 'y`(t)')
    draw_func(x, C[:, 0], 'u1(t)')
#    draw_func(x, C[:, 1], 'u2(t)')
    draw_func(x, C[:, 1], 'u3(t)')

    np.savetxt('y.txt', y(x))	
    np.savetxt('dy.txt', derivative(x))
    np.savetxt('u1.txt', C[:, 0])
#    np.savetxt('u2.txt', C[:, 1])
    np.savetxt('u3.txt', C[:, 1])

    z = calc_zhat(C, n)
    theta = calc_theta_ls(z, C)
    print (theta)
    ellipse = calc_ellipse(C)
    draw_func(ellipse[0, :], ellipse[1, :], 'mnk')


if __name__ == '__main__':
    main()
