from numpy import cos, sin, sqrt, e

def y(t):
    c1 = -17.801
    c2 = -50
    beta = 0.421
    alfa = -0.15
    yo = 50
    a = c1 * sin(beta * t)
    b = c2 * cos(beta * t)
    return (a + b) * (e ** (alfa * t)) + yo 

def derivative(t):
  # eps = 0.0001
  # return (f(x + eps) - f(x)) / eps
    c1 = -17.801
    c2 = -50
    c3 = 21.05
    c4 = -7.494221
    beta = 0.421
    alfa = -0.15
    a = c1 * sin(beta * t)
    b = c2 * cos(beta * t)
    c = c3 * sin(beta * t)
    d = c4 * cos(beta * t)
    return (c + d) * (e ** (alfa * t)) + alfa*(a + b) * (e ** (alfa * t)) 

def f1(t, u1, z1):
    return z1

def g1(t, u1, z1):
    return -0.3 * z1 - 0.2 * u1 - derivative(t)

def f2(t, u2, z2):
    return z2

def g2(t, u2, z2):
    return -0.3 * z2 - 0.2 * u2 - y(t)

def f3(t, u3, z3):
    return z3

def g3(t, u3, z3):
    return -0.3 * z3 - 0.2 * u3 + 1
