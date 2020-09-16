from numpy import cos, sin, sqrt, e

def y(t):
    c1 = -75
    c2 = -225/sqrt(71)
    beta = sqrt(71)/20
    alfa = -3/20
    yo = 75
    a = c1 * sin(beta * t)
    b = c2 * cos(beta * t)
    return (a + b) * (e ** (alfa * t)) + yo 

def x(t):
    return 1

def derivative(x, f):
   eps = 0.0001
   return (f(x + eps) - f(x)) / eps
  # c1 = -17.801
  #  c2 = -50
  #  c3 = 21.05
  #  c4 = -7.494221
  #  beta = 0.421
  #  alfa = -0.15
  #  a = c1 * sin(beta * t)
  #  b = c2 * cos(beta * t)
  #  c = c3 * sin(beta * t)
  #  d = c4 * cos(beta * t)
  #  return (c + d) * (e ** (alfa * t)) + alfa*(a + b) * (e ** (alfa * t)) 

def f1(t, u1, z1):
    return z1

def g1(t, u1, z1):
    Teta1 = 0.3
    Teta2 = 0.2
    Teta3 = 15
    return -Teta1 * z1 - Teta2 * u1 - derivative(t, y)

def f2(t, u2, z2):
    return z2

def g2(t, u2, z2):
    Teta1 = 0.3
    Teta2 = 0.2
    Teta3 = 15
    return -Teta1 * z2 - Teta2 * u2 - y(t)

def f3(t, u3, z3):
    return z3

def g3(t, u3, z3):
    Teta1 = 0.3
    Teta2 = 0.2
    Teta3 = 15
    return -Teta1 * z3 - Teta2 * u3 + x(t)
