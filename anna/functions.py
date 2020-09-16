from numpy import cos, sin, sqrt, e

def y(t):
    
    return (10/427) * (e ** (-0.2 * t))*(-28*(e ** (-0.2 * t))*cos(2*t) + 28*cos(sqrt(14)*t/5) - 119 * (e ** (-0.2 * t))*sin(2*t) + 87*sqrt(14)*sin(sqrt(14)*t/5)) 

def x(t):
    return sin(2*t)

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
    Teta1 = 0.4
    Teta2 = 0.6
    Teta3 = 10
    return -Teta2 * z1 - Teta1 * u1 - derivative(t, y)

def f2(t, u2, z2):
    return z2

def g2(t, u2, z2):
    Teta1 = 0.4
    Teta2 = 0.6
    Teta3 = 10
    return -Teta2 * z2 - Teta1 * u2 - y(t)

def f3(t, u3, z3):
    return z3

def g3(t, u3, z3):
    Teta1 = 0.4
    Teta2 = 0.6
    Teta3 = 10
    return -Teta2 * z3 - Teta1 * u3 + x(t)
