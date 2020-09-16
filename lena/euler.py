import functions
import numpy as np


def euler(f, g, y0, z0, a, b, h):
    n = int((b - a) / h) + 1
    x_prev = xi = a
    y_prev = yi = y0
    z_prev = zi = z0

    for i in range(n):
        yi = y_prev + h * f(x_prev, y_prev, z_prev)
        zi = z_prev + h * g(x_prev, y_prev, z_prev)
        xi = x_prev + h
        yield z_prev

        x_prev = xi
        y_prev = yi
        z_prev = zi
