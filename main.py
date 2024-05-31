import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math


def cubic_spline(x, y):
    n = len(x)
    h = [x[i + 1] - x[i] for i in range(n - 1)]
    a = y
    b = [0] * (n - 1)
    d = [0] * (n - 1)
    c = solve_c(x, y, h, n)
    compute_coefficients(a, b, c, d, h, n)

    return a, b, c, d, x


def solve_c(x, y, h, n):
    A = [[0] * n for _ in range(n)]
    b = [0] * n

    # Boundary conditions (natural spline)
    A[0][0] = 1
    A[n - 1][n - 1] = 1

    # Fill the system with equations
    for i in range(1, n - 1):
        A[i][i - 1] = h[i - 1]
        A[i][i] = 2 * (h[i - 1] + h[i])
        A[i][i + 1] = h[i]
        b[i] = 3 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1])

    # Solve the system using Gaussian elimination
    for i in range(1, n):
        m = A[i][i - 1] / A[i - 1][i - 1]
        A[i][i] -= m * A[i - 1][i]
        b[i] -= m * b[i - 1]

    c = [0] * n
    c[n - 1] = b[n - 1] / A[n - 1][n - 1]
    for i in range(n - 2, -1, -1):
        c[i] = (b[i] - A[i][i + 1] * c[i + 1]) / A[i][i]

    return c


def compute_coefficients(a, b, c, d, h, n):
    for i in range(n - 1):
        b[i] = (a[i + 1] - a[i]) / h[i] - h[i] * (2 * c[i] + c[i + 1]) / 3
        d[i] = (c[i + 1] - c[i]) / (3 * h[i])


def evaluate_spline(a, b, c, d, x, x_val):
    i = find_segment(x, x_val)
    dx = x_val - x[i]
    return a[i] + b[i] * dx + c[i] * dx**2 + d[i] * dx**3


def find_segment(x, x_val):
    for i in range(len(x) - 1):
        if x[i] <= x_val <= x[i + 1]:
            return i


def linearPoints(x, y, points):
    n = len(x)
    xPoints = [0 for _ in range(points)]
    yPoints = [0 for _ in range(points)]
    for i in range(points):
        xPoints[i] = x[int((n - 1) / (points - 1) * i)]
        yPoints[i] = y[int((n - 1) / (points - 1) * i)]

    return xPoints, yPoints


def chebyshevNodes(x, y, points):
    n = len(x) - 1
    xPoints = [0 for _ in range(points)]
    yPoints = [0 for _ in range(points)]
    for k in range(points):
        index = int((math.cos((k * math.pi) / (points - 1)) + 1) / 2 * n)
        xPoints[points - k - 1] = x[index]
        yPoints[points - k - 1] = y[index]

    return xPoints, yPoints


def lagrangeInterpolation(xPoints, yPoints, x):
    n = len(xPoints)
    result = 0.0

    for i in range(n):
        term = yPoints[i]
        for j in range(n):
            if j != i:
                term *= (x - xPoints[j]) / (xPoints[i] - xPoints[j])
        result += term

    return result


def showWithSpecialFormat(xPoints, yPoints, x_new, y_new, interpolationType, nodePointsType):
    plt.scatter(xPoints, yPoints, color='red', label=f'Punkty węzłowe {nodePointsType}')
    plt.plot(x_new, y_new, label='Wielomian interpolacyjny')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'Interpolacja {interpolationType}')
    plt.show()


def interpolationForSpecificData(path, nodePoints):
    data = pd.read_csv(path)
    x = data['Dystans'].to_numpy()
    y = data['Wysokość'].to_numpy()

    xPoints, yPoints = linearPoints(x, y, nodePoints)

    x_new = [i for i in range(int(min(xPoints)), int(max(xPoints)), int((max(xPoints) - min(xPoints)) / 1000))]
    y_new = [lagrangeInterpolation(xPoints, yPoints, xi) for xi in x_new]
    showWithSpecialFormat(xPoints, yPoints, x_new, y_new, 'Lagrange', 'liniowe')

    a, b, c, d, xPoints = cubic_spline(xPoints, yPoints)
    x_new = [i for i in range(int(min(xPoints)), int(max(xPoints)), int((max(xPoints) - min(xPoints)) / 1000))]
    y_new = [evaluate_spline(a, b, c, d, xPoints, val) for val in x_new]
    showWithSpecialFormat(xPoints, yPoints, x_new, y_new, 'funkcji sklejania trzeciego stopnia', 'liniowe')

    xPoints, yPoints = chebyshevNodes(x, y, nodePoints)

    x_new = [i for i in range(int(min(xPoints)), int(max(xPoints)), int((max(xPoints) - min(xPoints)) / 1000))]
    y_new = [lagrangeInterpolation(xPoints, yPoints, xi) for xi in x_new]
    showWithSpecialFormat(xPoints, yPoints, x_new, y_new, 'Lagrange', 'Czebyszewa')

    a, b, c, d, xPoints = cubic_spline(xPoints, yPoints)
    x_new = [i for i in range(int(min(xPoints)), int(max(xPoints)), int((max(xPoints) - min(xPoints)) / 1000))]
    y_new = [evaluate_spline(a, b, c, d, xPoints, val) for val in x_new]
    showWithSpecialFormat(xPoints, yPoints, x_new, y_new, 'funkcji sklejania trzeciego stopnia', 'Czebyszewa')


if __name__ == '__main__':
    interpolationForSpecificData('2018_paths/WielkiKanionKolorado.csv', 30)
    interpolationForSpecificData('2018_paths/SpacerniakGdansk.csv', 20)
    interpolationForSpecificData('2018_paths/MountEverest.csv', 10)
