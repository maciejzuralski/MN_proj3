import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math


def cubicSpline(x, y):
    n = len(x)
    h = [x[i + 1] - x[i] for i in range(n - 1)]
    a = y
    b = [0 for _ in range(n - 1)]
    d = [0 for _ in range(n - 1)]
    c = solveC(x, y, h, n)
    computeCoefficients(a, b, c, d, h, n)

    return a, b, c, d, x


def solveC(x, y, h, n):
    A = [[0 for _ in range(n)] for _ in range(n)]
    b = [0 for _ in range(n)]

    A[0][0] = 1
    A[n - 1][n - 1] = 1

    for i in range(1, n - 1):
        A[i][i - 1] = h[i - 1]
        A[i][i] = 2 * (h[i - 1] + h[i])
        A[i][i + 1] = h[i]
        b[i] = 3 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1])

    c = np.linalg.solve(np.array(A), np.array(b))
    return c


def computeCoefficients(a, b, c, d, h, n):
    for i in range(n - 1):
        b[i] = (a[i + 1] - a[i]) / h[i] - h[i] * (2 * c[i] + c[i + 1]) / 3
        d[i] = (c[i + 1] - c[i]) / (3 * h[i])


def evaluateSpline(a, b, c, d, x, xVal):
    i = findSegment(x, xVal)
    dx = xVal - x[i]
    return a[i] + b[i] * dx + c[i] * dx**2 + d[i] * dx**3


def findSegment(x, xVal):
    for i in range(len(x) - 1):
        if x[i] <= xVal <= x[i + 1]:
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


def showWithSpecialFormat(x, y, xPoints, yPoints, x_new, y_new, interpolationType, nodePointsType, path):
    plt.scatter(xPoints, yPoints, color='red', label=f'Punkty węzłowe {nodePointsType}')
    plt.plot(x_new, y_new, label='Wielomian interpolacyjny')
    plt.plot(x, y, label='oryginalne wartości')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'{path}\nInterpolacja {interpolationType}')
    plt.show()


def interpolationForSpecificData(path, nodePoints):
    data = pd.read_csv(path)
    x = data['Dystans'].to_numpy()
    y = data['Wysokość'].to_numpy()

    xPoints, yPoints = linearPoints(x, y, nodePoints)

    xNew = [i for i in range(int(min(xPoints)), int(max(xPoints)), int((max(xPoints) - min(xPoints)) / 1000))]
    yNew = [lagrangeInterpolation(xPoints, yPoints, xi) for xi in xNew]
    showWithSpecialFormat(x, y, xPoints, yPoints, xNew, yNew, 'Lagrange', 'liniowe', path)

    a, b, c, d, xPoints = cubicSpline(xPoints, yPoints)
    xNew = [i for i in range(int(min(xPoints)), int(max(xPoints)), int((max(xPoints) - min(xPoints)) / 1000))]
    yNew = [evaluateSpline(a, b, c, d, xPoints, val) for val in xNew]
    showWithSpecialFormat(x, y, xPoints, yPoints, xNew, yNew, 'funkcji sklejania trzeciego stopnia', 'liniowe', path)

    xPoints, yPoints = chebyshevNodes(x, y, nodePoints)

    xNew = [i for i in range(int(min(xPoints)), int(max(xPoints)), int((max(xPoints) - min(xPoints)) / 1000))]
    yNew = [lagrangeInterpolation(xPoints, yPoints, xi) for xi in xNew]
    showWithSpecialFormat(x, y, xPoints, yPoints, xNew, yNew, 'Lagrange', 'Czebyszewa', path)

    a, b, c, d, xPoints = cubicSpline(xPoints, yPoints)
    xNew = [i for i in range(int(min(xPoints)), int(max(xPoints)), int((max(xPoints) - min(xPoints)) / 1000))]
    yNew = [evaluateSpline(a, b, c, d, xPoints, val) for val in xNew]
    showWithSpecialFormat(x, y, xPoints, yPoints, xNew, yNew, 'funkcji sklejania trzeciego stopnia', 'Czebyszewa', path)


if __name__ == '__main__':
    interpolationForSpecificData('paths/MountEverest.csv', 5)
    interpolationForSpecificData('paths/WielkiKanionKolorado.csv', 20)
    interpolationForSpecificData('paths/SpacerniakGdansk.csv', 20)
