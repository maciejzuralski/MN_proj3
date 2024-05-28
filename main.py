import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math


def linearPoints(x, y, points):
    n = len(x)
    points = 10
    x_points = np.zeros([points])
    y_points = np.zeros([points])
    for i in range(points):
        x_points[i] = x[int((n - 1) / (points - 1) * i)]
        y_points[i] = y[int((n - 1) / (points - 1) * i)]

    return x_points, y_points


def chebyshevNodes(x, y, points):
    n = len(x) - 1
    x_points = np.zeros([points])
    y_points = np.zeros([points])
    for k in range(points):
        index = int((math.cos((k * math.pi) / (points - 1)) + 1) / 2 * n)
        x_points[k] = x[index]
        y_points[k] = y[index]

    return x_points, y_points


def lagrangeInterpolation(x_points, y_points, x):
    """
    Funkcja oblicza wartość wielomianu interpolacyjnego Lagrange'a w punkcie x.

    :param x_points: Lista współrzędnych x punktów węzłowych.
    :param y_points: Lista współrzędnych y punktów węzłowych.
    :param x: Punkt, w którym obliczamy wartość wielomianu interpolacyjnego.
    :return: Wartość wielomianu interpolacyjnego w punkcie x.
    """
    n = len(x_points)
    result = 0.0

    for i in range(n):
        term = y_points[i]
        for j in range(n):
            if j != i:
                term *= (x - x_points[j]) / (x_points[i] - x_points[j])
        result += term

    return result


def mountEverest():
    data = pd.read_csv('2018_paths/MountEverest.csv')
    x = data['Dystans (m)'].to_numpy()
    y = data['Wysokość (m)'].to_numpy()

    #x_points, y_points = linearPoints(x, y, 10)
    x_points, y_points = chebyshevNodes(x, y, 30)

    x_new = np.linspace(min(x_points), max(x_points), 1000)
    y_new = [lagrangeInterpolation(x_points, y_points, xi) for xi in x_new]

    plt.scatter(x_points, y_points, color='red', label='Punkty węzłowe')
    plt.plot(x_new, y_new, label='Wielomian interpolacyjny Lagrange\'a')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Interpolacja Lagrange\'a')
    plt.show()


if __name__ == '__main__':
    mountEverest()

    # Przykładowe punkty
    x_points = [1, 2, 3, 4]
    y_points = [1, 4, 9, 5]

    # Tworzenie nowych punktów do rysowania wielomianu interpolacyjnego
    x_new = np.linspace(min(x_points), max(x_points), 100)
    y_new = [lagrangeInterpolation(x_points, y_points, xi) for xi in x_new]

    # Rysowanie wyników
    plt.scatter(x_points, y_points, color='red', label='Punkty węzłowe')
    plt.plot(x_new, y_new, label='Wielomian interpolacyjny Lagrange\'a')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Interpolacja Lagrange\'a')
    #plt.show()