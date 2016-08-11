import numpy


def basic_linear_regression(x, y):

    """Basic linear regression. Returns the co-efficients needed to plot a least squares fit line

    Parameters
    ----------
    x : [] list of floats or integers
    y : [] list of floats or integers

    Returns
    -------
    a : float
        gradient of least squares fit line (y = ax + b)
    b : float
        intercept of least squares fit line (y = ax + b)
    r_squared : float
        value for r-squared

    Raises
    ------
    ValueError if the lists are not the same length

    """

    try:
        len(x) == len(y)

    except:

        raise ValueError('Lists of x and y should be the same length\n')

    length = len(x)
    sum_x = sum(x)
    sum_y = sum(y)

    sum_x_squared = sum(map(lambda g: g*g, x))
    sum_of_products = sum(x[i] * y[i] for i in range(length))

    a = (sum_of_products - (sum_x * sum_y) / length) / (sum_x_squared - ((sum_x ** 2) / length))
    b = (sum_y - a * sum_x) / length

    ss_tot = 0
    ss_res = 0

    mean_y = numpy.mean(y)

    for i in y:
        ss_tot += (i - mean_y) ** 2

    for i, j in zip(x, y):

        f = (i * a) + b
        ss_res += (j - f) ** 2

    r_squared = 1 - (ss_res / ss_tot)

    return a, b, r_squared
