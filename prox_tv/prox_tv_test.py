import numpy as np

from prox_tv import tv1_1d, tv1w_1d, tv2_1d, tv1_2d, tvp_1d, \
                    tv1w_2d, tvp_2d


def test_tv1w_1d():
    methods = ('tautstring', 'pn')
    for _ in range(20):
        dimension = np.random.randint(1e1, 1e3)
        x = 100*np.random.randn(dimension)
        w = 20*np.random.rand(dimension-1)
        solutions = [tv1w_1d(x, w, method=method) for method in methods]
        for i in range(len(solutions)-1):
            assert np.allclose(solutions[i], solutions[i+1])


def test_tv1_1d():
    methods = ('tautstring', 'pn', 'condat', 'dp')
    for _ in range(20):
        dimension = np.random.randint(1e1, 3e1)
        x = 100*np.random.randn(dimension)
        w = 20*np.random.rand()
        solutions = [tv1_1d(x, w, method=method) for method in methods]
        for i in range(1, len(solutions)):
            assert np.allclose(solutions[0], solutions[i], atol=1e-3)


def test_tvp_1d():
    methods = ('gp', 'fw', 'gpfw')
    for _ in range(20):
        dimension = np.random.randint(1e1, 3e1)
        x = 100*np.random.randn(dimension)
        w = 20*np.random.rand()
        p = 1 + 10 * np.random.rand()
        solutions = [tvp_1d(x, w, p, method=method, max_iters=100000)
                     for method in methods]
        for i in range(1, len(solutions)):
            assert np.allclose(solutions[0], solutions[i], atol=1e-3)


def test_tv2_1d():
    methods = ('ms', 'pg', 'mspg')
    for _ in range(20):
        dimension = np.random.randint(1e1, 3e1)
        x = 100*np.random.randn(dimension)
        w = 20*np.random.rand()
        solutions = [tv2_1d(x, w, method=method) for method in methods]
        for i in range(len(solutions)-1):
            assert np.allclose(solutions[i], solutions[i+1], atol=1e-3)


def test_tv1_2d():
    methods = ('yang', 'condat', 'chambolle-pock')
    for _ in range(20):
        rows = np.random.randint(1e1, 3e1)
        cols = np.random.randint(1e1, 3e1)
        x = 100*np.random.randn(rows, cols)
        w = 20*np.random.rand()
        solutions = [tv1_2d(x, w, method=method, max_iters=5000)
                     for method in methods]
        solutions.append([tvp_2d(x, w, w, 1, 1, max_iters=5000)])
        w_cols = w * np.ones((rows-1, cols))
        w_rows = w * np.ones((rows, cols-1))
        solutions.append(tv1w_2d(x, w_cols, w_rows, max_iters=5000))
        for i in range(1, len(solutions)):
            assert np.allclose(solutions[i], solutions[0], atol=1e-3)
