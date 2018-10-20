import numpy as np

from prox_tv import tv1_1d, tv1w_1d, tv2_1d, tv1_2d, tvp_1d, \
                    tv1w_2d, tvp_2d, tvgen


def test_tv1w_1d():
    methods = ('tautstring', 'pn')
    for _ in range(20):
        dimension = np.random.randint(1e1, 1e3)
        x = 100*np.random.randn(dimension)
        w = 20*np.random.rand(dimension-1)
        solutions = [tv1w_1d(x, w, method=method) for method in methods]
        for i in range(len(solutions)-1):
            assert np.allclose(solutions[i], solutions[i+1])


def test_tv1w_1d_uniform_weights():
    _test_tv1w_1d_uniform_weights(1e1, 1e4)


def test_tv1w_1d_uniform_weights_small_input():
    _test_tv1w_1d_uniform_weights(2, 4)


def _test_tv1w_1d_uniform_weights(min, max):
    for _ in range(1000):
        dimension = np.random.randint(min, max)
        x = 100*np.random.randn(dimension)
        w1 = np.random.rand()
        w = np.ones(dimension-1) * w1 
        solw = tv1w_1d(x, w)
        sol = tv1_1d(x, w1)
        assert np.allclose(solw, sol)


def test_tv1_1d():
    """Tests 1-dimensional TV methods"""
    for _ in range(20):
        dimension = np.random.randint(1e1, 3e1)
        x = 100*np.random.randn(dimension)
        w = 20*np.random.rand()
        _test_tv1_methods(x, w)


def test_tv1_1d_int():
    """Tests 1-dimensional TV methods for integer inputs"""
    for _ in range(20):
        dimension = np.random.randint(1e1, 3e1)
        x = (100*np.random.randn(dimension)).astype('int')
        w = 20*np.random.rand()
        _test_tv1_methods(x, w)


def _test_tv1_methods(x, w):
    """For given input signal and weight, all TV1 methods must be similar"""
    methods = ('classictautstring', 'linearizedtautstring', 'hybridtautstring',
               'pn', 'condat', 'dp', 'condattautstring', 'kolmogorov')
    solutions = [tv1_1d(x, w, method=method) for method in methods]
    solutions.append(tv1_1d(x, w, method='hybridtautstring', maxbacktracks=1.2))
    for i in range(1, len(solutions)):
        assert np.allclose(solutions[0], solutions[i], atol=1e-3)


def test_tvp_1d():
    """Test that all 1D-lp-TV methods produce equivalent results"""
    # Some of these methods are kind of unstable, so we ensure only that most
    # of the times the results are similar
    methods = ('gp', 'fw', 'gpfw')
    errors = 0
    for _ in range(20):
        dimension = np.random.randint(1e1, 3e1)
        x = 100*np.random.randn(dimension)
        w = 20*np.random.rand()
        p = 1 + 10 * np.random.rand()
        solutions = [tvp_1d(x, w, p, method=method, max_iters=100000)
                     for method in methods]
        for i in range(1, len(solutions)):
            try:
                assert np.allclose(solutions[0], solutions[i], atol=1e-3)
            except AssertionError:
                errors += 1

    assert(errors < 10)


def test_tv2_1d():
    # TODO: issues in osx with this test
    methods = ('ms', 'pg', 'mspg')
    for _ in range(20):
        dimension = np.random.randint(1e1, 3e1)
        x = 100*np.random.randn(dimension)
        w = 20*np.random.rand()
        solutions = [tv2_1d(x, w, method=method) for method in methods]
        for i in range(len(solutions)-1):
            assert np.allclose(solutions[i], solutions[i+1], atol=1e-3)


def _generate2d():
    """Generates a 2D array for the test"""
    rows = np.random.randint(1e1, 3e1)
    cols = np.random.randint(1e1, 3e1)
    return 100*np.random.randn(rows, cols) 


def test_tv1_2d():
    """Tests that all 2D-TV methods produce equivalent results"""
    methods = ('yang', 'condat', 'chambolle-pock', 'kolmogorov', 'pd', 'dr')
    for _ in range(20):
        x = _generate2d()
        w = 20*np.random.rand()
        solutions = [tv1_2d(x, w, method=method, max_iters=5000)
                     for method in methods]
        for i in range(1, len(solutions)):
            print(methods[i], ":", solutions[i])
            assert np.allclose(solutions[i], solutions[0], atol=1e-3)
            

def test_tv1_tvp_2d():
    """Tests that 2D-TVp == 2D-TV1 when p=1"""
    for _ in range(20):
        x = _generate2d()
        w = 20*np.random.rand()
        solution1 = tv1_2d(x, w, max_iters=5000)
        solutionp = tvp_2d(x, w, w, 1, 1, max_iters=5000)
        assert np.allclose(solution1, solutionp, atol=1e-3)
        

def test_tv1_tv1w_2d():
    """Tests that 2D-TV1w == 2D-TV1 for unit weights"""
    for _ in range(20):
        x = _generate2d()
        rows = len(x)
        cols = len(x[0])
        w = 20*np.random.rand()
        w_cols = w * np.ones((rows-1, cols))
        w_rows = w * np.ones((rows, cols-1))
        solution1 = tv1_2d(x, w, max_iters=5000)
        solutionp = tv1w_2d(x, w_cols, w_rows, max_iters=5000)
        assert np.allclose(solution1, solutionp, atol=1e-3)


def test_tv1w_2d_uniform_weights():
    for _ in range(20):
        x = _generate2d()
        rows = len(x)
        cols = len(x[0])
        w1 = np.random.rand()
        w_rows = np.ones([rows-1, cols]) * w1 
        w_cols = np.ones([rows, cols-1]) * w1 
        solw = tv1w_2d(x, w_rows, w_cols, max_iters=5000)
        solw1 = tv1_2d(x, w1, max_iters=5000)
        assert np.allclose(solw, solw1, atol=1e-3)


def test_tv1w_2d_uniform_weights_small_input():
    for _ in range(1000):
        rows = np.random.randint(2, 4)
        cols = np.random.randint(2, 4)
        x = 100*np.random.randn(rows, cols)
        w1 = np.random.rand()
        w_rows = np.ones([rows-1, cols]) * w1 
        w_cols = np.ones([rows, cols-1]) * w1 
        solw = tv1w_2d(x, w_rows, w_cols, max_iters=5000)
        solw1 = tv1_2d(x, w1, max_iters=5000)
        assert np.allclose(solw, solw1, atol=1e-3)


def test_tv1w_2d_emengd():
    r"""Issue reported by emengd
    
    Make the solver fail due to missing checks on integer arguments
    """        
    a = -np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])/10.
    sol1 = tv1w_2d(a, np.array([[1, 1, 1], [1, 1, 1]]),
                   np.array([[1, 1], [1, 1], [1, 1]]), max_iters=100)
    sol2 = tv1_2d(a, 1)
    assert np.allclose(sol1, sol2, atol=1e-3)


def test_tvgen_1d():
    """Tests that the general solver returns correct 1d solutions"""
    for _ in range(20):
        dimension = np.random.randint(1e1, 3e1)
        x = 100*np.random.randn(dimension)
        w = 20*np.random.rand()
        specific = tv1_1d(x, w)
        general = tvgen(x, [w], [1], [1])
        assert np.allclose(specific, general, atol=1e-3)


def test_tvgen_2d():
    """Tests that the general solver returns correct 2d solutions"""
    for _ in range(20):
        x = _generate2d()
        w = 20*np.random.rand()
        specific = tv1_2d(x, w, max_iters=1000)
        general = tvgen(x, [w, w], [1, 2], [1, 1], max_iters=1000)
        assert np.allclose(specific, general, atol=1e-2)


def test_tvgen_nd():
    """Test that the general solver does not crash for high-d tensors"""
    for _ in range(20):
        dims = np.random.randint(3, 5)
        shape = np.random.randint(2, 10, size=dims)
        x = np.random.randn(*shape)
        w = np.random.randn(dims)
        tvgen(x, w, list(range(1, dims+1)), np.ones(dims))


def test_tvgen_multireg():
    """Test applying several regularizers on same dimension"""
    for _ in range(20):
        x = _generate2d()
        w = 20*np.random.rand()
        specific = tv1_2d(x, w, max_iters=1000)
        general = tvgen(
                x,
                [w/2., w/2., w/3., w/3., w/3.],
                [1, 1, 2, 2, 2],
                [1, 1, 1, 1, 1],
                max_iters=1000
        )
        print("Max diff: " + str((specific-general).max()))
        assert np.allclose(specific, general, atol=1e-2)
