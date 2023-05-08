NEWTOWN_SIZE_LIMIT = 2000

def g_lam(lam, mu, Ainv, term3, term4, ones):
    '''Gradient of dg(lam, mu)/d_lam'''

    d_lam = 0.5 * Ainv.T @ lam
    d_lam -= Ainv @ lam
    d_lam -= term3
    d_lam -= term4
    d_lam += (0.25 * mu * Ainv @ ones).squeeze()
    d_lam += (0.5 * mu * Ainv.T @ ones).squeeze()

    return d_lam

def g_mu_non_const(lam, mu, Ainv, ones):
    '''Gradient of dg(lam, mu)/d_lam only variable dependant terms'''
    d = -0.25 * lam.T @ Ainv.T @ ones
    d -= 0.25 * ones.T @ Ainv @ lam
    d += 0.5 * (mu * ones.T @ Ainv.T @ ones).squeeze()
    d += 0.5 * (lam.T @ Ainv @ ones).squeeze()
    d += 0.5 * (ones.T @ Ainv @ lam).squeeze()
    return d


def d_mu_const(H, Ainv, y, ones):
    '''Gradient of dg(lam, mu)/d_lam only variable dependant terms'''
    dc = -0.5 * y.T @ H @ Ainv @ ones
    dc -= 0.5 * ones.T @ Ainv.T @ H.T @ y
    dc += y.T @ H @ Ainv @ ones
    dc += (ones.T @ Ainv @ H.T @ y).squeeze() - 0.5 * (ones.T @ Ainv @ ones).squeeze()
    dc -= 1
    return dc


def lg_BarrierMethod(H, y, Ainv, mu_init, lam_init, m, t, nu=0.01,
                     tol_barrier=1e-5,
                     tol_newton=1e-5, max_iter=100, eta=0.0001, n_netwon_steps=14):
    '''Main algorithm for Barrier Method'''

    lam = lam_init  # store initial value
    mu = mu_init
    duality_gap = [m / t]  # initialize tabulation of duality gap
    k = 0  # number of iterations

    ones = np.ones(shape=(H.shape[1], 1))
    g_term3 = 0.5 * Ainv @ H.T @ y
    g_term4 = 0.5 * Ainv.T @ H.T @ y
    Hess_lam = 0.5 * Ainv.T - Ainv
    Hess_mu = 0.5 * ones.T @ Ainv @ ones
    gmu_const = d_mu_const(H, Ainv, y, ones)
    x = np.hstack((lam, [mu]))

    Hlam_inv = np.linalg.solve(Hess_lam, np.eye(Hess_lam.shape[1]))
    Hmu_inv = np.linalg.solve(Hess_mu, np.eye(Hess_mu.shape[1]))

    # loop until stopping criterion is met
    while m / t > tol_barrier:

        # centering step: Newton Algorithm
        i = 0
        d = np.hstack((np.ones_like(lam), [1]))

        while np.linalg.norm(d) > tol_newton and i < max_iter:
            lam, mu = x[:-1], x[-1]
            i += 1

            # calculate step
            glam = g_lam(lam, mu, Ainv, g_term3, g_term4, ones)
            gmu_non_const = g_mu_non_const(lam, mu, Ainv, ones)
            gmu = gmu_const + gmu_non_const
            if i <= n_netwon_steps:  # Newton step
                d_lam = np.dot(Hlam_inv, glam)
                d_mu = np.dot(Hmu_inv, gmu)
            else:  # GD
                d_lam = eta * glam
                d_mu = eta * gmu
            d = np.hstack((d_lam, d_mu))

            # update x
            valid_inds = np.where(x + d > 0)[0]
            x[valid_inds] += d[valid_inds]
            x /= np.linalg.norm(x)

        # update parameter t
        t = (1 + 1 / (13 * np.sqrt(nu))) * t

        # update tabulations
        duality_gap.append(m / t)
        k += 1

    return x


def lagrangian(H, y):
    n = H.shape[1]
    x_init = np.random.random(n + 1)
    m = n
    t = 1
    Ainv = np.linalg.solve(H.T @ H, np.eye(H.shape[1]))
    lam_init, mu_init = x_init[:-1], x_init[-1]
    x = lg_BarrierMethod(H, y, Ainv, mu_init, lam_init, m, t, tol_newton=1e-6, max_iter=10000, eta=1e-8)
    lam, mu = x[:-1], x[-1]
    x_approx = Ainv @ (H.T @ y + 0.5 * lam - 0.5 * mu * np.ones_like(lam))
    x_approx[x_approx < 0] = 0
    x_approx /= np.sum(x_approx)
    return x_approx


def g(H, x, t, eps, g_term2):
    '''Gradient of f(x) + 1/t * ϕ(x)'''

    x = 2.0 * x.T @ H.T @ H
    x += g_term2
    x += - (1 / t) * (1 / (x + eps))

    return x


def Hess(x, t, H_term1, eps=1e-15):
    '''Hessian of f(x) + 1/t * ϕ(x)'''

    x = x.squeeze()

    # plug x into the following formula:
    # 2 * (H.T @ H) + diag(1 / x^-2)
    x = H_term1
    x += np.diag((1 / t) * (1 / (x ** 2 + eps)))
    return x


def gd_BarrierMethod(H, y, x_init, m, t, g_term2, H_term1, nu=0.01,
                     epsilon=1e-15, tol_barrier=1e-5,
                     tol_newton=1e-5, max_iter=100, eta=0.0001, n_netwon_steps=14):
    '''Main algorithm for Barrier Method'''
    x = x_init  # store initial value
    duality_gap = [m / t]  # initialize tabulation of duality gap
    k = 0  # number of iterations

    # loop until stopping criterion is met
    while m / t > tol_barrier:

        # centering step: Newton Algorithm
        i = 0
        d = np.ones_like(x)
        while np.linalg.norm(d) > tol_newton and i < max_iter:
            i += 1

            # calculate step
            gx = g(H, x, t, eps=epsilon, g_term2=g_term2)
            d = -(eta * gx)

            # update x
            valid_inds = np.where(x + d > 0)[0]
            x[valid_inds] += d[valid_inds]
            x /= np.sum(x)

        # update parameter t
        t = (1 + 1 / (13 * np.sqrt(nu))) * t

        # update tabulations
        duality_gap.append(m / t)
        k += 1

    return x


def simple_gd(H, y):
    n = H.shape[1]
    x_init = np.random.random(n)
    x_init = x_init / np.sum(x_init)
    g_term2 = - 2 * (y.T @ H)
    H_term1 = 2 * (H.T @ H)
    m = n
    t = 0.1
    x = gd_BarrierMethod(H, y, x_init, m, t,
                         max_iter=10,
                         eta=1e-8, g_term2=g_term2, H_term1=H_term1)
    return x


def solve(H, y):
    if np.max(H.shape) > NEWTOWN_SIZE_LIMIT:
        return simple_gd(H, y)
    else:
        return lagrangian(H, y)


if __name__ == '__main__':
    n = len(examples)
    for i in range(6, 8):
        H, y = examples[i].H, examples[i].y
        print(f'Example {i} (H shape: {H.shape})')
        x_approx = solve(H, y)
        print('solver score:', np.linalg.norm(examples[i].y - examples[i].H @ x_approx))
        print('validation:', np.isclose(np.sum(x_approx), 1.), 'sum=', np.sum(x_approx))
        print('validation xi >= 0:', np.all(x_approx >= 0))
        print("\n" + "#" * 20 + "\n")
