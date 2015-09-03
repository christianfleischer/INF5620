import sympy as sym
import numpy as np
import matplotlib.pyplot as plt

V, t, I, w, dt, a, b, c, d = sym.symbols('V t I w dt a b c d')  # global symbols
f = None  # global variable for the source term in the ODE

def ode_source_term(u):
    """Return the terms in the ODE that the source term
    must balance, here u'' + w**2*u.
    u is symbolic Python function of t."""
    return sym.diff(u(t), t, t) + w**2*u(t)

def residual_discrete_eq(u):
    """Return the residual of the discrete eq. with u inserted."""
    R = DtDt(u, dt) + w**2*u(t) - f
    return sym.simplify(R)

def residual_discrete_eq_step1(u):
    """Return the residual of the discrete eq. at the first
    step with u inserted."""
    R = (DtDt(u, dt) + w**2*u(t) - f).subs(t, 0)
    return sym.simplify(R)

def DtDt(u, dt):
    """Return 2nd-order finite difference for u_tt.
    u is a symbolic Python function of t.
    """
    return (u(t+dt)-2*u(t)+u(t-dt))/(dt**2)

def main(u):
    """
    Given some chosen solution u (as a function of t, implemented
    as a Python function), use the method of manufactured solutions
    to compute the source term f, and check if u also solves
    the discrete equations.
    """
    print '=== Testing exact solution: %s ===' % u(t)
    print "Initial conditions u(0)=%s, u'(0)=%s:" % \
          (u(t).subs(t, 0), sym.diff(u(t), t).subs(t, 0))

    # Method of manufactured solution requires fitting f
    global f  # source term in the ODE
    f = sym.simplify(ode_source_term(u))
    #f = sym.simplify(ode_lhs(u))

    # Residual in discrete equations (should be 0)
    print 'residual step1:', residual_discrete_eq_step1(u)
    print 'residual:', residual_discrete_eq(u)

def linear():
    main(lambda t: V*t + I)

def quadratic():
    main(lambda t: b*t**2 + c*t + d)

def cubic():
    main(lambda t: a*t**3 + b*t**2 + c*t + d)

def f_d(u_d, b_, w_, n):
    return 2*b_ + w_**2*u_d[n]

def solver(I, w_, dt, T, V_, b_):
    """
    Solve u'' + w**2*u = f for t in (0,T], u(0)=I and u'(0)=0,
    by a central finite difference method with time step dt.
    """
    dt = float(dt)
    Nt = int(round(T/dt))
    u_d = np.zeros(Nt+1)
    t = np.linspace(0, Nt*dt, Nt+1)
    #f_d = np.linspace(f(0), f(T), Nt+1)

    u_d[0] = I
    u_d[1] = u_d[0] - 0.5*dt**2*w_**2*u_d[0]+ dt*V_ + 0.5*dt**2*f_d(u_d, b_, w_, 0)
    #u_d[1] = u_d[0] - 0.5*dt**2*w_**2*u_d[0]+ dt*V_ + 0.5*dt**2*(2*b_ + w_**2*u_d[0])
    for n in range(1, Nt):
        u_d[n+1] = 2*u_d[n] - u_d[n-1] - dt**2*w_**2*u_d[n] + dt**2*f_d(u_d, b_, w_, n)
        #u_d[n+1] = 2*u_d[n] - u_d[n-1] - dt**2*w_**2*u_d[n] + dt**2*(2*b_ + w_**2*u_d[n])
    return u_d, t

def u_e(t, b_, c_, d_):
    return b_*t**2 + c_*t + d_

def nose_test(d_, w_, dt, T, c_, b_):
    Nt=T/dt
    t_e = np.linspace(0, T, Nt+1)

    [u_d, t] = solver(d_, w_, dt, T, c_, b_)
    plt.plot(t, u_d,
    t_e, u_e(t_e, b_, c_, d_))
    plt.show()

    tol = 1E-12
    error = abs(u_e(t_e, b_, c_, d_) - u_d)#/u_e(t, 1, 1, 1)
    assert error.max() < tol
    plt.plot(t_e, error)
    plt.show()

if __name__ == '__main__':
    linear()
    quadratic()
    cubic()
    nose_test(1, 1, 0.1, 10, 1, 1)

