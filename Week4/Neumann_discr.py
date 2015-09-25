from wave1D_dn_vc_edited import *
import matplotlib.pylab as mpl

def find_f(L_value, q):
        from sympy import symbols, cos, pi, diff, simplify, lambdify
        xs, ts, Ls = symbols('x t L')
        
        u_e = lambda x, t: cos(np.pi*x/Ls)*cos(t)
        #q   = lambda x: 1 + cos(np.pi*xs/Ls)#1. + (x-Ls/2)**4
        
        utt = diff(u_e(xs,ts),ts,ts)
        ux = diff(u_e(xs,ts),xs)   
        qux_x = diff(q(xs)*ux,xs)
        
        f_ = simplify(utt - qux_x).subs(Ls,L_value)
        return lambdify((xs,ts), f_, modules='numpy')

def ex(exercise):

    L = 1.
    u_exact = lambda x, t: np.cos(np.pi*x/L)*np.cos(t)
    
    if exercise == "a":
        q = lambda x: 1. + (x-L/2)**4
    elif exercise == "b":
        q = lambda x: 1. + np.cos(np.pi*x/L)
    else:
        q = lambda x: 1. + (x-L/2)**4
        
    c = lambda x: np.sqrt(q(x))#np.sqrt(1. + (x-L/2)**4)
    I = lambda x: u_exact(x, 0)
    V = None#lambda x: 0.5*u_exact(x, 0) 
    U_0 = None#lambda t: u_exact(0, t)
    U_L = None#lambda t: u_exact(L, t)#None
    C = 1.
    #Nx = 100.
    #dt = C*(L/Nx)/c(0)
    T = 3.    
    f = find_f(L, q)

    def plot_ex(Nx, approx):
        dt = C*(L/Nx)/c(0)
        solver(I, V, f, c, U_0, U_L, L, dt, C, T,
               user_action=PlotAndStoreSolution('moving_end'),
               version='vectorized',
               stability_safety_factor=1, approx=approx)
    
    def conv_rate(Nx_values, approx):
    
        def assert_no_error(u, x, t, n):
            u_e = u_exact(x, t[n])
            diff = np.abs(u - u_e).max()
            e = u_e - u
            dt = t[1]-t[0]
            E = np.sqrt(dt*sum(e**2))
            tol = 1E-13
            Errors_t.append(E)
            #print E
            #print diff
            #assert diff < tol
        
        Error_t_max = []
        dt_values = []
           
        for i in range(len(Nx_values)):
            dt = C*(L/Nx_values[i])/c(0)
            dt_values.append(dt)
        
            Errors_t = []
            t, cpu_time = solver(I, V, f, c, U_0, U_L, L, dt, C, T,
                                user_action=assert_no_error,
                                version='vectorized',
                                stability_safety_factor=1,
                                approx=approx)
            Error_t_max.append(max(Errors_t))

        r = compute_rates(dt_values, Error_t_max)
        for i in range(len(Nx_values)-1):
            print 'For Nx[i-1]=%s and Nx[i]=%s, r=%s (approximation used: %s)' % (Nx_values[i],Nx_values[i+1], r[i], approx)
        print ''
        #print dt_values
        #print len(t)
        mpl.plot(t, Errors_t)
        mpl.show()
        #print dt_values[-1]**r[-1]
           
    def compute_rates(dt_values, E_values):
        m = len(dt_values)
        r = [np.log(E_values[i-1]/E_values[i])/
        np.log(dt_values[i-1]/dt_values[i])
        for i in range(1, m, 1)]
        #Round to two decimals
        r = [round(r_, 5) for r_ in r]
        return r
    
    Nx_values = np.linspace(10, 100, 10)
    
    if exercise == "a" or exercise == "b":
        for i in ["Standard", "Approx"]:#, "OneSided"]:
            conv_rate(Nx_values, approx=i)
        
        rinp = raw_input("Plot solution for a given Nx (y/n)? ")
        if rinp == "y":
            Nx = raw_input("Enter Nx value: ")
            approx = raw_input("Enter desired scheme (Standard/Approx): ")
            plot_ex(int(Nx), approx)
            
    elif exercise == "c":
        for i in ["OneSided"]:
            conv_rate(Nx_values, approx=i)
        
        rinp = raw_input("Plot solution for a given Nx (y/n)? ")
        if rinp == "y":
            Nx = raw_input("Enter Nx value: ")
            plot_ex(int(Nx), approx="OneSided")

    
if __name__ == '__main__':    
    ex(raw_input("Which exercise (a/b/c/d)? "))
           




