from wave1D_dn_vc_edited import *
import matplotlib.pylab as mpl

def find_f(L_value, q, exercise):
        """Function for finding the source term 
        for the given exact solution, using sympy."""

        from sympy import symbols, cos, pi, diff, simplify, lambdify
        xs, ts, Ls = symbols('x t L')
        
        u_e = lambda x, t: cos(np.pi*x/Ls)*cos(t)   #Exact solution.
        if exercise == "b":
            q   = lambda x: 1 + cos(np.pi*xs/Ls)
        
        #Exact solution twice differentiated wrt. t:
        utt = diff(u_e(xs,ts),ts,ts)
        #Exact solution differentiated wrt. x:    
        ux = diff(u_e(xs,ts),xs)   
        #Dx(q*Dx(u)):
        qux_x = diff(q(xs)*ux,xs)
        
        #Finding a sympy expression for f:
        f_ = simplify(utt - qux_x).subs(Ls,L_value)
        #Returning a callable function for f:
        return lambdify((xs,ts), f_, modules='numpy')

def ex(exercise):

    L = 1.
    u_exact = lambda x, t: np.cos(np.pi*x/L)*np.cos(t)
    
    #Defining q for the different exercises, for c) and d) I use
    #the q from a).
    if exercise == "a":
        q = lambda x: 1. + (x-L/2)**4
    elif exercise == "b":
        q = lambda x: 1. + np.cos(np.pi*x/L)
    else:
        q = lambda x: 1. + (x-L/2)**4
        
    c = lambda x: np.sqrt(q(x)) #wave velocity.
    I = lambda x: u_exact(x, 0) #Initial condition.
    V = None 
    U_0 = None
    U_L = None
    C = 1.
    #Nx = 100.
    #dt = C*(L/Nx)/c(0)
    T = 3.    
    f = find_f(L, q, exercise=exercise)

    def plot_ex(Nx, approx):
        """Function for plotting the exact 
        and numerical solutions together."""
        dt = C*(L/Nx)/c(0)
        solver(I, V, f, c, U_0, U_L, L, dt, C, T,
               user_action=PlotAndStoreSolution('moving_end'),
               version='vectorized',
               stability_safety_factor=1, approx=approx)
    
    def conv_rate(Nx_values, approx):
        """Function for finding the convergence rates."""
        
        def assert_no_error(u, x, t, n):
            """Function for finding the errors needed 
            to calculate the convergence rate."""
            u_e = u_exact(x, t[n]) #Exact solution.
            #Maximum difference in exact and numerical:
            diff = np.abs(u - u_e).max() 
            #Difference in exact and numerical:
            e = u_e - u
            dt = t[1]-t[0]
            #The error on the form used in the equation 
            #for finding the convergence rate:
            E = np.sqrt(dt*sum(e**2))
            #Saving E for each t value:
            Errors_t.append(E)
        
        Error_t_max = []
        dt_values = []
           
        for i in range(len(Nx_values)):
            dt = float(C*(L/Nx_values[i])/c(0))
            dt_values.append(dt)
        
            Errors_t = []
            #Calling the solver with the user actions 
            #in assert_no_error:
            t, cpu_time = solver(I, V, f, c, U_0, U_L, L, dt, C, T,
                                user_action=assert_no_error,
                                version='vectorized',
                                stability_safety_factor=1,
                                approx=approx)
            #Saving the maximum error for each Nx value:
            Error_t_max.append(max(Errors_t))
        
        #r is the convergence rate:    
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
        """Function for computing the convergence rate."""
        m = len(dt_values)
        #Convergence rate:
        r = [np.log(E_values[i-1]/E_values[i])/
        np.log(dt_values[i-1]/dt_values[i])
        for i in range(1, m, 1)]
        #Round to two decimals
        r = [round(r_, 5) for r_ in r]
        return r
    
    #Run the program for multiple Nx values:
    Nx_values = np.linspace(100, 1000, 10)
    
    #Calling the necessary functions for the given exercise.
    #Standard is eq(57), Approx is eq(54), OneSided is for c) and 
    #GhostCell is for d).
    if exercise == "a" or exercise == "b":
        for i in ["Standard", "Approx"]:#, "OneSided"]:
            conv_rate(Nx_values, approx=i)
        
        rinp = raw_input("Plot solution for a given Nx (y/n)? ")
        if rinp == "y":
            Nx = raw_input("Enter Nx value: ")
            approx = raw_input("Enter desired scheme (Standard/Approx): ")
            plot_ex(int(Nx), approx)
            
    elif exercise == "c":
        conv_rate(Nx_values, approx="OneSided")
        
        rinp = raw_input("Plot solution for a given Nx (y/n)? ")
        if rinp == "y":
            Nx = raw_input("Enter Nx value: ")
            plot_ex(int(Nx), approx="OneSided")
            
    elif exercise == "d":
        conv_rate(Nx_values, approx="GhostCell")
        
        rinp = raw_input("Plot solution for a given Nx (y/n)? ")
        if rinp == "y":
            Nx = raw_input("Enter Nx value: ")
            plot_ex(int(Nx), approx="GhostCell")

    
if __name__ == '__main__':    
    ex(raw_input("Which exercise (a/b/c/d)? "))
           




