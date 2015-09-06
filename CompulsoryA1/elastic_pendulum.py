from numpy import *
import numpy as np
import matplotlib.pyplot as plt

def simulate(
    beta=0.9,                 # dimensionless parameter
    Theta=30,                 # initial angle in degrees
    epsilon=0,                # initial stretch of wire
    num_periods=6,            # simulate for num_periods
    time_steps_per_period=60, # time step resolution
    plot=True,                # make plots or not
    anim=True,                # animate pendulum movement or not
    ):

    P = 2.*pi
    DegToRad = pi/180.
    Theta = DegToRad*Theta
    T = num_periods*P
    #T = T/sqrt(beta/(1.-beta))
    dt = float(P/time_steps_per_period)
    #dt = dt/sqrt(beta/(1.-beta))
    Nt = int(round(T/dt))
    x = zeros(Nt+1)
    y = zeros(Nt+1)
    t = linspace(0, Nt*dt, Nt+1)
    Theta_Array = zeros(Nt+1)
    Const = dt**2*beta/(1.-beta)

    #Initial conditions
    x[0] = (1+epsilon)*sin(Theta)
    y[0] = 1 - (1+epsilon)*cos(Theta)
    L = sqrt(x[0]**2 + (y[0] - 1)**2)
    x[1] = (1 - Const/2.*(1-beta/L))*x[0]
    y[1] = y[0] - Const/2.*(1-beta/L)*(y[0]-1) - dt**2/2.*beta
    Theta_Array[0] = Theta
    Theta_Array[1] = arctan(x[1]/(1-y[1]))
    
    #Preparing animation
    if anim == True:
        plt.figure()
        #if Theta < 10.*2.*pi/360:
            #plt.axis([-.20,.20,-.01,.02])
        #elif Theta > 30.*2.*pi/360:
            #plt.axis([-1,1,-1,1])
        #else:
            #plt.axis([-.6,.6,-.1,.2])
        plt.plot(x[:1], y[:1])
        plt.ion()
        plt.show()

    for n in range(1, Nt):
        #Solving ODEs
        L = sqrt(x[n]**2 + (y[n] - 1)**2)
        x[n+1] = (2 - Const*(1-beta/L))*x[n] - x[n-1]
        y[n+1] = 2*y[n] - Const*(1-beta/L)*(y[n]-1) - y[n-1] - dt**2*beta
        Theta_Array[n+1] = arctan(x[n+1]/(1-y[n+1]))

        #Animating
        if anim == True:
            plt.clf()
            xl = [0, x[n]]
            yl = [1, y[n]]
            #if Theta < 10.*2.*pi/360:
                #plt.axis([-.20,.20,-.01,.02])
            #elif Theta > 30.*2.*pi/360:
                #plt.axis([-1,1,-1,1])
            #else:
                #plt.axis([-.6,.6,-.1,.2])
            plt.plot(x[:n+1], y[:n+1],
                    xl, yl, 'k-',
                    x[n], y[n], 'ro', markersize=15)
            plt.legend(['Drawn path', 'Elastic wire', 'Pendulum bob'])
            plt.ylabel("y [dimensionless]")
            plt.xlabel("x [dimensionless]")
            plt.draw()

    #Plotting solution for position (x,y)
    if plot == True:

        #Sorting out figures and legends for animation on/off
        if anim == False:
            plt.figure()
        plt.plot(x, y)
        plt.legend(['Drawn path', 'Elastic wire', 'Pendulum bob', 'Full path'])
        if anim == False:
            plt.legend(['Full path'])
        plt.ylabel("y [dimensionless]")
        plt.xlabel("x [dimensionless]")

        #Plotting theta(t) if it is non zero
        if Theta != 0:
            plt.figure()
            plt.plot(t, Theta_Array*1./DegToRad)
            plt.legend([r"$\theta$(t)"])
            plt.ylabel(r"$\theta$(t) [Deg]")
            plt.xlabel("t [dimensionless]")
            plt.gca().set_aspect('equal')

            #Plotting theta comparison for the elastic and inelastic case
            if Theta < 10.*DegToRad:# and Theta != 0:
                beta = 1.
                #x[0] = (1+epsilon)*sin(Theta)
                #y[0] = 1 - (1+epsilon)*cos(Theta)
                #L = sqrt(x[0]**2 + (y[0] - 1)**2)

                #Initial conditions for the inelastic case
                x[1] = (1 - dt**2/2.*1./L)*x[0]
                y[1] = y[0] - dt**2/2.*1./L*(y[0]-1) - dt**2/2.*beta
                Theta_Array2 = zeros(Nt+1)
                Theta_Array2[0] = Theta
                Theta_Array2[1] = arctan(x[1]/(1-y[1]))

                #Solving the ODEs for the inelastic case
                for n in range(1, Nt):
                    L = sqrt(x[n]**2 + (y[n] - 1)**2)
                    x[n+1] = (2 - dt**2*1./L)*x[n] - x[n-1]
                    y[n+1] = 2*y[n] - dt**2*1./L*(y[n]-1) - y[n-1] - dt**2*beta
                    Theta_Array2[n+1] = arctan(x[n+1]/(1-y[n+1]))

                #Plot of theta(t) comparison
                plt.figure()
                plt.plot(t, Theta_Array*1./DegToRad,
                    t, Theta_Array2*1./DegToRad)
                plt.legend([r"$\theta$(t) elastic", r"$\theta$(t) non-elastic"])
                plt.ylabel(r"$\theta$(t) [Deg]")
                plt.xlabel("t [dimensionless]")
                plt.show()


    return x, y, t, Theta_Array, Nt

def test_function_b(
    beta=0.9, 
    Theta=0, 
    epsilon=0, 
    num_periods=6, 
    time_steps_per_period=60, 
    plot=False, 
    anim=False):

    #Test function to verify that x=0 and y=0 for all times when
    #theta=0 and epsilon=0.

    x, y, t, Theta_Array, Nt = simulate(beta, Theta,  epsilon, num_periods, time_steps_per_period, plot, anim)
    #Checks if there are any non-zero elements in x and y:
    assert not x.any(), not y.any()

def test_function_c(
    beta=0.9, 
    Theta=0, 
    epsilon=1, 
    num_periods=1, 
    time_steps_per_period=6000, 
    plot=True, 
    anim=False):

    #Test function for pure vertical motion (theta=0, epsilon!=0)

    x, y, t, Theta_Array, Nt = simulate(beta, Theta,  epsilon, num_periods, time_steps_per_period, plot, anim)

    #Exact solution is y=A*cos(w*t) where:
    A = -epsilon
    #and:
    w = sqrt(beta/(1.-beta))
    
    #Finding the error in the numerical solution.
    y_e = A*cos(w*t)
    tol = 1E-5
    error = abs(y_e - y)
    #print error.max()
    assert error.max() < tol

def demo(beta, Theta):
    
    simulate(beta, Theta, num_periods=3, time_steps_per_period=600, plot=True , anim=False)

if __name__ == '__main__':

    #simulate(beta=0.9, Theta=9, num_periods=1, plot=True , anim=True)

    test_function_b()

    test_function_c()

    demo(0.9, 10)

    plt.show(block=True)













