from numpy import *
import matplotlib.pyplot as plt

def simulate(
    beta=0.9,                 # dimensionless parameter
    Theta=30,                 # initial angle in degrees
    epsilon=0,                # initial stretch of wire
    num_periods=6,            # simulate for num_periods
    time_steps_per_period=60, # time step resolution
    plot=True,                # make plots or not
    ):

    P = 2*pi
    T = num_periods*P
    dt = float(P/time_steps_per_period)
    Nt = int(round(T/dt))
    x = zeros(Nt+1)
    y = zeros(Nt+1)
    t = linspace(0, Nt*dt, Nt+1)
    Theta_Array = zeros(Nt+1)

    x[0] = (1+epsilon)*sin(Theta)
    y[0] = 1 - (1+epsilon)*cos(Theta)
    L = sqrt(x[0]**2 + (y[0] - 1)**2)
    x[1] = (1 - dt**2/2.*beta/(1-beta)*(1-beta/L))*x[0]
    y[1] = y[0] - dt**2/2.*beta/(1-beta)*(1-beta/L)*(y[0]-1) - dt**2/2.*beta
    Theta_Array[0] = Theta
    Theta_Array[1] = arctan(x[1]/(1-y[1]))

    for n in range(1, Nt):
        L = sqrt(x[n]**2 + (y[n] - 1)**2)
        x[n+1] = (2 - dt**2*beta/(1-beta)*(1-beta/L))*x[n] - x[n-1]
        y[n+1] = 2*y[n] - dt**2*beta/(1-beta)*(1-beta/L)*(y[n]-1) - y[n-1] - dt**2*beta
        Theta_Array[n+1] = arctan(x[n+1]/(1-y[n+1]))
    
    if plot == True:
        plt.plot(x, y)
        plt.figure()
        plt.plot(t, Theta_Array)
        plt.show()
        plt.gca().set_aspect('equal')
        if Theta < 10:
            beta = 1.
            #x[0] = (1+epsilon)*sin(Theta)
            #y[0] = 1 - (1+epsilon)*cos(Theta)
            #L = sqrt(x[0]**2 + (y[0] - 1)**2)
            x[1] = (1 - dt**2/2.*1./L)*x[0]
            y[1] = y[0] - dt**2/2.*1./L*(y[0]-1) - dt**2/2.*beta
            Theta_Array2 = zeros(Nt+1)
            Theta_Array2[0] = Theta
            Theta_Array2[1] = arctan(x[1]/(1-y[1]))

            for n in range(1, Nt):
                L = sqrt(x[n]**2 + (y[n] - 1)**2)
                x[n+1] = (2 - dt**2*1./L)*x[n] - x[n-1]
                y[n+1] = 2*y[n] - dt**2*1./L*(y[n]-1) - y[n-1] - dt**2*beta
                Theta_Array2[n+1] = arctan(x[n+1]/(1-y[n+1]))

            plt.plot(t, Theta_Array,
                    t, Theta_Array2)
            plt.show()

simulate(beta=0.9, Theta=9)

#def test_function()





#return x, y, t, Theta_Array








