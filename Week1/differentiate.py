from numpy import *
from pylab import *

def differentiate(u, dt):
	du = zeros(len(u))
	for i in range(1, (len(u)-1)):
		du[i] = (u[i+1]-u[i-1])/(2*dt)
	du[0] = (u[1]-u[0])/dt
	du[-1] = (u[-1] - u[-2])/dt
	return du

def vecdif(u, dt):
	du = zeros(len(u))
	du[1:-1] = (u[2:] - u[0:-2])/(2*dt)
	du[0] = (u[1]-u[0])/dt
	du[-1] = (u[-1] - u[-2])/dt
	return du

def test_differentiate():
	return t**2

print differentiate(test_differentiate(), dt)
print vecdif(test_differentiate(), dt)

