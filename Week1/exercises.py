from numpy import *
from pylab import *

def mesh_function(f, t):
	f_array = zeros(len(t))
	for i in range(len(t)):
		f_array[i] = f(t[i])
	return f_array

def f(t):
	if 0<=t<=3:
		return exp(-t)
	elif 3<t<=4:
		return exp(-3*t)

dt = 0.1
T = 4
t = linspace(0, T, T/dt+1)
mesh = mesh_function(f, t)
print mesh
plot(t, mesh)
show()

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


plot(t, differentiate(mesh, dt))
show()

