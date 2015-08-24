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
