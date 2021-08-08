import numpy as np 

def successive_approximations(x_start, iterations):
	'''
	Computes the successive approximations for the equation
	z = z^3 - 1 in the complex plane.
	'''
	x_now = x_start
	ls = []
	for i in range(iterations):
		f_now = x_now**5-x_now-1
		f_prime_now = 5*x_now**4 - 1 

		x_next = x_now - f_now / f_prime_now
		ls.append(x_now)
		x_now = x_next
	return ls

entry = [0.35140075+0.78236245j]
entry = [0.3514151 +0.7823557j]
# entry = [-0.5]

for i in entry:
	print (successive_approximations(i, 20))
m1 = [[1, 3, 5],
	[ 2, 4, 6]]

m2 = [[1 , 4, 3], 
	[0 , 5, 1], 
	[17, 6, 11]]

m2 = np.transpose(m2)
print (np.matmul(m1, m2))
