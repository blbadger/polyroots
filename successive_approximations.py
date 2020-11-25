def successive_approximations(x_start, iterations):
	'''
	Computes the successive approximations for the equation
	z = z^3 - 1 in the complex plane.
	'''
	x_now = x_start
	ls = []
	for i in range(iterations):
		f_now = x_now**3 - 1
		f_prime_now = 3*x_now**2 

		x_next = x_now - f_now / f_prime_now
		ls.append(x_now)
		x_now = x_next
	return ls

print (successive_approximations(-1j - 1, 30))