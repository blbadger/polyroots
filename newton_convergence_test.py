# libraries
import numpy as np 
import matplotlib.pyplot as plt 
plt.style.use('dark_background')
from Calculate import Calculate


def newton_raphson_map(equation, max_iterations, x_range, y_range):
	print (equation)
	# top left to bottom right
	y, x = np.ogrid[2: -2: y_range*1j, -2: 2: x_range*1j]
	z_array = x + y*1j

	iterations_until_rooted = max_iterations + np.zeros(z_array.shape)

	 # create a boolean table of all 'true'
	not_already_at_root = iterations_until_rooted < 10000

	for i in range(max_iterations):
		previous_z_array = z_array
		z = z_array
		f_now = Calculate(equation, z, differentiate=False).evaluate()
		f_prime_now = Calculate(equation, z, differentiate=True).evaluate()
		z_array = z_array - f_now / f_prime_now

		# the boolean map is tested for rooted values
		found_root = (abs(z_array - previous_z_array) < 0.0000001) & not_already_at_root
		iterations_until_rooted[found_root] = i
		not_already_at_root = np.invert(found_root) & not_already_at_root

	return iterations_until_rooted


def root_iterations(z, equation):
	z_now = Calculate(equation, z, differentiate=False).evaluate()
	z_prime = Calculate(equation, z, differentiate=True).evaluate()
	return z_now, z_prime
	
# number of iterations and array initialization
steps = 10000
y_range = 1400
x_range = 2000

Z = [0 for i in range(steps+1)]

# starting point
Z[0] = 1.8132 - 1.00405j

# add points to array
for i in range(steps):
	z_now, z_prime = root_iterations(Z[i], 'x^5-x-1')
	Z[i+1] = Z[i] - z_now/z_prime

# print (Z)
X = [i.real for i in Z]
Y = [i.imag for i in Z]

# plt.plot(X, Y, '-', color='white', alpha = 1, markersize=2)
# plt.imshow(newton_raphson_map('x^5-x-1', 30, 1558, 1558), extent=[-2, 2, -2, 2], cmap='inferno')
# plt.axis('off')
# plt.show()
# plt.close()

last, second_last, third_last = Z[-1], Z[-2], Z[-3] # assign periodic points

def convergence_to_period_3(equation, max_iterations, x_range, y_range, last, second_last, third_last):
	'''
	Determines which points converge on a period 3 trajectory
	defined by 'last, second_last, third_last' in the complex
	plane using Newton's root finiding method.  These points 
	need to be defined prior to calling this function.  Output 
	is a table corresponding to the first iteration coming close 
	to the specified orbit.
	'''
	print (equation)
	# top left to bottom right
	y, x = np.ogrid[2: -2: y_range*1j, -2: 2: x_range*1j]
	z_array = x + y*1j
	previous_z_array = z_array # initialize array for tracking previous values

	iterations_until_periodic = max_iterations + np.zeros(z_array.shape)

	 # create a boolean table of all 'true'
	not_already_periodic = iterations_until_periodic < 10000


	for i in range(max_iterations):
		double_previous_z_array = previous_z_array
		previous_z_array = z_array
		z = z_array
		f_now = Calculate(equation, z, differentiate=False).evaluate() 
		f_prime_now = Calculate(equation, z, differentiate=True).evaluate()
		z_array = z_array - f_now / f_prime_now

		# the boolean map is tested for rooted values
		periodic = (abs(z_array - last) < 0.0000001) & (abs(previous_z_array - second_last) < 0.0000001) & \
		(abs(double_previous_z_array - third_last) < 0.0000001) & not_already_periodic
		iterations_until_periodic[periodic] = i
		not_already_periodic = np.invert(periodic) & not_already_periodic

	return iterations_until_periodic

# plt.plot(X, Y, '^', color='white', alpha = 1, markersize=2)
# plt.imshow(convergence_to_period_3('x^5-x-1', 50, 1558, 1558, last, second_last, third_last), extent=[-2, 2, -2, 2], cmap='inferno')
# plt.axis('off')
# # plt.show()
# plt.savefig('convergence.png', bbox_inches='tight', dpi=420)
# plt.close()


# to overlay plots

# array_1 = convergence_to_period_3('x^5-x-1', 80, 1558, 1558, last, second_last, third_last)
# array_2 = newton_raphson_map('x^5-x-1', 80, 1558, 1558)

# plt.imshow(np.minimum(array_1, array_2), extent=[-2, 2, -2, 2], cmap='inferno')
# plt.axis('on')
# plt.show()
# # plt.savefig('convergence_overlayed.png', bbox_inches='tight', dpi=400)
# plt.close()


def traveling_together(equation, max_iterations, x_range, y_range):
	'''
	Returns points that stay near points nearby in future iteration
	'''
	y, x = np.ogrid[2: -2: y_range*1j, -2: 2: x_range*1j]
	z_array = x + y*1j
	z_array_2 = z_array + 0.0000001 # arbitrary change, can be any small amount
	iterations_until_together = max_iterations + np.zeros(z_array.shape)

	# create a boolean table of all 'true'
	not_already_together = iterations_until_together < 10000

	for i in range(max_iterations):
		f_now = Calculate(equation, z_array, differentiate=False).evaluate() 
		f_prime_now = Calculate(equation, z_array, differentiate=True).evaluate()
		z_array = z_array - f_now / f_prime_now

		f_2_now = Calculate(equation, z_array_2, differentiate=False).evaluate() 
		f_2_prime_now = Calculate(equation, z_array_2, differentiate=True).evaluate()		
		z_array_2 = z_array_2 - f_2_now / f_2_prime_now

		# the boolean map is tested for rooted values
		together = (abs(z_array - z_array_2) <= 0.0000001) & not_already_together
		iterations_until_together[together] = i
		not_already_together = np.invert(together) & not_already_together

	return iterations_until_together


plt.imshow(traveling_together('x^5-x-1', 40, 1558, 1558), extent=[-2, 2, -2, 2], cmap='inferno')
plt.axis('on')
plt.show()
# plt.savefig('together.png', bbox_inches='tight', dpi=420)
plt.close()
