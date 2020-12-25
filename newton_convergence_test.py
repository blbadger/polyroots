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

print (Z)
X = [i.real for i in Z]
Y = [i.imag for i in Z]

# plt.plot(X, Y, '-', color='white', alpha = 1, markersize=2)
# plt.imshow(newton_raphson_map('x^5-x-1', 30, 1558, 1558), extent=[-2, 2, -2, 2], cmap='inferno')
# plt.axis('off')
# plt.show()
# plt.close()

last, second_last, third_last = Z[-1], Z[-2], Z[-3] # assign periodic points

def convergence_to_period_3(equation, max_iterations, x_range, y_range, last, second_last, third_last):
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


# overlay plots

array_1 = convergence_to_period_3('x^5-x-1', 80, 1558, 1558, last, second_last, third_last)

print (array_1)

array_2 = newton_raphson_map('x^5-x-1', 80, 1558, 1558)

print (array_2)

plt.imshow(np.minimum(array_1, array_2), extent=[-2, 2, -2, 2], cmap='inferno')
plt.axis('off')
# plt.show()
plt.savefig('convergence_2.png', bbox_inches='tight', dpi=400)
plt.close()