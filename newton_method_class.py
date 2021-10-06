#! python3

# import third party libraries
import numpy as np 
import matplotlib.pyplot as plt 

from Calculate import Calculate 


def newton_raphson_map(equation, max_iterations, x_range, y_range, t):
	"""
	Generates a newton-raphson fractal.

	Args:
		equation: str, equation of interest
		max_iterations: int, number of iterations 
		x_range: int, number of real values per output
		y_range: int, number of imaginary values per output
		t: int

	Returns:
		iterations_until_rooted: np.arr (2D) of iterations until a root is found
			at each point in y_range and x_range

	"""

	print (equation)
	y, x = np.ogrid[0.7: -0.7: y_range*1j, -1.1: 1.1: x_range*1j]
	z_array = x + y*1j

	iterations_until_rooted = max_iterations + np.zeros(z_array.shape)

	 # create a boolean gridof all 'true'
	not_already_at_root = iterations_until_rooted < 10000

	nondiff = Calculate(equation, differentiate=False)
	diffed = Calculate(equation, differentiate=True)


	for i in range(max_iterations):
		print (i)
		previous_z_array = z_array
		z = z_array
		f_now = nondiff.evaluate(z)
		f_prime_now = diffed.evaluate(z)
		z_array = z_array - f_now / f_prime_now

		# the boolean map is tested for rooted values
		found_root = (abs(z_array - previous_z_array) < 0.0000001) & not_already_at_root
		iterations_until_rooted[found_root] = i
		not_already_at_root = np.invert(found_root) & not_already_at_root

	return iterations_until_rooted

t = 0
for t in range(1, 2000):
	plt.style.use('dark_background')
	plt.imshow(newton_raphson_map('x^(7.11+' + str(t*0.1/100000) + 'i)-x^(1+' + str(t*0.1/100000) + 'i)-1', 40, 2200, 1400, t), cmap='inferno')
	plt.axis('off')
	plt.savefig('Newton{0:03d}.png'.format(t), bbox_inches='tight', pad_inches=0, dpi=410)
	# plt.show()
	plt.close()

# No spaces are allowed in the input polynomial! Another example of a valid complex-valued input:
'(0.2-3i)x^(2j)-x-(1-2i)'



