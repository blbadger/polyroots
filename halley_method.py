# libraries
import numpy as np 
import matplotlib.pyplot as plt 
from Calculate import Calculate 
from optimized.CalculateFaster import OptiCalculate


def halley_method(equation, max_iterations, x_range, y_range, t):
	"""
	Returns an array of the number of iterations until a root is found
	using Halley's method in the complex plane.

	Args:
		equation: str, polynomial of interest
		max_iterations: int of iterations
		x_range: int, number of real values per output
		y_range: int, number of imaginary values per output
		t: int

	Returns:
		iterations_until_together: np.arr[int] (2D) 
		
	"""

	# top left to bottom right
	y, x = np.ogrid[1: -1: y_range*1j, -1: 1: x_range*1j]
	z_array = x + y*1j

	iterations_until_rooted = max_iterations + np.zeros(z_array.shape)

	# create a boolean table of all 'true'
	not_already_at_root = iterations_until_rooted < 10000

	# initialize calculate object
	nondiff = Calculate(equation, differentiate=False)
	diffed = Calculate(equation, differentiate=True)

	# second derivative calculation
	diff_string = diffed.to_string()
	double_diffed = Calculate(diff_string, differentiate=True)

	for i in range(max_iterations):
		previous_z_array = z_array
		z = z_array

		f_now = nondiff.evaluate(z)
		f_prime_now = diffed.evaluate(z) # first derivative evaluation

		f_double_prime_now = double_diffed.evaluate(z) # second derivative evaluation

		z_array = z - (2*f_now * f_prime_now / (2*(f_prime_now)**2 - f_now * f_double_prime_now))

		# test the boolean map for rooted values
		found_root = (abs(z_array - previous_z_array) < 0.000000001) & not_already_at_root
		iterations_until_rooted[found_root] = i
		not_already_at_root = np.invert(found_root) & not_already_at_root

	return iterations_until_rooted



plt.style.use('dark_background')
plt.imshow(halley_method('x^13-x-1', 30, 1558, 1558, 30), extent=[-1, 1, -1, 1], cmap='inferno')
plt.axis('on')
plt.show()
plt.close()

# for incrementation
# for i in range(1):
# 	t = i
# 	plt.imshow(halley_method('x^13-x-1', 25, 2000, 1400, t), cmap='inferno')
# 	plt.axis('off')
# 	plt.show()
# 	# plt.savefig('halley{0:03d}.png'.format(i), bbox_inches='tight', dpi=420)
# 	plt.close()
