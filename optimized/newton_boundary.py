# libraries
import numpy as np 
import matplotlib.pyplot as plt 
from Calculate import Calculate # real-value version of Calculate
from CalculateFaster import OptiCalculate
import numexpr as ne

def newton_boundary_optimized(equation, max_iterations, x_range, y_range, t):
	print (equation)
	xl = -1.4/(2**(t/60)) + 0.004559235
	xr = 1.4/(2**(t/60)) + 0.004559235
	yl = 1/(2**(t/60)) + 0.00113144
	yr = -1/(2**(t/60)) + 0.00113144
	y, x = np.ogrid[yl: yr: y_range*1j, xl: xr: x_range*1j]
	a_array = x + y*1j
	z_array = np.zeros(a_array.shape)

	iterations_until_rooted = max_iterations + np.zeros(z_array.shape)

	 # create a boolean table of all 'true'
	not_already_at_root = iterations_until_rooted < 10000

	nondiff = OptiCalculate(equation, differentiate=False)
	diff = OptiCalculate(equation, differentiate=True)

	for i in range(max_iterations):
		previous_z_array = z_array
		z = z_array
		f_now = nondiff.evaluate(z)
		f_prime_now = diff.evaluate(z)
		z_array = ne.evaluate('z_array - f_now / f_prime_now + a_array')

		# the boolean map is tested for rooted values
		found_root = (abs(z_array - previous_z_array) < 0.0000001) & not_already_at_root
		iterations_until_rooted[found_root] = i
		not_already_at_root = np.invert(found_root) & not_already_at_root

	return iterations_until_rooted


def newton_boundary(equation, max_iterations, x_range, y_range, t):
	"""
	Generates the boundary of a Newton fractal.

	Args:
		equation: str, equation of interest
		max_iterations: int, number of iterations 
		x_range: int, number of real values per output
		y_range: int, number of imaginary values per output
		t: int, identifier for image construction

	Returns:
		iterations_until_rooted: np.arr (2D) of iterations until a root is found
			at each point in y_range and x_range

	"""
	xl = -1.4/(2**(t/60)) + 0.004559235
	xr = 1.4/(2**(t/60)) + 0.004559235
	yl = 1/(2**(t/60)) + 0.00113144
	yr = -1/(2**(t/60)) + 0.00113144
	y, x = np.ogrid[yl: yr: y_range*1j, xl: xr: x_range*1j]
	a_array = x + y*1j
	z_array = np.zeros(a_array.shape)

	iterations_until_rooted = np.zeros(z_array.shape)

	 # create a boolean grid of all 'true'
	not_already_at_root = iterations_until_rooted < 10000

	nondiff = Calculate(equation, differentiate=False)
	diffed = Calculate(equation, differentiate=True)

	for i in range(max_iterations):
		iterations_until_rooted[not_already_at_root] = i
		previous_z_array = z_array
		z = z_array
		f_now = nondiff.evaluate(z)
		f_prime_now = diffed.evaluate(z)
		z_array = z_array - f_now / f_prime_now + a_array

		# the boolean map is tested for rooted values
		found_root = (abs(z_array - previous_z_array) < 0.0000001) & not_already_at_root
		iterations_until_rooted[found_root] = i
		not_already_at_root = np.invert(found_root) & not_already_at_root

	return iterations_until_rooted


if __name__ == '__main__':
	for t in range(1501):
		print (t)
		plt.style.use('dark_background')
		xl = -1.4/(2**(t/60)) + 0.004559235
		xr = 1.4/(2**(t/60)) + 0.004559235
		yl = 1/(2**(t/60)) + 0.00113144
		yr = -1/(2**(t/60)) + 0.00113144
		plt.imshow(newton_boundary_optimized('x^5-x-1', 50 + t // 30, 1920, 1080, t), cmap='inferno')
		plt.axis('off')
		plt.savefig('Newton_boundary{0:04d}.png'.format(t), bbox_inches='tight', pad_inches=0, dpi=500)
		# plt.show()
		plt.close()
