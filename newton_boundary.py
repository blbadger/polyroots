#! python3

# import third party libraries
import numpy as np 
import matplotlib.pyplot as plt 

from Calculate import Calculate 

plt.style.use('dark_background')

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

	xl = -1.4/(2**(t/60)) + 0.032
	xr = 1.4/(2**(t/60)) + 0.032
	yl = 1/(2**(t/60))
	yr = -1/(2**(t/60)) 
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
		# print (i)
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


for t in range(600):
	print (t)
	plt.style.use('dark_background')
	plt.imshow(newton_boundary('x^5-x-1', 60 + t // 10, 2400, 1600, t), cmap='inferno')
	plt.axis('off')
	plt.savefig('Newton_boundary{0:03d}.png'.format(t), bbox_inches='tight', pad_inches=0, dpi=500)
	# plt.show()
	plt.close()


def newton_boundary_iterations(equation, max_iterations, x_range, y_range, t):
	"""
	Generates a newton-raphson fractal.

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

	y, x = np.ogrid[0.5: -0.45: y_range*1j, -0.3: 0.95: x_range*1j]
	a_array = x + y*1j
	z_array = np.zeros(a_array.shape)

	iterations_until_rooted = np.zeros(z_array.shape)

	 # create a boolean grid of all 'true'
	not_already_at_root = iterations_until_rooted < 10000

	nondiff = Calculate(equation, differentiate=False)
	diffed = Calculate(equation, differentiate=True)

	for i in range(max_iterations):
		iterations_until_rooted[not_already_at_root] = i
		# print (i)
		previous_z_array = z_array
		z = z_array
		f_now = nondiff.evaluate(z)
		f_prime_now = diffed.evaluate(z)
		z_array = z_array - f_now / f_prime_now + a_array

		# the boolean map is tested for rooted values
		found_root = (abs(z_array - previous_z_array) < 0.0000001) & not_already_at_root
		iterations_until_rooted[found_root] = i
		not_already_at_root = np.invert(found_root) & not_already_at_root

		print ('i' + str(set([i for i in iterations_until_rooted[700]])))
		plt.plot()
		plt.imshow(iterations_until_rooted, cmap='inferno')
		plt.axis('off')
		plt.savefig('Newton_boundary{0:03d}.png'.format(i), bbox_inches='tight', pad_inches=0, dpi=410)
		# plt.show()
		plt.close()

	return iterations_until_rooted


# newton_boundary_iterations('x^5-x-1', 40, 2200, 1400, 0)
