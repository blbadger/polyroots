# rootmethod.py

# import standard libraries
import io
import base64

# impot third party libraries
import numpy as np 
import matplotlib.pyplot as plt 
from CalculateFaster import OptiCalculate
from CoCalculate import ComplexCalculate
import numexpr as ne


def convert_to_binary(arr, cmap):
	"""
	Encodes a numpy cmap as a binary array

	Args:
		arr: np.ogrid[int]
		cmap: string, color map choice for np.plt.imsave()

	Returns: 
		str, base64-encoded string for html loading

	"""
	buf = io.BytesIO()
	plt.style.use('dark_background')
	plt.imsave(buf, arr, cmap=cmap, format='png')
	data = base64.b64encode(buf.getbuffer()).decode("utf8") 

	return "data:image/png;base64,{}".format(data)


def return_roots(z_array, not_already_at_root):
	"""
	Returns the roots found in z_array

	Args:
		z_array: np.ogrid[complex] object
		not_already_at_root: np.opgrid[bool] object

	Returns:
		roots_arr: arr[complex] of found roots

	"""

	found_root = np.invert(not_already_at_root)
	z_array = np.around(z_array, 4)

	roots = set()
	for i in z_array[found_root]:
		roots.add(i)

	roots_arr = list(roots)
	roots_arr.sort()
	return roots_arr


def newton(equation, max_iterations, x_range, y_range, res_value, cmap):
	"""
	Newton's method, compatible with any complex numbered 'equation'

	Args:
		equation: str
		max_iterations: int, number of times Newton's method is applied
		x_range: arr[int], real axis bounds
		y_range: arr[int], imaginary axis bounds
		res_value: arr[int], resolution [x_resolution, y_resolution]
		cmap: str, color map choice for np.plt.imshow()

	Returns:
		bin_arr: str, base64-encoded binary array
		roots: arr[complex] of roots found

	"""

	y, x = np.ogrid[float(y_range[1]): float(y_range[0]): res_value[1]*1j, \
					float(x_range[0]): float(x_range[1]): res_value[0]*1j]
	z_array = x + y*1j

	iterations_until_rooted = max_iterations + np.zeros(z_array.shape)

	# create a boolean table of all 'true'
	not_already_at_root = iterations_until_rooted < 10000

	# initialization of Calculate object
	nondiff = ComplexCalculate(equation, differentiate=False)
	diff = ComplexCalculate(equation, differentiate=True)

	for i in range(max_iterations):
		previous_z_array = z_array
		z = z_array
		f_now = nondiff.evaluate(z)
		f_prime_now = diff.evaluate(z)
		z_array = z_array - f_now / f_prime_now

		# the boolean map is tested for rooted values
		found_root = (abs(z_array - previous_z_array) < 0.0000001) & not_already_at_root
		iterations_until_rooted[found_root] = i
		not_already_at_root = np.invert(found_root) & not_already_at_root
		
	arr = iterations_until_rooted
	bin_arr = convert_to_binary(arr, cmap)
	roots = return_roots(z_array, not_already_at_root)

	return bin_arr, roots


def newton_optimized(equation, max_iterations, x_range, y_range, res_value, cmap):
	"""
	Newton's method using a numexpr-optimized Calculate class.  Only for use
	with real-valued 'equation'

	Args:
		equation: str
		max_iterations: int, number of times Newton's method is applied
		x_range: arr[int], real axis bounds
		y_range: arr[int], imaginary axis bounds
		res_value: arr[int], resolution [x_resolution, y_resolution]
		cmap: str, color map choice for np.plt.imshow()

	Returns:
		bin_arr: str, base64-encoded binary array
		roots: arr[complex] of roots found

	"""

	y, x = np.ogrid[float(y_range[1]): float(y_range[0]): res_value[1]*1j, \
					float(x_range[0]): float(x_range[1]): res_value[0]*1j]
	z_array = x + y*1j

	iterations_until_rooted = max_iterations + np.zeros(z_array.shape)

	# create a boolean table of all 'true'
	not_already_at_root = iterations_until_rooted < 10000

	# initialization of OptiCalculate object
	nondiff = OptiCalculate(equation, differentiate=False)
	diff = OptiCalculate(equation, differentiate=True)

	for i in range(max_iterations):
		previous_z_array = z_array
		z = z_array
		f_now = nondiff.evaluate(z)
		f_prime_now = diff.evaluate(z)
		z_array = ne.evaluate('z_array - f_now / f_prime_now')

		# the boolean map is tested for rooted values
		found_root = (abs(z_array - previous_z_array) < 0.0000001) & not_already_at_root
		iterations_until_rooted[found_root] = i
		not_already_at_root = np.invert(found_root) & not_already_at_root
		
	arr = iterations_until_rooted
	bin_arr = convert_to_binary(arr, cmap)
	roots = return_roots(z_array, not_already_at_root)

	return bin_arr, roots


def halley(equation, max_iterations, x_range, y_range, res_value, cmap):
	"""
	Halley's method, compatible with any complex numbered 'equation' arg

	Args:
		equation: str
		max_iterations: int, number of times Newton's method is applied
		x_range: arr[int], real axis bounds
		y_range: arr[int], imaginary axis bounds
		res_value: arr[int], resolution [x_resolution, y_resolution]
		cmap: str, color map choice for np.plt.imshow()

	Returns:
		bin_arr: str, base64-encoded binary array
		roots: arr[complex] of roots found

	"""

	y, x = np.ogrid[float(y_range[1]): float(y_range[0]): res_value[1]*1j, \
					float(x_range[0]): float(x_range[1]): res_value[0]*1j]
	z_array = x + y*1j

	iterations_until_rooted = max_iterations + np.zeros(z_array.shape)

	# create a boolean table of all 'true'
	not_already_at_root = iterations_until_rooted < 10000

	nondiff = ComplexCalculate(equation, differentiate=False)
	diff = ComplexCalculate(equation, differentiate=True)

	# double derivative initialization
	diff_string = diff.to_string()
	double_diff = ComplexCalculate(diff_string, differentiate=True)

	for i in range(max_iterations):
		previous_z_array = z_array
		z = z_array
		
		f_now = nondiff.evaluate(z)
		f_prime_now = diff.evaluate(z) # first derivative evaluation
		f_double_prime_now = double_diff.evaluate(z) # second derivative evaluation
		z_array = z - (2*f_now * f_prime_now / (2*(f_prime_now)**2 - f_now * f_double_prime_now))

		# test the boolean map for rooted values
		found_root = (abs(z_array - previous_z_array) < 0.000000001) & not_already_at_root
		iterations_until_rooted[found_root] = i
		not_already_at_root = np.invert(found_root) & not_already_at_root

	arr = iterations_until_rooted
	bin_arr = convert_to_binary(arr, cmap)
	roots = return_roots(z_array, not_already_at_root)

	return bin_arr, roots


def halley_optimized(equation, max_iterations, x_range, y_range, res_value, cmap):
	"""
	Halley's method using an optimized Calculate class via numexpr.  Only for
	real-valued 'equation' args.

	Args:
		equation: str
		max_iterations: int, number of times Newton's method is applied
		x_range: arr[int], real axis bounds
		y_range: arr[int], imaginary axis bounds
		res_value: arr[int], resolution [x_resolution, y_resolution]
		cmap: str, color map choice for np.plt.imshow()

	Returns:
		bin_arr: str, base64-encoded binary array
		roots: arr[complex] of roots found

	"""

	y, x = np.ogrid[float(y_range[1]): float(y_range[0]): res_value[1]*1j, \
					float(x_range[0]): float(x_range[1]): res_value[0]*1j]
	z_array = x + y*1j

	iterations_until_rooted = max_iterations + np.zeros(z_array.shape)

	 # create a boolean table of all 'true'
	not_already_at_root = iterations_until_rooted < 10000

	#
	nondiff = OptiCalculate(equation, differentiate=False)
	diff = OptiCalculate(equation, differentiate=True)

	# double derivative initialization
	diff_string = diff.to_string()
	double_diff = OptiCalculate(diff_string, differentiate=True)


	for i in range(max_iterations):
		previous_z_array = z_array
		z = z_array
		
		f_now = nondiff.evaluate(z)
		f_prime_now = diff.evaluate(z) # first derivative evaluation
		f_double_prime_now = double_diff.evaluate(z) # second derivative evaluation
		z_array = ne.evaluate('z - (2*f_now * f_prime_now / (2*(f_prime_now)**2 - f_now * f_double_prime_now))')

		# test the boolean map for rooted values
		found_root = (abs(z_array - previous_z_array) < 0.000000001) & not_already_at_root
		iterations_until_rooted[found_root] = i
		not_already_at_root = np.invert(found_root) & not_already_at_root

	arr = iterations_until_rooted
	bin_arr = convert_to_binary(arr, cmap)
	roots = return_roots(z_array, not_already_at_root)

	return bin_arr, roots


def secant(equation, max_iterations, x_range, y_range, res_value, cmap):
	"""
	Secant method, compatible with any complex numbered 'equation' arguments

	Args:
		equation: str
		max_iterations: int, number of times Newton's method is applied
		x_range: arr[int], real axis bounds
		y_range: arr[int], imaginary axis bounds
		res_value: arr[int], resolution [x_resolution, y_resolution]
		cmap: str, color map choice for np.plt.imshow()

	Returns:
		bin_arr: str, base64-encoded binary array
		roots: arr[complex] of roots found

	"""

	y, x = np.ogrid[float(y_range[1]): float(y_range[0]): res_value[1]*1j, \
					float(x_range[0]): float(x_range[1]): res_value[0]*1j]
	z_array = x + y*1j
	iterations_until_rooted = max_iterations + np.zeros(z_array.shape)

	 # create a boolean table of all 'true'
	not_already_at_root = iterations_until_rooted < 10000
	zeros = np.zeros(z_array.shape) 
	z_0 = (z_array - zeros)/2 # setting the initial guess to half the distance to the origin from the second guess, which is plotted

	# initializatoin of ComplexCalculate object
	nondiff = ComplexCalculate(equation, differentiate=False)

	for i in range(max_iterations):
		previous_z_array = z_array
		z = z_array
		f_previous = nondiff.evaluate(z_0)
		f_now = nondiff.evaluate(z)
		z_array = z - f_now * (z - z_0)/(f_now - f_previous)

		# the boolean map is tested for rooted values
		found_root = (abs(z_array - previous_z_array) < 0.0000001) & not_already_at_root
		iterations_until_rooted[found_root] = i
		not_already_at_root = np.invert(found_root) & not_already_at_root
		z_0 = z 

	arr = iterations_until_rooted
	bin_arr = convert_to_binary(arr, cmap)

	return bin_arr, []


def secant_optimized(equation, max_iterations, x_range, y_range, res_value, cmap):
	"""
	Secant method with optimized Calculate class via numexpr. For use only with
	real-valued 'equation' arguments.

	Args:
		equation: str
		max_iterations: int, number of times Newton's method is applied
		x_range: arr[int], real axis bounds
		y_range: arr[int], imaginary axis bounds
		res_value: arr[int], resolution [x_resolution, y_resolution]
		cmap: str, color map choice for np.plt.imshow()

	Returns:
		bin_arr: str, base64-encoded binary array
		roots: arr[complex] of roots found

	"""
	y, x = np.ogrid[float(y_range[1]): float(y_range[0]): res_value[1]*1j, \
					float(x_range[0]): float(x_range[1]): res_value[0]*1j]
	z_array = x + y*1j
	iterations_until_rooted = max_iterations + np.zeros(z_array.shape)

	# create a boolean table of all 'true'
	not_already_at_root = iterations_until_rooted < 10000
	zeros = np.zeros(z_array.shape) 
	z_0 = ne.evaluate('(z_array - zeros)/2') # setting the initial guess to half the distance to the origin from the second guess, which is plotted

	# initialize OptiCalculate object
	nondiff = OptiCalculate(equation, differentiate=False)

	for i in range(max_iterations):
		previous_z_array = z_array
		z = z_array
		f_previous = nondiff.evaluate(z_0)
		f_now = nondiff.evaluate(z)
		z_array = ne.evaluate('z - f_now * (z - z_0)/(f_now - f_previous)')

		# the boolean map is tested for rooted values
		found_root = (abs(z_array - previous_z_array) < 0.0000001) & not_already_at_root
		iterations_until_rooted[found_root] = i
		not_already_at_root = np.invert(found_root) & not_already_at_root
		z_0 = z 

	arr = iterations_until_rooted
	bin_arr = convert_to_binary(arr, cmap)
	
	return bin_arr, []



