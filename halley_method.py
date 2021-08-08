# libraries
import numpy as np 
import matplotlib.pyplot as plt 
from Calculate import Calculate 
from optimized.CalculateFaster import OptiCalculate


def halley_method(equation, max_iterations, x_range, y_range, t):
	print (equation)
	# top left to bottom right
	y, x = np.ogrid[5: -5: y_range*1j, -5: 5: x_range*1j]
	z_array = x + y*1j

	iterations_until_rooted = max_iterations + np.zeros(z_array.shape)

	 # create a boolean table of all 'true'
	not_already_at_root = iterations_until_rooted < 10000

	for i in range(max_iterations):
		previous_z_array = z_array
		z = z_array
		
		diff = Calculate(equation, z, differentiate=True)
		nondiff = Calculate(equation, z, differentiate=False)
		f_now = nondiff.evaluate()
		f_prime_now = diff.evaluate() # first derivative evaluation
		diff_string = diff.to_string()

		double_diff = Calculate(diff_string, z, differentiate=True)
		f_double_prime_now = double_diff.evaluate() # second derivative evaluation
		z_array = z - (2*f_now * f_prime_now / (2*(f_prime_now)**2 - f_now * f_double_prime_now))

		# test the boolean map for rooted values
		found_root = (abs(z_array - previous_z_array) < 0.000000001) & not_already_at_root
		iterations_until_rooted[found_root] = i
		not_already_at_root = np.invert(found_root) & not_already_at_root

	return iterations_until_rooted

plt.style.use('dark_background')
# plt.imshow(halley_method('x^3-1', 30, 1558, 1558, 30), extent=[-5, 5, -5, 5], cmap='inferno')
# plt.axis('on')
# plt.show()
# plt.close()

# for incrementation
for i in range(5, 6):
	t = i
	plt.imshow(halley_method('x^' + str(3+t/500) + '-1', 25, 2000, 1400, t), cmap='inferno')
	plt.axis('off')
	# plt.show()
	plt.savefig('halley{0:03d}.png'.format(i), bbox_inches='tight', dpi=420)
	plt.close()
