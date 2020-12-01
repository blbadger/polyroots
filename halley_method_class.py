# libraries
import numpy as np 
import matplotlib.pyplot as plt 
plt.style.use('dark_background')
from Calculate import Calculate


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
		f_now = Calculate(equation, z, differentiate=False).evaluate()
		f_prime_now = Calculate(equation, z, differentiate=True).evaluate()
		diff_string = Calculate(equation, z, differentiate=True).to_string()
		f_double_prime_now = Calculate(diff_string, z, differentiate=True).evaluate()
		z_array = z - (2*f_now * f_prime_now / (2*(f_prime_now)**2 - f_now * f_double_prime_now))

		# the boolean map is tested for rooted values
		found_root = (abs(z_array - previous_z_array) < 0.0000001) & not_already_at_root
		iterations_until_rooted[found_root] = i
		not_already_at_root = np.invert(found_root) & not_already_at_root

	return iterations_until_rooted

# plt.imshow(halley_method('x^3-1', 30, 1558, 1558, 30), extent=[-5, 5, -5, 5], cmap='inferno')
# plt.axis('on')
# plt.show()
# plt.close()

# for incrementation
for i in range(501):
	t = i
	plt.imshow(halley_method('x^' + str(13+t/500) + '-x-1', 25, 2000, 1400, t), cmap='inferno')
	plt.axis('off')
	# plt.show()
	plt.savefig('halley{0:03d}.png'.format(i), bbox_inches='tight', dpi=420)
	plt.close()
