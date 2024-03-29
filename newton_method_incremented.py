# libraries
import numpy as np 
import matplotlib.pyplot as plt 
plt.style.use('dark_background')
from Calculate import Calculate


def newton_raphson_map(equation, max_iterations, x_range, y_range, t):
	print (equation)
	# top left to bottom right
	y, x = np.ogrid[0.7: -0.7: y_range*1j, -1: 1: x_range*1j]
	z_array = x + y*1j

	iterations_until_rooted = max_iterations + np.zeros(z_array.shape)

	 # create a boolean table of all 'true'
	not_already_at_root = iterations_until_rooted < 10000

	nondiff = Calculate(equation, differentiate=False)
	diff = Calculate(equation, differentiate=True)

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

	return iterations_until_rooted

# plt.imshow(newton_raphson_map('x^7-x+1', 40, 1558, 1558, 1), extent=[-5, 5, -5, 5], cmap='inferno')
# plt.axis('off')
# plt.show()
# plt.close()

for i in range(0, 10000):
	t = i
	plt.imshow(newton_raphson_map('x^' + str(7+i/90000) + '-x-1', 40, 2000, 1400, t), cmap='inferno')
	plt.axis('off')
	# plt.show()
	plt.savefig('Newton{0:03d}.png'.format(i), bbox_inches='tight', dpi=400)
	plt.close()
