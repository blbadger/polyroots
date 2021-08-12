#! python3

# import third party libraries
import numpy as np 
import matplotlib.pyplot as plt 

from Calculate import Calculate 


def newton_raphson_map(equation, max_iterations, x_range, y_range, t):
	print (equation)
	y, x = np.ogrid[1: -1: y_range*1j, -1: 1: x_range*1j]
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

t = 0
plt.style.use('dark_background')
plt.imshow(newton_raphson_map('x^5.14-x-1', 40, 2000, 2000, t), extent=[-1, 1, -1, 1], cmap='inferno')
plt.axis('off')
plt.savefig('Newton{0:03d}.png'.format(t), bbox_inches='tight', dpi=420)
# plt.show()
plt.close()






