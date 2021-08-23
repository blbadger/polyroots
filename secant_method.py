import numpy as np 
import matplotlib.pyplot as plt 
from Calculate import Calculate

plt.style.use('dark_background')


def secant_method(equation, max_iterations, x_range, y_range, t):
	print (equation)
	# top left to bottom right
	y, x = np.ogrid[1: -1: y_range*1j, -1: 1: x_range*1j]
	z_array = x + y*1j
	iterations_until_rooted = max_iterations + np.zeros(z_array.shape)

	 # create a boolean table of all 'true'
	not_already_at_root = iterations_until_rooted < 10000
	zeros = np.zeros(z_array.shape) 
	z_0 = (z_array - zeros)/2 # setting the initial guess to half the distance to the origin from the second guess, which is plotted

	for i in range(max_iterations):
		previous_z_array = z_array
		z = z_array
		f_previous = Calculate(equation, z_0, differentiate=False).evaluate()
		f_now = Calculate(equation, z, differentiate=False).evaluate()
		z_array = z - f_now * (z - z_0)/(f_now - f_previous)

		# the boolean map is tested for rooted values
		found_root = (abs(z_array - previous_z_array) < 0.0000001) & not_already_at_root
		iterations_until_rooted[found_root] = i
		not_already_at_root = np.invert(found_root) & not_already_at_root
		z_0 = z # set 

	return iterations_until_rooted

plt.imshow(secant_method('x^(7.3-i)-x-1', 50, 2000, 2000, 30), extent=[-1, 1, -1, 1], cmap='inferno')
plt.axis('off')
plt.show()
plt.close()
