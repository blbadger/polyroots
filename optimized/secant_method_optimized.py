import numpy as np 
import matplotlib.pyplot as plt 

from Calculate import Calculate
from CalculateFaster import OptiCalculate
import numexpr as ne



def secant_method(equation, max_iterations, x_range, y_range, t):
	print (equation)
	# top left to bottom right
	y, x = np.ogrid[1: -1: y_range*1j, -1: 1: x_range*1j]
	z_array = x + y*1j
	iterations_until_rooted = max_iterations + np.zeros(z_array.shape)

	 # create a boolean table of all 'true'
	not_already_at_root = iterations_until_rooted < 10000
	zeros = np.zeros(z_array.shape) 

	# setting the initial guess to half the distance to the origin from the second guess
	z_0 = ne.evaluate('(z_array - zeros)/2') 

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
		z_0 = z # set second previous value to nontransformed array

	return iterations_until_rooted


plt.style.use('dark_background')
plt.imshow(secant_method('x^3-x-x-1', 50, 2000, 2000, 30), extent=[-1, 1, -1, 1], cmap='inferno')
plt.axis('off')
plt.show()
plt.close()
