!# python3

# libraries
import numpy as np 
import matplotlib.pyplot as plt 
plt.style.use('dark_background')
from Calculate import Calculate # custom parser, differentiator, and evaluator class

def newton_raphson_map(equation, max_iterations, x_range, y_range, t):
	print (equation)
	xl = -10/(2**(t/30)) + 0.41187
	xr = 10/(2**(t/30)) + 0.41187
	yl = 10/(2**(t/30))
	yr = -10/(2**(t/30)) 
	y, x = np.ogrid[yl: yr: y_range*1j, xl: xr: x_range*1j]
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

plt.imshow(newton_raphson_map('x^5-x-1', 30, 2000, 2000, 60), extent=[-10/(2**(t/30)) + 0.41187,\
								      10/(2**(t/30)) + 0.41187, \
								      10/(2**(t/30)), \
								      -10/(2**(t/30))], cmap='inferno')
plt.axis('off')
plt.show()
plt.close()

# Use for zoom 

# for i in range(300):
#   t = i
#   plt.imshow(newton_raphson_map('x^5-x-1', 30, 2000, 2000, t), extent=[-10/(2**(t/30)) + 0.41187, 10/(2**(t/30)) + 0.41187, 10/(2**(t/30)), -10/(2**(t/30))], cmap='inferno')
#   plt.axis('off')
#   # plt.show()
#   plt.savefig('Newton_Raphson{0:03d}.png'.format(i), bbox_inches='tight', dpi=420)
#   plt.close()
