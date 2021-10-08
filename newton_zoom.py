# libraries
import numpy as np 
import matplotlib.pyplot as plt 
plt.style.use('dark_background')
from Calculate import Calculate

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

# check image
# plt.imshow(newton_raphson_map('x^5-x-1', 30, 1558, 1558), extent=[0.4, 0.42, -0.01, 0.01], cmap='inferno')
# plt.axis('off')
# plt.show()
# plt.close()

# zoom in by a factor of 2^10

for i in range(300):
	t = i
	plt.imshow(newton_raphson_map('x^5-x-1', 30, 1558, 1558, t), extent=[-10/(2**(t/30)) + 0.41187, 10/(2**(t/30)) + 0.41187, 10/(2**(t/30)), -10/(2**(t/30))], cmap='inferno')
	plt.axis('off')
	# plt.show()
	plt.savefig('Newton_Raphson{0:03d}.png'.format(i), bbox_inches='tight', dpi=420)
	plt.close()
