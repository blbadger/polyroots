#! python3 Newton_method_optimized.py

# libraries
import numpy as np 
import matplotlib.pyplot as plt 
from Calculate import Calculate # real value version of Calculate
from CalculateFaster import OptiCalculate
import numexpr as ne

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
		f_now = OptiCalculate(equation, z, differentiate=False).evaluate()
		f_prime_now = OptiCalculate(equation, z, differentiate=True).evaluate()
		z_array = ne.evaluate('z_array - f_now / f_prime_now')

		# the boolean map is tested for rooted values
		found_root = (abs(z_array - previous_z_array) < 0.0000001) & not_already_at_root
		iterations_until_rooted[found_root] = i
		not_already_at_root = np.invert(found_root) & not_already_at_root

	return iterations_until_rooted

plt.style.use('dark_background')

t = 0 # incrementor 
plt.imshow(newton_raphson_map('x^5.14-x-1', 40, 2000, 2000, t), extent=[-1, 1, -1, 1], cmap='inferno')
plt.axis('off')
# plt.savefig('Newton{0:03d}.png'.format(t), bbox_inches='tight', dpi=420)
plt.show()
plt.close()


# def period_2_finder(equation, max_iterations):
# 	"""
# 	Turns out there is no period-2 trajectory for the netwon map of x^5-x-1,
# 	even though there is a period-3 trajectory.  This does not violate the
# 	Li-Yorke or Sharkovskii's theorem, as the function in question is discontinuous
# 	at points where f'(x) = 0, and thus cannot be applied to either theorem.

# 	"""
# 	x_range, y_range = 5000, 5000
# 	y, x = np.ogrid[0.7823: 0.7824: y_range*1j, 0.3514: 0.3515: x_range*1j]
# 	z_array = x + y*1j
# 	initial_arr = z_array
# 	p2_set = set()
# 	z = z_array
# 	previous_z_array = z_array
# 	f_now = OptiCalculate(equation, z, differentiate=False).evaluate()
# 	f_prime_now = OptiCalculate(equation, z, differentiate=True).evaluate()
# 	z_array = ne.evaluate('z_array - f_now / f_prime_now')
# 	p2s_final = []

# 	for i in range(max_iterations):
# 		second_previous_z_array = previous_z_array
# 		previous_z_array = z_array
# 		z = z_array
# 		# compute next z
# 		f_now = OptiCalculate(equation, z, differentiate=False).evaluate()
# 		f_prime_now = OptiCalculate(equation, z, differentiate=True).evaluate()
# 		z_array = ne.evaluate('z_array - f_now / f_prime_now')
# 		found_p2 = (abs(z_array - second_previous_z_array) < 0.00001) & abs(z_array - previous_z_array > 0.00001)
# 		print (z_array[found_p2], second_previous_z_array[found_p2])

# 	for i, bool_list in enumerate(found_p2):
# 		for j, boolv in enumerate(bool_list):
# 			if boolv:
# 				p2s_final.append(initial_arr[i][j])
		
# 	return p2s_final

p2 = period_2_finder('x^5-x-1', 20)

print (p2[:10])