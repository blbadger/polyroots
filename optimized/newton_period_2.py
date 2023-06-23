# libraries
import numpy as np 
import matplotlib.pyplot as plt 
from Calculate import Calculate # real-value version of Calculate
from CalculateFaster import OptiCalculate
import numexpr as ne

def period_2_finder(equation, max_iterations):
	"""
	Turns out there is no period-2 trajectory for the netwon map of x^5-x-1,
	even though there is a period-3 trajectory.  This does not violate the
	Li-Yorke or Sharkovskii's theorem, as the function in question is discontinuous
	at points where f'(x) = 0, and thus cannot be applied to either theorem.

	"""
	x_range, y_range = 5000, 5000
	y, x = np.ogrid[0.7823: 0.7824: y_range*1j, 0.3514: 0.3515: x_range*1j]
	z_array = x + y*1j
	initial_arr = z_array
	p2_set = set()
	z = z_array
	previous_z_array = z_array
	nondiff = OptiCalculate(equation, differentiate=False)
	diff = OptiCalculate(equation, differentiate=True)
	f_now = nondiff.evaluate(z)
	f_prime_now = diff.evaluate(z)
	z_array = ne.evaluate('z_array - f_now / f_prime_now')
	p2s_final = []

	for i in range(max_iterations):
		second_previous_z_array = previous_z_array
		previous_z_array = z_array
		z = z_array
		# compute next z
		f_now = nondiff.evaluate(z)
		f_prime_now = diff.evaluate(z)
		z_array = ne.evaluate('z_array - f_now / f_prime_now')
		found_p2 = (abs(z_array - second_previous_z_array) < 0.00001) & abs(z_array - previous_z_array > 0.00001)
		print (z_array[found_p2], second_previous_z_array[found_p2])

	for i, bool_list in enumerate(found_p2):
		for j, boolv in enumerate(bool_list):
			if boolv:
				p2s_final.append(initial_arr[i][j])
		
	return p2s_final


# p2 = period_2_finder('x^5-x-1', 20)