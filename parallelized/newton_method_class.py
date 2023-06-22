#! python3

# import third party libraries
import numpy as np 
import matplotlib.pyplot as plt 
import torch

from Calculate import Calculate 

# specify device
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print (device)
print(torch.version.cuda) 

def newton_raphson_map(equation, max_iterations, x_range, y_range, t):
	"""
	Generates a newton-raphson fractal.

	Args:
		equation: str, equation of interest
		max_iterations: int, number of iterations 
		x_range: int, number of real values per output
		y_range: int, number of imaginary values per output
		t: int

	Returns:
		iterations_until_rooted: np.arr (2D) of iterations until a root is found
			at each point in y_range and x_range

	"""
	y, x = np.ogrid[1.1: -1.1: y_range*1j, -1.9: 1.9: x_range*1j]
	z_array = torch.tensor(x + y*1j).to(device)

	iterations_until_rooted = torch.tensor(max_iterations + np.zeros(z_array.shape)).to(device)

	 # create a boolean grid of all 'true'
	not_already_at_root = torch.tensor(iterations_until_rooted < 10000).to(device)

	nondiff = Calculate(equation, differentiate=False)
	diffed = Calculate(equation, differentiate=True)

	for i in range(max_iterations):
		previous_z_array = z_array
		z = z_array
		f_now = nondiff.evaluate(z)
		f_prime_now = diffed.evaluate(z)
		z_array = z_array - f_now / f_prime_now

		# the boolean map is tested for rooted values
		found_root = torch.logical_and((torch.abs(z_array - previous_z_array) < 0.0000001), not_already_at_root)
		iterations_until_rooted[found_root] = i
		not_already_at_root = torch.logical_and(~found_root, not_already_at_root)

	return iterations_until_rooted.cpu().numpy()


if __name__ == '__main__':

	plt.style.use('dark_background')
	output = newton_raphson_map('x^5-x-1', 45, 2700, 1500, 0)
	plt.imshow(output, cmap='inferno')
	plt.axis('off')
	plt.savefig('Newton.png', bbox_inches='tight', pad_inches=0, dpi=400)
	plt.close()

	# No spaces are allowed in the input polynomial! Another example of a valid complex-valued input:
	# '(0.2-3i)x^(2j)-x-(1-2i)'



