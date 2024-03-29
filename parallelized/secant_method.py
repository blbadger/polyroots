import numpy as np 
import matplotlib.pyplot as plt 
from Calculate import Calculate
import torch

plt.style.use('dark_background')
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print (f'Device: {device}')


def secant_method(equation, max_iterations, x_range, y_range, t):
	"""
	Returns an array of the number of iterations until a root is found
	using the Secant method in the complex plane.

	Args:
		equation: str, polynomial of interest
		max_iterations: int of iterations
		x_range: int, number of real values per output
		y_range: int, number of imaginary values per output
		t: int

	Returns:
		iterations_until_together: np.arr[int] (2D) 
		
	"""

	# top left to bottom right
	y, x = np.ogrid[1: -1: y_range*1j, -1: 1: x_range*1j]
	z_array = x + y*1j
	z_array = torch.tensor(z_array).to(device)
	iterations_until_rooted = torch.tensor(max_iterations + torch.zeros(z_array.shape)).to(device)

	 # create a boolean table of all 'true'
	not_already_at_root = torch.ones(iterations_until_rooted.shape).to(device)
	zeros = torch.zeros(z_array.shape).to(device)

	# set the initial guess to half the distance to the origin from the second guess
	z_0 = (z_array - zeros)/2 

	# initialize calculate object
	nondiff = Calculate(equation, differentiate=False)

	for i in range(max_iterations):
		previous_z_array = z_array
		z = z_array
		f_previous = nondiff.evaluate(z_0)
		f_now = nondiff.evaluate(z)
		z_array = z - f_now * (z - z_0)/(f_now - f_previous)

		# the boolean map is tested for rooted values
		found_root = torch.logical_and(abs(z_array - previous_z_array) < 1e-8, not_already_at_root)
		iterations_until_rooted[found_root] = i
		not_already_at_root = torch.logical_and(~found_root, not_already_at_root)

		# set previous array to current values
		z_0 = z 

	return iterations_until_rooted.cpu()


plt.imshow(secant_method('x^3-x-1', 50, 2000, 2000, 30), extent=[-1, 1, -1, 1], cmap='inferno')
plt.axis('off')
plt.show()
plt.close()
