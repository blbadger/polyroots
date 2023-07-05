#! python3

# import third party libraries
import numpy as np 
import matplotlib.pyplot as plt 
import torch
from torchvision.io import write_video
import time

from Calculate import Calculate 

# specify device
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print (device)
print(torch.version.cuda) 

def newton_raphson_map(equation, max_iterations, x_range, y_range, tensor_output=False):
	"""
	Generates a newton-raphson fractal.

	Args:
		equation: str, equation of interest
		max_iterations: int, number of iterations 
		x_range: int, number of real values per output
		y_range: int, number of imaginary values per output
		tensor_output: bool, if True then returns a torch.tensor object

	Returns:
		iterations_until_rooted: np.arr (2D) of iterations until a root is found
			at each point in y_range and x_range

	"""
	y, x = np.ogrid[1.1: -1.1: y_range*1j, -1.8: 1.8: x_range*1j]
	z_array = torch.tensor(x + y*1j).to(device)

	iterations_until_rooted = torch.tensor(max_iterations + np.zeros(z_array.shape)).to(device)

	 # create a boolean grid of all 'true'
	not_already_at_root = torch.ones(iterations_until_rooted.shape).to(device)

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

	if not tensor_output:
		iterations_until_rooted = iterations_until_rooted.cpu().numpy()

	return iterations_until_rooted


def newton_raphson_video(equation, max_iterations, z_array, iterations_until_rooted, not_already_at_root, tensor_output=False):
	"""
	Generates a newton-raphson fractal.

	Args:
		equation: str, equation of interest
		max_iterations: int, number of iterations 
		x_range: int, number of real values per output
		y_range: int, number of imaginary values per output
		tensor_output: bool, if True then returns a torch.tensor object

	Returns:
		iterations_until_rooted: np.arr (2D) of iterations until a root is found
			at each point in y_range and x_range

	"""
	nondiff = Calculate(equation, differentiate=False)
	diffed = Calculate(equation, differentiate=True)

	for i in range(max_iterations):
		# print (i)
		previous_z_array = z_array
		z = z_array
		f_now = nondiff.evaluate(z)
		f_prime_now = diffed.evaluate(z)
		z_array = z_array - f_now / f_prime_now

		# the boolean map is tested for rooted values
		found_root = torch.logical_and((torch.abs(z_array - previous_z_array) < 0.0000001), not_already_at_root)
		iterations_until_rooted[found_root] = i
		not_already_at_root = torch.logical_and(~found_root, not_already_at_root)

	if not tensor_output:
		iterations_until_rooted = iterations_until_rooted.cpu().numpy()

	return iterations_until_rooted

def increment_powers():
	n_frames = 100
	x_range, y_range = 2200, 1400
	max_iterations = 40
	y, x = np.ogrid[1.1: -1.1: y_range*1j, -1.8: 1.8: x_range*1j]
	z_array = torch.tensor(x + y*1j).to(device)
	iterations_until_rooted = torch.tensor(max_iterations + np.zeros(z_array.shape)).to(device)
	not_already_at_root = torch.ones(iterations_until_rooted.shape).to(device)
	video = torch.empty(n_frames, y_range, x_range)
	start_time = time.time()
	for t in range(100):
		plt.style.use('dark_background')
		frame = newton_raphson_video('x^(7.11+' + str(t*0.1/100000) + 'i)-x^(1+' + str(t*0.1/100000) + 'i)-1', 
			max_iterations, 
			z_array.clone(),
			iterations_until_rooted.clone(),
			not_already_at_root.clone(),
			tensor_output=False)
		# video[t] = frame
		plt.imshow(frame, cmap='inferno')
		plt.axis('off')
		plt.savefig('Newton.png', bbox_inches='tight', pad_inches=0, dpi=400)
		plt.close()
		print (time.time() - start_time)

	return

if __name__ == '__main__':
	plt.style.use('dark_background')
	# NB: no spaces are allowed in the input polynomial! Another example of a valid complex-valued input: '(0.2-3i)x^(2j)-x-(1-2i)'
	output = newton_raphson_map('x^5-x-1', 45, 2700, 1500, 0)
	plt.imshow(output, cmap='inferno')
	plt.axis('off')
	plt.savefig('Newton.png', bbox_inches='tight', pad_inches=0, dpi=400)
	plt.close()




