#!python3 Calculate.py

class Calculate:
	'''
	Parses, differentiates, and evalutes polynomial expression for
	root finding algorithms.  Single-threaded version.
	'''

	def __init__(self, equation, differentiate=False):
		self.equation = equation
		self.diff = differentiate

	def parse(self):
		"""
		Iterative parser to prepare a polynomial string for evaluation or 
		differentiation. Assumes that negative-value or complex exponents 
		are contained in parentheses.

		Args:
			None

		Returns:
			ls: arr[str or int or np.complex], parsed expression
		
		"""

		equation = self.equation
		digits = '0123456789e.i'
		characters_ls = [i for i in equation]

		# add start and end tags
		characters_ls = ['start'] + characters_ls
		characters_ls.append('end')

		# add constant a (ie in ax+b) if not supplied
		for i in range(len(characters_ls)-1):
			if (characters_ls[i] not in digits and characters_ls[i] != ')') and characters_ls[i+1] == 'x':
				characters_ls.insert(i+1, '1')

		# parse expression into list
		ls, i = [], 0
		while i in range(len(characters_ls)):

			# if parentheses exist, compile number contained and append to ls
			if characters_ls[i] in '()':
				number = ''
				j = 1
				while characters_ls[i+j] != ')':
					if characters_ls[i+j] in 'e.0123456789-+':
						number += characters_ls[i+j] # append char to number string
					else:
						number += 'j' # convert 'i' to 'j' for python
					j += 1

				# remove closing paren if necessary
				if characters_ls[i+j] == ')':
					j += 1
				ls.append(complex(number))
				i += j

			if characters_ls[i] in digits:
				number = ''
				j = 0
				while characters_ls[i+j] in digits:
					if characters_ls[i+j] in 'e.0123456789':
						number += characters_ls[i+j] # append char to number string
					else:
						number += 'j' # convert 'i' to 'j' for python
					j += 1

				ls.append(complex(number))
				i += j

			else:
				ls.append(characters_ls[i])
				i += 1

		return ls

	def to_string(self):
		"""
		Converts a list of components of a differentiated
		expression into a string.

		Args:
			None

		Returns:
			expression: str representing equation in ls

		"""

		if self.diff:
			ls = self.differentiate()
		else:
			ls = self.parse()

		expression = ''.join([str(i) for i in ls])

		return expression


	def differentiate(self):
		"""
		Finds the derivative of a given function 'equation'. 
		Accepts any polynomial with positive exponent values.

		Args:
			None (accesses self.equation via self.parse())

		Returns:
			final_ls: arr[str or float or np.complex]
		
		"""

		parsed_exp = self.parse()
		ls = parsed_exp

		# differentiate polynomial
		final_ls = []
		for i in range(len(ls)-1):

			if isinstance(ls[i], complex) and ls[i+1] == 'x' or (ls[i-1] == '^' and ls[i-2] == 'x'):
				final_ls.append(ls[i])

			if ls[i] == 'x':
				if ls[i+1] == '^':
					final_ls[-1] *= ls[i+2]

					# prevent divide by 0 error if exponent is 0
					if ls[i+2].real > 0 or ls[i+2].real < 0: 
						ls[i+2] -= 1

				final_ls.append(ls[i])

			if ls[i] == 'x' and ls[i+1] != '^':
				final_ls.append('^')
				final_ls.append(complex(0))

			if ls[i] in ['+', '-', '^']:
				final_ls.append(ls[i])


		# remove trailing symbols
		while True:
			# check if last list value is numeric, and if so stop look
			if isinstance(final_ls[-1], complex) or isinstance(final_ls[-1], float):
				break
			# otherwise remove
			final_ls.pop()

		return final_ls


	def evaluate(self, point):
		"""
		Computes the value of an equation for array `point`. Usually 
		a numpy ogrid, in which case the derivative is computed at each
		point on the grid. 

		Args:
			point: np.ogrid(np.complex), 2D ogrid of interest

		Returns:
			total: np.ogrid(complex) resulting from evaluating every
				   value in point

		"""

		if self.diff:
			final_ls = self.differentiate()
		else:
			final_ls = self.parse()

		if final_ls[0] != 'start':
			final_ls = ['start'] + final_ls

		if final_ls[-1] != 'end':
			final_ls.append('end')
		
		# change 'start' and 'end' to appropriate markers
		final_ls[0], final_ls[-1] = '+', '+' 


		# split parsed polynomial into blocks
		i = 0
		final_blocks = [[]]
		while i in range(len(final_ls)):
			ls = []
			j = 0
			while final_ls[i+j] not in ['+', '-']:
				ls.append(final_ls[i+j])
				j += 1

			if final_ls[i-1] == '-':
				if ls:
					ls[0] = -1 * ls[0]

			final_blocks.append(ls)
			i += j + 1

		# evaluate blocks sequentially
		total = 0
		for block in final_blocks:
			if block:
				if '^' not in block:
					if 'x' not in block:
						block += ['x', '^', 0]
					else:
						block += ['^', 1]

				start = block[0] * point ** block[-1]
				total += start

		return total



