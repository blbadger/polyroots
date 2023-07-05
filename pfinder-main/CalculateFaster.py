#!python3 CalculateFaster.py

# import third party libraries
import numexpr as ne
from Calculate import Calculate
from RealCalculate import RealCalculate

class OptiCalculate(RealCalculate):
	'''
	Parses, differentiates, and evalutes polynomial expression for
	root finding algorithms.  Multithreading accross CPU cores
	and bytecode optimization via numexpr.  Inherits from 'Calculate'.
	'''

	def __init__(self, equation, differentiate):
		super(OptiCalculate, self).__init__(equation, differentiate)

	def evaluate(self, point):
		'''
		Evaluate expression using numexpr
		'''

		if self.diff:
			ls = self.differentiate()
		else:
			ls = self.parse()

		final_ls = []

		# convert to expression appropriate for numexpr evaluation
		for i, val in enumerate(ls):
			if val == '^':
				final_ls.append('**')

			elif val == 'x' and isinstance(ls[i-1], float):
				final_ls.append('*')
				final_ls.append('point')

			elif val not in ['start', 'end']:
				final_ls.append(val)

		if final_ls[-1] == '+':
			final_ls.pop()
		expression_string = ''.join([str(i) for i in final_ls])

		return ne.evaluate(expression_string)








