# Implementation of a wrapper for objects returned by Sage commands.
# Author: R. Mark Prosser

from sage.all import * # Sage libraries

# The following script imports relevant constants from C++ which may be used here
execfile("ParseC++Constants.py")

# Sage expression operators
import operator as _operator
relationOperators = {\
_operator.eq:'=', _operator.ne:'neq', _operator.lt:'<', _operator.le:'leq', _operator.gt:'>', _operator.ge:'geq' }

# Wrapper for the simplify operation
def DoSimplify(expr):
	if isinstance(expr, sage.symbolic.expression.Expression):
		try:
			return expr.simplify_full()
		except:
			pass
	return simplify(expr)

# Wrapper for the factor operation
def DoFactor(expr):
	decimal_classes = (sage.rings.real_mpfr.RealNumber, float)
	if isinstance(expr, decimal_classes):
		return factor(Rational(expr))
	return factor(expr)

# Wrapper for the solve operation
def DoSolve(expr,v):
	# todo: potentially improve results by setting explicit_solutions=True or to_poly_solve=True	
	return solve(expr,v)

#def DoSolve(expr):
	# e.g. solve([eq1,eq2,eq3,p==1],p,q,x,y)
#	return solve(expr)

# ExpressionTree is a Python class which represents the semantics of an expression tree
# and mimics the SCG expression tree type in C++.
# The ExpressionTree constructor takes as input a Sage object and initializes itself by 
# examining the Sage object and determining the expression tree structure.
class ExpressionTree(object):

	def __init__(self, obj=None, sid=None):
		
		self.sid = sid
		self.children = []

		# List of recognized classes for numeric types in sage
		numeric_classes = (sage.rings.integer.Integer, int, long, \
		sage.rings.real_mpfr.RealNumber, float, \
		#sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic, \
		sage.symbolic.constants.Constant, sage.symbolic.constants_c.E) # includes constants (e,pi)
		# Note: in sage, e is an instance of sage.symbolic.constants_c.E, which derives from sage.symbolic.expression.Expression
		# Therefore, isinstance(obj, numeric_classes) must be checked BEFORE isinstance(obj, sage.symbolic.expression.Expression)

		if obj==None: 
			return # nothing more to do

		# Complex number
		elif isinstance(obj, sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_quadratic):
			if obj.real()==0:
				if obj.imag()==1:
					self.__init__(repr(I)) # TODO: perhaps this isn't the best way to deal with i
				else:	
					self.__initMultExpr([obj.imag(),I]) # or ([obj.imag(),I],'*')
			elif obj.imag()==0:
				self.__init__(obj.real())
			else:
				self.__initAddExpr([obj.real(), obj.imag()*I])
		
		# Numbers
		elif isinstance(obj, numeric_classes):
			leaf = ExpressionTree()
			leaf.sid = TERMINAL_EXPR
			leaf.terms = repr(obj)
			self.sid = NUM_EXPR
			self.children.append(leaf)
			
		# Strings
		elif isinstance(obj, str):			
			if obj.isalpha() and sid!=TERMINAL_EXPR:
				leaf = ExpressionTree()
				leaf.sid = TERMINAL_EXPR
				leaf.terms = obj
				if sid==None: # allows the ID to be defined previously (e.g. NAME_EXPR)
					self.sid = VAR_EXPR # default
				self.children.append(leaf)
			else: # e.g. symbols +,-
				self.sid = TERMINAL_EXPR
				self.terms = obj

		# Expressions enclosed in parentheses
		elif sid==PAREN_EXPR:
			self.children.append(ExpressionTree('('))
			if isinstance(obj, ExpressionTree):
				self.children.append(obj)
			else:
				self.children.append(ExpressionTree(obj))
			self.children.append(ExpressionTree(')'))

		# General symbolic expressions
		elif isinstance(obj, sage.symbolic.expression.Expression):
		
			operator = obj.operator()
			operands = obj.operands()
			
			if operator is _operator.pow:
				self.__initExponentExpr(operands[0],operands[1])

			elif operator is _operator.neg:
				self.sid = NEG_EXPR
				self.children.append(ExpressionTree(operands[0]))
				
			elif operator is _operator.add:
				self.__initAddExpr(operands)#,'+')
				
			# Note: apparently the subtraction operation is never used in Sage's expression tree
			# It is represented by addition with a negative operand. 
			# So this code may never have a chance to execute.
			elif operator is _operator.sub:
				# Convert to addition: negate all operands after the first one
				for i in range(1,len(operands)):
					operands[i] = -operands[i]
				self.__initAddExpr(operands)#,'-')

			elif operator is _operator.mul:
				# "sort" the operands to ensure that coefficients appear on the left
				def comparator(a,b):
#					if not isinstance(b, numeric_classes) and not b.is_numeric():
					if (isinstance(b, sage.symbolic.expression.Expression) and \
					b.is_symbol() or (b.operator() is _operator.pow and b.op[0].is_symbol())):
						if isinstance(a, numeric_classes) or (isinstance(a, sage.symbolic.expression.Expression) and a.is_numeric()):
							return -1 # treat a < b (a is coefficient, b is variable)
					return 0
				operands.sort(comparator)
				self.__initMultExpr(operands)
				
			# Note: apparently the division operation is never used in Sage's expression tree
			# It is represented using multiplication and an exponent of -1.
			# So this code may never have a chance to execute.
			elif operator is _operator.div:
				self.__initFracExpr(operands[0],operands[1])

			elif operator in relationOperators.keys():
				self.sid = REL_EXPR
				self.children.append(ExpressionTree(operands[0]))
				self.children.append(ExpressionTree(operands[1]))
				self.children.append(ExpressionTree(relationOperators[operator],TERMINAL_EXPR)) 
				
			elif isinstance(operator, sage.symbolic.function.Function): # e.g. sin, exp, erf
				# Note: Euler's number can appear as a constant or a function
				if type(operator) == sage.functions.log.Function_exp:
					self.__initExponentExpr(e, operands[0])
				elif type(operator) == sage.symbolic.integration.integral.IndefiniteIntegral:
					self.__initIntegralExpr(operands[0],operands[1])
				elif type(operator) == sage.symbolic.integration.integral.DefiniteIntegral:
					self.__initIntegralExpr(operands[0],operands[1],operands[2],operands[3])
				else:
					self.sid = FN_EXPR
					self.children.append(ExpressionTree(repr(operator), NAME_EXPR))
					self.children.append(ExpressionTree(operands[0], PAREN_EXPR))  # operand (in parentheses)
				
			elif operator is None: # terminal
				if obj.is_numeric():
					self.__init__(obj.pyobject())
				elif obj.is_symbol():
					self.__init__(str(obj.variables()[0]))
				elif obj.is_constant(): # is_constant means "is a constant" (not "is constant")
					self.__init__(obj.pyobject())
					
			else: 
				raise TypeError("Unsupported sage operator for symbolic expression")

		elif sid==ADD_EXPR:
			# Note: we get here from __initAddExpr
			self.__initAddExpr(obj)
		
		# Lists of expressions
		elif isinstance(obj, sage.structure.sequence.Sequence_generic) or (type(obj)==list and sid==None):
			# (possibly nested) list of expressions
			if len(obj)==1:
				self.__init__(obj[0])
			elif len(obj)>1:
				self.sid = LIST_EXPR
				for o in obj:
					self.children.append(ExpressionTree(o))
	
		# Factorization expression (numeric)
		# Note: polynomial factorizations are represented as symbolic expressions (handled above)
		elif isinstance(obj, sage.structure.factorization.Factorization):
			# Note: use asterisk multiplication operator for factorizations
			f = list(obj) # get list of factors (ignoring unit factor)
			if obj.unit() == -1:
				self.sid = NEG_EXPR
				child = ExpressionTree()
				child.__initMultExpr(f,'*')
				child.__PostProcess()
				self.children.append(child)
			else:
				self.__initMultExpr(f,'*')

		# Tuples are assumed to represent (base, exponent) of a factor, in a numeric factorization
		elif type(obj)==tuple: 
			self.__initExponentExpr(obj[0],obj[1])

		# Fractions (rationals)
		elif isinstance(obj, sage.rings.rational.Rational):
			self.__initFracExpr(obj.numerator(), obj.denominator())

		# Matrices
		elif isinstance(obj, sage.matrix.matrix_space.matrix.Matrix):
			self.__initMatrixExpr(obj)

		# Expression representing the nullspace (right kernel)
		elif isinstance(obj, sage.modules.free_module.FreeModule_generic):
			if len(obj.basis()) > 0:
				# TODO: List of expressions. For now: first expression.
				self.__initMatrixExpr(obj.basis()[0].column()) # pass vector as 1-column matrix
			else:
				# zero vector is the only solution
				self.__initMatrixExpr(zero_vector(obj.degree()).column())
		
		# Apply post-processing to the expression for pretty-printing
		self.__PostProcess()
		
	# Initialize addition expression
	def __initAddExpr(self, operands):
		
		self.sid = ADD_EXPR
		self.children.append(ExpressionTree(operands[0]))
		opsymbol = '+'				
		if self.__IsNegated(operands[1]): 
			operands[1] = -operands[1]
			opsymbol = '-'
		if len(operands)==2:
			self.children.append(ExpressionTree(operands[1]))
		else:
			self.children.append(ExpressionTree(operands[1:], ADD_EXPR))
		self.children.append(ExpressionTree(opsymbol)) 

	# Initialize multiplication expression
	def __initMultExpr(self, factors, opsymbol=None):

		self.sid = MULT_EXPR
		if len(factors)==0:
			self.__init__(1)
		elif len(factors)==1:
			self.__init__(factors[0])
		else:
			left = ExpressionTree(factors.pop(0))
			right = ExpressionTree()
			right.__initMultExpr(factors,opsymbol) # recursion
			right.__PostProcess()
			self.children.append(left)
			self.children.append(right) 
			if opsymbol:
				self.children.append(ExpressionTree(opsymbol))				
			else:
				self.children.append(ExpressionTree(None, BLANK_EXPR))
	
	# Initialize fraction expression
	def __initFracExpr(self, numerator, denominator):

		if denominator==1:
			self.__init__(numerator)
		else:
			self.sid = FRAC_EXPR
			self.children.append(ExpressionTree(numerator))
			self.children.append(ExpressionTree(denominator))

	# Initialize exponent expression
	def __initExponentExpr(self, base, exponent):
		
		if exponent == -1:
			self.sid = FRAC_EXPR
			self.children.append(ExpressionTree(1))
			self.children.append(ExpressionTree(base))
		elif exponent == 1:
			self.__init__(base)
		elif exponent == 0.5:
			self.sid = ROOT_EXPR
			self.children.append(ExpressionTree(base))
		else:
			self.sid = SUP_EXPR
			self.children.append(ExpressionTree(base))
			self.children.append(ExpressionTree(exponent))		

	# Initialize integral expression
	def __initIntegralExpr(self, integrand, variable, a=None, b=None):
		
		self.sid = INTEGRAL_EXPR
		if (a is not None) and (b is not None):
			self.children.append(ExpressionTree(a))
			self.children.append(ExpressionTree(b))
		else:
			self.children.append(ExpressionTree(None, BLANK_EXPR))
			self.children.append(ExpressionTree(None, BLANK_EXPR))
		self.children.append(ExpressionTree(integrand))
		self.children.append(ExpressionTree(repr(variable)[0])) # only one-character variables
		
		

	# Initialize matrix expression
	def __initMatrixExpr(self, obj):
		
		d = obj.dimensions() 

		self.sid = PAREN_EXPR
		self.children.append(ExpressionTree('['))
		
		matrix = ExpressionTree()
		matrix.sid = MATRIX_EXPR
		
		numrows = ExpressionTree()
		numrows.sid = TERMINAL_EXPR
		numrows.terms = repr(d[0])
		matrix.children.append(numrows)

		numcols = ExpressionTree()
		numcols.sid = TERMINAL_EXPR
		numcols.terms = repr(d[1])
		matrix.children.append(numcols)
		
		rows = ExpressionTree()
		rows.sid = MATRIXROWS_EXPR
		
		for r in range(d[0]):
			row = ExpressionTree() 
			row.sid = MATRIXROW_EXPR
			
			for c in range(d[1]):
				row.children.append(ExpressionTree(obj[r,c]))
			
			rows.children.append(row)
			
		matrix.children.append(rows)
		self.children.append(matrix)
		self.children.append(ExpressionTree(']'))

	# Determine if a Sage object is negated (meaning it has a negative sign in front)
	@staticmethod
	def __IsNegated(obj):
		
		if isinstance(obj,sage.symbolic.expression.Expression): 
			if obj.is_negative():
				return true
			if obj.operator() is _operator.mul:
				for o in obj.operands():
					if o.is_numeric() and o.pyobject()==-1:
						return true
		
		return false
	
	# Determine if the first term of a non-parentheses, non-exponent expression tree is negated, 
	# AND optionally undo the negation (mutates the expression in that case).
	# If remove is set (optional), then remove the negative (if applicable)
	def __LeftSideIsNegated(self, remove=0):
		
		if self.sid==TERMINAL_EXPR and self.terms[0]=='-':
			if remove:
				self.terms = self.terms[1:]
			return true
		
		if self.sid==NUM_EXPR or self.sid==ADD_EXPR or self.sid==MULT_EXPR or self.sid==FRAC_EXPR:			
			return self.children[0].__LeftSideIsNegated(remove)
		
		if self.sid==NEG_EXPR:
			if remove:
				#self = self.children[0]
				# copy child to self (note: self-assignment won't work here)
				self.sid = self.children[0].sid
				self.children = self.children[0].children
				if self.sid==TERMINAL_EXPR:
					self.terms = self.children[0].terms
			return true
	
		return false
	
	# Determine if the first term (on the left) of a non-parentheses, non-exponent expression tree is numeric
	def __LeftSideIsNumeric(self):
		
		if self.sid==NUM_EXPR:
			return true
			
		if self.sid==ADD_EXPR or self.sid==MULT_EXPR or self.sid==FRAC_EXPR or self.sid==NEG_EXPR or self.sid==SUP_EXPR:
			return self.children[0].__LeftSideIsNumeric()
		
		return false
	
	# Post-processing: apply additional rules where necessary for presentation purposes (pretty-printing)
	def __PostProcess(self):
		
		if self.sid == MULT_EXPR:
			
			left = self.children[0]
			right = self.children[1]
			
			# add parentheses if necessary
			if left.sid == ADD_EXPR:
				self.children[0] = ExpressionTree(left, PAREN_EXPR)
				left = self.children[0]
			if right.sid == ADD_EXPR:
				self.children[1] = ExpressionTree(right, PAREN_EXPR)	
				right = self.children[1]

			# remove unit factors
			if left.sid == NUM_EXPR and left.children[0].terms == repr(-1): 
				self.sid = NEG_EXPR
				self.children = [right]				
			elif right.sid == NUM_EXPR and right.children[0].terms == repr(-1): 
				self.sid = NEG_EXPR
				self.children = [left]				
			elif left.sid == NUM_EXPR and left.children[0].terms == repr(1): 
				self = right # note: self-assignment works here because it follows from __init__
			elif right.sid == NUM_EXPR and right.children[0].terms == repr(1): 
				self = left	# note: self-assignment works here because it follows from __init__	
				
			# combine fractions
			elif left.sid==FRAC_EXPR and right.sid==FRAC_EXPR:
				
				# combine the product of fractions into a fraction of products
				left.sid = MULT_EXPR
				right.sid = MULT_EXPR
				self.sid = FRAC_EXPR
				
				# swap children appropriately
				temp = left.children[1]
				left.children[1] = right.children[0]
				right.children[0] = temp
				
				# multiplication operator (symbol)
				op = self.children.pop()
				left.children.append(op)
				right.children.append(op)
				
				# must post process the newly formed expression trees
				left.__PostProcess()
				right.__PostProcess()
				
			else:				
	
				# if both factors are numbers, use asterisk
				if left.sid == NUM_EXPR and right.sid == NUM_EXPR:
					self.children[2] = ExpressionTree('*')
					
				# if both factors are numbers, use asterisk
#				elif left.sid==NUM_EXPR and (right.sid==MULT_EXPR or right.sid==FRAC_EXPR or right.sid==SUP_EXPR):
				elif left.sid==NUM_EXPR and right.__LeftSideIsNumeric():
#					if right.children[0].sid == NUM_EXPR:
					self.children[2] = ExpressionTree('*')
		
				# if both factors are numbers, use asterisk
				elif right.sid==NUM_EXPR and (left.sid==MULT_EXPR or left.sid==FRAC_EXPR or left.sid==SUP_EXPR):
					if left.children[0].sid == NUM_EXPR:
						self.children[2] = ExpressionTree('*')
		
				# put numbers on the left
				elif right.sid == NUM_EXPR and left.sid != NUM_EXPR:
					# switch left and right
					temp = self.children[1]
					self.children[1] = self.children[0]
					self.children[0] = temp
					
				# put numbers on the left
				elif right.sid == MULT_EXPR and right.children[0].sid == NUM_EXPR and left.sid != NUM_EXPR:
					# switch expressions
					temp = right.children[0]
					right.children[0] = left
					self.children[0] = temp

					# must post process the newly formed expression tree
					right.__PostProcess()
					
			# if we're still a product at this point
			if self.sid == MULT_EXPR:
				left = self.children[0]
				right = self.children[1]
				
				# keep negative terms on left
				if right.__LeftSideIsNegated() and not left.__LeftSideIsNegated():
					# switch
					temp = left
					self.children[0] = right
					self.children[1] = temp			
		
		if self.sid == ADD_EXPR:
			left = self.children[0]
			right = self.children[1]
			op = self.children[2].terms
			
			# if there's a negative sign on the left of the right expression, absorb it
			if right.__LeftSideIsNegated(True):
				self.children[2].terms = '-' if op=='+' else '+' # flip the symbol	

			# if we're adding a negative to a positive, make it a subtraction
			# (this is subjective and unnecessary)
#			if op=='+' and left.__LeftSideIsNegated(True):
#				temp = left
#				self.children[0] = right
#				self.children[1] = temp
#				self.children[2].terms = '-'				
		
		if self.sid == SUP_EXPR:
			
			base = self.children[0]
			exponent = self.children[1]

			# add parentheses if necessary
			if base.sid == ADD_EXPR or base.sid == MULT_EXPR:
				self.children[0] = ExpressionTree(base, PAREN_EXPR)				
		
	# Get a string representation of the expression tree in prefix notation (like Scheme)
	def PrefixNotation(self):
		s = "(%s" % GetExprString(self.sid)
		if self.children == [] and hasattr(self,'terms'):
			s += " %s" % self.terms
		else:
			for child in self.children:
				if isinstance(child, ExpressionTree):
					s += " %s" % child.PrefixNotation()
				elif isinstance(child, str):
					s += " %s" % child
		s += ")"
		return s

	# Returns a string representation of the semantic structure of the expression tree
	def ToString(self):
		
		s = ""
		
		if self.sid == TERMINAL_EXPR:
			s = self.terms
		elif self.sid == BLANK_EXPR:
			s = ""
		elif self.sid == VAR_EXPR or self.sid == NUM_EXPR or self.sid == NAME_EXPR:
			s = self.children[0].ToString()
		elif self.sid == REL_EXPR or self.sid == ADD_EXPR:
			s = "{{{0}}} {1} {{{2}}}".format(self.children[0].ToString(), self.children[2].ToString(), self.children[1].ToString())
		elif self.sid == PAREN_EXPR:
			s = "{0}{1}{2}".format(self.children[0].ToString(), self.children[1].ToString(), self.children[2].ToString())
		elif self.sid == ROOT_EXPR:
			s = "sqrt(%s)" % self.children[0].ToString()
		elif self.sid == NEG_EXPR:
			s = "-{%s}" % self.children[0].ToString()
		elif self.sid == MULT_EXPR:
			s = "{{{0}}}{1}{{{2}}}".format(self.children[0].ToString(), self.children[2].ToString(), self.children[1].ToString())
		elif self.sid == FRAC_EXPR:
			s = "{{{0}}}/{{{1}}}".format(self.children[0].ToString(), self.children[1].ToString())
		elif self.sid == SUP_EXPR:
			s = "{{{0}}}^{{{1}}}".format(self.children[0].ToString(), self.children[1].ToString())
		elif self.sid == SUBSCR_EXPR:
			s = "{{{0}}}_{{{1}}}".format(self.children[0].ToString(), self.children[1].ToString())
		elif self.sid == FN_EXPR:
			s = "{{{0}}}{{{1}}}".format(self.children[0].ToString(), self.children[1].ToString())
		elif self.sid == INTEGRAL_EXPR:
			s = "integral({0},{1},{2},{3})".format(self.children[0].ToString(), self.children[1].ToString(), self.children[2].ToString(), self.children[3].ToString())
		elif self.sid == LIST_EXPR:
			s = "["
			for child in self.children:
				s += "{0},".format(child.ToString())
			s += "]\n"
		elif self.sid == MATRIX_EXPR:
			rows = self.children[2]
			for r in range(int(self.children[0].terms)):
				row = rows.children[r]
				s += "["
				for c in range(int(self.children[1].terms)):
					s += "%s, " % row.children[c].ToString()
				s += "]\n"

		return s
