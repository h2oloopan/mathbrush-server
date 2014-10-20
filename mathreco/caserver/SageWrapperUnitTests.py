from SageWrapper import *
import unittest

# Unit testing for SageWrapper. Run "sage -python SageWrapperUnitTests.py"
#
# See http://docs.python.org/2/library/unittest.html

# See http://www.sagemath.org/doc/reference/sage/symbolic/expression.html#sage.symbolic.expression.Expression.solve
class TestSolve(unittest.TestCase):

    def test1(self):
        z = var('z')
        ex = z**5 - 1
        obj = solve(ex, z)
        self.assertEqual(ExpressionTree(obj[0:2]).ToString(), \
        '(LIST (REL (VAR (TERM z)) (FN (VAR (TERM exp)) (MULT (VAR (TERM pi)) (NUM (TERM 2/5*I)) (BLANK))) (TERM =)) (REL (VAR (TERM z)) (FN (VAR (TERM exp)) (MULT (VAR (TERM pi)) (NUM (TERM 4/5*I)) (BLANK))) (TERM =)))')

    def test2(self):
        x = var('x')
        ex = x**2>1
        obj = solve(ex, x)
        self.assertEqual(ExpressionTree(obj).ToString(), \
        '(LIST (REL (VAR (TERM x)) (NUM (TERM -1)) (TERM <)) (REL (VAR (TERM x)) (NUM (TERM 1)) (TERM >)))')
        
    # inconsistent results 
#    def test3(self):
#        x,y = var('x,y')
#        ex = ln(x)>ln(y)
#        obj = solve(ex, x)
#        self.assertEqual(ExpressionTree(obj).ToString(), \
#        '(LIST (REL (VAR (TERM y)) (VAR (TERM x)) (TERM <)) (REL (NUM (TERM 0)) (VAR (TERM y)) (TERM <)))')

class TestRelations(unittest.TestCase):

    def test1(self):
        x = var('x')
        obj = x == 1
        self.assertEqual(ExpressionTree(obj).ToString(), \
        '(REL (VAR (TERM x)) (NUM (TERM 1)) (TERM =))')

    def test1(self):
        x = var('x')
        obj = x != 1
        self.assertEqual(ExpressionTree(obj).ToString(), \
        '(REL (VAR (TERM x)) (NUM (TERM 1)) (TERM !=))')

    def test1(self):
        x = var('x')
        obj = x > 1
        self.assertEqual(ExpressionTree(obj).ToString(), \
        '(REL (VAR (TERM x)) (NUM (TERM 1)) (TERM >))')

    def test1(self):
        x = var('x')
        obj = x >= 1
        self.assertEqual(ExpressionTree(obj).ToString(), \
        '(REL (VAR (TERM x)) (NUM (TERM 1)) (TERM >=))')

    def test1(self):
        x = var('x')
        obj = x < 1
        self.assertEqual(ExpressionTree(obj).ToString(), \
        '(REL (VAR (TERM x)) (NUM (TERM 1)) (TERM <))')

    def test1(self):
        x = var('x')
        obj = x <= 1
        self.assertEqual(ExpressionTree(obj).ToString(), \
        '(REL (VAR (TERM x)) (NUM (TERM 1)) (TERM <=))')

class TestFactor(unittest.TestCase):

    def test1(self):
        n = -46
        obj = factor(n)
        self.assertEqual(ExpressionTree(obj).ToString(), \
        '(NEG (MULT (NUM (TERM 2)) (NUM (TERM 23)) (BLANK)))')

    def test2(self):
        x = var('x')
        ex = x**2-1
        obj = factor(ex)
        self.assertEqual(ExpressionTree(obj).ToString(), \
        '(MULT (ADD (VAR (TERM x)) (NUM (TERM 1)) (TERM -)) (ADD (VAR (TERM x)) (NUM (TERM 1)) (TERM +)) (BLANK))')
        
    def test3(self):
        n = 2012
        obj = factor(n)
        self.assertEqual(ExpressionTree(obj).ToString(), \
        '(MULT (SUP (NUM (TERM 2)) (NUM (TERM 2))) (NUM (TERM 503)) (BLANK))')

class TestNum(unittest.TestCase):

    def test1(self):
        obj = 12        
        self.assertEqual(ExpressionTree(obj).ToString(), \
        '(NUM (TERM 12))')
        
    def test2(self):
        obj = 12345678901234567890
        self.assertEqual(ExpressionTree(obj).ToString(), \
        '(NUM (TERM 12345678901234567890L))')
        
    def test3(self):
        obj = 2**0.5
        self.assertEqual(ExpressionTree(obj).ToString(), \
        '(NUM (TERM 1.4142135623730951))')
        
class TestOther(unittest.TestCase):

    def test_sqrt(self):
        obj = sqrt(2)
        self.assertEqual(ExpressionTree(obj).ToString(), \
        '(SUP (NUM (TERM 2)) (FRAC (NUM (TERM 1)) (NUM (TERM 2))))')
        
    def test_var(self):
        obj = var('x')        
        self.assertEqual(ExpressionTree(obj).ToString(), \
        '(VAR (TERM x))')
        
    def test_sin(self):
        x = var('x')
        obj = sin(x)
        self.assertEqual(ExpressionTree(obj).ToString(), \
        '(FN (VAR (TERM sin)) (VAR (TERM x)))')

if __name__ == '__main__':
    unittest.main()    