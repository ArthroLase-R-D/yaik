import sympy as sp
from sympy.printing.codeprinter import CodePrinter
from sympy.printing.rust import RustCodePrinter
from sympy.printing.precedence import PRECEDENCE


class YaikRustCodePrinter(RustCodePrinter):
    """
    Rust printer where every number is an f64. The sympy RustCodePrinter emits integer
    literals (``2*x``) and relies on assumptions (is_real) that the yaik symbols do not carry,
    which produces code that does not type-check. Here all numeric literals are f64 and
    the only integer literals are powi exponents.
    """

    def __init__(self, use_safe_operator: bool = False, settings=None):
        super().__init__(settings or {})
        self._use_safe_operator = use_safe_operator

    # Numbers
    def _print_Integer(self, expr, _type=False):
        return '{}.0'.format(int(expr))

    def _print_Float(self, expr, _type=False):
        text = CodePrinter._print_Float(self, expr)
        if not any(ch in text for ch in '.eE') and text not in ('inf', 'nan'):
            text += '.0'
        return text

    def _print_Rational(self, expr):
        return '({}.0/{}.0)'.format(int(expr.p), int(expr.q))

    def _print_Half(self, expr):
        return self._print_Rational(expr)

    def _print_NegativeOne(self, expr):
        return '-1.0'

    def _print_One(self, expr):
        return '1.0'

    def _print_Zero(self, expr):
        return '0.0'

    def _print_Pi(self, expr, _type=False):
        return 'PI'

    def _print_Exp1(self, expr, _type=False):
        return 'E'

    def _print_Infinity(self, expr, _type=False):
        return 'f64::INFINITY'

    def _print_NegativeInfinity(self, expr, _type=False):
        return 'f64::NEG_INFINITY'

    def _print_NaN(self, expr, _type=False):
        return 'f64::NAN'

    # Arithmetic: use the plain C-like printing, without RustCodePrinter's casting logic
    def _print_Mul(self, expr):
        return CodePrinter._print_Mul(self, expr)

    def _print_Add(self, expr, order=None):
        return CodePrinter._print_Add(self, expr, order)

    def _method_receiver(self, expr) -> str:
        """
        Print an expression used as the receiver of a f64 method call, i.e. ``(<expr>).sqrt()``.
        """
        if expr.is_Symbol or expr.is_NumberSymbol or isinstance(expr, sp.Function):
            # Names, constants and calls (x.sin(), safe_asin(x)) bind tighter than a method call
            return self._print(expr)
        if expr.is_Rational and not expr.is_Integer:
            return self._print(expr)
        if expr.is_Number:
            # A bare literal has an ambiguous float type, so it needs a suffix
            return '({}_f64)'.format(self._print(expr))
        return '(' + self._print(expr) + ')'

    def _print_Pow(self, expr):
        base, exp = expr.base, expr.exp
        if exp == sp.S.Half:
            if self._use_safe_operator:
                return 'safe_sqrt({})'.format(self._print(base))
            return '{}.sqrt()'.format(self._method_receiver(base))
        if exp == -sp.S.Half:
            if self._use_safe_operator:
                return '(1.0/safe_sqrt({}))'.format(self._print(base))
            return '(1.0/{}.sqrt())'.format(self._method_receiver(base))
        if exp == -1:
            return '(1.0/{})'.format(self.parenthesize(base, PRECEDENCE['Mul'] + 1))
        if exp.is_Integer:
            return '{}.powi({})'.format(self._method_receiver(base), int(exp))
        if exp == sp.Rational(1, 3):
            return '{}.cbrt()'.format(self._method_receiver(base))
        return '{}.powf({})'.format(self._method_receiver(base), self._print(exp))

    # Functions
    def _print_Function(self, expr):
        name = expr.func.__name__
        if self._use_safe_operator and name in ('asin', 'acos'):
            return 'safe_{}({})'.format(name, self._print(expr.args[0]))
        unary_methods = {
            'sin': 'sin', 'cos': 'cos', 'tan': 'tan',
            'asin': 'asin', 'acos': 'acos', 'atan': 'atan',
            'sinh': 'sinh', 'cosh': 'cosh', 'tanh': 'tanh',
            'exp': 'exp', 'log': 'ln', 'Abs': 'abs', 'sign': 'signum',
            'floor': 'floor', 'ceiling': 'ceil',
        }
        if name in unary_methods and len(expr.args) == 1:
            return '{}.{}()'.format(self._method_receiver(expr.args[0]), unary_methods[name])
        if name == 'atan2':
            y, x = expr.args
            return '{}.atan2({})'.format(self._method_receiver(y), self._print(x))
        if name in ('Max', 'Min'):
            method = 'max' if name == 'Max' else 'min'
            text = self._print(expr.args[0])
            for arg in expr.args[1:]:
                text = '({}).{}({})'.format(text, method, self._print(arg))
            return text
        return super()._print_Function(expr)

    def _print_sign(self, expr):
        return self._print_Function(expr)

    # Logic
    def _print_Relational(self, expr):
        lhs = self._print(expr.lhs)
        rhs = self._print(expr.rhs)
        return '({} {} {})'.format(lhs, expr.rel_op, rhs)

    def _print_And(self, expr):
        return '(' + ' && '.join(self._print(arg) for arg in expr.args) + ')'

    def _print_Or(self, expr):
        return '(' + ' || '.join(self._print(arg) for arg in expr.args) + ')'

    def _print_Not(self, expr):
        return '!(' + self._print(expr.args[0]) + ')'

    def _print_BooleanTrue(self, expr, _type=False):
        return 'true'

    def _print_BooleanFalse(self, expr, _type=False):
        return 'false'


def rust_code(expr, use_safe_operator: bool = False) -> str:
    return YaikRustCodePrinter(use_safe_operator=use_safe_operator).doprint(expr)
