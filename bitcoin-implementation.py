import random
from typing import Union, Callable, Tuple


class EllipticPoint(object):
    def __init__(
        self,
        x: int,
        y: int = None,
        a: int = 486662,
        b: int = 1,
        c: int = 0,
        prime: int = ((2 ** 255) - 19),
    ):
        """
        For storing Elliptic Points corresponding to elements:

        y^2 = x^3 + ax^2 + bx + c

        Parameters a through d default to the ed25519 curve:

        y^2 = x^3 + 486662x^2 + x

        where a=486662, b = 1, & c = 0

        :param x: the integer for the x coordinate of the point
        :param y: the integer for the y cooirdinate. If none, it will calculate
            this based on the curve parameters provided (a - c)
        :param a: corresponds to the elliptic function (see above)
        :param b: corresponds to the elliptic function (see above)
        :param c: corresponds to the elliptic function (see above)
        :param prime: prime to use for modulo operations, defaults to ((2**255) - 19)
            from the ed25519 curve
        """
        self.x = x
        self._curve = EllipticPoint._curve_factory(a, b, c, prime)
        if y is not None:
            self.y = y
        else:
            self.y = self._curve(x)
        self.a, self.b, self.c = a, b, c
        self.prime = prime

    def to_tuple(self):
        return (self.x, self.y)

    def __repr__(self):
        return "<EllipticPoint: %r,%r>" % (self.to_tuple())

    def __eq__(self, anotherObj):
        if isinstance(anotherObj, EllipticPoint):
            return (
                self.x == anotherObj.x
                and self.y == anotherObj.y
                and self._same_curve(anotherObj)
            )
        elif isinstance(anotherObj, tuple):
            return self.x == anotherObj[0] and self.y == anotherObj[1]

    def __rmul__(self, scalar):
        return self.__mul__(scalar)

    def __lmul__(self, scalar):
        return self.__mul__(scalar)

    def __mul__(self, scalar):
        """
        Elliptic Multiply two points together. Both must have the same parameters.

        >>> G = EllipticPoint(x=5, y=1, b=2, prime=17)
        >>> 1*G == G
        True
        >>> 2*G == G+G
        True
        >>> 3*G == G+G+G
        True
        >>> (2*G)+G == G+(2*G)
        True
        >>> 4*G == G+G+G+G
        True
        >>> 19*G
        <EllipticPoint: inf,inf>
        >>> 20*G == G
        True
        >>> 21*G == 2*G
        True
        """
        if not type(scalar) in [int]:
            try:
                scalar = int(scalar)
                if scalar == 0:
                    raise
            except:
                raise ValueError(
                    "Can only multiply against scalar values -- value of type %r passed in: %r"
                    % (type(scalar), scalar)
                )
        assert scalar >= 0
        n = scalar
        R, Q = self.copy(), self.copy()
        R.x, R.y = float("inf"), float("inf")
        while n:
            if (n % 2) == 1:
                R += Q
            Q += Q
            n //= 2
        return R

    def __add__(self, anotherObj):
        """
        Given another point, it will add the two together and return the
        corresponding point in an EllipticPoint object. Both must have the
        same parameters

        >>> P = EllipticPoint(x=5, y=1, a=0, b=2, c=2, prime=17)
        >>> Q = EllipticPoint(x=5, y=1, a=0, b=2, c=2, prime=17)
        >>> P + Q
        <EllipticPoint: 6,3>

        >>> P = EllipticPoint(x=7, y=6, a=0, b=2, c=2, prime=17)
        >>> Q = EllipticPoint(x=5, y=1, a=0, b=2, c=2, prime=17)
        >>> P + Q
        <EllipticPoint: 7,11>
        >>> Q + P == P + Q
        True
        """
        if not self._same_curve(anotherObj):
            raise ValueError("Both points do not correspond to the same curve")

        if self == (float("inf"), float("inf")):
            return anotherObj
        elif anotherObj == (float("inf"), float("inf")):
            return self

        x_p, y_p = self.x, self.y
        x_q, y_q = anotherObj.x, anotherObj.y

        # point doubling
        if x_p == x_q and y_p == y_q:
            s_num = 3 * (x_p ** 2) + self.b
            s_denom = 2 * y_p
            s = s_num * EllipticPoint._multiplicative_inverse(s_denom, self.prime)
            x_r = (s ** 2) - (2 * x_p)
            y_r = (s * (x_p - x_r)) - y_p

        # adding vertical points -- point at infinity
        elif x_p == x_q and self.y != anotherObj.y:
            tmp = self.copy()
            tmp.x, tmp.y = float("inf"), float("inf")
            return tmp

        # standard addition
        else:
            s_num = y_p - y_q
            s_denom = x_p - x_q
            s = s_num * EllipticPoint._multiplicative_inverse(s_denom, self.prime)
            x_r = (s ** 2) - (x_p + x_q)
            y_r = (s * (x_p - x_r)) - y_p

        to_ret = self.copy()
        to_ret.x = x_r % self.prime
        to_ret.y = y_r % self.prime
        return to_ret

    def copy(self):
        """
        Return a perfect copy of this elliptic point.
        """
        return EllipticPoint(
            x=self.x, y=self.y, a=self.a, b=self.b, c=self.c, prime=self.prime
        )

    def _same_curve(self, anotherObj) -> bool:
        """
        Given a passed function _curve, this will return true if this curve
        has the same parameters as the curve represented within the object. Not
        meant for use outside of internal functions.
        """
        return (
            self.a == anotherObj.a
            and self.b == anotherObj.b
            and self.c == anotherObj.c
            and self.prime == anotherObj.prime
        )

    @staticmethod
    def _curve_factory(a: int, b: int, c: int, prime: int) -> Callable[[int], int]:
        """
        Returns a lambda function corresponding to the curve:

        y^2 = ( x^3 + ax^2 + bx + c ) % prime

        All parameters passed in correspond

        :returns: function corresponding to the curve
        """
        return lambda x: ((x ** 3) + (a * (x ** 2)) + (b * x) + c) % prime

    @staticmethod
    def _multiplicative_inverse(a: int, m: int) -> int:
        """
        Find the multiplicative inverse using the GCD & Extended
        Euclidean Algorith -- taken from
        https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Pseudocode

        :param x: number
        :param m: the modulo for the inverse

        :returns: the multiplicative inverse

        >>> EllipticPoint._multiplicative_inverse(3, 11)
        4

        >>> EllipticPoint._multiplicative_inverse(3, 17)
        6

        >>> EllipticPoint._multiplicative_inverse(2, 17)
        9
        """
        _r, r = a, m
        _s, s = 1, 0
        _t, t = 0, 1
        while r != 0:
            q = _r // r
            _r, r = r, _r - q * r
            _s, s = s, _s - q * s
            _t, t = t, _t - q * t

        return _s % m


if __name__ == "__main__":

    G = EllipticPoint(x=5, y=1, b=2, prime=17)
    # print(2 * G)
    # print((20 * G) + G)
    # print(G + (20 * G))
    # print((20 * G) + G)
    # print(21 * G)

    import doctest

    doctest.testmod()
