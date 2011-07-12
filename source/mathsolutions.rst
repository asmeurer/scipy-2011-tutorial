
.. include:: definitions.def

=====================
Solutions for Part II
=====================

Partial fraction decomposition
==============================

1.
  * `\frac{3 x + 5}{(2 x + 1)^2}`

Script style
------------

::

    f = (3*x + 5)/(2*x + 1)**2;f
    var('A:B')
    p1 = A/(2*x + 1)
    p2 = B/(2*x + 1)**2
    p1, p2
    h = sum((p1, p2));h
    hcombined = factor(together(h));hcombined
    eq = Eq(hcombined, f);eq
    coeffeq = Eq(numer(eq.lhs), numer(eq.rhs));eq
    sol = solve_undetermined_coeffs(coeffeq, [A, B], x);sol
    solution = h.subs(sol);solution
    apart(f, x)
    solution == apart(f, x)

Doctest Style
-------------

::

    >>> f = (3*x + 5)/(2*x + 1)**2;f
     3⋅x + 5  
    ──────────
             2
    (2⋅x + 1) 
    >>> var('A:B')
    (A, B)
    >>> p1 = A/(2*x + 1)
    >>> p2 = B/(2*x + 1)**2
    >>> p1, p2
    ⎛   A         B     ⎞
    ⎜───────, ──────────⎟
    ⎜2⋅x + 1           2⎟
    ⎝         (2⋅x + 1) ⎠
    >>> h = sum((p1, p2));h
       A          B     
    ─────── + ──────────
    2⋅x + 1            2
              (2⋅x + 1) 
    >>> hcombined = factor(together(h));hcombined
    2⋅A⋅x + A + B
    ─────────────
               2 
      (2⋅x + 1)  
    >>> eq = Eq(hcombined, f);eq
    2⋅A⋅x + A + B    3⋅x + 5  
    ───────────── = ──────────
               2             2
      (2⋅x + 1)     (2⋅x + 1) 
    >>> coeffeq = Eq(numer(eq.lhs), numer(eq.rhs));eq
    2⋅A⋅x + A + B = 3⋅x + 5
    >>> sol = solve_undetermined_coeffs(coeffeq, [A, B], x);sol
    {A: 3/2, B: 7/2}
    >>> solution = h.subs(sol);solution
         3             7      
    ─────────── + ────────────
    2⋅(2⋅x + 1)              2
                  2⋅(2⋅x + 1) 
    >>> apart(f, x)
         3             7      
    ─────────── + ────────────
    2⋅(2⋅x + 1)              2
                  2⋅(2⋅x + 1) 
    >>> solution == apart(f, x)
    True

  * `\frac{3 x + 5}{(u x + v)^2}`
  

Script Style
------------

::

    var('u v')
    f = (3*x + 5)/(u*x + v)**2;f
    var('A:B')
    p1 = A/(u*x + v)
    p2 = B/(u*x + v)**2
    p1, p2
    h = sum((p1, p2));h
    hcombined = factor(together(h));hcombined
    eq = Eq(hcombined, f);eq
    coeffeq = Eq(numer(eq.lhs), numer(eq.rhs));eq
    lhs, rhs = Poly(coeffeq.lhs, x), Poly(coeffeq.rhs, x)
    lhs
    rhs
    coeffeqs = [ Eq(lhs.nth(i), rhs.nth(i)) for i in xrange(2) ];coeffeqs
    solution = h.subs(sol);solution
    apart(f, x) # Note, you must include x!
    solution == apart(f, x)

Doctest Style
-------------

::

    >>> var('u v')
    (u, v)
    >>> f = (3*x + 5)/(u*x + v)**2;f
     3⋅x + 5  
    ──────────
             2
    (u⋅x + v) 
    >>> var('A:B')
    (A, B)
    >>> p1 = A/(u*x + v)
    >>> p2 = B/(u*x + v)**2
    >>> p1, p2
    ⎛   A         B     ⎞
    ⎜───────, ──────────⎟
    ⎜u⋅x + v           2⎟
    ⎝         (u⋅x + v) ⎠
    >>> h = sum((p1, p2));h
       A          B     
    ─────── + ──────────
    u⋅x + v            2
              (u⋅x + v) 
    >>> hcombined = factor(together(h));hcombined
    A⋅u⋅x + A⋅v + B
    ───────────────
                2  
       (u⋅x + v)   
    >>> eq = Eq(hcombined, f);eq
    A⋅u⋅x + A⋅v + B    3⋅x + 5  
    ─────────────── = ──────────
                2              2
       (u⋅x + v)      (u⋅x + v) 
    >>> coeffeq = Eq(numer(eq.lhs), numer(eq.rhs));eq
    A⋅u⋅x + A⋅v + B    3⋅x + 5  
    ─────────────── = ──────────
                2              2
       (u⋅x + v)      (u⋅x + v) 
    >>> lhs, rhs = Poly(coeffeq.lhs, x), Poly(coeffeq.rhs, x)
    >>> lhs
    Poly(A*u*x + A*v + B, x, domain='ZZ[u,v,A,B]')
    >>> rhs
    Poly(3*x + 5, x, domain='ZZ')
    >>> coeffeqs = [ Eq(lhs.nth(i), rhs.nth(i)) for i in xrange(2) ];coeffeqs
    [A⋅v + B = 5, A⋅u = 3]
    >>> solution = h.subs(sol);solution
     5⋅u - 3⋅v          3     
    ──────────── + ───────────
               2   u⋅(u⋅x + v)
    u⋅(u⋅x + v)               
    >>> apart(f, x) # Note, you must include x!
     5⋅u - 3⋅v          3     
    ──────────── + ───────────
               2   u⋅(u⋅x + v)
    u⋅(u⋅x + v)               
    >>> solution == apart(f, x)
    True

  * `\frac{(3 x + 5)^2}{(2 x + 1)^2}`

Script Style
------------

::

    f = (3*x + 5)**2/(2*x + 1)**2;f
    # Here, deg(p) == deg(q), so we have to use div()
    p, q = numer(f), denom(f)
    d, r = div(p, q) # p = d*q + r, i.e., f = p/q = d + r/q
    d, r
    solution = d # This is the polynomial part of the solution
    var('A:B')
    p1 = A/(2*x + 1)
    p2 = B/(2*x + 1)**2
    p1, p2
    h = sum((p1, p2));h
    hcombined = factor(together(h));hcombined
    eq = Eq(hcombined, r/q);eq
    coeffeq = Eq(numer(eq.lhs), r);eq # No need to call numer() here; we know it's r
    coeffeqs = collect(coeffeq.lhs - coeffeq.rhs, x, evaluate=False);sol
    # Note, since we're just subtracting coeffeq.lhs from coeffeq.rhs, we could have
    # done it without using Eq() in the first place.
    sol = solve(coeffeqs.values())
    solution += h.subs(sol);solution # Note the +=
    apart(f, x)
    solution == apart(f, x)

Doctest Style
-------------

::

    >>> f = (3*x + 5)**2/(2*x + 1)**2;f
             2
    (3⋅x + 5) 
    ──────────
             2
    (2⋅x + 1) 
    >>> # Here, deg(p) == deg(q), so we have to use div()
    >>> p, q = numer(f), denom(f)
    >>> d, r = div(p, q) # p = d*q + r, i.e., f = p/q = d + r/q
    >>> d, r
    (9/4, 21⋅x + 91/4)
    >>> solution = d # This is the polynomial part of the solution
    >>> var('A:B')
    (A, B)
    >>> p1 = A/(2*x + 1)
    >>> p2 = B/(2*x + 1)**2
    >>> p1, p2
    ⎛   A         B     ⎞
    ⎜───────, ──────────⎟
    ⎜2⋅x + 1           2⎟
    ⎝         (2⋅x + 1) ⎠
    >>> h = sum((p1, p2));h
       A          B     
    ─────── + ──────────
    2⋅x + 1            2
              (2⋅x + 1) 
    >>> hcombined = factor(together(h));hcombined
    2⋅A⋅x + A + B
    ─────────────
               2 
      (2⋅x + 1)  
    >>> eq = Eq(hcombined, r/q);eq
    2⋅A⋅x + A + B   21⋅x + 91/4
    ───────────── = ───────────
               2              2
      (2⋅x + 1)      (2⋅x + 1) 
    >>> coeffeq = Eq(numer(eq.lhs), r);eq # No need to call numer() here; we know it's r
    2⋅A⋅x + A + B   21⋅x + 91/4
    ───────────── = ───────────
               2              2
      (2⋅x + 1)      (2⋅x + 1) 
    >>> coeffeqs = collect(coeffeq.lhs - coeffeq.rhs, x, evaluate=False);sol
    {A: 21/2, B: 49/4}
    >>> # Note, since we're just subtracting coeffeq.lhs from coeffeq.rhs, we could have
    >>> # done it without using Eq() in the first place.
    >>> sol = solve(coeffeqs.values())
    >>> solution += h.subs(sol);solution # Note the +=
    9        21            49     
    ─ + ─────────── + ────────────
    4   2⋅(2⋅x + 1)              2
                      4⋅(2⋅x + 1) 
    >>> apart(f, x)
    9        21            49     
    ─ + ─────────── + ────────────
    4   2⋅(2⋅x + 1)              2
                      4⋅(2⋅x + 1) 
    >>> solution == apart(f, x)
    True

2. Even though you can use :func:`Expr.coeff` to get the coefficient of
`x^n` in a polynomial for `n > 0`, it will not work for `n = 0`.  This
is because :func:`Expr.coeff` considers any term with no numerical
coefficient to be the coefficient of 1.  From the docstring::

    You can select terms with no rational coefficient:
    >>> (x+2*y).coeff(1)
    x
    >>> (3+2*x+4*x**2).coeff(1)

There is an issue open (`2558
<http://code.google.com/p/sympy/issues/detail?id=2558>`_) that suggests
a way to fix this, but it has not yet been implemented.

Deriving Trigonometric Identities
====================================

1.

Script Style
------------

::

    var('a,b')
    f = sin(a + b);f
    g = f.series(a, 0, 10)
    g
    g = collect(g, [sin(b), cos(b)]).removeO()
    g
    sina = sin(a).series(a, 0, 10).removeO()
    sina
    h = g.subs(sina, sin(a))
    h
    cosa = cos(a).series(a, 0, 10).removeO()
    cosa
    h = h.subs(cosa, cos(a))
    h
    eq = Eq(f, h)
    eq
    eq.lhs.expand(trig=True) == eq.rhs

Doctest Style
-------------

::

    >>> var('a,b')
    (a, b)
    >>> f = sin(a + b);f
    sin(a + b)
    >>> g = f.series(a, 0, 10)
    >>> g
                         2           3           4           5           6           7           8           9
                        a ⋅sin(b)   a ⋅cos(b)   a ⋅sin(b)   a ⋅cos(b)   a ⋅sin(b)   a ⋅cos(b)   a ⋅sin(b)   a ⋅cos(b)
    sin(b) + a⋅cos(b) - ───────── - ───────── + ───────── + ───────── - ───────── - ───────── + ───────── + ───────── + O(a**10)
                            2           6           24         120         720         5040       40320       362880
    >>> g = collect(g, [sin(b), cos(b)]).removeO()
    >>> g
    ⎛   8      6    4    2    ⎞          ⎛   9       7      5    3    ⎞
    ⎜  a      a    a    a     ⎟          ⎜  a       a      a    a     ⎟
    ⎜───── - ─── + ── - ── + 1⎟⋅sin(b) + ⎜────── - ──── + ─── - ── + a⎟⋅cos(b)
    ⎝40320   720   24   2     ⎠          ⎝362880   5040   120   6     ⎠
    >>> sina = sin(a).series(a, 0, 10).removeO()
    >>> sina
       9       7      5    3
      a       a      a    a
    ────── - ──── + ─── - ── + a
    362880   5040   120   6
    >>> h = g.subs(sina, sin(a))
    >>> h
    ⎛   8      6    4    2    ⎞
    ⎜  a      a    a    a     ⎟
    ⎜───── - ─── + ── - ── + 1⎟⋅sin(b) + sin(a)⋅cos(b)
    ⎝40320   720   24   2     ⎠
    >>> cosa = cos(a).series(a, 0, 10).removeO()
    >>> cosa
       8      6    4    2
      a      a    a    a
    ───── - ─── + ── - ── + 1
    40320   720   24   2
    >>> h = h.subs(cosa, cos(a))
    >>> h
    sin(a)⋅cos(b) + sin(b)⋅cos(a)
    >>> eq = Eq(f, h)
    >>> eq
    sin(a + b) = sin(a)⋅cos(b) + sin(b)⋅cos(a)
    >>> eq.lhs.expand(trig=True) == eq.rhs
    True

2.

Script Style
------------

::

    var('a,b')
    f = cos(a + b);f
    g = f.series(a, 0, 10)
    g
    g = collect(g, [sin(b), cos(b)]).removeO()
    g
    sina = sin(a).series(a, 0, 10).removeO()
    sina
    # Note that subs will not work directly with sina, because -sina appears
    # in g.  We get around it by substituting -sina.
    h = g.subs(-sina, -sin(a))
    h
    cosa = cos(a).series(a, 0, 10).removeO()
    cosa
    h = h.subs(cosa, cos(a))
    h
    eq = Eq(f, h)
    eq
    eq.lhs.expand(trig=True) == eq.rhs

Doctest Style
-------------

::

    >>> var('a,b')
    (a, b)
    >>> f = cos(a + b);f
    cos(a + b)
    >>> g = f.series(a, 0, 10)
    >>> g
                         2           3           4           5           6           7           8           9
                        a ⋅cos(b)   a ⋅sin(b)   a ⋅cos(b)   a ⋅sin(b)   a ⋅cos(b)   a ⋅sin(b)   a ⋅cos(b)   a ⋅sin(b)
    cos(b) - a⋅sin(b) - ───────── + ───────── + ───────── - ───────── - ───────── + ───────── + ───────── - ───────── + O(a**10)
                            2           6           24         120         720         5040       40320       362880
    >>> g = collect(g, [sin(b), cos(b)]).removeO()
    >>> g
    ⎛   8      6    4    2    ⎞          ⎛     9       7      5    3    ⎞
    ⎜  a      a    a    a     ⎟          ⎜    a       a      a    a     ⎟
    ⎜───── - ─── + ── - ── + 1⎟⋅cos(b) + ⎜- ────── + ──── - ─── + ── - a⎟⋅sin(b)
    ⎝40320   720   24   2     ⎠          ⎝  362880   5040   120   6     ⎠
    >>> sina = sin(a).series(a, 0, 10).removeO()
    >>> sina
       9       7      5    3
      a       a      a    a
    ────── - ──── + ─── - ── + a
    362880   5040   120   6
    >>> # Note that subs will not work directly with sina, because -sina appears
    >>> # in g.  We get around it by substituting -sina.
    >>> h = g.subs(-sina, -sin(a))
    >>> h
    ⎛   8      6    4    2    ⎞
    ⎜  a      a    a    a     ⎟
    ⎜───── - ─── + ── - ── + 1⎟⋅cos(b) - sin(a)⋅sin(b)
    ⎝40320   720   24   2     ⎠
    >>> cosa = cos(a).series(a, 0, 10).removeO()
    >>> cosa
       8      6    4    2
      a      a    a    a
    ───── - ─── + ── - ── + 1
    40320   720   24   2
    >>> h = h.subs(cosa, cos(a))
    >>> h
    -sin(a)⋅sin(b) + cos(a)⋅cos(b)
    >>> eq = Eq(f, h)
    >>> eq
    cos(a + b) = -sin(a)⋅sin(b) + cos(a)⋅cos(b)
    >>> eq.lhs.expand(trig=True) == eq.rhs
    True

Not only Symbolics: Numerical Computing
=======================================

1.

Script Style
-------------

::

    f = x**(1 - log(log(log(log(1/x)))))
    f.subs(x, pi).evalf()

Doctest Style
-------------

::

    >>> f = x**(1 - log(log(log(log(1/x)))))
    >>> f.subs(x, pi).evalf()
    0.730290019485505 - 1.30803618190213⋅ⅈ

2.

Script Style
------------

::

    E_million_digits = str(E.evalf(n=1000000))
    pi_million_digits = str(pi.evalf(n=1000000))
    # We can use the .find() method of str.
    # First, let's see what we should offset the value
    # (it's not immediately obvious because of 0 indexing
    # and because of the '.' in the string.
    E_million_digits[:5]
    E_million_digits.find('71')
    # So, counting the 2 as the first digit, we see that
    # there is no need for any offset!
    E999999 = E_million_digits.find('999999')
    E999999
    "...%s..." % E_million_digits[E999999:E999999+10]
    # We see that there are actually 8 consecutive 9's here
    pi789 = pi_million_digits.find('789')
    pi789
    "...%s..." % pi_million_digits[pi789:pi789+10]
    # Now, double check the values from Dinosaur Comics.
    E789 = E_million_digits.find('789')
    E789
    "...%s..." % E_million_digits[E789:E789+10]
    pi999999 = pi_million_digits.find('999999')
    pi999999
    # Hey, this one was wrong in the webcomic (T-Rex said "762nd")!
    "...%s..." % pi_million_digits[pi999999:pi999999+10]

Doctest Style
-------------

::

    >>> E_million_digits = str(E.evalf(n=1000000))
    >>> pi_million_digits = str(pi.evalf(n=1000000))
    >>> # We can use the .find() method of str.
    >>> # First, let's see what we should offset the value
    >>> # (it's not immediately obvious because of 0 indexing
    >>> # and because of the '.' in the string.
    >>> E_million_digits[:5]
    2.718
    >>> E_million_digits.find('71')
    2
    >>> # So, counting the 2 as the first digit, we see that
    >>> # there is no need for any offset!
    >>> E999999 = E_million_digits.find('999999')
    >>> E999999
    384341
    >>> "...%s..." % E_million_digits[E999999:E999999+10]
    ...9999999951...
    >>> # We see that there are actually 8 consecutive 9's here
    >>> pi789 = pi_million_digits.find('789')
    >>> pi789
    352
    >>> "...%s..." % pi_million_digits[pi789:pi789+10]
    ...7892590360...
    >>> # Now, double check the values from Dinosaur Comics.
    >>> E789 = E_million_digits.find('789')
    >>> E789
    2501
    >>> "...%s..." % E_million_digits[E789:E789+10]
    ...7893652605...
    >>> pi999999 = pi_million_digits.find('999999')
    >>> pi999999
    763
    >>> # Hey, this one was wrong in the webcomic (T-Rex said "762nd")!
    >>> "...%s..." % pi_million_digits[pi999999:pi999999+10]
    ...9999998372...

3. According to Gruntz's thesis (pg. 7), `n`, the number of digits of
   precision in base 10 required to numerically evaluate the function at
   the point `x`, must satisfy `n \log(10) \geq e^x + \log(x)`,
   otherwise the floating point approximation of the difference of the
   error functions will be 0.  In other words, the number of *digits*
   required grows exponentially with `x` (!)  We can see that for very
   small (compared to a numerical conception of `\infty`) values of `x`,
   the number of digits is more than can be computed in a reasonable
   amount of time (mpmath can compute with arbitrary precision, but it
   of course requires a large amount of time and memory to compute a
   large number)::

    >>> [(x, int(((exp(x) + log(x))/log(10)).evalf())) for x in xrange(1, 15)]
    [(1, 1), (2, 3), (3, 9), (4, 24), (5, 65), (6, 175), (7, 477), (8, 1295), (9, 3520), (10, 9566), (11, 26004), (12, 70684), (13, 192138), (14, 522285)]

A little experimentation will reveal that around `x = 11` or so, the
required value of `n` to pass to :func:`Expr.evalf` becomes too large to
reasonably compute.  On the other hand, SymPy's :func:`limit` function
can compute the limit symbolically exactly and very quickly.

Script Style
------------

::

    limit((erf(x - exp(-exp(x))) - erf(x))*exp(exp(x))*exp(x**2), x, oo)

Doctest Style
-------------

::

    >>> limit((erf(x - exp(-exp(x))) - erf(x))*exp(exp(x))*exp(x**2), x, oo)
      -2
    ─────
      ⎽⎽⎽
    ╲╱ π


Summing Roots of Polynomials
============================

1.

* `f = z^5 + z + a` and `g = \frac{1}{z + 1}`

Script Style
------------

::

    var('a')
    f = z**5 + z + a
    f
    g = 1/(z + 1)
    g
    R = var('r:5')
    G = together(sum( [ g.subs(z, r) for r in R] ))
    G
    V = viete(f, R, z)
    for lhs, rhs in V:
        pprint(Eq(lhs, rhs))

    p = expand(numer(G))
    q = expand(denom(G))
    (P, Q), mapping = symmetrize((p, q), R, formal=True)
    P
    Q
    for s, poly in mapping:
        pprint(Eq(s, poly))

    subslist = [ (s, c) for (s, _), (_, c) in zip(mapping, V) ]
    subslist
    P[0].subs(subslist)/Q[0].subs(subslist)
    # Note, we must give a variable to RootSum, as f has more than one Symbol
    cancel(P[0].subs(subslist)/Q[0].subs(subslist) - RootSum(f, Lambda(z, g), z)) == 0

Doctest Style
-------------

::

    >>> var('a')
    a
    >>> f = z**5 + z + a
    >>> f
         5
    a + z  + z
    >>> g = 1/(z + 1)
    >>> g
      1
    ─────
    z + 1
    >>> R = var('r:5')
    >>> G = together(sum( [ g.subs(z, r) for r in R] ))
    >>> G
    (r₀ + 1)⋅(r₁ + 1)⋅(r₂ + 1)⋅(r₃ + 1) + (r₀ + 1)⋅(r₁ + 1)⋅(r₂ + 1)⋅(r₄ + 1) + (r₀ + 1)⋅(r₁ + 1)⋅(r₃ + 1)⋅(r₄ + 1) + (r₀ + 1)⋅(r₂ + 1)⋅(r₃ + 1)⋅(r₄ + 1) + (r₁ + 1)⋅(r₂ + 1)⋅(r₃ + 1)⋅(r₄ + 1)
    ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                                                                            (r₀ + 1)⋅(r₁ + 1)⋅(r₂ + 1)⋅(r₃ + 1)⋅(r₄ + 1)
    >>> V = viete(f, R, z)
    >>> for lhs, rhs in V:
    ...     pprint(Eq(lhs, rhs))
    ...
    r₀ + r₁ + r₂ + r₃ + r₄ = 0
    r₀⋅r₁ + r₀⋅r₂ + r₀⋅r₃ + r₀⋅r₄ + r₁⋅r₂ + r₁⋅r₃ + r₁⋅r₄ + r₂⋅r₃ + r₂⋅r₄ + r₃⋅r₄ = 0
    r₀⋅r₁⋅r₂ + r₀⋅r₁⋅r₃ + r₀⋅r₁⋅r₄ + r₀⋅r₂⋅r₃ + r₀⋅r₂⋅r₄ + r₀⋅r₃⋅r₄ + r₁⋅r₂⋅r₃ + r₁⋅r₂⋅r₄ + r₁⋅r₃⋅r₄ + r₂⋅r₃⋅r₄ = 0
    r₀⋅r₁⋅r₂⋅r₃ + r₀⋅r₁⋅r₂⋅r₄ + r₀⋅r₁⋅r₃⋅r₄ + r₀⋅r₂⋅r₃⋅r₄ + r₁⋅r₂⋅r₃⋅r₄ = 1
    r₀⋅r₁⋅r₂⋅r₃⋅r₄ = -a
    >>> p = expand(numer(G))
    >>> q = expand(denom(G))
    >>> (P, Q), mapping = symmetrize((p, q), R, formal=True)
    >>> P
    (4⋅s₁ + 3⋅s₂ + 2⋅s₃ + s₄ + 5, 0)
    >>> Q
    (s₁ + s₂ + s₃ + s₄ + s₅ + 1, 0)
    >>> for s, poly in mapping:
    ...     pprint(Eq(s, poly))
    ...
    s₁ = r₀ + r₁ + r₂ + r₃ + r₄
    s₂ = r₀⋅r₁ + r₀⋅r₂ + r₀⋅r₃ + r₀⋅r₄ + r₁⋅r₂ + r₁⋅r₃ + r₁⋅r₄ + r₂⋅r₃ + r₂⋅r₄ + r₃⋅r₄
    s₃ = r₀⋅r₁⋅r₂ + r₀⋅r₁⋅r₃ + r₀⋅r₁⋅r₄ + r₀⋅r₂⋅r₃ + r₀⋅r₂⋅r₄ + r₀⋅r₃⋅r₄ + r₁⋅r₂⋅r₃ + r₁⋅r₂⋅r₄ + r₁⋅r₃⋅r₄ + r₂⋅r₃⋅r₄
    s₄ = r₀⋅r₁⋅r₂⋅r₃ + r₀⋅r₁⋅r₂⋅r₄ + r₀⋅r₁⋅r₃⋅r₄ + r₀⋅r₂⋅r₃⋅r₄ + r₁⋅r₂⋅r₃⋅r₄
    s₅ = r₀⋅r₁⋅r₂⋅r₃⋅r₄
    >>> subslist = [ (s, c) for (s, _), (_, c) in zip(mapping, V) ]
    >>> subslist
    [(s₁, 0), (s₂, 0), (s₃, 0), (s₄, 1), (s₅, -a)]
    >>> P[0].subs(subslist)/Q[0].subs(subslist)
      6
    ──────
    -a + 2
    >>> # Note, we must give a variable to RootSum, as f has more than one Symbol
    >>> cancel(P[0].subs(subslist)/Q[0].subs(subslist) - RootSum(f, Lambda(z, g), z)) == 0
    True


* `f = z^5 + z + a` and `g = \frac{1}{z + b}`

Script Style
------------

::

    var('a b')
    f = z**5 + z + a
    f
    g = 1/(z + b)
    g
    R = var('r:5')
    G = together(sum( [ g.subs(z, r) for r in R] ))
    G
    V = viete(f, R, z)
    for lhs, rhs in V:
        pprint(Eq(lhs, rhs))

    p = expand(numer(G))
    q = expand(denom(G))
    (P, Q), mapping = symmetrize((p, q), R, formal=True)
    P
    Q
    for s, poly in mapping:
        pprint(Eq(s, poly))

    subslist = [ (s, c) for (s, _), (_, c) in zip(mapping, V) ]
    subslist
    P[0].subs(subslist)/Q[0].subs(subslist)
    # Note, we must give a variable to RootSum, as f has more than one Symbol
    cancel(P[0].subs(subslist)/Q[0].subs(subslist) - RootSum(f, Lambda(z, g), z)) == 0

Doctest Style
-------------

::

    >>> var('a b')
    (a, b)
    >>> f = z**5 + z + a
    >>> f
         5
    a + z  + z
    >>> g = 1/(z + b)
    >>> g
      1
    ─────
    b + z
    >>> R = var('r:5')
    >>> G = together(sum( [ g.subs(z, r) for r in R] ))
    >>> G
    (b + r₀)⋅(b + r₁)⋅(b + r₂)⋅(b + r₃) + (b + r₀)⋅(b + r₁)⋅(b + r₂)⋅(b + r₄) + (b + r₀)⋅(b + r₁)⋅(b + r₃)⋅(b + r₄) + (b + r₀)⋅(b + r₂)⋅(b + r₃)⋅(b + r₄) + (b + r₁)⋅(b + r₂)⋅(b + r₃)⋅(b + r₄)
    ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                                                                            (b + r₀)⋅(b + r₁)⋅(b + r₂)⋅(b + r₃)⋅(b + r₄)
    >>> V = viete(f, R, z)
    >>> for lhs, rhs in V:
    ...     pprint(Eq(lhs, rhs))
    ...
    r₀ + r₁ + r₂ + r₃ + r₄ = 0
    r₀⋅r₁ + r₀⋅r₂ + r₀⋅r₃ + r₀⋅r₄ + r₁⋅r₂ + r₁⋅r₃ + r₁⋅r₄ + r₂⋅r₃ + r₂⋅r₄ + r₃⋅r₄ = 0
    r₀⋅r₁⋅r₂ + r₀⋅r₁⋅r₃ + r₀⋅r₁⋅r₄ + r₀⋅r₂⋅r₃ + r₀⋅r₂⋅r₄ + r₀⋅r₃⋅r₄ + r₁⋅r₂⋅r₃ + r₁⋅r₂⋅r₄ + r₁⋅r₃⋅r₄ + r₂⋅r₃⋅r₄ = 0
    r₀⋅r₁⋅r₂⋅r₃ + r₀⋅r₁⋅r₂⋅r₄ + r₀⋅r₁⋅r₃⋅r₄ + r₀⋅r₂⋅r₃⋅r₄ + r₁⋅r₂⋅r₃⋅r₄ = 1
    r₀⋅r₁⋅r₂⋅r₃⋅r₄ = -a
    >>> p = expand(numer(G))
    >>> q = expand(denom(G))
    >>> (P, Q), mapping = symmetrize((p, q), R, formal=True)
    >>> P
    ⎛   4      3         2                    ⎞
    ⎝5⋅b  + 4⋅b ⋅s₁ + 3⋅b ⋅s₂ + 2⋅b⋅s₃ + s₄, 0⎠
    >>> Q
    ⎛ 5    4       3       2                  ⎞
    ⎝b  + b ⋅s₁ + b ⋅s₂ + b ⋅s₃ + b⋅s₄ + s₅, 0⎠
    >>> for s, poly in mapping:
    ...     pprint(Eq(s, poly))
    ...
    s₁ = r₀ + r₁ + r₂ + r₃ + r₄
    s₂ = r₀⋅r₁ + r₀⋅r₂ + r₀⋅r₃ + r₀⋅r₄ + r₁⋅r₂ + r₁⋅r₃ + r₁⋅r₄ + r₂⋅r₃ + r₂⋅r₄ + r₃⋅r₄
    s₃ = r₀⋅r₁⋅r₂ + r₀⋅r₁⋅r₃ + r₀⋅r₁⋅r₄ + r₀⋅r₂⋅r₃ + r₀⋅r₂⋅r₄ + r₀⋅r₃⋅r₄ + r₁⋅r₂⋅r₃ + r₁⋅r₂⋅r₄ + r₁⋅r₃⋅r₄ + r₂⋅r₃⋅r₄
    s₄ = r₀⋅r₁⋅r₂⋅r₃ + r₀⋅r₁⋅r₂⋅r₄ + r₀⋅r₁⋅r₃⋅r₄ + r₀⋅r₂⋅r₃⋅r₄ + r₁⋅r₂⋅r₃⋅r₄
    s₅ = r₀⋅r₁⋅r₂⋅r₃⋅r₄
    >>> subslist = [ (s, c) for (s, _), (_, c) in zip(mapping, V) ]
    >>> subslist
    [(s₁, 0), (s₂, 0), (s₃, 0), (s₄, 1), (s₅, -a)]
    >>> P[0].subs(subslist)/Q[0].subs(subslist)
         4
      5⋅b  + 1
    ───────────
          5
    -a + b  + b
    >>> # Note, we must give a variable to RootSum, as f has more than one Symbol
    >>> cancel(P[0].subs(subslist)/Q[0].subs(subslist) - RootSum(f, Lambda(z, g), z)) == 0
    True
