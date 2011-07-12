
.. include:: definitions.def

=====================
Solutions for Part II
=====================

Partial fraction decomposition
==============================

1.
  * `\frac{3 x + 5}{(2 x + 1)^2}`::

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
