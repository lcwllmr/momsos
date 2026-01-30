---
title: A Practical Introduction to Polynomial Optimization using the Moment-SOS Hierarchy
author: Luca Wellmeier
date: January 29, 2026
toc: true
abstract: |
  These notes are mostly based on the great survey article [(Laurent, 2010)](https://homepages.cwi.nl/~monique/files/moment-ima-update-new.pdf).
  This text and the code used to perform the experiments are hosted at [`gh:lcwllmr/momsos`](https://github.com/lcwllmr/momsos).
  Read the instructions on GitHub if you would like to experiment with the code on your own.
  These notes have been used for the following talks:

  - Jan 28, 2026 at IIT Bombay (17:15 in Ramanujan Hall) covering sections 1, 2 and 3
  - Feb 5, 2026 at IIT Bombay (??? in ???) covering sections 4 and 5
macros:
  '\R': '\mathbb{R}'
  '\deg': '\operatorname{deg}'
  '\set': '\left\{#1\;\colon\;#2\right\}'
  '\NN': '\mathcal{P}^+'
  '\SOS': '\Sigma'
  '\pmin': 'p^{\text{min}}'
  '\psos': 'p^{\text{sos}}'
---

## Non-negativity and sums of squares

Let $\R[x]$ be the polynomial ring in variables $x = (x_1, \dots, x_n)$ with real coefficients.
We start with the following lead problem:
decide whether a given polynomial $p \in \R[x]$ defines
a globally non-negative function (in short, $p \geq 0$).
While we expect this to be very hard to solve in general,
we may consider some special cases of polynomials.

> **Example (Univariate case)**:
> Let $n = 1$ and assume that $p \in \R[x]$ is non-negative.
> Necessarily, we must have that $\deg(p) = 2d$ is even,
> and that the leading coefficient of $p$ is positive
> (we might as well take $p$ monic).
> We also note that any real root $p$ might have must appear with
> even multiplicity as $p$ otherwise flip sign there.
> Consider a factorization of $p$ into real roots and complex
> conjugate root pairs:
> $$\begin{aligned}
>   p(x) &= \prod_i (x - r_i)^{2 m_i}
>   \prod_j (x - z_j)^{n_j} (x - \overline{z_j})^{n_j} \\
>   &= \underbrace{\left(\prod_i (x - r_i)^{m_i} \prod_j (x - z_j)^{n_j}\right)}_{g(x)}
>      \underbrace{\left(\prod_i (x - r_i)^{m_i} \prod_j (x - \overline{z_j})^{n_j}\right).}_{\overline{g(x)}}
> \end{aligned}$$
> Thus, with $g(x) = u(x) + i v(x)$, we obtain the expression
> $p(x) = g(x) \overline{g(x)} = u(x)^2 + v(x)^2$.

We have proven that every non-negative univariate polynomial is a sum of
squares of other polynomials.
Of course the reverse, i.e. being a sum of squares, trivially implies non-negativity.
Is it possible that non-negativity can be characterized by the sum-of-squares
(or abbreviated, *sos*) criterion? Let's see another case.

> **Example (Quadratic case)**: We now let $n$ be arbitrary and pick
> $p \in \R[x]$ such that $p \geq 0$ and $\deg(p) = 2$.
> We may write $$ p(x) = x^T A x + b^T x + c, $$
> where $A$ is a symmetric $(n \times n)$-matrix, be an $n$-vector and $c$ a constant.
> First, observe that the non-negativity of $p$ implies that $A$ is positive semidefinite
> (or *psd* for brevity):
> indeed, if this were not the case, there would be a direction $u$ s.t. $u^T A u < 0$
> and thus $p(tu) = t^2 (u^T A u) + O(t) \to -\infty$ as $t \to \infty$.
> Next, consider the homogenization $\tilde{p}$ of $p$ in vectorized form:
> $$ \tilde{p}(x, t) = x^T A x + b^T x t + c t^2
> = \begin{bmatrix} x & t \end{bmatrix}
> \underbrace{\begin{bmatrix} A & \frac 12 b \\ \frac 12 b^T & c \end{bmatrix}}_{Q}
> \begin{bmatrix} x \\ t \end{bmatrix}.$$
> Note that $\tilde p (x, 0) = x^T A x \geq 0$ since $A$ is *psd* and
> $\tilde p (x, t) = t^2 p(x / t) \geq 0$  since $p \geq 0$ whenever $t \neq 0$,
> so that $\tilde p \geq 0$ globally on $\R^n \times \R$.
> It follows that its defining matrix $Q$ is *psd* itself and we may consider a
> [Cholesky factorization](https://en.wikipedia.org/wiki/Cholesky_decomposition)
> $Q = U^T U$ to find
> $$ p(x) = \tilde p(x, 1) = \left( U \begin{bmatrix} x \\ 1 \end{bmatrix} \right)^T
> \left( U \begin{bmatrix} x \\ 1 \end{bmatrix} \right) = \left\| U \begin{bmatrix} x \\ 1 \end{bmatrix} \right\|^2.$$
> This last expression for $p$ is again a sum of squares of affine functions.

As for the univariate case, it turns out that non-negativity and being a *sos* is
equivalent in case of quadratic polynomials as well.
This should provide enough motivation to introduce dedicated symbols.

> **Definition**:
> The set of non-negative polynomials and the set of *sos* polynomials will be
> denoted respectively by
> $ \NN_n = \set{p \in \R[x]}{p \geq 0} $
> and $\SOS_n = \set{\sum_i q_i^2}{q_i \in \R[x]}$.
> Their respective truncations at degree $2d$ are denoted by by
> $\NN_{n,2d}$ and $\SOS_{n,2d}$.
> The subscript $n$ may be dropped if the number of variables is clear from the context.

Both $\NN$ and $\SOS$ are [convex cones](https://en.wikipedia.org/wiki/Convex_cone):
they are closed under addition and multiplication by non-negative scalars and clearly
$\SOS \subseteq \NN$.
We have seen two families for which we even have equality - how far can we go with this?
Not much further unfortunately - a result which motivated [Problem 17](https://en.wikipedia.org/wiki/Hilbert%27s_seventeenth_problem) in Hilbert's 1900 list of 23 problems.

> **Theorem (Hilbert, 1888)**:
> The equality $\NN_{n,2d} = \SOS_{n,2d}$ holds if and only if $(n,2d)$ is either of 
> $(1,2d)$, $(n,2)$ or $(2,4)$.

In fact, the situation looks rather hopeless:

> **Theorem [(Blekherman, 2006)](https://doi.org/10.1007/BF02771790)**:
> Fix an even degree $d$ and let $H$ be the compact set of all $p \in \NN_{n,d}$
> such that $\int_{\mathbb{S}^{n-1}} p \, d\sigma = 1$.
> Then, as $n \to \infty$
> $$ \frac{\operatorname{vol}(\SOS_{n,d} \cap H)}{\operatorname{vol}(\NN_{n,d} \cap H)} \to 0. $$

The paper gives exact asymptotic rates (think $n^{-d}$; very fast), but the point is clear:
the density of $\SOS$ in $\NN$ is vanishingly small in many variables or high degree.
Nonetheless, *sos*-based methods have proven very practical for a wide variety of applications.
The next section will justify this claim.

## Automated proofs

Even though Hilbert proved that the equality $\NN = \SOS$ holds only in some special cases,
the first explicit counterexample was considered much later by Motzkin in 1967:
$$ M(x,y) = x^4 y^2 + x^2 y^4 - 3 x^2 y^2 + 1. $$
When plotted as a function on $\R^2$ its surface looks like this:

![](./motzkin-surface.light.png) ![](./motzkin-surface.dark.png)

Proving that it is globally non-negative is a simple application of the
[AM-GM inequality](https://en.wikipedia.org/wiki/AM%E2%80%93GM_inequality):
$$ x^2 y^2 = \sqrt[3]{x^4 y^2 \; x^2 y^4 \; 1} \leq \frac{x^4 y^2 + x^2 y^4 + 1}{3}. $$
Rearranging yields $M \geq 0$.
The next result is the key to algorithmically deciding whether a given polynomial is *sos*.

> **Proposition (Parameterization of $\SOS$)**:
> Let $v_{n,d}(x) = (1, x_1, \dots, x_n, x_1^2, x_1 x_2, \dots, x_n^d)^T$
> be a monomial basis of the truncated polynomial ring $\R[x]_d$ written as column vector.
> A polynomial $p \in \R[x]_{2d}$ is *sos* if and only if there exists a
> *psd* matrix $Q$ such that
> $$ p(x) = v_{n,d}(x)^T Q v_{n,d}(x). $$

The exact choice of basis does not matter of course.
This matrix $Q$ will be refered to as the Gram matrix.

**Proof**:
If there exists such a $Q$ we may proceed as in the quadratic polynomial case above:
using a Cholesky factorization $Q = U^T U$ and plugging this in 
yields a sum-of-squares expression $p(x) = \| U v_{n,d}(x) \|^2$.
Conversely, if $p = \sum_i q_i^2$ is *sos*, we may write each $q_i$
(which is necessarily of degree at most $d$) in terms of the basis:
$q_i(x) = c_i^T v_{n,d}(x)$. Thus, $q_i(x)^2 = v_{n,d}(x)^T (c_i c_i^T) v_{n,d}(x)$
and by summing these up we get $p(x) = v_{n,d}(x)^T (\sum_i c_i c_i^T) v_{n,d}(x)$.
Finally, we observe that $Q = \sum_i c_i c_i^T$ is *psd* as a sum of rank-1 *psd* matrices.
$\square$

Checking whether the Motzkin polynomial $M(x,y)$ (or any polynomial for that matter)
is *sos* is therefore equivalent to checking whether the linear system
$$
  M(x,y) = v_{2,3}(x,y)^T Q v_{2,3} (x,y)
$$
obtained via coefficient matching has a *psd* solution $Q$.
In other words, it is equivalent to checking whether the intersection of an affine subspace
with the convex cone of *psd* matrices is non-empty.
A problem of this type can be solved in practice by
[semidefinite programming](https://en.wikipedia.org/wiki/Semidefinite_programming)
(short, *sdp*) solvers.
The experiment [`motzkin_not_sos.py`](https://github.com/lcwllmr/momsos/blob/main/code/momsos/experiments/motzkin_not_sos.py) computes exactly this - see the first of the two tests.

```
(.venv) $ motzkin-not-sos

Starting test: Is M(x,y) a sum of squares?  ...  Computation finished.
  Solution status: infeasible
  Solve time: 0.00112448s
  Iterations: 9

Starting test: Is (x^2 + y^2) * M(x,y) a sum of squares?  ...  Computation finished.
  Solution status: optimal
  Solve time: 0.005551222s
  Iterations: 12
```

This tells us that the the decomposition failed and that the Motzkin polynomial is not *sos*.
Some implementations of the so-called [interior point method](https://en.wikipedia.org/wiki/Interior-point_method) will also return a "dual infeasibility certificate":
To understand what that means exactly we need to study the conic dual of $\SOS$.
For now it is enough to think about such certificates as
a concrete hyperplane in some closely related cone that separates the compact feasible region
from the point of interest. More on this later.

So Motzkin is now verfiably not *sos*.
However, according to the second test that is performed in the experiment,
it turns out that the product of Motzkin with some other *sos* polynomial
results in something *sos*.
This is not just a coincidence.
Hilbert's result from the last section motivated a new question (which became Problem 17):
Is every non-negative polynomial a sum of squares of rational functions?
This has been answered positively by [(Artin, 1927)](https://doi.org/10.1007%2FBF02952513).
In our case, $M(x,y)$ is not *sos*, but it is fraction of two *sos* polynomials:
the denominator is $x^2 + y^2$ and the numerator was provided by the Gram matrix $Q$,
which was computed by the *sdp* solver.

In fact, one does not even have to guess the denominator as we did here.
The documentation of the Julia package [`SumsOfSquares.jl`](https://jump.dev/SumOfSquares.jl/stable/generated/Getting%20started/motzkin/) demonstrates how to build a program that not only
searches the decomposition of the product but also the denominator in one run given some
upper bound on the degree.
Let's rephrase this:
Any polynomial inequality can - in principle - be automatically proven or disproven
by means of semidefinite programming.
If it is true, then the *sdp* solver will provide a proof in the shape of a Gram matrix,
which is equivalent to an explicit *sos* decomposition, which in turn, proves non-negativity.
On the other hand, if the inequality is false, then the solver may return an infeasibility
certificate as we have seen before.

As powerful as this sounds, there are two important caveats.
For one, it is always important to be aware of possible issues with accuracy since
interior point methods are implemented using numerical linear algebra routines.
Moreover, the dimensions of such *sdps* explode very quickly if constructed naively.
The monomial vector $v_{n,d}(x)$, for instance, will have length $\binom{n + d}{n}$,
and an interior point solver itself will have rather high runtime costs being a
second-order method working on large matrices.
The research communities around these methods have found various tricks to circumvent cost
or trade it off with accuracy
(e.g. by exploiting more structure of problem instances or using more specialized solvers).


## Optimization over semialgebraic sets

In this section we will localize the non-negativity to sets described
by polynomial inequalities and formulate a effective procedure for solving
optimization problems over them.

Let $g_1(x), \dots, g_m(x) \in \R[x]$ be arbitrary polynomials and let them define
the basic semialgebraic set
$$
  K = \set{x \in \R^n}{g_1(x) \geq 0, \dots, g_m(x) \geq 0}.
$$
Such sets $K$ may take a wide variety of shapes. Here are just a few basic compact examples:

- closed unit ball: $g_1(x) = 1 - \| x \|^2$
- spherical shell: $g_1(x) = 2 - \| x \|^2$ and $g_2(x) = \| x \|^2 - 1$
- hypercube: $g_i (x) = 1 - x_i^2$ for $i = 1, \dots, n$
- simplex: $g_i(x) = x_i$ for $i = 1, \dots, n$ and $g_{n+1}(x) = 1 - \sum_{i=1}^n x_i$
- binary hypercube: $g_i(x) = x_i (x_i - 1)$ and $g_{n + i}(x) = -g_i (x)$ for $i = 1, \dots, n$

Note that this last example illustrates how to achieve equality constaints as well: simply add the negative of the polynomial to the mix.
We will assume throughout that $K$ is compact.

Next, consider a polynomial objective function $p \in \R[x]$ and consider the polynomial optimization problem (or *pop* for short)
$$ \pmin = \inf_{x \in K} p(x). $$

Problems of this rather innocent looking form encode a large variety of problems appearing in practice:
for starters, note that this already includes the whole classes of linear and quadratically constrained programming (themselves classes of broad application) as well as that of non-linear programming in binary variables which is of high interest to the field of combinatorial optimization.
Research have applied the method we are going to discuss in domains like power systems design, dynamical systems analysis, quantum information theory and to find sharp approximations to many famous NP-hard problems like [MaxCut](https://en.wikipedia.org/wiki/Maximum_cut) or [StableSet](https://en.wikipedia.org/wiki/Maximal_independent_set).

While *pops* are hard, non-convex problems in general, we can reformulate them into a convex problem (albeit an infinite dimensional one) by using the language of non-negativity.
Define the new convex cone $\NN(K) = \set{p \in \R[x]}{p(x) \geq 0 \; \forall x \in K}$
(exercise: verify that this is really a convex cone).
With it, reformulate
$$ \pmin = \sup \set{\lambda \in \R}{p - \lambda \in \NN(K)}. $$
Indeed, let $x^* \in K$ be a minimizer in the first formulation (remember: $K$ was assumed to be compact) and $\lambda^*$ be the maximizer of the second formulation.
Then, the polynomial $p - \lambda^*$ is non-negative on $K$.
If $p(x^*) - \lambda^*$ were strictly positive on $K$, it would contradict the optimality of $\lambda^*$, and therefore $\pmin = p(x^*) = \lambda^*$.

How lucky for us, that we already know very well how to relax non-negativity - let's go back to sums of squares.

> **Definition**: The quadratic module associated to the polynomials $g_1, \dots, g_m$ is defined as
> $$ Q(g_1, \dots, g_m) = \set{\sigma_0 + \sum_{j = 1}^{m} \sigma_j g_j }{\sigma_0, \sigma_1 \dots, \sigma_m \in \SOS}. $$
> For a given degree $2d$ we also define the truncated quadratic module
> $$ Q(g_1, \dots, g_m)_{2d} = \set{p \in Q(g_1, \dots, g_m)}{\deg(p) \leq 2d}. $$

The purpose is clear: any $p \in Q(g_1, \dots, g_m)$ is non-negative on $K$ since $g_i(x) \geq 0$ on $K$ for all $i = 1, \dots, m$.
Thus, similar to how $\SOS$ relaxes $\NN$, the new convex cone $Q(g_1, \dots, g_m)$ relaxes $\NN(K)$.

Note that membership in the truncated quadratic module will only certify non-negativity on $K$ if all of the summands appear. The truncation degree $2d$ must be chosen to be greater or equal $\max(\deg(g_1), \dots, \deg(g_m))$.

Moreover, observe that we declared explicit dependence on the defining inequalities instead of $K$ itself.
Different descriptions may result in different quadratic modules, but we can only write down an effective certificate if we know a description.

Now we are ready to define the sums-of-squares hierarchy:
$$ \psos_d = \sup \set{\lambda \in \R}{p - \lambda \in Q(g_1, \dots, g_m)_{2d}}. $$
These define a sequence of convex, finite-dimensional optimization problems.
As we increase the hierarchy level $d$, the search space grows and we find ourselves in a chain of inequalities
$$ \psos_1 \leq \psos_2 \leq \psos_3 \leq \dots \leq \pmin. $$

The first immediate question is: will this converge? That is, will $\psos_d \to \pmin$ as $d \to \infty$? Lucky us again:

> **Theorem (Putinar's Positivstellensatz, 1993)**:
> Assume that $K$ is compact and the quadratic module $Q = Q(g_1, \dots, g_m)$ is Archimedean, i.e. if $R - \|x\|^2 \in Q$ for some $R > 0$.
> If $p \in \NN(K)$, then $p + \varepsilon \in Q$ for all $\varepsilon > 0$.

The word Positivstellensatz was chosen in analogy to the German Nullstellensatz in algebraic geometry.
A modern, elementary proof of this theorem can be found in [(Schweighofer, 2005)](https://doi.org/10.1137/S1052623403431779).

Let us first have a look at the Archimedean condition, which is a technical assumption only and closely related to the discussion above: the quadratic module depends on the description of the set and not the set itself.
All it says is that you need to have such a ball constraint in the module for guaranteed convergence.
Practically, this is not a deal-breaker.
If you already know (maybe approximately) diameter and location of $K$, or if you have an a-priori bound on where the minimizers should be, it is enough to add this inequality without changing anything.

So how does this help us?
Think of $\varepsilon$ as the sharpness of the bound.
For any desired sharpness, Putinar tells us that we will be able to decompose a polynomial $p - \lambda$ into the form required for membership in the quadratic module as long as we are willing to wiggle $\lambda$ by $\varepsilon$.
Therefore, as we go up in hierarchy level, i.e. in truncation degree, we will be able to tighten our bounds, i.e. to let $\varepsilon \to 0$, because we gradually enlarge the search space.
This proves asymptotic convergence of the *sos* hierarchy.

Since we did not really live the framework of the last section, we are again able to use semidefinite programming to solve such optimization problems.
Let's take again the Motzkin polynomial $M(x,y)$ from the previous section and try to minimize it on a non-centered ball.
In other words, we attempt to solve the problem
of maximizing a scalar $\lambda$ such that the linear coefficient matching system 
$$M(x,y) + \lambda = \sigma_0 + \sigma_1 (r^2 - (x - x_0)^2 - (y - y_0)^2)$$
is satisfied and such that the the two polynomials $\sigma_0$ and $\sigma_1$ are *sos*, i.e. can be described by two *psd* matrix variables.
Before, we were dealing only with feasibility and now we also need to optimize a linear function at the same time.
Also, note that the first hierarchy level that makes sense here is $d=3$ - anything below and the coefficient matching system would be over-determined.
This problem is implemented in the experiment [`motzkin_minimize.py`](https://github.com/lcwllmr/momsos/blob/main/code/momsos/experiments/motzkin_minimize.py); check out the code to see how to hand these data to a solver.
Here is the result:

![](./motzkin-minimize.light.png) ![](./motzkin-minimize.dark.png)

Note that this is converging to the minimum (or rather to numerical precision) very fast even though the Archimedean condition is not satisfied.
