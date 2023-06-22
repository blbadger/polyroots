# Polynomial roots

Find which regions of the complex plane converge on a polynomial root, which do not, and what roots are found.  Evaluate trajectories and generate images for Newton's method, Halley's method, and the secant method in the complex plane.  

If you do have access to a graphic processing unit, see `/parallelized` for modules that implement root-finding methods using the `torch` library.  This is by far the fastest option, with hundreds of root finding iterations on a 2k by 2k array taking less than a second on an RTX 3060.

If you do not have a GPU, see `/optimized` for modules that have been optimized using `numexpr`, capable of iterating through tens of iterations on a 2k by 2k array per second on a mid-range CPU.

As an example, here is a map where light pixels take longer to settle on a root than dark, of x^13-x-1 using Halley's method:
![halley's](https://github.com/blbadger/blbadger.github.io/blob/master/newton-method/halley_x%5E13-x-1.png)

and Newton's method for x^7.11-x-1
![halley's](https://github.com/blbadger/blbadger.github.io/blob/master/newton-method/newton_z%5E7.11-z-1.png)
