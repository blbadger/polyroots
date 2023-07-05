## *** Sunsetted 11/2022 with the removal of free tier Heroku services ***

## Pfinder: Polynomial Root Finder

### Overview

The Pfinder uses one of three analytic methods to find roots of algebraic polynomial equations.  The root values found are displayed along with a map showing how many iterations each initial point in the complex plane took to arrive at a root: using the defualt color map, lighter color indicates more iterations until convergence.  In the general case, these maps are self-similar fractals and some can be exquisitely detailed.  

### How to use 

After following the link in the upper right hand corner of this page, a short server activation time ensues (this app runs on a free server which sleeps if not in use) and the following interface will appear:

![screenshot](/assets/pfinder_screenshot.png)

Specifically, this is what will be seen for Firefox: use of Chrome or Safari will lead to slightly different aesthetics. The Pfinder was designed for use by Firefox or Chrome and is not optimized for Safari, although it is still usable for that browser.

`Real Bounds` and `Imaginary Bounds` specify the set of initial points in the complex plane.  Real and imaginary components are mapped on the horizontal and vertical axes, respectively: for example, imaginary bounds `-1, 1` results in an image generated for initial values between `-i` (at the bottom of the map) and `i` (top).  Inputs may or may not contain spaces, but must have a single comma `,` between boundary elements, and the first element must be less than the second.

The `Choose color map:` dropdown specify the colors assigned to the array of the number of iterations before arriving at a root.  The`Choose resolution:` dropdown menu specifies two integers, which are the number of values calculated between the `Real Bounds` and the `Imaginary Bounds`, respectively. There are approximately equal to the resolution of the image produced.

Any of Newton's method, Halley's method, or the Secant method (with the initial guess plotted and the second guess being half way between the origin `0 + 0i` and the first guess) can be selected in the `Choose method:` dropdown.

The `Maximum iterations` input selection field records the maximum number of iterations of the root method before the programs halts, and can be fed any value between `0` and `300`.  Roughly speaking, more iterations are required for higher resolution images.  Note that the different root methods require different iterations to make interesting maps, with the Secant method requiring the most and Halley's the fewest.

The last input is the `Specify Equation` field, where the user specifies the equation that is then fed to the selected root method in order to generate a map of quickly- versus slowly-converging initial values. Inputs should contain no spaces, and only the characters `0123456789.ie^+-` should be entered.  The equation is parsed and rendered in MathJax to the right of the input field, and checking to make sure that the rendered equation matches the intended input is helpful to make sure that the parser recognizes the input properly.  

The example input equation `x^7.14-x-1` results in

![cover](/assets/pfinder_example.png)

Complex-valued constants and exponents are currrently supported as well, for example:

`x^(7.11+0.1i)-x-(1-0.2i)`

Complex numbers may be entered with or without real values, ie `(2i)` is acceptable alternative to `(0+2i)`, but all complex or negative contants and exponents must exist in single parentheses.  Note also complex-numbered equation inputs will lead to longer computations than real-valued equations.

Once all desired inputs have been made, press the `CLICK TO RUN` button to activate the program.  Expect to wait anywhere from a couple seconds for integer-powered polynomials at low resolution to over a minute for fractional or complex-valued polynomials at 4k resolution. Upon completion of the computation, root values are displayed above the map of convergence of initial values.  Note that root values are not currently displayed for the Secant method, and note also that this field allows for horizontal scrolling but that the scroll bar is practically invisible for iOS versions of Safari.

### Behind the scenes

The Pfinder is a Flask web app running on a gunicorn server using a Plotly Dash interface for front end layout and callbacks, with styling in CSS and HTML.  Array-based computation necessary for root finding and image generation is started by a worker and sent to a Redis server via a Redis Queue message broker, and the Redis server is pinged every two seconds by the app to see if computation is complete.  The worker computation is performed using Numpy, and the resulting array signifying rate of convergence is transformed into an image via Matplotlib and saved to a temporary memory buffer as a bytestring using a base64 binary encoding.  


When the Redis job is complete, it is fetched and the bytestring is decoded into a PNG that can be opened in a separate page for maximum resolution.  At the same time, the roots values are fetched and converted to mathematical notation via MathJax and displayed.

For clarity, the callback graph is as follows:

![cover](/assets/pfinder_graph.png)

This system of background Redis processes circumvents the problem of long computation times faced by complex-valued operations at very high resolution.  Heroku, Azure, and most other cloud PaaS providers have hard time limits (30s in this case) for safety and efficiency concerns, meaning that the long computations required to generate high-resolution images of complex-number arrays would simply time out without this system in place.  Timeouts do not occur even when computations run for more than 5 minutes with the current configuration because the app is continually 'active' as it pings the Redis server.  

The generation of a polynomial root map only occurs when the `CLICK TO RUN` button is pressed, but real-time callbacks occurs for rendering the specified equation with MathJax (with a delay of ~100ms).  This process does not employ redis but is computed by  helper functions, which allows for new equations to be rendered whilst a prior equation map is generated.

The Pfinder is hosted by Heroku, using a free version of a platform as a service in which a linux (specifically ubuntu) operating system is used as an app container, also called a dyno, in order to run the app worker.  The version of Redis used is Heroku Redis, which supports up to 20 clients and offers 25MB of memory, more than enough even for 4k image resolution of intricate fractals. This Redis service is free but generally superior to other free Redis instances with regards to memory capacity.


### Learn More

To learn more, see [this page](https://blbadger.github.io/polynomial-roots.html).  If this sort of thing is of interest to you, see also the [Jenerator](https://github.com/blbadger/jenerator)

