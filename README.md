# QUIJIBO: Quaternion Julia Set Shape Optimization

![](http://www.tkim.graphics/JULIA/images/bunnies_small.png)

This is an updated version of Dr. Theodore Kim's [QUIJIBO](http://www.tkim.graphics/JULIA/index.html) source code, modified to use [bempp-cl](http://bempp.com/) instead of [BEM++](https://github.com/bempp/bempp-legacy), which is now rather difficult to build and run. 

For information about usage of this program, see [the QUIJIBO source code page](http://www.tkim.graphics/JULIA/source.html). This modified version still runs the BEM solve via a system call, so to run the full pipeline (see "Usage: Advanced" on the source code page) you will need to do the following:
1. Install [`bempp-cl`](https://github.com/bempp/bempp-cl) and its dependencies (including Gmsh capabilities). You can do this using the `requirements.txt` file in `projects/bemWrapper` or manually, following the instructions on the [Bempp GitHub page](https://github.com/bempp/bempp-cl) or [bempp.org](bempp.org). `solve_potential.py` is written to be compatible with Bempp 0.2.3.
2. Modify [the system calls](https://github.com/alexaschor/QUIJIBO/blob/main/projects/bemWrapper/bemWrapper.cpp#L64-L65) in `bemWrapper.cpp` to activate the proper Python environment (if you installed Bempp in your main Python environment you can skip this part) and then specify the location of `solve_potential.py` on your disk. You can use a relative path for this, but note that it will be relative to your current working directory when you execute `./bemWrapper`, not relative to where you build QUIJIBO.
