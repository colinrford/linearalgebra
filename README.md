# linearalgebra

just another math library..

Built with unique\_ptr, compile with c++23. Currently only 'trivial' 
classes and methods (and all that) are implemented. Not a high performance 
library in its current state. 

It's a work-in-progress. Vector has iterator. Matrix uses c++23 
multidimensional [0, ..., k] index operator. Very incomplete still but 
writing some converters for Eigen and lapack++, unnecessary for everyone but 
myself I imagine... but will lazily yet greatly enhance its utility for the 
aesthetically challenged. Hopefully I can get around to setting up modules,
which would also greatly improve the quality of life of this single digit user 
library :0) Toying with the idea of doing constexpr unique_ptr things, as well 
as investigating shared_ptr and atomic<shared_ptr> or things like that. TBD.

2024

