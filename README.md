# linearalgebra

The goal behind this repository is to create from the ground up a robust, powerful linear algebra C++ library - it is a test of my understanding of linear algebra and the C++ language. My vision for this is different than that of, e.g., Matlab's; I want to eventually introduce this as both a tool for professional mathematicians as well as a learning source for the young mathematician. I realize that there are existing libraries or softwares out there that share (some of) this vision, however, I hope to provide perhaps a more intuitive interface that reduces the time involved with various computations.

-- update (Feb 2020)
the repo has been untouched for some time, though I have been trying to make some improvements lately, and will continue to do so. Part of this motivation is to refresh C++ for some numerical work I am going to be doing soon (will share more details later, likely in a few months). Another aspect is that Metal is based off of C++, and so my hope is that I can familiarize myself a bit with memory management, and additionally that these skills are portable from C++ to Metal. Let's see.

-- update (Mar 2020)
currently classes vector and matrix are functioning. Trying to implement newer C++ practices. Matrix does not use the vector class written here and instead uses std::vector. The vector class here as of this time is not robust enough to build matrices with, though hopefully it will be used (by me) for matrices someday! Have not done extensive testing on either class. Up next will be adding methods for finding inverse, transpose, determinant of a matrix. Also hope to improve error handling and stick to various common practices. By May I would like this to be ready for concepts / c++20.  
// been compiling with clang++ 17
