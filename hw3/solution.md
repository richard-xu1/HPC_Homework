# Question 1


We vectorized the code by simply adding the next few terms in the Taylor Expansion for sine into the function for vectorized sine evaluation. Since the class already defined the vectorized operators for addition and multiplication, we didn't have to change the code for that. In order to extend sine to the rest of the real line, we used the identities sin(x) = cos(pi/2 - x) and sin(x) = -cos(pi/2 + x). Then we used the fmod function in C++ to convert anything outside of that range to within the [-pi/2, pi/2] interval. I wasn't sure how to vectorize fmod so I couldn't vectorize this version of sine. We see that the vectorized version of sin_taylor is much faster but the version extended to the real line is almost 100 times slower.

 # Question 2
We were interested in seeing how the scan function performs when parallelized over multiple threads.  This was done on an Intel i5-6300U.
  
| Threads | Run Time |
|---------|----------|
| 1       | 0.111440s|
| 2       | 0.156272s|
| 3       | 0.190271s|
| 4       | 0.178406s|

There was no increase in time for parallelizing the code. In fact, the algorithm actually slowed down despite increasing the number of cores used. This is likely because the problem has low computational intensity, since for every two read and 1 write, we are only doing 1 flop. The problem is bandwidth bound very quickly. 