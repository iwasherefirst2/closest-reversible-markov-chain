# getRevMarkov
A matlab script to compute the nearest reversible Markov chain.

# What?
For any Markov chain exists according to a given norm and probability vector a unique closest reversible Markov chain.

![image](https://user-images.githubusercontent.com/1765602/68806665-9259a200-0666-11ea-943b-dc52b16fd2f1.png)

For the full theory check out the paper: https://www.researchgate.net/publication/271449135_Computing_the_nearest_reversible_Markov_chain or the [Wikipedia article](https://en.wikipedia.org/wiki/Markov_chain#Closest_reversible_Markov_chain)

# How to use it?
Go with Matlab in the folder where you have saved the file **getClosestSparse.m**.

Assume you have a NxN dimensonal Markov chain **A** and a N-dimensional probability vector **m** then you can compute the closest reversible NxN dimensonal Markov chain **U** according to the Frobenius-Norm with the command
```
U = getClosestSparse(A,m)
```

If you want that **A** and **U** have similar rare events (this implies that they have similar eigenvalues, see the above mentioned paper) then you can compute the closest reversible NxN dimensonal Markov chain **U** according to a special weighted norm with the command
```
U = getClosestSparse(A,m,1)
```

Note that the computation of the Markov chain **U** according to a special weighted norm is only usable if **m** has no zero entries. For the closest reversible Markov chain **U** according to the Frobenius-Norm zero entries are allowed.

