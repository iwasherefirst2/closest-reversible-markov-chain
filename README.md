# What is this?
A matlab script to compute the nearest reversible Markov chain.

# What?
For any Markov chain exists according to a given norm and probability vector a unique closest reversible Markov chain.
For example, consider the followin Markov chain:

<p align="center">
  <img src="https://user-images.githubusercontent.com/1765602/68806995-4b1fe100-0667-11ea-9356-fe204e537432.png">
</p>

This Markov chain is not reversible. According to the Frobenius norm the closest reversible Markov chain according to 

<p align="center">
  <img src="https://user-images.githubusercontent.com/1765602/68806912-21ff5080-0667-11ea-94e3-c56ec44956a9.png">
</p>

can be computed as 

<p align="center">
  <img src="https://user-images.githubusercontent.com/1765602/68807042-62f76500-0667-11ea-945b-1dbdd1221088.png">
</p>

If we choose the probability vector as 

<p align="center">
  <img src="https://user-images.githubusercontent.com/1765602/68807075-730f4480-0667-11ea-8fe8-566536cdc48e.png">
</p>

then the closest reversible Markov chain according to the Frobenius norm is approximately given by 

<p align="center">
  <img src="https://user-images.githubusercontent.com/1765602/68807109-81f5f700-0667-11ea-9f2f-d529c0cb487e.png">
</p>

For the full theory check out the [paper](https://www.researchgate.net/publication/271449135_Computing_the_nearest_reversible_Markov_chain ) or the [Wikipedia article](https://en.wikipedia.org/wiki/Markov_chain#Closest_reversible_Markov_chain)

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

