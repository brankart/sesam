# Preliminary tests of efficiency on GPUs 

These test codes reproduces the structure of the code of the MCMC sampler:
same memory structure, same loops, but without input data
and without some intermediate computation (negligible in cost).

In the following:
 - C is the execution time obtained without GPU accelerators (CPU only).
 - G is the execution time obtained with GPU accelerators (CPU and GPU).

The performance is mainly driven by the size of the problem (N, M, Niter):
 - larger vector size N -> larger memory transfer from CPU to GPU and larger cost per iteration.
 - larger M -> larger cost per iteration and per kernel.
 - larger number of iteration Niter -> larger cost for given memory transfer.

## Performance on Jean-Zay (partition V100) for test5_mpi

This test is performed with MPI on 2 CPUs and 2 GPUs (1 GPU per CPU),
with Nvidia compiler (option mpigpu_nvi in compilation script).

 -  N=1e5, Niter=1e4: C=  7.35s, G= 0.72s, ratio=10
 -  N=1e5, Niter=1e5: C=  74.1s, G= 6.09s, ratio=13
 -  N=1e5, Niter=5e5: C= 373.5s, G=30.03s, ratio=12
 -  N=1e6, Niter=1e4: C=  74.0s, G= 1.39s, ratio=53
 -  N=1e6, Niter=5e4: C= 365.5s, G= 6.14s, ratio=60
 -  N=5e6, Niter=1e4: C= 370.0s, G= 4.34s, ratio=85
 -  N=5e6, Niter=1e5: C=3700.0s, G=40.54s, ratio=91
 -  N=1e7, Niter=1e4: C= 740.0s, G= 7.82s, ratio=95
 -  N=1e7, Niter=1e5: C=7400.0s, G=73.75s, ratio=100
 -  N=1e7, Niter=1e4: C= 740.0s, G= 7.58s, ratio=98
 -  N=1e7, Niter=1e3: C=  74.0s, G= 0.95s, ratio=78
 -  N=2e7, Niter=1e4: C=1480.0s, G=13.83s, ratio=107
 -  N=5e7, Niter=1e4: C=3700.0s, G=31.07s, ratio=119

 -  N=1e8, Niter=1e4: C=7400.0s, G=60.70s, ratio=122
 -  N=2e8, Niter=1e4: OUT OF MEMORY

From these test, we conclude:
 - A large size N is necessary for sufficient efficiency.
 - Minimum size of observation vector per core: N = 5e6
 - Maximum size of ensemble vectors for a 16Gb node
(for a 30-member ensemble): Nmax(node) = 4e7

If we use 4 cores per node, the max size of vector size per core is:
        Nmax(core) = 1e7

The max ratio between full vector and observation vector is thus:
	rho = Nmax(core)/N = 2

If rho=10 (N=1e6), the efficiency ratio will drop from ratio=85 to ratio=53

Typical run time at optimum (use number of nodes and cores to have that):
 - for 1e4 iterations: between T=  5s and T=  10s
 - for 1e5 iterations: between T=  1m and T=   2m
 - for 1e6 iterations: between T= 10m and T=  20m
 - for 1e7 iterations: between T=  1h and T=   2h

## Performance on Adastra (partition MI250) for test5_mpi

This test is performed with MPI on 2 CPUs and 2 GPUs (1 GPU per CPU),
with Cray compiler (option mpiomp_crayftn in compilation script).

  - N=2e6, Niter=1e4: C=  27s, G= 3.44s, ratio=8
  - N=5e6, Niter=1e4: C=  83s, G= 6.04s, ratio=14
  - N=1e7, Niter=1e4: C= 303s, G= 9.54s, ratio=32
  - N=1e7, Niter=1e4: C= 190s, G= 9.65s, ratio=20
  - N=2e7, Niter=1e4: C= 403s, G= 15.7s, ratio=26
  - N=1e8, Niter=1e3: C= 220s, G= 5.37s, ratio=41
  - N=1e8, Niter=1e3: C= 385s, G= 5.40s, ratio=71

From a discussion with ChatGPT about the comparison between the V100 and the MI250,
we provisionally conclude that:

 - The CPU is much less efficient on the V100 (as compared to the MI250),
which makes it easier to have a good ratio on the V100.

 - The MI250 would require more intensive computation within each kernel.
Ok, it is fine that there is a lot of computation on the same 'target data' (large Niter),
and this was key to have a good performance on the V100,
but there is still a latency time between the running of the kernels,
which may penalize the MI250. This might be tested by increasing M.

 - There is more memory on the MI250, so we can expect to divide
the big problem in larger parts, and thus just use less processors.
Futher test needed.

 - There may be a difficulty in adressing the big ensemble array (ens),
because of the use of a vector of indices (sample),
which is here necessary because it changes at each iteration in the true code.

 - It is normal that there are more fluctuations on C than on G (?). Make an average to have a reliable value of the ratio. But predictability of cost is also an asset in the applications.

 - This is just preliminary, there is still hope for improvement.
