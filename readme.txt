Simulation of STDP learning during homeostatic control, as described in [1].
- The 'decorrelation_learning.m' script simulates learning under decorrelation homeostasis
- The 'rate_control_learning.m' script simulates learning under rate-control homeostasis

Upon simulation completion, both scripts plot the real and imaginary parts of all eigenvalues of connectivity, across time.

The 'eigenshuffle' function is used within the simulation scripts, and should therefore be included the execution folder.

Pre-initialized weights and state-vector are stored in 'data.mat' and are automatically loaded by the scripts.




[1] Stable memory with unstable synapses, L. Susman, N. Brenner and O. Barak.