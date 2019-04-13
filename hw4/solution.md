# Question 1

| CUDA Server| CPU Bandwidth | GPU Bandwidth|
|---------|----------|--------------|
| 1       | 86.184295 GB/s|203.522325 GB/s |
| 2       | 84.506014 GB/s| 70.390337 GB/s |
| 3       | 43.548502 GB/s|87.935884 GB/s |
| 4       | 43.4347711 GB/s|80.211482 GB/s |
|5        | 63.407023 GB/s| 123.629778 GB/s |

We can see that the GPU Bandwidth is higher for all the cuda servers except for CUDA2.

# Question 2

We use a matrix size of 6402 x 6402 where the extra 2 entries are the 0 boundary conditions. The GPU kernel operates on the inner 6400 x 6400 matrix and solves a system iteratively where u0 = 0, and f = 1 .

| CUDA Server| GPU Time for 1000 iterations|
|---------|-----------|
|1 | 0.003321s |
|2 | 0.010319s |
|3 | 0.005153s|
|4 | 0.010976s|
|5 | 0.004444s|
