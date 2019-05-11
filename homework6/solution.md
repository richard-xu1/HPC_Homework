# Question 1

Weak Scaling Study

1000 iterations on prince. Nl = 100. We ran our jacobi-mpi code on 1,4 and 16 processors. We did not do 64 because we could not get the 
64 nodes from the scheduler

|number of nodes| 1 | 4 | 16|
|---------------|--|----|---|
|time|0.024622s |0.033896s |0.052960s|


Strong Scaling Study

N = 1600. 1000 iterations.

|number of nodes| 1 |4 |16 |
|-|-|-|-|
|time| 4.253198s| 1.170796s| 0.352467s|
|linear time| 4.253198s|1.0632295s|0.2658249s|

# Question 2

Timing with 8 nodes and 8 tasks per node on prince0. 

|N|Time|
|-|-|
|1e4| 1.967086s|
|1e5|2.270425s |
|1e6| 3.02500s|
