# Benchmarking

Runtime benchmarking was performed on the Setonix HPC system at the Pawsey Supercomputing Research Centre

## CoverM
### Using local
coral: CPUs = 8 Memory = 20GB
Charonnay: CPUs=8, Mem 32GB

| Genome name | Genome size | size of reads files | walltime | peak mem GB |
--------------|-------------|---------------------|----------|-----------
| Coral       | XXXX        | 9.1GB               | 16m55s   | 4.89     |
| Chardonnay  | XXXX        | 1.6TB               | 12h04m   | 10.31     |

### Using SLURM profile to allow jobs to spawn in parallel

| Genome name | Genome size | size of reads files | walltime | peak mem GB |
--------------|-------------|---------------------|----------|-----------
| Coral       | XXXX        | 9.1GB               | 2m03s   | 4.88     |
| Chardonnay  | XXXX        | 1.6TB               | 1h02m   | 10.98     |

## Map
### Using local
CPUs = 8 Memory = 20GB for coral
| Genome name | Genome size | size of reads files | walltime | peak mem GB |
--------------|-------------|---------------------|----------|-----------
| Coral       | XXXX        | 9.1GB               | 25m41s   | 2.33     |
| Chardonnay  | XXXX        | 1.6TB               |      |  | ##REMEMBER TO ALSO ADD IN JOB 3528859 THAT GOT KILLED

Job ID: 3528859
Cluster: setonix
User/Group: sbeecroft/sbeecroft
State: CANCELLED (exit code 0)
Nodes: 1
Cores per node: 32
CPU Utilized: 3-10:34:46
CPU Efficiency: 3.72% of 92-09:36:00 core-walltime
Job Wall-clock time: 2-21:18:00
Memory Utilized: 30.02 GB
Memory Efficiency: 75.06% of 40.00 GB


### Using SLURM profile to allow jobs to spawn in parallel
| Genome name | Genome size | size of reads files | walltime | peak mem GB |
--------------|-------------|---------------------|----------|-----------
| Coral       | XXXX        | 9.1GB               | 2m55s   | 2.29     |
| Chardonnay  | XXXX        | 1.6TB                | XXXX     | XXXX     |

## Kmer
### Using local
for coral: #SBATCH --cpus-per-task=8 --time=1:00:00 --mem=10GB
| Genome name | Genome size | size of reads files | walltime | peak mem GB |
--------------|-------------|---------------------|----------|-----------
| Coral       | XXXX        | 9.1GB               |  37m16s  |  4.16    |
| Chardonnay  | XXXX        | 1.6TB                | 14h50m  | 53.52    |


### Using SLURM profile to allow jobs to spawn in parallel
| Genome name | Genome size | size of reads files | walltime | peak mem GB |
--------------|-------------|---------------------|----------|-----------
| Coral       | XXXX        | 9.1GB               | 4m23s   | 1.69     |
| Chardonnay  | XXXX        | 1.6TB               | XXXX     | XXXX     |
