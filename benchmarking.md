# Benchmarking

Runtime benchmarking was performed on the Setonix HPC system at the Pawsey Supercomputing Research Centre

## CoverM
### Using local
CPUs = 8 
Memory = 20GB

| Genome name | Genome size | size of reads files | walltime | peak mem GB |
--------------|-------------|---------------------|----------|-----------
| Coral       | XXXX        | 9.1GB               | 16m55s   | 4.89     |
| Chardonnay  | XXXX        | 1.6TB                | XXXX     | XXXX     |

### Using SLURM profile to allow jobs to spawn in parallel

| Genome name | Genome size | size of reads files | walltime | peak mem GB |
--------------|-------------|---------------------|----------|-----------
| Coral       | XXXX        | 9.1GB               | 2m03s   | 4.88     |
| Chardonnay  | XXXX        | 1.6TB               | 1h02m   | 10.98     |

## Map
### Using local
CPUs = 8 
Memory = 20GB
| Genome name | Genome size | size of reads files | walltime | peak mem GB |
--------------|-------------|---------------------|----------|-----------
| Coral       | XXXX        | 9.1GB               | 25m41s   | 2.33     |
| Chardonnay  | XXXX        | 1.6TB               | XXXX     | XXXX     |


### Using SLURM profile to allow jobs to spawn in parallel
| Genome name | Genome size | size of reads files | walltime | peak mem GB |
--------------|-------------|---------------------|----------|-----------
| Coral       | XXXX        | 9.1GB               | 2m55s   | 2.29     |
| Chardonnay  | XXXX        | 1.6TB                | XXXX     | XXXX     |

## Kmer
### Using local
#SBATCH --cpus-per-task=8
#SBATCH --time=1:00:00
#SBATCH --mem=10GB
| Genome name | Genome size | size of reads files | walltime | peak mem GB |
--------------|-------------|---------------------|----------|-----------
| Coral       | XXXX        | 9.1GB               |  37m16s  |  4.16    |
| Chardonnay  | XXXX        | 1.6TB                | XXXX     | XXXX     |


### Using SLURM profile to allow jobs to spawn in parallel
| Genome name | Genome size | size of reads files | walltime | peak mem GB |
--------------|-------------|---------------------|----------|-----------
| Coral       | XXXX        | 9.1GB               | 4m23s   | 1.69     |
| Chardonnay  | XXXX        | 1.6TB               | XXXX     | XXXX     |