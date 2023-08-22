# Benchmarking

Runtime benchmarking was performed on the Setonix HPC system at the Pawsey Supercomputing Research Centre. 


'Local' mode refers to running Koverage on a single node without using the scheduler to run jobs in parallel. SLURM mode refers to using a Snakemake config file to allow jobs to be sent to the SLURM scheduler individually.  The 'Speedup' column shows the times speedup from using the SLURM mode (i.e. parallelism) rather than running jobs in serial. 


Max RSS is the peak memory used in MB.

Walltime is HH:MM:SS

## Coral Genome
9.1GB of read files
<!DOCTYPE html>
<html>
<body>

<table>
  <tr>
    <th></th>
    <th colspan="2">Local</th>
    <th colspan="2">SLURM</th>
    <th></th>
  </tr>
  <tr>
    <th>Module</th>
    <th>Max RSS</th>
    <th>Walltime</th>
    <th>Max RSS</th>
    <th>Walltime</th>
    <th>Speedup</th>
  </tr>
  <tr>
    <td>CoverM</td>
    <td>4602.67</td>
    <td>0:18:41</td>
    <td>4652.13</td>
    <td>0:10:12</td>
    <td>1.8x</td>
  </tr>
  <tr>
    <td>Kmer</td>
    <td>3966.83</td>
    <td>0:36:39</td>
    <td>3967.35</td>
    <td>0:16:57</td>
    <td>2.2x</td>
  </tr>
  <tr>
    <td>Map</td>
    <td>2066.34</td>
    <td>1:46:11</td>
    <td>2116.61</td>
    <td>0:17:07</td>
    <td>6.2x</td>
  </tr>
</table>

</body>
</html>

## Chardonnay Genome
1.6TB of read files


<!DOCTYPE html>
<html>
<body>

<table>
  <tr>
    <th></th>
    <th colspan="2">Local</th>
    <th colspan="2">SLURM</th>
    <th></th>
  </tr>
  <tr>
    <th>Module</th>
    <th>Max RSS</th>
    <th>Walltime</th>
    <th>Max RSS</th>
    <th>Walltime</th>
    <th>Speedup</th>
  </tr>
  <tr>
    <td>CoverM</td>
    <td>10534.47</td>
    <td>12:12:44</td>
    <td>10418.07</td>
    <td>1:10:21</td>
    <td>10.4x</td>
  </tr>
  <tr>
    <td>Kmer</td>
    <td>54756.71</td>
    <td>7:24:45</td>
    <td>54756.38</td>
    <td>2:08:15</td>
    <td>3.5x</td>
  </tr>
  <tr>
    <td>Map</td>
    <td>2441.05</td>
    <td>17:28:40</td>
    <td>2428.76</td>
    <td>0:58:55</td>
    <td>17.8x</td>
  </tr>
</table>

</body>
</html>
