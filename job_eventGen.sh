#!/bin/bash
SBATCH --job-name=2p2hEventGen
SBATCH --output=output_%j.log
SBATCH --error=error_%j.log
SBATCH --ntasks=1
SBATCH --time=12:00:00
SBATCH --partition private-dpnc-cpu,shared-cpu
SBATCH --mem=10G

module load GCC/6.4.0-2.28
module load OpenMPI/2.1.2
module load ROOT/6.26.10


cd /path/to/your/executable

./build/program 0 0 0 0 2000 200 0 0 0 0 2000 100