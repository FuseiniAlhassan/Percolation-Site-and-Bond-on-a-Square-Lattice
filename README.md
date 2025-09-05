# Percolation-Site-and-Bond-on-a-Square-Lattice

**Abstract:**  
We present a Python‑based simulation framework for studying **site** and **bond percolation** on a two‑dimensional square lattice, with emphasis on cluster identification, spanning detection, and statistical characterization. The model implements both **site occupation** and **bond activation** using independent Bernoulli trials with adjustable probability \( p \), and employs an efficient **Union–Find (disjoint set)** data structure to label connected clusters and compute their sizes.

For site percolation, the code detects spanning clusters under open boundary conditions, identifying whether any cluster connects opposite edges of the lattice. For bond percolation, horizontal and vertical bonds are generated independently, and connectivity is determined via the same Union–Find approach. The framework supports **Monte Carlo sweeps** over \( p \) to estimate the percolation probability curve, revealing the critical threshold \( p_c \) characteristic of the square lattice.

Visualization tools include:
- Cluster‑colored lattice plots
- Log–log cluster size distributions
- Percolation probability curves
- An **animated growth sequence** showing the emergence of large‑scale connectivity as $$\( p \)$$ increases from 0 to 1

The implementation is self‑contained, requires only standard scientific Python libraries, and is designed for **educational use**, **research prototyping**, and **statistical physics outreach**. It provides a reproducible platform for exploring universality, scaling, and critical phenomena in percolation theory.


# Site and Bond Percolation on a Square Lattice

## Overview
This project simulates **site** and **bond percolation** on a 2D square lattice, detects spanning clusters, and visualizes cluster statistics and growth.  
It is implemented in **pure Python** with `NumPy` and `Matplotlib`, using an efficient **Union–Find** algorithm for cluster labeling.

The code can:
- Generate random site or bond configurations
- Identify connected clusters and their sizes
- Detect percolation (spanning) events
- Estimate percolation probability vs. occupation probability
- Visualize configurations, cluster size distributions, and percolation curves
- Animate the growth of occupied sites from $$\( p = 0 \)$$ to $$\( p = 1 \)$$


## Features
- **Site Percolation**: Occupied sites with probability $$\( p \)$$
- **Bond Percolation**: Activated horizontal/vertical bonds with probability $$\( p \)$$
- **Union–Find Cluster Detection**: Fast connectivity analysis
- **Spanning Cluster Detection**: Left–right and top–bottom percolation checks
- **Monte Carlo Sweeps**: Estimate percolation probability curves
- **Visualization**:
  - Cluster‑colored lattice plots
  - Log–log cluster size distributions
  - Percolation probability curves
  - Animated site occupation growth


## Requirements
- Python 3.x
- NumPy
- Matplotlib
- FFmpeg or ImageMagick (for saving animations)

Install dependencies:

pip install numpy matplotlib


## Usage

### Run the Demo
Percolation. Site and Bond on a Square Lattice.py

The demo will:
1. Generate a site percolation configuration near $$\( p_c \approx 0.5927 \)$$
2. Plot the configuration with clusters colored
3. Compute and plot the cluster size distribution
4. Sweep $$\( p \)$$ to estimate the percolation probability curve
5. Run a bond percolation example
6. Optionally animate site growth

## Key Functions

| Function | Purpose |
|----------|---------|
| `generate_site_configuration(L, p)` | Random site occupation |
| `generate_bond_configuration(L, p)` | Random bond activation |
| `site_percolation_clusters()` | Label site clusters |
| `bond_percolation_clusters()` | Label bond clusters |
| `has_spanning_cluster_site()` | Detect spanning in site percolation |
| `has_spanning_cluster_bond()` | Detect spanning in bond percolation |
| `percolation_probability_site()` | Monte Carlo percolation probability (site) |
| `percolation_probability_bond()` | Monte Carlo percolation probability (bond) |
| `plot_configuration_site()` | Visualize site configuration and clusters |
| `plot_cluster_sizes()` | Plot cluster size distribution |
| `plot_percolation_curve()` | Plot percolation probability vs. \( p \) |
| `animate_site_growth()` | Animate site occupation from \( p=0 \) to \( p=1 \) |

## Example Output
- **Cluster plot**: Sites colored by cluster ID
- **Cluster size distribution**: Power‑law behavior near $$\( p_c \)$$
- **Percolation curve**: S‑shaped transition near critical $$\( p \)$$
- **Animation**: Growth of occupied sites and emergence of spanning cluster


## Customization
- `L`: Lattice size
- `p`: Occupation probability
- `trials`: Number of Monte Carlo trials
- `p_steps`: Number of animation frames
- `seed`: Random seed for reproducibility

## References
- D. Stauffer & A. Aharony, *Introduction to Percolation Theory*, Taylor & Francis (1994)
- M. E. J. Newman & R. M. Ziff, "Efficient Monte Carlo algorithm and high-precision results for percolation," *Phys. Rev. Lett.*, 85, 4104 (2000)


## License
MIT License: free to use, modify, and distribute with attribution.
