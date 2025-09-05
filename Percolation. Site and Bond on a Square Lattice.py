#Percolation: site and bond on a square lattice
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from collections import Counter
import random
import math

# Union-Find (Disjoint Set)

class UnionFind:
    def __init__(self, n):
        self.parent = np.arange(n, dtype=int)
        self.rank = np.zeros(n, dtype=int)
        self.size = np.ones(n, dtype=int)  # track cluster sizes

    def find(self, a):
        # path compression
        while self.parent[a] != a:
            self.parent[a] = self.parent[self.parent[a]]
            a = self.parent[a]
        return a

    def union(self, a, b):
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return ra
        # union by rank
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        self.size[ra] += self.size[rb]
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1
        return ra

    def get_size(self, a):
        return self.size[self.find(a)]


# Helper indexing utilities

def idx(i, j, L):
    return i * L + j

def neighbors_4(i, j, L):
    # periodic? by default we use open boundaries for spanning checks; toggling is easy
    out = []
    if i > 0:        out.append((i-1, j))
    if i < L - 1:    out.append((i+1, j))
    if j > 0:        out.append((i, j-1))
    if j < L - 1:    out.append((i, j+1))
    return out


# Site percolation cluster finder

def site_percolation_clusters(occupied, L):

   # occupied: boolean array shape (L,L) indicating occupied sites.
    #Returns:
      #uf: UnionFind object
      #root_map: array shape (L,L) with root id or -1 for empty

    N = L*L
    uf = UnionFind(N)
    root_map = -np.ones((L, L), dtype=int)

    for i in range(L):
        for j in range(L):
            if not occupied[i, j]:
                continue
            root_map[i, j] = idx(i, j, L)
            # join with previously-seen neighbors (up and left) to avoid double work
            if i > 0 and occupied[i-1, j]:
                uf.union(idx(i, j, L), idx(i-1, j, L))
            if j > 0 and occupied[i, j-1]:
                uf.union(idx(i, j, L), idx(i, j-1, L))
    # compress roots and produce final map
    for i in range(L):
        for j in range(L):
            if root_map[i, j] >= 0:
                root_map[i, j] = uf.find(root_map[i, j])
    return uf, root_map


# Bond percolation cluster finder

def bond_percolation_clusters(bonds_h, bonds_v, L):
    
    #bonds_h: horizontal bonds present (L x (L-1)) connecting (i,j) to (i,j+1)
    #bonds_v: vertical bonds present ((L-1) x L) connecting (i,j) to (i+1,j)
    #Returns:
      #uf: UnionFind object
      #root_map: array shape (L,L) with root id (always all sites)
    
    N = L*L
    uf = UnionFind(N)
    # iterate bonds and union connected sites
    for i in range(L):
        for j in range(L-1):
            if bonds_h[i, j]:
                uf.union(idx(i, j, L), idx(i, j+1, L))
    for i in range(L-1):
        for j in range(L):
            if bonds_v[i, j]:
                uf.union(idx(i, j, L), idx(i+1, j, L))
    root_map = np.zeros((L, L), dtype=int)
    for i in range(L):
        for j in range(L):
            root_map[i, j] = uf.find(idx(i, j, L))
    return uf, root_map


# Percolation test (spanning cluster)

def has_spanning_cluster_site(occupied, root_map, L, wrap=False):

   # Check whether there exists a cluster connecting top to bottom or left to right.
   # If wrap=True, periodic boundary conditions are used conceptually for detection (wrap detection),
   # otherwise we check spanning in open boundaries.
   # Returns:
      #percolates_lr, percolates_tb, root_ids (set of roots that span)
    
    if wrap:
        # For wrap detection in site percolation, more complex algorithms used.
        # Here: we detect whether a cluster connects leftmost column to rightmost column (wrap)
        left_roots = set(root_map[:, 0][occupied[:, 0]])
        right_roots = set(root_map[:, -1][occupied[:, -1]])
        top_roots = set(root_map[0, :][occupied[0, :]])
        bottom_roots = set(root_map[-1, :][occupied[-1, :]])
    else:
        # open boundaries: collect roots touching each side
        left_roots = set(root_map[:, 0][occupied[:, 0]]) if np.any(occupied[:,0]) else set()
        right_roots = set(root_map[:, -1][occupied[:, -1]]) if np.any(occupied[:,-1]) else set()
        top_roots = set(root_map[0, :][occupied[0, :]]) if np.any(occupied[0,:]) else set()
        bottom_roots = set(root_map[-1, :][occupied[-1, :]]) if np.any(occupied[-1,:]) else set()

    percolates_lr = len(left_roots & right_roots) > 0
    percolates_tb = len(top_roots & bottom_roots) > 0
    percolating_roots = (left_roots & right_roots) | (top_roots & bottom_roots)
    return percolates_lr, percolates_tb, percolating_roots

def has_spanning_cluster_bond(bonds_h, bonds_v, L):

    #For bond percolation: generate uf and root map and then use same side-checks.

    uf, root_map = bond_percolation_clusters(bonds_h, bonds_v, L)
    # every site exists; identify roots touching sides
    left_roots = set(root_map[:, 0])
    right_roots = set(root_map[:, -1])
    top_roots = set(root_map[0, :])
    bottom_roots = set(root_map[-1, :])
    percolates_lr = len(left_roots & right_roots) > 0
    percolates_tb = len(top_roots & bottom_roots) > 0
    percolating_roots = (left_roots & right_roots) | (top_roots & bottom_roots)
    return percolates_lr, percolates_tb, percolating_roots, root_map, uf


# Utilities to generate random configurations

def generate_site_configuration(L, p, seed=None):
    rng = np.random.RandomState(seed)
    occupied = rng.rand(L, L) < p
    return occupied

def generate_bond_configuration(L, p, seed=None):
    rng = np.random.RandomState(seed)
    bonds_h = rng.rand(L, L-1) < p
    bonds_v = rng.rand(L-1, L) < p
    return bonds_h, bonds_v

# Analyze cluster sizes

def cluster_size_distribution(root_map, occupied_mask=None):
    
    #Returns Counter of cluster sizes (root id -> size)
    #If occupied_mask is provided (site percolation), only include occupied sites; otherwise include all.
    
    if occupied_mask is not None:
        roots = root_map[occupied_mask]
    else:
        roots = root_map.flatten()
    counts = Counter(roots)
    # remove invalid root id -1 if present (empty sites)
    if -1 in counts:
        del counts[-1]
    sizes = sorted(counts.values(), reverse=True)
    return sizes, counts

# Sweep over p to estimate percolation probability

def percolation_probability_site(L, p_vals, trials=50, seed=None):
    rng = np.random.RandomState(seed)
    perc_probs = []
    for p in p_vals:
        count = 0
        for t in range(trials):
            occupied = generate_site_configuration(L, p, seed=rng.randint(0,2**31-1))
            uf, root_map = site_percolation_clusters(occupied, L)
            lr, tb, roots = has_spanning_cluster_site(occupied, root_map, L)
            if lr or tb:
                count += 1
        perc_probs.append(count / trials)
    return np.array(perc_probs)

def percolation_probability_bond(L, p_vals, trials=50, seed=None):
    rng = np.random.RandomState(seed)
    perc_probs = []
    for p in p_vals:
        count = 0
        for t in range(trials):
            bonds_h, bonds_v = generate_bond_configuration(L, p, seed=rng.randint(0,2**31-1))
            lr, tb, roots, root_map, uf = has_spanning_cluster_bond(bonds_h, bonds_v, L)
            if lr or tb:
                count += 1
        perc_probs.append(count / trials)
    return np.array(perc_probs)


# Visualization helpers

def plot_configuration_site(occupied, root_map=None, title=None, cmap='tab20'):
    L = occupied.shape[0]
    plt.figure(figsize=(5,5))
    if root_map is None:
        plt.imshow(occupied, cmap='gray_r', origin='lower')
    else:
        # Map root ids to dense labels for color mapping
        unique_roots = np.unique(root_map[root_map >= 0])
        mapping = {r: i for i, r in enumerate(unique_roots)}
        colored = -np.ones_like(root_map)
        for i in range(L):
            for j in range(L):
                if root_map[i,j] >= 0:
                    colored[i,j] = mapping[root_map[i,j]]
        plt.imshow(colored, cmap=cmap, origin='lower')
    plt.title(title or "Site configuration")
    plt.axis('off')
    plt.savefig("site configuration")
    plt.show()

def plot_cluster_sizes(sizes, title="Cluster sizes"):
    if len(sizes) == 0:
        print("No clusters to plot.")
        return
    plt.figure(figsize=(6,4))
    counts = Counter(sizes)
    xs = np.array(sorted(counts.keys()))
    ys = np.array([counts[x] for x in xs])
    plt.loglog(xs, ys, 'o-')
    plt.xlabel("Cluster size (s)")
    plt.ylabel("Number of clusters of size s")
    plt.title(title)
    plt.grid(True, which='both', ls='--', alpha=0.3)
    plt.savefig("clusters size pot")
    plt.show()

def plot_percolation_curve(p_vals, perc_probs, label=None):
    plt.figure(figsize=(6,4))
    plt.plot(p_vals, perc_probs, 'o-')
    plt.xlabel("Occupation probability p")
    plt.ylabel("Percolation probability")
    if label:
        plt.title(label)
    plt.grid(True)
    plt.savefig('percolavation curve')
    plt.show()
    
# Animation: occupation growth for site percolation

def animate_site_growth(L=100, p_steps=50, seed=None, savefile=None):
    
    #Animate site occupation by gradually increasing p from 0 to 1 in p_steps increments.
    #If savefile ends with .mp4 or .gif, the animation will be saved (requires ffmpeg or imagemagick).
    
    rng = np.random.RandomState(seed)
    fig, ax = plt.subplots(figsize=(5,5))
    im = ax.imshow(np.zeros((L,L)), cmap='viridis', origin='lower')
    ax.axis('off')

    # Pre-generate random values and then threshold with growing p to keep correspondence
    rand_grid = rng.rand(L, L)
    p_list = np.linspace(0.0, 1.0, p_steps)

    def update(k):
        p = p_list[k]
        occupied = rand_grid < p
        im.set_data(occupied.astype(float))
        ax.set_title(f"Site occupation p={p:.3f}")
        return [im]

    ani = animation.FuncAnimation(fig, update, frames=p_steps, blit=True, interval=200, repeat=False)

    if savefile:
        if savefile.endswith('.mp4'):
            ani.save(savefile, writer='ffmpeg', fps=10)
        elif savefile.endswith('.gif'):
            ani.save(savefile, writer='imagemagick', fps=10)
        print(f"Saved animation to {savefile}")
        plt.close(fig)
    else:
        plt.show()


# Small demonstration and tests

def demo():
    L = 60
    seed = 12345
    print("Demo: site percolation L=", L)

    # Single instance near critical p ~ 0.5927 for site percolation on square lattice
    p = 0.59
    occupied = generate_site_configuration(L, p, seed=seed)
    uf, root_map = site_percolation_clusters(occupied, L)
    lr, tb, roots = has_spanning_cluster_site(occupied, root_map, L)
    print(f"Percolates (LR or TB)? {lr or tb}")

    # Plot occupied sites and cluster coloring
    plot_configuration_site(occupied, root_map=root_map, title=f"Site configuration p={p:.3f}")

    # Cluster sizes
    sizes, counts = cluster_size_distribution(root_map, occupied_mask=occupied)
    print(f"Top cluster sizes (site): {sizes[:10]}")
    plot_cluster_sizes(sizes, title="Cluster size distribution (site)")

    # Sweep p to estimate percolation probability (coarse)
    p_vals = np.linspace(0.4, 0.75, 12)
    perc_probs = percolation_probability_site(L=40, p_vals=p_vals, trials=40, seed=seed)
    plot_percolation_curve(p_vals, perc_probs, label="Site percolation probability (L=40)")

    # Bond percolation example
    p_bond = 0.5  # critical is ~0.5 for bond percolation on square lattice
    bonds_h, bonds_v = generate_bond_configuration(L, p_bond, seed=seed)
    lr_b, tb_b, roots_b, root_map_b, uf_b = has_spanning_cluster_bond(bonds_h, bonds_v, L)
    print(f"Bond percolation (p={p_bond}), percolates: {lr_b or tb_b}")
    # visualize bond occupancy not implemented as image; show basic textual info
    print("Bond percolation demo done.")

    # Optional: animate site growth and save
    print("Generating growth animation (small L=80 to keep file reasonable)...")
    animate_site_growth(L=80, p_steps=60, seed=seed, savefile=None)  # set filename to 'growth.mp4' to save

if __name__ == "__main__":
    demo()
