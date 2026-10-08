"""
Finite-element mesh geometry and finite-volume transport on the element graph.

Elements are Abaqus eight-node quadrilaterals (CPE8RP): node ids 1-4 are the corners
(counter-clockwise), 5-8 the mid-side nodes. Two elements are neighbours when they share
an edge (two corner nodes); elements that touch only at a corner are not neighbours.

Transport (Methods 5.4.3): cell-centred finite volumes with a two-point flux,

    dc_e/dt = (1/A_e) * sum_{j in N(e)} D * (L_ej / delta_ej) * (c_j - c_e),

with zero-flux outer boundaries. The explicit update is sub-stepped so that every
sub-step satisfies the positivity/stability bound dt_sub <= min_e A_e / (D * sum_j L_ej/delta_ej).
The scheme conserves sum_e A_e c_e exactly.
"""
import math

import numpy as np


class MeshGeometry:
    def __init__(self, nodes, elements):
        self.ids = list(elements.keys())
        self.index = {eid: i for i, eid in enumerate(self.ids)}
        n = len(self.ids)

        corners = np.array([[nodes[nid] for nid in elements[eid][:4]] for eid in self.ids])  # (n, 4, 2)
        x, y = corners[..., 0], corners[..., 1]
        self.area = 0.5 * np.abs(np.sum(x * np.roll(y, -1, axis=1) - np.roll(x, -1, axis=1) * y, axis=1))
        self.centroid = np.array([np.mean([nodes[nid] for nid in elements[eid]], axis=0) for eid in self.ids])
        self.mean_area = float(self.area.mean())
        # area relative to the mean element; capacities and sources are expressed per mean-sized element
        self.area_ratio = self.area / self.mean_area

        edge_owner = {}
        pairs, lengths = [], []
        for i, eid in enumerate(self.ids):
            c = elements[eid][:4]
            for k in range(4):
                a, b = c[k], c[(k + 1) % 4]
                key = (min(a, b), max(a, b))
                if key in edge_owner:
                    j = edge_owner.pop(key)
                    pairs.append((j, i))
                    lengths.append(math.dist(nodes[a], nodes[b]))
                else:
                    edge_owner[key] = i
        self.pairs = np.array(pairs, dtype=int)                 # interior edges (i, j)
        self.edge_length = np.array(lengths)
        self.delta = np.linalg.norm(self.centroid[self.pairs[:, 0]] - self.centroid[self.pairs[:, 1]], axis=1)
        self.n_boundary_edges = len(edge_owner)

        self.neighbors = {eid: [] for eid in self.ids}
        for i, j in self.pairs:
            self.neighbors[self.ids[i]].append(self.ids[j])
            self.neighbors[self.ids[j]].append(self.ids[i])

        # sum of conductances per element (stability bound)
        g = self.edge_length / self.delta
        self._g = g
        self._g_sum = np.zeros(n)
        np.add.at(self._g_sum, self.pairs[:, 0], g)
        np.add.at(self._g_sum, self.pairs[:, 1], g)

    def diffuse(self, c, D, dt):
        """Advance concentration array c (n,) by dt with diffusion coefficient D; returns new array."""
        if D <= 0.0:
            return c
        dt_max = float(np.min(self.area / (D * self._g_sum)))
        n_sub = max(1, math.ceil(dt / (0.9 * dt_max)))
        h = dt / n_sub
        i, j = self.pairs[:, 0], self.pairs[:, 1]
        coef = D * self._g
        c = c.copy()
        for _ in range(n_sub):
            flux = coef * (c[j] - c[i])            # amount per unit time from j into i
            dc = np.zeros_like(c)
            np.add.at(dc, i, flux)
            np.add.at(dc, j, -flux)
            c += h * dc / self.area
        return c
