from __future__ import annotations

from pathlib import Path

import networkx as nx
import numpy as np

from data_io import DATA_DIR, LAYER_ORDER, SEED


def _randomize_layer_once(A: np.ndarray, seed: int) -> np.ndarray:
    if A.sum() == 0:
        return A.copy()
    G = nx.from_numpy_array(A.T, create_using=nx.DiGraph)
    G.remove_edges_from(nx.selfloop_edges(G))
    if G.number_of_edges() < 4:
        return A.copy()
    H = G.copy()
    edges = max(H.number_of_edges(), 1)
    nswap = max(edges, 10)
    max_tries = max(edges * 10, 200)
    success = False
    for extra in [1, 2, 4, 8]:
        trial = H.copy()
        try:
            nx.algorithms.swap.directed_edge_swap(
                trial,
                nswap=max(1, nswap // extra),
                max_tries=max_tries * extra,
                seed=seed + extra,
            )
            H = trial
            success = True
            break
        except (nx.NetworkXAlgorithmError, nx.NetworkXError):
            continue
    if not success:
        return A.copy()
    R = nx.to_numpy_array(H, nodelist=range(A.shape[0]), dtype=np.uint8).T
    np.fill_diagonal(R, 0)
    if not np.array_equal(A.sum(axis=0), R.sum(axis=0)):
        raise RuntimeError("Out-degree sequence not preserved")
    if not np.array_equal(A.sum(axis=1), R.sum(axis=1)):
        raise RuntimeError("In-degree sequence not preserved")
    return R


def generate_degree_preserving_nulls(
    tensor: np.ndarray,
    n_draws: int = 1000,
    seed: int = SEED,
    prefix: str = "null",
) -> dict[str, np.ndarray]:
    """Degree-preserving null draws per layer, cached on disk.

    ``prefix`` namespaces the cache. Any analysis run on a *subset* of strains
    (e.g. the n=56 taxon sensitivity in p2_13) must pass its own prefix, or it
    would silently load the full-cohort nulls whose node count no longer matches.
    """
    rng = np.random.default_rng(seed)
    # Every layer's seeds are drawn up front, in LAYER_ORDER, before any cache
    # lookup. Consuming the stream inside the loop instead would let a cached
    # layer shift the seeds of the layers after it, so a partially warm cache
    # would silently produce different draws from a cold run.
    layer_seeds = {
        layer: rng.integers(0, 2**31 - 1, size=n_draws) for layer in LAYER_ORDER
    }
    out = {}
    for idx, layer in enumerate(LAYER_ORDER):
        A = tensor[idx].astype(np.uint8)
        layer_path = DATA_DIR / f"{prefix}_{layer}_{n_draws}.npy"
        if layer_path.exists():
            out[layer] = np.load(layer_path, allow_pickle=False)
            continue
        draws = np.zeros((n_draws, A.shape[0], A.shape[1]), dtype=np.uint8)
        for d in range(n_draws):
            draws[d] = _randomize_layer_once(A, int(layer_seeds[layer][d]))
            if (d + 1) % 100 == 0:
                print(f"[null_models] {layer}: generated {d + 1}/{n_draws} nulls", flush=True)
        out[layer] = draws
        np.save(layer_path, draws, allow_pickle=False)
    return out
