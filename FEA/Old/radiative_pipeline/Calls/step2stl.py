"""step2stl.py - tessellate a simple B-spline STEP solid into a binary STL.

No CAD kernel needed (numpy + scipy only). Handles what the waverider
generator writes: ADVANCED_FACEs on B_SPLINE_SURFACE_WITH_KNOTS (untrimmed,
sampled over their full parameter domain) and PLANE faces (triangulated from
their boundary edges). Checks the result is watertight and reports it.

usage:  python step2stl.py in.step out.stl [target_edge_mm]
"""
import re
import sys

import numpy as np
from scipy.interpolate import BSpline
from scipy.spatial import cKDTree


# ---------------------------------------------------------------- parsing
def parse_step(path):
    txt = open(path, encoding='latin-1').read()
    data = txt[txt.index('DATA;') + 5:txt.index('ENDSEC;', txt.index('DATA;'))]
    ents = {}
    for m in re.finditer(r'#(\d+)\s*=\s*(.*?);\s*(?=#\d+\s*=|$)', data, re.S):
        body = ' '.join(m.group(2).split())
        km = re.match(r'([A-Z0-9_]+)\s*\((.*)\)$', body)
        if km:
            ents[int(m.group(1))] = (km.group(1), parse_args(km.group(2)))
        else:
            ents[int(m.group(1))] = ('COMPLEX', body)
    return ents


def parse_args(s):
    """Parse a STEP argument list into nested python lists."""
    out, stack, tok, i = [], [], '', 0
    cur = out
    while i < len(s):
        c = s[i]
        if c == "'":
            j = s.index("'", i + 1)
            tok += s[i:j + 1]
            i = j + 1
            continue
        if c == '(':
            new = []
            cur.append(new)
            stack.append(cur)
            cur = new
        elif c == ')':
            if tok.strip():
                cur.append(conv(tok.strip()))
            tok = ''
            cur = stack.pop()
        elif c == ',':
            if tok.strip():
                cur.append(conv(tok.strip()))
            tok = ''
        else:
            tok += c
        i += 1
    if tok.strip():
        cur.append(conv(tok.strip()))
    return out


def conv(t):
    if t.startswith('#'):
        return ('ref', int(t[1:]))
    try:
        return float(t)
    except ValueError:
        return t


def pt(ents, r):
    return np.array(ents[r[1]][1][1], float)


# ---------------------------------------------------------------- geometry
def knot_vector(mults, knots):
    return np.repeat(np.array(knots, float), np.array(mults, int).astype(int))


def bspline_surface(ents, sid):
    a = ents[sid][1]
    p, q = int(a[1]), int(a[2])
    P = np.array([[pt(ents, r) for r in row] for row in a[3]])      # (nu, nv, 3)
    U = knot_vector(a[8], a[10])
    V = knot_vector(a[9], a[11])
    return p, q, P, U, V


def eval_surface(p, q, P, U, V, us, vs):
    nu, nv, _ = P.shape
    Bu = BSpline.design_matrix(us, U, p).toarray()                    # (len(us), nu)
    Bv = BSpline.design_matrix(vs, V, q).toarray()
    return np.einsum('ai,ijk,bj->abk', Bu, P, Bv)                    # (len(us), len(vs), 3)


def eval_edge(ents, eid, n):
    """Sample an EDGE_CURVE from its start vertex to its end vertex."""
    a = ents[eid][1]
    v0, v1 = pt(ents, ents[a[1][1]][1][1]), pt(ents, ents[a[2][1]][1][1])
    crv = ents[a[3][1]]
    if crv[0] == 'SURFACE_CURVE':
        crv = ents[crv[1][1][1]]
    if crv[0] == 'LINE':
        return np.linspace(v0, v1, 2)
    if crv[0] == 'B_SPLINE_CURVE_WITH_KNOTS':
        c = crv[1]
        deg = int(c[1])
        C = np.array([pt(ents, r) for r in c[2]])
        T = knot_vector(c[6], c[7])
        ts = np.linspace(T[deg], T[-deg - 1], n)
        xs = BSpline.design_matrix(ts, T, deg).toarray() @ C
        if np.linalg.norm(xs[0] - v0) > np.linalg.norm(xs[-1] - v0):
            xs = xs[::-1]
        return xs
    raise NotImplementedError(f'edge curve type {crv[0]}')


def grid_counts(X, h):
    du = np.linalg.norm(np.diff(X, axis=0), axis=2).sum(axis=0).max()
    dv = np.linalg.norm(np.diff(X, axis=1), axis=2).sum(axis=1).max()
    return max(int(np.ceil(du / h)) + 1, 9), max(int(np.ceil(dv / h)) + 1, 9)   # >=9 across thin leading-edge strips


def tess_bspline(ents, sid, h):
    p, q, P, U, V = bspline_surface(ents, sid)
    u0, u1, v0, v1 = U[p], U[-p - 1], V[q], V[-q - 1]
    coarse = eval_surface(p, q, P, U, V, np.linspace(u0, u1, 41), np.linspace(v0, v1, 41))
    nu, nv = grid_counts(coarse, h)
    # cluster samples toward the patch edges (leading edge lives there)
    s = lambda n: 0.5 - 0.5 * np.cos(np.linspace(0, np.pi, n))
    X = eval_surface(p, q, P, U, V, u0 + (u1 - u0) * s(nu), v0 + (v1 - v0) * s(nv))
    idx = np.arange(nu * nv).reshape(nu, nv)
    a, b, c, d = idx[:-1, :-1], idx[1:, :-1], idx[1:, 1:], idx[:-1, 1:]
    tris = np.vstack([np.stack([a, b, c], -1).reshape(-1, 3), np.stack([a, c, d], -1).reshape(-1, 3)])
    return X.reshape(-1, 3), tris


def tess_plane(ents, face, h, other_boundary):
    """Triangulate a planar face from its boundary loop (ear clipping in 2D)."""
    a = ents[face][1]
    loop = ents[ents[a[1][0][1]][1][1][1]][1][1]
    pts = []
    for oe in loop:
        e = ents[oe[1]][1]
        xs = eval_edge(ents, e[3][1], 400)
        if e[4] == '.F.':
            xs = xs[::-1]
        if xs.shape[0] == 2:                                          # straight edge: subdivide
            n = max(int(np.linalg.norm(xs[1] - xs[0]) / h), 1) + 1
            xs = np.linspace(xs[0], xs[1], n)
        pts.append(xs[:-1])
    poly = np.vstack(pts)
    # snap to the neighbouring surfaces' boundary vertices where they coincide, for watertightness
    if other_boundary is not None and len(other_boundary):
        d, j = cKDTree(other_boundary).query(poly)
        # resample the loop at the neighbours' boundary points that lie on this plane
        pl = ents[a[2][1]][1]
        ax = ents[pl[1][1]][1]
        o, n = pt(ents, ax[1]), np.array(ents[ax[2][1]][1][1], float)
        on = other_boundary[np.abs((other_boundary - o) @ n) < 1e-6 * max(1, np.abs(o).max())]
        if len(on) > 3:
            poly = order_on_loop(poly, on)
    tris = ear_clip(poly)
    # orient along the face's outward normal (plane normal, flipped if face sense is .F.)
    pl = ents[ents[a[2][1]][1][1][1]][1]
    nout = np.array(ents[pl[2][1]][1][1], float) * (1 if a[3] == '.T.' else -1)
    tn = np.cross(poly[tris[:, 1]] - poly[tris[:, 0]], poly[tris[:, 2]] - poly[tris[:, 0]])
    if (tn @ nout).sum() < 0:
        tris = tris[:, ::-1]
    return poly, tris


def order_on_loop(loop, cand):
    """Order candidate points by arc-length position along a closed polyline."""
    seg = np.diff(np.vstack([loop, loop[:1]]), axis=0)
    cum = np.concatenate([[0], np.cumsum(np.linalg.norm(seg, axis=1))])
    pos = []
    for c in cand:
        t = np.clip(((c - loop) * seg).sum(1) / (seg * seg).sum(1).clip(1e-30), 0, 1)
        dist = np.linalg.norm(loop + t[:, None] * seg - c, axis=1)
        k = dist.argmin()
        pos.append(cum[k] + t[k] * np.linalg.norm(seg[k]))
    cand = cand[np.argsort(pos)]
    keep = np.r_[True, np.linalg.norm(np.diff(cand, axis=0), axis=1) > 1e-9]
    return cand[keep]


def ear_clip(poly3):
    c = poly3.mean(0)
    _, _, vt = np.linalg.svd(poly3 - c)
    xy = (poly3 - c) @ vt[:2].T
    idx = list(range(len(xy)))
    area = 0.5 * np.sum(xy[:, 0] * np.roll(xy[:, 1], -1) - np.roll(xy[:, 0], -1) * xy[:, 1])
    if area < 0:
        idx.reverse()
    tris, guard = [], 0
    while len(idx) > 3 and guard < 10 * len(xy) ** 2:
        guard += 1
        n = len(idx)
        best = None
        for k in range(n):
            i0, i1, i2 = idx[k - 1], idx[k], idx[(k + 1) % n]
            A, B, C = xy[i0], xy[i1], xy[i2]
            cr = (B[0] - A[0]) * (C[1] - A[1]) - (B[1] - A[1]) * (C[0] - A[0])
            if cr <= 1e-14:
                continue
            others = xy[[j for j in idx if j not in (i0, i1, i2)]]
            if len(others) and inside(others, A, B, C).any():
                continue
            # prefer the fattest ear
            e = np.array([np.linalg.norm(B - A), np.linalg.norm(C - B), np.linalg.norm(A - C)])
            qual = cr / (e ** 2).sum()
            if best is None or qual > best[0]:
                best = (qual, k)
        if best is None:
            break
        k = best[1]
        tris.append([idx[k - 1], idx[k], idx[(k + 1) % n]])
        idx.pop(k)
    if len(idx) == 3:
        tris.append(idx)
    return np.array(tris)


def inside(P, A, B, C):
    def s(p1, p2, p3):
        return (p1[:, 0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[:, 1] - p3[1])
    d1, d2, d3 = s(P, A, B), s(P, B, C), s(P, C, A)
    return ~(((d1 < 0) | (d2 < 0) | (d3 < 0)) & ((d1 > 0) | (d2 > 0) | (d3 > 0)))


# ---------------------------------------------------------------- mesh utils
def weld(V, F, tol):
    key = np.round(V / tol).astype(np.int64)
    _, first, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    V2, F2 = V[first], inv.ravel()[F]
    good = (F2[:, 0] != F2[:, 1]) & (F2[:, 1] != F2[:, 2]) & (F2[:, 0] != F2[:, 2])
    return V2, F2[good]


def edge_stats(F):
    e = np.sort(np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]), axis=1)
    _, cnt = np.unique(e, axis=0, return_counts=True)
    return (cnt == 1).sum(), (cnt > 2).sum()


def signed_volume(V, F):
    return np.einsum('ij,ij->i', V[F[:, 0]], np.cross(V[F[:, 1]], V[F[:, 2]])).sum() / 6


def write_stl(path, V, F):
    n = np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]])
    n /= np.linalg.norm(n, axis=1, keepdims=True).clip(1e-30)
    rec = np.zeros(len(F), dtype=[('n', '<f4', 3), ('v', '<f4', (3, 3)), ('a', '<u2')])
    rec['n'], rec['v'] = n, V[F]
    with open(path, 'wb') as f:
        f.write(b'step2stl'.ljust(80, b' '))
        f.write(np.uint32(len(F)).tobytes())
        f.write(rec.tobytes())


# ---------------------------------------------------------------- main
def main(inp, out, h):
    ents = parse_step(inp)
    faces = [k for k, v in ents.items() if v[0] == 'ADVANCED_FACE']
    Vs, Fs, planes, off = [], [], [], 0
    for f in faces:
        a = ents[f][1]
        surf = ents[a[2][1]][0]
        if surf == 'B_SPLINE_SURFACE_WITH_KNOTS':
            V, F = tess_bspline(ents, a[2][1], h)
            if a[3] == '.F.':
                F = F[:, ::-1]
            Vs.append(V); Fs.append(F + off); off += len(V)
            print(f'face #{f}: B-spline, {len(F)} tris')
        elif surf == 'PLANE':
            planes.append(f)
        else:
            raise NotImplementedError(f'face #{f}: surface type {surf}')
    V, F = weld(np.vstack(Vs), np.vstack(Fs), 1e-4)
    # boundary vertices of the B-spline shell -> the planar face must stitch to these
    e = np.sort(np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]), axis=1)
    u, cnt = np.unique(e, axis=0, return_counts=True)
    bverts = V[np.unique(u[cnt == 1])]
    for f in planes:
        pv, pf = tess_plane(ents, f, h, bverts)
        V = np.vstack([V, pv]); F = np.vstack([F, pf + len(V) - len(pv)])
        print(f'face #{f}: plane, {len(pf)} tris')
    V, F = weld(V, F, 1e-4)
    if signed_volume(V, F) < 0:
        F = F[:, ::-1]
    nb, nm = edge_stats(F)
    lo, hi = V.min(0), V.max(0)
    print(f'total {len(F)} triangles, {len(V)} vertices')
    print(f'extent x[{lo[0]:.1f},{hi[0]:.1f}] y[{lo[1]:.1f},{hi[1]:.1f}] z[{lo[2]:.1f},{hi[2]:.1f}] (STEP units)')
    print(f'open edges {nb}, non-manifold edges {nm}  ->  {"WATERTIGHT" if nb == 0 and nm == 0 else "NOT watertight"}')
    print(f'enclosed volume {signed_volume(V, F):.4g} (STEP units^3)')
    write_stl(out, V, F)
    print(f'wrote {out}')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]) if len(sys.argv) > 3 else 5.0)
