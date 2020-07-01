import math, numpy, const

## MATH FUNCTIONS ##

# compare if two values are the same within specified threshold
def are_same(n1, n2, tol, minval):
    same = False
    nmax = max(abs(n1), abs(n2), abs(minval))
    comp = abs((n2 - n1) / nmax)
    if (comp <= 10**(-tol)):
        same = True
    return same

# calculate distance between two 3-d cartesian coordinates
def get_r12(coords1, coords2):
    r2 = 0.0
    for p in range(3):
        r2 += (coords2[p] - coords1[p])**2
    r = math.sqrt(r2)
    return r

# calculate unit vector between to 3-d cartesian coordinates
def get_u12(coords1, coords2):
    r12 = get_r12(coords1, coords2)
    u12 = [0.0 for p in range(3)]
    for p in range(3):
        u12[p] = (coords2[p] - coords1[p]) / r12
    return u12

# calculate dot product between two unit vectors
def get_udp(uvec1, uvec2):
    udp = 0.0
    for p in range(3):
        udp += uvec1[p] * uvec2[p]
    udp = max(min(udp, 1.0), -1.0)
    return udp

# calculate unit cross product between two unit vectors
def get_ucp(uvec1, uvec2):
    ucp = [0.0 for p in range(3)]
    cos_12 = get_udp(uvec1, uvec2)
    sin_12 = math.sqrt(1 - cos_12**2)
    ucp[0] = (uvec1[1]*uvec2[2] - uvec1[2]*uvec2[1]) / sin_12
    ucp[1] = (uvec1[2]*uvec2[0] - uvec1[0]*uvec2[2]) / sin_12
    ucp[2] = (uvec1[0]*uvec2[1] - uvec1[1]*uvec2[0]) / sin_12
    return ucp

# calculate angle between three 3-d cartesian coordinates
def get_a123(coords1, coords2, coords3):
    u21 = get_u12(coords2, coords1)
    u23 = get_u12(coords2, coords3)
    dp2123 = get_udp(u21, u23)
    a123 = const.rad2deg * math.acos(dp2123)
    return a123

# calculate torsion angle between four 3-d cartesian coordinates
def get_t1234(coords1, coords2, coords3, coords4):
    u21 = get_u12(coords2, coords1)
    u23 = get_u12(coords2, coords3)
    u32 = get_u12(coords3, coords2)
    u34 = get_u12(coords3, coords4)
    u21c23 = get_ucp(u21, u23)
    u32c34 = get_ucp(u32, u34)
    dp = get_udp(u21c23, u32c34)
    sign = 2 * float(get_udp(u21c23, u34) < 0) - 1
    t1234 = const.rad2deg * sign * math.acos(dp)
    return t1234

# calculate out-of-plane (improper torsion) angle between four 3-d cartesian coordinates
def get_o1234(coords1, coords2, coords3, coords4):
    u42 = get_u12(coords4, coords2)
    u43 = get_u12(coords4, coords3)
    u41 = get_u12(coords4, coords1)
    u42c43 = get_ucp(u42, u43)
    dp = get_udp(u42c43, u41)
    o1234 = const.rad2deg * math.asin(dp)
    return o1234

# translate coordinates by a defined vector and scale factor
def translate_coords(geom, vector, scale):
    coords = geom[1]
    n_atoms = len(coords)
    for i in range(n_atoms):
        for j in range(3):
            coords[i][j] += scale * vector[j]
    return geom

# calculate center of mass of a set of atoms
def get_com(geom):
    at_types, coords = geom[0:2]
    n_atoms = len(at_types)
    com = [0.0 for p in range(3)]
    mass = 0.0
    for i in range(n_atoms):
        at_mass = const.at_masses[at_types[i]]
        mass += at_mass
        for j in range(3):
            com[j] += at_mass * coords[i][j]
    for p in range(3):
        com[p] /= mass
    return com

# calculate moment of inertia tensor for a set of atoms
def get_moi(geom):
    at_types, coords = geom[0:2]
    n_atoms = len(at_types)
    moi = [[0.0 for q in range(3)] for p in range(3)]
    for i in range(n_atoms):
        at_mass = const.at_masses[at_types[i]]
        for p in range(3):
            for q in range(3):
                if (p == q):
                    r = (p+1) % 3
                    s = (p+2) % 3
                    moi[p][p] += at_mass * (coords[i][r]**2 + coords[i][s]**2)
                else:
                    moi[p][q] += -at_mass * coords[i][p] * coords[i][q]
    moi = numpy.matrix(moi)
    return moi

# calculate principal moments of inertia (eigenvalues of tensor)
def get_prinmom(moi):
    prinmom = numpy.linalg.eigvalsh(moi)
    return prinmom

# calculate rotational frequencies in MHz and wavenumbers (cm^-1)
def get_rotfreq(prinmom):
    freqcm1, freqmhz = [], []
    for p in range(3):
        iszero = are_same(prinmom[p], 0.0, const.mom_thresh, const.mom_min)
        degen = (p>0 and are_same(prinmom[p-1], prinmom[p], const.mom_thresh, const.mom_min))
        if (iszero or degen):
            continue
        mhz = const.h / (8 * math.pi**2 * prinmom[p])
        mhz *= (10**10)**2 * const.na * 10**3 * 10**(-6)
        cm1 = mhz / const.c * 10**6
        freqmhz.append(mhz)
        freqcm1.append(cm1)
    return freqmhz, freqcm1

# rotate molecule to inertial frame of principal moments
def get_inertial_coords(geom, moi):
    moi_eigvals, moi_eigvecs = numpy.linalg.eig(moi)
    coords = geom[1]
    coords = numpy.array(numpy.dot(coords, moi_eigvecs))
    geom[1] = coords
    moi = get_moi(geom)
    order = [0, 1, 2]
    for p in range(3):
        for q in range(p+1, 3):
            if (moi.item(p, p) < moi.item(q, q)):
                temp = order[p]
                order[p] = order[q]
                order[q] = temp
    moveaxes = numpy.zeros((3, 3))
    for p in range(3):
        moveaxes[p][order[p]] = 1.0
    coords = numpy.dot(coords, moveaxes)
    geom[1] = coords
    return geom

## TOPOLOGY FUNCTIONS ##

# build graph of which atoms are covalently bonded
def get_bond_graph(geom):
    at_types, coords = geom[0:2]
    n_atoms = len(at_types)
    bond_graph = [[] for i in range(n_atoms)]
    for i in range(n_atoms):
        covrad1 = const.cov_rads[at_types[i]]
        for j in range(i+1, n_atoms):
            covrad2 = const.cov_rads[at_types[j]]
            thresh = const.bond_thresh * (covrad1 + covrad2)
            r12 = get_r12(coords[i], coords[j])
            if (r12 < thresh):
                bond_graph[i].append(j)
                bond_graph[j].append(i)
# check for H-bonds
    nhb = 0
    for i in range(n_atoms):
        if at_types[i] in const.hB_atoms:
            for j in filter(lambda x: x != i, range(n_atoms)):
                if at_types[j] == 'H' and i not in bond_graph[j]:
                    if any(hB_atom in list(at_types[x] for x in bond_graph[j]) for hB_atom in const.hB_atoms):
                        r12 = get_r12(coords[i], coords[j])
                        if r12 < 2.2:
                            nhb += 1
                            bond_graph[i].append(j)
                            bond_graph[j].append(i)
    return bond_graph, nhb

# determine atoms which are covalently bonded from bond graph
def get_bonds(geom, bond_graph):
    at_types, coords = geom[0:2]
    n_atoms = len(at_types)
    bonds = []
    for i in range(n_atoms):
        for a in range(len(bond_graph[i])):
            j = bond_graph[i][a]
            if (i < j):
                r12 = get_r12(coords[i], coords[j])
                bonds.append([i, j, r12])
    return bonds

# determine atoms which form a bond angle from bond graph
def get_angles(geom, bond_graph):
    at_types, coords = geom[0:2]
    n_atoms = len(at_types)
    angles = []
    for j in range(n_atoms):
        n_jbonds = len(bond_graph[j])
        for a in range(n_jbonds):
            i = bond_graph[j][a]
            for b in range(a+1, n_jbonds):
                k = bond_graph[j][b]
                a123 = get_a123(coords[i], coords[j], coords[k])
                angles.append([i, j, k, a123])
    angles = sorted(angles, key=lambda angle:angle[0])
    return angles

# determine atoms which form torsion angles from bond graph
def get_torsions(geom, bond_graph):
    at_types, coords = geom[0:2]
    n_atoms = len(at_types)
    torsions = []
    for j in range(n_atoms):
        n_jbonds = len(bond_graph[j])
        for a in range(n_jbonds):
            k = bond_graph[j][a]
            if (k < j):
                continue
            n_kbonds = len(bond_graph[k])
            for b in range(n_jbonds):
                i = bond_graph[j][b]
                if (i == k):
                    continue
                for c in range(n_kbonds):
                    l = bond_graph[k][c]
                    if (l == j or l == i):
                        continue
                    t1234 = get_t1234(coords[i], coords[j], coords[k], coords[l])
                    torsions.append([i, j, k, l, t1234])
    torsions = sorted(torsions, key=lambda torsion:torsion[0])
    return torsions

# determine atoms which form out-of-plane angles from bond graph
def get_outofplanes(geom, bond_graph):
    at_types, coords = geom[0:2]
    n_atoms = len(at_types)
    outofplanes = []
    for l in range(n_atoms):
        n_lbonds = len(bond_graph[l])
        for a in range(n_lbonds):
            i = bond_graph[l][a]
            for b in range(n_lbonds):
                j = bond_graph[l][b]
                if (i == j):
                    continue
                for c in range(b+1, n_lbonds):
                    k = bond_graph[l][c]
                    if (i == k):
                        continue
                    o1234 = get_o1234(coords[i], coords[j], coords[k], coords[l])
                    outofplanes.append([i, j, k, l, o1234])
    outofplanes = sorted(outofplanes, key=lambda outofplane:outofplane[0])
    return outofplanes

# determine molecule type based on principal moments of inertia
def get_moltype(geom, pm):
    same12  = are_same(pm[0], pm[1], const.mom_thresh, const.mom_min)
    same13  = are_same(pm[0], pm[2], const.mom_thresh, const.mom_min)
    same23  = are_same(pm[1], pm[2], const.mom_thresh, const.mom_min)
    onezero = are_same(pm[0], 0.0,   const.mom_thresh, const.mom_min)
    allzero = are_same(pm[2], 0.0,   const.mom_thresh, const.mom_min)
    if (allzero):
        moltype = 'monatomic'
    elif (onezero):
        moltype = 'linear'
    elif (same13):
        moltype = 'a spherical top'
    elif (same12 or same23):
        moltype = 'a symmetric top'
    else:
        moltype = 'an asymmetric top'
    return moltype
