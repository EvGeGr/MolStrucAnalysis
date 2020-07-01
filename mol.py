import numpy, functions, const

class Mol:

    #Read in molecular Cartesian coordinates from a file
    #First line of the file - number of atoms (n_atoms)
    #Seconf line of the file - comment
    #Each next of n_atoms lines: Atom type (at_type) X-coordinate  Y-coordinate  Z-coordinate
    def read_geom(self, geom_file=None):
        if self.__geom:
            print("Geometry has already been read in!\n")
            return
        self.__geom_file = geom_file
        lines = None
        if geom_file:
            try:
                f = open(self.__geom_file, 'r')
                lines = f.readlines()
                f.close()
            except (FileNotFoundError, IOError):
                print("File {} not found".format(self.__geom_file))
        else:
            print("No geometry file has been provided!")
            self.read_in_geom_err()
            self.__geom = None

        if lines:
            xyz_array = []
            for line in lines:
                xyz_array.append(line.split())

            n_atoms = 0
            at_types = []

            try:
                n_atoms = int(xyz_array[0][0])
            except(TypeError):
                print("A number of atoms (integer) is expected in the first line of the geometry file")
                return

            comment = ' '.join(xyz_array[1][:])

            coords = numpy.zeros((n_atoms, 3))
            for i in range(n_atoms):
                at_types.append(xyz_array[i + 2][0])
                for j in range(3):
                    try:
                        coords[i][j] = float(xyz_array[i + 2][j + 1])
                    except(TypeError):
                        print("Error in reading in coordinates:\n"
                              "x, y, z coordinates are expected for each atom after the atom symbol!\n"
                              "A coordinate is assigned to '0' if you see this error")
            self.__geom = [comment, at_types, coords]

    # build graph of covalently and hydrogen bonded atoms (bond_graph)
    def get_bond_graph(self):
        if self.__geom:
            at_types, coords = self.__geom[1:3]
            n_atoms = len(at_types)
            self.__bond_graph = [[] for i in range(n_atoms)]
            for i in range(n_atoms):
                covrad1 = const.cov_rads[at_types[i]]
                for j in range(i + 1, n_atoms):
                    covrad2 = const.cov_rads[at_types[j]]
                    thresh = const.bond_thresh * (covrad1 + covrad2)
                    r12 = functions.get_r12(coords[i], coords[j])
                    if (r12 < thresh):
                        self.__bond_graph[i].append(j)
                        self.__bond_graph[j].append(i)

            # check for H-bonds and if there are H-bonds building a graph for them (hbond_graph)
            # add H-bonds to bond_graph
            self.__hbond_graph = [[] for i in range(n_atoms)]
            for i in range(n_atoms):
                if at_types[i] in const.hB_atoms:
                    for j in filter(lambda x: x != i, range(n_atoms)):
                        if at_types[j] == 'H' and i not in self.__bond_graph[j]:
                            if any(hB_atom in list(at_types[x] for x in self.__bond_graph[j]) for hB_atom in
                                   const.hB_atoms):
                                r12 = functions.get_r12(coords[i], coords[j])
                                if r12 < const.hB_thresh:
                                    self.__nhb+=1
                                    self.__hbond_graph[i].append(j)
                                    self.__hbond_graph[j].append(i)
                                    self.__bond_graph[i].append(j)
                                    self.__bond_graph[j].append(i)
        else:
            self.__read_in_geom_err()

    # determine atoms which are covalently and hydrogen bonded from bond graph
    def get_bonds(self):
        self.bonds=None
        if self.__bond_graph:
            at_types, coords = self.__geom[1:3]
            n_atoms = len(at_types)
            self.bonds = []
            for i in range(n_atoms):
                for a in range(len(self.__bond_graph[i])):
                    j = self.__bond_graph[i][a]
                    if (i < j):
                        r12 = functions.get_r12(coords[i], coords[j])
                        self.bonds.append([i, j, r12])
        else:
            self.read_in_geom_err()

    # determine atoms which form a bond angle from bond graph
    def get_angles(self):
        self.angles=None
        if self.__bond_graph:
            at_types, coords = self.__geom[1:3]
            n_atoms = len(at_types)
            self.angles = []
            for j in range(n_atoms):
                n_jbonds = len(self.__bond_graph[j])
                for a in range(n_jbonds):
                    i = self.__bond_graph[j][a]
                    for b in range(a + 1, n_jbonds):
                        k = self.__bond_graph[j][b]
                        a123 = functions.get_a123(coords[i], coords[j], coords[k])
                        self.angles.append([i, j, k, a123])
            self.angles = sorted(self.angles, key=lambda angle: angle[0])
            self.angle_names=[]
            self.angle_values=[]
            for a in self.angles:
                self.angle_names.append('%i-%i-%i' % (a[0] + 1, a[1] + 1, a[2] + 1))
                self.angle_values.append(a[3])
        else:
            self.read_in_geom_err()

    # determine atoms which form torsion angles from bond graph
    def get_torsions(self):
        self.torsions=None
        if self.__bond_graph:
            print("XXX")
            at_types, coords = self.__geom[1:3]
            n_atoms = len(at_types)
            self.torsions = []
            for j in range(n_atoms):
                n_jbonds = len(self.__bond_graph[j])
                for a in range(n_jbonds):
                    k = self.__bond_graph[j][a]
                    if (k < j):
                        continue
                    n_kbonds = len(self.__bond_graph[k])
                    for b in range(n_jbonds):
                        i = self.__bond_graph[j][b]
                        if (i == k):
                            continue
                        for c in range(n_kbonds):
                            l = self.__bond_graph[k][c]
                            if (l == j or l == i):
                                continue
                            t1234 = functions.get_t1234(coords[i], coords[j], coords[k], coords[l])
                            self.torsions.append([i, j, k, l, t1234])
            self.torsions = sorted(self.torsions, key=lambda torsion: torsion[0])
        else:
            self.__no_bond_graph_err()

    # determine atoms which form out-of-plane angles from bond graph
    def get_outofplanes(self):
        self.outofplanes = None
        if self.__bond_graph:
            at_types, coords = self.__geom[1:3]
            n_atoms = len(at_types)
            self.outofplanes = []
            for l in range(n_atoms):
                n_lbonds = len(self.__bond_graph[l])
                for a in range(n_lbonds):
                    i = self.__bond_graph[l][a]
                    for b in range(n_lbonds):
                        j = self.__bond_graph[l][b]
                        if (i == j):
                            continue
                        for c in range(b + 1, n_lbonds):
                            k = self.__bond_graph[l][c]
                            if (i == k):
                                continue
                            o1234 = functions.get_o1234(coords[i], coords[j], coords[k], coords[l])
                            self.outofplanes.append([i, j, k, l, o1234])
            self.outofplanes = sorted(self.outofplanes, key=lambda outofplane: outofplane[0])
        else:
            self.__no_bond_graph_err()


    #Processing error "no geomtry has been read in"
    def __read_in_geom_err(self):
        print("No geometry has been read in!\n\rUse '.read_geom(file)' method to read a geometry from a file!\n")

    #Processing error "no bond graph has been constructed"
    def __no_bond_graph_err(self):
        print("Bond graph has not been constructed\n\rUse '.get_bond_graph()' method to construct a bond graph!\n")

    # Processing error "no geometry file has been specified"
    def __no_geom_file(self):
        print("No geometry file has been given!\n\rUse '.read_geom(file)' method to read a geometry from a file!\n")

    #Initialize class Mol ...
    def __init__(self,geom_file=None):
        self.geom_file=geom_file
        self.__geom=None
        self.__bond_graph=None
        self.__hbond_graph=None
        self.__nhb=0
        self.bonds=None
        self.angles=None
        self.torsions=None
        self.outofplanes=None
        if self.geom_file:
            self.read_geom(self.geom_file)
        else:
            self.__no_geom_file()
        if self.__geom:
            self.get_bond_graph()
        else:
            self.__no_bond_graph_err()
        if self.__bond_graph:
            self.get_bond_graph()
            self.get_bonds()
            self.get_angles()
            self.get_torsions()
            self.get_outofplanes()
        else:
            self.__no_bond_graph_err()

    # print coordinates
    def print_geom(self):
        if self.__geom:
            comment, at_types, coords = self.__geom[0:3]
            n_atoms = len(at_types)
            print("Molecular geometry '{}':".format(comment))
            print("Number of atoms - {}".format(n_atoms))
            for i in range(n_atoms):
                print('%-2s' % (at_types[i]), end='')
                for j in range(3):
                    print(' %12.6f' % (coords[i][j]), end='')
                print('\n', end='')
            print('\n', end='')
        else:
            self.__read_in_geom_err()

    # print bond graph
    def print_bond_graph(self):
        if self.__bond_graph:
            comment, at_types, coords = self.__geom[0:3]
            n_atoms = len(at_types)
            print("Bond graph '{}':".format(comment))
            for i in range(n_atoms):
                print("{:<2} {:<2}- ".format(i + 1, at_types[i]),end='')
                for j in range(len(self.__bond_graph[i])):
                    print("{}".format(self.__bond_graph[i][j] + 1),end='')
                    if j < len(self.__bond_graph[i])-1:
                        print(', ',end='')
                    else:
                        print('\n', end='')
            print('\n', end='')

        else:
            self.__no_bond_graph_err()
    #print hbond_graph
    def print_hbond_graph(self):
        if self.__hbond_graph:
            if self.__nhb == 0:
                print("No H-bonds have been identified!\n")
                return
            comment, at_types, coords = self.__geom[0:3]
            n_atoms = len(at_types)
            print("H-Bonds graph:")
            for i in range(n_atoms):
                for j in range(len(self.__hbond_graph[i])):
                    if len(self.__hbond_graph[i] != 0):
                        print("{:<2} {:<2}- ".format(i + 1, at_types[i]), end='')
                        print("{}".format(self.__hbond_graph[i][j] + 1), end='')
                        if j < len(self.__hbond_graph[i]) - 1:
                            print(', ', end='')
            else:
                print()
        else:
            self.__no_bond_graph_err()

    # print list of bond lengths
    def print_bonds(self):
        if self.bonds:
            at_types = self.__geom[1]
            n_bonds = len(self.bonds)
            print('%i bond length(s) found (Angstrom)' % (n_bonds))
            if (n_bonds > 0):
                print(' atoms            elements         values')
            for q in range(n_bonds):
                n1, n2 = self.bonds[q][0:2]
                r12 = self.bonds[q][2]
                nstr = '%i-%i' % (n1 + 1, n2 + 1)
                tstr = '(%s-%s) ' % (at_types[n1], at_types[n2])
                print(' %-15s  %-13s    %6.4f\n' % (nstr, tstr, r12), end='')
            print('\n', end='')
        else:
            print("There are no bonds for this molecule or they have not (yet) been determined")

    # print list of bond angles
    def print_angles(self):
        if self.angles:
            at_types = self.__geom[1]
            n_angles = len(self.angles)
            print('%i bond angle(s) found (degrees)' % (n_angles))
            if (n_angles > 0):
                print(' atoms            elements         values')
            for q in range(n_angles):
                n1, n2, n3 = self.angles[q][0:3]
                a123 = self.angles[q][3]
                nstr = '%i-%i-%i' % (n1 + 1, n2 + 1, n3 + 1)
                tstr = '(%s-%s-%s) ' % (at_types[n1], at_types[n2], at_types[n3])
                print(' %-15s  %-13s   %7.3f\n' % (nstr, tstr, a123), end='')
            print('\n', end='')
        else:
            print("There are no angles for this molecule or they have not (yet) been determined")

    # print list of torsion angles
    def print_torsions(self):
        if self.torsions==0:
            at_types = self.__geom[1]
            n_torsions = len(self.torsions)
            print('%i torsion angle(s) found (degrees)' % (n_torsions))
            if (n_torsions > 0):
                print(' atoms            elements         values')
            for q in range(n_torsions):
                n1, n2, n3, n4 = self.torsions[q][0:4]
                t1234 = self.torsions[q][4]
                nstr = '%i-%i-%i-%i' % (n1 + 1, n2 + 1, n3 + 1, n4 + 1)
                tstr = '(%s-%s-%s-%s) ' % (at_types[n1], at_types[n2], at_types[n3], at_types[n4])
                print(' %-15s  %-13s  %8.3f\n' % (nstr, tstr, t1234), end='')
            print('\n', end='')
        else:
            print("There are no dihedrals for this molecule or they have not (yet) been determined")

    # print list of out-of-plane angles to screen
    def print_outofplanes(self):
        if self.outofplanes==None:
            at_types = self.__geom[1]
            n_outofplanes = len(self.outofplanes)
            print('%i out-of-plane angle(s) found (degrees)' % (n_outofplanes))
            if (n_outofplanes > 0):
                print(' atoms            elements         values')
            for q in range(n_outofplanes):
                n1, n2, n3, n4 = self.outofplanes[q][0:4]
                o1234 = self.outofplanes[q][4]
                nstr = '%i-%i-%i-%i' % (n1 + 1, n2 + 1, n3 + 1, n4 + 1)
                tstr = '(%s-%s-%s-%s) ' % (at_types[n1], at_types[n2], at_types[n3], at_types[n4])
                print(' %-15s  %-13s  %8.3f\n' % (nstr, tstr, o1234), end='')
            print('\n', end='')
        else:
            print("There are no out-of-plane angles for this molecule or they have not (yet) been determined")
## Further