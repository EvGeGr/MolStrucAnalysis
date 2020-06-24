import numpy, functions, const

class Mol:

    #Reading in molecular Cartesians from a file
    #First line of the file - number of atoms (n_atoms)
    #Seconf line of the file - comment
    #Each next of n_atoms lines: Atom type (at_type) X-coordinate  Y-coordinate  Z-coordinate
    def read_geom(self, geom_file=None):
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
            print("No geometry file is given!\n"
                  "Geometry is empty!")
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

    # build graph of which atoms are covalently bonded
    def get_bond_graph(self,geom):
        self.__bond_graph=None
        comment, at_types, coords = geom[0:3]
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
        # check for H-bonds
        nhb = 0
        for i in range(n_atoms):
            if at_types[i] in const.hB_atoms:
                for j in filter(lambda x: x != i, range(n_atoms)):
                    if at_types[j] == 'H' and i not in self.__bond_graph[j]:
                        if any(hB_atom in list(at_types[x] for x in self.__bond_graph[j]) for hB_atom in const.hB_atoms):
                            r12 = functions.get_r12(coords[i], coords[j])
                            if r12 < const.hB_thresh:
                                nhb += 1
                                self.__bond_graph[i].append(j)
                                self.__bond_graph[j].append(i)
#        return bond_graph, nhb

    def __init__(self,geom_file=None):
        self.geom_file=geom_file
        self.__geom=None
        self.__bond_graph=None
        if self.geom_file:
            self.read_geom(self.geom_file)
        if self.__geom:
            self.get_bond_graph(self.__geom)

    # print coordinates to screen
    def print_geom(self):
        if self.__geom:
            comment, at_types, coords = self.__geom[0:3]
            n_atoms = len(at_types)
            print("Molecular geometry:\t'{}'".format(comment))
            print("Number of atoms:\t{}".format(n_atoms))
            for i in range(n_atoms):
                print('%-2s' % (at_types[i]), end='')
                for j in range(3):
                    print(' %12.6f' % (coords[i][j]), end='')
                print('\n', end='')
            print('\n', end='')
        else:
            print("No geometry has been read in")

    # print bond graph to screen
    def print_bond_graph(self):
        if self.__bond_graph:
            comment, at_types, coords = self.__geom[0:3]
            n_atoms = len(at_types)
            print("Bond graph:\t'{}'".format(comment))
            for i in range(n_atoms):
                print("{:<2} {:<2}- ".format(i + 1, at_types[i]),end='')
                for j in range(len(self.__bond_graph[i])):
                    print("{}".format(self.__bond_graph[i][j] + 1),end='')
                    if j < len(self.__bond_graph[i])-1:
                        print(', ',end='')
                    else:
                        print()
        else:
            print("Bond graph has not yet been created!\n\rUse '{}.get_bond_graph' to creat it!".format
                  (str(self.__class__).split('.')[0].split("'")[1]))


## Further