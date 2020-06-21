import numpy

class Mol:

    #Reading in molecular Cartesians from a file
    #First line of the file - number of atoms (n_atoms)
    #Seconf line of the file - comment
    #Each next of n_atoms lines: Atom type (at_type) X-coordinate  Y-coordinate  Z-coordinate
    def read_geom(self, geom_file=None):
        self.geom_file = geom_file
        if self.geom_file:
            try:
                f = open(self.geom_file, 'r')
                lines = f.readlines()
                f.close()
            except (FileNotFoundError, IOError):
                print("Wrong file or file path")
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

    def __init__(self,geom_file=None):
        self.geom_file=geom_file
        self.__geom=None
        if self.geom_file:
            self.read_geom(self.geom_file)

    # print coordinates to screen
    def print_geom(self):
        if self.__geom:
            comment, at_types, coords = self.__geom[0:3]
            n_atoms = len(at_types)
            print('%i\n%s\n' % (n_atoms, comment), end='')
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
        if self.__geom:
            comment, at_types, coords = self.__geom[0:3]
            n_atoms = len(at_types)
            print('%s\n' % (comment), end='')
            for i in range(n_atoms):
                print(' %4i %-2s -' % (i + 1, at_types[i]), end='')
                for j in range(len(bond_graph[i])):
                    print(' %i' % (bond_graph[i][j] + 1), end='')
                print('\n', end='')
            print('\n', end='')
        else:
            print("No geometry has been read in")

## Further