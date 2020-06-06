class Mol:
    def __init__(self,geom_file):
        self.geom_file=geom_file
        if geom_file:
            try:
                f=open(self.geom_file, 'r')
                lines = f.readlines()
                f.close()
            except (FileNotFoundError, IOError):
                print("Wrong file or file path")
        else:
            print("No geometry file is given!\n"
                  "Geometry is empty!")
            self.geom=None

        xyz_array = []
        for line in lines:
            xyz_array.append(line.split())

        try:
            n_atoms = int(xyz_array[0][0])
        except(TypeError):
            print("A number of atoms (integer) is expected in the first line of the geometry file")

        at_types=[]
        coords = numpy.zeros((n_atoms, 3))
        for i in range(n_atoms):
            at_types.append(xyz_array[i + 2][0])
            for j in range(3):
                coords[i][j] = float(xyz_array[i + 2][j + 1])
        self.geom = [at_types, coords]




    pass

## Further