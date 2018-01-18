class Matrix:
    def __init__(self, m, n, rws=[]):
        self.m = m
        self.n = n
        self.rrs = rws
        if not self.rrs: self.rrs = [[int(input()) for _ in range(self.n)] for _ in range(self.m)]
        self.cols = [[i[j] for i in self.rrs] for j in range(self.n)]
        assert len(self.rows) == self.m and len(self.columns) == self.n

    @property
    def rows(self):
        if len(self.cols) < len(self.rrs) and self.m == self.n:
            for i, j in zip(range(len(self.rrs)), range(len(self.cols))):
                if self.rrs[i][j] != self.cols[j][i]:
                    for a in self.rrs:
                        a.remove(a[j])
                    self.n -= 1
        return self.rrs
    @property
    def columns(self):
        if len(self.cols) > len(self.rrs) and self.m == self.n:
            for i, j in zip(range(len(self.rrs)), range(len(self.cols))):
                if self.rrs[i][j] != self.cols[j][i]:
                    for a in self.cols:
                        a.remove(a[j])
                    self.m -= 1
        return self.cols

    # Operations and Properties
    @property
    def det(self):
        assert self.m == self.n, "Only square matrices have a well-defined determinant."
        if self.m == 2:
            return (self.rrs[0][0]*self.rrs[1][1]) - (self.rrs[0][1]*self.rrs[1][0])
        elif self.m == 1:
            return self.rrs[0][0]
        def deter_help(self, total):
            for i in range(self.m):
                rest = self.rows[1:]
                helper_mat = Matrix(self.m - 1, self.n, rest)
                helper_mat.columns.remove(helper_mat.columns[i])
                sub_mat = Matrix(self.m - 1, self.n - 1, helper_mat.cols)
                if i % 2 == 0:
                    total += self.rrs[0][i] * sub_mat.det
                else:
                    total -= self.rrs[0][i] * sub_mat.det
            return total
        return deter_help(self, 0)
    def transpose(self):
        return Matrix(self.n, self.m, self.columns)
    def inverse(self):
        assert self.m == self.n, "Only square matrices can be inverted."
        assert self.is_invertible(), "This matrix is not invertible."
        if self.m == 2:
            changed = [[self.rows[1][1], -1 * self.rows[0][1]], \
                        [-1 * self.rows[1][0], self.rows[0][0]]]
            new_matrix = Matrix(2, 2, changed)
            return (1 / self.det) * new_matrix
        def minors(self, entries=[]):
            for i in range(self.m):
                row = []
                copy = self.rows[:]
                copy.remove(self.rows[i])
                rest = copy
                for j in range(self.n):
                    helper_mat = Matrix(self.m - 1, self.n, rest)
                    helper_mat.cols.remove(helper_mat.cols[j])
                    sub_mat = Matrix(self.m - 1, self.n - 1, helper_mat.columns)
                    row.append(sub_mat.det)
                entries.append(row)
            return Matrix(self.m, self.n, entries)
        def cofactors(self, entries=[]):
            minor = minors(self)
            for i in range(self.m):
                row = []
                for j in range(self.n):
                    ent = minor.rows[i][j]
                    if (i + j) % 2 != 0:
                        ent = -1 * ent
                    row.append(ent)
                entries.append(row)
            return 1/self.det * Matrix(self.m, self.n, entries).transpose()
        return cofactors(self)
    def is_invertible(self):
        return self.det != 0
    def is_symmetric(self):
        return self.rows == self.columns

    # Arithmetic and Display methods
    def __add__(self, B):
        assert self.m == B.m and self.n == B.n, "Matrices must be the same size to be added."
        added = [[self.cols[i][j] + B.cols[i][j] for i in range(self.m)] for j in range(self.n)]
        return Matrix(self.m, self.n, added)
    def __sub__(self, B):
        assert self.m == B.m and self.n == B.n, "Matrices must be the same size to be subtracted."
        subbed = [[self.cols[i][j] - B.cols[i][j] for i in range(self.m)] for j in range(self.n)]
        return Matrix(self.m, self.n, subbed)
    def __mul__(self, B):
        assert self.n == B.m, "These matrices are not compatible for multiplication."
        ents = [[sum([self.rows[k][i] * B.columns[j][i] for i in range(self.n)]) \
                    for j in range(B.n)] for k in range(self.m)]
        rounded = []
        for row in ents:
            z = []
            for entry in row:
                if abs(entry - (round(entry))) != 0 and abs(entry - (round(entry))) < 0.00000001:
                    entry = round(entry)
                z.append(float(entry))
            rounded.append(z)
        return Matrix(self.m, B.n, rounded)
    def __rmul__(self, x):
        rounded = []
        scaled = [[self.rows[i][j] * x for j in range(self.m)] for i in range(self.n)]
        for row in scaled:
            z = []
            for entry in row:
                if abs(entry - (round(entry))) != 0 and abs(entry - (round(entry))) < 0.00000001:
                    entry = round(entry)
                z.append(float(entry))
            rounded.append(z)

        return Matrix(self.m, self.n, rounded)
    def __str__(self):
        return "{}".format(self.rows)
    def __repr__(self):
        for i in self.rrs:
            print("| ", end=" ")
            for entry in i:
                print(entry, " ", end =" ")
            print("|")
        return ""

class Vector(Matrix):
    def __init__(self, m, n=1):
        self.m = m
        self.n = n
        assert self.n == 1, "Column vectors have only one dimension."
        return Matrix.__init__(self, self.m, 1)
    def __mul__(self, v):
        return sum([self.rows[i][0] * v.rows[i][0] for i in range(self.m)])
    def __rmul__(self, x):
        scaled = [[self.rows[i][0] * x] for i in range(self.m)]
        return Matrix(self.m, 1, scaled)

# Test Matrices and Vectors
three = Matrix(3, 3, [[1, 2, 3], [4, 5, 6], [7, 8, 9]]) # 0
identity = Matrix(3, 3, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]) # 1
up = Matrix(3, 3, [[1, 2, 3], [0, 9, 5], [0, 0, 2]]) # 18
low = Matrix(3, 3, [[9, 0, 0], [-6, 2, 5], [1, -2, -3]]) # 36
four = Matrix(4, 4, [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]) # 1
a = Matrix(3, 3, [[0, 0, 0], [1, 2, 3], [5, 2, 6]]) # 0
b = Matrix(3, 3, [[1, 0, 4], [3, 0, 0], [5, 5, 2]]) # 60
c = Matrix(4, 4, [[1, 2, 3, 4], [5, 6, 7, 8], [3, 9, 3, 6], [4, 8, 1, 5]]) # 132
d = Matrix(1, 1, [[5]]) # 5
e = Matrix(5, 5, [[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], \
                    [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]]) # 1
f = Matrix(5, 5, [[1, 2, 3, 4, 5], [2, 7, 1, 8, 2], [3, 1, 4, 1, 5], \
                    [0, 1, 1, 2, 3], [5, 8, 5, 3, 2]]) # -280
g = Matrix(2, 3, [[1, 2, 3], [4, 5, 6]])
h = Matrix(3, 4, [[2, 7, 1, 8], [2, 8, 1, 8], [1, 7, 2, 9]])
