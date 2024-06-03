from Vec import Vec
from Poly import Poly, Monome
        
class Line(Vec):
    """Classe pour les lignes des matrices en tant que vecteur.
       Les méthodes modifiant le vecteur vont modifier la matrice.
       Vous ne devriez pas avoir besoin de créer d'objet de cette classe,
       cette classe est utilisée en interne par `Mat`."""
  
    def __init__(self,mat,i):
        self.i   = i
        self.mat = mat
        self.dim = self.mat.dim2
    
    def __getitem__(self,j):
        return self.mat.__coefs__[self.i][j]

    def __setitem__(self,j,x):
        self.mat.__coefs__[self.i][j] = x

class Col(Vec):
    """Classe pour les colonnes des matrices en tant que vecteur.
       Les méthodes modifiant le vecteur vont modifier la matrice.
       ous ne devriez pas avoir besoin de créer d'objet de cette classe,
       cette classe est utilisée en interne par `Mat`."""
    def __init__(self,mat,j):
        self.j   = j
        self.mat = mat
        self.dim = self.mat.dim

    def __getitem__(self,i):
        return self.mat.__coefs__[i][self.j]

    def __setitem__(self,i,x):
        self.mat.__coefs__[i][self.j] = x
       
class Mat:
    """Classe pour les matrices. Les matrices sont représentées en interne
       par des tableaux python, et peuvent être utilisée avec la même syntaxe.
       Un petit changement pour les boucles:
         `for (i, l) in mat:`
       permet d'itérer sur une matrice avec l'index et la chacune des lignes.
       
       Rien dans l'implémentation ne suppose que les coefficients sont
       des nombres. On peut faire des matrices de matrices ou de polynômes.
       Biensur la division par un scalaire ne sera supportée (et les autres
       opérations) que si elle existe.
       
       Les matrices supportent les opérateurs suivants:
       - `M + N` addition de matrice
       - `M += N` addition en place
       - `M - N` soustraction
       - `M -= N` soustraction en place
       - `-M` opposée d'un vecteur
       - `M * N` le produit
       - `M * v` le produit matrice vecteur
       - `v * M` le produit ̀M.transpose() * v`
       - `x * M` et `M * x` le produit par un scalaire.
       - `M *= x` multiplication en place par un scalaire
       - `M / x` division par un scalaire
       - `M /= x` divition en place par un scalaire
       - `M ** n` puissance entière d'une matrice carrée
       - `v == w` et `v != w` test d'égalité/inégalité
       
       Les fonctions habituelles suivantes sont supportées
       - abs(v): carré de la norme euclidienne (somme des carrées des
         coefficients
       - len(v), str(v) et repr(v)
       """
    def __init__(self,*args,coefs=None, dim=None, dim2=None,
                 init=lambda i, j: 0):
        """Constructeur pour les matrices. On peut utiliser:
           `Mat(coefs=tableaux2d)` 
           `Mat(ligne1, ligne2, ...)`
           `Mat(dim=3, init)` pour initialiser une matrice carrée.
           `Mat(dim=3, dim2=4, init` pour initialiser une matrice
            qui peut être rectangulaire.
           Dans les deux derniers cas, `init(i,j)` sera le
           coefficient, ligne i, colonne j. Si `init` n'est pas fournie,
           la matrice nulle est construite.
        """
        if coefs == None and len(args) > 0: coefs = args
        if coefs != None and dim == None: dim = len(coefs)
        if coefs != None and len(coefs) > 0 and dim2 == None: dim2 = len(coefs[0])
        if dim2 == None and dim != None: dim2 = dim
        if coefs == None:
            coefs = [[init(i,j) for j in range(dim2)] for i in range(dim)]
        if len(coefs) != dim or dim == 0: raise ValueError("Mat: Bad dim")
        if len(coefs[0]) != dim2 or dim2 == 0: raise ValueError("Mat: Bad dim2")
        self.__coefs__ = coefs
        self.dim  = dim
        self.dim2  = dim2
        self.lines = [Line(self,i) for i in range(dim)]
        self.cols  = [Col(self,j) for j in range(dim2)]

    def __len__(self):
        return self.dim
    
    def __getitem__(self,i):
        return self.lines[i]

    def __setitem__(self,i,l):
        self.lines[i] = l
        
    def __str__(self):
        return str([[self[i][j] for j in range(self.dim2)]
                    for i in range(self.dim)])

    def __repr__(self):
        return str(self)

    def __iter__(self):
        class Iter:
            def __init__(self, mat):
                self.mat = mat
                self.i = 0
            def __next__(self):
                i = self.i
                if i >= self.mat.dim: raise StopIteration
                self.i += 1
                return (i, self.mat[i])
        return Iter(self)
      
    def transpose(self):
        """`M.transpose()` transpose la matrice M. La matrice n'est pas
        copiée. Toute modification en place de la transposée modifie la
        matrice d'origine."""
        from copy import copy
        self = copy(self)
        self.lines, self.cols = self.cols, self.lines
        self.dim, self.dim2 = self.dim2, self.dim
        return self

    def swap(self,i1,i2):
        """̀M.swap(i1,i2)` ou `M.swap_lines(i1,i2)` échange en place
        les lignes i1 et i2 de M."""
        for j in range(self.dim2):
            self.cols[j].swap(i1,i2)

    swap_lines = swap
    
    def swap_cols(self,j1,j2):
        """̀M.swap_cols(i1,i2)` échange en place les coplonnes i1 et
        i2 de M."""
        for i in range(self.dim):
            self[i].swap(j1,j2)

    def __add__(self, m):
        if not isinstance(m,Mat): m = m * Idt(self.dim)
        if self.dim != m.dim or self.dim2 != m.dim2: raise ValueError("bad dim in Mat.__add__")
        def init(i,j): return self[i][j] + m[i][j]      
        return Mat(dim = self.dim, dim2 = self.dim2, init = init)

    def __radd__(self, m): return self + m
    
    def __iadd__(self, m):
        if not isinstance(m,Mat): m = m * Idt(self.dim)
        if self.dim != m.dim or self.dim2 != m.dim2: raise ValueError("bad dim in Mat.__additem__")
        for i in range(self.dim):
            self[i] += m[i]
        return self

    def __sub__(self, m):
        if not isinstance(m,Mat): m = m * Idt(self.dim)
        if self.dim != m.dim or self.dim2 != m.dim2: raise ValueError("bad dim in Mat.__sub__")
        def init(i,j): return self[i][j] - m[i][j]      
        return Mat(dim = self.dim, dim2 = self.dim2, init = init)

    def __rsub__(self, m): return - self + m
    
    def __isub__(self, m):
        if not isinstance(m,Mat): m = m * Idt(self.dim)
        if self.dim != m.dim or self.dim2 != m.dim2: raise ValueError("bad dim in Mat.__subitem__")
        for i in range(self.dim):
            self[i] -= m[i]
        return self

    def __neg__(self):
        return -1 * self
        
    def __mul__(self, x):
        if isinstance(x, Mat):
            if self.dim2 != x.dim: raise ValueError("bad dim in Mat.__mul__")
            def init(i,j): return self[i] * x.cols[j]
            return Mat(dim=self.dim, dim2 = x.dim2, init = init)
        elif isinstance(x, Vec) and not isinstance(x, Poly):
            if self.dim2 != x.dim: raise ValueError("bad dim in Mat.__mul__")
            def init(i): return self[i] * x
            return Vec(dim=self.dim, init = init)
        else:
            def init(i,j):
                return self[i][j] *x
            return Mat(dim=self.dim, dim2 = self.dim2, init = init)

    def __rmul__(self,x):
        if isinstance(x, Mat):
            if x.dim2 != self.dim: raise ValueError("bad dim in Mat.__mul__")
            def init(i,j): return x[i] * self.cols[j]
            return Mat(dim=x.dim, dim2 = self.dim2, init = init)
        elif isinstance(x, Vec):
            if x.dim != self.dim: raise ValueError("bad dim in Mat.__mul__")
            def init(j): return x * self.cols[j]
            return Vec(dim=self.dim2, init = init)
        else:
            def init(i,j):
                return x * self[i][j]
            return Mat(dim=self.dim, dim2 = self.dim2, init = init)
        
    def __truediv__(self, x):
        def init(i,j):
            return self[i][j] / x
        return Mat(dim=self.dim, dim2 = self.dim2, init = init)

    def __eq__(self, m):
        if not isinstance(m, Mat):
            if m == 0: m = Zero(dim=self.dim, dim2 = self.dim2)
            else: m = m * Idt(self.dim)
        if self.dim != m.dim or self.dim2 != m.dim2: return False
        for i in range(self.dim):
            for j in range(self.dim2):
                if self[i][j] != m[i][j]: return False
        return True
    
    def copy(self):
        """M.copy(): copie la matrice"""
        return Mat(dim=self.dim,dim2=self.dim2,
                   init = lambda i,j: self[i][j])

    def __abs__(self):
        s = 0
        for (_,l) in self:
            s += abs(l)
        return s

    def gauss(self,b,abs=abs):
        """`M.gauss(B)̀: résoud l'équation M X = B. B peut être aussi
           bien une matrice qu'un vecteur. L'algorithme utilisé est
           le pivot de gauss avec un pivot partiel. Le paramètre optionnel
           ̀abs` permet de modifier le choix du pivot (on choisit le pivot
           qui maximise `abs` sur la colonne."""
        b = b.copy()
        a = self.copy()
        if a.dim != b.dim: raise ValueError("Mat.gausss: bad dim")
        for j in range(a.dim2):
            p = j
            for i in range(j+1, a.dim):
                if abs(a[i][j]) > abs(a[p][j]): p = i
            if p != j:
                a.swap(p,j)
                b.swap(p,j)
            b[j] /= a[j][j]
            a[j] /= a[j][j]
            for i in range(j+1, a.dim):
                b[i] -= b[j] * a[i][j]
                a[i] -= a[j] * a[i][j]
        for j in range(a.dim2-1,0,-1):
            for i in range(0,j):
                b[i] -= b[j] * a[i][j]
                a[i] -= a[j] * a[i][j]                
        return b
    
    def det(self,abs=abs,div=None):
        """`M.det()̀: calcule le détermiinant de M. L'algorithme utilisé est
           le pivot de gauss avec un pivot partiel. Le paramètre optionnel
           ̀abs` permet de modifier le choix du pivot (on choisit le pivot
           qui maximise `abs` sur la colonne). Le paramètre optionnel
           `div` est l'algorithme de division utilisée. Dans le cas de
           matrice à coefficient entier ou dans une anneau Euclidien, toutes
           les divisions tomberons juste."""
        a = self.copy()
        if a.dim != a.dim2: raise ValueError("Mat.gausss: bad dim")
        n = 1
        d = 1
        for j in range(a.dim2):
            p = j
            for i in range(j+1, a.dim):
                if abs(a[i][j]) > abs(a[p][j]): p = i
            if p != j:
                a[p].swap(a[j]); n = -n
            x = a[j][j]
            did_mul = False
            if div == None: a[j] /= x
            for i in range(j+1, a.dim):
                q = a[i][j]
                if q != 0:
                    if div != None:
                        a[i] *= x
                        if did_mul: d *= x
                        else: did_mul = True
                    a[i] -= q * a[j]
            if not did_mul: n *= x
        return (n if div == None else div(n,d))

    def carPoly(self):
        """`M.carPoly()` calcule le polynôme caractéristique de M"""
        X = Monome(1, 1)
        M = X * Idt(self.dim) - self
        return M.det(abs = lambda p: p.deg(),div = lambda p, q: p // q)
        
    def inv(self):
        """`M.inv()` calcule l'inverse de M"""
        if self.dim != self.dim2: raise ValueError("Mat.inv: bad dim")
        return self.gauss(Idt(self.dim))

    def __pow__(self,n):
        if n <  0: return self.inv() ** (-n)
        if n == 0: return Idt(self.dim)
        if n == 1: return self
        M = self ** (n // 2)
        M = M * M
        if n %  2 == 0:
            return M
        else: return self * M
        
    
class Zero(Mat):
    """Classe avec un constructeur de la classe Mat pour la matrice nulle"""
    def __init__(self,dim,dim2=None):
        """`Zero(dim)` ou `Zero(dim,dim2)` construit une matrice nulle
        avec les dimensions spécifiés.""" 
        if dim2 == None: dim2 = dim
        return Mat.__init__(self,dim=dim,dim2=dim2)

class Idt(Mat):
    """Classe avec un constructeur de la classe Mat pour la matrice idéntité"""
    def __init__(self,dim):
        """`Idt(dim)` construit la matrice identité de taille dim."""
        return Mat.__init__(self,dim=dim,dim2=dim,init=lambda i,j: 1 if i==j else 0)
    
    
