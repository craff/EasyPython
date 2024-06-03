
class Vec:
    """Classe pour les vecteurs. Les vecteurs sont représentés en interne
       par des tableaux python, et peuvent être utilisé avec la même syntaxe.
       Un petit changement pour les boucles:
         `for (i, x) in vec:`
       permet d'itérer sur un vecteur avec l'index est la valeur de chaque
       coefficient.

       Rien dans l'implémentation ne suppose que les coefficients sont
       des nombres. On peut faire des vecteurs de matrices ou de polynômes.
       Biensur la division par un scalaire ne sera supportée (et les autres
       opérations) que si elle existe.

       Par contre, les vecteurs de vecteurs ne sont pas supportés, il faut
       utiliser la classe Mat fournie.
       
       Les vecteurs supportent les opérateurs suivants:
       - `v + w` addition de vecteur
       - `v += w` addition en place
       - `v - w` soustraction
       - `v -= w` soustraction en place
       - `-v` opposée d'un vecteur
       - `v * w` le produit scalaire
       - `v * x` et `x * v` le produit par un scalaire
       - `v *= x` multiplication en place par un scalaire
       - `v / x` division par un scalaire
       - `v /= x` divition en place par un scalaire
       - `v == w` et `v != w` test d'égalité/inégalité
       
       Les fonctions habituelles suivantes sont supportées
       - abs(v): carré de la norme euclidienne
       - len(v), str(v) et repr(v)
    """
    def __init__(self,*args,dim=None,coefs=None,init=lambda i: 0):
        """Constructeur des vecteurs.

           paramètres:
           - `dim`: la dimension (taille) du vecteur
           - `coefs` ou l'ensemble des arguments non nommés:
             les coefficients du vecteur.
           - `init`: une fonction pour initialiser le vecteur.

           Ainsi les trois formes suivantes sont équivalents:
           - `Vec(0,1,2,3)`
           - `Vec(coefs=[0,1,2,3])`
           - `Vec(dim=4, init=lambda i: i)

           raise `ValueError` si les paramètres sont incohérents ou
           que la dimension est <= 0"""
        if coefs == None and len(args) > 0: coefs = args
        if coefs != None and dim == None: dim = len(coefs)
        if coefs == None:
            coefs = [init(i) for i in range(dim)] 
        if len(coefs) != dim or dim <= 0: raise ValueError("Vec: Bad dim")
        self.__coefs__ = coefs
        self.dim  = dim

    def __len__(self):
        return self.dim
    
    def __getitem__(self,i):
        return self.__coefs__[i]

    def __setitem__(self,i,x):
        self.__coefs__[i] = x

    def __iter__(self):
        class Iter:
            def __init__(self, vec):
                self.vec = vec
                self.i = 0
            def __next__(self):
                i = self.i
                if i >= self.vec.dim: raise StopIteration
                self.i += 1
                return (i, self.vec[i])
        return Iter(self)
 
    def __str__(self):
        return str([x for (_,x) in self])

    def __repr__(self):
        return str(self)
    
    def __add__(self, vec):
        return [x + vec[i] for (i, x) in self]

    def __iadd__(self, vec):
        for (i, _) in self: self[i] += vec[i]
        return self
            
    def __sub__(self, vec):
        return [x - vec[i] for (i, x) in self]

    def __isub__(self, vec):
        for (i, _) in self: self[i] -= vec[i]
        return self

    def __mul__(self, x):
        if isinstance(x,Vec) and not isinstance(x, Poly):
            if len(self) != len(x): raise ValueError("Vec.Mul: bad dim")
            s = 0
            for (i, y) in self: s += y * x[i]
            return s
        if isinstance(x,Mat):
            return x.transpose() * self
        else: return Vec(coefs=[self[j] * x for j in range(len(self))])

    def __neg__(self):
        return -1 * self

    def __rmul__(self,x):
        if isinstance(x,Mat):
            return x * self
        else: return self.__mul__(x)
    
    def __imul__(self, x):
        for (i, _) in self: self[i] *= x
        return self

    def __truediv__(self, x):
        return Vec(coefs=[y / x for (_,y) in self])
        
    def __itruediv__(self, x):
        for (i, _) in self: self[i] /= x
        return self

    def __abs__(self):
        from math import sqrt
        s = 0
        for (_, x) in self: s += x * x 
        return s

    def swap(self,i,j):
        """`v.swap(i,j)` échange les coefficients i et j en place"""
        self[i], self[j] = self[j], self[i]

    def copy(self):
        """`v.copy()` retourne une copie de v"""
        return Vec(dim=self.dim,init = lambda i: self[i])

    def __eq__(self, vec):
        if len(self) != len(vec): return False
        for (i, x) in self:
            if x != vec[i]: return False
        return True       

class Poly(Vec):
    """La classe Poly, qui hérite de Vec, est une représentation
       des polynômes à une variable.

       En plus des opérateurs sur les vecteurs, on a:
       - P * Q: change de sens, c'est le produit de polynômes,
         par le produit scalaire.
       - P // Q: quotient de la division euclidienne
       - P % Q: reste de la divition euclidienne

       Les polynômes sont "callable". On peut écrire P(x) pour évaluer
       le polynôme en x. On peut évaluer un polynôme sur une matrice carrée.
       """
    def deg(self):
        """`P.deg()` donne le degré de P"""
        i = len(self) - 1
        while i > 0 and self[i] == self.zero: i -= 1
        return i

    def __init__(self,*args,deg=None,coefs=None,init=None, zero = 0):
        """Constructeur des polynômes
           paramètres:
           - `deg`: le degré 
           - `coefs` ou l'ensemble des arguments non nommés:
             les coefficients du polynôme.
           - `init`: une fonction pour initialiser le polynôme.
           - `zero`: la représentation des coefficients nulles.
             Remarque: on peut laisser la valeur par défaut 0
             pour les matrices car 0 + M, 0 * M, etc
             fonctionne.

           Ainsi les trois formes suivantes sont équivalentes:
           - `Poly(0,1,2,3)`
           - `Poly(coefs=[0,1,2,3])`
           - `Poly(deg=3, init=lambda i: i)
        """ 
        if coefs == None and len(args) > 0: coefs = args
        if init == None: init = lambda i: zero
        if deg == None and coefs != None:
            deg = len(coefs) - 1
        Vec.__init__(self,dim=deg+1,coefs=coefs,init=init)
        self.zero = zero

    def __str__(self):       
        for (i,x) in self:
            if i == 0: s = str(x)
            else:
                s += ("+" + str(x) if x > 0 else "-" + str(-x)) + ("X^" + str(i) if i > 1 else "X")
        return(s)
    
    def __getitem__(self,i):
        if i >= 0 and i < len(self): return self.__coefs__[i]
        else: return self.zero
        
    def cl(self, a, p, shift = 0):
        """`P.cl(a,Q)` renvoie P + a Q"""
        return Poly(coefs=[self[i] + a * p[i-shift]
                for i in range(max(len(self), len(p) + shift))])
    
    def __add__(self, p):
        if isinstance(p,Poly): return self.cl(1,p)
        else: return self + Monome(0, p)

    __iadd__ = __add__
        
    def __radd__(self, p):
        if isinstance(p,Poly): return self.cl(1,p)
        else: return self + Monome(0, p)
        
    def __sub__(self, p):
        if isinstance(p,Poly): return self.cl(-1,p)
        else: return self - Monome(0, p)

    __isub__ = __sub__ 

    def __rsub__(self, p):
        if isinstance(p,Poly): return self.cl(-1,p)
        else: return Monome(0, p) - self
    
    def __mul__(self, p):
        if (not isinstance(p,Poly) and
            (isinstance(p, Vec) or isinstance(p,Mat))):
            return p * self
        if isinstance(p, Poly):
            r = Poly(deg=p.deg() + self.deg(), zero = self.zero)
            for i in range(0,len(self)):
                r = r.cl(self[i], p, shift = i)
            return(r)
        M = Monome(0,p,self.zero)
        return self * M

    __imul__ = __mul__
    
    def euclidian_div(self, p):
        """`D.euclidian_div(P)` retourne un couple `(Q,R)` tel que
           `P = Q D + R et R.deg() < D.deg()`. C'est bien sur la
           division euclidienne des polynômes à une variable."""
        d = self.deg()
        b = self[d]
        if d == 0:
            if b == 1: return (p, 0)
            return (Poly(deg = p.deg(), init = lambda i: p[i] / b), 0)
        r = p.copy()
        q = Poly(deg = 0)
        
        while r.deg() >= d:
            dr = r.deg()
            a = r[dr] / b
            q += Monome(dr - d, a, zero = self.zero)
            
            r = r.cl(-a,self,shift=dr - d)
            r[dr] = self.zero
        return (q, r)

    def __floordiv__(self, p):
        if not isinstance(p, Poly): p = Monome(0,p,zero=self.zero)
        return p.euclidian_div(self)[0]

    def __mod__(self,p):
        if not isinstance(p, Poly): p = Monome(0,p,zero=self.zero)
        return p.euclidian_div(self)[1]

    def copy(self):
        """`P.copy()` copie le polynéôme P"""
        return Poly(deg=self.deg(),init = lambda i: self[i])

    def __eq__(self, vec):
        if not isinstance(vec, Poly): vec = Monome(0, vec, zero=self.zero)
        if self.deg() != vec.deg(): return False
        for i in range(self.deg() + 1):
            if self[i] != vec[i]: return False
        return True

    def __call__(self,x):
        r = self[self.deg()]
        for i in range(self.deg()-1, -1, -1):
            r *= x
            r += self[i]
        return(r)
    
class Monome(Poly):
    """Contructeur de la classe `Polynome pour créer un monome:
       `Poly(d,c)` crée le polynôme c X^d. L'argument optionel zero, permet
       de préciser la représentation du coefficient zero."""
    def __init__(self,deg,c, zero = 0):
        Poly.__init__(self,deg=deg, init =  lambda i: c if i == deg else 0)
        
class Line(Vec):
    def __init__(self,mat,i):
        self.i   = i
        self.mat = mat
        self.dim = self.mat.dim2
    
    def __getitem__(self,j):
        return self.mat.__coefs__[self.i][j]

    def __setitem__(self,j,x):
        self.mat.__coefs__[self.i][j] = x

class Col(Vec):
    def __init__(self,mat,j):
        self.j   = j
        self.mat = mat
        self.dim = self.mat.dim

    def __getitem__(self,i):
        return self.mat.__coefs__[i][self.j]

    def __setitem__(self,i,x):
        self.mat.__coefs__[i][self.j] = x
       
class Mat:
    def __init__(self,coefs=None, dim=None, dim2=None, init=lambda i, j: 0):
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

    def __str__(self):
        return str([[self[i][j] for j in range(self.dim2)]
                    for i in range(self.dim)])

    def __repr__(self):
        return str(self)
        
    def get(self,i,j):
        return self[i][j]

    def transpose(self):
        from copy import copy
        self = copy(self)
        self.lines, self.cols = self.cols, self.lines
        self.dim, self.dim2 = self.dim2, self.dim
        return self

    def __len__(self):
        return self.dim
    
    def __getitem__(self,i):
        return self.lines[i]

    def __setitem__(self,i,l):
        self.lines[i] = l
        
    def swap(self,i1,i2):
        for j in range(self.dim2):
            self.cols[j].swap(i1,i2)

    swap_lines = swap
    
    def swap_cols(self,j1,j2):
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
        return Mat(dim=self.dim,dim2=self.dim2,
                   init = lambda i,j: self[i][j])

    def gauss(self,b,abs=abs):
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

    def dot(self,m):
        s = 0
        for i in range(self.dim):
            s += self[i] * m[i]
        return s

    def norm2(self):
        return self.dot(self)
    
    def det(self,abs=abs,div=None):
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
        X = Monome(1, 1)
        M = X * Idt(self.dim) - self
        return M.det(abs = lambda p: p.deg(),div = lambda p, q: p // q)
        
    def inv(self):
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
    def __init__(self,dim,dim2=None):
        if dim2 == None: dim2 = dim
        return Mat.__init__(self,dim=dim,dim2=dim2)

class Idt(Mat):
    def __init__(self,dim):
        return Mat.__init__(self,dim=dim,dim2=dim,init=lambda i,j: 1 if i==j else 0)
    
    
if __name__ == "__main__" or True:
    from unittest import main as startTests, TestCase

    class Test(TestCase):
        def testPoly(self):
            P = Poly(coefs=[1,2,3]) 
            Q = Poly(coefs=[1,-1,2,-2,0])
            self.assertEqual(P.deg(), 2)
            self.assertEqual(Q.deg(), 3)
            self.assertEqual(P + Q, Poly(coefs=[2,1,5,-2]))
            self.assertEqual(P - Q, Poly(coefs=[0,3,1,2]))
            self.assertEqual(P * Q, Poly(coefs=[1,1,3,-1,2,-6]))
            self.assertEqual(Q * P, Poly(coefs=[1,1,3,-1,2,-6]))
            (D,R) = P.euclidian_div(Q)
            self.assertAlmostEqual(abs(P * D + R - Q), 0)
            self.assertEqual(P(1), 6)
            self.assertEqual(P(2), 17)
            self.assertEqual(P(Idt(2)), 6 * Idt(2))
            self.assertEqual(P(2 * Idt(2)), 17 * Idt(2))
        def testMat1(self):
            M = Idt(3)
            for i in range(3):
                for j in range(3):
                    x = 1 if i == j else 0
                    self.assertEqual(M[i][j], x)
                    self.assertEqual(M.lines[i][j], x)
                    self.assertEqual(M.cols[j][i], x)
        def testMat2(self):
            M = Idt(3)
            self.assertTrue(M == (M + M) - M)
            self.assertTrue(M == M * 2 - M)
            self.assertTrue(M == M*M)
        def testMatInv(self):
            M = Mat(coefs=[[2,1],[-1,0]])
            self.assertAlmostEqual((M * M.inv() - Idt(2)).norm2(), 0)
            self.assertAlmostEqual((M.inv() * M - Idt(2)).norm2(), 0)
            self.assertEqual(M.carPoly()(M), Zero(2))
            M = Mat(coefs=[[1,2,3],[4,5,6],[7,8,1]])
            self.assertAlmostEqual((M * M.inv() - Idt(3)).norm2(), 0)
            self.assertAlmostEqual((M.inv() * M - Idt(3)).norm2(), 0)
            self.assertEqual(M.carPoly()(M), Zero(3))
    startTests()
    
