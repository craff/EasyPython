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
        import Poly
        import Mat
        if isinstance(x,Vec) and not isinstance(x, Poly.Poly):
            if len(self) != len(x): raise ValueError("Vec.Mul: bad dim")
            s = 0
            for (i, y) in self: s += y * x[i]
            return s
        if isinstance(x,Mat.Mat):
            return x.transpose() * self
        else: return Vec(coefs=[self[j] * x for j in range(len(self))])

    def __neg__(self):
        return -1 * self

    def __rmul__(self,x):
        import Mat
        if isinstance(x,Mat.Mat):
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
