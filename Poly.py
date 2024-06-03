from Vec import Vec

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
        from Mat import Mat
        
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
