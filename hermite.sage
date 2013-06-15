class Hermite(object):
    def __init__(self,n,k,field,alpha):
        """
        EXAMPLES::

            sage: l.<a> = GF(16) # not tested
            sage: RS = ReedSolomon(5, 3, l,[a^(3*i) for i in range(5)]) # not tested
            sage: RS.G # not tested
            sage: RS.phi # not tested
            sage: RS.V([a^(3*i) for i in range(5)]) in RS.C # not tested
            True

        """
        assert type(n) is sage.rings.integer.Integer
        assert type(k) is sage.rings.integer.Integer
        assert field.is_field()
        self.n=n # test bla
        self.k=k
        self.field=field
        self.V=field^n
        self.W=field^k
        self.d=n-k+1
        self.t=floor((self.d-1)/2)
        self.alpha=alpha
        assert len(alpha)==n
        assert [el.parent()==field for el in alpha] == [True]*n
        basis = [[el^i for el in alpha] for i in range(k)]
        self.G = matrix(basis)
        self.phi=self.W.hom([x*self.G for x in self.W.basis()])
        self.C=self.V.subspace_with_basis(basis)


    def __repr__(self):
        """
        EXAMPLES::

            sage: l.<a> = GF(16) # not tested
            sage: RS = ReedSolomon(5, 3, l,[a^(3*i) for i in range(5)]) # not tested
            sage: RS # not tested
            (5,3)-RS-Code
        """
        return "("+str(self.n)+","+str(self.k)+")-RS-Code"


    def w2pk(self,w):
        """
        EXAMPLES::

            sage: l.<a> = GF(16) # not tested
            sage: RS = ReedSolomon(5, 3, l,[a^(3*i) for i in range(5)]) # not tested
            sage: RS.w2pk(RS.W([0,0,1])) # not tested
            x^2
        """
        R.<x>=PolynomialRing(self.field)
        return R(w.list())

    def encode(self,w):
        """
        EXAMPLES::

            sage: l.<a> = GF(16) # not tested
            sage: RS = ReedSolomon(5, 3, l,[a^(3*i) for i in range(5)]) # not tested
            sage: RS.phi(RS.W([0,0,1])) # not tested
            (1, a^3 + a^2, a^3 + a^2 + a + 1, a^3, a^3 + a)

        """

        return self.phi(w)

    def find_codeword(self,r):
        """
        EXAMPLES::

            sage: l.<a> = GF(16) # not tested
            sage: RS = ReedSolomon(5, 3, l,[a^(3*i) for i in range(5)]) # not tested
            sage: RS.find_codeword(RS.V([0,1,0,a^3,1])) # not tested
            (0, 1, a^3 + a^2 + a + 1, a^3, 1)
            sage: RS.find_codeword(RS.V([a,a^13,a^11,a^14,a^7])) # not tested
            (0, a^3 + a^2 + 1, a^3 + a^2 + a, a^3 + 1, a^3 + a + 1)
        """

        for c in self.C:

            if (c-r).hamming_weight() < self.t+1:
                return c


    def decode(self,r):
        """
        EXAMPLES::

            sage: l.<a> = GF(16) # not tested
            sage: RS = ReedSolomon(5, 3, l,[a^(3*i) for i in range(5)]) # not tested
            sage: RS.decode(RS.V([a,a^13,a^11,a^14,a^7])) # not tested
            (1, 0, 1)
        """
        return column_matrix(self.G)\self.find_codeword(RS.find_codeword(r))
