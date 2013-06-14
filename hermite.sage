class ReedSolomon(object):
    def __init__(self,n,k,field,alpha):
        """
        EXAMPLES::

            sage: l.<a> = GF(16)
            sage: RS = ReedSolomon(5, 3, l,[a^(3*i) for i in range(5)])
            sage: RS.G
            [                1                 1
1                 1                 1]
            [                1               a^3         a^3 +
a^2           a^3 + a a^3 + a^2 + a + 1]
            [                1         a^3 + a^2 a^3 + a^2 + a +
1               a^3           a^3 + a]
            sage: RS.phi
            Vector space morphism represented by the matrix:
            [                1                 1
1                 1                 1]
            [                1               a^3         a^3 +
a^2           a^3 + a a^3 + a^2 + a + 1]
            [                1         a^3 + a^2 a^3 + a^2 + a +
1               a^3           a^3 + a]
            Domain: Vector space of dimension 3 over Finite Field in a of
size 2^4
            Codomain: Vector space of dimension 5 over Finite Field in a of
size 2^4
            sage: RS.V([a^(3*i) for i in range(5)]) in RS.C
            True

        """
        assert type(n) is sage.rings.integer.Integer
        assert type(k) is sage.rings.integer.Integer
        assert field.is_field()
        self.n=n
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

            sage: l.<a> = GF(16)
            sage: RS = ReedSolomon(5, 3, l,[a^(3*i) for i in range(5)])
            sage: RS
            (5,3)-RS-Code
        """
        return "("+str(self.n)+","+str(self.k)+")-RS-Code"


    def w2pk(self,w):
        """
        EXAMPLES::

            sage: l.<a> = GF(16)
            sage: RS = ReedSolomon(5, 3, l,[a^(3*i) for i in range(5)])
            sage: RS.w2pk(RS.W([0,0,1]))
            x^2
        """
        R.<x>=PolynomialRing(self.field)
        return R(w.list())

    def encode(self,w):
        """
        EXAMPLES::

            sage: l.<a> = GF(16)
            sage: RS = ReedSolomon(5, 3, l,[a^(3*i) for i in range(5)])
            sage: RS.phi(RS.W([0,0,1]))
            (1, a^3 + a^2, a^3 + a^2 + a + 1, a^3, a^3 + a)

        """

        return self.phi(w)

    def find_codeword(self,r):
        """
        EXAMPLES::

            sage: l.<a> = GF(16)
            sage: RS = ReedSolomon(5, 3, l,[a^(3*i) for i in range(5)])
            sage:
RS.find_codeword(RS.V([0,1,0,a^3,1]))

            (0, 1, a^3 + a^2 + a + 1, a^3, 1)
            sage:
RS.find_codeword(RS.V([a,a^13,a^11,a^14,a^7]))

            (0, a^3 + a^2 + 1, a^3 + a^2 + a, a^3 + 1, a^3 + a + 1)
        """

        for c in self.C:

            if (c-r).hamming_weight() < self.t+1:
                return c


    def decode(self,r):
        """
        EXAMPLES::

            sage: l.<a> = GF(16)
            sage: RS = ReedSolomon(5, 3, l,[a^(3*i) for i in range(5)])
            sage: RS.decode(RS.V([a,a^13,a^11,a^14,a^7]))
            (1, 0, 1)
        """
        return column_matrix(self.G)\self.find_codeword(RS.find_codeword(r))



l.<a> = GF(16)
RS = ReedSolomon(5, 3, l,[a^(3*i) for i in range(5)])

vec =  RS.V([a,a^13,a^11,a^14,a^7])
print 'Codeword c:             ', vec
print 'Encoded from decoded c: ', RS.encode(RS.decode(  vec))

vec =  RS.V([0,1,0,a^3,1])
print 'Codeword c:             ', vec
print 'Encoded from decoded c: ', RS.encode(RS.decode(  vec))

vec =  RS.W([1,0,1])
print 'Plaintext w:            ', vec
print 'Decoded from encoded w: ', RS.decode(RS.encode(  vec))
