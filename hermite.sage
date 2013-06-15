# /home/sage/sage/sage
# %attach "/home/boerner/sage/hermite.sage"
# run_doctests("/home/boerner/sage/hermite.sage")
class Hermite(object):
    def __init__(self,q,m):
        """
        EXAMPLES::

            sage: Her = Hermite(3,4)  
            Finite Field in a of size 3^2
            [[0, 0], [0, a + 1], [0, 2*a + 2], [a, 2*a], [a, a + 2], [a, 1], [a + 1, a], [a + 1, 2*a + 1], [a + 1, 2], [2*a + 1, 2*a], [2*a + 1, a + 2], [2*a + 1, 1], [2, a], [2, 2*a + 1], [2, 2], [2*a, 2*a], [2*a, a + 2], [2*a, 1], [2*a + 2, a], [2*a + 2, 2*a + 1], [2*a + 2, 2], [a + 2, 2*a], [a + 2, a + 2], [a + 2, 1], [1, a], [1, 2*a + 1], [1, 2]]

            sage: Her.G # not tested
            sage: Her.phiHer # not tested
            sage: Her.V([a^(3*i) for i in range(5)]) in Her.C # not tested
            True

        """
        assert type(q) is sage.rings.integer.Integer
        assert type(m) is sage.rings.integer.Integer
        ## Parameters, cf. page 30
        self.n=q^3    # Length of Code
        self.k=m+1-(q^2-1)/2    # Dimension of Code
        self.m_perp=self.n-self.m+q^2+q-2
        self.k_perp=self.n-self.k    # dimension of L(m_perp P)
        self.d_ast=self.n-self.m    # virtual minimal Distance
        self.a_m=q
        self.b_m=-666    # Maximum of b in x^ay^b in L(mP)
        if mod(q,2) == 0:
            self.sigma_q=(q-2)^2/8+1/2    # q = 2^k
        else
            self.sigma_q=(q-1)^2/8+1/2    # q=p^k for p prime, p>2
        self.ell=q*self.a_m+(q+1)*self.b_m-m   
        ##        
        self.q=q
        field.<a>=GF(q^2)
        self.field=field
        print field
        self.V=field^n
        self.W=field^(self.k)
        self.d=n-self.k+1
        self.t=floor((self.d-1)/2)
        R.<x,y>=PolynomialRing(self.field)
        self.hermite_curve=x^(self.q+1)-y^(self.q)-y
        self.points=[[c,d] for c in field for d in field if self.hermite_curve(c,d)==0] 
        print self.points
        #basis = [[el^i for el in alpha] for i in range(self.k)]
        #self.G = matrix(basis)
        #self.phiHer=self.W.hom([x*self.G for x in self.W.basis()])
        #self.C=self.V.subspace_with_basis(basis)


    def __repr__(self):
        """
        EXAMPLES::

            sage: Her = Hermite(5, 3, l,[a^(3*i) for i in range(5)]) # not tested
            sage: Her # not tested
            (5,3)-Her-Code
        """
        return "("+str(self.n)+","+str(self.k)+")-Her-Code"



        

        
        
        

    def w2pk(self,w):
        """
        EXAMPLES::

            sage: l.<a> = GF(16) # not tested
            sage: Her = Hermite(5, 3, l,[a^(3*i) for i in range(5)]) # not tested
            sage: Her.w2pk(Her.W([0,0,1])) # not tested
            x^2
        """
        R.<x,y>=PolynomialRing(self.field)
        return R(w.list())

    def encode(self,w):
        """
        EXAMPLES::

            sage: l.<a> = GF(16) # not tested
            sage: Her = Hermite(5, 3, l,[a^(3*i) for i in range(5)]) # not tested
            sage: Her.phiHer(Her.W([0,0,1])) # not tested
            (1, a^3 + a^2, a^3 + a^2 + a + 1, a^3, a^3 + a)

        """

        return self.phiHer(w)

    def find_codeword(self,r):
        """
        EXAMPLES::

            sage: l.<a> = GF(16) # not tested
            sage: Her = Hermite(5, 3, l,[a^(3*i) for i in range(5)]) # not tested
            sage: Her.find_codeword(Her.V([0,1,0,a^3,1])) # not tested
            (0, 1, a^3 + a^2 + a + 1, a^3, 1)
            sage: Her.find_codeword(Her.V([a,a^13,a^11,a^14,a^7])) # not tested
            (0, a^3 + a^2 + 1, a^3 + a^2 + a, a^3 + 1, a^3 + a + 1)
        """

        for c in self.C:

            if (c-r).hamming_weight() < self.t+1:
                return c


    def decode(self,r):
        """
        EXAMPLES::

            sage: l.<a> = GF(16) # not tested
            sage: Her = Hermite(5, 3, l,[a^(3*i) for i in range(5)]) # not tested
            sage: Her.decode(Her.V([a,a^13,a^11,a^14,a^7])) # not tested
            (1, 0, 1)
        """
        return column_matrix(self.G)\self.find_codeword(Her.find_codeword(r))
