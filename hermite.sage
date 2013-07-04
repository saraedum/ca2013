# /home/sage/sage/sage
# %attach "/home/boerner/sage/hermite.sage"
# run_doctests("/home/boerner/sage/hermite.sage")
class Hermite(object):
    def __init__(self,q,m):
        """
        EXAMPLES::
            sage: her = Hermite(4,51) 
            sage: her.n                                                                                         
            64
            sage: her.q                                                                                         
            4
            sage: her.k                                                                                         
            46
            sage: her.d_ast                                                                                     
            13
            sage: her.ell                                                                            
            12
            sage: her.b_m                                                                            
            4
            sage: her.G 
            51 x 64 dense matrix over Finite Field in a of size 2^4
            sage: her.phi_her # not tested
            sage: her.V([a^(3*i) for i in range(5)]) in her.C # not tested
            True

        """
        assert type(q) is sage.rings.integer.Integer
        assert type(m) is sage.rings.integer.Integer
        assert (q*q-q-1 <= m)    # Def. 3.7
        self.m=m
        self.q=q
        field.<a>=GF(q^2)
        self.field=field
        R.<x,y>=self.field[]
        ## Parameters, cf. page 30
        self.n=q^3    # Length of Code
        self.k=m+1-(q^2-q)/2    # Dimension of Code,    TYPO on PAGE 30
        self.m_perp=self.n-self.m+q^2-q-2 
        self.k_perp=self.n-self.k    # dimension of L(m_perp P)
        self.d_ast=self.n-self.m    # virtual minimal Distance
        self.a_m=q
        self.basis_poly = [x^a*y^b for a in range(q+1) for b in range(m) if ((q*a+(q+1)*b) <= self.m )  ]     # S. 19 
        basis_poly_powers_of_y = [b for a in range(q+1) for b in range(m) if ((q*a+(q+1)*b) <= self.m_perp )]   # S. 23  
        self.b_m=max(basis_poly_powers_of_y)   # Maximum of b in x^ay^b in L(m_perp P)
        if mod(q,2) == 0:
            self.sigma_q=(q-2)^2/8+1/2    # q = 2^k
        else:
            self.sigma_q=(q-1)^2/8+1/2    # q=p^k for p prime, p>2
        self.ell=q*self.a_m+(q+1)*self.b_m-self.m_perp-1   
        ##        

        self.ring=R
        self.V=field^self.n
        self.W=field^(self.k)
        self.d=self.n-self.k+1
        #self.t=floor((self.d-1)/2)
        R.<x,y>=PolynomialRing(self.field)
        self.hermite_curve=x^(self.q+1)-y^(self.q)-y
        self.points=[[c,d] for c in field for d in field if self.hermite_curve(c,d)==0] 
        self.G=matrix([[self.create_element(i)(self.points[j]) for j in range(self.n)] for i in range(self.ord_element(self.m)+1) if self.element_exists(i)])
        self.H=matrix([[self.create_element(i)(self.points[j]) for j in range(self.n)] for i in range(self.ord_element(self.n-self.k)+1) if self.element_exists(i)])

        #self.G = matrix(basis)
        #self.phi_her=self.W.hom([x*self.G for x in self.W.basis()])
        #self.C=self.V.subspace_with_basis(basis)

    def ord(self,f):
        """
        EXAMPLES::
            sage: her=Hermite(3,5)
            sage: R.<x,y>=her.ring
            sage: f=x^2*y
            sage: her.ord(f)
            10
        """
        exps=f.exponents()
        ret=0;
        for exp in exps:
            ret=max(ret,self.q*exp[0]+(self.q+1)*exp[1])
        return ret

    def create_element(self,n):
        """
        EXAMPLES::
            sage: her=Hermite(3,5)
            sage: her.create_element(3)
            x
            sage: her.create_element(17)
            x^3*y^2
        """
        assert type(n) is sage.rings.integer.Integer or n in ZZ
        R.<x,y>=self.ring
        if self.element_exists(n):
            if n==0:
                return R.one()
            else:
                a=-n%(self.q+1)
                b=n%self.q
                if n<(self.q*(self.q+1)):
                    return x^a*y^b
                else:
                    b=(n-self.q*a)/(self.q+1)
                    return x^a*y^b
        else:
            return R.zero()

    def element_exists(self,n):
        """
        EXAMPLES::
            sage: her=Hermite(3,5)
            sage: her.element_exists(2)
            False
            sage: her.element_exists(13)
            True
        """
        assert type(n) is sage.rings.integer.Integer or n in ZZ
        if n>=(self.q*(self.q+1)):
            return True
        else:
            R.<x,y>=self.ring
            a=-n%(self.q+1)
            b=n%self.q
            if self.ord(x^a*y^b)==n:
                return True
            else:
                return False

    def ord_element(self,n):
        """
        EXAMPLES::
            sage: her=Hermite(3,5)
            sage: her.ord_element(4)
            6
        """
        assert type(n) is sage.rings.integer.Integer or n in ZZ
        ret=0
        i=0
        while i<n:
            if self.element_exists(ret):
                i+=1
            ret+=1
        return ret-1

    def L(self,j):
        """
        EXAMPLES::
            sage: her=Hermite(3,5)
            sage: her.L(13)
            [1, x, y, x^2, x*y, y^2, x^3, x^2*y, x*y^2, y^3, x^3*y]
        """
        assert type(j) is sage.rings.integer.Integer or j in ZZ
        return [self.create_element(i) for i in range(j+1) if self.element_exists(i)]

    def __repr__(self):
        """
        EXAMPLES::

            sage: her = Hermite(3,5)
            sage: her 
            (3,5)-Hermite-Code over Finite Field in a of size 3^2
        """
        return "("+str(self.q)+","+str(self.m)+")-Hermite-Code"+" over "+str(self.field)

    # def w2pk(self,w):
        # """
        # EXAMPLES::

            # sage: l.<a> = GF(16) # not tested
            # sage: her = Hermite(5, 3, l,[a^(3*i) for i in range(5)]) # not tested
            # sage: her.w2pk(her.W([0,0,1])) # not tested
            # x^2
        # """
        # R.<x,y>=PolynomialRing(self.field)
        # return R(w.list())

    def encode(self,w):
        """
        EXAMPLES::

            sage: l.<a> = GF(16) # not tested
            sage: her = Hermite(5, 3, l,[a^(3*i) for i in range(5)]) # not tested
            sage: her.phi_her(her.W([0,0,1])) # not tested
            (1, a^3 + a^2, a^3 + a^2 + a + 1, a^3, a^3 + a)

        """

        return self.phi_her(w)
        
        
        
        
    def error_values(self,error_loc):
        """
        EXAMPLES::
            sage: her = Hermite(2,3) 
            sage: her.error_values(1234) # not tested
            got:  [[5, a, a^2], [6, a^2, a]]



        """
        a=self.field.gen()
        error_loc = [[5,a,a^2],[6,a^2,a]]
        print 'got: ', error_loc
        
        Ring.<x,y>=self.ring
        S=x*y+x*x+y
        LL=y+x+1
        R=x^2+y+1
        
        u=[[0,0]]*len(self.points)
        e=[0]*len(self.points)

        print e
        print error_loc[0][0]
        list= [error_loc[j][0] for j in range(len(error_loc))]
        for i in  range(len(error_loc)):
            # #u[error_loc[i][0]][0] = y^(self.q)+y-(error_loc[i][1])^((self.q)+1)
            # #u[error_loc[i][0]][1] = (x-error_loc[i][1])*(y-error_loc[i][2])
            # print 'Position: ', error_loc[i][0]
             f=R/error_loc[i][2]^(self.b_m+1)
             print 'f= ', f
             print 'f P', f(error_loc[i][1],error_loc[i][2])
             print 'f R', R(error_loc[i][1],error_loc[i][2])
             print error_loc[i]
             print e
             g=(x-error_loc[i][1])/LL
             print 'LL(x): ', LL(x,error_loc[i][2])
             print '(x-alph)/LL: ', (x-error_loc[i][1])/LL
             gg=g(x,error_loc[i][2])
             
             print 'gg: ', gg
             e[list[i]]=f(error_loc[i][1],error_loc[i][2])*gg
             print e
             print R(error_loc[i][1],error_loc[i][2])
        #print u
        print e
        

        
        
        
        
        

    def find_codeword(self,r):
        """
        EXAMPLES::

        """

        for c in self.C:

            if (c-r).hamming_weight() < self.t+1:
                return c


    def decode(self,r):
        """
        EXAMPLES::

            sage: l.<a> = GF(16) # not tested
            sage: her = Hermite(5, 3, l,[a^(3*i) for i in range(5)]) # not tested
            sage: her.decode(her.V([a,a^13,a^11,a^14,a^7])) # not tested
            (1, 0, 1)
        """
        return column_matrix(self.G)\self.find_codeword(her.find_codeword(r))
