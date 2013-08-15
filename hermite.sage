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
            46 x 64 dense matrix over Finite Field in a of size 2^4
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
        self.t=floor((self.d_ast-1)/2)
        self.a_m=q
        self.basis_poly = [x^a*y^b for a in range(q+1) for b in range(m) if ((q*a+(q+1)*b) <= self.m )  ]     # S. 19 
        basis_poly_powers_of_y = [b for a in range(q+1) for b in range(m) if ((q*a+(q+1)*b) <= self.m_perp )]   # S. 23  
        self.b_m=max(basis_poly_powers_of_y)   # Maximum of b in x^ay^b in L(m_perp P)
        if mod(q,2) == 0:
            self.sigma_q=(q-2)^2/8+1/2    # q = 2^k
        else:
            self.sigma_q=(q-1)^2/8+1/2    # q=p^k for p prime, p>2
        self.ell=q*self.a_m+(q+1)*self.b_m-self.m_perp-1   
        self.ring=R
        self.V=field^self.n
        self.W=field^(self.k)
        self.d=self.n-self.k+1
        R.<x,y>=PolynomialRing(self.field)
        self.hermite_curve=x^(self.q+1)-y^(self.q)-y
        self.points=[[c,d] for c in field for d in field if self.hermite_curve(c,d)==0] 
        assert self.n==len(self.points)
        self.G=matrix([[self.basis_polynomial(i)(self.points[j]) for j in range(self.n)] for i in range(self.k)])
        self.H=matrix([[self.basis_polynomial(i)(self.points[j]) for j in range(self.n)] for i in range(self.n-self.k)])
        self.phi=self.W.hom(self.G)
        assert self.G*self.H.transpose()==0
        self.C=self.V.subspace_with_basis(self.G)

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

    def basis_polynomial(self,i):
        """
        EXAMPLES::

            sage: her=Hermite(2,3)
            sage: her.basis_polynomial(0)
            1
            sage: her.basis_polynomial(1)
            x
            sage: her.basis_polynomial(2)
            y

        """
        assert type(i) is sage.rings.integer.Integer or i in ZZ
        return self.create_element(self.ord_element(i))

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
            7

        """
        assert type(n) is sage.rings.integer.Integer or n in ZZ
        ret=0
        i=0
        if n==0:
            return 0
        else:
            while i<n:
                if self.element_exists(ret+1):
                    i+=1
                ret+=1
            return ret

    def L(self,j):
        """
        EXAMPLES::

            sage: her=Hermite(3,5)
            sage: her.L(13)
            [1, x, y, x^2, x*y, y^2, x^3, x^2*y, x*y^2, y^3, x^3*y]

        """
        assert type(j) is sage.rings.integer.Integer or j in ZZ
        return [self.create_element(i) for i in range(j+1) if self.element_exists(i)]

    def syndrome(self,r):
        """
        EXAMPLES::

            sage: her=Hermite(2,3)
            sage: r=her.V([0,0,0,1,1,0,0,0])
            sage: her.H*r
            (0, 1, 1, 1, 0)
            sage: TODO

        """
        assert r in self.V
        return self.H*r

    def syndrome_polynomial(self,r):
        """
        EXAMPLES::

            sage: her=Hermite(2,3)
            sage: r=her.V([0,0,0,1,1,0,0,0])
            sage: her.syndrome_polynomial(r)
            x^2 + x*y + y
            
        """
        s=self.syndrome(r)
        R.<x,y>=self.ring
        return sum([s[i]*x**self.a_m*y**self.b_m//self.basis_polynomial(i) for i in range(self.k_perp)])

    def division_algorithm(self,r,mu):
        """
        EXAMPLES::

            sage: her=Hermite(4,51)
            sage: TODO

        """
        #Sieh Skript Seite 27 Paragraf 3.5
        Ring.<x,y>=self.ring
        S=Ring(self.syndrome_polynomial(r))
        #Schritt 1
        f=y^(self.b_m+1)
        R_0=[S%f]
        gamma=[[self.ring.zero() for j in range(mu+1)] for k in range(mu+1)]
        Delta=[self.ring.zero() for j in range(mu+1)]
        Delta[0]=Ring.one()
        #Schritt 2
        for i in range(1,mu+1):
            phi_i=self.basis_polynomial(i)
            if phi_i.exponents()[0][1]>0:
                psi=y
            else:
                psi=x
            nu=[self.basis_polynomial(k) for k in range(i+1)].index(phi_i/psi)
            theta=self.reduce(psi*R_0[nu])%f
            R=theta
            for j in range(i):
                gamma[i][i-j-1],R=self.euclidian(R,R_0[i-j-1],i,i-j-1)
                assert self.ord(self.basis_polynomial(j)) + self.ord(gamma[i][j]) < self.ord(self.basis_polynomial(i))
            R_0.append(self.reduce(theta-sum([gamma[i][j]*R_0[j] for j in range(i)])))
            assert sum([gamma[i][j]*R_0[j]]) + R_0[i] == theta
            Delta[i]=psi*Delta[nu]-sum([gamma[i][j]*Delta[i] for j in range(i)])
            assert self.ord(Delta[i] - self.basis_polynomial(i)) < self.ord(Delta[i])
            assert (self.reduce(Delta[i] * S - R_0[i])%f).is_zero()
        #Schritt 3
        Lambda=Delta[mu]+sum([gamma[mu][j]*Delta[j] for j in range(mu+1) if self.ord(R_0[j])-self.ord(Delta[mu])<self.ell])
        R=R_0[mu]+sum([gamma[mu][j]*R_0[j] for j in range(mu+1) if self.ord(R_0[j])-self.ord(Delta[mu])<self.ell])
        if (self.reduce(Lambda * S - R)%f).is_zero():
            return Lambda, R
        return self.division_algorithm(r,mu+1)
    
    def euclidian(self,A,B,i,j):
        """
        
        EXAMPLES::

                sage: TODO

        """
        #In Schritt 2 "theta durch R_i-1,...,R_0 in dieser Ordnung teilen"
        if A==0:
            return self.ring.zero(),self.ring.zero()
        A_monos=A.monomials()
        B_monos=B.monomials()
        #Führende Monome holen
        a=A_monos[[self.ord(mono) for mono in A_monos].index(max([self.ord(mono) for mono in A_monos]))]
        a=A.monomial_coefficient(a)*a
        b=B_monos[[self.ord(mono) for mono in B_monos].index(max([self.ord(mono) for mono in B_monos]))]
        b=B.monomial_coefficient(b)*b
        #abgewandelter euklidischer Algorithmus mit Zusatzbedingung. Siehe Gleichung (16)
        if not b.divides(a):
            gamma,R=self.euclidian(A-a,B,i,j)
            R+=a
            return gamma,R
        else:
            c=a//b
            if self.ord(c)<self.ord(self.basis_polynomial(i))-self.ord(self.basis_polynomial(j)):
                A-=B*c
                gamma, R=self.euclidian(A,B,i,j)
                gamma+=c
                return gamma,R
            else:
                gamma,R=self.euclidian(A-a,B,i,j)
                R+=a
                return gamma,R

    def reduce(self,f):
        r"""
        Reduce ``f`` modulo the equation `x^(q+1) = y^q + y` such that it does
        not contain powers of `x` exceeding `q`.

        EXAMPLES::

            sage: her = Hermite(3,5)
            sage: R.<x,y> = her.ring
            sage: her.reduce(x^10)
            x^2*y^6 - x^2*y^4 + x^2*y^2
            sage: her.reduce(y^10)
            y^10
                
        """
        #Magic made by Julian
        R.<x,y>=self.ring
        g=x^(self.q+1)-y^self.q-y
        S.<y>=R.base_ring()[]
        T.<x>=S[]
        U=T.quo(g)
        return R(U(f).lift()(R.gen(0)))

    def __repr__(self):
        """
        EXAMPLES::

            sage: her = Hermite(3,5)
            sage: her 
            (3,5)-Hermite-Code over Finite Field in a of size 3^2

        """
        return "("+str(self.q)+","+str(self.m)+")-Hermite-Code"+" over "+str(self.field)

    def encode(self,w):
        """
        EXAMPLES::

            sage: her = Hermite(3,5)
            sage: her.phi_her(her.W([0,0,1])) # not tested
            (0, a + 1, 2*a + 2, 2*a, a + 2, 1, a, 2*a + 1, 2, 2*a, a + 2, 1, a, 2*a + 1, 2, 2*a, a + 2, 1, a, 2*a + 1, 2, 2*a, a + 2, 1, a, 2*a + 1, 2)

        """
        return self.phi(w)
        
    def error_locations(self,LL):
        """
        EXAMPLES::

            sage: her = Hermite(2,3) 
            sage: field.<a>=GF(2^2)
            sage: R.<x,y>=field[]
            sage: her.error_locations(y+x+1) 
            [[1, 0, 1], [3, a, a + 1], [4, a + 1, a]]

        """     
        a=self.field.gen()
        return[[P[0],P[1][0],P[1][1]] for  P in enumerate(self.points) if LL(P[1][0],P[1][1])==0 ]

    def uniformizer(self,x0,y0):
        """
        EXAMPLES::

            sage: her = Hermite(2,3) 
            sage: a=her.field.gen()
            sage: her.uniformizer(1,1)
            [x - 1, True]
            sage: her.uniformizer(a,a^2)
            [x + a, True]

        """
        a=self.field.gen()
        y=var('y')

        if ((self.q+1)*x0^(self.q)  != 0):
            return [x-x0,True]
        if (self.q*y0^(2*self.q-2)+1 != 0):
            return [y-y0,False]    
        
    def error_values(self,LL,R,S):
        """
        EXAMPLES::

            sage: her = Hermite(2,3) 
            sage: field.<a>=GF(2^2)
            sage: R.<x,y>=field[]
            sage: her.error_values(y+x+1,x*x,x*y+x*x+y) 
            [0, 0, 0, 1, 1, 0, 0, 0]

        """
        a=self.field.gen()
        Ring.<x,y>=self.ring

        error_loc=self.error_locations(LL)
        e=[0]*len(self.points)

        lst= [error_loc[j][0] for j in range(len(error_loc))]
        for i in  range(len(error_loc)):
            if R(error_loc[i][1],error_loc[i][2]) == 0:
                e[lst[i]]=0
                continue
            if (error_loc[i][1]==0) and(error_loc[i][2]) == 0:
                continue
            
            [t,x_is_uniformizer]=self.uniformizer(error_loc[i][1],error_loc[i][2])
            uring.<t>=self.field[[]]
            Hring.<H>=uring[]     

            
            if x_is_uniformizer:
                g=self.hermite_curve(t+error_loc[i][1],H+error_loc[i][2])
            else:
                g=self.hermite_curve(H+error_loc[i][1],t+error_loc[i][2])
            rk=0
            g_prime=diff(g)
            if (g_prime(rk) != 0): 
               for k in range(4):
                     rk=rk-g(rk)/(g_prime(rk))
            if x_is_uniformizer:
                Quo=t/LL(t+error_loc[i][1],rk+error_loc[i][2])
            else:
                Quo=rk/LL(rk+error_loc[i][1],t+error_loc[i][2])
            
            
            e[lst[i]]=-R(error_loc[i][1],error_loc[i][2])/((error_loc[i][2])^(self.b_m+1))*Quo(0)


        e[0]=S.coefficient(x^(self.a_m)*y^(self.b_m))-sum(e)
        return self.V(e)

    def find_codeword(self,r):
        """
        EXAMPLES::

            sage: her = Hermite(3,5)
            sage: a = her.field.gen()
            sage: code_ = her.V((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, a + 1, 0, 0, 0, 0, 0, 0, 0))
            sage: code = her.find_codeword(code_)
            sage: code.is_zero()
            True
            
        """
        Ring.<x,y>=self.ring
        LL,R=self.division_algorithm(r,1)
        S=self.syndrome_polynomial(r)
        e=self.error_values(LL,R,S)
        c=r-e
        return c

    def decode(self,r):
        """
        EXAMPLES::

            sage: TODO

        """
        return self.G.transpose()\self.find_codeword(r)

    def noise(self, max_amount):
        r"""
        EXAMPLES::

            sage: her = Hermite(3,5)
            sage: her.noise(0)
            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
            sage: her.noise(3) # random output
            (0, 2*a + 2, 0, a + 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0)

        """
        pos = range(self.V.dimension())
        from random import shuffle
        shuffle(pos)
        ret = copy(self.V.zero())
        for p in pos[:max_amount]:
            ret[p] = self.V.base_field().random_element()
        return ret

    def test(self, max_noise):
        r"""
        EXAMPLES::

            sage: her = Hermite(3,5)
            sage: her.test(1)

        """
        msg = self.W.random_element()
        code = self.encode(msg)
        code_ = code+self.noise(max_noise)
        msg_ = self.decode(code_)
        assert msg == msg_
