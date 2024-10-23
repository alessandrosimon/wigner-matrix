import numpy as np
import numba    


'''
define necessary matrices first
'''

@numba.jit
def a(l,m,mp):
    if abs(m) == l:
        return 0.
    elif (abs(mp) == l):
        return 0.
    else:
        return np.sqrt(((l+m)*(l-m))/((l+mp)*(l-mp)))
    
@numba.jit
def b(l,m,mp):
    if (m == -l or m == -l+1) :
        return 0.
    elif (abs(mp) == l):
        return 0.
    else:
        return np.sqrt(((l+m)*(l+m-1))/(2.*(l+mp)*(l-mp)))
    
@numba.jit
def c(l,m,mp):
    if (abs(m) == l or mp == -l or mp-1 == -l):
        return 0.
    else:
        return np.sqrt((2.*(l+m)*(l-m))/((l+mp)*(l+mp-1)))
    
@numba.jit
def d(l,m,mp):
    if (m == -l or m ==(-l+1) or mp == -l or mp-1 == -l):
        return 0.
    else:
        return np.sqrt(((l+m)*(l+m-1))/((l+mp)*(l+mp-1)))


def get_a(l):
    dim = 2*l+1
    al = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            I = i-l
            J = j-l
            al[i,j] = a(l, I, J)
    return al

def get_b(l):
    dim = 2*l+1
    bl = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            I = i-l
            J = j-l
            bl[i,j] = b(l, I, J)
    return bl


def get_c(l):
    dim = 2*l+1
    cl = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            I = i-l
            J = j-l
            cl[i,j] = c(l, I, J)
    return cl


def get_d(l):
    dim = 2*l+1
    dl = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            I = i-l
            J = j-l
            dl[i,j] = d(l, I, J)
    return dl

def checkerboard(shape):
    return np.indices(shape).sum(axis=0) % 2

class cwigner:
    def __init__(self, lmax:int) -> None:
        self.lmax = lmax
        self.cached_al = []
        self.cached_bl = []
        self.cached_cl = []
        self.cached_dl = []

        for l in range(1, lmax+1):
            self.cached_al.append(get_a(l))
            self.cached_bl.append(get_b(l))
            self.cached_cl.append(get_c(l))
            self.cached_dl.append(get_d(l))

 
        
    def get_FG(self, R):
        sqrt2 = np.sqrt(2)
        D = np.zeros((3, 3), dtype=complex)
        F = np.array([
            [(R[1, 1]+R[0,0])/2.0, R[0,2]/sqrt2, (R[1,1]-R[0,0])/2.],
            [R[2,0]/sqrt2, R[2,2], -R[2,0]/sqrt2],
            [(R[1,1]-R[0,0])/2.0, -R[0,2]/sqrt2, (R[1,1]+R[0,0])/2.0]
        ])
        G = np.array([
            [(R[1, 0]-R[0,1])/2.0, R[1,2]/sqrt2, -(R[1,0]+R[0,1])/2.],
            [-R[2,1]/sqrt2, 0.0, -R[2,1]/sqrt2],
            [(R[1,0]+R[0,1])/2.0, R[1,2]/sqrt2, (R[0,1]-R[1,0])/2.0]
        ])
        return F, G
   


    def calc_D(self, R):
        lmax = self.lmax
        F1, G1 = self.get_FG(R)
        
        D1 = F1+1.j*G1
        if lmax == 1:
            return D1
        
        Dtot = [D1]
        l = 2
        F = F1
        G = G1
        while l <= lmax:
            # al = get_a(l)
            # bl = get_b(l)
            # cl = get_c(l)
            # dl = get_d(l)

            cl = self.cached_cl[l-1]
            dl = self.cached_dl[l-1]
                    
            def H(i,j):
                return F1[i,j]*F - G1[i,j]*G

            def K(i,j):
                return F1[i,j]*G + G1[i,j]*F
            
    #         Flmid = al * np.pad(H(1, 1), 1)
    #         Flmid += bl * np.pad(H(2, 1), ((0, 2), (1, 1)))
    #         Flmid += bl[::-1,:] * np.pad(H(0, 1), ((2, 0), (1, 1)))

    #         Glmid = al * np.pad(K(1, 1), 1)
    #         Glmid += bl * np.pad(K(2, 1), ((0, 2), (1, 1)))
    #         Glmid += bl[::-1,:] * np.pad(K(0, 1), ((2, 0), (1, 1)))
            

            Flright = cl * np.pad(H(1, 2), ((1, 1), (2, 0)))
            Flright += dl * np.pad(H(2, 2), ((2, 0), (2, 0)))
            Flright += dl[::-1, :] * np.pad(H(0, 2), ((0, 2), (2, 0)))
            
            Glright = cl * np.pad(K(1, 2), ((1, 1), (2, 0)))
            Glright += dl * np.pad(K(2, 2), ((2, 0), (2, 0)))
            Glright += dl[::-1, :] * np.pad(K(0, 2), ((0, 2), (2, 0)))
            
            
            Flleft = cl[:,::-1] * np.pad(H(1,0), ((1, 1), (0, 2)))
            Flleft += dl[:,::-1] * np.pad(H(2,0), ((2, 0), (0, 2)))
            Flleft += dl[::-1, ::-1] * np.pad(H(0, 0), ((0, 2), (0, 2)))

            Glleft = cl[:,::-1] * np.pad(K(1, 0), ((1, 1), (0, 2)))
            Glleft += dl[:,::-1] * np.pad(K(2, 0), ((2, 0), (0, 2)))
            Glleft += dl[::-1, ::-1] * np.pad(K(0, 0), ((0, 2), (0, 2)))
            
            
            Dlright = Flright + 1.j*Glright
            Dlleft = Flleft + 1.j*Glleft
            
            
            Dl = np.zeros_like(Dlleft)
            Dl[:, :l] = Dlleft[:,:l]
            Dl[:,l:] = Dlright[:,l:]
            

            Dtot.append(Dl)
            F = np.zeros_like(Flleft)
            F[:, :l] = Flleft[:,:l]
            F[:, l:] = Flright[:, l:]
            G = np.zeros_like(Glleft)
            G[:, :l] = Glleft[:,:l]
            G[:, l:] = Glright[:, l:]
            l += 1
        return Dtot

