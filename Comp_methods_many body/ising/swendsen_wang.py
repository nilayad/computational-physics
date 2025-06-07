import numpy as np
import time
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components


rng = np.random.default_rng()

#____________________________________________________________________________________________

try:
    from numba import jit
    from numba import njit
except ImportError:
    def jit(func, **kwargs):
        import warnings
        warnings.warn("The Python `numba` package is not installed.\n"
                      "It's recommended to install it and restart python/jupyter to enable "
                      "just-in-time compilation with `jit`, "
                      "which makes the code run significantly faster.")
        return func


#____________________________________________________________________________________________

#@jit(nopython=True)
def init_lat(lx,ly):
    rng = np.random.default_rng()
    lattice = 2 * (0.5 - rng.integers(0, 2, size=(lx, ly)))
    return lattice

@jit(nopython=True)
def gen_bonds(lattice):
    bonds=[]
    lx, ly = lattice.shape
    N = lx * ly
    def red_dim(x, y):
        return x*ly+y
    def ret_dim(n):
        return n//ly, np.mod(n,ly)
        
    for i in range(lx):
        for j in range(ly):
            n= red_dim(i,j)
            b1= red_dim((i+1)%lx, j)
            b2 = red_dim(i,(j+1)% ly )
            bonds.append([n,b1])
            bonds.append([n,b2])
    return np.array(bonds)

@jit(nopython=True)
def update_bond_config(lattice,bonds,T,J=None):
    lx,ly = lattice.shape
    if J == None:
        J = 1
    weights = np.zeros(len(bonds))
    p_w = 1-np.exp(-(2*J)/T)
    for b in range(len(bonds)):
        s1 = (bonds[b,0]//ly, np.mod(bonds[b,0],ly))
        s2 = (bonds[b,1]//ly, np.mod(bonds[b,1],ly))
       #s1 = np.unravel_index(bonds[b,0],shape=(lx,ly))
        #s2 = np.unravel_index(bonds[b,1],shape=(lx,ly))
        if (lattice[*s1] == lattice[*s2] 
            and (np.random.rand() < p_w) ):
            weights[b] = 1

    return weights


@njit
def energy(lattice,bonds):
    sys =np.ravel(lattice)
    return -1*np.sum(sys[bonds[:,0]]*sys[bonds[:,1]])
@njit
def magnetization(lattice):
    return np.sum(lattice)
    
@njit
def abs_magnetization(lattice):
    return np.abs(np.mean(lattice))

#@@jit(nopython=True)
def find_cluster (bonds,weights):
    N = bonds.max()+1
    adj = csr_matrix((weights,(bonds[:,0],bonds[:,1])),shape=(N,N))
    adj = adj + adj.T
    connected_clusters , labels = connected_components(csgraph= adj)
    return adj, connected_clusters,labels
@jit(nopython=True)
def flip_cluster(lattice,cluster,label):
    lx,ly = lattice.shape
    flip_prob = np.random.random(cluster) <0.5
    #sys = lattice.copy()
    for i in range(len(label)):
        if flip_prob[label[i]] :
            lattice[i//ly,np.mod(i,ly)] = lattice[i//ly,np.mod(i,ly)]*-1
    
#@jit
def swendsen_wang_update(lattice,bonds,T):
    
    bond_config = update_bond_config(lattice,bonds,T)
    adj_matrix,con_clust,labls = find_cluster(bonds,bond_config)
    flip_cluster(lattice,con_clust,labls)
#@jit
def swendsen_simulation(lattice,bonds,T,N_sweeps, N_eq):
    if N_eq is None:
        N_eq = N_sweeps//10

    for i in range(N_eq):
        swendsen_wang_update(lattice,bonds,T)
    E= []
    Ms = []
    
    for n in range(N_sweeps):
        swendsen_wang_update(lattice,bonds,T)
        E.append(energy(lattice,bonds))
        Ms.append(magnetization(lattice))
        
    return np.array(E),np.array(Ms)



#--------------------------------------------------------------------------------------------
def a_cor(E, delta):
    if delta ==0:
        return 1
    N = len(E)
    mean_E = np.mean(E)
    var_E = np.var(E)
    prod = [E[i] * E[(i + delta) % N] for i in range(N)]
    return (np.mean(prod) - mean_E**2) / var_E
@njit
def auto_corr(E,delta):
    if delta == 0:
        return 1
        
    de = E-np.mean(E)
    corr = np.mean(de[delta:]*de[:-delta])/np.mean(de**2)
    return corr

#--------------------------------------------------------------------------------------------


#____________________________________________________________________________________________



#@jit(nopython=True)
def run(T_range, L, N_sweeps, N_eq=None, N_bins=10):
    if N_eq is None:
        N_eq = N_sweeps//10
    system = init_lat(L,L)
    bonds = gen_bonds(system)

    
    C_list = []
    E_list = []
    M_list =[]
    M_abslist =[]

    for T in T_range:
        C_list_bin = []
        E_list_bin = []
        M_list_bin = []
        M_abs_list_bin = []
        for k in range(N_bins):
            # run N_sweeps swendsen_simulation and store result
            Es,Ms  = swendsen_simulation(system, bonds,T, N_sweeps,N_eq)
            # evaluate observables

            C_list_bin.append(1. / T**2. / L**2. * np.var(Es))
            E_list_bin.append(1. / L**2. * np.mean(Es))
            M_list_bin.append((np.mean(Ms))/L**2)
            M_abs_list_bin.append((np.mean(np.abs(Ms)))/L**2)
    
        C_list.append([np.mean(C_list_bin), np.std(C_list_bin) / np.sqrt(N_bins)])
        E_list.append([np.mean(E_list_bin), np.std(E_list_bin) / np.sqrt(N_bins)])
        M_list.append([np.mean(M_list_bin), np.std(M_list_bin)/np.sqrt(N_bins)])
        M_abslist.append([np.mean(M_abs_list_bin), np.std(M_abs_list_bin)/np.sqrt(N_bins)])



    
    return np.array(E_list), np.array(C_list),np.array(M_list), np.array(M_abslist)


#____________________________________________________________________________________________

if __name__ == "__main__":
   # T_range = np.arange(1.5, 3.1, 0.1)  # temperatures of interest
    #el,cl,ml,mal = run(T_range,16,100,100,10)

    L_list =[4,8,16,32]
    T_range = np.linspace(1.5,3. , 20)
    Tc = 2. / np.log(1. + np.sqrt(2))
    plt.figure(figsize=(12,6))
    for i in range(len(L_list)):
        l =L_list[i]
        t0=time.time()
        #system= metropolis.prepare_system(l,rng)
        E,c,ml,mal = run(T_range,l,100,10,10)
        print(f'time for L={l:} => {(time.time()-t0): .4f} sec')
        plt.subplot(2,2,1)
        plt.errorbar(T_range, c[:, 0], c[:, 1], fmt='-', capsize=2, label=f'L={l}')
        plt.subplot(2,2,2)
        plt.errorbar(T_range, E[:, 0], E[:, 1], fmt='-', capsize=2, label=f'L={l}')
        plt.subplot(2,2,3)
        plt.errorbar(T_range, ml[:, 0], ml[:, 1], fmt='-', capsize=2, label=f'L={l}')
        plt.subplot(2,2,4)
        plt.errorbar(T_range, mal[:, 0], mal[:, 1], fmt='-', capsize=2, label=f'L={l}')

    plt.subplot(2,2,1)
    plt.xlabel('$T$',fontsize = 12)
    plt.ylabel('$c$',fontsize = 12)
    plt.title(r'$C_v$ Vs T')
    plt.axvline(Tc, color='r', linestyle='--', label=r'$T_{c}$ (critical temp)')
    plt.legend()

    plt.subplot(2,2,2)
    plt.xlabel('$T$',fontsize = 12)
    plt.ylabel('$E$',fontsize = 12)
    plt.title(r'E Vs T')
    plt.axvline(Tc, color='r', linestyle='--', label=r'$T_{c}$ (critical temp)')
    plt.legend()
    plt.subplot(2,2,3)
    plt.xlabel('$T$',fontsize = 12)
    plt.ylabel('$M$',fontsize = 12)
    plt.title(r'M Vs T')
    plt.axvline(Tc, color='r', linestyle='--', label=r'$T_{c}$ (critical temp)')
    plt.legend()
    plt.subplot(2,2,4)
    plt.xlabel('$T$',fontsize = 12)
    plt.ylabel('$M_{abs}$',fontsize = 12)
    plt.title(r'$M_{abs}$ Vs T')
    plt.axvline(Tc, color='r', linestyle='--', label=r'$T_{c}$ (critical temp)')
    plt.legend()


    plt.tight_layout()
    plt.show()
    
