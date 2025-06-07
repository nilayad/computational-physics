import numpy as np
import matplotlib.pyplot as plt
import time

try:
    from numba import jit
except ImportError:
    def jit(func, **kwargs):
        import warnings
        warnings.warn("The Python `numba` package is not installed.\n"
                      "It's recommended to install it and restart python/jupyter to enable "
                      "just-in-time compilation with `jit`, "
                      "which makes the code run significantly faster.")
        return func


@jit(nopython=True)
def energy(system, i, j, L):
    """Energy function of spins connected to site (i, j)."""
    return -1. * system[i, j] * (system[np.mod(i - 1, L), j] + system[np.mod(i + 1, L), j] +
                                 system[i, np.mod(j - 1, L)] + system[i, np.mod(j + 1, L)])


@jit
def prepare_system(L, rng):
    """Initialize the system."""
    system = 2 * (0.5 - rng.integers(0, 2, size=(L, L)))
    return system


@jit(nopython=True)
def measure_energy(system):
    L = system.shape[0]
    E = 0
    for i in range(L):
        for j in range(L):
            E += energy(system, i, j, L) / 2.
    return E


@jit(nopython=True)
def metropolis_loop(system, T, N_sweeps, N_eq, N_flips, rng):
    """ Main loop doing the Metropolis algorithm."""
    E = measure_energy(system)
    L = system.shape[0]
    E_list = []
    for step in range(N_sweeps + N_eq):
        i = rng.integers(0, L)
        j = rng.integers(0, L)

        dE = -2. * energy(system, i, j, L)
        if dE <= 0.:
            system[i, j] *= -1
            E += dE
        elif np.exp(-1. / T * dE) > rng.random():
            system[i, j] *= -1
            E += dE

        if step >= N_eq and np.mod(step, N_flips) == 0:
            # measurement
            E_list.append(E)
    return np.array(E_list)


def run(T_range, L, N_sweeps, N_eq=None, N_bins=10, N_flips=10, rng=None):
    """run a monte carlo simulation for given system size.

    Paramaters
    ----------
    T_range : temperatures
    L : system size
    N_sweeps : number of steps for the measurements
    N_eq : number of equilibration steps. Defaults to N_sweeps/10.
    N_bins : number of bins used for the error analysis
    N_flips : number of updates between measurements.
    rng : Numpy (pseudo) random number generator.

    Returns
    -------
    E_list : energies (and errors) at the different temperatures
    C_list : specific heat (and errors) at the different temperatures
    """
    if N_eq is None:
        N_eq = N_sweeps//10  # integer division
    if rng is None:
        rng = np.random.default_rng()
    C_list = []
    E_list = []
    system = prepare_system(L, rng)
    for T in T_range:
        C_list_bin = []
        E_list_bin = []
        M_list_bin = []
        M_abs_list_bin = []
        for k in range(N_bins):
            # run N_sweeps metropolis updates and store results the array Es
            Es = metropolis_loop(system, T, N_sweeps, N_eq, N_flips, rng)
            # evaluate observables
            mean_E = np.mean(Es)
            mean_E2 = np.mean(Es**2)
            C_list_bin.append(1. / T**2. / L**2. * (mean_E2 - mean_E**2))
            E_list_bin.append(1. / L**2. * mean_E)
        # estimate error from binning analysis
        C_list.append([np.mean(C_list_bin), np.std(C_list_bin) / np.sqrt(N_bins)])
        E_list.append([np.mean(E_list_bin), np.std(E_list_bin) / np.sqrt(N_bins)])
    return np.array(E_list), np.array(C_list)


#--------------------------------------------------------------------------------------------
@jit(nopython=True)
def metropolis_loop_mod(system, T, N_sweeps, N_eq, N_flips, rng):
    """ Main loop doing the Metropolis algorithm."""
    E = measure_energy(system)
    M =np.sum(system)
    L = system.shape[0]
    E_list = []
    M_list=[]
    for step in range(N_sweeps + N_eq):
        i = rng.integers(0, L)
        j = rng.integers(0, L)

        dE = -2. * energy(system, i, j, L)
        if dE <= 0.:
            system[i, j] *= -1
            E += dE
            M+= 2*system[i,j]
        elif np.exp(-1. / T * dE) > rng.random():
            system[i, j] *= -1
            E += dE
            M+= 2*system[i,j]
            
        if step >= N_eq and np.mod(step, N_flips) == 0:
            # measurement
            E_list.append(E)
            M_list.append(M)
    return np.array(E_list),np.array(M_list)


def run_mod(T_range, L, N_sweeps, N_eq=None, N_bins=10, N_flips=10, rng=None):
    """run a monte carlo simulation for given system size.

    Paramaters
    ----------
    T_range : temperatures
    L : system size
    N_sweeps : number of steps for the measurements
    N_eq : number of equilibration steps. Defaults to N_sweeps/10.
    N_bins : number of bins used for the error analysis
    N_flips : number of updates between measurements.
    rng : Numpy (pseudo) random number generator.

    Returns
    -------
    E_list : energies (and errors) at the different temperatures
    C_list : specific heat (and errors) at the different temperatures
    """
    if N_eq is None:
        N_eq = N_sweeps//10  # integer division
    if rng is None:
        rng = np.random.default_rng()
    C_list = []
    E_list = []
    M_list =[]
    M_abslist =[]
    system = prepare_system(L, rng)
    M= np.sum(system)
    for T in T_range:
        C_list_bin = []
        E_list_bin = []
        M_list_bin = []
        M_abs_list_bin = []
        for k in range(N_bins):
            # run N_sweeps metropolis updates and store results the array Es
            Es,Ms  = metropolis_loop_mod(system, T, N_sweeps, N_eq, N_flips, rng)
            # evaluate observables
            mean_E = np.mean(Es)
            mean_E2 = np.mean(Es**2)
            
            C_list_bin.append(1. / T**2. / L**2. * (mean_E2 - mean_E**2))
            E_list_bin.append(1. / L**2. * mean_E)
            M_list_bin.append((np.mean(Ms))/L**2)
            M_abs_list_bin.append((np.mean(np.abs(Ms)))/L**2)
        # estimate error from binning analysis
        C_list.append([np.mean(C_list_bin), np.std(C_list_bin) / np.sqrt(N_bins)])
        E_list.append([np.mean(E_list_bin), np.std(E_list_bin) / np.sqrt(N_bins)])
        M_list.append([np.mean(M_list_bin), np.std(M_list_bin)/np.sqrt(N_bins)])
        M_abslist.append([np.mean(M_abs_list_bin), np.std(M_abs_list_bin)/np.sqrt(N_bins)])
    return np.array(E_list), np.array(C_list),np.array(M_list), np.array(M_abslist)
#--------------------------------------------------------------------------------------------
#____________________________________________________________________________________________
def metropolis_loop_mod_timed(system, T, N_sweeps, N_eq, N_flips, rng):
    """ Main loop doing the Metropolis algorithm."""
    E = measure_energy(system)
    M =np.sum(system)
    L = system.shape[0]
    E_list = []
    M_list=[]
    time_list=[]
    for step in range(N_sweeps + N_eq):
        i = rng.integers(0, L)
        j = rng.integers(0, L)

        dE = -2. * energy(system, i, j, L)
        if dE <= 0.:
            system[i, j] *= -1
            E += dE
            M+= 2*system[i,j]
        elif np.exp(-1. / T * dE) > rng.random():
            system[i, j] *= -1
            E += dE
            M+= 2*system[i,j]
            
        if step >= N_eq and np.mod(step, N_flips) == 0:
            # measurement
            E_list.append(E)
            M_list.append(M)
            time_list.append(time.time())

    t0 = time_list[0]
    time_list =[t-t0 for t in time_list]
    return np.array(E_list),np.array(M_list),np.array(time_list)
#--------------------------------------------------------------------------------------------
#____________________________________________________________________________________________

#e)

@jit(nopython=True)
def measure_energy_e(system):
    L = system.shape[0]
    E = 0
    for i in range(L):
        for j in range(L):
            E += energy(system, i, j, L) / 2.
    return E
    
@jit(nopython=True)
def metropolis_loop_e(system, T, N_sweeps, N_eq, N_flips, h,rng,):
    """ Main loop doing the Metropolis algorithm."""
    M =np.sum(system)
    E = measure_energy(system) -h*M
    
    L = system.shape[0]
    E_list = []
    M_list=[]
    for step in range(N_sweeps + N_eq):
        i = rng.integers(0, L)
        j = rng.integers(0, L)

        dE = -2. * energy(system, i, j, L) + 2*h*system[i,j]
        if dE <= 0.:
            system[i, j] *= -1
            E += dE
            M+= 2*system[i,j]
        elif np.exp(-1. / T * dE) > rng.random():
            system[i, j] *= -1
            E += dE
            M+= 2*system[i,j]
            
        if step >= N_eq and np.mod(step, N_flips) == 0:
            # measurement
            E_list.append(E)
            M_list.append(M)
    return np.array(E_list),np.array(M_list)


def run_e(T_range, L, N_sweeps, N_eq=None, N_bins=10, N_flips=10, h=None,rng=None):
    """run a monte carlo simulation for given system size.

    Paramaters
    ----------
    T_range : temperatures
    L : system size
    N_sweeps : number of steps for the measurements
    N_eq : number of equilibration steps. Defaults to N_sweeps/10.
    N_bins : number of bins used for the error analysis
    N_flips : number of updates between measurements.
    rng : Numpy (pseudo) random number generator.

    Returns
    -------
    E_list : energies (and errors) at the different temperatures
    C_list : specific heat (and errors) at the different temperatures
    """
    if N_eq is None:
        N_eq = N_sweeps//10  # integer division
    if rng is None:
        rng = np.random.default_rng()

    if h is None:
        h=0
    C_list = []
    E_list = []
    M_list =[]
    M_abslist =[]
    system = prepare_system(L, rng)
    M= np.sum(system)
    for T in T_range:
        C_list_bin = []
        E_list_bin = []
        M_list_bin = []
        M_abs_list_bin = []
        for k in range(N_bins):
            # run N_sweeps metropolis updates and store results the array Es
            Es,Ms  = metropolis_loop_e(system, T, N_sweeps, N_eq, N_flips, h,rng)
            # evaluate observables
            mean_E = np.mean(Es)
            mean_E2 = np.mean(Es**2)
            
            C_list_bin.append(1. / T**2. / L**2. * (mean_E2 - mean_E**2))
            E_list_bin.append(1. / L**2. * mean_E)
            M_list_bin.append((np.mean(Ms))/L**2)
            M_abs_list_bin.append((np.mean(np.abs(Ms)))/L**2)
        # estimate error from binning analysis
        C_list.append([np.mean(C_list_bin), np.std(C_list_bin) / np.sqrt(N_bins)])
        E_list.append([np.mean(E_list_bin), np.std(E_list_bin) / np.sqrt(N_bins)])
        M_list.append([np.mean(M_list_bin), np.std(M_list_bin)/np.sqrt(N_bins)])
        M_abslist.append([np.mean(M_abs_list_bin), np.std(M_abs_list_bin)/np.sqrt(N_bins)])
    return np.array(E_list), np.array(C_list),np.array(M_list), np.array(M_abslist)

#____________________________________________________________________________________________

if __name__ == "__main__":
    T_range = np.arange(1.5, 3.1, 0.1)  # temperatures of interest
    C_list = run(
        T_range=T_range,
        L=10,                      # system size
        N_sweeps=10_000,           # number of sweeps
        N_eq=1_000,                # Number of equilibration steps before the measurements start
        N_flips=10,                # Number of steps between measurements
        N_bins=10,                 # Number of bins used for the error analysis
    )[1]

    # Plot the results
    plt.errorbar(T_range, C_list[:, 0], C_list[:, 1])
    Tc = 2. / np.log(1. + np.sqrt(2))
    print(Tc)
    plt.axvline(Tc, color='r', linestyle='--')
    plt.xlabel('$T$')
  
    plt.ylabel('$c$')
    plt.show()
