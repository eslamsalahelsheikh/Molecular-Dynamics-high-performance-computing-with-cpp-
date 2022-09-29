import pandas as pd
from matplotlib import pyplot as plt

atoms_number=923

plt.rcParams["figure.figsize"] = [7.00, 3.50]
plt.rcParams["figure.autolayout"] = True
columns = ["cores","real","user","sys"]
df = pd.read_csv(str(atoms_number)+"/equilibrium_on_923_atoms.csv", usecols=columns)
plt.plot(df.cores, df.user)
plt.xlabel("Cores number")
plt.ylabel("Time to demonstrate energy conservation (fs)")
plt.title("Time to demonstrate energy conservation vs Cores number")
plt.grid(0.1)
plt.savefig(str(atoms_number)+"/time_to_demonstrate_energy_conservation_vs_cores_number.png")
plt.show()

columns = ["iteration","total_energy","potential_energy"]
cores8 = pd.read_csv(str(atoms_number)+"/8_energies.csv", usecols=columns)
cores4 = pd.read_csv(str(atoms_number)+"/4_energies.csv", usecols=columns)
cores2 = pd.read_csv(str(atoms_number)+"/2_energies.csv", usecols=columns)
cores1 = pd.read_csv(str(atoms_number)+"/1_energies.csv", usecols=columns)

plt.plot(cores8.iteration, cores8.total_energy)
plt.plot(cores4.iteration, cores4.total_energy)
plt.plot(cores2.iteration, cores2.total_energy)
plt.plot(cores1.iteration, cores1.total_energy)
plt.legend(["8 cores", "4 cores", "2 cores", "1 core"])
plt.xlabel("time steps (fs)")
plt.ylabel("Total energy (eV)")
plt.title("Total energy vs Time steps")
plt.grid(0.1)
plt.savefig(str(atoms_number)+"/total_energy_vs_time_steps.png")
plt.show()