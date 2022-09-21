import pandas as pd
from matplotlib import pyplot as plt
plt.rcParams["figure.figsize"] = [7.00, 3.50]
plt.rcParams["figure.autolayout"] = True
columns = ["atoms_number","simulation_time"]
df = pd.read_csv("simulation_time_vs_atoms_number.csv", usecols=columns)
plt.plot(df.atoms_number, df.simulation_time)
plt.xlabel("Atoms number")
plt.ylabel("Simulation time (s)")
plt.title("Simulation time vs Atoms number")
plt.grid(0.1)
plt.savefig("simulation_time_vs_atoms_number.png")
plt.show()

