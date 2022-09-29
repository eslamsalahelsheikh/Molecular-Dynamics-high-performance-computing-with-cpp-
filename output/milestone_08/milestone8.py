import pandas as pd
from matplotlib import pyplot as plt
plt.rcParams["figure.figsize"] = [7.00, 3.50]
plt.rcParams["figure.autolayout"] = True
columns = ["cores","real","user","sys"]
df = pd.read_csv("equilibrium_on_923_atoms.csv", usecols=columns)
plt.plot(df.cores, df.user)
plt.xlabel("Cores number")
plt.ylabel("Time to demonstrate energy conservation")
plt.title("Time to demonstrate energy conservation vs Cores number")
plt.grid(0.1)
plt.savefig("time_to_demonstrate_energy_conservation_vs_cores_number.png")
plt.show()
