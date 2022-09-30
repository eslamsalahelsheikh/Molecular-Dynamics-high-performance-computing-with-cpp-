import pandas as pd
from matplotlib import pyplot as plt
plt.rcParams["figure.figsize"] = [7.00, 3.50]
plt.rcParams["figure.autolayout"] = True
columns = ["time_step", "total_energy", "kinetic_energy", "potential_energy"]
df = pd.read_csv("energies_vs_time_step.csv", usecols=columns)
plt.plot(df.time_step, df.total_energy)
plt.xlabel("Time Step")
plt.ylabel("Total Energy")
plt.grid()
plt.title("Total Energy vs Time Step")
plt.savefig("total_energy_vs_time_step.png")
plt.show()

