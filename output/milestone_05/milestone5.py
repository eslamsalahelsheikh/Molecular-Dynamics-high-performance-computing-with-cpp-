import pandas as pd
from matplotlib import pyplot as plt
plt.rcParams["figure.figsize"] = [7.00, 3.50]
plt.rcParams["figure.autolayout"] = True
columns = ["cluster_size","simulation_time"]
df = pd.read_csv("simulation_time_vs_cluster_size.csv", usecols=columns)
plt.plot(df.cluster_size, df.simulation_time)
plt.xlabel("Cluster size")
plt.ylabel("Simulation time (s)")
plt.title("Simulation Time vs Cluster Size")
plt.grid(0.1)
plt.savefig("simulation_time_vs_cluster_size.png")
plt.show()

