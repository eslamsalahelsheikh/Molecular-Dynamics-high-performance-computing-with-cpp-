import pandas as pd
from matplotlib import pyplot as plt
import os
import glob
import re

plt.rcParams["figure.figsize"] = [7.00, 3.50]
plt.rcParams["figure.autolayout"] = True
columns = ["iteration","total_energy","average_temp","potential_energy"]

directory = '/home/eslam/Desktop/Molecular-Dynamics/output/milestone_07'
extension = 'csv'

def num_sort(test_string):
    return list(map(int, re.findall(r'\d+', test_string)))[0]

os.chdir(directory)
energy_files = glob.glob('**/*.{}'.format(extension))
energy_files.sort(key=num_sort)
legend = []

heat_capacity = []
start_energies = []
for file in energy_files:
    energy_file = pd.read_csv(file, usecols=columns)
    plt.plot(energy_file.average_temp, energy_file.total_energy)
    legend.append(file.split("/")[0])
    delta_energy = energy_file.total_energy[-1:]-energy_file.total_energy[0]
    delta_temp = energy_file.average_temp[-1:]-energy_file.average_temp[0]
    heat_capacity.append(delta_energy/delta_temp)
    start_energies.append(energy_file.total_energy[0])

plt.xlabel("average_temp")
plt.ylabel("total_energy (s)")
plt.title("total_energy vs average_temp")
plt.grid(0.1)
plt.legend(legend)
plt.savefig("total_energy_vs_average_temp.png")
plt.show()

columns2 = ["cluster_size","melting_point","number_of_iterations_till_melting","added_energy","atoms_number","latent_energy"]
melting_file = pd.read_csv("/home/eslam/Desktop/Molecular-Dynamics/output/milestone_07/melting_points_vs_cluster_size.csv", usecols=columns2)
plt.plot(melting_file.cluster_size, melting_file.melting_point)
plt.xlabel("cluster_size")
plt.ylabel("melting_point (K)")
plt.title("Melting point vs Cluster size")
plt.grid(0.1)
plt.savefig("melting_point_vs_cluster_size.png")
plt.show()

plt.plot(melting_file.cluster_size, heat_capacity)
plt.xlabel("cluster_size")
plt.ylabel("heat_capacity (J/K)")
plt.title("Heat capacity vs Cluster size")
plt.grid(0.1)
plt.savefig("heat_capacity_vs_cluster_size.png")
plt.show()


plt.plot(melting_file.cluster_size, melting_file.latent_energy)
plt.xlabel("cluster_size")
plt.ylabel("latent_heat (J)")
plt.title("Latent heat vs Cluster size")
plt.grid(0.1)
plt.savefig("latent_heat_vs_cluster_size.png")
plt.show()
