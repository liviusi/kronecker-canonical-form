from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import os


dirname = "./profiles"
xlines = []

for file in os.listdir(dirname):
    if file.endswith(".prof"):
        with open(f'{dirname}/{file}', "r") as f:
            lines = f.readlines()
            # TODO: parse
    else:
        continue

fig = plt.figure()
ax = plt.axes(projection='3d')
plt.show()