import matplotlib.pyplot as plt
import numpy as np
func_names = ["xy_euler", "xy_rk4", "sol_analythique"]
for f in func_names:
    lecture = np.loadtxt(open(f, 'rt').readlines())
    plt.plot(lecture[:,0], lecture[:,1], marker = "+", alpha = 0.5, label = f)
plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
plt.show()
plt.close()


names = ["err_euler", "err_rk2heun", "err_rk2euler", "err_rk4"]

#names = ["err_euler","err_rk2heun"]

for n in names:
    file = np.loadtxt(open(n, 'rt').readlines())
    plt.scatter(file[:,0], file[:,1], label = n)

plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
plt.xscale("log")
plt.yscale("log")
plt.xlabel("h")
plt.ylabel("Erreur")
plt.show()
plt.close()