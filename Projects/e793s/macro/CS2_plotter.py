
import matplotlib.pyplot as plt
import numpy as np



xExp, yExp = np.loadtxt('./CS2_TextFiles/18May23_Exp_0000.txt', unpack=True)





x = np.linspace(0, 2 * np.pi, 400)
y = np.sin(x ** 2)

fig = plt.figure()
gs = fig.add_gridspec(2, 2, hspace=0, wspace=0)
(ax1, ax2), (ax3, ax4) = gs.subplots(sharex='col', sharey='row')







fig.suptitle('Sharing x per column, y per row')
ax1.plot(xExp, yExp)
ax2.plot(x, y**2, 'tab:orange')
ax3.plot(x + 1, -y, 'tab:green')
ax4.plot(x + 2, -y**2, 'tab:red')

for ax in fig.get_axes():
    ax.label_outer()

plt.show()
