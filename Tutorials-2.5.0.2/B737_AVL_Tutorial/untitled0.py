x = [0, 4, 25, 35]
y = [1, 2, 2, 4]

import matplotlib.pyplot as plt

fig, ax = plt.subplots()
ax.plot(x, y)
ax.set_aspect('equal')
plt.show()