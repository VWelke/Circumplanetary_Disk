import matplotlib.pyplot as plt
import numpy as np
from mockmodel import r_i, theta_i, phi_i, r_c, theta_c, phi_c, rr, tt, zr

# Plot radial coordinates
plt.figure()
plt.plot(r_i, np.zeros_like(r_i), 'o-')
plt.xlabel('Radial Distance (cm)')
plt.title('Radial Coordinates')
plt.grid(True)
plt.show()

# Plot theta coordinates
plt.figure()
plt.plot(theta_i, np.zeros_like(theta_i), 'o-')
plt.xlabel('Theta (radians)')
plt.title('Theta Coordinates')
plt.grid(True)
plt.show()

# Plot zr coordinates
plt.figure()
plt.imshow(zr[:, :, 0], origin='lower', aspect='auto', cmap='viridis')
plt.colorbar(label='Vertical Coordinate (zr)')
plt.xlabel('Theta Index')
plt.ylabel('Radial Index')
plt.title('Vertical Coordinate (zr) Slice')
plt.grid(True)
plt.show()

