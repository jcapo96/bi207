import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
x = sp.symbols('x')
channelNumber = -8

# Define the inclined line as a piecewise function
vertical_line = sp.Piecewise((5, x == 5), (0, True))
# intersection_x = sp.solve(vertical_line - inclined_line, x)
deg = 30*np.pi/180
xs = np.linspace(-48*5, 48*5, 100)
for channelNumber in range(0, 41):
    m = np.tan(deg)
    m2 = -np.tan(deg)
    # x0 = -(channelNumber-7)*(5/np.sin(deg))/np.tan(deg) 
    x0 = (8-channelNumber)*(7.5/np.sin(deg))
    x02 = (8+channelNumber)*(7.5/np.sin(deg))
    inclined_line = m * (x - x0)
    inclined_line2 = m2 * (x - x02)
    ys, ys2 = [], []
    for i in xs:
        result = inclined_line.subs({x: i})
        result2 = inclined_line2.subs({x: i})
        ys.append(result)
        ys2.append(result2)
    plt.plot(xs, ys)
    plt.plot(xs, ys2)
plt.ylim(0, 33*(7.5/np.sin(np.pi/2 - deg)))
plt.xlim(0, 48*5)
plt.savefig("mesh.png")











