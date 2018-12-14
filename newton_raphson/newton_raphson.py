import numpy as np
import pandas as pd
import scipy
import scipy.optimize as opt
import matplotlib.pyplot as plt

def function(x):
    out = [0,0,0]
    out[0] = 4*x[0]**3 + 8
    out[1] = -5*x[1]**1 + 4
    out[2] = 4*x[2]**2 - 10 
    # print(out)   
    return out

def load_flow(M_interco, V_list, theta, G, B):
    """
    function to calculate load flow

    input:
    M_interco : matrix with the following shape
    [
        [0,1,0,0,0],
        [1,0,1,0,0],
        [0,1,0,1,0],
        [0,0,1,0,1],
        [0,0,0,1,0],
    ]
    allow to know which bus is connected to another
    V_list: Tension of each bus. Has the following shape
    [
        [V_phase_1, V_phase_2, V_phase_3],
        [V_phase_1, V_phase_2, V_phase_3],
        [V_phase_1, V_phase_2, V_phase_3],
        [V_phase_1, V_phase_2, V_phase_3],
    ]
    theta: Has the following shape
    [
        [delta_phase_1, delta_phase_2, delta_phase_3],
        [delta_phase_1, delta_phase_2, delta_phase_3],
        [delta_phase_1, delta_phase_2, delta_phase_3],
        [delta_phase_1, delta_phase_2, delta_phase_3],
    ]
    G: real part admitance between two buses
        hypotesis:
        only three phases for now
        the matrix admitance is the same between every buses
    symetric and positive
    [
        [g11, g12, g13],
        [g21, g22, g23],
        [g31, g32, g33],
    ]
    B: imaginary part admitance between two buses
        hypotesis:
        only three phases for now
        the matrix admitance is the same between every buses
    symetric and positive
    [
        [g11, g12, g13],
        [g21, g22, g23],
        [g31, g32, g33],
    ]
    """
    P_load = [[0,0,0] for k in M_interco]
    Q_load = [[0,0,0] for k in M_interco]
    for i, bus in enumerate(M_interco):
        for phase in range(3):
            P_load[i][phase] = sum([ 
                sum([
                    (
                        z
                        * V_list[i][phase]
                        * V_list[voisin][k]
                        * (
                            G[phase][k] * np.cos(theta[i][phase] - theta[voisin][k])
                            + B[phase][k] * np.sin(theta[i][phase] - theta[voisin][k])
                        )
                    ) for k in range(3)
                ]) for voisin, z in enumerate(bus)
            ]) * 0.001
            Q_load[i][phase] = sum([ 
                sum([
                    (
                        z
                        * V_list[i][phase]
                        * V_list[voisin][k]
                        * (
                            G[phase][k] * np.sin(theta[i][phase] - theta[voisin][k])
                            + B[phase][k] * np.cos(theta[i][phase] - theta[voisin][k])
                        )
                    ) for k in range(3)
                ]) for voisin, z in enumerate(bus)
            ]) * 0.001
    return P_load, Q_load

def plotting_load_flow(P_load, Q_load):
    
    class MyFigure(plt.Figure):
        def __init__(self, *args, figtitle=None, **kwargs):
            """
            custom kwarg figtitle is a figure title
            """
            super().__init__(*args, **kwargs)
            if not figtitle:
                pass
            self.text(
                0.5,
                0.95,
                figtitle,
                ha='center',
                fontsize=18
            )
    
    plt.figure(
        figsize=(10,12),
        FigureClass=MyFigure,
        figtitle="LOAD FLOWS REGARDING THE BUSES NUMBERS"
    )
    plt.subplot(211)
    plt.grid(True)
    plt.title("ACTIVE POWER")
    plt.xlabel("bus number")
    plt.ylabel("power (kW)")
    plt.plot(P_load, "-^")
    plt.legend(["phase 1", "phase 2", "phase 3"])
    plt.hold()
    plt.subplot(212)
    plt.title("REACTIVE  POWER")
    plt.grid(True)
    plt.xlabel("bus number")
    plt.ylabel("power (kVA)")
    plt.plot(Q_load, "-^")
    plt.legend(["phase 1", "phase 2", "phase 3"])
    plt.show()

if __name__ == "__main__":

    # Here is an example of newton raphson use
    # on a function described earlier
    solution = opt.fsolve(
        function,
        x0=[-1.3,0.78,1.2],
    )

    # Here is an example of interconnextion matrix
    # 0 -- 1 -- 2 -- 3 -- 4
    M_interco = [
        [0,1,0,0,0],
        [1,0,1,0,0],
        [0,1,0,1,0],
        [0,0,1,0,1],
        [0,0,0,1,0],
    ]
    V_list = [
        [120,120,120],
        [110,120,120],
        [120,110,120],
        [120,120,110],
        [100,100,100],
    ]
    theta = [
        [0,120,240],
        [0,120,240],
        [0,120,240],
        [0,120,240],
        [0,120,240],
    ]
    G = [
        [50,0,0],
        [0,50,0],
        [0,0,50],
    ]
    B = [
        [150,0,0],
        [0,150,0],
        [0,0,150],
    ]
    P, Q = load_flow(M_interco, V_list, theta, G, B)
    print(P,Q)
    plotting_load_flow(P,Q)