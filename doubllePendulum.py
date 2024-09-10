import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from collections import deque
import numpy as np
from scipy.integrate import odeint

VECTOR_NUM = 64
ROUND_FIX = 10
DT = 0.001


def Ek(V, m):
    return V**2 * m * 0.5


def U(h, m, g):
    return h * m * g


# The function to pass to the odeint
def module(params, t, u, L1, L2, g, k=None, func=None):
    # params init
    O1, O2, w1, w2 = params
    if func != None:  # call the special condition function
        u, L1, L2, O1, O2, w1, w2 = func(t, k, u, L1, L2, O1, O2, w1, w2)
    c, s = np.cos(O1 - O2), np.sin(O1 - O2)  # shortcut because this is used a lot
    w1dot = (
        (
            u * g * np.sin(O2) * c
            - u * s * (L1 * w1**2 * c + L2 * w2**2)
            - (1 + u) * g * np.sin(O1)
        )
        / L1
        / (1 + u * s**2)
    )
    w2dot = (
        (
            (1 + u) * (L1 * w1**2 * s - g * np.sin(O2) + g * np.sin(O1) * c)
            + u * L2 * w2**2 * s * c
        )
        / L2
        / (1 + u * s**2)
    )
    return (
        round(w1, ROUND_FIX),
        round(w2, ROUND_FIX),
        round(w1dot, ROUND_FIX),
        round(w2dot, ROUND_FIX),
    )


def cartesian_convert_pos(theta, L):
    return np.vstack((np.sin(theta), -np.cos(theta))) * L


def cartesian_convert_vel(theta, dthetadt, L):
    return np.vstack((np.cos(theta) * dthetadt, np.sin(theta) * dthetadt)) * L


def vector_size(x, y):
    return (x**2 + y**2) ** 0.5


def sim(
    u, L1, L2, O10, O20, w10=0, w20=0, g=9.8, t_max=10, delta_t=DT, k=None, func=None
):
    """simulate a double pendulum accation with the given parameters

    Parameters
    ----------
    u : float
        masses' ratio (m2/m1)
    L1 : float
        length of the first pendulum
    L2 : float
        length of the second pendulum
    O10 : float
        starting angle of the first pendulum
    O20 : float
        starting angle of the second pendulum
    w10 : float
        starting angular velocity of the first pendulum
    w20 : float
        starting angular velocity of the second pendulum
    g : float
        gravitational acceleration
    t_max : float
        maximum time to simulate
    delta_t : float
        time step
    k : float
        anouther parameter that could be passed to func
    func : function
        a function that could be passed to change the parameters dynamically at runtime

    Return
    -------
    t : list
        the time arrat
    ang : list
        array of the anngular properties of the results ( angle1, angle2, angular velocity1, angular velocity2)
    cart : list
        array of the cartesian properties of the results all of them are vectors containers ( pos1, pos2, velocity1, velocity2)
    E : list
        array of the energy properties of the results ( kinetic energy1, kinetic energy2, potential energy1, potential energy2)
    """
    t = np.arange(0, t_max, delta_t)  # set a the range of the calculations to 0-10
    q_not_sorted = odeint(
        module, [O10, O20, w10, w20], t, args=(u, L1, L2, g, k, func)
    )  # run the simulation
    q = q_not_sorted.T
    pos1 = cartesian_convert_pos(q[0], L1)
    pos2 = pos1 + cartesian_convert_pos(q[1], L2)
    vel1 = cartesian_convert_vel(q[0], q[2], L1)
    vel2 = vel1 + cartesian_convert_vel(q[1], q[3], L1)
    Ek1 = Ek(vector_size(*vel1), 1)
    Ek2 = Ek(vector_size(*vel2), u)
    U1 = U(pos1[1] + L1, 1, g)
    U2 = U(pos2[1] + L1 + L2, u, g)
    return t, q, [pos1, pos2, vel1, vel2], [Ek1, Ek2, U1, U2]


def plott(
    t,
    ang,
    cart,
    energy,
    graphs=None,
    ext_label="",
    c1="r",
    c2="b",
    c3="#7f007f",
    vectors=True,
):
    """plot the results of the simulation

    Parameters
    ----------
    ~ all the return values of sim ~
    graphs : list
        optional for appending of existing graphs instead of creating new ones
    ext_label : str
        optional label to add for the graphs
    c[1-3] : str
        color of the graphs
    vector : bool
        show the velocity vectors on the path or not
    """
    if graphs == None:  # init the subplots
        graphs = [plt.subplot(2, 2, index) for index in range(1, 5)]
    # plot the x by time
    graphs[0].set_title("X By Time", color="r")
    graphs[0].plot(t, cart[0][0], label=ext_label + " X1", color=c1)
    graphs[0].plot(t, cart[1][0], label=ext_label + " X2", color=c2)
    graphs[0].set_xlabel("Time[sec]", color="r")
    graphs[0].set_ylabel("X[m]", color="r")
    # plot the path of the penduloms in the space
    graphs[1].set_title("Path In XY Space", color="r")
    graphs[1].plot(*cart[0], label=ext_label + " Path 1", color=c1)
    graphs[1].plot(*cart[1], label=ext_label + " Path 2", color=c2)
    # plot velocity vectors on the path
    if vectors:
        graphs[1].quiver(
            *[i[:: len(i) // VECTOR_NUM] for i in cart[0]],
            *[i[:: len(i) // VECTOR_NUM] for i in cart[2]],
            units="xy",
            color=c1,
            width=0.025
        )
        graphs[1].quiver(
            *[i[:: len(i) // VECTOR_NUM] for i in cart[1]],
            *[i[:: len(i) // VECTOR_NUM] for i in cart[3]],
            units="xy",
            color=c2,
            width=0.025
        )
    graphs[1].set_xlabel("X[m]", color="r")
    graphs[1].set_ylabel("Y[m]", color="r")
    # plot the angle by time
    graphs[2].set_title("Angles By Time", color="r")
    graphs[2].plot(t, ang[0], color=c1, label=ext_label + " θ1")
    graphs[2].plot(t, ang[1], color=c2, label=ext_label + " θ2")
    graphs[2].set_xlabel("Time[sec]", color="r")
    graphs[2].set_ylabel("Angle[rad]", color="r")
    # plot the kinetick energy by time
    graphs[3].set_title("Kinetic Energies By Time", color="r")
    graphs[3].plot(t, energy[0], lw=0.6, label=ext_label + " Ek 1", color=c1)
    graphs[3].plot(t, energy[1], lw=0.6, label=ext_label + " Ek 2", color=c2)
    graphs[3].plot(
        t, energy[0] + energy[1], label=ext_label + " Total Ek", color=c3, lw=0.8
    )
    graphs[3].set_xlabel("t[s]", color="r")
    graphs[3].set_ylabel("E[j]", color="r")
    for i in graphs:  # Turn on the lengend and grid in each graph and
        i.legend(loc="lower right", shadow=True)
        i.grid(visible=True)
    return graphs


def animation(t, pos1, pos2, size, name, fps=30):
    """Create an animation of the simulation and save it as an mp4 file

    Parameters
    ----------
    t : list
        the array of the time
    pos[1-2] : list
        the positions of the pendulums
    size : int
        half of the size of the plot (recommended L1 + L2)
    name : str
        the name of the animation mp4 file
    fps : int
        the frames per second of the animation
    """
    delta_t = t[1] - t[0]
    douration = int(
        round(t[-1], 0)
    )  # auto set the douration of the video to be the sim time
    # init the plots
    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(autoscale_on=False, xlim=(-size, size), ylim=(-size, size))
    ax.set_aspect("equal")
    ax.grid()
    (line,) = ax.plot([], [], "o-", lw=5)
    (trace,) = ax.plot([], [], ".-", lw=1, ms=2)
    time_template = "time = %.1fs"
    time_text = ax.text(0.05, 0.9, "", transform=ax.transAxes)
    history_x, history_y = deque(maxlen=douration * fps), deque(maxlen=douration * fps)

    # a function that sets every frame of the animation.
    def animate(i):
        thisx = [
            0,
            pos1[0][:: int(1 / fps / delta_t)][i],
            pos2[0][:: int(1 / fps / delta_t)][i],
        ]
        thisy = [
            0,
            pos1[1][:: int(1 / fps / delta_t)][i],
            pos2[1][:: int(1 / fps / delta_t)][i],
        ]
        if i == 0:
            history_x.clear()
            history_y.clear()
        history_x.appendleft(thisx[2])
        history_y.appendleft(thisy[2])
        # sets the data of the frame's plots
        line.set_data(thisx, thisy)
        trace.set_data(history_x, history_y)
        time_text.set_text(time_template % (t[:: int(1 / fps / delta_t)][i]))
        return line, trace, time_text

    # create the animation
    ani = FuncAnimation(fig, animate, fps * douration, interval=1 / fps, blit=True)
    # save the animation

    ani.save(filename=("%s.gif" % name), writer="ffmpeg", fps=fps)
    plt.show()


def changeMass(_t, _k, _u, _L1, _L2, _O1, _O2, _w1, _w2):
    return _u + _k * _t, _L1, _L2, _O1, _O2, _w1, _w2


def changeLengh(_t, _k, _u, _L1, _L2, _O1, _O2, _w1, _w2):
    return _u, _L1 + _k * _t, _L2 + _k * _t, _O1, _O2, _w1, _w2


if __name__ == "__main__":
    u = 2
    l1 = 1
    l2 = 2
    o1 = 1.5
    o2 = 0.5
    output = sim(u, l1, l2, o1, o2)
    plott(*output)
    plt.show()
    animation(output[0], output[2][0], output[2][1], l1 + l2, "anime/ere")
