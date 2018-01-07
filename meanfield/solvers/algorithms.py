
import matplotlib.pyplot as plt
import numpy as np


def custom_gradient_solver(mfstate, p_0, plotting: bool=False, dt: float = .1, tmax: float = 30.):
    """Simple gradient descent along the error."""

    t = 0.
    state = np.array(p_0)
    states = [state]

    while t < tmax:
        state -= [dt * v for v in mfstate(state)]
        t += dt
        states.append(list(state))

    states = np.array(states).T

    if plotting:
        nspl = states.shape[0]
        for i in range(nspl):
            plt.subplot(nspl, 1, i + 1)
            plt.plot(states[i, :])
        plt.show()

    # minimal solution object to return
    class sol:
        x = state
        fun = mfstate.error

    return sol()
