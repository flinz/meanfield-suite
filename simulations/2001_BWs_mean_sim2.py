from brian2 import *

from MFParams import MFParams
from MFPop import MFLinearPop
from MFSolver import MFSolverRatesVoltages
from MFSource import MFSource
from MFSystem import MFSystem

BrianLogger.log_level_debug()

# neurons
N = 1000
N_E = int(N * 0.8)  # pyramidal neurons
N_I = int(N * 0.2)  # interneurons

# voltage
V_L = -70.  # resting
V_thr = -50.
V_reset = -55.
V_E = 0.
V_I = -70.

# membrane capacitance
C_m_E = 0.5
C_m_I = 0.2

# membrane leak
g_m_E = 25. * nS
g_m_I = 20. * nS

# refactorty period
tau_rp_E = 2.
tau_rp_I = 1.

# external stimuli
rate = 0.003
C_ext = 800

# synapses
C_E = N_E
C_I = N_I

# AMPA (excitatory)
g_AMPA_ext_E = 2.08
g_AMPA_ext_I = 1.62
g_AMPA_rec_E =  0 * 0.104 * 800. / N_E
g_AMPA_rec_I =  0 * 0.081 * 800. / N_E
tau_AMPA = 2.

# NMDA (excitatory)
g_NMDA_E = 0.01 * 0.327 * 800. / N_E
g_NMDA_I = 0.01 * 0.258 * 800. / N_E
tau_NMDA_rise = 2.
tau_NMDA_decay = 100.
alpha = 0. # 0.5 / ms
Mg2 = 0. # 1.

# GABAergic (inhibitory)
g_GABA_E = 0.5 * 1.25 * 200. / N_I
g_GABA_I = 0.1 * 0.973 * 200. / N_I
tau_GABA = 10.

# subpopulations
f = 0.1
p = 1
N_sub = int(N_E * f)
N_non = int(N_E * (1. - f * p))
w_plus = 1
w_minus = 1. - f * (w_plus - 1.) / (1. - f)

def mean():

    E_params = MFParams(
        gamma=0.280112,
        beta=0.062,
        g_L=g_m_E,
        C_m=C_m_E * 1e3,
        E_L=V_L,
        V_th=V_thr,
        V_reset=V_reset,
        t_ref=tau_rp_E
    )

    I_params = MFParams(
        gamma=0.280112,
        beta=0.062,
        g_L=g_m_I,
        C_m=C_m_I * 1e3,
        E_L=V_L,
        V_th=V_thr,
        V_reset=V_reset,
        t_ref=tau_rp_I
    )

    nu_e = 0.01
    nu_i = 0.01

    system = MFSystem("Brunel Wang simplified")

    pop_e1 = MFLinearPop("E", N_non, E_params)
    pop_e1.rate_ms = nu_e
    pop_e1.v_mean = -52.

    pop_e2 = MFLinearPop("Edown", N_sub, E_params)
    pop_e2.rate_ms = nu_e
    pop_e2.v_mean = -52.

    pop_i = MFLinearPop("I", N_I, I_params)
    pop_i.rate_ms = nu_i
    pop_i.v_mean = -52.

    system.pops += [pop_e1, pop_e2, pop_i]

    # noise pops
    source_e_noise1 = MFSource("E_noise1", pop_e1)
    source_e_noise1.noise_tau = tau_AMPA
    source_e_noise1.g_base = g_AMPA_ext_E
    source_e_noise1.g_dyn = lambda: rate * C_ext * tau_AMPA
    pop_e1.noise = source_e_noise1

    source_e_noise2 = MFSource("E_noise2", pop_e2)
    source_e_noise2.noise_tau = tau_AMPA
    source_e_noise2.g_base = g_AMPA_ext_E
    source_e_noise2.g_dyn = lambda: rate * C_ext * tau_AMPA
    pop_e2.noise = source_e_noise2

    source_i_noise = MFSource("I_noise", pop_i)
    source_i_noise.noise_tau = tau_AMPA
    source_i_noise.g_base = g_AMPA_ext_I
    source_i_noise.g_dyn = lambda: rate * C_ext * tau_AMPA
    pop_i.noise = source_i_noise

    # E->E NMDA
    source_ee_nmda1 = MFSource('EE NMDA 1', pop_e1)
    #source_ee_nmda1.is_nmda = True
    source_ee_nmda1.g_base = g_NMDA_E
    source_ee_nmda1.g_dyn = lambda: (
        pop_e1.n * f * w_plus * pop_e1.rate_ms * tau_NMDA_decay +
        pop_e2.n * (1. - f) * w_minus * pop_e2.rate_ms * tau_NMDA_decay
    )

    source_ee_nmda2 = MFSource('EE NMDA 2', pop_e2)
    #source_ee_nmda2.is_nmda = True
    source_ee_nmda2.g_base = g_NMDA_E
    source_ee_nmda2.g_dyn = lambda: (
        pop_e1.n * f * w_minus * pop_e1.rate_ms * tau_NMDA_decay +
        pop_e2.n * (1. - f) * w_plus * pop_e2.rate_ms * tau_NMDA_decay
    )

    # E->E AMPA
    source_ee_ampa1 = MFSource('EE AMPA 1', pop_e1)
    source_ee_ampa1.g_base = g_AMPA_rec_E
    source_ee_ampa1.g_dyn = lambda: (
        pop_e1.n * f * w_plus * pop_e1.rate_ms * tau_AMPA +
        pop_e2.n * (1. - f) * w_minus * pop_e2.rate_ms * tau_AMPA
    )

    source_ee_ampa2 = MFSource('EE AMPA 2', pop_e2)
    source_ee_ampa2.g_base = g_AMPA_rec_E
    source_ee_ampa2.g_dyn = lambda: (
        pop_e1.n * f * w_minus * pop_e1.rate_ms * tau_AMPA +
        pop_e2.n * (1. - f) * w_plus * pop_e2.rate_ms * tau_AMPA
    )

    # E->I NMDA
    source_ie_nmda = MFSource('EI NMDA', pop_i)
    #source_ie_nmda.is_nmda = True
    source_ie_nmda.g_base = g_NMDA_I
    source_ie_nmda.g_dyn = lambda: (
        pop_e1.n * f * pop_e1.rate_ms * tau_NMDA_decay +
        pop_e2.n * (1. - f) * pop_e2.rate_ms * tau_NMDA_decay
    )

    # E->I AMPA
    source_ie_ampa = MFSource('EI AMPA', pop_i)
    source_ie_ampa.g_base = g_AMPA_rec_E
    source_ie_ampa.g_dyn = lambda: (
        pop_e1.n * f * pop_e1.rate_ms * tau_AMPA +
        pop_e2.n * (1. - f) * pop_e2.rate_ms * tau_AMPA
    )

    # I->I GABA
    source_ii_gaba = MFSource('II GABA', pop_i)
    source_ii_gaba.E_rev = -70.
    source_ii_gaba.g_base = g_GABA_I
    source_ii_gaba.g_dyn = lambda: pop_i.n * pop_i.rate_ms * tau_GABA

    # I->E GABA
    source_ie_gaba1 = MFSource('IE GABA 1', pop_e1)
    source_ie_gaba1.E_rev = -70.
    source_ie_gaba1.g_base = g_GABA_E
    source_ie_gaba1.g_dyn = lambda: pop_i.n * pop_i.rate_ms * tau_GABA

    source_ie_gaba2 = MFSource('IE GABA 2', pop_e2)
    source_ie_gaba2.E_rev = -70.
    source_ie_gaba2.g_base = g_GABA_E
    source_ie_gaba2.g_dyn = lambda: pop_i.n * pop_i.rate_ms * tau_GABA

    solver = MFSolverRatesVoltages(system, solver='gradient')
    solver.run()


def sim():

    eqs_E = '''
    dv / dt = (- g_m_E * nS * (v - V_L * mV) - I_syn) / (C_m_E * nF) : volt (unless refractory)

    I_syn = I_AMPA_ext + I_AMPA_rec + I_NMDA_rec + I_GABA_rec : amp

    I_AMPA_ext = g_AMPA_ext_E * nS * (v - V_E * mV) * s_AMPA_ext : amp
    I_AMPA_rec = g_AMPA_rec_E * nS * (v - V_E * mV) * s_AMPA : amp
    ds_AMPA_ext / dt = - s_AMPA_ext / (tau_AMPA * ms) : 1
    ds_AMPA / dt = - s_AMPA / (tau_AMPA * ms) : 1

    I_NMDA_rec = g_NMDA_E * nS * (v - V_E * mV) * s_NMDA : amp
    ds_NMDA / dt = - s_NMDA / (tau_NMDA_decay * ms) : 1

    I_GABA_rec = g_GABA_E * nS * (v - V_I * mV) * s_GABA : amp
    ds_GABA / dt = - s_GABA / (tau_GABA * ms) : 1
    '''

    eqs_I = '''
    dv / dt = (- g_m_I * nS * (v - V_L * mV) - I_syn) / (C_m_I * nF) : volt (unless refractory)

    I_syn = I_AMPA_ext + I_AMPA_rec + I_NMDA_rec + I_GABA_rec : amp

    I_AMPA_ext = g_AMPA_ext_I * nS * (v - V_E * mV) * s_AMPA_ext : amp
    I_AMPA_rec = g_AMPA_rec_I * nS * (v - V_E * mV) * s_AMPA : amp
    ds_AMPA_ext / dt = - s_AMPA_ext / (tau_AMPA * ms) : 1
    ds_AMPA / dt = - s_AMPA / (tau_AMPA * ms) : 1

    I_NMDA_rec = g_NMDA_I * nS * (v - V_E * mV) * s_NMDA : amp
    ds_NMDA / dt = - s_NMDA / (tau_NMDA_decay * ms) : 1

    I_GABA_rec = g_GABA_I * nS * (v - V_I * mV) * s_GABA : amp
    ds_GABA / dt = - s_GABA / (tau_GABA * ms) : 1
    '''

    P_E = NeuronGroup(N_E, eqs_E, method='euler', threshold='v > V_thr * mV', reset='v = V_reset * mV', refractory=tau_rp_E * ms)
    P_E.v = V_L * mV
    P_I = NeuronGroup(N_I, eqs_I, method='euler', threshold='v > V_thr * mV', reset='v = V_reset * mV', refractory=tau_rp_I * ms)
    P_I.v = V_L * mV

    eqs_glut = '''
    w : 1
    '''

    eqs_pre_glut = '''
    s_AMPA += w
    s_NMDA += w
    '''

    eqs_pre_gaba = '''
    s_GABA += 1
    '''

    eqs_pre_ext = '''
    s_AMPA_ext += 1
    '''

    # recurrent E to E
    C_E_E = Synapses(P_E, P_E, method='euler', model=eqs_glut, on_pre=eqs_pre_glut)
    C_E_E.connect('i != j')
    C_E_E.w[:] = 1

    for pi in range(N_non, N_non + p * N_sub, N_sub):

        # internal other subpopulation to current nonselective
        C_E_E.w[C_E_E.indices[:, pi:pi + N_sub]] = w_minus

        # internal current subpopulation to current subpopulation
        C_E_E.w[C_E_E.indices[pi:pi + N_sub, pi:pi + N_sub]] = w_plus

    # E to I
    C_E_I = Synapses(P_E, P_I, method='euler', model=eqs_glut, on_pre=eqs_pre_glut)
    C_E_I.connect()
    C_E_I.w[:] = 1

    # recurrent I to I
    C_I_I = Synapses(P_I, P_I, on_pre=eqs_pre_gaba)
    C_I_I.connect('i != j')

    # I to E
    C_I_E = Synapses(P_I, P_E, on_pre=eqs_pre_gaba)
    C_I_E.connect()

    # external
    C_P_E = PoissonInput(P_E, 's_AMPA_ext', C_ext, rate * Hz, 500)
    C_P_I = PoissonInput(P_I, 's_AMPA_ext', C_ext, rate * Hz, 500)

    # monitors
    sp_E = SpikeMonitor(P_E[:40])
    sp_I = SpikeMonitor(P_I[:40])
    r_E = PopulationRateMonitor(P_E)
    r_I = PopulationRateMonitor(P_I)

    run(3000 * ms, report='text', report_period=0.5 * second)

    subplot(311)
    plot(r_E.t / ms, r_E.smooth_rate(width=25 * ms) / Hz, label='pyramidal neuron')
    plot(r_I.t / ms, r_I.smooth_rate(width=25 * ms) / Hz, label='interneuron')
    legend()

    subplot(312)
    plot(sp_E.t / ms, sp_E.i, '.', markersize=5, label='nonselective')

    subplot(313)
    plot(sp_I.t / ms, sp_I.i, '.', markersize=5)

    show()

mean()
sim()
