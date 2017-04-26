from brian2 import units

from MFParams import MFParams
from MFPop import MFLinearPop
from MFSolver import MFSolver, MFSolverRatesVoltages
from MFSource import MFSource, MFDynamicSource, MFStaticSource
from MFSystem import MFSystem
from brian2 import *

from params import NP, SP

BrianLogger.log_level_debug()

# neurons
N = 1000
N_E = int(N * 0.8)  # pyramidal neurons
N_I = int(N * 0.2)  # interneurons

# voltage
V_L = -70. * mV # resting
V_thr = -50. * mV
V_reset = -55. * mV
V_E = 0. * mV
V_I = -70.

# membrane capacitance
C_m_E = 0.5 * nF
C_m_I = 0.2 * nF

# membrane leak
g_m_E = 25. * nS
g_m_I = 20. * nS

# refactorty period
tau_rp_E = 2. * ms
tau_rp_I = 1. * ms

# external stimuli
rate = 3 * Hz
C_ext = 800

# synapses
C_E = N_E
C_I = N_I

# AMPA (excitatory)
g_AMPA_ext_E = 2.08 * nS
g_AMPA_ext_I = 1.62 * nS
g_AMPA_rec_E =  0 * 0.104 * 800. / N_E * nS
g_AMPA_rec_I =  0 * 0.081 * 800. / N_E * nS
tau_AMPA = 2. * ms

# NMDA (excitatory)
g_NMDA_E = 0 * 0.01 * 0.327 * 800. / N_E * nS
g_NMDA_I = 0 * 0.01 * 0.258 * 800. / N_E * nS
tau_NMDA_rise = 2. * ms
tau_NMDA_decay = 100. * ms
alpha = 0. # 0.5 / ms
Mg2 = 0. # 1.

# GABAergic (inhibitory)
g_GABA_E = 0 * 0.5 * 1.25 * 200. / N_I * nS
g_GABA_I = 0 * 0.1 * 0.973 * 200. / N_I * nS
tau_GABA = 10. * ms

# subpopulations
f = 0.1
p = 1
N_sub = int(N_E * f)
N_non = int(N_E * (1. - f * p))
w_plus = 1
w_minus = 1. - f * (w_plus - 1.) / (1. - f)

def mean():

    nu_e = 10. * Hz
    nu_i = 10. * Hz

    system = MFSystem("Brunel Wang simplified")

    e_params = {
        NP.GAMMA: 0.280112,
        NP.BETA: 0.062,
        NP.GM: g_m_E,
        NP.CM: C_m_E,# * 1e3,
        NP.VL: V_L,
        NP.VTHR: V_thr,
        NP.VRES: V_reset,
        NP.TAU_RP: tau_rp_E
    }

    pop_e1 = MFLinearPop("E", N_non, e_params)
    pop_e1.rate = nu_e
    pop_e1.v_mean = -52. * mV

    pop_e2 = MFLinearPop("Edown", N_sub, e_params)
    pop_e2.rate = nu_e
    pop_e2.v_mean = -52. * mV

    i_params = {
        NP.GAMMA: 0.280112,
        NP.BETA: 0.062,
        NP.GM: g_m_I,
        NP.CM: C_m_I,# * 1e3,
        NP.VL: V_L,
        NP.VTHR: V_thr,
        NP.VRES: V_reset,
        NP.TAU_RP: tau_rp_I
    }

    pop_i = MFLinearPop("I", N_I, i_params)
    pop_i.rate = nu_i
    pop_i.v_mean = -52. * mV

    system.pops += [pop_e1, pop_e2, pop_i]

    # noise pops
    source_e_noise1 = MFStaticSource("E_noise1", pop_e1, {
        SP.GM: g_AMPA_ext_E,
        SP.VE: 0. * mV,
        SP.TAU: tau_AMPA
    }, rate, C_ext)
    pop_e1.noise = source_e_noise1

    source_e_noise2 = MFStaticSource("E_noise2", pop_e2, {
        SP.GM: g_AMPA_ext_E,
        SP.VE: 0. * mV,
        SP.TAU: tau_AMPA
    }, rate, C_ext)
    pop_e2.noise = source_e_noise2

    source_i_noise = MFStaticSource("I_noise", pop_i, {
        SP.GM: g_AMPA_ext_I,
        SP.VE: 0. * mV,
        SP.TAU: tau_AMPA
    }, rate, C_ext)
    pop_i.noise = source_i_noise

    # E->E NMDA 1.1
    source_ee_nmda11 = MFDynamicSource('EE NMDA 11', pop_e1, {
        SP.GM: g_NMDA_E,
        SP.VE: 0. * mV,
        SP.TAU: tau_NMDA_decay,
        SP.W: w_plus,
        SP.FRAC: f
    }, from_pop=pop_e1)
    # source_ee_nmda1.is_nmda = True
    # E->E NMDA 1.2
    source_ee_nmda12 = MFDynamicSource('EE NMDA 12', pop_e1, {
        SP.GM: g_NMDA_E,
        SP.VE: 0. * mV,
        SP.TAU: tau_NMDA_decay,
        SP.W: w_minus,
        SP.FRAC: 1 - f
    }, from_pop=pop_e2)
    # source_ee_nmda1.is_nmda = True

    # E->E NDMA 2.1
    source_ee_nmda21 = MFDynamicSource('EE NMDA 21', pop_e2, {
        SP.GM: g_NMDA_E,
        SP.VE: 0. * mV,
        SP.TAU: tau_NMDA_decay,
        SP.W: w_minus,
        SP.FRAC: f
    }, from_pop=pop_e1)
    # source_ee_nmda2.is_nmda = True
    # E->E NDMA 2.2
    source_ee_nmda22 = MFDynamicSource('EE NMDA 22', pop_e2, {
        SP.GM: g_NMDA_E,
        SP.VE: 0. * mV,
        SP.TAU: tau_NMDA_decay,
        SP.W: w_plus,
        SP.FRAC: 1 - f
    }, from_pop=pop_e2)
    # source_ee_nmda2.is_nmda = True

    # E->E AMPA 1
    source_ee_ampa11 = MFDynamicSource('EE AMPA 11', pop_e1, {
        SP.GM: g_AMPA_rec_E,
        SP.VE: 0. * mV,
        SP.TAU: tau_AMPA,
        SP.W: w_plus,
        SP.FRAC: f
    }, from_pop=pop_e1)

    # E->E AMPA 2
    source_ee_ampa12 = MFDynamicSource('EE AMPA 12', pop_e1, {
        SP.GM: g_AMPA_rec_E,
        SP.VE: 0. * mV,
        SP.TAU: tau_AMPA,
        SP.W: w_minus,
        SP.FRAC: 1 - f
    }, from_pop=pop_e2)

    source_ee_ampa21 = MFDynamicSource('EE AMPA 21', pop_e2, {
        SP.GM: g_AMPA_rec_E,
        SP.VE: 0. * mV,
        SP.TAU: tau_AMPA,
        SP.W: w_minus,
        SP.FRAC: f
    }, from_pop=pop_e1)
    source_ee_ampa22 = MFDynamicSource('EE AMPA 22', pop_e2, {
        SP.GM: g_AMPA_rec_E,
        SP.VE: 0. * mV,
        SP.TAU: tau_AMPA,
        SP.W: w_plus,
        SP.FRAC: 1 - f
    }, from_pop=pop_e2)

    # E->I NMDA
    source_ie_nmda1 = MFDynamicSource('EI NMDA 1', pop_i, {
        SP.GM: g_NMDA_I,
        SP.VE: 0. * mV,
        SP.TAU: tau_NMDA_decay,
        SP.W: 1,
        SP.FRAC: f
    }, from_pop=pop_e1)
    # source_ie_nmda.is_nmda = True

    source_ie_nmda2 = MFDynamicSource('EI NMDA 2', pop_i, {
        SP.GM: g_NMDA_I,
        SP.VE: 0. * mV,
        SP.TAU: tau_NMDA_decay,
        SP.W: 1,
        SP.FRAC: 1 - f
    }, from_pop=pop_e2)
    # source_ie_nmda.is_nmda = True

    # E->I AMPA
    source_ie_ampa1 = MFDynamicSource('EI AMPA 1', pop_i, {
        SP.GM: g_AMPA_rec_E,
        SP.VE: 0. * mV,
        SP.TAU: tau_AMPA,
        SP.W: 1,
        SP.FRAC: f
    }, from_pop=pop_e1)
    source_ie_ampa2 = MFDynamicSource('EI AMPA 2', pop_i, {
        SP.GM: g_AMPA_rec_E,
        SP.VE: 0. * mV,
        SP.TAU: tau_AMPA,
        SP.W: 1,
        SP.FRAC: 1 - f
    }, from_pop=pop_e2)

    # I->I GABA
    source_ii_gaba = MFDynamicSource('II GABA', pop_i, {
        SP.GM: g_GABA_I,
        SP.VE: -70. * mV,
        SP.TAU: tau_GABA,
        SP.W: 1,
        SP.FRAC: 1
    }, from_pop=pop_i)

    # I->E GABA
    source_ie_gaba1 = MFDynamicSource('IE GABA 1', pop_e1, {
        SP.GM: g_GABA_E,
        SP.VE: -70. * mV,
        SP.TAU: tau_GABA,
        SP.W: 1,
        SP.FRAC: 1
    }, from_pop=pop_i)

    source_ie_gaba2 = MFDynamicSource('IE GABA 2', pop_e2, {
        SP.GM: g_GABA_E,
        SP.VE: -70. * mV,
        SP.TAU: tau_GABA,
        SP.W: 1,
        SP.FRAC: 1
    }, from_pop=pop_i)

    solver = MFSolverRatesVoltages(system, solver='mse')
    solver.run()

    print(pop_e1.brian2_model())
    print(pop_e2.brian2_model())
    print(pop_i.brian2_model())

    net = Network()
    net.add(pop_e1.brian2)
    net.add(pop_e2.brian2)
    net.add(pop_i.brian2)

    for p in system.pops:
        for i, s in enumerate(p.sources):
            print(s.brian2)
            if s.brian2:
                net.add(s.brian2)

    n1 = PoissonInput(pop_e1.brian2, 's_E_noise1', C_ext, rate, 1)
    n2 = PoissonInput(pop_e2.brian2, 's_E_noise2', C_ext, rate, 1)
    n3 = PoissonInput(pop_i.brian2, 's_I_noise', C_ext, rate, 1)
    net.add(n1)
    net.add(n2)
    net.add(n3)

    sp_E = SpikeMonitor(pop_e1.brian2[:40])
    sp_I = SpikeMonitor(pop_i.brian2[:40])
    r_E = PopulationRateMonitor(pop_e1.brian2)
    r_I = PopulationRateMonitor(pop_i.brian2)

    net.add(sp_E)
    net.add(sp_I)
    net.add(r_E)
    net.add(r_I)

    net.run(3000 * ms, report='text', report_period=0.5 * second)

    subplot(311)
    plot(r_E.t / ms, r_E.smooth_rate(width=25 * ms) / Hz, label='pyramidal neuron')
    plot(r_I.t / ms, r_I.smooth_rate(width=25 * ms) / Hz, label='interneuron')
    legend()

    subplot(312)
    plot(sp_E.t / ms, sp_E.i, '.', markersize=5, label='nonselective')

    subplot(313)
    plot(sp_I.t / ms, sp_I.i, '.', markersize=5)

    show()

# def sim():
#
#     eqs_E = '''
#     dv / dt = (- g_m_E * nS * (v - V_L * mV) - I_syn) / (C_m_E * nF) : volt (unless refractory)
#
#     I_syn = I_AMPA_ext + I_AMPA_rec + I_NMDA_rec + I_GABA_rec : amp
#
#     I_AMPA_ext = g_AMPA_ext_E * nS * (v - V_E * mV) * s_AMPA_ext : amp
#     I_AMPA_rec = g_AMPA_rec_E * nS * (v - V_E * mV) * s_AMPA : amp
#     ds_AMPA_ext / dt = - s_AMPA_ext / (tau_AMPA * ms) : 1
#     ds_AMPA / dt = - s_AMPA / (tau_AMPA * ms) : 1
#
#     I_NMDA_rec = g_NMDA_E * nS * (v - V_E * mV) * s_NMDA : amp
#     ds_NMDA / dt = - s_NMDA / (tau_NMDA_decay * ms) : 1
#
#     I_GABA_rec = g_GABA_E * nS * (v - V_I * mV) * s_GABA : amp
#     ds_GABA / dt = - s_GABA / (tau_GABA * ms) : 1
#     '''
#
#     eqs_I = '''
#     dv / dt = (- g_m_I * nS * (v - V_L * mV) - I_syn) / (C_m_I * nF) : volt (unless refractory)
#
#     I_syn = I_AMPA_ext + I_AMPA_rec + I_NMDA_rec + I_GABA_rec : amp
#
#     I_AMPA_ext = g_AMPA_ext_I * nS * (v - V_E * mV) * s_AMPA_ext : amp
#     I_AMPA_rec = g_AMPA_rec_I * nS * (v - V_E * mV) * s_AMPA : amp
#     ds_AMPA_ext / dt = - s_AMPA_ext / (tau_AMPA * ms) : 1
#     ds_AMPA / dt = - s_AMPA / (tau_AMPA * ms) : 1
#
#     I_NMDA_rec = g_NMDA_I * nS * (v - V_E * mV) * s_NMDA : amp
#     ds_NMDA / dt = - s_NMDA / (tau_NMDA_decay * ms) : 1
#
#     I_GABA_rec = g_GABA_I * nS * (v - V_I * mV) * s_GABA : amp
#     ds_GABA / dt = - s_GABA / (tau_GABA * ms) : 1
#     '''
#
#     P_E = NeuronGroup(N_E, eqs_E, method='euler', threshold='v > V_thr * mV', reset='v = V_reset * mV', refractory=tau_rp_E * ms)
#     P_E.v = V_L * mV
#     P_I = NeuronGroup(N_I, eqs_I, method='euler', threshold='v > V_thr * mV', reset='v = V_reset * mV', refractory=tau_rp_I * ms)
#     P_I.v = V_L * mV
#
#     eqs_glut = '''
#     w : 1
#     '''
#
#     eqs_pre_glut = '''
#     s_AMPA += w
#     s_NMDA += w
#     '''
#
#     eqs_pre_gaba = '''
#     s_GABA += 1
#     '''
#
#     eqs_pre_ext = '''
#     s_AMPA_ext += 1
#     '''
#
#     # recurrent E to E
#     C_E_E = Synapses(P_E, P_E, method='euler', model=eqs_glut, on_pre=eqs_pre_glut)
#     C_E_E.connect('i != j')
#     C_E_E.w[:] = 1
#
#     for pi in range(N_non, N_non + p * N_sub, N_sub):
#
#         # internal other subpopulation to current nonselective
#         C_E_E.w[C_E_E.indices[:, pi:pi + N_sub]] = w_minus
#
#         # internal current subpopulation to current subpopulation
#         C_E_E.w[C_E_E.indices[pi:pi + N_sub, pi:pi + N_sub]] = w_plus
#
#     # E to I
#     C_E_I = Synapses(P_E, P_I, method='euler', model=eqs_glut, on_pre=eqs_pre_glut)
#     C_E_I.connect()
#     C_E_I.w[:] = 1
#
#     # recurrent I to I
#     C_I_I = Synapses(P_I, P_I, on_pre=eqs_pre_gaba)
#     C_I_I.connect('i != j')
#
#     # I to E
#     C_I_E = Synapses(P_I, P_E, on_pre=eqs_pre_gaba)
#     C_I_E.connect()
#
#     # external
#     C_P_E = PoissonInput(P_E, 's_AMPA_ext', C_ext, rate * Hz, 500)
#     C_P_I = PoissonInput(P_I, 's_AMPA_ext', C_ext, rate * Hz, 500)
#
#     # monitors
#     sp_E = SpikeMonitor(P_E[:40])
#     sp_I = SpikeMonitor(P_I[:40])
#     r_E = PopulationRateMonitor(P_E)
#     r_I = PopulationRateMonitor(P_I)
#
#     run(3000 * ms, report='text', report_period=0.5 * second)
#
#     subplot(311)
#     #plot(r_E.t / ms, r_E.smooth_rate(width=25 * ms) / Hz, label='pyramidal neuron')
#     #plot(r_I.t / ms, r_I.smooth_rate(width=25 * ms) / Hz, label='interneuron')
#     legend()
#
#     subplot(312)
#     plot(sp_E.t / ms, sp_E.i, '.', markersize=5, label='nonselective')
#
#     subplot(313)
#     plot(sp_I.t / ms, sp_I.i, '.', markersize=5)
#
#     show()

mean()
