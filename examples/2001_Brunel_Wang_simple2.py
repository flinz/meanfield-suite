from brian2 import *
from meanfield.parameters import PP
from meanfield.parameters import IP

from meanfield.populations.MFLinearPop import MFLinearPop
from meanfield.solvers.MFSolver import MFSolverRatesVoltages
from meanfield.MFSystem import MFSystem
from meanfield.parameters import Connection
from meanfield.inputs.MFLinearInput import MFLinearInput
from meanfield.inputs.MFStaticInput import MFStaticInput
from meanfield.utils import brian2_introspect

BrianLogger.log_level_debug()
set_device('cpp_standalone')

# neurons
N = 1000
N_I = int(N * 0.2)  # interneurons

# voltage
V_L = -70. * mV  # resting
V_thr = -50. * mV
V_reset = -55. * mV
V_E = 0 * mV
V_I = -70. * mV

# membrane capacitance
C_m_I = 0.2 * nF

# membrane leak
g_m_I = 20. * nS

# refactorty period
tau_rp_E = 2. * ms
tau_rp_I = 1. * ms

# external stimuli
rate = 3 * Hz
C_ext = 800

# synapses
C_I = N_I

# AMPA (excitatory)
g_AMPA_ext_E = 2.08 * nS
g_AMPA_ext_I = 1.62 * nS
tau_AMPA = 2. * ms

# GABAergic (inhibitory)
g_GABA_E = 0*1.25 * nS * 200. / N_I
g_GABA_I = 0.01*0.973 * nS * 200. / N_I
tau_GABA = 10. * ms


I_params = {
    PP.GM: g_m_I,
    PP.CM: C_m_I,
    PP.VL: V_L,
    PP.VTHR: V_thr,
    PP.VRES: V_reset,
    PP.TAU_RP: tau_rp_I,
}

# GABA
pop_i = MFLinearPop("I", N_I, I_params)


source_i_noise = MFStaticInput("I_noise", pop_i, C_ext, rate, {
    IP.GM: g_AMPA_ext_I,
    IP.VREV: 0 * mvolt,
    IP.TAU: tau_AMPA,
})
pop_i.add_noise(source_i_noise)

# I->I GABA
source_ii_gaba = MFLinearInput('II GABA', pop_i, {
    IP.GM: g_GABA_I,
    IP.VREV: -70 * mvolt,
    IP.TAU: tau_GABA,
}, pop_i, Connection.all_to_others())


system = MFSystem("Brunel Wang simplified")
system.pops += [pop_i]

solver = MFSolverRatesVoltages(system, solver='mse')
#solver.run()

sp3 = SpikeMonitor(pop_i.brian2[:10])
rate3 = PopulationRateMonitor(pop_i.brian2)
s = StateMonitor(pop_i.brian2, ['v'], record=[100])

net = Network()
net.add(pop_i.brian2)
net.add(source_i_noise.brian2)
net.add(source_ii_gaba.brian2)
net.add(sp3)
net.add(rate3)
net.add(s)

net.run(2000 * ms)

brian2_introspect(net, globals())
#print()
#system.print_introspect()


if True:

    suptitle('new')

    subplot(211)
    title('rates')
    plot(rate3.t / ms, rate3.smooth_rate(width=10 * ms) / Hz, label='interneuron')
    legend()

    subplot(212)
    title('spikes')
    plot(sp3.t / ms, sp3.i, '.', markersize=5, label='interneuron')
    legend()

    show()


