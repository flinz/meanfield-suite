from brian2 import StateMonitor
from brian2.units import *

from meanfield.parameters import NP
from meanfield.parameters import SP
from meanfield.populations.MFLinearPop import MFLinearPop
from meanfield.populations.MFPoissonSource import MFPoissonSource
from meanfield.sources.MFLinear3TSSource import MFLinear3TSSource

params_pop = {
    NP.GAMMA: 0.280112,
    NP.BETA: 0.062,
    NP.GM: 25. * nS,
    NP.CM: 0.5 * nF,  # * 1e3,
    NP.VL: -70. * mV,
    NP.VTHR: -50. * mV,
    NP.VRES: -55. * mV,
    NP.TAU_RP: 2. * ms
}

params_source = {
    SP.GM: 0 * siemens,
    SP.VE: 0 * volt,
    SP.TAU: 10 * ms,
}

n = 100

poisson = MFPoissonSource('poisson', n, n * 10 * Hz)
pop = MFLinearPop('pop', n, {
    NP.GM: 10 * nsiemens,
    NP.VL: 0 * mV,
    NP.CM: 5 * nfarad,
    NP.VTHR: 0 * mV,
    NP.VRES: 0 * mV,
    NP.TAU_RP: 15 * ms
})
syn = MFLinear3TSSource('syn', pop, {
    SP.GM: 10 * nsiemens,
    SP.VREV: 0 * volt,
    SP.TAU: 0 * ms,
    SP.TAU_RISE: 10 * ms,
    SP.TAU_D1: 20 * ms,
    SP.TAU_D2: 30 * ms,
    SP.ALPHA: 1
}, poisson)

print(syn.b2_syn)
print([syn.post_variable_name_1, syn.post_variable_name_2, syn.post_variable_name_3])

#m = StateMonitor(syn.b2_syn, [syn.post_variable_name_1, syn.post_variable_name_2, syn.post_variable_name_3], record=True)
m1 = StateMonitor(syn.b2_syn, syn.post_variable_name_1, record=True)
m2 = StateMonitor(syn.b2_syn, syn.post_variable_name_2, record=True)
m3 = StateMonitor(syn.b2_syn, syn.post_variable_name_3, record=True)

