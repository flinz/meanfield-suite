from brian2 import *

 20000 # neurons
C_i = int(C * 0.2) # inibitory
C_e = int(C * 0.2) # excitatory
tau_0 = 2*ms

tau = 10*ms
V_reset = -60

eqs_i = '''
dv/dt = (i - v) / tau
di/dt = - i : amp
i += 1 : amp (event-driven)
'''

eqs_j = '''
dv/dt = (i - v) / tau
di/dt = - i : amp
i += 1 : amp (event-driven)
'''

#dv/dt  = (ge + gi - (v + 49*mV)) / (20*ms) : volt
#dge/dt = - ge / (5 * ms)                : volt
#dgi/dt = - gi / (10 * ms)               : volt

P = NeuronGroup(4000, eqs, method='linear', threshold='v>-50*mV', reset='v=-60*mV')
P.v = -60*mV
Pe = P[:3200]
Pi = P[3200:]

Ce = Synapses(Pe, P, on_pre='ge+=1.62*mV')
Ce.connect(p=0.02)

Ci = Synapses(Pi, P, on_pre='gi-=9*mV')
Ci.connect(p=0.02)

M = SpikeMonitor(P)

run(1*second)

plot(M.t/ms, M.i, '.')
show()
