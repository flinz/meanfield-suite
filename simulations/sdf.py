from brian2 import *

gm = 25 * nS
VL = -70 * mV
Cm = 0.5 * nF
VE = 0 * mV
g_nmda = 0.327 * nS
t_nmda_decay = 100 * ms
alpha_nmda = 0.5 * ms**-1
t_nmda_rise = 2 * ms


set_global_preferences(useweave=True, usecodegenweave=True,
                       usenewpropagate=True,
                       gcc_options=['-march=native'],
                       usecodegenstateupdate=True, usecodegen=True)

# The neuron
eqs = '''
dv/dt = (-gm*(v-VL) - I_syn)/Cm : volt
I_syn = g_nmda*(v-VE)*s_nmda/(1 + 0.001 * exp(-0.062 * v / volt) / 3.57) : amp
s_nmda : 1
'''

G = NeuronGroup(1000, model=eqs, threshold=-50 * mV, reset=-55 * mV,
                refractory= 2 * ms)
G.v = -70 * mV
G[0].v = -49 * mV # set the first neuron's potential to a suprathreshold value

# The NMDA synapses
nmda_eqs = '''
s_nmda_ in = WW*s_nmda : 1
ds_nmda/dt = -s_nmda/t_nmda_decay + alpha_nmda * x * (1 - s_nmda) : 1
dx/dt = -x / (t_nmda_rise) : 1
WW : 1
'''

syn = Synapses(G, G, model=nmda_eqs, pre='x += w')
G.s_nmda = syn.s_nmda_in

syn[:, :] = 'i != j' #all to all but not self-connections
syn.w[:, :] = 1
syn.WW[:, :] = 2.1

syn.delay[:, :] = 0.5*ms
mon = StateMonitor(G, 'v', record=True)
net = Network(G, syn, mon)

# Run and plot results
net.run(50 * msecond, report='text')
plot(mon.times / ms, mon[0] / mV)
plot(mon.times / ms, mon[1] / mV)

show()
