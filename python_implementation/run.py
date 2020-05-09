from crn import *
a, a1, a2, b, c, t, z = species("A A1 A2 B C T Z")

sys = CRN(
    a >> a1 + a2,
    a1 + b >> t,
    c >> z,
    z + t >> 0,
    name="presentation")

sim = sys.stoch_simulate({a: 25, b: 20, c: 15}, t=5)
sim.plot("simulation_biltin_function.png",title="Simulation - Built-in Function")





