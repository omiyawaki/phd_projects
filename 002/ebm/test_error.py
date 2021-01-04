import climlab

state = climlab.column_state(num_lev=3, num_lat=2)
model = climlab.TimeDependentProcess(state=state)
conv = climlab.convection.ConvectiveAdjustment(state=state, adj_lapse_rate='MALR')
model.add_subprocess('convection', conv)
model.integrate_days(1)
