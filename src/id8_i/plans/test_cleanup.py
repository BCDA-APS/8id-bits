from apsbits.core.instrument_init import oregistry

pv_registers = oregistry["pv_registers"]

print(pv_registers.sample0_pos.get())