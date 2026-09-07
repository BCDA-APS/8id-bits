
from apsbits.core.instrument_init import oregistry
from bluesky import plan_stubs as bps

pv_registers = oregistry["pv_registers"]

def rheo_xpcs_quiescent():

    for ii in range(1):
         yield from wait_for_mcr()

    pv_registers.sample_name.put('P50_quiescent')
    yield from run_measurement_info('measurement_info_quiescent.json') 


def rheo_xpcs_attenuation_test():

    # for ii in range():
    #     yield from wait_for_mcr()

    pv_registers.sample_name.put('G48_attenuation_test')
    yield from run_measurement_info('measurement_info_attenuation_test.json') 





def rheo_xpcs_steadystate():

    for ii in range(2):
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_P52_Load1_RampUpRampDn')
    yield from run_measurement_info('measurement_info_steadystate.json') 
    
    #for ii in range(1):
     #   yield from wait_for_mcr()

    #pv_registers.sample_name.put('SIO2_P_Load2_RampDown')
    #yield from run_measurement_info('measurement_info_steadystate.json') 

def rheo_xpcs_3ITT():

    for ii in range(3):
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_G53_Load1_3ITT_90_1')
    yield from run_measurement_info('measurement_info_3ITT_1.json') 
    
    for ii in range(1):
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_PEG50_Load3_3ITT_005_25s_2')
    yield from run_measurement_info('measurement_info_3ITT_2.json') 


# 3iTT Measurement 20th june
def rheo_xpcs_3ITTrepitition():
    for ii in range(3):
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_PEG48_Load2_3ITT_1_25s')
    yield from run_measurement_info('measurement_info_3ITT_1.json') 
    
    for ii in range(1):
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_PEG48_Load2_3ITT_2_25s')
    yield from run_measurement_info('measurement_info_3ITT_2.json') 

    for ii in range(4):
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_PEG48_Load2_3ITT_1_96s')
    yield from run_measurement_info('measurement_info_3ITT_1.json') 
    
    for ii in range(1):
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_PEG48_Load2_3ITT_2_96s')
    yield from run_measurement_info('measurement_info_3ITT_2.json') 

    for ii in range(4):
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_PEG48_Load2_3ITT_1_173s')
    yield from run_measurement_info('measurement_info_3ITT_1.json') 
    
    for ii in range(1):
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_PEG48_Load2_3ITT_2_173s')
    yield from run_measurement_info('measurement_info_3ITT_2.json') 

    for ii in range(4):
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_PEG48_Load2_3ITT_1_416s')
    yield from run_measurement_info('measurement_info_3ITT_1.json') 
    
    for ii in range(1):
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_PEG48_Load2_3ITT_2_416s')
    yield from run_measurement_info('measurement_info_3ITT_2.json') 

#3ITT Long Run Procedure

def rheo_xpcs_3ITT_longrun():

    for ii in range(1):
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_P50_Load4_3ITT_Ref_LongRun_RepeatLong_')
    yield from run_measurement_info('measurement_info_3ITT_1.json') 
    
    for ii in range(1):
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_P50_Load4_3ITT_10_LongRun_RepeatLong')
    yield from run_measurement_info('measurement_info_3ITT_2.json') 

    for ii in range():
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_P50_Load4_3ITT_30_LongRun_RepeatLong')
    yield from run_measurement_info('measurement_info_3ITT_2.json') 
    
    for ii in range(2):
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_P50_Load4_3ITT_60_LongRun_RepeatLong')
    yield from run_measurement_info('measurement_info_3ITT_2.json') 
    
    for ii in range(2):
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_P50_Load4_3ITT_170_LongRun_Repeat_2')
    yield from run_measurement_info('measurement_info_3ITT_2.json') 


#Shear Stop Procedure

def rheo_xpcs_shearstop():

    for ii in range(5):
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_PEG50_Load1_01_120s_shearstop')
    yield from run_measurement_info('measurement_info_shearstop.json') 

    for ii in range(6):
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_PEG50_Load1_01_230s_shearstop')
    yield from run_measurement_info('measurement_info_shearstop.json') 


def rheo_xpcs_overnight():

#Quiescent Measurement

    for ii in range(2):
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_PEG48_Load1_Quiescent')
    yield from run_measurement_info('measurement_info_quiescent.json') 

#Steady Shear Measurement

    for ii in range(4):
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_PEG48_Load1_SteadyState_RampUp')
    yield from run_measurement_info('measurement_info_steadystate.json') 
    
    for ii in range(2):
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_PEG48_Load1_SteadyState_RampDown')
    yield from run_measurement_info('measurement_info_steadystate.json') 

#3iTT Measurement

    for loop_idx in range(4):
    
        for ii in range(4):
            yield from wait_for_mcr()

        pv_registers.sample_name.put('SIO2_PEG48_Load1_3ITT_1')
        yield from run_measurement_info('measurement_info_3ITT_1.json') 
        
        for ii in range(1):
            yield from wait_for_mcr()

        pv_registers.sample_name.put('SIO2_PEG48_Load1_3ITT_2')
        yield from run_measurement_info('measurement_info_3ITT_2.json') 


#Shear Stop Measurement

    for loop_idx in range(6):
    
        for ii in range(4):
            yield from wait_for_mcr()

        pv_registers.sample_name.put('SIO2_PEG48_Load1_shearstop')
        yield from run_measurement_info('measurement_info_shearstop.json') 

#Full Measurement Set (For P50)

def rheo_xpcs_fulltestset():
    
    #pv_registers.sample_name.put('SIO2_P50_Load4_quiescent_test')
    #yield from run_measurement_info('measurement_info_quiescent.json')
            
    for ii in range(2):
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_P50_Load4_RampUpRampDn_02')
    yield from run_measurement_info('measurement_info_steadystate.json')

    #for ii in range(2):
    #    yield from wait_for_mcr()

    #pv_registers.sample_name.put('SIO2_P50_Load4_3ITT_Ref')
    #yield from run_measurement_info('measurement_info_3ITT_1.json') 
    
    #for ii in range(1):
    #    yield from wait_for_mcr()

    #pv_registers.sample_name.put('SIO2_P50_Load4_3ITT_60s')
    #yield from run_measurement_info('measurement_info_3ITT_2.json') 

    #for ii in range(2):
    #    yield from wait_for_mcr()

    #pv_registers.sample_name.put('SIO2_P50_Load1_3ITT_120s')
    #yield from run_measurement_info('measurement_info_3ITT_2.json') 
    
    #for ii in range(2):
    #    yield from wait_for_mcr()

    #pv_registers.sample_name.put('SIO2_P50_Load1_3ITT_230s')
    #yield from run_measurement_info('measurement_info_3ITT_2.json') 
    
    #for ii in range(2):
    #    yield from wait_for_mcr()

    #pv_registers.sample_name.put('SIO2_P50_Load1_3ITT_415s')
    #yield from run_measurement_info('measurement_info_3ITT_2.json') 

    
    
def rheo_xpcs_3ITT_continuous():
   # for ii in range(1):
       # yield from wait_for_mcr()

    pv_registers.sample_name.put('3ITT_continuous_repeat3')
    yield from run_measurement_info('measurement_info_3ITT_1.json') 

#Catch Up

def rheo_xpcs_3ITT_1():

    for ii in range(1):
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_P50_Load4_3ITT_10_LongRun_Repeat_2')
    yield from run_measurement_info('measurement_info_3ITT_2.json') 

def rheo_xpcs_3ITT_2():

    for ii in range(1):
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_P50_Load4_3ITT_30_LongRun_Repeat_2')
    yield from run_measurement_info('measurement_info_3ITT_2.json') 
    
def rheo_xpcs_3ITT_3():

    for ii in range(1):
        yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_P50_Load4_3ITT_60_LongRun_Repeat_2')
    yield from run_measurement_info('measurement_info_3ITT_2.json') 
    
def rheo_xpcs_3ITT_4():

    #for ii in range(0):
    #    yield from wait_for_mcr()

    pv_registers.sample_name.put('SIO2_P50_Load4_3ITT_170_LongRun_Repeat_2')
    yield from run_measurement_info('measurement_info_3ITT_2.json') 
