import sys, getopt
import json
from math import sqrt,log

# global constants
R = 287.06 # Specific gas constant
Rho = 1.225 # Air density at sea level, kg/m^3
LHV = 4.31E7 # Kerosene lower heating value, J/kg
CRITICAL_NPR = 1.8929 # Chocked Nozzle Mach No

CP_VALUES = {
     200 : 1001, 550 : 1040, 1000 : 1142, 1700 : 1229,
     250 : 1003, 600 : 1051, 1100 : 1155, 1800 : 1237,
     300 : 1005, 650 : 1063, 1200 : 1173, 1900 : 1244,
     350 : 1008, 700 : 1075, 1300 : 1190, 2000 : 1250,
     400 : 1013, 750 : 1087, 1400 : 1204, 2100 : 1255,
     450 : 1020, 800 : 1099, 1500 : 1216, 2200 : 1260,
     500 : 1029, 900 : 1121, 1600 : 1221, 2300 : 1265
}

GAMMA_VALUES = {
     200 : 1.402, 550 : 1.381, 1000 : 1.366, 1700 : 1.305,
     250 : 1.401, 600 : 1.376, 1100 : 1.331, 1800 : 1.302,
     300 : 1.400, 650 : 1.370, 1200 : 1.324, 1900 : 1.300,
     350 : 1.398, 700 : 1.364, 1300 : 1.318, 2000 : 1.298,
     400 : 1.395, 750 : 1.359, 1400 : 1.313, 2100 : 1.295,
     450 : 1.391, 800 : 1.354, 1500 : 1.309, 2200 : 1.294,
     500 : 1.387, 900 : 1.344, 1600 : 1.308, 2300 : 1.293
}

# computed values
engine_data = {
    "stations": [],
    "thrusts": [],
    "energies": [],
    "works": {},
    "p_ambient": 0.0,
    "t_ambient": 0.0,
    "mass_flow": 0.0,
    "total_flow": 0.0,
    "fuel_flow": 0.0,
    "net_work": 0.0,
    "net_thrust": 0.0,
    "net_heat": 0.0,
    "engine_class": "TURBOJET",
    "comments":""
}

def interpolate(t,TABLE):
    if t in TABLE:
        return TABLE[t]
    else:
        lower=int(round(t/100,0)*100)
        upper=int(round(t/100,0)+1)*100
        x = (TABLE[lower]+(t-lower)*(TABLE[upper]-TABLE[lower])/(upper-lower))
        return x

def cp_calc(t):
    return interpolate(t,CP_VALUES)

def gamma_calc(t):
    return interpolate(t,GAMMA_VALUES)

def heat_exchanger(t_in, p_in, t_exit, p_loss, n_is):
    t_out = t_in - n_is*(t_in-t_exit)
    p_out = p_in * (1-p_loss)
    return t_out, p_out

def converging_nozzle_thrust(t,p,p_amb,mflow):
    npr = p/p_amb
    y = gamma_calc(t)
    if npr>CRITICAL_NPR:
        npr = CRITICAL_NPR # Mach No = 1.0
    t_out = t/(npr**((y-1)/y))
    p_out = p/npr
    density = p_out/(R*t_out)
    if npr<CRITICAL_NPR:
        mach_no = sqrt((t/t_out-1)*2/(y-1))
    else:
        mach_no = 1.0
    u = sqrt(y*t_out*R)
    area = mflow/(density*u*mach_no)
    f = mflow*u*mach_no+(p_out-p_amb)*area
    return {'AREA_NOZZLE':area,'EXIT_P':p_out,'MACH_NO':mach_no, 'V_EXIT':u,'THRUST':f }

def compressor_stage(T_in, P_in, PR, N_eff):
    y = gamma_calc(T_in)
    T_out_is = T_in*PR**((y-1)/y)
    T = T_in + (T_out_is-T_in)/N_eff
    P = P_in*PR
    return (T,P)

def compute_intake(item,engine_data):
    stations = engine_data["stations"]
    engine_data["t_ambient"] = item['T_INTAKE']
    engine_data["p_ambient"] = item['P_INTAKE']
    engine_data["mass_flow"] = item['MFLOW_TOTAL']
    engine_data["total_flow"] = item['MFLOW_TOTAL']
    stations.extend([
        {
            'ID':'ENTRY',
            'T':item['T_INTAKE'],
            'P':item['P_INTAKE']
        },
        {
            'ID':item['ID'],
            'T':item['T_INTAKE'],
            'P':item['P_INTAKE']*(1-item['P_LOSS'])
        }
    ])
    return engine_data

def compute_fan(item,engine_data):
    engine_data["engine_class"] = "TURBOFAN"
    stations = engine_data["stations"]
    p_ambient = engine_data["p_ambient"]
    T_in = stations[-1]['T']
    P_in = stations[-1]['P']
    T_out, P_out = compressor_stage(T_in,P_in,item['FPR'],item['EFF'])
    engine_data["works"][item['ID']]=cp_calc(T_in)*(T_out - T_in)*engine_data["mass_flow"]
    #split the air flow after fan exit
    engine_data["mass_flow"] = engine_data["mass_flow"]/(item['BPR']+1)
    fan_flow = engine_data["mass_flow"]*item['BPR']
    engine_data["comments"]="%s%s" % (engine_data["comments"],"FAN:\n")
    fan_thrust = converging_nozzle_thrust(T_out,P_out,p_ambient,fan_flow)
    area = fan_thrust['AREA_NOZZLE']
    mach_no = fan_thrust['MACH_NO']
    engine_data["comments"]="%s%s" % (engine_data["comments"],("\tNozzle area is %3.2f m^2\n" % area))
    engine_data["comments"]="%s%s" % (engine_data["comments"],("\tExit Mach number is %3.2f\n" % mach_no))
    engine_data["thrusts"].append(fan_thrust['THRUST'])
    engine_data["energies"].append((fan_flow*fan_thrust['V_EXIT']**2)/2)
    stations.append({'ID':item['ID'],'T':T_out,'P':P_out})
    return engine_data

def compute_compressor(item,engine_data):
    stations=engine_data["stations"]
    T_in = stations[-1]['T']
    P_in = stations[-1]['P']
    T_out, P_out = compressor_stage(T_in,P_in,item['PR'],item['EFF'])
    engine_data["works"][item['ID']]=cp_calc(T_in)*(T_out - T_in)*engine_data["mass_flow"]
    stations.append({'ID':item['ID'],'T':T_out,'P':P_out})
    return engine_data

def compute_intercooler(item,engine_data):
    stations = engine_data["stations"]
    T_in = stations[-1]['T']
    P_in = stations[-1]['P']
    T_out, P_out = heat_exchanger(T_in, P_in, item['T_COLD'], item['P_LOSS'], item['EFF'])
    stations.append({'ID':item['ID'],'T':T_out,'P':P_out})
    return engine_data

def compute_recuperator_ploss(item,engine_data):
    stations = engine_data["stations"]
    T_in = stations[-1]['T']
    P_in = stations[-1]['P']
    p_loss =  item['P_LOSS']
    T_out = T_in
    P_out = P_in*(1-p_loss)
    stations.append({'ID':item['ID'],'T':T_out,'P':P_out})
    return engine_data

def compute_recuperator(item,engine_data):
    stations = engine_data["stations"]
    n_eff = item['EFF']
    p_loss = item['P_LOSS']
    index = -1
    for x in stations:
        if x['ID'] == item['INPUT_SOURCE']:
            index = stations.index(x)
        if x['ID'] == item['HEAT_SOURCE']:
            T_source = x['T']
    T_in = stations[index]['T']
    P_in = stations[index]['P']
    T_out, P_out = heat_exchanger(T_in, P_in, T_source, p_loss, n_eff)
    index = index + 1
    stations = stations[:index]
    stations.append({'ID':item['ID'],'T':T_out,'P':P_out})
    engine_data["stations"] = stations
    return engine_data

def compute_combustor(item,engine_data):
    stations = engine_data["stations"]
    P_in = stations[-1]['P']
    T_in = stations[-1]['T']
    T_out = item['TET']
    Q_in = cp_calc(T_out)*(T_out - T_in)*engine_data["mass_flow"]
    fuel_flow = Q_in/LHV
    engine_data["fuel_flow"] = fuel_flow
    engine_data["mass_flow"] = engine_data["mass_flow"] + fuel_flow
    engine_data["total_flow"] = engine_data["total_flow"] + fuel_flow
    engine_data["net_heat"] = Q_in
    stations.append({'ID':item['ID'],'T':T_out,'P':P_in*(1-item['P_LOSS'])})
    return engine_data

def compute_turbine(item,engine_data):
    stations = engine_data["stations"]
    T_in = stations[-1]['T']
    P_in = stations[-1]['P']
    cp_turb = cp_calc(T_in)
    y = gamma_calc(T_in)
    # match the turbine to corresponding compressor
    W_t = 0
    for x in item['W_MATCH']:
        W_t = W_t + engine_data["works"][x]
    T_out = T_in - W_t/(cp_turb*engine_data["mass_flow"])
    T_out_is = T_in-(T_in-T_out)/item['EFF']
    P_out = P_in / ((T_in/T_out_is)**(y/(y-1)))
    stations.append({'ID':item['ID'],'T':T_out,'P':P_out})
    return engine_data

def compute_power_turbine(item,engine_data):
    stations = engine_data["stations"]
    engine_data["engine_class"] = "TURBOSHAFT"
    T_in = stations[-1]['T']
    P_in = stations[-1]['P']
    cp_turb = cp_calc(T_in)
    M_exit = item['M_EXIT']
    y = gamma_calc(engine_data["t_ambient"])
    P_out = engine_data["p_ambient"]*(1+(y-1)/2*M_exit**2)**(y/(y-1));
    T_out_is = T_in*(P_in/P_out)**((1-y)/y);
    T_out = T_in - (T_in-T_out_is)*item['EFF'];
    # subtract the residual work from net work
    W_net = (T_in-T_out)*cp_turb*engine_data["mass_flow"];
    for x in item['W_MATCH']:
        W_net = W_net - engine_data["works"][x]
    stations.append({'ID':item['ID'],'T':T_out,'P':P_out})
    engine_data["net_work"] = W_net
    return engine_data

def compute_nozzle(item,engine_data):
    stations = engine_data["stations"]
    T_in = stations[-1]['T']
    P_in = stations[-1]['P']
    T_out = T_in
    P_out = P_in*(1-item['P_LOSS'])
    engine_data["comments"]="%s%s\n" % (engine_data["comments"],'CORE:')
    core_thrust = converging_nozzle_thrust(T_out,P_out,engine_data["p_ambient"],engine_data["mass_flow"])
    area = core_thrust['AREA_NOZZLE']
    mach_no = core_thrust['MACH_NO']
    engine_data["comments"]="%s%s" % (engine_data["comments"],("\tNozzle area is %3.2f m^2\n" % area))
    engine_data["comments"]="%s%s" % (engine_data["comments"],("\tExit Mach number is %3.2f" % mach_no))
    engine_data["energies"].append((engine_data["mass_flow"]*core_thrust['V_EXIT']**2)/2)
    engine_data["thrusts"].append(core_thrust['THRUST'])
    return engine_data

def produce_report(engine_data):
    print(engine_data["comments"])
    print("STATIONS:")
    for item in engine_data["stations"]:
        print("\t%s" % item)
    print ('SUMMARY:')
    #for item in stations:
    #    print ('\t %s' % item)
    print ("\tEngine class: %s" % engine_data["engine_class"])
    print ("\tTotal mass flow: %3.2f kg/s" % engine_data["total_flow"])
    print ("\tFuel flow: %3.2f kg/s" %(engine_data["fuel_flow"]))
    W_net = engine_data["net_work"]
    Q_in = engine_data["net_heat"]
    energies = engine_data["energies"]
    thrusts = engine_data["thrusts"]
    fuel_flow = engine_data["fuel_flow"]
    if engine_data["engine_class"] == 'TURBOSHAFT':
        print ("\tThermal efficiency: %3.2f%%" % (W_net*100/Q_in))
        print ("\tAvailable shaft power: %3.2f MW" % (W_net/1000000))
        print ("\tSFC is %3.3f kg/kW/hr" % (3600*fuel_flow/(W_net/1000)))
    else:
        print ("\tThermal efficiency: %3.2f%%" % (sum(energies)*100/Q_in))
        print ("\tNet thrust: %3.2f kN" % (sum(thrusts)/1000))
        print ("\tSFC is %3.3f kg/kN-s" % (engine_data["fuel_flow"]/(sum(thrusts)/1000)))

def compute_engine(component_data, engine_data):
    second_pass = False
    for item in component_data:
        if item['TYPE'] == 'INTAKE':
            engine_data=compute_intake(item,engine_data)
        if item['TYPE'] == 'FAN':
            engine_data = compute_fan(item,engine_data)
        if item['TYPE'] == 'COMPRESSOR':
            engine_data = compute_compressor(item,engine_data)
        if item['TYPE'] == 'INTERCOOLER':
            engine_data = compute_intercooler(item,engine_data)
        if item['TYPE'] == 'RECUPERATOR':
            # first pass - pressure loss only
            engine_data = compute_recuperator_ploss(item,engine_data)
            second_pass = True
        if item['TYPE'] == 'COMBUSTOR':
            engine_data = compute_combustor(item,engine_data)
        if item['TYPE'] == 'TURBINE':
            engine_data = compute_turbine(item,engine_data)
        if item['TYPE'] == 'POWER_TURBINE':
            engine_data = compute_power_turbine(item,engine_data)
        if item['TYPE'] == 'NOZZLE':
            engine_data = compute_nozzle(item,engine_data)
    # second pass for recuperator calculation
    if second_pass:
        for item in component_data:
            if item['TYPE'] == 'RECUPERATOR':
                engine_data = compute_recuperator(item,engine_data)
            if item['TYPE'] == 'COMBUSTOR':
                engine_data = compute_combustor(item,engine_data)
            if item['TYPE'] == 'POWER_TURBINE':
                engine_data = compute_power_turbine(item,engine_data)
    # print report
    produce_report(engine_data)

def load_json_file(filename):
    f = open(filename, "r")
    str_value  = f.read()
    json_object = json.loads(str_value)
    return json_object

def main(argv):
    inputfile =""
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        print ('turbocalc.py -i <inputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('turbocalc.py -i <inputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
            component_list = load_json_file(inputfile)
            compute_engine(component_list,engine_data)

main(sys.argv[1:])
