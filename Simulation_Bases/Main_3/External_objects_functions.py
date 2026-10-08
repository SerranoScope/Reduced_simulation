#!/usr/bin/env python
# coding: utf-8

# In[13]:


# In this notebook each one of the exterior objetc's class will be gathered. As some of them need some iterations according to the net conditions, there will be two network actors: 
# - Static actors: Only describe boundary conditions and make calculations after the convergence of the network. 
# - Network actor Objects that iterates according to initial suppositions until convergence. Their main characteristics are that the boundary conditions imposed are dynamic and dependent of exterior factors. 
#   


# In[14]:


from abc import ABC, abstractmethod
import numpy as np
from dataclasses import dataclass,field
import pandapipes as pp
from pandapipes.pf.pipeflow_setup import get_fluid
from tespy.components import Compressor, HeatExchanger, CycleCloser, Valve, Condenser, Sink, Source
from tespy.connections import Connection
from tespy.networks import Network
from tespy.tools import UserDefinedEquation


# ## Static actors: 
# 
# Normally is composed by the slack plant and heat consumers. They only imposed boundary conditions and read results after convergence


# In[ ]:


class StaticActor(ABC): 
    name: str

    @abstractmethod
    def Def_input_variables(self,input_variables):
        raise NotImplementedError
    def Boundary_condition(self,t,):
        raise NotImplementedError
@dataclass
class Pipe_class(StaticActor):
    name:str
    net: object
    pipe_id:int
    T_amb: list=field(default_factory=list)
    def Def_input_variables(self,input_variables):
        self.T_amb=input_variables["T_amb"]
    def Boundary_condition(self,t):
        self.net.pipe.at[self.pipe_id, "text_k"] = self.T_amb[t]

    

@dataclass
class ConsumerActor(StaticActor):
    name:str
    net:object
    heat_consumer_id:int
    T_return_design= float #Kelvin
    Q_building_profile: list=field(default_factory=list)

    def Def_input_variables(self,input_variables):
        self.Q_building_profile=input_variables["Q_building_profile"]
        self.T_return_design=input_variables["T_return_design"]


    def Boundary_condition(self,t):
        Q=self.Q_building_profile[t]
        self.net.heat_consumer.at[self.heat_consumer_id, "qext_w"] = Q
        self.net.heat_consumer.at[self.heat_consumer_id, "treturn_k"] = self.T_return_design

@dataclass
class Slack_Central_production(StaticActor): 
    name:str
    net:object
    Circ_pump_id : int
    Fuel_type: str
    Nominal_efficiency: float 
    Emission_factor: float 
    cost_constant: float 
    T_out_design: float=305 #Kelvin
    P_out_design: float=3 #bar
    P_lift_design: float=0.5 #bar
    type: str='pt' 

    def Def_input_variables(self,input_variables):
        self.T_out_design=input_variables["T_out_design"]
        self.P_out_design=input_variables["P_out_design"]
        self.P_lift_design= input_variables["P_lift_design"]
    def Boundary_condition(self,t):
        self.net.circ_pump_pressure.at[self.Circ_pump_id, "t_flow_k"] = self.T_out_design
        self.net.circ_pump_pressure.at[self.Circ_pump_id, "p_flow_bar"] = self.P_out_design
        self.net.circ_pump_pressure.at[self.Circ_pump_id, "plift_bar"] = self.P_lift_design
    def Recover_information(self):
        self.Q_net_delivered=self.net.res_circ_pump_pressure.at[self.Circ_pump_id, "qext_w"]
    def Consumption_Fuel(self):
        self.Q_cons=self.Q_net_delivered/self.Nominal_efficiency
    def CO2_emissions(self):
        CO2=self.Emission_factor*self.Q_cons
        return CO2
    def Operation_cost(self):
        Euros=self.cost_constant*self.Q_cons
        return Euros    


# ## Network actors: 
# 
# This class will change boundary conditions and everything according to the convergence method defined by the class 


# In[ ]:


class NetworkActor(ABC):
    name:str
    @abstractmethod
    def Def_input_variables(self,input_variables):
        raise NotImplementedError
    def Boundary_condition(self,t_guess,t,input_variable):
        raise NotImplementedError
    def read_T_result(self,t):
        raise NotImplementedError
    def log_hour(self,t,t_real_k):
        raise NotImplementedError 
@dataclass
class Bidirectional_W_to_WHeatPump(NetworkActor):
    name: str
    refrigerant: str
    HP_heat_consumer_id: int
    net: object
    Q_consumer : list=field(default_factory=list)
    mode : list=field(default_factory=list)
    T_cons: list=field(default_factory=list)
    log: dict=field(default_factory=dict)
    dT_water: float= 5
    eta_s : float=0.95
    def _post_init(self):
        self._build_tespy_Cycles()
    def _build_tespy_Cycles(self):
        def my_ude(ude):
            dt_DHN=ude.params['dt']
            return ude.conns[0].calc_T()+dt_DHN-ude.conns[1].calc_T()
        def my_ude_dependents(ude):
            c1, c2 = ude.conns
            return [c1.p,c1.h, c2.p,c2.h]
        self.nw_Cooling_net = Network()
        self.nw_Cooling_net.units.set_defaults(temperature="K", pressure="bar",pressure_difference="bar",  enthalpy="J/kg", heat="W", power="W")
        self.Cooling_net_compressor = Compressor("compresor")
        self.Cooling_net_condenser = Condenser("condensador")
        self.Cooling_net_valve = Valve("valvula_expansion")
        self.Cooling_net_evaporator = HeatExchanger("evaporador")
        self.Cooling_net_cc=CycleCloser('CycleCloser')
        self.Cooling_net_source_consumer=Source("Source_Consumer ")
        self.Cooling_net_sink_consumer=Sink("Sink_consumer")
        self.Cooling_net_source_reseau=Source("Source_Reseau")
        self.Cooling_net_sink_reseau=Sink("Sink_Reseau")
        self.Cooling_net_c0=Connection(self.Cooling_net_valve, 'out1', self.Cooling_net_cc, 'in1', label='0')
        self.Cooling_net_c1 = Connection(self.Cooling_net_cc, 'out1', self.Cooling_net_evaporator, 'in2', label='1')
        self.Cooling_net_c2 = Connection(self.Cooling_net_evaporator, 'out2', self.Cooling_net_compressor, 'in1', label='2')
        self.Cooling_net_c3 = Connection(self.Cooling_net_compressor, 'out1', self.Cooling_net_condenser, 'in1', label='3')
        self.Cooling_net_c4 = Connection(self.Cooling_net_condenser, 'out1', self.Cooling_net_valve, 'in1', label='4')
        self.Cooling_net_c5=Connection(self.Cooling_net_condenser, 'out2',self.Cooling_net_sink_consumer, 'in1', label='5')
        self.Cooling_net_c6=Connection(self.Cooling_net_source_consumer, 'out1',self.Cooling_net_condenser , 'in2', label='6')
        self.Cooling_net_c7=Connection(self.Cooling_net_source_reseau, 'out1',self.Cooling_net_evaporator , 'in1', label='7')
        self.Cooling_net_c8=Connection(self.Cooling_net_evaporator, 'out1',self.Cooling_net_sink_reseau , 'in1', label='8')
        self.nw_Cooling_net.add_conns(self.Cooling_net_c0,  self.Cooling_net_c1,  self.Cooling_net_c2,  self.Cooling_net_c3,  self.Cooling_net_c4 , self.Cooling_net_c5,  self.Cooling_net_c6, self.Cooling_net_c7, self.Cooling_net_c8)
        self.Cooling_ude = UserDefinedEquation(
                    'my_ude', my_ude, my_ude_dependents, conns=[self.Cooling_net_c8, self.Cooling_net_c7],params={'dt': 5})
        self.nw_Cooling_net.add_ude(self.Cooling_ude)

        #Heating net 
        self.nw_Heating_net = Network()
        self.nw_Heating_net.units.set_defaults(temperature="K", pressure="bar",pressure_difference="bar",  enthalpy="J/kg", heat="W", power="W")
        self.Heating_net_compressor = Compressor("compresor")
        self.Heating_net_condenser = Condenser("condensador")
        self.Heating_net_valve = Valve("valvula_expansion")
        self.Heating_net_evaporator = HeatExchanger("evaporador")
        self.Heating_net_cc=CycleCloser('CycleCloser')
        self.Heating_net_source_consumer=Source("Source_Consumer ")
        self.Heating_net_sink_consumer=Sink("Sink_consumer")
        self.Heating_net_source_reseau=Source("Source_Reseau")
        self.Heating_net_sink_reseau=Sink("Sink_Reseau")
        self.Heating_net_c0=Connection(self.Heating_net_valve, 'out1', self.Heating_net_cc, 'in1', label='0')
        self.Heating_net_c1 = Connection(self.Heating_net_cc, 'out1', self.Heating_net_evaporator, 'in2', label='1')
        self.Heating_net_c2 = Connection(self.Heating_net_evaporator, 'out2', self.Heating_net_compressor, 'in1', label='2')
        self.Heating_net_c3 = Connection(self.Heating_net_compressor, 'out1', self.Heating_net_condenser, 'in1', label='3')
        self.Heating_net_c4 = Connection(self.Heating_net_condenser, 'out1', self.Heating_net_valve, 'in1', label='4')
        self.Heating_net_c5 = Connection(self.Heating_net_evaporator, 'out1', self.Heating_net_sink_consumer, 'in1', label='5')
        self.Heating_net_c6 = Connection(self.Heating_net_source_consumer, 'out1', self.Heating_net_evaporator, 'in1', label='6')  
        self.Heating_net_c7 = Connection(self.Heating_net_source_reseau, 'out1', self.Heating_net_condenser, 'in2', label='7')
        self.Heating_net_c8 = Connection(self.Heating_net_condenser, 'out2', self.Heating_net_sink_reseau, 'in1', label='8')
        self.nw_Heating_net.add_conns(self.Heating_net_c0,  self.Heating_net_c1,  self.Heating_net_c2,  self.Heating_net_c3,  self.Heating_net_c4 , self.Heating_net_c5,  self.Heating_net_c6, self.Heating_net_c7, self.Heating_net_c8)
        self.Heating_ude = UserDefinedEquation(
                        'my_ude_4', my_ude, my_ude_dependents, conns=[self.Heating_net_c7, self.Heating_net_c8],params={'dt': 5})
        self.nw_Heating_net.add_ude(self.Heating_ude)

    def solve_cycle(self,mode,t,T_estimate_in,T_cons):
        #The pressure values for the consumer and district heating side in the heating pump are arbitrary values
        T_cons_in = T_cons[0]
        T_cons_out = T_cons[1]
        T_DH_in = T_estimate_in
        if mode=="COOLING_NET":
            self.Cooling_net_evaporator.set_attr(pr1=1, pr2=1, ttd_l=5)
            self.Cooling_net_condenser.set_attr(pr1=1, pr2=1, ttd_u=5,Q=self.Q_consumer[t])
            self.Cooling_net_compressor.set_attr(eta_s=self.eta_s)
            self.Cooling_net_c2.set_attr(fluid={self.refrigerant: 1},td_dew=5)
            # 6. Parámetros del Consumidor 
            self.Cooling_net_c5.set_attr(T=T_cons_out, p=3, fluid={"water": 1})
            self.Cooling_net_c6.set_attr(T=T_cons_in)
            # 7. Parámetros de la Red de Distrito (Recibe calor, se calienta de 50 °C a 60 °C)
            self.Cooling_net_c7.set_attr(T=T_DH_in, p=2.5, fluid={"water": 1})
            import CoolProp.CoolProp as CP
            T_triple = CP.Props1SI("Ttriple", self.refrigerant)        # Triple point temperature (K)
            p_triple = CP.Props1SI("ptriple", self.refrigerant)        # Triple point pressure (Pa)
            T_critical = CP.Props1SI("T_critical", self.refrigerant)  # Critical temperature (K)
            p_critical = CP.Props1SI("p_critical", self.refrigerant)  # Critical pressure (Pa)
            h_min = CP.PropsSI("H", "T", T_triple + 0.1, "Q", 0, self.refrigerant)
            T_max_K =  T_critical*0.9
            p_high = min(p_critical * 0.9, 30e5) 
            h_max = CP.PropsSI("H", "T", T_max_K, "P", p_high, self.refrigerant)
            self.nw_Cooling_net._set_p_range([p_triple, p_high])
            self.nw_Cooling_net._set_h_range([h_min,h_max])
            self.Cooling_ude.params['dt'] = abs(self.dT_water[t])
            self.nw_Cooling_net.solve('design')        
        elif self.mode[t]=="HEATING_NET":
            self.Heating_net_evaporator.set_attr(pr1=1, pr2=1, ttd_l=5,Q=self.Q_consumer[t])
            self.Heating_net_condenser.set_attr(pr1=1, pr2=1, ttd_u=5)
            self.Heating_net_compressor.set_attr(eta_s=self.eta_s)
            self.Heating_net_c2.set_attr(fluid={self.refrigerant: 1},td_dew=5)
            # 6. Parámetros del Consumidor 
            self.Heating_net_c5.set_attr(T=T_cons_out, p=3, fluid={"water": 1})
            self.Heating_net_c6.set_attr(T=T_cons_in)
            # 7. Parámetros de la Red de Distrito (Recibe calor, se calienta de 50 °C a 60 °C)
            self.Heating_net_c7.set_attr(T=T_DH_in, p=2.5, fluid={"water": 1})
            import CoolProp.CoolProp as CP
            T_triple = CP.Props1SI("Ttriple", self.refrigerant)        # Triple point temperature (K)
            p_triple = CP.Props1SI("ptriple", self.refrigerant)        # Triple point pressure (Pa)
            T_critical = CP.Props1SI("T_critical", self.refrigerant)  # Critical temperature (K)
            p_critical = CP.Props1SI("p_critical", self.refrigerant)  # Critical pressure (Pa)
            h_min = CP.PropsSI("H", "T", T_triple + 0.1, "Q", 0, self.refrigerant)
            T_max_K =  T_critical*0.9
            p_high = min(p_critical * 0.9, 30e5) 
            h_max = CP.PropsSI("H", "T", T_max_K, "P", p_high, self.refrigerant)
            self.nw_Heating_net._set_p_range([p_triple, p_high])
            self.nw_Heating_net._set_h_range([h_min,h_max])
            self.Heating_ude.params['dt'] = abs(self.dT_water[t])
            self.nw_Heating_net.solve('design')
    def Def_input_variables(self,input_variables):
            self.mode=input_variables["mode"]
            self.Q_consumer=input_variables["Q_consumer"]
            self.T_cons=input_variables["T_cons"]
            self.dT_water=input_variables["dT_water"]
    def Boundary_condition(self,t_guess,t):
        self.T_estimate=t_guess
        self.solve_cycle(self.mode[t],t,self.T_estimate,self.T_cons[t])
        if self.mode[t]=="COOLING_NET":
            self.net.heat_consumer.at[self.HP_heat_consumer_id, "qext_w"] =abs(self.Cooling_net_evaporator.Q.val)
            self.net.heat_consumer.at[self.HP_heat_consumer_id,"deltat_k"]=self.dT_water[t]

        elif self.mode[t]=="HEATING_NET":
            self.net.heat_consumer.at[self.HP_heat_consumer_id, "qext_w"] = self.Heating_net_condenser.Q.val
            self.net.heat_consumer.at[self.HP_heat_consumer_id,"deltat_k"]=self.dT_water[t]

        else: 
            self.net.heat_consumer.at[self.HP_heat_consumer_id, "in_service"] =False
            self.net.heat_consumer.at[self.HP_heat_consumer_id, "in_service"] =False

    def get_main_results(self):
        if self.mode=="COOLING_NET":
            results = {
                "Q_consumer": abs(self.Cooling_net_condenser.Q.val),
                "Q_network": abs(self.Cooling_net_evaporator.Q.val),
                "W_compressor":self.Cooling_net_compressor.P.val ,
                "COP": abs(self.Cooling_net_condenser.Q.val)/self.Cooling_net_compressor.P.val 
            }
        elif self.mode=="HEATING_NET":
            results = {
                "Q_consumer": abs(self.Heating_net_evaporator.Q.val),
                "Q_network": abs(self.Heating_net_condenser.Q.val),
                "W_compressor":self.Heating_net_compressor.P.val ,
                "COP": abs(self.Heating_net_evaporator.Q.val)/self.Heating_net_compressor.P.val 
            }
        else:
            results = {}
        return results

    def get_results_solver(self):
        if self.mode=="COOLING_NET":
            self.nw_Cooling_net.print_results()
        elif self.mode=="HEATING_NET":
            self.nw_Heating_net.print_results()
        else:
            return None   
    def read_T_result(self,t):
        return self.net.res_heat_consumer.at[self.HP_heat_consumer_id, "t_from_k"]
    def log_hour(self, t, t_real_k):
        self.log[t] = {"T_source_k": t_real_k, **self.get_main_results()}



# In[ ]:


@dataclass
class ThermoclineTwoLayer(NetworkActor):
    "fc= flow control and cp= circ pump mass"
    net: object 
    name: str
    fc_charge : int 
    fc_bypass_charge : int 
    cp_charge_storage : int 
    fc_charge_storage : int 
    fc_discharge : int 
    fc_bypass_discharge : int 
    cp_discharge_storage : int 
    fc_discharge_storage : int  
    V_tot: float
    T_hot: float       # Kelvin
    T_cold: float     # Kelvin
    v_hot_fraction: float
    UA_loss: float
    UA_interface: float
    T_amb: list=field(default_factory=list)
    log: dict=field(default_factory=dict)
    mass_flow : list=field(default_factory=list)
    bypass: bool=False
    direction : str | bool=False
    mdot_entering: float=0.0
    mdot_bypass: float= 0.0
    v_in =float =0.0
    dt_s: float=3600

    def _post_init(self):
        self.V_hot= self.V_tot * self.v_hot_fraction
        self.V_cold=self.V_tot - self.V_hot
        self.V_MIN=self.V_tot * 0.01
        self.fluid = get_fluid(self.net)
        self.density_average_heat_capacity()

    def density_average_heat_capacity(self):
        T_min_K, T_max_K = 273.15 + 20, 273.15 + 100
        T_samples = np.linspace(T_min_K, T_max_K, 100)
        self.density = np.mean([self.fluid.get_density(T) for T in T_samples])
        self.heat_capacity = np.mean([self.fluid.get_heat_capacity(T) for T in T_samples]) 

    def evaluate_bypass(self, mass_flow,dt_s,T_net): 
        bypass,mdot_entering,mdot_bypass,v_in=False,0.0,0.0,0.0
        if mass_flow > 0:  
            if T_net < self.T_hot:
                self.direction = "to_Tcold"
                mdot_entering=mass_flow
            else:
                self.direction = "to_Thot"
                v_available = max(0.0, self.V_cold - self.V_MIN)
                Dv_requested = (mass_flow * dt_s) / self.density
                v_in = min(Dv_requested, v_available)
                mdot_entering = (v_in * self.density) / dt_s if dt_s > 0 else 0.0
                mdot_bypass = mass_flow - mdot_entering
                if mdot_bypass > 1e-6 or v_available <= 0:
                    bypass = True
                    print(f"  ⚠ [{self.name}] Tanque saturado de calor: "f"{Dv_requested - v_in:.2f} m3 no se pudieron cargar")
                    print(f"  ⚠ [{self.name}] Tanque saturado de calor: "f"{mdot_bypass:.2f} kg/s no se pudieron cargar")
        elif mass_flow < 0: # Discharge
            self.direction = "discharge"
            abs_mass_flow = abs(mass_flow)
            v_available = max(0.0, self.V_hot - self.V_MIN)
            Dv_requested = (abs_mass_flow * dt_s) / self.density
            v_in = min(Dv_requested, v_available)
            mdot_entering = (v_in * self.density) / dt_s if dt_s > 0 else 0.0
            mdot_bypass = abs_mass_flow - mdot_entering
            if mdot_bypass > 1e-6 or v_available <= 0:
                bypass = True
                print(f"  ⚠ [{self.name}] Tanque saturado de frio: "f"{Dv_requested - v_in:.2f} m3 no se pudieron descargar")
                print(f"  ⚠ [{self.name}] Tanque saturado de frio: "f"{mdot_bypass:.2f} kg/s no se pudieron descargar")
        elif mass_flow == 0:
            self.direction = "static"
        self.bypass,self.mdot_entering,self.mdot_bypass,self.v_in= bypass,mdot_entering,mdot_bypass,v_in

    def Def_input_variables(self,input_variables):
        self.mass_flow=input_variables["mass_flow"]
        self.T_amb=input_variables["T_amb"]
    def Boundary_condition(self,t_guess,t):
        self.T_estimate=t_guess
        mass_flow=self.mass_flow[t]
        self.evaluate_bypass(mass_flow,dt_s=3600,T_net=t_guess)
        if mass_flow >= 0:
            self.net.flow_control.at[self.fc_discharge_storage, "in_service"] = False
            self.net.flow_control.at[self.fc_discharge, "in_service"] = False
            self.net.flow_control.at[self.fc_bypass_discharge, "in_service"] = False
            self.net.flow_control.at[self.fc_charge_storage, "in_service"] = True
            self.net.flow_control.at[self.fc_charge, "in_service"] = True
            self.net.flow_control.at[self.fc_bypass_charge, "in_service"] = True
            self.net.circ_pump_mass.at[self.cp_discharge_storage, "mdot_flow_kg_per_s"] = 0
            self.net.flow_control.at[self.fc_bypass_charge, "controlled_mdot_kg_per_s"] = self.mdot_bypass
            self.net.circ_pump_mass.at[self.cp_charge_storage, "mdot_flow_kg_per_s"] = self.mdot_entering
            self.net.flow_control.at[self.fc_charge_storage, "controlled_mdot_kg_per_s"] = self.mdot_entering
            self.net.flow_control.at[self.fc_charge, "controlled_mdot_kg_per_s"] = mass_flow
            self.net.circ_pump_mass.at[self.cp_charge_storage, "t_flow_k"] = self.T_cold
        else:
            self.net.flow_control.at[self.fc_charge_storage, "in_service"] = False
            self.net.flow_control.at[self.fc_charge, "in_service"] = False
            self.net.flow_control.at[self.fc_bypass_charge, "in_service"] = False
            self.net.flow_control.at[self.fc_discharge_storage, "in_service"] = True
            self.net.flow_control.at[self.fc_discharge, "in_service"] = True
            self.net.flow_control.at[self.fc_bypass_discharge, "in_service"] = True
            self.net.circ_pump_mass.at[self.cp_charge_storage, "mdot_flow_kg_per_s"] = 0
            self.net.flow_control.at[self.fc_bypass_discharge, "controlled_mdot_kg_per_s"] = self.mdot_bypass
            self.net.circ_pump_mass.at[self.cp_discharge_storage, "mdot_flow_kg_per_s"] = self.mdot_entering
            self.net.flow_control.at[self.fc_discharge_storage, "controlled_mdot_kg_per_s"] = self.mdot_entering
            self.net.flow_control.at[self.fc_discharge, "controlled_mdot_kg_per_s"] = abs(mass_flow)
            self.net.circ_pump_mass.at[self.cp_discharge_storage, "t_flow_k"] = self.T_hot

    def read_T_result(self,t):
        fc = self.fc_charge if self.mass_flow[t] >= 0 else self.fc_discharge
        return self.net.res_flow_control.at[fc, "t_from_k"]
    def log_hour(self, t, t_real_k):
        self.log[t] = {"T_source_k": t_real_k}
    def V_dis(self, dt_s, T_net,t):
        "The resolution of the equation will be done with euler implicit"
        mass_flow=self.mass_flow[t]
        rho, cp = self.density, self.heat_capacity
        b=(dt_s*self.UA_loss)/(rho*cp*self.V_tot)
        #Euler implicit volume
        if mass_flow >= 0:
            if self.direction=='to_Thot':
                self.V_cold -= self.v_in
                self.V_hot += self.v_in
                self.v_hot_fraction = self.V_hot / self.V_tot
                self.UA_hot_layer = self.UA_loss * (self.V_hot / self.V_tot)
                self.UA_cold_layer = self.UA_loss * (self.V_cold / self.V_tot)
                a_h = dt_s * self.mdot_entering / (rho * self.V_hot)
                a_c = dt_s * self.mdot_entering / (rho * self.V_cold)
                self.T_hot  = (self.T_hot + a_h * T_net + b * self.T_amb[t]) / (1 + a_h + b)
                self.T_cold = (self.T_cold + b * self.T_amb[t]) / (1 + b)
            elif self.direction=='to_Tcold':
                self.UA_hot_layer = self.UA_loss * (self.V_hot / self.V_tot)
                self.UA_cold_layer = self.UA_loss * (self.V_cold / self.V_tot)
                a_h = dt_s * self.mdot_entering / (rho * self.V_hot)
                a_c = dt_s * self.mdot_entering / (rho * self.V_cold)
                self.T_cold  = (self.T_cold + a_c * T_net + b * self.T_amb[t]) / (1 + a_c + b)
                self.T_hot = (self.T_hot + b * self.T_amb[t]) / (1 + b)
            else:
                self.UA_hot_layer = self.UA_loss * (self.V_hot / self.V_tot)
                self.UA_cold_layer = self.UA_loss * (self.V_cold / self.V_tot)
                self.T_hot = (self.T_hot + b * self.T_amb[t]) / (1 + b)
                self.T_cold = (self.T_cold + b * self.T_amb[t]) / (1 + b)
        else:
            self.V_cold += self.v_in
            self.V_hot -= self.v_in
            self.UA_hot_layer = self.UA_loss * (self.V_hot / self.V_tot)
            self.UA_cold_layer = self.UA_loss * (self.V_cold / self.V_tot)
            a_h = dt_s * self.mdot_entering / (rho * self.V_hot)
            a_c = dt_s * self.mdot_entering / (rho * self.V_cold)
            self.T_cold  = (self.T_cold + a_c * T_net + b * self.T_amb[t]) / (1 + a_c + b)
            self.T_hot = (self.T_hot + b * self.T_amb[t]) / (1 + b)

    def finalize_hour(self,t_real_k,t):
        mass_flow=self.mass_flow[t]
        self.evaluate_bypass(mass_flow,dt_s=self.dt_s,T_net=t_real_k)
        self.V_dis(dt_s=self.dt_s,T_net=t_real_k,t=t)


# In[ ):

@dataclass
class Central_production (NetworkActor):
    name:str
    net:object
    Circ_mass_id : int
    flow_control_id : int
    Fuel_type: str
    Nominal_efficiency: float 
    Emission_factor: float 
    cost_constant: float 
    Q_nominal: list=field(default_factory=list) #W
    T_out_design: float=305
    log: dict=field(default_factory=dict)
    def _post_init(self):
        self.fluid = get_fluid(self.net)
    def Def_input_variables(self,input_variables):
        self.T_out_design=input_variables['T_out_design']
        self.Q_nominal=input_variables["Q_nominal"]
    def log_hour(self, t, t_real_k):
            self.log[t] = {"T_source_k": t_real_k}
    def Boundary_condition(self,T_guess,t):
        self.net.circ_pump_mass.at[self.Circ_mass_id, "t_flow_k"] = self.T_out_design
        self.T_estimate=T_guess
        delta_t=self.T_out_design- self.T_estimate
        cp_supply = self.fluid.get_heat_capacity(self.T_out_design)
        cp_return = self.fluid.get_heat_capacity(self.T_estimate)
        M_dot_estimative=self.Q_nominal[t]/(((cp_supply+cp_return)/2)*delta_t)
        self.net.flow_control.at[self.flow_control_id,'controlled_mdot_kg_per_s']= M_dot_estimative
        self.net.circ_pump_mass.at[self.Circ_mass_id,'mdot_flow_kg_per_s']=M_dot_estimative
    def read_T_result(self, t):
        T_result=self.net.res_circ_pump_mass.at[self.Circ_mass_id, "t_from_k"]
        return T_result
    def Recover_information(self):
        self.Q_net_delivered=self.net.res_circ_pump_mass.at[self.Circ_mass_id, "qext_w"]
    def Consumption_Fuel(self):
        self.Q_cons=self.Q_net_delivered/self.Nominal_efficiency
    def CO2_emissions(self):
        CO2=self.Emission_factor*self.Q_cons
        return CO2
    def Operation_cost(self):
        Euros=self.cost_constant*self.Q_cons
        return Euros    


