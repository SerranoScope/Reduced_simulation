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
    """
    Almacenamiento termoclina de dos capas con UNA sola bomba reversible
    (circ_pump_mass), sin flow controls.
 
    Convención de signo de mass_flow[t]:
        > 0 : carga     (toma de j_supply, inyecta en j_return a T_cold)
        < 0 : descarga  (toma de j_return, inyecta en j_supply a T_hot)
        = 0 : estático
 
    Protocolo con el orquestador (igual que los demás NetworkActor):
        Def_input_variables(...)            una vez
        Boundary_condition(T_guess, t)      en cada iteración
        pp.pipeflow(net, mode="bidirectional")
        read_T_result(t)                    -> nuevo T_guess hasta converger
        log_hour(t, T_real) / finalize_hour(T_real, t)   tras converger
    """
    net: object
    name: str
    Circ_pump_id: int
    j_supply: int
    j_return: int
    V_tot: float
    T_hot: float              # K
    T_cold: float             # K
    v_hot_fraction: float
    UA_loss: float
    UA_interface: float = 0.0
    T_amb: list = field(default_factory=list)
    mass_flow: list = field(default_factory=list)   # perfil horario (kg/s)
    log: dict = field(default_factory=dict)
    dt_s: float = 3600
    # --- estado operativo (no son parámetros de entrada) ---
    mass_flow_t: float = 0.0
    bypass: bool = False
    direction: str = "static"
    mdot_entering: float = 0.0
    mdot_bypass: float = 0.0
    v_in: float = 0.0
 
    def _post_init(self):
        self.V_hot = self.V_tot * self.v_hot_fraction
        self.V_cold = self.V_tot - self.V_hot
        self.V_MIN = self.V_tot * 0.01
        self.fluid = get_fluid(self.net)
        self.density_average_heat_capacity()
 
    # ------------------------------------------------------------------
    def density_average_heat_capacity(self):
        T_samples = np.linspace(273.15 + 20, 273.15 + 100, 100)
        self.density = float(np.mean([self.fluid.get_density(T) for T in T_samples]))
        self.heat_capacity = float(np.mean([self.fluid.get_heat_capacity(T) for T in T_samples]))
 
    def Def_input_variables(self, input_variables):
        self.mass_flow = input_variables["mass_flow"]   # ya con el signo carga(+)/descarga(-)
        self.T_amb = input_variables["T_amb"]
 
    def _set_pump(self, mdot):
        n, p = self.net, self.Circ_pump_id
        mdot = max(abs(float(mdot)), 1e-4)
        if self.mass_flow_t >= 0:   # carga
            n.circ_pump_mass.at[p, "return_junction"] = self.j_supply
            n.circ_pump_mass.at[p, "flow_junction"] = self.j_return
            n.circ_pump_mass.at[p, "t_flow_k"] = self.T_cold
        else:                       # descarga
            n.circ_pump_mass.at[p, "return_junction"] = self.j_return
            n.circ_pump_mass.at[p, "flow_junction"] = self.j_supply
            n.circ_pump_mass.at[p, "t_flow_k"] = self.T_hot
        n.circ_pump_mass.at[p, "mdot_flow_kg_per_s"] = mdot
        n.circ_pump_mass.at[p, "in_service"] = True
 
    # ------------------------------------------------------------------
    def evaluate_bypass(self, mass_flow, dt_s, T_net, verbose=False):
        bypass, mdot_entering, mdot_bypass, v_in = False, 0.0, 0.0, 0.0
        if mass_flow > 0:                                   # carga
            if T_net < self.T_hot:
                self.direction = "to_Tcold"
                mdot_entering = mass_flow
            else:
                self.direction = "to_Thot"
                v_available = max(0.0, self.V_cold - self.V_MIN)
                Dv_requested = (mass_flow * dt_s) / self.density
                v_in = min(Dv_requested, v_available)
                mdot_entering = (v_in * self.density) / dt_s if dt_s > 0 else 0.0
                mdot_bypass = mass_flow - mdot_entering
                if mdot_bypass > 1e-6 or v_available <= 0:
                    bypass = True
                    if verbose:
                        print(f"  ⚠ [{self.name}] Tanque saturado de calor: "
                              f"{Dv_requested - v_in:.2f} m3 / {mdot_bypass:.2f} kg/s no se pudieron cargar")
        elif mass_flow < 0:                                 # descarga
            self.direction = "discharge"
            abs_mf = abs(mass_flow)
            v_available = max(0.0, self.V_hot - self.V_MIN)
            Dv_requested = (abs_mf * dt_s) / self.density
            v_in = min(Dv_requested, v_available)
            mdot_entering = (v_in * self.density) / dt_s if dt_s > 0 else 0.0
            mdot_bypass = abs_mf - mdot_entering
            if mdot_bypass > 1e-6 or v_available <= 0:
                bypass = True
                if verbose:
                    print(f"  ⚠ [{self.name}] Tanque saturado de frío: "
                          f"{Dv_requested - v_in:.2f} m3 / {mdot_bypass:.2f} kg/s no se pudieron descargar")
        else:
            self.direction = "static"
        self.bypass, self.mdot_entering = bypass, mdot_entering
        self.mdot_bypass, self.v_in = mdot_bypass, v_in
 
    # ------------------------------------------------------------------
    # Interfaz NetworkActor
    # ------------------------------------------------------------------
    def Boundary_condition(self, t_guess, t):
        mf = float(self.mass_flow[t])
        self.mass_flow_t = 0.0 if abs(mf) < 1e-6 else mf    # elimina el -0.0
        self.T_estimate = t_guess
        self.evaluate_bypass(self.mass_flow_t, self.dt_s, t_guess)
        self._set_pump(self.mdot_entering)
 
    def read_T_result(self, t):
        return float(self.net.res_circ_pump_mass.at[self.Circ_pump_id, "t_from_k"])
 
    def Recover_information(self):
        self.Q_net_delivered = float(self.net.res_circ_pump_mass.at[self.Circ_pump_id, "qext_w"])
 
    def Estimate_loss_ambient(self, T_amb):
        UA_hot = self.UA_loss * (self.V_hot / self.V_tot)
        UA_cold = self.UA_loss * (self.V_cold / self.V_tot)
        return UA_hot * (self.T_hot - T_amb), UA_cold * (self.T_cold - T_amb)
 
    def log_hour(self, t, t_real_k):
        self.Recover_information()
        self.log.setdefault(t, {}).update({
            "T_source_k": t_real_k,
            "Q_storage_w": self.Q_net_delivered,
            "direction": self.direction,
            "bypass": self.bypass,
            "mdot_entering": self.mdot_entering,
            "mdot_bypass": self.mdot_bypass,
        })
 
    # ------------------------------------------------------------------
    def V_dis(self, dt_s, T_net, t):
        "Evolución del tanque con Euler implícito"
        rho, cp = self.density, self.heat_capacity
        T_amb = self.T_amb[t]
        b = (dt_s * self.UA_loss) / (rho * cp * self.V_tot)
        if self.mass_flow_t >= 0:
            if self.direction == "to_Thot":
                self.V_cold -= self.v_in
                self.V_hot += self.v_in
                self.v_hot_fraction = self.V_hot / self.V_tot
                a_h = dt_s * self.mdot_entering / (rho * self.V_hot)
                self.T_hot = (self.T_hot + a_h * T_net + b * T_amb) / (1 + a_h + b)
                self.T_cold = (self.T_cold + b * T_amb) / (1 + b)
            elif self.direction == "to_Tcold":
                a_c = dt_s * self.mdot_entering / (rho * self.V_cold)
                self.T_cold = (self.T_cold + a_c * T_net + b * T_amb) / (1 + a_c + b)
                self.T_hot = (self.T_hot + b * T_amb) / (1 + b)
            else:  # estático
                self.T_hot = (self.T_hot + b * T_amb) / (1 + b)
                self.T_cold = (self.T_cold + b * T_amb) / (1 + b)
        else:  # descarga
            self.V_cold += self.v_in
            self.V_hot -= self.v_in
            self.v_hot_fraction = self.V_hot / self.V_tot
            a_c = dt_s * self.mdot_entering / (rho * self.V_cold)
            self.T_cold = (self.T_cold + a_c * T_net + b * T_amb) / (1 + a_c + b)
            self.T_hot = (self.T_hot + b * T_amb) / (1 + b)
 
    def finalize_hour(self, t_real_k, t):
        # Recalcula con la T real convergida (misma lógica que el script suelto)
        self.evaluate_bypass(self.mass_flow_t, self.dt_s, t_real_k, verbose=True)
        Q_loss_hot, Q_loss_cold = self.Estimate_loss_ambient(self.T_amb[t])
        self.V_dis(self.dt_s, t_real_k, t)
        # setdefault: funciona aunque log_hour se llame antes o después
        self.log.setdefault(t, {}).update({
            "Q_loss_storage": Q_loss_hot + Q_loss_cold,
            "T_hot_storage": self.T_hot,
            "T_cold_storage": self.T_cold,
            "V_hot": self.V_hot,
            "V_cold": self.V_cold,
        })

# ---------------------------------------------------------------------------
# REGLAS DE RECORTE: reparten el calor disponible entre los consumidores cuando no alcanza.
# Firma:  regla(demand, budget_w, t, net) -> {nombre_consumidor: W entregados}
#   demand   : {nombre: W pedidos}      budget_w : W totales que se pueden repartir
#   t        : hora                     net      : red de pandapipes (por si la regla la necesita)
# El codigo recorta cada valor al rango [0, demanda], asi que la regla no tiene que vigilarlo.
# ---------------------------------------------------------------------------
def proportional_cut(demand, budget_w, t=None, net=None):
    """Todos reciben el mismo porcentaje de lo que piden."""
    tot = sum(demand.values())
    f = 1.0 if tot <= 0 else min(1.0, max(budget_w, 0.0) / tot)
    return {k: v * f for k, v in demand.items()}

def priority_cut(order):
    """Se sirve completo en el orden dado hasta agotar el calor; el que no cabe recibe lo que queda y el resto 0.
       Los consumidores que no aparezcan en `order` van al final (en su orden original)."""
    def rule(demand, budget_w, t=None, net=None):
        out = {k: 0.0 for k in demand}
        left = max(budget_w, 0.0)
        names = [n for n in order if n in demand] + [n for n in demand if n not in order]
        for n in names:
            give = min(demand[n], left)
            out[n] = give
            left -= give
        return out
    return rule

class MeritOrderDispatcher:
    def __init__(self, plants, order=None, ref=None, tol_w=1.0, order_fn=None, dead_band=1e-3):
        self.plants = plants
        self.tol_w = tol_w
        # Banda muerta (fraccion de Q_max de cada planta): un exceso menor que esto NO obliga a
        # arrancar la planta siguiente. Evita pedirle a una planta 5 W (caudal ~0 -> pipeflow no converge)
        # y evita falsos "infactible" por ruido numerico cuando tol_w es minusculo.
        self.dead_band = dead_band
        self.snapshots = {}      # {hora: estado verdadero vs estado posible}
        self.by_node = {}
        for i, p in enumerate(plants):
            self.by_node[p.node] = i
        self.order_fixed = None
        if order is not None:
            self.order_fixed = [self._idx(x) for x in order]
        self.order_fn = order_fn
        self.ref = None
        if ref is not None:
            self.ref = self._idx(ref)
        self.balance_i = None
        self.slack_i = None
        self.marginal_i = None
        self.stage_k = 0
        self.n_stages = 0
        self.feasible = True
        self.deficit_w = 0.0

    def _idx(self, x):
        if isinstance(x, str):
            return self.by_node[x]
        return int(x)

    def _order_at(self, t):
        if self.order_fn is not None:
            return [self._idx(x) for x in self.order_fn(t)]
        if self.order_fixed is not None:
            return self.order_fixed
        pares = []
        for i, p in enumerate(self.plants):
            pares.append((p.merit_key(t), i))
        pares.sort()
        return [i for clave, i in pares]

    def start_hour(self, t):
        for p in self.plants:
            p.update_Q_max(t)
            p.set_role("off")
        order = self._order_at(t)
        self.stage_order = [i for i in order if self.plants[i].Q_max_t > 0]
        self.n_stages = len(self.stage_order)
        if self.n_stages == 0:
            raise ValueError(f"Hora {t}: ninguna chaufferie disponible")
        ref = self.ref
        if ref is None:
            refs = [i for i in order if self.plants[i].is_reference]
            if refs:
                ref = refs[0]
            else:
                ref = order[0]
        if ref in self.stage_order:
            self.ref_t = ref
        else:
            self.ref_t = self.stage_order[0]
        self.feasible = True
        self.deficit_w = 0.0

    def set_stage(self, k, demand_w):
        P = self.plants
        so = self.stage_order
        active = so[:k + 1]
        sat = so[:k]
        M = so[k]
        if self.ref_t in active:
            slack = self.ref_t
        else:
            slack = M
        for p in P:
            p.set_role("off")
        P[slack].set_role("slack")
        for i in sat:
            if i != slack:
                P[i].set_role("fixed", P[i].Q_max_t)
        self.slack_i = slack
        self.marginal_i = M
        if slack in sat:
            seed = demand_w - sum(P[i].Q_max_t for i in sat)
            P[M].set_role("balance", max(seed, 0.01 * P[M].Q_max_t))
            self.balance_i = M
        else:
            self.balance_i = None

    def set_cut_stage(self):
        """Estado recortado: TODAS las plantas disponibles a su Q_max (rol 'fixed'), salvo la de referencia
        que queda como slack (fija la presion). No hay planta 'balance': la potencia del slack se
        ajusta recortando a los consumidores (ver DHN._solve_cut_state)."""
        P = self.plants
        slack = self.ref_t
        for p in P:
            p.set_role("off")
        P[slack].set_role("slack")
        for i in self.stage_order:
            if i != slack:
                P[i].set_role("fixed", P[i].Q_max_t)
        self.slack_i = slack
        self.marginal_i = slack
        self.balance_i = None

    def update_balance(self):
        if self.balance_i is None:
            return 0.0
        slack = self.plants[self.slack_i]
        b = self.plants[self.balance_i]
        err = slack.Q_delivered() - slack.Q_max_t
        b.Q_target = max(b.Q_target + err, 1.0)
        return abs(err)

    def band(self, i):
        """Exceso admisible de la planta i antes de considerarla fuera de limite."""
        return max(self.tol_w, self.dead_band * self.plants[i].Q_max_t)

    def stage_ok(self):
        self.deficit_by_plant = {}
        for i, p in enumerate(self.plants):
            excess = p.Q_delivered() - p.Q_max_t
            if excess > self.tol_w:
                self.deficit_by_plant[i] = excess
        self.deficit_w = sum(self.deficit_by_plant.values())      # exceso real (0.0 si ninguno)
        ok = all(e <= self.band(i) for i, e in self.deficit_by_plant.items())
        self.feasible = ok
        return ok

    def snapshot(self, t, demand_w=None, converged=True, note=""):
        """Guarda las DOS fotos de la hora:
           - Q_required : lo que la red necesita de cada planta (estado hidraulico convergido, puede pasar de Q_max)
           - Q_feasible : lo que cada planta puede dar de verdad = min(Q_required, Q_max)
           - Q_deficit  : lo que falta = Q_required - Q_feasible
        """
        rows = {}
        for i, p in enumerate(self.plants):
            q_req = p.Q_delivered()
            q_max = p.Q_max_t
            q_fea = min(q_req, q_max)
            rows[p.name] = dict(role=p.role, Q_required=q_req, Q_max=q_max,
                                Q_feasible=q_fea, Q_deficit=max(q_req - q_max, 0.0))
        req = sum(r["Q_required"] for r in rows.values())
        cap = sum(r["Q_max"] for r in rows.values())
        deficit = max(req - cap, 0.0)       # calor que falta en TOTAL (incluye perdidas de red)
        snap = dict(feasible=self.feasible, converged=converged, note=note,
                    demand_w=demand_w, Q_required_total=req, capacity_total=cap,
                    deficit_w=deficit, coverage=(min(1.0, cap / req) if req > 0 else 1.0),
                    overload_by_plant={self.plants[i].name: e for i, e in self.deficit_by_plant.items()},
                    plants=rows)
        self.snapshots[t] = snap
        return snap

@dataclass
class Chaufferie:
    name: str
    node: str
    Fuel_type: str
    Q_nom: float                      # W térmicos
    Nominal_efficiency: float
    Emission_factor: float
    cost_constant: float
    Q_min: float = 0.0                # W (hoy solo informativo, ver Q_available)
    Q_max_profile: list = field(default_factory=list)   # W por hora (disponibilidad externa), opcional
    load: float = 0.0                 # W asignados en la hora (se rellena en allocate)
    Q_cons: float = 0.0               # W de combustible
    priority: object = None                 # None | float | list per hour | callable(t)->float

    def Def_input_variables(self, input_variables):
        self.Q_max_profile = input_variables.get("Q_max_profile", [])

    # ---- ÚNICO sitio donde se decide el Q_max de la hora -------------------
    def Q_max(self, t):
        """Here is the space to put the description of times and ramps or normal ."""
        q = self.Q_nom
        if len(self.Q_max_profile) > 0:
            q = min(q, self.Q_max_profile[t])
        return max(float(q), 0.0)
    def marginal_cost(self):
        """Marginal cost of the chaufferie in €/MWh. It is a constant value for now."""
        return self.cost_constant / self.Nominal_efficiency #€/W

    def merit_key(self, t):
        """Here the merit order is defined. It can be a constant, a callable or a list/array of values."""
        p = self.priority
        if p is None or t is None:
            return self.marginal_cost()
        if callable(p):
            return float(p(t))
        if isinstance(p, (list, tuple, np.ndarray)):
            return float(p[t])
        return float(p)

@dataclass
class Central_production (NetworkActor):
    name:str
    net:object
    node: str
    Circ_mass_id : int
    Circ_pressure_id : int
    Chaufferies: list=field(default_factory=list)
    is_reference: bool=False # Is reference in pressure 
    T_out_design: float=377
    P_out_design: float=3 
    P_lift_design: float=0.5 
    log: dict=field(default_factory=dict)
    #Estado que cambia el dispatchyer en cada etapa
    role: str = "off"                 # off | slack | fixed | balance
    iterates: bool = False 
    Q_target: float = 0.0             # W (solo fixed / balance)
    Q_max_t: float = 0.0              # W (suma de chaufferies, hora actual)
    T_estimate: float = 305.0
    Q_overload: float = 0.0           # W que la red pide por encima de lo que las chaufferies pueden dar (hora actual)
    

    def _post_init(self):
        self.fluid = get_fluid(self.net)
    def Def_input_variables(self,input_variables):
        self.T_out_design=input_variables["T_out_design"]
        # Here if the balue is not described in the input variables, the default value will be used.
        self.P_out_design = input_variables.get("P_out_design", self.P_out_design)
        self.P_lift_design = input_variables.get("P_lift_design", self.P_lift_design)
        per_chaufferies = input_variables.get("per_chaufferies", {})
        for c in self.Chaufferies:
            c.Def_input_variables(per_chaufferies.get(c.name, {}))
    def update_Q_max(self, t):
            self.Q_max_t = sum(c.Q_max(t) for c in self.Chaufferies)
    def set_role(self, role, Q_target=0.0):
        self.role = role
        self.Q_target = float(Q_target)
        # iterates se actualiza aquí, cada vez que cambia el rol
        if role == "fixed" or role == "balance":
            self.iterates = True
        else:
            self.iterates = False       
    def m_dot(self,Q_nominal):
        cp_supply = self.fluid.get_heat_capacity(self.T_out_design)
        cp_return = self.fluid.get_heat_capacity(self.T_estimate)
        delta_t=self.T_out_design- self.T_estimate
        M_dot_estimative=Q_nominal/(((cp_supply+cp_return)/2)*delta_t)
        return M_dot_estimative
    def Boundary_condition(self,T_guess,t):
        self.T_estimate=T_guess
        if self.role=="slack":
            self.net.circ_pump_pressure.at[self.Circ_pressure_id, "in_service"] =True
            self.net.circ_pump_pressure.at[self.Circ_pressure_id, "t_flow_k"] = self.T_out_design
            self.net.circ_pump_pressure.at[self.Circ_pressure_id, "p_flow_bar"] = self.P_out_design
            self.net.circ_pump_pressure.at[self.Circ_pressure_id, "plift_bar"] = self.P_lift_design
            self.net.circ_pump_mass.at[self.Circ_mass_id, "in_service"] =False
        elif self.role in ("fixed","balance"):
            self.net.circ_pump_pressure.at[self.Circ_pressure_id, "in_service"] =False
            self.net.circ_pump_mass.at[self.Circ_mass_id, "in_service"] =True
            M_dot_estimative=max(self.m_dot(self.Q_target), 1e-4)   # suelo: caudal ~0 rompe pipeflow
            self.net.circ_pump_mass.at[self.Circ_mass_id,'mdot_flow_kg_per_s']=M_dot_estimative
            self.net.circ_pump_mass.at[self.Circ_mass_id, "t_flow_k"] = self.T_out_design      
        else:   # "off": apagar AMBAS bombas (si no, una bomba activa de otra etapa/hora sigue inyectando)
            self.net.circ_pump_pressure.at[self.Circ_pressure_id, "in_service"] =False
            self.net.circ_pump_mass.at[self.Circ_mass_id, "in_service"] =False
    def read_T_result(self, t):
        if self.role == "slack":
            return float(self.net.res_circ_pump_pressure.loc[self.Circ_pressure_id, "t_from_k"])
        if self.role in ("fixed", "balance"):
            return float(self.net.res_circ_pump_mass.loc[self.Circ_mass_id, "t_from_k"])
        return self.T_estimate
    def Q_delivered(self):
        if self.role == "slack":
            return float(self.net.res_circ_pump_pressure.at[self.Circ_pressure_id, "qext_w"])
        if self.role in ("fixed", "balance"):
            return float(self.net.res_circ_pump_mass.at[self.Circ_mass_id, "qext_w"])
        return 0.0
    def get_results(self):
        if self.role == "slack":
            tabla, idx = self.net.res_circ_pump_pressure, self.Circ_pressure_id
        elif self.role in ("fixed", "balance"):
            tabla, idx = self.net.res_circ_pump_mass, self.Circ_mass_id
        else:
            return dict(role="off", mdot_kg_s=0.0, Q_net=0.0, Delta_k=0.0,
                        Pin=np.nan, Pout=np.nan, Q_max=self.Q_max_t,
                        Q_feasible=0.0, Q_deficit=0.0)
        return dict(
            role=self.role,
            mdot_kg_s=float(tabla.at[idx, "mdot_from_kg_per_s"]),
            Q_net=float(tabla.at[idx, "qext_w"]),
            Delta_k=float(tabla.at[idx, "deltat_k"]),
            Pin=float(tabla.at[idx, "p_from_bar"]),
            Pout=float(tabla.at[idx, "p_to_bar"]),
            Q_max=self.Q_max_t,
            Q_feasible=min(float(tabla.at[idx, "qext_w"]), self.Q_max_t),
            Q_deficit=max(float(tabla.at[idx, "qext_w"]) - self.Q_max_t, 0.0),
        )
    def log_hour(self, t, t_real_k):
        self.log[t] = {"T_source_k": t_real_k, "role": self.role, **self.get_results()}
    def Recover_information(self,t=None):
        #Distribution between boiler after convergence.
        #Q_net_delivered = lo que la red REALMENTE pide a esta planta (estado verdadero, puede pasar de Q_max).
        #Las chaufferias solo reciben lo que pueden dar (estado posible); el resto queda en Q_overload.
        self.Q_net_delivered = self.Q_delivered()
        remaining = self.Q_net_delivered
        pares = []
        for i, c in enumerate(self.Chaufferies):
            if t is not None:
                clave = c.merit_key(t)
            else:
                clave = c.marginal_cost()
            pares.append((clave, i, c))
        pares.sort()
        ordered = [c for clave, i, c in pares]
        for c in ordered:
            c.load = min(remaining, c.Q_max(t)) if t is not None else remaining
            remaining -= c.load
        self.Q_overload = max(remaining, 0.0) if t is not None else 0.0
    def Consumption_Fuel(self):
        for c in self.Chaufferies:
            c.Q_cons = c.load / c.Nominal_efficiency
        self.Q_cons = sum(c.Q_cons for c in self.Chaufferies)
    def CO2_emissions(self):
        return sum(c.Emission_factor * c.Q_cons for c in self.Chaufferies)
    def Operation_cost(self):
        return sum(c.cost_constant * c.Q_cons for c in self.Chaufferies)