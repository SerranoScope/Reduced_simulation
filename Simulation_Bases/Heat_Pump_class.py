# %%
import pandapipes as pp
from tespy.components import Compressor,SimpleHeatExchanger,CycleCloser, Valve, HeatExchanger,Condenser, Sink, Source
from tespy.connections import Connection
from tespy.networks import Network
import numpy as np
from tespy.tools import UserDefinedEquation


class Bidirectional_W_to_WHeatPump:

    def __init__( self,name,net,refrigerant,HC_ext_id,HC_inj_id):
        self.name = name
        self.refrigerant = refrigerant
        self.HC_ext_id=HC_ext_id
        self.HC_inj_id=HC_inj_id
        self.net=net

    def _build_tespy_Cycles(self):

        
        #Cooling net
        self.nw_Cooling_net = Network()
        self.nw_Cooling_net.units.set_defaults(temperature="degC", pressure="bar",pressure_difference="bar",  enthalpy="J/kg", heat="W", power="W")
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
        
                
        #Heating net 
        self.nw_Heating_net = Network()
        self.nw_Heating_net.units.set_defaults(temperature="degC", pressure="bar",pressure_difference="bar",  enthalpy="J/kg", heat="W", power="W")
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
        
    def solve_cycle(self,mode,Q_consumer,eta_s,T_nework_in,T_cons,dt_DHN):
        #The pressure values for the consumer and district heating side in the heating pump are arbitrary values
        self.mode=mode
        T_cons_in = T_cons[0]
        T_cons_out = T_cons[1]
        T_DH_in = T_nework_in
        def my_ude(ude):
            return ude.conns[0].calc_T_dew()+5-ude.conns[1].calc_T()
        def my_ude_dependents(ude):
            c1, c2 = ude.conns
            return [c1.p,c1.h, c2.p,c2.h]
        def my_ude_2(ude):
            dt_DHN=ude.params['dt']
            return ude.conns[0].calc_T()+dt_DHN-ude.conns[1].calc_T()
        def my_ude_dependents_2(ude):
            c1, c2 = ude.conns
            return [c1.p,c1.h, c2.p,c2.h]
        if self.mode=="COOLING_NET":
            self.Cooling_net_evaporator.set_attr(pr1=1, pr2=1, ttd_l=5)
            self.Cooling_net_condenser.set_attr(pr1=1, pr2=1, ttd_u=5,Q=Q_consumer)
            self.Cooling_net_compressor.set_attr(eta_s=eta_s)
            self.Cooling_net_c2.set_attr(fluid={self.refrigerant: 1})
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
            self.Cooling_ude = UserDefinedEquation(
            'my_ude', my_ude, my_ude_dependents, conns=[self.Cooling_net_c1, self.Cooling_net_c2])
            self.Cooling_ude_2 = UserDefinedEquation(
                            'my_ude_2', my_ude_2, my_ude_dependents_2, conns=[self.Cooling_net_c8, self.Cooling_net_c7],params={'dt':dt_DHN})
            try:
                self.nw_Cooling_net.add_ude(self.Cooling_ude)
                self.nw_Cooling_net.add_ude(self.Cooling_ude_2)
                self.nw_Cooling_net.solve('design')
            except ValueError: 
                self.nw_Cooling_net.del_ude(self.Cooling_ude)
                self.nw_Cooling_net.del_ude(self.Cooling_ude_2)
                self.nw_Cooling_net.add_ude(self.Cooling_ude)
                self.nw_Cooling_net.add_ude(self.Cooling_ude_2)
                self.nw_Cooling_net.solve('design')


            
        elif self.mode=="HEATING_NET":
            self.Heating_net_evaporator.set_attr(pr1=1, pr2=1, ttd_l=5,Q=Q_consumer)
            self.Heating_net_condenser.set_attr(pr1=1, pr2=1, ttd_u=5)
            self.Heating_net_compressor.set_attr(eta_s=eta_s)
            self.Heating_net_c2.set_attr(fluid={self.refrigerant: 1})
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
            self.Heating_ude = UserDefinedEquation(
            'my_ude_3', my_ude, my_ude_dependents, conns=[self.Heating_net_c1, self.Heating_net_c2])
            self.Heating_ude_2 = UserDefinedEquation(
                            'my_ude_4', my_ude_2, my_ude_dependents_2, conns=[self.Heating_net_c7, self.Heating_net_c8],params={'dt':dt_DHN})
            try:
                self.nw_Heating_net.add_ude(self.Heating_ude)
                self.nw_Heating_net.add_ude(self.Heating_ude_2)
                self.nw_Heating_net.solve('design')
            except ValueError: 
                self.nw_Heating_net.del_ude(self.Heating_ude)
                self.nw_Heating_net.del_ude(self.Heating_ude_2)
                self.nw_Heating_net.add_ude(self.Heating_ude)
                self.nw_Heating_net.add_ude(self.Heating_ude_2)
                self.nw_Heating_net.solve('design')
    def interaction_simulation(self,dT_water):
        if self.mode=="COOLING_NET":
            self.net.heat_consumer.at[self.HC_ext_id, "qext_w"] =abs( self.Cooling_net_evaporator.Q.val)
            self.net.heat_consumer.at[self.HC_ext_id,"deltat_k"]=dT_water
            self.net.heat_consumer.at[self.HC_inj_id, "in_service"] =False
            
        elif self.mode=="HEATING_NET":
            self.net.heat_consumer.at[self.HC_inj_id, "qext_w"] = self.Heating_net_condenser.Q.val
            self.net.heat_consumer.at[self.HC_inj_id,"deltat_k"]=dT_water
            self.net.heat_consumer.at[self.HC_ext_id, "in_service"] =False

        else: 
            self.net.heat_consumer.at[self.HC_ext_id, "in_service"] =False
            self.net.heat_consumer.at[self.HC_inj_id, "in_service"] =False

    def get_main_results(self):
        if self.mode=="COOLING_NET":
            results = {
                "Q_consumer": self.Cooling_net_condenser.Q.val,
                "Q_network": self.Cooling_net_evaporator.Q.val,
                "COP": abs(self.Cooling_net_condenser.Q.val)/self.Cooling_net_compressor.P.val 
            }
        elif self.mode=="HEATING_NET":
            results = {
                "Q_consumer": self.Heating_net_evaporator.Q.val,
                "Q_network": self.Heating_net_condenser.Q.val,
                "COP": self.Heating_net_evaporator.Q.val/self.Heating_net_compressor.P.val 
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



