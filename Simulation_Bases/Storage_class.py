# %%
import pandapipes as pp
from pandapipes.pf.pipeflow_setup import get_fluid
import pandas as pd
import numpy as np
class ThermoclineTwoLayer():
    def __init__(self, net, name, Flow_control_charge_id,Flow_control_bypass_charge_id,Circ_pump_charge_storage_id,Flow_control_charge_storage_id,Flow_control_discharge_id,
                 Flow_control_bypass_discharge_id,Circ_pump_discharge_storage_id,Flow_control_discharge_storage_id,
                 volume_m3, t_hot_init, t_cold_init, v_hot_fraction_init, UA_loss, UA_interface=0.0):
        self.net = net
        self.fluid = get_fluid(net)
        self.name = name
        self.Flow_control_charge_id = Flow_control_charge_id
        self.Flow_control_bypass_charge_id = Flow_control_bypass_charge_id
        self.Circ_pump_charge_storage_id= Circ_pump_charge_storage_id
        self.Flow_control_charge_storage_id = Flow_control_charge_storage_id
        self.Flow_control_discharge_id = Flow_control_discharge_id
        self.Flow_control_bypass_discharge_id = Flow_control_bypass_discharge_id
        self.Circ_pump_discharge_storage_id = Circ_pump_discharge_storage_id
        self.Flow_control_discharge_storage_id = Flow_control_discharge_storage_id
        self.V_tot = volume_m3
        self.T_hot = t_hot_init
        self.T_cold = t_cold_init
        self.v_hot_fraction = v_hot_fraction_init
        self.UA_loss = UA_loss
        self.UA_interface = UA_interface
        self.V_hot = volume_m3 * self.v_hot_fraction
        self.V_cold = volume_m3 - self.V_hot
        self.V_MIN = volume_m3 * 0.01
        # Inicialización de variables operativas
        self.bypass = False
        self.direction = None
        self.mass_flow = 0.0
        self.mdot_entering = 0.0
        self.mdot_bypass = 0.0
        self.v_in=0
        self.density_average_heat_capacity()

    

    def density_average_heat_capacity(self):
        T_min_K = 273.15 + 20   
        T_max_K = 273.15 + 100  

        # Muestreo fino del rango (cuantos más puntos, más precisa la integral)
        T_samples = np.linspace(T_min_K, T_max_K, 100)

        rho_samples = np.array([self.fluid.get_density(T) for T in T_samples])
        cp_samples = np.array([self.fluid.get_heat_capacity(T) for T in T_samples])

        rho_avg = np.mean(rho_samples)
        cp_avg = np.mean(cp_samples)
        self.density = rho_avg
        self.heat_capacity = cp_avg

    def evaluate_T_network_in(self):
        #Charging 
        if self.mass_flow >= 0:
            self.net.flow_control.at[self.Flow_control_discharge_storage_id, "in_service"] = False  
            self.net.flow_control.at[self.Flow_control_discharge_id, "in_service"] = False
            self.net.flow_control.at[self.Flow_control_bypass_discharge_id, "in_service"] = False
            self.net.flow_control.at[self.Flow_control_bypass_charge_id, "in_service"] = False  
            self.net.flow_control.at[self.Flow_control_charge_storage_id, "in_service"] = True  
            self.net.flow_control.at[self.Flow_control_charge_id, "in_service"] = True
            self.net.circ_pump_mass.at[self.Circ_pump_discharge_storage_id, "mdot_flow_kg_per_s"] = 0
            self.net.circ_pump_mass.at[self.Circ_pump_charge_storage_id, "mdot_flow_kg_per_s"] = self.mass_flow
            self.net.flow_control.at[self.Flow_control_charge_storage_id, "controlled_mdot_kg_per_s"] =self.mass_flow
            self.net.flow_control.at[self.Flow_control_charge_id, "controlled_mdot_kg_per_s"] =self.mass_flow
            self.net.circ_pump_mass.at[self.Circ_pump_charge_storage_id, "t_flow_k"] = self.T_cold
            pp.pipeflow(self.net, mode="bidirectional")
            Tnet=self.net.res_flow_control.at[self.Flow_control_charge_id, "t_from_k"]
        else: #Discharging
            self.net.flow_control.at[self.Flow_control_charge_storage_id, "in_service"] = False  
            self.net.flow_control.at[self.Flow_control_charge_id, "in_service"] = False
            self.net.flow_control.at[self.Flow_control_bypass_charge_id, "in_service"] = False
            self.net.flow_control.at[self.Flow_control_discharge_storage_id, "in_service"] = True  
            self.net.flow_control.at[self.Flow_control_discharge_id, "in_service"] = True
            self.net.flow_control.at[self.Flow_control_bypass_discharge_id, "in_service"] = False
            self.net.circ_pump_mass.at[self.Circ_pump_charge_storage_id, "mdot_flow_kg_per_s"] = 0
            self.net.flow_control.at[self.Flow_control_bypass_charge_id, "controlled_mdot_kg_per_s"] = 0
            self.net.circ_pump_mass.at[self.Circ_pump_discharge_storage_id, "mdot_flow_kg_per_s"] = abs(self.mass_flow)
            self.net.flow_control.at[self.Flow_control_discharge_storage_id, "controlled_mdot_kg_per_s"] = abs(self.mass_flow)
            self.net.flow_control.at[self.Flow_control_discharge_id, "controlled_mdot_kg_per_s"] =abs(self.mass_flow)
            self.net.circ_pump_mass.at[self.Circ_pump_discharge_storage_id, "t_flow_k"] = self.T_hot
            pp.pipeflow(self.net, mode="bidirectional")
            Tnet=self.net.res_flow_control.at[self.Flow_control_discharge_id, "t_from_k"]
        return Tnet
    def evaluate_bypass(self, dt_s,T_net): 
        bypass = False
        mdot_entering = 0.0
        mdot_bypass = 0.0
        v_in=0.0
        if self.mass_flow > 0: # Charge  
            if T_net < self.T_hot:
                self.direction = "to_Tcold"
                mdot_entering=self.mass_flow
            else:
                self.direction = "to_Thot"
                v_available = max(0.0, self.V_cold - self.V_MIN)
                Dv_requested = (self.mass_flow * dt_s) / self.density
                v_in = min(Dv_requested, v_available)
                mdot_entering = (v_in * self.density) / dt_s if dt_s > 0 else 0.0
                mdot_bypass = self.mass_flow - mdot_entering
                if mdot_bypass > 1e-6 or v_available <= 0:
                    bypass = True
                    print(f"  ⚠ [{self.name}] Tanque saturado de calor: "f"{Dv_requested - v_in:.2f} m3 no se pudieron cargar")
                    print(f"  ⚠ [{self.name}] Tanque saturado de calor: "f"{mdot_bypass:.2f} kg/s no se pudieron cargar")
        elif self.mass_flow < 0: # Discharge
            self.direction = "discharge"
            abs_mass_flow = abs(self.mass_flow)
            v_available = max(0.0, self.V_hot - self.V_MIN)
            Dv_requested = (abs_mass_flow * dt_s) / self.density
            v_in = min(Dv_requested, v_available)
            mdot_entering = (v_in * self.density) / dt_s if dt_s > 0 else 0.0
            mdot_bypass = abs_mass_flow - mdot_entering
            if mdot_bypass > 1e-6 or v_available <= 0:
                bypass = True
                print(f"  ⚠ [{self.name}] Tanque saturado de frio: "f"{Dv_requested - v_in:.2f} m3 no se pudieron descargar")
                print(f"  ⚠ [{self.name}] Tanque saturado de frio: "f"{mdot_bypass:.2f} kg/s no se pudieron descargar")
        elif self.mass_flow == 0:
            self.direction = "static"
        self.bypass = bypass
        self.mdot_entering = mdot_entering
        self.mdot_bypass = mdot_bypass
        self.v_in=v_in
    def apply_control_settings(self):
        if self.mass_flow > 0:
            #Deactivate and activate flow controls
            self.net.flow_control.at[self.Flow_control_discharge_storage_id, "in_service"] = False  
            self.net.flow_control.at[self.Flow_control_discharge_id, "in_service"] = False
            self.net.flow_control.at[self.Flow_control_bypass_discharge_id, "in_service"] = False  
            self.net.flow_control.at[self.Flow_control_charge_storage_id, "in_service"] = True  
            self.net.flow_control.at[self.Flow_control_charge_id, "in_service"] = True
            self.net.flow_control.at[self.Flow_control_bypass_charge_id, "in_service"] = True
            self.net.circ_pump_mass.at[self.Circ_pump_discharge_storage_id, "mdot_flow_kg_per_s"] = 0
            self.net.flow_control.at[self.Flow_control_bypass_charge_id, "controlled_mdot_kg_per_s"] = self.mdot_bypass
            self.net.circ_pump_mass.at[self.Circ_pump_charge_storage_id, "mdot_flow_kg_per_s"] = self.mdot_entering
            self.net.flow_control.at[self.Flow_control_charge_storage_id, "controlled_mdot_kg_per_s"] =self.mdot_entering
            self.net.flow_control.at[self.Flow_control_charge_id, "controlled_mdot_kg_per_s"] =self.mass_flow
            self.net.circ_pump_mass.at[self.Circ_pump_charge_storage_id, "t_flow_k"] = self.T_cold
        else: #Discharging 
            self.net.flow_control.at[self.Flow_control_charge_storage_id, "in_service"] = False  
            self.net.flow_control.at[self.Flow_control_charge_id, "in_service"] = False
            self.net.flow_control.at[self.Flow_control_bypass_charge_id, "in_service"] = False
            self.net.flow_control.at[self.Flow_control_discharge_storage_id, "in_service"] = True  
            self.net.flow_control.at[self.Flow_control_discharge_id, "in_service"] = True
            self.net.flow_control.at[self.Flow_control_bypass_discharge_id, "in_service"] = True 
            self.net.circ_pump_mass.at[self.Circ_pump_charge_storage_id, "mdot_flow_kg_per_s"] = 0
            self.net.flow_control.at[self.Flow_control_bypass_charge_id, "controlled_mdot_kg_per_s"] = 0
            self.net.flow_control.at[self.Flow_control_bypass_discharge_id, "controlled_mdot_kg_per_s"] = self.mdot_bypass
            self.net.circ_pump_mass.at[self.Circ_pump_discharge_storage_id, "mdot_flow_kg_per_s"] = self.mdot_entering
            self.net.flow_control.at[self.Flow_control_discharge_storage_id, "controlled_mdot_kg_per_s"] = self.mdot_entering
            self.net.flow_control.at[self.Flow_control_discharge_id, "controlled_mdot_kg_per_s"] =abs(self.mass_flow)
            self.net.circ_pump_mass.at[self.Circ_pump_discharge_storage_id, "t_flow_k"] = self.T_hot

    def Estimate_loss_ambient(self, T_amb):
        self.UA_hot_layer = self.UA_loss * (self.V_hot / self.V_tot)
        self.UA_cold_layer = self.UA_loss * (self.V_cold / self.V_tot)
        Q_loss_hot = self.UA_hot_layer * (self.T_hot - T_amb)
        Q_loss_cold = self.UA_cold_layer * (self.T_cold - T_amb)
        return Q_loss_hot, Q_loss_cold

    def V_dis(self, dt_s, T_amb,T_net):
        self.UA_hot_layer = self.UA_loss * (self.V_hot / self.V_tot)
        self.UA_cold_layer = self.UA_loss * (self.V_cold / self.V_tot)
        if self.mass_flow >= 0:
            if self.direction=='to_Tcold':
                self.T_hot += dt_s * (- ((self.UA_hot_layer * (self.T_hot - T_amb)) / (self.density * self.V_hot * self.heat_capacity)))
                self.T_cold += dt_s * (((self.mdot_entering * (T_net - self.T_cold)) / (self.density * self.V_cold)) - ((self.UA_cold_layer * (self.T_cold - T_amb)) / (self.density * self.V_cold * self.heat_capacity)))     
            elif self.direction=='to_Thot':
                self.T_hot += dt_s * (((self.mdot_entering * (T_net - self.T_hot)) / (self.density * self.V_hot)) - ((self.UA_hot_layer * (self.T_hot - T_amb)) / (self.density * self.V_hot * self.heat_capacity)))
                self.T_cold += dt_s * (-self.UA_cold_layer * (self.T_cold - T_amb) / (self.density * self.V_cold * self.heat_capacity))
                self.V_hot += self.v_in
                self.V_cold -= self.v_in
            else:
                self.T_hot += dt_s * (-self.UA_hot_layer * (self.T_hot - T_amb) / (self.density * self.V_hot * self.heat_capacity))
                self.T_cold += dt_s * (-self.UA_cold_layer * (self.T_cold - T_amb) / (self.density * self.V_cold * self.heat_capacity))   
        else:          
            self.T_cold += dt_s * ((self.mdot_entering * (T_net - self.T_cold) / (self.density * self.V_cold)) - (self.UA_cold_layer * (self.T_cold - T_amb) / (self.density * self.V_cold * self.heat_capacity)))
            self.T_hot += dt_s * (-(self.UA_hot_layer * (self.T_hot - T_amb) / (self.density * self.V_hot * self.heat_capacity)))
            self.V_cold += self.v_in
            self.V_hot -= self.v_in
        self.v_hot_fraction = self.V_hot / self.V_tot



# %%



