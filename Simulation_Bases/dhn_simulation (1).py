"""
Simulación horaria de una red 5GDHC (pandapipes + tespy):
    - N HeatPumpActor      (bombas de calor reversibles, ciclo real en TESPy)
    - 1 StorageActor       (tanque de estratificación / thermocline)
    - Coordinador          Gauss-Seidel con relajación, 1 pipeflow por iteración

Fusiona:
    - network_simulation.py       -> arquitectura de actores
    - multi_plant_coordinator.py  -> física real (ThermoclineStorage + TESPy HP)

Sin CapacityManager: se asume que la red, por diseño (dimensionado,
consignas de Q_need_profile y mass_flow_profile), no puede sufrir una
inversión de sentido. Las HP inyectan/extraen directamente lo que les pide
su perfil horario, sin limitar ni redirigir excedente hacia el storage.

===========================================================================
CORRECCIONES FÍSICAS respecto a los originales (importante leerlas)
===========================================================================

1) BYPASS DEL STORAGE CALCULADO DENTRO DE LA ITERACIÓN, NO DESPUÉS.
   En multi_plant_coordinator.py, `set_boundary()` asumía durante TODA la
   iteración de Gauss-Seidel que el 100% del caudal entraba al tanque
   (sin bypass), y solo AL FINAL, una vez convergido, se llamaba a
   `evaluate_bypass()` + `apply_control_settings()` para calcular el
   reparto real entrada/bypass -- pero sin volver a resolver pandapipes.
   Resultado: si el tanque estaba cerca de saturarse, la red quedaba en un
   estado (temperaturas, caudales) que NUNCA fue verificado por un
   pipeflow convergido; la "convergencia" reportada era la de un problema
   ligeramente distinto del que realmente se aplicaba a la red.
   Aquí, `evaluate_bypass()` + la escritura de fronteras se hacen en
   CADA iteración, con la T adivinada de esa iteración. Esto es correcto
   porque V_hot/V_cold (lo único que limita el bypass) no cambia dentro de
   la misma hora -- solo cambia entre horas, tras `finalize_hour()`.

2) UNIDADES DE TEMPERATURA MEZCLADAS ENTRE PANDAPIPES/TESPy.
   `Bidirectional_W_to_WHeatPump` usaba `T_estimate = 280.15` como valor
   por defecto y lo pasaba directo a `solve_cycle()` -> TESPy, pero la red
   TESPy está configurada en grados Celsius (`temperature="degC"`).
   280.15 interpretado como °C es una temperatura de red disparatada
   (debía ser ~7 °C, no 280 °C). Es decir, se estaba mezclando un valor en
   Kelvin (convención de pandapipes/ThermoclineStorage) con una API que
   espera Celsius (TESPy). Aquí el coordinador trabaja SIEMPRE en Kelvin
   (para casar con pandapipes) y se convierte explícitamente a Celsius
   justo antes de llamar a TESPy (`k_to_c`), y no antes.

3) EL ESTADO FÍSICO DEL TANQUE (volúmenes, temperaturas) SOLO SE ACTUALIZA
   UNA VEZ, AL CONVERGER LA HORA. Esto ya era así en el original y sigue
   siendo correcto: son estados que evolucionan en el tiempo, no
   condiciones de frontera que deban re-evaluarse en cada iteración.
   Aquí, antes de aplicar `V_dis()`, se hace una última llamada a
   `evaluate_bypass()` con la T ya convergida (no la última T adivinada),
   por consistencia -- la diferencia es < tol_k pero es gratis y evita
   arrastrar un pequeño sesgo sistemático.

4) SIN GESTIÓN DE CAPACIDAD: si el `mass_flow_profile` del storage y los
   `Q_need_profile` de las HP piden, en algún instante, más de lo que la
   red/tanque pueden dar físicamente (p. ej. tanque ya saturado, caudal
   pedido negativo imposible), el `pipeflow` puede no converger o devolver
   un resultado sin sentido físico (T fuera de rango, caudal invertido).
   Como no hay curtailment ni redirección automática, esa responsabilidad
   queda enteramente en cómo se construyen los perfiles horarios de entrada.
"""

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Optional
import numpy as np
import pandapipes as pp
from pandapipes.pf.pipeflow_setup import get_fluid
from tespy.components import Compressor, HeatExchanger, CycleCloser, Valve, Condenser, Sink, Source
from tespy.connections import Connection
from tespy.networks import Network
from tespy.tools import UserDefinedEquation


def k_to_c(t_k: float) -> float:
    return t_k - 273.15


def c_to_k(t_c: float) -> float:
    return t_c + 273.15


# ===========================================================================
# 1. INTERFAZ COMÚN DE ACTOR
#    El coordinador solo conoce prepare() / read_result(); nunca sabe si
#    detrás hay una bomba de calor TESPy o un tanque de estratificación.
#    Todas las temperaturas que cruzan esta interfaz están en KELVIN.
# ===========================================================================
class NetworkActor(ABC):
    name: str

    @abstractmethod
    def prepare(self, t_guess_k: float, hour: int) -> None:
        """Decide y escribe su condición de frontera en `net`. NO resuelve pipeflow."""
        raise NotImplementedError

    @abstractmethod
    def read_result(self, net) -> float:
        """Lee, tras el pipeflow, la T real (K) en su nodo de conexión con la red."""
        raise NotImplementedError

    def log_hour(self, hour: int, t_real_k: float) -> None:
        pass


# ===========================================================================
# 2. TANQUE DE ESTRATIFICACIÓN (thermocline) -- física original conservada
# ===========================================================================
class ThermoclineStorage:
    def __init__(self, net, name, flow_control_charge_id, flow_control_bypass_charge_id,
                 circ_pump_charge_storage_id, flow_control_charge_storage_id,
                 flow_control_discharge_id, flow_control_bypass_discharge_id,
                 circ_pump_discharge_storage_id, flow_control_discharge_storage_id,
                 volume_m3, t_hot_init_k, t_cold_init_k, v_hot_fraction_init, ua_loss, ua_interface=0.0):
        self.net = net
        self.fluid = get_fluid(net)
        self.name = name
        self.fc_charge = flow_control_charge_id
        self.fc_bypass_charge = flow_control_bypass_charge_id
        self.cp_charge_storage = circ_pump_charge_storage_id
        self.fc_charge_storage = flow_control_charge_storage_id
        self.fc_discharge = flow_control_discharge_id
        self.fc_bypass_discharge = flow_control_bypass_discharge_id
        self.cp_discharge_storage = circ_pump_discharge_storage_id
        self.fc_discharge_storage = flow_control_discharge_storage_id
        self.V_tot = volume_m3
        self.T_hot = t_hot_init_k       # Kelvin
        self.T_cold = t_cold_init_k     # Kelvin
        self.v_hot_fraction = v_hot_fraction_init
        self.UA_loss = ua_loss
        self.UA_interface = ua_interface
        self.V_hot = volume_m3 * v_hot_fraction_init
        self.V_cold = volume_m3 - self.V_hot
        self.V_MIN = volume_m3 * 0.01
        self.bypass = False
        self.direction = None
        self.mass_flow = 0.0
        self.mdot_entering = 0.0
        self.mdot_bypass = 0.0
        self.v_in = 0.0
        self.T_estimate = t_cold_init_k  # grado de libertad del coordinador (Kelvin)
        self._density_average_heat_capacity()

    def _density_average_heat_capacity(self):
        T_min_K, T_max_K = 273.15 + 20, 273.15 + 100
        T_samples = np.linspace(T_min_K, T_max_K, 100)
        self.density = np.mean([self.fluid.get_density(T) for T in T_samples])
        self.heat_capacity = np.mean([self.fluid.get_heat_capacity(T) for T in T_samples])

    def evaluate_bypass(self, dt_s: float, T_net_k: float) -> None:
        """Calcula cuánto caudal entra realmente al tanque y cuánto hace bypass,
        según la capacidad de volumen disponible EN ESTE INSTANTE (no cambia
        dentro de la misma hora). Puede llamarse en cada iteración de G-S."""
        bypass, mdot_entering, mdot_bypass, v_in = False, 0.0, 0.0, 0.0
        if self.mass_flow > 0:  # Charge
            if T_net_k < self.T_hot:
                self.direction = "to_Tcold"
                mdot_entering = self.mass_flow
            else:
                self.direction = "to_Thot"
                v_available = max(0.0, self.V_cold - self.V_MIN)
                Dv_requested = (self.mass_flow * dt_s) / self.density
                v_in = min(Dv_requested, v_available)
                mdot_entering = (v_in * self.density) / dt_s if dt_s > 0 else 0.0
                mdot_bypass = self.mass_flow - mdot_entering
                bypass = mdot_bypass > 1e-6 or v_available <= 0
        elif self.mass_flow < 0:  # Discharge
            self.direction = "discharge"
            abs_mf = abs(self.mass_flow)
            v_available = max(0.0, self.V_hot - self.V_MIN)
            Dv_requested = (abs_mf * dt_s) / self.density
            v_in = min(Dv_requested, v_available)
            mdot_entering = (v_in * self.density) / dt_s if dt_s > 0 else 0.0
            mdot_bypass = abs_mf - mdot_entering
            bypass = mdot_bypass > 1e-6 or v_available <= 0
        else:
            self.direction = "static"

        self.bypass, self.mdot_entering, self.mdot_bypass, self.v_in = bypass, mdot_entering, mdot_bypass, v_in

    def write_boundary(self) -> None:
        """Escribe en `net` el reparto entrada/bypass YA calculado por
        evaluate_bypass(). Reemplaza los antiguos set_boundary()+
        apply_control_settings() separados: ahora es una única función,
        llamada en cada iteración con el reparto real, nunca aproximado."""
        if self.mass_flow >= 0:
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
            self.net.flow_control.at[self.fc_charge, "controlled_mdot_kg_per_s"] = self.mass_flow
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
            self.net.flow_control.at[self.fc_discharge, "controlled_mdot_kg_per_s"] = abs(self.mass_flow)
            self.net.circ_pump_mass.at[self.cp_discharge_storage, "t_flow_k"] = self.T_hot

    def read_result(self) -> float:
        fc = self.fc_charge if self.mass_flow >= 0 else self.fc_discharge
        return self.net.res_flow_control.at[fc, "t_from_k"]

    def V_dis(self, dt_s: float, T_amb_k: float, T_net_k: float) -> None:
        UA_hot = self.UA_loss * (self.V_hot / self.V_tot)
        UA_cold = self.UA_loss * (self.V_cold / self.V_tot)
        if self.mass_flow >= 0:
            if self.direction == "to_Tcold":
                self.T_hot += dt_s * (-((UA_hot * (self.T_hot - T_amb_k)) / (self.density * self.V_hot * self.heat_capacity)))
                self.T_cold += dt_s * (((self.mdot_entering * (T_net_k - self.T_cold)) / (self.density * self.V_cold))
                                        - ((UA_cold * (self.T_cold - T_amb_k)) / (self.density * self.V_cold * self.heat_capacity)))
            elif self.direction == "to_Thot":
                self.T_hot += dt_s * (((self.mdot_entering * (T_net_k - self.T_hot)) / (self.density * self.V_hot))
                                       - ((UA_hot * (self.T_hot - T_amb_k)) / (self.density * self.V_hot * self.heat_capacity)))
                self.T_cold += dt_s * (-UA_cold * (self.T_cold - T_amb_k) / (self.density * self.V_cold * self.heat_capacity))
                self.V_hot += self.v_in
                self.V_cold -= self.v_in
            else:
                self.T_hot += dt_s * (-UA_hot * (self.T_hot - T_amb_k) / (self.density * self.V_hot * self.heat_capacity))
                self.T_cold += dt_s * (-UA_cold * (self.T_cold - T_amb_k) / (self.density * self.V_cold * self.heat_capacity))
        else:
            self.T_cold += dt_s * ((self.mdot_entering * (T_net_k - self.T_cold) / (self.density * self.V_cold))
                                    - (UA_cold * (self.T_cold - T_amb_k) / (self.density * self.V_cold * self.heat_capacity)))
            self.T_hot += dt_s * (-(UA_hot * (self.T_hot - T_amb_k) / (self.density * self.V_hot * self.heat_capacity)))
            self.V_cold += self.v_in
            self.V_hot -= self.v_in
        self.v_hot_fraction = self.V_hot / self.V_tot


# ===========================================================================
# 3. BOMBA DE CALOR REVERSIBLE (TESPy) -- física original conservada,
#    la única API nueva es set_boundary()/read_result() (ya existía) y
#    la conversión K<->C explícita en el punto de contacto con TESPy.
# ===========================================================================
class BidirectionalWtoWHeatPump:
    def __init__(self, name, net, refrigerant, hc_ext_id, hc_inj_id):
        self.name = name
        self.refrigerant = refrigerant
        self.HC_ext_id = hc_ext_id
        self.HC_inj_id = hc_inj_id
        self.net = net
        self._build_tespy_cycles()
        self.T_estimate_c = 7.0        # °C -- SOLO para TESPy, nunca cruza al coordinador
        self.mode = None
        self.Q_consumer = 0.0
        self.eta_s = 0.7
        self.T_cons = [8.0, 3.0]
        self.dt_DHN = 3.0
        self.dT_water = 4.0

    def _build_tespy_cycles(self):
        def my_ude(ude):
            return ude.conns[0].calc_T() + ude.params["dt"] - ude.conns[1].calc_T()

        def my_ude_dependents(ude):
            c1, c2 = ude.conns
            return [c1.p, c1.h, c2.p, c2.h]

        self.nw_cooling = Network()
        self.nw_cooling.units.set_defaults(temperature="degC", pressure="bar", pressure_difference="bar",
                                            enthalpy="J/kg", heat="W", power="W")
        comp, cond, valve, evap = Compressor("compresor"), Condenser("condensador"), Valve("valvula"), HeatExchanger("evaporador")
        cc = CycleCloser("CycleCloser")
        src_cons, snk_cons = Source("Source_Consumer"), Sink("Sink_consumer")
        src_red, snk_red = Source("Source_Reseau"), Sink("Sink_Reseau")
        self._cool = dict(comp=comp, cond=cond, valve=valve, evap=evap)
        c0 = Connection(valve, "out1", cc, "in1", label="0")
        c1 = Connection(cc, "out1", evap, "in2", label="1")
        c2 = Connection(evap, "out2", comp, "in1", label="2")
        c3 = Connection(comp, "out1", cond, "in1", label="3")
        c4 = Connection(cond, "out1", valve, "in1", label="4")
        c5 = Connection(cond, "out2", snk_cons, "in1", label="5")
        c6 = Connection(src_cons, "out1", cond, "in2", label="6")
        c7 = Connection(src_red, "out1", evap, "in1", label="7")
        c8 = Connection(evap, "out1", snk_red, "in1", label="8")
        self.nw_cooling.add_conns(c0, c1, c2, c3, c4, c5, c6, c7, c8)
        self._cool.update(c2=c2, c5=c5, c6=c6, c7=c7)
        ude = UserDefinedEquation("ude_cool", my_ude, my_ude_dependents, conns=[c8, c7], params={"dt": 5})
        self.nw_cooling.add_ude(ude)
        self._cool["ude"] = ude

        self.nw_heating = Network()
        self.nw_heating.units.set_defaults(temperature="degC", pressure="bar", pressure_difference="bar",
                                            enthalpy="J/kg", heat="W", power="W")
        comp2, cond2, valve2, evap2 = Compressor("compresor"), Condenser("condensador"), Valve("valvula"), HeatExchanger("evaporador")
        cc2 = CycleCloser("CycleCloser")
        src_cons2, snk_cons2 = Source("Source_Consumer"), Sink("Sink_consumer")
        src_red2, snk_red2 = Source("Source_Reseau"), Sink("Sink_Reseau")
        self._heat = dict(comp=comp2, cond=cond2, valve=valve2, evap=evap2)
        h0 = Connection(valve2, "out1", cc2, "in1", label="0")
        h1 = Connection(cc2, "out1", evap2, "in2", label="1")
        h2 = Connection(evap2, "out2", comp2, "in1", label="2")
        h3 = Connection(comp2, "out1", cond2, "in1", label="3")
        h4 = Connection(cond2, "out1", valve2, "in1", label="4")
        h5 = Connection(evap2, "out1", snk_cons2, "in1", label="5")
        h6 = Connection(src_cons2, "out1", evap2, "in1", label="6")
        h7 = Connection(src_red2, "out1", cond2, "in2", label="7")
        h8 = Connection(cond2, "out2", snk_red2, "in1", label="8")
        self.nw_heating.add_conns(h0, h1, h2, h3, h4, h5, h6, h7, h8)
        self._heat.update(h2=h2, h5=h5, h6=h6, h7=h7)
        ude2 = UserDefinedEquation("ude_heat", my_ude, my_ude_dependents, conns=[h7, h8], params={"dt": 5})
        self.nw_heating.add_ude(ude2)
        self._heat["ude"] = ude2

    def solve_cycle(self, mode, Q_consumer, eta_s, T_network_in_c, T_cons, dt_DHN):
        """T_network_in_c: temperatura de la red en °C (¡no Kelvin!)."""
        self.mode = mode
        T_cons_in, T_cons_out = T_cons
        import CoolProp.CoolProp as CP
        T_triple = CP.Props1SI("Ttriple", self.refrigerant)
        p_triple = CP.Props1SI("ptriple", self.refrigerant)
        T_critical = CP.Props1SI("T_critical", self.refrigerant)
        p_critical = CP.Props1SI("p_critical", self.refrigerant)
        h_min = CP.PropsSI("H", "T", T_triple + 0.1, "Q", 0, self.refrigerant)
        p_high = min(p_critical * 0.9, 30e5)
        h_max = CP.PropsSI("H", "T", T_critical * 0.9, "P", p_high, self.refrigerant)

        if mode == "COOLING_NET":
            g = self._cool
            g["evap"].set_attr(pr1=1, pr2=1, ttd_l=5)
            g["cond"].set_attr(pr1=1, pr2=1, ttd_u=5, Q=Q_consumer)
            g["comp"].set_attr(eta_s=eta_s)
            g["c2"].set_attr(fluid={self.refrigerant: 1}, td_dew=5)
            g["c5"].set_attr(T=T_cons_out, p=3, fluid={"water": 1})
            g["c6"].set_attr(T=T_cons_in)
            g["c7"].set_attr(T=T_network_in_c, p=2.5, fluid={"water": 1})
            self.nw_cooling._set_p_range([p_triple, p_high])
            self.nw_cooling._set_h_range([h_min, h_max])
            g["ude"].params["dt"] = dt_DHN
            self.nw_cooling.solve("design")
        elif mode == "HEATING_NET":
            g = self._heat
            g["evap"].set_attr(pr1=1, pr2=1, ttd_l=5, Q=Q_consumer)
            g["cond"].set_attr(pr1=1, pr2=1, ttd_u=5)
            g["comp"].set_attr(eta_s=eta_s)
            g["c2"].set_attr(fluid={self.refrigerant: 1}, td_dew=5)
            g["c5"].set_attr(T=T_cons_out, p=3, fluid={"water": 1})
            g["c6"].set_attr(T=T_cons_in)
            g["c7"].set_attr(T=T_network_in_c, p=2.5, fluid={"water": 1})
            self.nw_heating._set_p_range([p_triple, p_high])
            self.nw_heating._set_h_range([h_min, h_max])
            g["ude"].params["dt"] = dt_DHN
            self.nw_heating.solve("design")

    def interaction_simulation(self, dT_water):
        if self.mode == "COOLING_NET":
            self.net.heat_consumer.at[self.HC_ext_id, "qext_w"] = abs(self._cool["evap"].Q.val)
            self.net.heat_consumer.at[self.HC_ext_id, "deltat_k"] = dT_water
            self.net.heat_consumer.at[self.HC_ext_id, "in_service"] = True
            self.net.heat_consumer.at[self.HC_inj_id, "in_service"] = False
        elif self.mode == "HEATING_NET":
            self.net.heat_consumer.at[self.HC_inj_id, "qext_w"] = self._heat["cond"].Q.val
            self.net.heat_consumer.at[self.HC_inj_id, "deltat_k"] = dT_water
            self.net.heat_consumer.at[self.HC_inj_id, "in_service"] = True
            self.net.heat_consumer.at[self.HC_ext_id, "in_service"] = False
        else:
            self.net.heat_consumer.at[self.HC_ext_id, "in_service"] = False
            self.net.heat_consumer.at[self.HC_inj_id, "in_service"] = False

    def get_main_results(self) -> dict:
        if self.mode == "COOLING_NET":
            g = self._cool
            return {"Q_consumer": abs(g["cond"].Q.val), "Q_network": abs(g["evap"].Q.val),
                    "W_compressor": g["comp"].P.val, "COP": abs(g["cond"].Q.val) / g["comp"].P.val}
        if self.mode == "HEATING_NET":
            g = self._heat
            return {"Q_consumer": abs(g["evap"].Q.val), "Q_network": abs(g["cond"].Q.val),
                    "W_compressor": g["comp"].P.val, "COP": abs(g["evap"].Q.val) / g["comp"].P.val}
        return {}

    def set_boundary(self) -> None:
        if self.mode is None:
            self.net.heat_consumer.at[self.HC_ext_id, "in_service"] = False
            self.net.heat_consumer.at[self.HC_inj_id, "in_service"] = False
            return
        self.solve_cycle(self.mode, self.Q_consumer, self.eta_s, self.T_estimate_c, self.T_cons, self.dt_DHN)
        self.interaction_simulation(self.dT_water)

    def read_result(self) -> float:
        """Devuelve la T en KELVIN (convención pandapipes), no en °C."""
        if self.mode is None:
            return c_to_k(self.T_estimate_c)
        active_id = self.HC_inj_id if self.mode == "HEATING_NET" else self.HC_ext_id
        return self.net.res_heat_consumer.at[active_id, "t_from_k"]


# ===========================================================================
# 4. ACTORES (envoltorio fino sobre la física de arriba para el coordinador)
# ===========================================================================
@dataclass
class StorageActor(NetworkActor):
    name: str
    storage: ThermoclineStorage
    mass_flow_profile: list          # perfil horario (kg/s, +carga/-descarga)
    dt_s: float = 3600.0
    log: dict = field(default_factory=dict)

    def prepare(self, t_guess_k: float, hour: int) -> None:
        self.storage.mass_flow = self.mass_flow_profile[hour]
        self.storage.T_estimate = t_guess_k
        self.storage.evaluate_bypass(dt_s=self.dt_s, T_net_k=t_guess_k)
        self.storage.write_boundary()

    def read_result(self, net) -> float:
        return self.storage.read_result()

    def log_hour(self, hour: int, t_real_k: float) -> None:
        self.log[hour] = {
            "T_source_k": t_real_k, "mode": self.storage.direction,
            "mdot_bypass": self.storage.mdot_bypass,
        }

    def finalize_hour(self, t_amb_k: float, t_real_k: float) -> None:
        """Se llama UNA vez, tras converger la hora, para avanzar el estado
        físico del tanque (volúmenes y temperaturas). Recalcula el reparto
        con la T ya convergida (t_real_k), no con la última T adivinada."""
        self.storage.evaluate_bypass(dt_s=self.dt_s, T_net_k=t_real_k)
        self.storage.V_dis(dt_s=self.dt_s, T_amb_k=t_amb_k, T_net_k=t_real_k)


@dataclass
class HeatPumpActor(NetworkActor):
    name: str
    hp: BidirectionalWtoWHeatPump
    Q_need_profile: list             # perfil horario [W], + calefacción / - refrigeración
    eta_s: float = 0.7
    T_cons: list = field(default_factory=lambda: [8.0, 3.0])
    dt_DHN: float = 3.0
    dT_water: float = 4.0
    log: dict = field(default_factory=dict)

    def prepare(self, t_guess_k: float, hour: int) -> None:
        Q_need = self.Q_need_profile[hour]

        self.hp.mode = "HEATING_NET" if Q_need >= 0 else "COOLING_NET"
        self.hp.Q_consumer = abs(Q_need)
        self.hp.eta_s = self.eta_s
        self.hp.T_cons = self.T_cons
        self.hp.dt_DHN = self.dt_DHN
        self.hp.dT_water = self.dT_water
        self.hp.T_estimate_c = k_to_c(t_guess_k)   # única conversión K->C del sistema
        self.hp.set_boundary()

    def read_result(self, net) -> float:
        return self.hp.read_result()

    def log_hour(self, hour: int, t_real_k: float) -> None:
        self.log[hour] = {"T_source_k": t_real_k, **self.hp.get_main_results()}


# ===========================================================================
# 5. COORDINADOR: Gauss-Seidel con relajación, 1 solo pipeflow/iteración
# ===========================================================================
class ConvergenceError(RuntimeError):
    pass


def solve_hour(net, actors: list[NetworkActor], T_prev: dict, hour: int,
               tol_k: float = 0.05, max_iter: int = 25, relax: float = 0.5,
               verbose: bool = False) -> tuple[dict, int]:
    """Todas las T en Kelvin. El orden de `actors` en la lista es el orden
    de barrido de Gauss-Seidel (cada actor ve, dentro de la misma iteración,
    la T adivinada -- no depende de otros actores, ya que no hay reparto
    de excedente entre ellos)."""
    T_guess = dict(T_prev)

    for it in range(1, max_iter + 1):
        for actor in actors:
            actor.prepare(T_guess[actor.name], hour)

        pp.pipeflow(net, mode="bidirectional")

        T_real = {a.name: a.read_result(net) for a in actors}
        error = max(abs(T_real[a.name] - T_guess[a.name]) for a in actors)

        if verbose:
            detalle = ", ".join(f"{a.name}={k_to_c(T_real[a.name]):.2f}°C" for a in actors)
            print(f"  hora {hour} iter {it}: error={error:.4f} K | {detalle}")

        if error < tol_k:
            for a in actors:
                a.log_hour(hour, T_real[a.name])
            return T_real, it

        T_guess = {a.name: relax * T_real[a.name] + (1 - relax) * T_guess[a.name] for a in actors}

    raise ConvergenceError(f"Hora {hour}: no convergió en {max_iter} iteraciones (error={error:.4f} K, tol={tol_k}).")


def simulate_year(net, actors: list[NetworkActor], storage_actor: Optional[StorageActor],
                   T_amb_profile_k: list, n_hours: int = 8760, T_init_k: float = 288.15,
                   tol_k: float = 0.05, max_iter: int = 25, relax: float = 0.5,
                   verbose_hours: bool = False) -> dict:
    T_prev = {a.name: T_init_k for a in actors}
    history = {a.name: {} for a in actors}
    n_iter_log = []

    for hour in range(n_hours):
        T_real, n_iter = solve_hour(net, actors, T_prev, hour, tol_k, max_iter, relax, verbose_hours)

        for a in actors:
            history[a.name][hour] = T_real[a.name]

        if storage_actor is not None:
            storage_actor.finalize_hour(t_amb_k=T_amb_profile_k[hour], t_real_k=T_real[storage_actor.name])

        n_iter_log.append(n_iter)
        T_prev = T_real

        if hour % 500 == 0:
            avg = sum(n_iter_log[-500:]) / len(n_iter_log[-500:])
            resumen = " | ".join(f"{k}: {k_to_c(v):5.2f}°C" for k, v in T_real.items())
            print(f"Hora {hour:5d}/{n_hours} | iter={n_iter} (media 500h={avg:.1f}) | {resumen}")

    print(f"\nMedia global: {sum(n_iter_log) / len(n_iter_log):.2f} iteraciones/hora en {n_hours} horas.")
    return history


# ===========================================================================
# 6. EJEMPLO DE USO (esqueleto -- ajustar IDs a tu red real)
# ===========================================================================
if __name__ == "__main__":
    # net = pp.create_empty_network(fluid="water")
    # ... construir junctions, pipes, heat_consumer, flow_control, circ_pump_mass, ext_grid ...
    #
    # storage_phys = ThermoclineStorage(
    #     net, "storage_1",
    #     flow_control_charge_id=..., flow_control_bypass_charge_id=...,
    #     circ_pump_charge_storage_id=..., flow_control_charge_storage_id=...,
    #     flow_control_discharge_id=..., flow_control_bypass_discharge_id=...,
    #     circ_pump_discharge_storage_id=..., flow_control_discharge_storage_id=...,
    #     volume_m3=50.0, t_hot_init_k=318.15, t_cold_init_k=283.15,
    #     v_hot_fraction_init=0.3, ua_loss=15.0,
    # )
    # storage = StorageActor(name="Storage", storage=storage_phys,
    #                         mass_flow_profile=mdot_storage_profile, dt_s=3600.0)
    #
    # hp1_phys = BidirectionalWtoWHeatPump("HP1", net, "R410A", hc_ext_id=0, hc_inj_id=1)
    # hp1 = HeatPumpActor(name="HP1", hp=hp1_phys, Q_need_profile=Q_need_hp1)
    #
    # actors = [hp1, storage]
    #
    # history = simulate_year(net, actors, storage_actor=storage,
    #                          T_amb_profile_k=[278.15] * 8760, n_hours=8760, T_init_k=288.15)
    print("Esqueleto listo. Rellena la sección 6 con tu red real.")
