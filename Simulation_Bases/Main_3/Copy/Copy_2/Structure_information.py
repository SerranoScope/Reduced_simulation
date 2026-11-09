#!/usr/bin/env python
# coding: utf-8

# This file will give the structure to build entirelyh a district heating network from a graph-representation. 
# 
# The physical connection of the entire district heating and cooling networks is described by a incidence matrix. Each node represent one element of a the district heating. 
# 
# For each element linked in the graph, there is a pair od pipes. One named high pressure line and the other low pressure line. The incidence matrix will described only the high pressure and the low pressure is that matrix *-1. 
# 
# There is going to be the option to build a a network from the shp file or from a incidence matrix. 
# 
# 

# In[393]:


import networkx as nx
import pandas as pd
import numpy as np
import pandapipes as pp 
from External_objects_functions import Bidirectional_W_to_WHeatPump,ThermoclineTwoLayer,ConsumerActor,Pipe_class,Chaufferie,MeritOrderDispatcher,Central_production,proportional_cut,priority_cut
import pandapipes_patch   # parche CirculationPump (in_service mixto) 
import time
try:
    from pandapipes.pf.pipeflow_setup import PipeflowNotConverged
except ImportError:                      # versiones de pandapipes donde vive en otro sitio
    PipeflowNotConverged = RuntimeError


# In[394]:


class ConvergenceError(Exception):
    """El bucle iterativo de temperaturas no convergió."""
    pass


# In[395]:


def bold(text):
    out = []
    for ch in text:
        if "A" <= ch <= "Z":
            out.append(chr(0x1D400 + ord(ch) - ord("A")))
        elif "a" <= ch <= "z":
            out.append(chr(0x1D41A + ord(ch) - ord("a")))
        elif "0" <= ch <= "9":
            out.append(chr(0x1D7CE + ord(ch) - ord("0")))
        else:
            out.append(ch)
    return "".join(out)


# In[396]:


class DHN ():
    "This class will deal with the mathematical graph structure , the construction of the network in pandapipes and the instantiate of the "
    "external objects "
    def __init__(self,type_Graph,layout):
        self.type_Graph=type_Graph
        self.layout=layout
        # --- estados de cada hora (ver results / compare / export_states) ---
        self.cut_rule=proportional_cut      # regla de reparto cuando el calor no alcanza
        self.solve_cut=True                 # resolver el estado recortado en horas infactibles
        self.capture_states="infeasible"    # "none" | "infeasible" | "all"
        self.states={}                      # {hora: {"required": ..., "cut": ...}}
        self.cut_notes={}                   # {hora: por que no hay estado recortado}
        self.dispatcher=None 
    def Graph_From_Excel(self,excel_path,sheet=0,Node_Column=0): 
        """This function will create a graph representation from a incidenct matrix coming from an excel descrption"""
        Sheet_1=pd.read_excel(excel_path,sheet)
        Sheet_2=pd.read_excel(excel_path,sheet+1)
        Sheet_3=pd.read_excel(excel_path,sheet+2)
        Sheet_4=pd.read_excel(excel_path,sheet+3)
        Sheet_6=pd.read_excel(excel_path,sheet+5)
        # Graph characteristics
        self.Nodes=Sheet_1.iloc[:,Node_Column].tolist()
        self.Edges=list(Sheet_1.columns.values[1:])
        self.Type_Nodes=Sheet_2.iloc[0,1:].tolist()
        self.IM=Sheet_1.iloc[:,Node_Column+1:].to_numpy()
        #Pipe characteristics
        self.Pipes_charact=Sheet_3.iloc[:,Node_Column+1:]
        if self.Pipes_charact.empty: 
            pass
        else:
            lengths=self.Pipes_charact.iloc[0, :].tolist()
            inner_diameters=self.Pipes_charact.iloc[1, :].tolist()
            u_w_per_m2k=self.Pipes_charact.iloc[2, :].tolist()
            K_mm=self.Pipes_charact.iloc[3, :].tolist()
        #Chaufferies (una columna por chaufferie; varias pueden colgar del mismo nodo)
        ch = pd.read_excel(excel_path, sheet_name="Chaufferies", index_col=0)
        self.chaufferies_by_node = {}
        for cname in ch.columns:
            p = ch[cname]
            q_min = 0.0 if pd.isna(p["Q_min"]) else float(p["Q_min"])
            self.chaufferies_by_node.setdefault(p["Node"], []).append(dict(
                name=cname, node=p["Node"], Fuel_type=p["Fuel_type"], Q_nom=float(p["Q_nom"]), Q_min=q_min,
                Nominal_efficiency=float(p["Nominal_efficiency"]), Emission_factor=float(p["Emission_factor"]),
                cost_constant=float(p["cost_constant"])))
        #Storage characteristics
        self.Storage_characteristics={"node":[],"V_tot":[],"T_hot_init":[],"T_cold_init":[],"v_hot_fraction":[],"UA_loss":[],"UA_interface":[]}
        self.Storage_characteristics_excel=Sheet_4.iloc[:,Node_Column+1:]
        if self.Storage_characteristics_excel.empty:
            pass
        else:
            self.Storage_characteristics["node"]=self.Storage_characteristics_excel.columns.tolist()
            self.Storage_characteristics["V_tot"]=self.Storage_characteristics_excel.iloc[0,:].tolist()
            self.Storage_characteristics["T_hot_init"]=self.Storage_characteristics_excel.iloc[1, :].tolist()
            self.Storage_characteristics["T_cold_init"]=self.Storage_characteristics_excel.iloc[2, :].tolist()
            self.Storage_characteristics["v_hot_fraction"]=self.Storage_characteristics_excel.iloc[3, :].tolist()
            self.Storage_characteristics["UA_loss"]=self.Storage_characteristics_excel.iloc[4, :].tolist()
            self.Storage_characteristics["UA_interface"]=self.Storage_characteristics_excel.iloc[5, :].tolist()
        #HP characteristics 
        self.HP_characteristics={"node":[],"refrigerant ":[],"dT_water":[],"eta_s":[]}
        self.HP_characteristics_excel=Sheet_6.iloc[:,Node_Column+1:]
        if  self.HP_characteristics_excel.empty:
            pass
        else:
            self.HP_characteristics["node"]=self.HP_characteristics_excel.columns.tolist()
            self.HP_characteristics["refrigerant"]=self.HP_characteristics_excel.iloc[0,:].tolist()
            self.HP_characteristics["dT_water"]=self.HP_characteristics_excel.iloc[1, :].tolist()
            self.HP_characteristics["eta_s"]=self.HP_characteristics_excel.iloc[2, :].tolist()

        if self.type_Graph== "Undirected":
            self.graph=nx.Graph()
            for i in range(len(self.Nodes)):
                self.graph.add_node(self.Nodes[i],type=self.Type_Nodes[i])
            edges=[]
            for i in range(self.IM.shape[1]): 
                list_nodes=[]
                col=self.IM[:,i]
                for j in range (self.IM.shape[0]): 
                    if col[j] == 1:
                        list_nodes.append(self.Nodes[j])
                edges.append((list_nodes[0],list_nodes[1],{"length": lengths[i],"name":self.Edges[i],'inner_diameter':inner_diameters[i],'u_w_per_m2k':u_w_per_m2k[i],'K_mm:':K_mm[i]}))
        elif self.type_Graph== "Directed":
            self.graph=nx.DiGraph()
            for i in range(len(self.Nodes)):
                self.graph.add_node(self.Nodes[i],type=self.Type_Nodes[i])
            edges=[]
            for i in range(self.IM.shape[1]): 
                list_nodes=[]
                col=self.IM[:,i]
                tails=int(np.where(col == 1)[0][0])
                heads=int(np.where(col == -1)[0][0])
                edges.append((self.Nodes[tails],self.Nodes[heads],{"length": lengths[i],"name":self.Edges[i],'inner_diameter':inner_diameters[i],'u_w_per_m2k':u_w_per_m2k[i],'K_mm:':K_mm[i]}))
        self.graph.add_edges_from(edges)
        self.draw_graph_graphviz_from_excel()

    def graph_from_SHP(self,shp_path):
        """This function will create a graph from a shp file"""
        pass

    def draw_graph_graphviz_from_excel(self, filename="network.png"):
            color_map = {
                "plant": "#e74c3c",
                "node": "#95a5a6",
                "storage": "#f39c12",
                "prosumer": "#5dade2",
                "consumer": "#2ecc71",
                "plant_support": "#FF00F7"
            }

            A = nx.nx_agraph.to_agraph(self.graph)

            A.graph_attr.update(
                mode="KK",
                overlap="scale",
                splines="true",
                sep="+10",
                fontsize="12",
                bgcolor="white",
                dpi="150",
                pack="true",
            )

            A.node_attr.update(
                shape="ellipse",
                style="filled",
                fontname="Helvetica",

                fontsize="11",
                width="0.7",
                height="0.5",
            )

            A.edge_attr.update(
                fontname="Helvetica",
                fontsize="9",
                color="#34495e",
                penwidth="1.3",
                arrowsize="0.8",
                dectorator='True'
            )

            for node in self.graph.nodes():
                node_type = str(self.graph.nodes[node]["type"]).strip().lower()
                color = color_map.get(node_type, "#f1c40f")
                n = A.get_node(node)
                n.attr["fillcolor"] = color
                n.attr["fontcolor"] = "white" if node_type == "plant" else "black"

            for u, v in self.graph.edges():
                edge_name = self.graph.edges[u, v].get("name", "")
                real_length = self.graph.edges[u, v].get("length", 1)
                e = A.get_edge(u, v)
                e.attr["label"] = edge_name
                e.attr["len"] = str(float(real_length) / 20)

            # --- Leyenda como UN SOLO nodo con tabla HTML ---
            rows = "".join(
                f'<TR>'
                f'<TD BGCOLOR="{color}" WIDTH="18" HEIGHT="18"></TD>'
                f'<TD ALIGN="LEFT"><FONT POINT-SIZE="11">{type_name.capitalize()}</FONT></TD>'
                f'</TR>'
                for type_name, color in color_map.items()
            )
            legend_html = (
                '<<TABLE BORDER="1" CELLBORDER="0" CELLSPACING="4" CELLPADDING="4" '
                'COLOR="gray">'
                '<TR><TD COLSPAN="2"><FONT POINT-SIZE="13"><B>Leyenda</B></FONT></TD></TR>'
                f'{rows}'
                '</TABLE>>'
            )
            A.add_node(
                "legend_box",
                label=legend_html,
                shape="plain",
                margin="0",
                pin="true",
                pos="0,-2!",   # ajusta según el tamaño de tu red; ver nota abajo
            )

            A.layout(prog="neato")
            A.draw(filename)

    def create_pumps_centrals(self,Node):
        pp.create_circ_pump_const_pressure(
            self.net, return_junction=self.net.junction[self.net.junction["name"]==Node+ '_lowp'].index[0], 
            flow_junction=self.net.junction[self.net.junction["name"]==Node + '_highp'].index[0],
            p_flow_bar=5, plift_bar=0.5,
            t_flow_k=333, type="pt", in_service=False,
            name="Circ_pump_const_pressure_plant_"+Node)
        pp.create_circ_pump_const_mass_flow(
            self.net, return_junction=self.net.junction[self.net.junction["name"]==Node+ '_lowp'].index[0], 
            flow_junction=self.net.junction[self.net.junction["name"]==Node + '_highp'].index[0],
            mdot_flow_kg_per_s=0.1, t_flow_k=333,
            p_flow_bar=None, type="t", in_service=False,
            name="Circ_pump_mass_plant_"+Node)

    def create_pandapipes_plus_initialization(self,fluid="water"):
        self.net=pp.create_empty_network(fluid=fluid)
        # Creation of pairs of nodes
        for i in range(len(self.Nodes)):
            pp.create_junction(self.net,pn_bar=5,tfluid_k=353.15, name=self.Nodes[i] + '_highp')
            pp.create_junction(self.net,pn_bar=5,tfluid_k=353.15, name=self.Nodes[i] + '_lowp')
        #Creation of pipes 
        for u, v in self.graph.edges():
            name=self.graph.edges[u, v].get("name", "")
            length=self.graph.edges[u, v].get("length", 1)
            inner_diameters=self.graph.edges[u, v].get("inner_diameter", 1)
            u_w_per_m2k=self.graph.edges[u, v].get("u_w_per_m2k", 1)
            K_mm=self.graph.edges[u, v].get("K_mm", 0.1)
            pp.create_pipe_from_parameters(self.net, from_junction=self.net.junction[self.net.junction["name"]==u+ '_highp'].index[0], 
                                           to_junction=self.net.junction[self.net.junction["name"]==v+ '_highp'].index[0],length_km=length, 
                                           inner_diameter_mm=inner_diameters, u_w_per_m2k=u_w_per_m2k, name=name+"_supply",k_mm=K_mm,text_k=283.15)
            pp.create_pipe_from_parameters(self.net, from_junction=self.net.junction[self.net.junction["name"]==v+ '_lowp'].index[0], 
                                           to_junction=self.net.junction[self.net.junction["name"]==u+ '_lowp'].index[0],length_km=length
                                           , inner_diameter_mm=inner_diameters, u_w_per_m2k=u_w_per_m2k, name=name+"_return",k_mm=K_mm,text_k=283.15)

        #Creation of elements in pandapipes and elements
        for i in range(len(self.Nodes)):
            if self.Type_Nodes[i]=="Storage":
                pp.create_circ_pump_const_mass_flow(
                    self.net, return_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_highp'].index[0], 
                    flow_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_lowp'].index[0],
                    mdot_flow_kg_per_s=0.1, t_flow_k=300, p_flow_bar=None,type="t",name="Circ_pump_Storage_"+ self.Nodes[i])
            elif self.Type_Nodes[i]=="Plant_support" or self.Type_Nodes[i]=="Plant":
                self.create_pumps_centrals(self.Nodes[i])
            elif self.Type_Nodes[i]=="Prosumer":
                pp.create_heat_consumer(
                                    self.net,from_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_highp'].index[0],
                                    to_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_lowp'].index[0],
                                    qext_w=150000,deltat_k=10,name=self.Nodes[i]+ "_Prosumer")

            elif self.Type_Nodes[i]=="Consumer":
                pp.create_heat_consumer(
                    self.net,from_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_highp'].index[0],
                    to_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_lowp'].index[0],
                    qext_w=150000,treturn_k=305,name=self.Nodes[i]+ "_Consumer")

    def Creation_objects(self):
        self.Pipe_objects={"name":[],"Object":[]}
        self.Consumer_actors={"Node":[],"name":[],"Object":[]}
        self.HP_actors={"Node":[],"name":[],"Object":[]}
        self.Storage_actors={"Node":[],"name":[],"Object":[]}
        self.plant_actors={"Node":[],"name":[],"Object":[]}
        for i in range(len(self.Nodes)):
            if self.Type_Nodes[i]=="Storage":
                Storage=ThermoclineTwoLayer(net=self.net,name="Storage_"+ self.Nodes[i],
                                            Circ_pump_id=self.net.circ_pump_mass[self.net.circ_pump_mass["name"]=="Circ_pump_Storage_"+ self.Nodes[i]].index[0],
                                            j_supply=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_highp'].index[0],
                                            j_return=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_lowp'].index[0],
                                            V_tot=self.Storage_characteristics["V_tot"][self.Storage_characteristics["node"].index(self.Nodes[i])],
                                            T_hot=self.Storage_characteristics["T_hot_init"][self.Storage_characteristics["node"].index(self.Nodes[i])],
                                            T_cold=self.Storage_characteristics["T_cold_init"][self.Storage_characteristics["node"].index(self.Nodes[i])],
                                            v_hot_fraction=self.Storage_characteristics["v_hot_fraction"][self.Storage_characteristics["node"].index(self.Nodes[i])],
                                            UA_loss=self.Storage_characteristics["UA_loss"][self.Storage_characteristics["node"].index(self.Nodes[i])],
                                            UA_interface=self.Storage_characteristics["UA_interface"][self.Storage_characteristics["node"].index(self.Nodes[i])]
                                            ) 
                Storage._post_init()
                self.Storage_actors["Node"].append(self.Nodes[i])
                self.Storage_actors["name"].append(Storage.name)
                self.Storage_actors["Object"].append(Storage)
            elif self.Type_Nodes[i]=="Prosumer":
                HP=Bidirectional_W_to_WHeatPump(net=self.net,name="HP_"+ self.Nodes[i],refrigerant=self.HP_characteristics["refrigerant"][self.HP_characteristics["node"].index(self.Nodes[i])],
                                                HP_heat_consumer_id=self.net.heat_consumer[self.net.heat_consumer["name"]==self.Nodes[i]+ "_Prosumer"].index[0])
                HP._post_init()
                self.HP_actors["Node"].append(self.Nodes[i])
                self.HP_actors["name"].append(HP.name)
                self.HP_actors["Object"].append(HP)

            elif self.Type_Nodes[i] in ("Plant","Plant_support"):
                node=self.Nodes[i]
                specs=self.chaufferies_by_node.get(node,[])
                if not specs:
                    raise ValueError(f"Nodo {node} es {self.Type_Nodes[i]} pero no tiene columnas en la hoja 'Chaufferies'")
                plant=Central_production(net=self.net,name="Plant_"+node,node=node,
                                         Circ_mass_id=self.net.circ_pump_mass[self.net.circ_pump_mass["name"]=="Circ_pump_mass_plant_"+node].index[0],
                                         Circ_pressure_id=self.net.circ_pump_pressure[self.net.circ_pump_pressure["name"]=="Circ_pump_const_pressure_plant_"+node].index[0],
                                         Chaufferies=[Chaufferie(**sp) for sp in specs],
                                         is_reference=(self.Type_Nodes[i]=="Plant"))
                plant._post_init()
                self.plant_actors["Node"].append(node)
                self.plant_actors["name"].append(plant.name)
                self.plant_actors["Object"].append(plant)

            elif self.Type_Nodes[i]=="Consumer":
                Consumer=ConsumerActor(net=self.net,name=self.Nodes[i]+ "_Consumer",
                heat_consumer_id=self.net.heat_consumer[self.net.heat_consumer["name"]==self.Nodes[i]+ "_Consumer"].index[0])
                self.Consumer_actors["Node"].append(self.Nodes[i])
                self.Consumer_actors["name"].append(Consumer.name)
                self.Consumer_actors["Object"].append(Consumer)
        for u, v in self.graph.edges():
            name=self.graph.edges[u, v].get("name", "")
            Pipe_supply=Pipe_class(name=name+"_supply",net=self.net,
                             pipe_id=self.net.pipe[self.net.pipe["name"]==name+"_supply"].index[0])
            Pipe_return=Pipe_class(name=name+"_return",net=self.net,
                             pipe_id=self.net.pipe[self.net.pipe["name"]==name+"_return"].index[0])
            self.Pipe_objects["name"].append(Pipe_supply.name)
            self.Pipe_objects["name"].append(Pipe_return.name)
            self.Pipe_objects["Object"].append(Pipe_supply)
            self.Pipe_objects["Object"].append(Pipe_return)      
    def def_input_variables_simulation(self,input_variables,consumer_actors,HP_actors,Storage_actors,plant_actors,pipe_objects):
        for c, inp in zip(consumer_actors["Object"], input_variables["Consumer"]):
            c.Def_input_variables(inp)
        for hp, inp in zip(HP_actors["Object"], input_variables["HP_actors"]):
            hp.Def_input_variables(inp)

        for s, inp in zip(Storage_actors["Object"], input_variables["Storage_actors"]):
            s.Def_input_variables(inp)

        for m, inp in zip(plant_actors["Object"], input_variables["plants"]):
            m.Def_input_variables(inp)
        for p in pipe_objects["Object"]:
            p.Def_input_variables(input_variables["pipe_objects"])
        disp=input_variables.get("dispatch",{})      # {"order":["N18","N1","N12"],"ref":"N18","tol_w":1.0}
        self.dispatcher=MeritOrderDispatcher(plant_actors["Object"],order=disp.get("order"),
                                             ref=disp.get("ref"),tol_w=disp.get("tol_w",1.0),dead_band=disp.get("dead_band",1e-3))
        self.cut_rule=disp.get("cut_rule",proportional_cut)
        self.solve_cut=disp.get("cut_state",True)
        self.capture_states=disp.get("capture_states","infeasible")

    def _fixed_point(self,net,actors,T_guess,t,tol_k,tol_w,max_iter,relax,print_iteration):
        """UN pipeflow por iteración. Converge a la vez:
           (i)  las T de retorno de los actores a caudal impuesto (Gauss-Seidel con relax)
           (ii) el balance de potencia del despacho (planta marginal ajusta su Q_target)."""
        disp=self.dispatcher
        for it in range(1,max_iter+1):
            for actor in actors:
                actor.Boundary_condition(T_guess[actor.name],t)
            try: 
                pp.pipeflow(net, mode="bidirectional", tol_T=tol_k)
            except Exception: 
                print(net.circ_pump_mass[["mdot_flow_kg_per_s", "in_service"]])
                raise
            T_real={a.name: a.read_T_result(t) for a in actors}
            for a in actors:
                if T_real[a.name] is None or T_guess[a.name] is None:
                    print(a.name, type(a).__name__, "T_real:", T_real[a.name], "T_guess:", T_guess[a.name], "role:", getattr(a, "role", None))
            error_T=max((abs(T_real[a.name]-T_guess[a.name]) for a in actors if getattr(a,"iterates",True)),default=0.0)
            error_Q=disp.update_balance()
            if print_iteration:
                print(f"  hour {t} stage {disp.stage_k} iter {it}: dT={error_T:.2e} K | dQ_slack={error_Q:.2e} W")
            if error_T<tol_k and error_Q<tol_w:
                return T_real,it
            T_guess={a.name: relax*T_real[a.name]+(1-relax)*T_guess[a.name] for a in actors}
        raise ConvergenceError(f"Hour{t}: stage {disp.stage_k} sin convergencia tras {it} iteraciones (dT={error_T}, dQ={error_Q})")

    def solve_hour(self,net,consumer_actors,plant_actors,HP_actors,Storage_actors,pipe_objects,T_prev,t,tol_k=1E-6,tol_w=1.0,max_iter=1000,relax=0.5,print_iteration=False):
        for c in consumer_actors["Object"]:
            c.Boundary_condition(t)
        for p in pipe_objects["Object"]:
            p.Boundary_condition(t)
        disp=self.dispatcher
        disp.start_hour(t)                      # Q_max de cada planta (chaufferies) ANTES de iterar
        actors=plant_actors["Object"]+HP_actors["Object"]+Storage_actors["Object"]
        T_guess=dict(T_prev)
        beginning=time.perf_counter()
        n_iter_total=0
        last_ok=None                            # (k, T_real) de la ultima etapa que SI convergio
        converged=True
        note=""
        demand=float(np.nansum(net.heat_consumer.loc[net.heat_consumer["in_service"],"qext_w"]))
        for k in range(disp.n_stages):          # orden de mérito: una etapa = una planta más
            disp.stage_k=k
            disp.set_stage(k,demand)
            try:
                T_real,it=self._fixed_point(net,actors,T_guess,t,tol_k,tol_w,max_iter,relax,print_iteration)
            except (ConvergenceError,PipeflowNotConverged) as e:
                if last_ok is None:
                    raise                       # ni la primera etapa converge: no hay nada que rescatar
                # Rescate: la etapa k no converge, pero la etapa k-1 SI. Se re-resuelve y se conserva como
                # estado "verdadero" de la hora (la red cubre toda la demanda, las plantas se pasan de Q_max).
                converged=False
                note=f"etapa {k+1} no convergio ({type(e).__name__}); se conserva el estado de la etapa {k}"
                k_ok,T_ok=last_ok
                disp.stage_k=k_ok
                disp.set_stage(k_ok,demand)
                T_real,it=self._fixed_point(net,actors,T_ok,t,tol_k,tol_w,max_iter,relax,print_iteration)
                n_iter_total+=it
                disp.stage_ok()
                disp.feasible=False
                break
            n_iter_total+=it
            T_guess=T_real                      # arranque en caliente para la etapa siguiente
            last_ok=(k,T_real)
            if disp.stage_ok():
                break
        snap=disp.snapshot(t,demand_w=demand,converged=converged,note=note)
        if not disp.feasible:
            print(f"INFACTIBLE hora {t}: faltan {snap['deficit_w']/1e3:.3f} kW "
                  f"(la red pide {snap['Q_required_total']/1e3:.1f} kW, las plantas dan como maximo "
                  f"{snap['capacity_total']/1e3:.1f} kW, cobertura {100*snap['coverage']:.2f} %)"
                  +(f" [{note}]" if note else ""))
        self.dispatch_status[t]=dict(feasible=disp.feasible,deficit_w=snap["deficit_w"],coverage=snap["coverage"],
                                     converged=converged,stages=disp.stage_k+1,slack=disp.plants[disp.slack_i].name,note=note)
        if self.capture_states=="all" or (self.capture_states=="infeasible" and not disp.feasible):
            self._capture_hour(net,t,actors,consumer_actors,T_real,snap,demand,tol_k,tol_w,max_iter,relax)
        for a in actors:
            a.log_hour(t,T_real[a.name])
        return T_real,n_iter_total,time.perf_counter()-beginning

    def simulate_period_time(self,net,consumer_actors,plant_actors,HP_actors,Storage_actors,pipe_objects,n_hours,T_prev=None,tol_k=1E-6,tol_w=1.0,max_iter=1000,relax=0.5,print_iteration=False):
        self.time_per_step_hour = []
        self.dispatch_status = {}
        self.states = {}
        self.cut_notes = {}
        plants=plant_actors["Object"]
        actors=plants+HP_actors["Object"]+Storage_actors["Object"]
        T_prev=dict(T_prev or {})
        for a in actors:
            T_prev.setdefault(a.name,305.0)
        self.history = {a.name: {}  for a in actors}
        self.CO2={c.name: {} for c in plants}
        self.Cost={c.name: {} for c in plants}
        n_iter_log = []
        self.create_output_variables(n_consumers=len(consumer_actors["Object"]),n_HP_actors=len(HP_actors["Object"]),n_Storage_actors=len(Storage_actors["Object"]),n_plant_actors=len(plants),n_pipe_objects=len(pipe_objects["Object"]))
        self.organize_output_variables_objects(consumer_actors,HP_actors,Storage_actors,plant_actors,pipe_objects)
        for hour in range(n_hours):
            T_real, n_iter, duration = self.solve_hour(net,consumer_actors,plant_actors,HP_actors,Storage_actors,pipe_objects,T_prev,tol_k=tol_k,tol_w=tol_w,max_iter=max_iter,relax=relax,print_iteration=print_iteration,t=hour)
            self.time_per_step_hour.append(duration)
            for a in actors:
                self.history[a.name][hour] = T_real[a.name]
            for storage in Storage_actors["Object"]:
                storage.finalize_hour(t_real_k=T_real[storage.name],t=hour)
            for p in plants:
                p.Recover_information(hour)
                p.Consumption_Fuel()
                self.CO2[p.name][hour]=p.CO2_emissions()
                self.Cost[p.name][hour]=p.Operation_cost()
            n_iter_log.append(n_iter)
            T_prev = T_real
            self.attribute_output_variables(hour,consumer_actors,HP_actors,Storage_actors,plant_actors,pipe_objects)

    def generate_report(self, input_variables, Detailed_report=True,generate_data=True,write_detailed_report_name="simulation_report.txt",data_report_name="simulation_data_report.xlsx"):
        # This function will generate a report of the simulation in format 
        tables = self.transform_output_variables_to_dataframe()
        if generate_data:
            with pd.ExcelWriter(data_report_name, engine='openpyxl') as writer:
                for keys,values in tables.items():
                    values.to_excel(writer, sheet_name=keys)
        if Detailed_report:
            with open(write_detailed_report_name, "w", encoding="utf-8") as f:
                f.write(bold("Simulation Report")+"\n")
                f.write("==================================================\n\n")
                f.write('The simulation was performed using the following objects:\n')
                for objects_names in tables.keys():
                        f.write(f"  - {objects_names}\n")
                f.write("\n\n")
                f.write(bold("Their input variables were defined as follows:")+"\n\n")
                for key,values in input_variables.items():
                    if isinstance(values, dict):
                        if "pipe_objects" in key:
                            f.write("==================================================\n")
                            f.write("Pipe Objects Input Variables:\n")
                            f.write("==================================================\n")
                            for key2, value in values.items():
                                f.write(f"{key2}: \n{value}\n")
                                f.write("\n")
                        if "Slack_Central_production" in key:
                            f.write("==================================================\n")
                            f.write("Slack Central production_time_serie_Input Variables:\n")
                            f.write("==================================================\n\n")
                            for key2, value in values.items():
                                f.write(f"{key2}: \n {value}\n")
                                f.write("\n")
                    if isinstance(values, list):
                        if "Consumer" in key:
                            f.write("==================================================\n")
                            f.write("Consumer_time serie Input Variables:\n")
                            f.write("==================================================\n")              
                            for idx,obj in enumerate(values):
                                f.write("  -"+self.Consumer_actors["name"][idx] +":\n")
                                for key2, value in obj.items():
                                    f.write(f"{key2}: \n{value}\n")
                                    f.write("\n")
                        if "HP_actors" in key:
                            f.write("==================================================\n")
                            f.write("HP Actors_time_serie_Input Variables:\n")
                            f.write("==================================================\n")
                            if len(self.HP_actors["name"])==0:
                                f.write("There is not existence of Bidirectional heat pump actors\n") 
                            else:
                                for idx,obj in enumerate(values):
                                    f.write("  -"+self.HP_actors["name"][idx]+":\n")
                                    for key2, value in obj.items():
                                        f.write(f"{key2}: \n{value}\n")
                                        f.write("\n")
                        if "Storage_actors" in key:
                            f.write("==================================================\n")
                            f.write("Storage Actors Input Variables:\n")
                            f.write("==================================================\n")
                            if len(self.Storage_actors["name"])==0:
                                    f.write("There is not existence of storage actors\n") 
                            else:
                                for idx,obj in enumerate(values):
                                    f.write("  -"+ self.Storage_actors["name"][idx]+ ":\n")
                                    for key2, value in obj.items():
                                        f.write(f"{key2}: \n{value}\n")
                                        f.write("\n")
                        if "plants" in key:
                            f.write("==================================================\n")
                            f.write("Multi Plant Input Variables:\n")
                            f.write("==================================================\n")
                            if len(self.plant_actors["name"])==0:
                                f.write("There is not existence of multiplant actors\n\n") 
                            else:
                                for idx,obj in enumerate(values):
                                        f.write("  -"+ self.plant_actors["name"][idx]+":\n")
                                        for key2, value in obj.items(): 
                                            f.write(f"{key2}: \n{value}\n")
                                            f.write("\n")

                f.write("==================================================\n\n")
                f.write(bold('Initial state of tank energy storage devices:\n'))
                for i in range(len(self.Storage_characteristics["node"])):
                    f.write("  -"+ self.Storage_actors["name"][i]+ "\n")
                    f.write("T_hot_init="+ str(self.Storage_characteristics["T_hot_init"][i])+ "\n")
                    f.write("T_cold_init="+str(self.Storage_characteristics["T_cold_init"][i])+ "\n")
                    f.write("v_hot_fraction="+ str(self.Storage_characteristics["v_hot_fraction"][i])+ "\n")
                f.write("==================================================\n\n")   
                f.write(f"Total simulation time: {sum(self.time_per_step_hour):.2f} seconds\n")
                f.write(f"Average time per step: {np.mean(self.time_per_step_hour):.2f} seconds\n")
                f.write(f"Number of hours simulated: {len(self.time_per_step_hour)} hours\n")
                f.write("==================================================\n\n")
                f.write("Time per step (seconds):\n")
                f.write("==================================================\n")
                for hour, duration in enumerate(self.time_per_step_hour):
                    f.write(f"Hour {hour}: {duration:.5f} seconds\n")
                f.write("==================================================\n\n")
                f.write(bold("The results of the simulation for each object are:")+"\n\n")
                f.write("==================================================\n")
                for keys,values in tables.items():
                    f.write(bold(f"  - {keys}")+"\n\n")
                    f.write(values.to_string(float_format="{:,.5f}".format))
                    f.write("\n\n")
                    f.write("==================================================\n\n")

    def create_output_variables(self,n_consumers,n_HP_actors,n_Storage_actors,n_plant_actors,n_pipe_objects):
        self.output_variables = {"Consumers":[], "HP_actors":[], "Storage_actors":[], "plants":[], "pipe_objects":[]}
        for n in range(n_consumers):
            self.output_variables["Consumers"].append({"name": [], 'hour': [], "mdot_kg_s": [], "Pin": [], "Pout": []})
        for n in range(n_HP_actors):
            self.output_variables["HP_actors"].append({"name": [], 'hour': [], "Q_net": [], "W_elec": [], "COP": []})
        for n in range(n_Storage_actors):
            self.output_variables["Storage_actors"].append({"name": [], 'hour': [], "T_hot": [], "T_cold": [], "V_hot": [], "V_cold": [], "by_pass": [],"direction": []})
        for n in range(n_plant_actors):
            self.output_variables["plants"].append({"name": [], 'hour': [], "role": [], "Q_max": [], "mdot_kg_s":[],"Q_net": [],"Q_feasible":[],"Q_deficit":[],"Delta_k":[],"Pin":[],"Pout":[]})
        for n in range(n_pipe_objects):
            self.output_variables["pipe_objects"].append({"name": [], 'hour': [], "mdot_kg_s": [], "Tin": [], "Tout": [],"P_in": [],"P_out": []})
    def organize_output_variables_objects(self, consumer_actors, HP_actors, Storage_actors, plant_actors, pipe_objects):
        for i in range(len(consumer_actors["Object"])):
            c = consumer_actors["Object"][i]
            self.output_variables["Consumers"][i]["name"].append(c.name)
        for i in range(len(HP_actors["Object"])):
            hp = HP_actors["Object"][i]
            self.output_variables["HP_actors"][i]["name"].append(hp.name)
        for i in range(len(Storage_actors["Object"])):
            s = Storage_actors["Object"][i]
            self.output_variables["Storage_actors"][i]["name"].append(s.name)
        for i in range(len(plant_actors["Object"])):
            m = plant_actors["Object"][i]
            self.output_variables["plants"][i]["name"].append(m.name)
        for i in range(len(pipe_objects["Object"])):
            p = pipe_objects["Object"][i]
            self.output_variables["pipe_objects"][i]["name"].append(p.name)
    def attribute_output_variables(self,hour,consumer_actors,HP_actors,Storage_actors,plant_actors,pipe_objects):
        for i in range(len(consumer_actors["Object"])):
            c = consumer_actors["Object"][i]
            self.output_variables["Consumers"][i]["hour"].append(hour+1)
            self.output_variables["Consumers"][i]["mdot_kg_s"].append(self.net.res_heat_consumer.at[c.heat_consumer_id, "mdot_from_kg_per_s"])
            self.output_variables["Consumers"][i]["Pin"].append(self.net.res_heat_consumer.at[c.heat_consumer_id, "p_from_bar"])
            self.output_variables["Consumers"][i]["Pout"].append(self.net.res_heat_consumer.at[c.heat_consumer_id, "p_to_bar"])
        for i in range(len(HP_actors["Object"])):
            hp = HP_actors["Object"][i]
            self.output_variables["HP_actors"][i]["hour"].append(hour+1)
            if hp.mode[hour]=="COOLING_NET":
                self.output_variables["HP_actors"][i]["Q_net"].append(abs(self.net.res_heat_consumer.at[hp.HP_heat_consumer_id, "qext_w"]))
                self.output_variables["HP_actors"][i]["W_elec"].append(hp.Cooling_net_compressor.P.val)
                self.output_variables["HP_actors"][i]["COP"].append(abs(hp.Cooling_net_condenser.Q.val)/hp.Cooling_net_compressor.P.val)
            elif hp.mode[hour]=="HEATING_NET":
                self.output_variables["HP_actors"][i]["Q_net"].append(self.net.res_heat_consumer.at[hp.HP_heat_consumer_id, "qext_w"])
                self.output_variables["HP_actors"][i]["W_elec"].append(hp.Heating_net_compressor.P.val)
                self.output_variables["HP_actors"][i]["COP"].append(abs(hp.Heating_net_evaporator.Q.val)/hp.Heating_net_compressor.P.val)
        for i in range(len(Storage_actors["Object"])):
            s = Storage_actors["Object"][i]
            self.output_variables["Storage_actors"][i]["hour"].append(hour + 1)
            self.output_variables["Storage_actors"][i]["T_hot"].append(s.T_hot)
            self.output_variables["Storage_actors"][i]["T_cold"].append(s.T_cold)
            self.output_variables["Storage_actors"][i]["V_hot"].append(s.V_hot)
            self.output_variables["Storage_actors"][i]["V_cold"].append(s.V_cold)
            self.output_variables["Storage_actors"][i]["by_pass"].append(s.mdot_bypass)
            self.output_variables["Storage_actors"][i]["direction"].append(s.direction)
        for i in range(len(plant_actors["Object"])):
            m = plant_actors["Object"][i]
            self.output_variables["plants"][i]["hour"].append(hour+1)
            r = m.get_results()
            for key in ("role","Q_max","mdot_kg_s","Q_net","Q_feasible","Q_deficit","Delta_k","Pin","Pout"):
                self.output_variables["plants"][i][key].append(r[key])
        for i in range(len(pipe_objects["Object"])):
            p = pipe_objects["Object"][i]
            self.output_variables["pipe_objects"][i]["hour"].append(hour+ 1)
            self.output_variables["pipe_objects"][i]["mdot_kg_s"].append(self.net.res_pipe.at[p.pipe_id, "mdot_from_kg_per_s"])
            self.output_variables["pipe_objects"][i]["Tin"].append(self.net.res_pipe.at[p.pipe_id, "t_from_k"])
            self.output_variables["pipe_objects"][i]["Tout"].append(self.net.res_pipe.at[p.pipe_id, "t_to_k"])
            self.output_variables["pipe_objects"][i]["P_in"].append(self.net.res_pipe.at[p.pipe_id, "p_from_bar"])
            self.output_variables["pipe_objects"][i]["P_out"].append(self.net.res_pipe.at[p.pipe_id, "p_to_bar"])
    # ======================================================================
    # ESTADOS DE UNA HORA
    #   "required" : estado convergido con TODA la demanda servida (las plantas pueden pasarse de Q_max)
    #   "cut"      : estado recortado: plantas a Q_max y consumidores recortados segun self.cut_rule
    #                (en una hora factible no hay recorte: "cut" devuelve el mismo estado "required")
    # ======================================================================
    _RES={"pipe":"res_pipe","junction":"res_junction","consumer":"res_heat_consumer",
          "pump_pressure":"res_circ_pump_pressure","pump_mass":"res_circ_pump_mass"}
    _COMP={"pipe":"pipe","junction":"junction","consumer":"heat_consumer",
           "pump_pressure":"circ_pump_pressure","pump_mass":"circ_pump_mass"}

    def _raw_results(self,net):
        return {res:net[res].copy() for res in self._RES.values()}

    def _restore_results(self,net,raw):
        for res,df in raw.items():
            net[res]=df.copy()

    def _named_tables(self,net):
        out={}
        for key,res in self._RES.items():
            df=net[res].copy()
            comp=net[self._COMP[key]]
            if "name" in comp.columns:
                df.insert(0,"name",comp["name"].reindex(df.index))
            out[key]=df
        return out

    def _plants_table(self,disp):
        rows={}
        for p in disp.plants:
            q=p.Q_delivered()
            rows[p.name]=dict(role=p.role,Q_kW=q/1e3,Q_max_kW=p.Q_max_t/1e3,
                              Q_feasible_kW=min(q,p.Q_max_t)/1e3,Q_deficit_kW=max(q-p.Q_max_t,0.0)/1e3)
        return pd.DataFrame.from_dict(rows,orient="index")

    def _make_state(self,net,disp,demand,delivered,info):
        service=pd.DataFrame({"demand_kW":{k:v/1e3 for k,v in demand.items()},
                              "delivered_kW":{k:delivered[k]/1e3 for k in demand}})
        service["served_pct"]=100*service["delivered_kW"]/service["demand_kW"].where(service["demand_kW"]>0,np.nan)
        return dict(tables=self._named_tables(net),plants=self._plants_table(disp),service=service,info=info)

    def _solve_cut_state(self,net,actors,consumer_actors,T_start,t,demand_all,snap,tol_k,tol_w,max_iter,relax):
        """Plantas a Q_max (la de referencia es slack) y consumidores recortados segun self.cut_rule.
        El presupuesto de calor para los consumidores se corrige hasta que el slack entregue su Q_max."""
        disp=self.dispatcher
        cons=consumer_actors["Object"]
        demand={c.name:float(net.heat_consumer.at[c.heat_consumer_id,"qext_w"]) for c in cons}
        non_cuttable=demand_all-sum(demand.values())          # prosumers/HP: no se recortan
        losses=max(snap["Q_required_total"]-demand_all,0.0)   # perdidas del estado "required" (primera estimacion)
        budget=max(snap["capacity_total"]-losses-non_cuttable,0.0)
        disp.set_cut_stage()
        slack=disp.plants[disp.slack_i]
        band=disp.band(disp.slack_i)
        T_guess=dict(T_start)
        note=""
        for outer in range(1,21):
            raw=self.cut_rule(dict(demand),budget,t,net)
            delivered={k:min(max(float(raw.get(k,0.0)),0.0),demand[k]) for k in demand}
            for c in cons:
                q=delivered[c.name]
                net.heat_consumer.at[c.heat_consumer_id,"qext_w"]=max(q,1.0)
                net.heat_consumer.at[c.heat_consumer_id,"in_service"]=q>1.0     # sin calor: consumidor fuera de servicio
            T_cut,_=self._fixed_point(net,actors,T_guess,t,tol_k,tol_w,max_iter,relax,False)
            T_guess=T_cut
            err=slack.Q_delivered()-slack.Q_max_t
            if abs(err)<=band:
                break
            budget=max(budget-err,0.0)
        else:
            note=f"el slack no llego a Q_max tras 20 ajustes (error {err:.1f} W)"
        info=dict(budget_w=budget,outer_iterations=outer,served_pct=100*sum(delivered.values())/max(sum(demand.values()),1e-9),note=note)
        return self._make_state(net,disp,demand,delivered,info)

    def _capture_hour(self,net,t,actors,consumer_actors,T_real,snap,demand_all,tol_k,tol_w,max_iter,relax):
        disp=self.dispatcher
        cons=consumer_actors["Object"]
        ids=[c.heat_consumer_id for c in cons]
        dem={c.name:float(net.heat_consumer.at[c.heat_consumer_id,"qext_w"]) for c in cons}
        raw_req=self._raw_results(net)
        self.states[t]={"required":self._make_state(net,disp,dem,dem,
                        dict(feasible=disp.feasible,converged=snap["converged"],note=snap["note"]))}
        if disp.feasible or not self.solve_cut:
            return
        saved_roles=[(p.role,p.Q_target) for p in disp.plants]
        saved_disp=(disp.slack_i,disp.marginal_i,disp.balance_i)
        saved_q=net.heat_consumer.loc[ids,"qext_w"].copy()
        saved_srv=net.heat_consumer.loc[ids,"in_service"].copy()
        try:
            self.states[t]["cut"]=self._solve_cut_state(net,actors,consumer_actors,T_real,t,demand_all,snap,tol_k,tol_w,max_iter,relax)
        except (ConvergenceError,PipeflowNotConverged) as e:
            self.cut_notes[t]=f"el estado recortado no convergio ({type(e).__name__})"
            print(f"hora {t}: {self.cut_notes[t]}; solo esta disponible el estado 'required'")
        finally:      # dejar la red y las plantas EXACTAMENTE como estaban (estado "required")
            for p,(role,q) in zip(disp.plants,saved_roles):
                p.set_role(role,q)
            disp.slack_i,disp.marginal_i,disp.balance_i=saved_disp
            net.heat_consumer.loc[ids,"qext_w"]=saved_q
            net.heat_consumer.loc[ids,"in_service"]=saved_srv
            self._restore_results(net,raw_req)

    def results(self,hour,state="required",table="pipe"):
        """DataFrame de una hora y un estado.
           state : "required" (convergido, demanda completa) | "cut" (recortado)
           table : "pipe" | "junction" | "consumer" | "pump_pressure" | "pump_mass" | "plants" | "service"
        """
        st=self.states.get(hour)
        if not st:
            raise KeyError(f"Hora {hour}: no hay estados guardados (capture_states={self.capture_states!r}). "
                           "Usa capture_states='all' para guardar tambien las horas factibles.")
        if state not in st:
            if state=="cut" and st["required"]["info"]["feasible"]:
                obj=st["required"]                      # hora factible: no hubo recorte
            else:
                raise KeyError(f"Hora {hour}: no hay estado {state!r}. {self.cut_notes.get(hour,'')}")
        else:
            obj=st[state]
        if table in ("plants","service"):
            return obj[table].copy()
        return obj["tables"][table].copy()

    def compare(self,hour,table="pipe",columns=None):
        """required | cut | diff (cut - required) lado a lado."""
        a=self.results(hour,"required",table)
        b=self.results(hour,"cut",table)
        num=[c for c in a.columns if pd.api.types.is_numeric_dtype(a[c])]
        if columns:
            num=[c for c in columns if c in num]
        out=pd.concat({"required":a[num],"cut":b[num],"diff":b[num]-a[num]},axis=1)
        if "name" in a.columns:
            out.insert(0,("","name"),a["name"])
        return out

    def export_states(self,path="simulation_states.xlsx"):
        """Una hoja por estado y tabla (p. ej. required_pipe, cut_pipe), con todas las horas guardadas."""
        with pd.ExcelWriter(path,engine="openpyxl") as w:
            for state in ("required","cut"):
                for table in list(self._RES)+["plants","service"]:
                    frames=[]
                    for h,st in self.states.items():
                        if state in st:
                            df=(st[state][table] if table in ("plants","service") else st[state]["tables"][table]).copy()
                            df.insert(0,"hour",h)
                            frames.append(df)
                    if frames:
                        pd.concat(frames).to_excel(w,sheet_name=f"{state}_{table}"[:31])

    def deficit_report(self):
        """Una fila por hora: estado verdadero (lo que pide la red) vs estado posible (lo que dan las plantas)
        y cuanto calor falta. Columnas Deficit_<planta> = cuanto se pasa cada planta de su Q_max."""
        rows={}
        for t,snap in self.dispatcher.snapshots.items():
            row=dict(feasible=snap["feasible"],converged=snap["converged"],
                     demand_kW=snap["demand_w"]/1e3,
                     Q_required_kW=snap["Q_required_total"]/1e3,
                     Q_capacity_kW=snap["capacity_total"]/1e3,
                     Deficit_kW=snap["deficit_w"]/1e3,
                     coverage_pct=100*snap["coverage"],
                     stages=self.dispatch_status[t]["stages"],note=snap["note"])
            for name,r in snap["plants"].items():
                row["Required_kW_"+name]=r["Q_required"]/1e3
                row["Feasible_kW_"+name]=r["Q_feasible"]/1e3
                row["Deficit_kW_"+name]=r["Q_deficit"]/1e3
            rows[t]=row
        df=pd.DataFrame.from_dict(rows,orient="index")
        df.index.name="hour"
        return df

    def transform_output_variables_to_dataframe(self):
        tables={}
        for key, value in self.output_variables.items():
            for actor in value:
                name = actor["name"][0]
                datos={}
                for var, vals in actor.items():
                    if var not in ("name", "hour"):
                        datos[var]=vals
                df=pd.DataFrame(datos, index=actor["hour"])
                df.index.name = "hour"
                tables[name]= df
        if getattr(self,"dispatcher",None) is not None and self.dispatcher.snapshots:
            tables["Dispatch_deficit"]=self.deficit_report()
        return tables