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
from External_objects_functions import Bidirectional_W_to_WHeatPump,ThermoclineTwoLayer,ConsumerActor,Slack_Central_production,Central_production,Pipe_class
import time


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
    def Graph_From_Excel(self,excel_path,sheet=0,Node_Column=0,Edge_row=0): 
        """This function will create a graph representation from a incidenct matrix coming from an excel descrption"""
        Sheet_1=pd.read_excel(excel_path,sheet)
        Sheet_2=pd.read_excel(excel_path,sheet+1)
        Sheet_3=pd.read_excel(excel_path,sheet+2)
        Sheet_4=pd.read_excel(excel_path,sheet+3)
        Sheet_5=pd.read_excel(excel_path,sheet+4)
        Sheet_6=pd.read_excel(excel_path,sheet+5)
        Sheet_7=pd.read_excel(excel_path,sheet+6)
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
        #Central plant characteristics 
        self.central_plant_characteristicst={"node":[],"Fuel_type":[],"Nominal_efficiency":[],"Emission_factor":[],"cost_constant":[],'T_out_design':[],'P_out_design':[],'P_lift_design':[]}
        self.central_plant_characteristicst_excel=Sheet_5.iloc[:,Node_Column+1:]
        self.central_plant_characteristicst["node"]=self.central_plant_characteristicst_excel.columns.tolist()
        self.central_plant_characteristicst["Fuel_type"]=self.central_plant_characteristicst_excel.iloc[0,:].tolist()
        self.central_plant_characteristicst["Nominal_efficiency"]=self.central_plant_characteristicst_excel.iloc[1, :].tolist()
        self.central_plant_characteristicst["Emission_factor"]=self.central_plant_characteristicst_excel.iloc[2, :].tolist()
        self.central_plant_characteristicst["cost_constant"]=self.central_plant_characteristicst_excel.iloc[3, :].tolist()

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
        #Plant_support characteristics 
        self.plant_support_charateristics={"node":[],"Fuel_type":[],"Nominal_efficiency ":[],"Emission_factor":[],"cost_constant":[],'T_out_design':[]}
        self.plant_support_charateristics_excel=Sheet_7.iloc[:,Node_Column+1:]
        if self.plant_support_charateristics_excel.empty:
            pass
        else:
            self.plant_support_charateristics["node"]=self.plant_support_charateristics_excel.columns.tolist()
            self.plant_support_charateristics["Fuel_type"]=self.plant_support_charateristics_excel.iloc[0,:].tolist()
            self.plant_support_charateristics["Nominal_efficiency"]=self.plant_support_charateristics_excel.iloc[1,:].tolist()
            self.plant_support_charateristics["Emission_factor"]=self.plant_support_charateristics_excel.iloc[2, :].tolist()
            self.plant_support_charateristics["cost_constant"]=self.plant_support_charateristics_excel.iloc[3, :].tolist()

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

    def create_pandapipes_plus_initialization(self,fluid="water"):
        self.net=pp.create_empty_network(fluid=fluid)
        # Creation of pairs of nodes
        for i in range(len(self.Nodes)):
            if self.Type_Nodes[i]=="Storage":
                pp.create_junction(self.net,pn_bar=5,tfluid_k=353.15, name=self.Nodes[i] + '_highp')
                pp.create_junction(self.net,pn_bar=5,tfluid_k=353.15, name=self.Nodes[i] + '_inter_1')
                pp.create_junction(self.net,pn_bar=5,tfluid_k=353.15, name=self.Nodes[i] + '_inter_2')
                pp.create_junction(self.net,pn_bar=5,tfluid_k=353.15, name=self.Nodes[i] + '_lowp')

            elif self.Type_Nodes[i]=="Plant_support":
                pp.create_junction(self.net,pn_bar=5,tfluid_k=353.15, name=self.Nodes[i] + '_highp')
                pp.create_junction(self.net,pn_bar=5,tfluid_k=353.15, name=self.Nodes[i] + '_inter')
                pp.create_junction(self.net,pn_bar=5,tfluid_k=353.15, name=self.Nodes[i] + '_lowp')
            else:
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
                #Charging part
                pp.create_flow_control(
                    self.net, from_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_highp'].index[0], 
                    to_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_inter_1'].index[0],
                    controlled_mdot_kg_per_s=0.1,name="Flow_control_charge__input_Storage_"+self.Nodes[i])
                pp.create_flow_control(
                    self.net, from_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_inter_1'].index[0], 
                    to_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_highp'].index[0],
                    controlled_mdot_kg_per_s=0.1,name="Flow_control_bypass_charge_Storage_"+self.Nodes[i])
                pp.create_circ_pump_const_mass_flow(
                    self.net, return_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_inter_1'].index[0],
                    flow_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_inter_2'].index[0],
                    mdot_flow_kg_per_s=0.1, t_flow_k=300, p_flow_bar=5 ,name="Circ_pump_charge_Storage_"+ self.Nodes[i])
                pp.create_flow_control(
                    self.net, from_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_inter_2'].index[0],
                    to_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_lowp'].index[0],
                    controlled_mdot_kg_per_s=0.1,name="Flow_control_charge_output_Storage_"+ self.Nodes[i])
                #discharging part 
                pp.create_flow_control(
                    self.net, from_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_lowp'].index[0], 
                    to_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_inter_2'].index[0],
                    controlled_mdot_kg_per_s=0.1, name="Flow_control_discharge_input_Storage_"+ self.Nodes[i])
                pp.create_flow_control(
                    self.net, from_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_inter_2'].index[0], 
                    to_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_lowp'].index[0],
                    controlled_mdot_kg_per_s=0.1,name="Flow_control_bypass_discharge_Storage"+ self.Nodes[i])
                pp.create_circ_pump_const_mass_flow(
                    self.net, return_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_inter_2'].index[0], 
                    flow_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_inter_1'].index[0],
                    mdot_flow_kg_per_s=0.1, t_flow_k=300, p_flow_bar=5 ,name="Circ_pump_discharge_Storage_"+ self.Nodes[i])
                pp.create_flow_control(
                    self.net, from_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_inter_1'].index[0], 
                    to_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_highp'].index[0],
                    controlled_mdot_kg_per_s=0.1,name="Flow_control_discharge_output_Storage"+ self.Nodes[i])
            elif self.Type_Nodes[i]=="Plant_support":
                pp.create_circ_pump_const_mass_flow(
                    self.net, return_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_lowp'].index[0], 
                    flow_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_inter'].index[0],
                    mdot_flow_kg_per_s=0.5, t_flow_k=305, p_flow_bar=2.5,name="Plant_support_circ_pump_"+self.Nodes[i])
                pp.create_flow_control(self.net, from_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_inter'].index[0], 
                                    to_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_highp'].index[0],
                                    controlled_mdot_kg_per_s=0.1,name="Flow_control_support_plant"+ self.Nodes[i])

            elif self.Type_Nodes[i]=="Plant":
                print("im here ")
                pp.create_circ_pump_const_pressure(
                                    self.net,flow_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_highp'].index[0],
                                    return_junction=self.net.junction[self.net.junction["name"]==self.Nodes[i] + '_lowp'].index[0],
                                    p_flow_bar=5,plift_bar=5,t_flow_k=333 ,name='Central_plant')

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
        self.slack_plant={"Node":0,"name":0,"Object":0}
        self.HP_actors={"Node":[],"name":[],"Object":[]}
        self.Storage_actors={"Node":[],"name":[],"Object":[]}
        self.multi_plant_actors={"Node":[],"name":[],"Object":[]}
        for i in range(len(self.Nodes)):
            if self.Type_Nodes[i]=="Storage":
                Storage=ThermoclineTwoLayer(net=self.net,name="Storage_"+ self.Nodes[i],
                                            fc_charge=self.net.flow_control[self.net.flow_control["name"]=="Flow_control_charge__input_Storage_"+self.Nodes[i]].index[0],
                                            fc_bypass_charge=self.net.flow_control[self.net.flow_control["name"]=="Flow_control_bypass_charge_Storage_"+ self.Nodes[i]].index[0],
                                            cp_charge_storage=self.net.circ_pump_mass[self.net.circ_pump_mass["name"]=="Circ_pump_charge_Storage_"+ self.Nodes[i]].index[0],
                                            fc_charge_storage=self.net.flow_control[self.net.flow_control["name"]=="Flow_control_charge_output_Storage_"+ self.Nodes[i]].index[0],
                                            fc_discharge=self.net.flow_control[self.net.flow_control["name"]=="Flow_control_discharge_input_Storage_"+ self.Nodes[i]].index[0],
                                            fc_bypass_discharge=self.net.flow_control[self.net.flow_control["name"]=="Flow_control_bypass_discharge_Storage"+ self.Nodes[i]].index[0],
                                            cp_discharge_storage=self.net.circ_pump_mass[self.net.circ_pump_mass["name"]=="Circ_pump_discharge_Storage_"+ self.Nodes[i]].index[0],
                                            fc_discharge_storage=self.net.flow_control[self.net.flow_control["name"]=="Flow_control_discharge_output_Storage"+ self.Nodes[i]].index[0],
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

            elif self.Type_Nodes[i]=="Plant_support":
                Plant_support=Central_production(net=self.net,name='Central_Production'+self.Nodes[i],
                                      flow_control_id=self.net.flow_control[self.net.flow_control["name"]=="Flow_control_support_plant"+ self.Nodes[i]].index[0],
                                      Circ_mass_id=self.net.circ_pump_mass[self.net.circ_pump_mass["name"]=="Plant_support_circ_pump_"+self.Nodes[i]].index[0],
                                      Fuel_type=self.plant_support_charateristics["Fuel_type"][self.plant_support_charateristics["node"].index(self.Nodes[i])],
                                      Nominal_efficiency=self.plant_support_charateristics["Nominal_efficiency"][self.plant_support_charateristics["node"].index(self.Nodes[i])],
                                      Emission_factor=self.plant_support_charateristics["Emission_factor"][self.plant_support_charateristics["node"].index(self.Nodes[i])],
                                      cost_constant=self.plant_support_charateristics["cost_constant"][self.plant_support_charateristics["node"].index(self.Nodes[i])])
                Plant_support._post_init()
                self.multi_plant_actors["Node"].append(self.Nodes[i])
                self.multi_plant_actors["name"].append(Plant_support.name)
                self.multi_plant_actors["Object"].append(Plant_support)

            elif self.Type_Nodes[i]=="Plant":
                slack_plant=Slack_Central_production(net=self.net,
                                                     name='Central_plant',
                                                     Circ_pump_id=self.net.circ_pump_pressure[self.net.circ_pump_pressure["name"]=='Central_plant'].index[0],
                                                     Fuel_type=self.central_plant_characteristicst["Fuel_type"][self.central_plant_characteristicst["node"].index(self.Nodes[i])],
                                                     Nominal_efficiency=self.central_plant_characteristicst["Nominal_efficiency"][self.central_plant_characteristicst["node"].index(self.Nodes[i])],
                                                     Emission_factor=self.central_plant_characteristicst["Emission_factor"][self.central_plant_characteristicst["node"].index(self.Nodes[i])],
                                                     cost_constant=self.central_plant_characteristicst["cost_constant"][self.central_plant_characteristicst["node"].index(self.Nodes[i])])
                self.slack_plant["Node"]=self.Nodes[i]
                self.slack_plant["name"]=slack_plant.name
                self.slack_plant["Object"]=slack_plant

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
    def def_input_variables_simulation(self,input_variables,consumer_actors,slack_plant,HP_actors,Storage_actors,multi_plant_actors,pipe_objects):
        for c, inp in zip(consumer_actors["Object"], input_variables["Consumer"]):
            c.Def_input_variables(inp)

        slack_plant["Object"].Def_input_variables(input_variables["Slack_Central_production"])

        for hp, inp in zip(HP_actors["Object"], input_variables["HP_actors"]):
            hp.Def_input_variables(inp)

        for s, inp in zip(Storage_actors["Object"], input_variables["Storage_actors"]):
            s.Def_input_variables(inp)

        for m, inp in zip(multi_plant_actors["Object"], input_variables["multi_plant"]):
            m.Def_input_variables(inp)
        for p in pipe_objects["Object"]:
            p.Def_input_variables(input_variables["pipe_objects"])

    def solve_hour(self,net,consumer_actors,slack_plant,HP_actors,Storage_actors,multi_plant_actors,pipe_objects,T_prev,t,tol_k=1E-6,max_iter=1000,relax=0.5,print_iteration=False):
        for c in consumer_actors["Object"]:
            c.Boundary_condition(t)
        slack_plant["Object"].Boundary_condition(t)
        for p in pipe_objects["Object"]:
            p.Boundary_condition(t)
        actors=HP_actors["Object"]+Storage_actors["Object"]+multi_plant_actors["Object"]
        T_guess=T_prev
        # 1. Guardar el tiempo de inicio
        beginning = time.perf_counter()

        for it in range(1, max_iter + 1):
            for actor in actors:
                actor.Boundary_condition(T_guess[actor.name],t)
            pp.pipeflow(net,mode="bidirectional")
            T_real={a.name: a.read_T_result(t) for a in actors}
            error=max(abs(T_real[a.name] - T_guess[a.name]) for a in actors)

            if print_iteration:
                detail= ", ".join(f"{a.name}={T_real[a.name]:.2f}°C" for a in actors)
                print(f"  hour {t} iter {it}: error={error:.4f} K | {detail}")

            if error<tol_k:
                for a in actors:
                    a.log_hour(t, T_real[a.name])

                fin = time.perf_counter()
                duracion = fin - beginning
                return T_real, it,duracion
            T_guess={a.name: relax * T_real[a.name] + (1 - relax) * T_guess[a.name] for a in actors}
        fin = time.perf_counter()
        duracion = fin - beginning
        print(f"El código tardó {duracion:.5f} segundos.")
        raise ConvergenceError(f"Hour{t}:Not convergence, after N_iterations{it}, tolerance{error}")

    def simulate_period_time(self,net,consumer_actors,slack_plant,HP_actors,Storage_actors,multi_plant_actors,pipe_objects,n_hours,T_prev,tol_k=1E-6,max_iter=1000,relax=0.5,print_iteration=False):
        self.time_per_step_hour = []
        plants=[slack_plant["Object"]]+multi_plant_actors["Object"]
        actors=HP_actors["Object"]+Storage_actors["Object"]+multi_plant_actors["Object"]
        self.history = {a.name: {}  for a in actors}
        self.CO2={c.name: {} for c in plants}
        self.Cost={c.name: {} for c in plants}
        n_iter_log = []
        self.create_output_variables(n_consumers=len(consumer_actors["Object"]),n_HP_actors=len(HP_actors["Object"]),n_Storage_actors=len(Storage_actors["Object"]),n_multi_plant_actors=len(multi_plant_actors["Object"]),n_pipe_objects=len(pipe_objects["Object"]))
        self.organize_output_variables_objects(consumer_actors,slack_plant,HP_actors,Storage_actors,multi_plant_actors,pipe_objects)
        for hour in range(n_hours):
            T_real, n_iter, duration = self.solve_hour(net,consumer_actors,slack_plant,HP_actors,Storage_actors,multi_plant_actors,pipe_objects,T_prev,tol_k=tol_k,max_iter=max_iter,relax=relax,print_iteration=print_iteration,t=hour)
            self.time_per_step_hour.append(duration)
            for a in actors:
                self.history[a.name][hour] = T_real[a.name]
            if Storage_actors is not None:
                for storage in Storage_actors["Object"]:
                    storage.finalize_hour(t_real_k=T_real[storage.name],t=hour)
            if plants is not None:
                for p in plants:
                    p.Recover_information()
                    p.Consumption_Fuel()
                    CO2_emmision=p.CO2_emissions()
                    cost=p.Operation_cost()
                    self.CO2[p.name][hour]=CO2_emmision
                    self.Cost[p.name][hour]=cost
            n_iter_log.append(n_iter)
            T_prev = T_real
            self.attribute_output_variables(hour,consumer_actors,slack_plant,HP_actors,Storage_actors,multi_plant_actors,pipe_objects)

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
                        if "multi_plant" in key:
                            f.write("==================================================\n")
                            f.write("Multi Plant Input Variables:\n")
                            f.write("==================================================\n")
                            if len(self.multi_plant_actors["name"])==0:
                                f.write("There is not existence of multiplant actors\n\n") 
                            else:
                                for idx,obj in enumerate(values):
                                        f.write("  -"+ self.multi_plant_actors["name"][idx]+":\n")
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

    def create_output_variables(self,n_consumers,n_HP_actors,n_Storage_actors,n_multi_plant_actors,n_pipe_objects,n_slack_plant=1):
        self.output_variables = {"Consumers":[], "Slack_Central_production":[], "HP_actors":[], "Storage_actors":[], "multi_plant":[], "pipe_objects":[]}
        for n in range(n_consumers):
            self.output_variables["Consumers"].append({"name": [], 'hour': [], "mdot_kg_s": [], "Pin": [], "Pout": []})
        for n in range(n_slack_plant):
            self.output_variables["Slack_Central_production"].append({"name": [], 'hour': [], "mdot_kg_s":[],"Q_net": [],"Delta_k":[]})
        for n in range(n_HP_actors):
            self.output_variables["HP_actors"].append({"name": [], 'hour': [], "Q_net": [], "W_elec": [], "COP": []})
        for n in range(n_Storage_actors):
            self.output_variables["Storage_actors"].append({"name": [], 'hour': [], "T_hot": [], "T_cold": [], "V_hot": [], "V_cold": [], "by_pass": [],"direction": []})
        for n in range(n_multi_plant_actors):
            self.output_variables["multi_plant"].append({"name": [], 'hour': [], "mdot_kg_s":[],"Q_net": [],"Delta_k":[],"Pin":[],"Pout":[]})
        for n in range(n_pipe_objects):
            self.output_variables["pipe_objects"].append({"name": [], 'hour': [], "mdot_kg_s": [], "Tin": [], "Tout": [],"P_in": [],"P_out": []})
    def organize_output_variables_objects(self, consumer_actors, slack_plant, HP_actors, Storage_actors, multi_plant_actors, pipe_objects):
        for i in range(len(consumer_actors["Object"])):
            c = consumer_actors["Object"][i]
            self.output_variables["Consumers"][i]["name"].append(c.name)
        slk=slack_plant["Object"]
        self.output_variables["Slack_Central_production"][0]["name"].append(slk.name)
        for i in range(len(HP_actors["Object"])):
            hp = HP_actors["Object"][i]
            self.output_variables["HP_actors"][i]["name"].append(hp.name)
        for i in range(len(Storage_actors["Object"])):
            s = Storage_actors["Object"][i]
            self.output_variables["Storage_actors"][i]["name"].append(s.name)
        for i in range(len(multi_plant_actors["Object"])):
            m = multi_plant_actors["Object"][i]
            self.output_variables["multi_plant"][i]["name"].append(m.name)
        for i in range(len(pipe_objects["Object"])):
            p = pipe_objects["Object"][i]
            self.output_variables["pipe_objects"][i]["name"].append(p.name)
    def attribute_output_variables(self,hour,consumer_actors,slack_plant,HP_actors,Storage_actors,multi_plant_actors,pipe_objects):
        for i in range(len(consumer_actors["Object"])):
            c = consumer_actors["Object"][i]
            self.output_variables["Consumers"][i]["hour"].append(hour+1)
            self.output_variables["Consumers"][i]["mdot_kg_s"].append(self.net.res_heat_consumer.at[c.heat_consumer_id, "mdot_from_kg_per_s"])
            self.output_variables["Consumers"][i]["Pin"].append(self.net.res_heat_consumer.at[c.heat_consumer_id, "p_from_bar"])
            self.output_variables["Consumers"][i]["Pout"].append(self.net.res_heat_consumer.at[c.heat_consumer_id, "p_to_bar"])
        slk=slack_plant["Object"]
        self.output_variables["Slack_Central_production"][0]["hour"].append(hour+1)
        self.output_variables["Slack_Central_production"][0]["mdot_kg_s"].append(self.net.res_circ_pump_pressure.at[slk.Circ_pump_id, "mdot_from_kg_per_s"])
        self.output_variables["Slack_Central_production"][0]["Q_net"].append(self.net.res_circ_pump_pressure.at[slk.Circ_pump_id, "qext_w"])
        self.output_variables["Slack_Central_production"][0]["Delta_k"].append(self.net.res_circ_pump_pressure.at[slk.Circ_pump_id, "deltat_k"])
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
        for i in range(len(multi_plant_actors["Object"])):
            m = multi_plant_actors["Object"][i]
            self.output_variables["multi_plant"][i]["hour"].append(hour+1)
            self.output_variables["multi_plant"][i]["mdot_kg_s"].append(self.net.res_flow_control.at[m.flow_control_id, "mdot_from_kg_per_s"])
            self.output_variables["multi_plant"][i]["Q_net"].append(self.net.res_circ_pump_mass.at[m.Circ_mass_id, "qext_w"])
            self.output_variables["multi_plant"][i]["Delta_k"].append(self.net.res_circ_pump_mass.at[m.Circ_mass_id, "deltat_k"])
            self.output_variables["multi_plant"][i]["Pin"].append(self.net.res_circ_pump_mass.at[m.Circ_mass_id, "p_from_bar"])
            self.output_variables["multi_plant"][i]["Pout"].append(self.net.res_flow_control.at[m.flow_control_id, "p_to_bar"])
        for i in range(len(pipe_objects["Object"])):
            p = pipe_objects["Object"][i]
            self.output_variables["pipe_objects"][i]["hour"].append(hour+ 1)
            self.output_variables["pipe_objects"][i]["mdot_kg_s"].append(self.net.res_pipe.at[p.pipe_id, "mdot_from_kg_per_s"])
            self.output_variables["pipe_objects"][i]["Tin"].append(self.net.res_pipe.at[p.pipe_id, "t_from_k"])
            self.output_variables["pipe_objects"][i]["Tout"].append(self.net.res_pipe.at[p.pipe_id, "t_to_k"])
            self.output_variables["pipe_objects"][i]["P_in"].append(self.net.res_pipe.at[p.pipe_id, "p_from_bar"])
            self.output_variables["pipe_objects"][i]["P_out"].append(self.net.res_pipe.at[p.pipe_id, "p_to_bar"])
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
        return tables

