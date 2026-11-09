"""Parche pandapipes 0.14.0: permite circ_pumps con in_service mixto (True/False) en la misma tabla.
Importar UNA vez antes de pipeflow (Structure_information ya lo hace)."""
import numpy as np
from pandapipes.component_models.abstract_models.circulation_pump import CirculationPump
from pandapipes.idx_branch import D, AREA, TOUTINIT
from pandapipes.idx_node import TINIT
from pandapipes.pf.pipeflow_setup import get_lookup

def _create_pit_branch_entries_fixed(cls, net, branch_pit):
    node_pit = net["_pit"]["node"]
    junction_idx_lookups = get_lookup(net, "node", "index")[cls.get_connected_node_type().table_name()]
    circ_pump_tbl = net[cls.table_name()]          # tabla COMPLETA (no filtrar por in_service)
    circ_pump_pit = super(CirculationPump, cls).create_pit_branch_entries(net, branch_pit)
    circ_pump_pit[:, D] = 0.1
    circ_pump_pit[:, AREA] = circ_pump_pit[:, D] ** 2 * np.pi / 4
    types = circ_pump_tbl.type.values
    mask_t = np.isin(types, ["pt", "t"])
    juncts = circ_pump_tbl[cls.from_to_node_cols()[1]].values
    to_nodes = junction_idx_lookups[juncts]
    circ_pump_pit[mask_t, TOUTINIT] = circ_pump_tbl.t_flow_k.values[mask_t]
    circ_pump_pit[~mask_t, TOUTINIT] = node_pit[to_nodes[~mask_t], TINIT]
    return circ_pump_pit

CirculationPump.create_pit_branch_entries = classmethod(_create_pit_branch_entries_fixed)
