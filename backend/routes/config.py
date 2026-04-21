"""
VaxElan 2.0 — Configuration endpoints.

Provides the frontend with available MHC methods, strategies,
and tool lists so the UI can render dynamic forms.
"""
from fastapi import APIRouter
from backend.models import MethodOption, StrategyInfo

router = APIRouter(prefix="/api/config", tags=["Configuration"])


# ── MHC-I Methods ────────────────────────────────────────────────

MHCI_METHODS = [
    MethodOption(key="a", name="Artificial Neural Network (ANN)"),
    MethodOption(key="b", name="Comblib Sidney 2008"),
    MethodOption(key="c", name="Consensus"),
    MethodOption(key="d", name="NetMHCcons"),
    MethodOption(key="e", name="NetMHCpan BA"),
    MethodOption(key="f", name="NetMHCpan EL"),
    MethodOption(key="g", name="NetMHCstabpan"),
    MethodOption(key="h", name="PickPocket"),
    MethodOption(key="i", name="SMM"),
    MethodOption(key="j", name="SMMPMBEC"),
]

@router.get("/mhci-methods", response_model=list[MethodOption])
async def get_mhci_methods():
    """Return available MHC-I prediction methods."""
    return MHCI_METHODS


# ── MHC-II Methods ───────────────────────────────────────────────

MHCII_METHODS = [
    MethodOption(key="1", name="NetMHCIIpan EL"),
    MethodOption(key="2", name="NetMHCIIpan BA"),
    MethodOption(key="3", name="Consensus3"),
    MethodOption(key="4", name="NN_align"),
    MethodOption(key="5", name="SMM_align"),
    MethodOption(key="6", name="Combinatorial library"),
    MethodOption(key="7", name="Sturniolo"),
    MethodOption(key="8", name="NetMHCIIpan EL 4.2"),
    MethodOption(key="9", name="NetMHCIIpan BA 4.2"),
    MethodOption(key="10", name="NetMHCIIpan EL 4.3"),
    MethodOption(key="11", name="NetMHCIIpan BA 4.3"),
]

@router.get("/mhcii-methods", response_model=list[MethodOption])
async def get_mhcii_methods():
    """Return available MHC-II prediction methods."""
    return MHCII_METHODS


# ── Strategies ───────────────────────────────────────────────────

STRATEGIES = [
    StrategyInfo(
        number=1,
        name="Epitope Prediction",
        description="Predict B-cell, MHC-I, MHC-II epitopes using IEDB tools, NetCTL, NetChop, and PSORTb.",
        available_tools=["mhc1", "mhc2", "netctl", "netchop", "bcell", "psortb"],
    ),
    StrategyInfo(
        number=2,
        name="Protein Prioritization",
        description="Assess antigenicity (IAPred), allergenicity (AlgPred), instability, molecular weight, and subcellular localization (WoLF PSORT).",
        available_tools=["iapred", "algpred", "instability", "molwt", "wolfpsort"],
    ),
    StrategyInfo(
        number=3,
        name="Virulence, Toxicity & Host Compatibility",
        description="Evaluate toxicity (ToxinPred), transmembrane topology (DeepTMHMM), solubility (NetSol), signal peptides (SignalP), and human homology.",
        available_tools=["toxinpred", "gutflora", "virulence", "deeptmhmm", "netsol", "signalp", "human"],
    ),
    StrategyInfo(
        number=4,
        name="Comprehensive Multi-Tool Analysis",
        description="Run strategies 1, 2, and 3 sequentially for a full analysis.",
        available_tools=[],
    ),
    StrategyInfo(
        number=5,
        name="Predicted Epitope Analysis",
        description="Pre-filter proteins, then predict and analyze epitopes with downstream tools (AlgPred, IAPred, ToxinPred3, PepMatch).",
        available_tools=["algpred", "iapred", "toxinpred3", "pepmatch", "clbtope"],
    ),
    StrategyInfo(
        number=6,
        name="Multi-Epitope Vaccine Assembly & 3D Modeling",
        description="Assemble multi-epitope vaccine constructs from predicted epitopes, run ESMFold 3D structural validation, and optional SASA analysis.",
        available_tools=["esmfold", "sasa"],
    ),
]

@router.get("/strategies", response_model=list[StrategyInfo])
async def get_strategies():
    """Return all 6 VaxElan strategies with descriptions and tool lists."""
    return STRATEGIES
