import aiohttp
import asyncio
import os
from typing import List, Optional
from langchain_core.tools import tool
from pydantic import BaseModel
from rdkit import Chem
from rdkit.Chem import Descriptors
from dotenv import load_dotenv

load_dotenv()

SERPAPI_API_KEY = os.getenv("SERPAPI_API_KEY")

### 1. PubChem Tool
class PubChemResponse(BaseModel):
    smiles: str
    molecular_formula: str

@tool
async def get_smiles_from_pubchem(name: str) -> PubChemResponse:
    """Retrieve SMILES and molecular formula for a given compound name using PubChem."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES,MolecularFormula/JSON"
    async with aiohttp.ClientSession() as session:
        async with session.get(url) as response:
            if response.status == 200:
                data = await response.json()
                props = data["PropertyTable"]["Properties"][0]
                return PubChemResponse(smiles=props["ConnectivitySMILES"], molecular_formula=props["MolecularFormula"])
            raise ValueError("Compound not found or PubChem API error")

### 2. RDKit Property Tool
class RDKitProperties(BaseModel):
    molecular_weight: float
    tpsa: float
    logp: float

@tool
async def calculate_rdkit_properties(smiles: str) -> RDKitProperties:
    """Calculate MW, TPSA, and logP from SMILES using RDKit."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    mw = Descriptors.MolWt(mol)
    tpsa = Descriptors.TPSA(mol)
    logp = Descriptors.MolLogP(mol)
    return RDKitProperties(molecular_weight=mw, tpsa=tpsa, logp=logp)

### 3. PDB Metadata Tool
class PDBMetadata(BaseModel):
    resolution: Optional[str]
    organism: Optional[str]
    expression_system: Optional[str]
    deposition_date: Optional[str]

@tool
async def get_pdb_metadata(pdb_id: str) -> PDBMetadata:
    """Query metadata for a protein from RCSB PDB."""
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    async with aiohttp.ClientSession() as session:
        async with session.get(url) as response:
            if response.status == 200:
                data = await response.json()
                struct = data.get("struct", {})
                src = data.get("rcsb_entry_container_identifiers", {})
                exp = data.get("rcsb_accession_info", {})
                return PDBMetadata(
                    resolution=data.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0],
                    organism=struct.get("title"),
                    expression_system=data.get("exptl", [{}])[0].get("expression_system"),
                    deposition_date=data.get("rcsb_accession_info", {}).get("initial_release_date")
                )
            raise ValueError("PDB ID not found or API error")

### 4. SmallWorld Similarity Tool
class SimilarMolecule(BaseModel):
    smiles: str
    name: Optional[str]
    similarity: float

@tool
async def search_similar_molecules(smiles: str, db: str = "chembl_31", dist: int = 4, top: int = 5) -> List[SimilarMolecule]:
    """Use SmallWorld API to find structurally similar molecules."""
    submit_url = "https://sw.docking.org/api/submit"
    query_url = "https://sw.docking.org/api/query"
    async with aiohttp.ClientSession() as session:
        async with session.get(submit_url, params={"db": db, "smiles": smiles, "dist": dist, "top": top}) as resp:
            if resp.status != 200:
                raise ValueError("SmallWorld submit error")
            hlid = (await resp.json())["hlid"]
        for _ in range(10):
            await asyncio.sleep(2)
            async with session.get(query_url, params={"hlid": hlid}) as qresp:
                if qresp.status == 200:
                    result = await qresp.json()
                    hits = result.get("hits", [])[:top]
                    return [SimilarMolecule(smiles=hit["smiles"], name=hit.get("name"), similarity=hit["dist"]) for hit in hits]
        raise ValueError("SmallWorld search timed out")

### 5. Web Search (SerpAPI)
@tool
async def search_chemistry_web(query: str) -> List[str]:
    """Perform a general web search for chemistry queries."""
    url = f"https://serpapi.com/search.json?q={query}&api_key={SERPAPI_API_KEY}&engine=google"
    async with aiohttp.ClientSession() as session:
        async with session.get(url) as response:
            if response.status == 200:
                data = await response.json()
                return [r.get("title") for r in data.get("organic_results", [])[:5]]
            raise ValueError("Web search failed")

### 6. Patent Search Tool
@tool
async def search_patents(target: str) -> List[str]:
    """Search patents related to a target (e.g., GLP-1, BACE) via SerpAPI."""
    url = f"https://serpapi.com/search.json?engine=google_patents&q={target}&api_key={SERPAPI_API_KEY}"
    async with aiohttp.ClientSession() as session:
        async with session.get(url) as response:
            if response.status == 200:
                data = await response.json()
                return [r.get("title", "No title") for r in data.get("organic_results", [])[:5]]
            raise ValueError("Patent search failed")
