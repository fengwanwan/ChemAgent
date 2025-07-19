import os
import json
import asyncio
import requests
import aiohttp
import time
from typing import Optional

from langchain.tools import tool
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, AllChem, DataStructs
from serpapi import GoogleSearch
from dotenv import load_dotenv
from PIL import Image

load_dotenv()

db = 'REAL-Database-22Q1'
serpapi_key = os.getenv("SERPAPI_KEY")


@tool
async def get_smiles_from_pubchem(name: str) -> str:
    """Looks up a chemical compound by name on PubChem and returns its SMILES string and molecular formula."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES,MolecularFormula/JSON"
    print(f"[PubChem][GET] {url}")
    try:
        async with aiohttp.ClientSession() as session:
            await asyncio.sleep(1)  # Respect API rate limits
            async with session.get(url) as response:
                print(f"[PubChem][Status] {response.status}")
                if response.status == 200:
                    data = await response.json()
                    props = data["PropertyTable"]["Properties"][0]
                    smiles = props["ConnectivitySMILES"]
                    formula = props["MolecularFormula"]
                    return f"{name.capitalize()} -> SMILES: {smiles}, Formula: {formula}"
                else:
                    return f"Error: Could not find {name} in PubChem. Status: {response.status}"
    except Exception as e:
        return f"An error occurred while querying PubChem: {e}"


@tool
def get_molecular_properties(smiles: str) -> str:
    """
    Calculates molecular properties for a given SMILES string.
    
    Args:
        smiles (str): The SMILES string of the molecule.
        
    Returns:
        str: A JSON string of the properties (Molecular Weight, LogP, TPSA).
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "Error: Invalid SMILES string provided."
        properties = {
            "MolecularWeight": Descriptors.MolWt(mol),
            "LogP": Descriptors.MolLogP(mol),
            "TPSA": Descriptors.TPSA(mol)
        }
        return json.dumps(properties, indent=2)
    except Exception as e:
        return f"An error occurred during property calculation: {e}"


@tool
def get_pdb_metadata(pdb_id: str) -> str:
    """
    Fetches metadata for a given PDB ID from the RCSB PDB database.
    
    Args:
        pdb_id (str): The PDB ID (e.g., 1HVR).
        
    Returns:
        str: A JSON string of the PDB metadata.
    """
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.lower()}"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            info = {
                "Title": data.get("struct", {}).get("title", "N/A"),
                "Organism": data.get("rcsb_entry_container_identifiers", {}).get("entity_ids", ["N/A"]),
                "Resolution": data.get("rcsb_entry_info", {}).get("resolution_combined", ["N/A"]),
                "Deposition Date": data.get("rcsb_accession_info", {}).get("deposit_date", "N/A"),
                "Release Date": data.get("rcsb_accession_info", {}).get("initial_release_date", "N/A")
            }
            return json.dumps(info, indent=2)
        else:
            return f"Error: {response.status_code} – invalid PDB ID or not found."
    except Exception as e:
        return f"An error occurred while fetching PDB data: {e}"


async def _smallworld_submit_search(smiles: str, db: str = db, dist: int = 4, top: int = 5):
    url = "https://sw.docking.org/search/submit"
    params = {"smi": [smiles], "db": db, "dist": str(dist), "top": str(top)}
    async with aiohttp.ClientSession() as session:
        await asyncio.sleep(1) # Respect API rate limits
        async with session.get(url, params=params) as resp:
            if resp.status != 200:
                raise Exception(f"Submit error {resp.status}: {await resp.text()}")
            text = await resp.text()
            for line in text.splitlines():
                if line.startswith("data:"):
                    data = json.loads(line[5:])
                    if "hlid" in data and data["status"] == "END":
                        return data["hlid"]
    raise Exception("HLID not found in SmallWorld response.")


async def _smallworld_view_results(hlid: int, state_smiles: str):
    url = "https://sw.docking.org/search/view"
    params = {"hlid": hlid, "fmt": "tsv", "scores": "graph"}
    async with aiohttp.ClientSession() as session:
        await asyncio.sleep(1) # Respect API rate limits
        async with session.get(url, params=params) as resp:
            if resp.status != 200:
                raise Exception(f"View error {resp.status}: {await resp.text()}")
            lines = (await resp.text()).strip().splitlines()
            results = []
            for line in lines:
                if line.strip().lower().startswith("smiles"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue
                smiles_and_id = parts[0]
                try:
                    smiles, zid = smiles_and_id.rsplit(" ", 1)
                except ValueError:
                    smiles, zid = smiles_and_id, "N/A"
                try:
                    query_mol = Chem.MolFromSmiles(state_smiles)
                    result_mol = Chem.MolFromSmiles(smiles)
                    query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=2048)
                    result_fp = AllChem.GetMorganFingerprintAsBitVect(result_mol, 2, nBits=2048)
                    similarity = DataStructs.TanimotoSimilarity(query_fp, result_fp)
                except Exception:
                    similarity = 0.0
                results.append({"id": zid, "smiles": smiles, "similarity": similarity})
            results.sort(key=lambda x: x["similarity"], reverse=True)
            return results[:5]


@tool
async def find_similar_molecules(smiles: str) -> str:
    """Finds molecules similar to the given SMILES string using the SmallWorld API."""
    try:
        hlid = await _smallworld_submit_search(smiles)
        results = await _smallworld_view_results(hlid, smiles)
        if not results:
            return "No similar molecules found."
        
        top_lines = [
            f"{i+1}. SMILES: {item.get('smiles', 'N/A')}, ID: {item.get('id', 'N/A')}, Similarity: {item.get('similarity', 0.0):.4f}"
            for i, item in enumerate(results)
        ]
        return "Top similar molecules:\n" + "\n".join(top_lines)
    except Exception as e:
        return f"SmallWorld search failed: {e}"


@tool
def web_search(query: str) -> str:
    """Performs a web search using SerpAPI for general science or chemistry-related queries."""
    
    if not serpapi_key:
        return "Error: Missing SERPAPI_KEY."
    params = {"q": query, "api_key": serpapi_key, "engine": "google", "num": 5}
    try:
        search = GoogleSearch(params)
        result = search.get_dict()
        organic = result.get("organic_results", [])[:5]
        if not organic:
            return "No web results found."
        
        lines = [
            f"{i+1}. {item.get('title', 'N/A')}\n{item.get('snippet', 'N/A')}\n{item.get('link', '')}"
            for i, item in enumerate(organic)
        ]
        return "Top web search results:\n" + "\n\n".join(lines)
    except Exception as e:
        return f"Web search error: {e}"


def _format_patent_results(results: dict) -> str:
    lines = []
    for i, p in enumerate(results.get("organic_results", []), 1):
        lines.append(f"{i}. {p.get('title', 'N/A')}")
        lines.append(f"  Inventor: {p.get('inventor', 'N/A')}")
        lines.append(f"  Publication Date: {p.get('publication_date', 'N/A')}")
        lines.append(f"  Link: {p.get('patent_link', 'N/A')}")
    return "\n".join(lines)


@tool
def search_patents(query: str) -> str:
    """Searches for patents related to a chemical compound or target using Google Patents via SerpAPI."""
    if not serpapi_key:
        return "Error: SERPAPI_KEY not set."
    
    params = {"engine": "google_patents", "q": query, "api_key": serpapi_key, "num": 5}
    try:
        response = requests.get("https://serpapi.com/search.json", params=params)
        response.raise_for_status()
        data = response.json()
        return _format_patent_results(data)
    except Exception as e:
        return f"Patent search error: {e}"


@tool
def visualize_molecule(smiles: str) -> str:
    """Generates a 2D image of a molecule from a SMILES string and saves it to 'output.png'."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            img = Draw.MolToImage(mol, size=(300, 300))
            output_path = os.path.join(os.getcwd(), "output.png")
            img.save(output_path)
            return "Molecule image saved to output.png. It will be displayed in the chat."
        else:
            return "Failed to generate molecule from SMILES."
    except Exception as e:
        return f"An error occurred during visualization: {e}"
