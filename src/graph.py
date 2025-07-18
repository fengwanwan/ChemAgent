import os
from typing import TypedDict, Annotated, Sequence
from langgraph.graph import StateGraph
from langchain_openai import ChatOpenAI
from langchain_core.messages import BaseMessage, HumanMessage, AIMessage
from langchain_core.runnables.history import RunnableWithMessageHistory
from langgraph.prebuilt import ToolNode,create_react_agent
from prompts import SYSTEM_PROMPT



from tools import (
    get_smiles_from_pubchem,
    calculate_rdkit_properties,
    get_pdb_metadata,
    search_similar_molecules,
    search_chemistry_web,
    search_patents
)


from dotenv import load_dotenv

# Load environment variables
load_dotenv()

from prompts import SYSTEM_PROMPT

# 1. Define the state for our agent
class AgentState(TypedDict):
    messages: Annotated[Sequence[BaseMessage], lambda x, y: x + y]

# Define the list of tools
tools = [
    get_smiles_from_pubchem,
    calculate_rdkit_properties,
    get_pdb_metadata,
    search_similar_molecules,
    search_chemistry_web,
    search_patents
]


llm = ChatOpenAI(
    model="deepseek-ai/DeepSeek-R1",
    #model="Qwen/Qwen2.5-72B-Instruct",
    api_key=os.getenv("OPENAI_API_KEY"),
    base_url=os.getenv("OPENAI_API_BASE", "https://api.openai.com/v1"),
)

graph = create_react_agent(
    model = llm, 
    tools = tools, 
    prompt = SYSTEM_PROMPT)

# Compile it into an executable graph
app_graph = StateGraph(AgentState)
app_graph.add_node("agent", graph)
app_graph.set_entry_point("agent")
app_graph.set_finish_point("agent")
compiled_graph = app_graph.compile()

# Export a run function
async def run_agent(messages: Sequence[BaseMessage]) -> AgentState:
    return await compiled_graph.ainvoke({"messages": messages})