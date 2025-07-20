import os
from typing import List, TypedDict,  Sequence, Annotated
from langchain_openai import ChatOpenAI
from langchain_core.messages import BaseMessage
from langgraph.graph import StateGraph, END
from langgraph.prebuilt import ToolNode,create_react_agent
from langchain_core.messages import HumanMessage
from langgraph.graph.message import add_messages
from prompts import SYSTEM_PROMPT
from dotenv import load_dotenv

from tools import (
    get_smiles_from_pubchem,
    get_molecular_properties,
    get_pdb_metadata,
    find_similar_molecules,
    web_search,
    search_patents,
    visualize_molecule
)

load_dotenv()

# ------------------
# Agent Setup
# ------------------

# Define the tools for the agent
tools = [
    get_smiles_from_pubchem,
    get_molecular_properties,
    get_pdb_metadata,
    find_similar_molecules,
    web_search,
    search_patents,
    visualize_molecule
]

# Set up the model
llm = ChatOpenAI(
    model="deepseek-ai/DeepSeek-R1",
    #model="Qwen/Qwen2.5-72B-Instruct",
    api_key=os.getenv("OPENAI_API_KEY"),
    streaming=True,
    temperature=0.0,
    base_url=os.getenv("OPENAI_API_BASE", "https://api.openai.com/v1"),
)


# ------------------
# Agent State and Graph
# ------------------

# Define the state to be passed around in LangGraph
class AgentState(TypedDict):
    messages: Annotated[List[BaseMessage], add_messages]

# Create the ReAct agent with both tools
agent_node = create_react_agent(
    model=llm,
    tools=tools,
    prompt=SYSTEM_PROMPT
)

# Build the state graph
workflow = StateGraph(AgentState)
workflow.add_node("agent", agent_node)
workflow.set_entry_point("agent")
workflow.set_finish_point("agent")
app_graph = workflow.compile()

# Limit how many messages are passed back to the agent
MAX_HISTORY = 8


# Wrapper for the Streamlit app
async def run_agent(messages: list[BaseMessage]) -> AgentState:
    """Invoke the graph with just the latest message history."""
    history = messages[-MAX_HISTORY:]
    return await app_graph.ainvoke({"messages": history})