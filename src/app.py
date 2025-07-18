# app.py

import streamlit as st
from langchain_core.messages import HumanMessage, AIMessage, SystemMessage
from graph import run_agent
from prompts import SYSTEM_PROMPT

st.set_page_config(page_title="ChemBot", page_icon="🧪")
st.title("🧪 ChemBot: Chemistry Assistant")
st.write("Example: *What is the molecular weight of aspirin?*")






import asyncio




if "chat_history" not in st.session_state:
    st.session_state.chat_history = []

user_input = st.chat_input("Ask me about a molecule (e.g., aspirin)")
if user_input:
    st.session_state.chat_history.append(HumanMessage(content=user_input))

if st.session_state.chat_history:
    for msg in st.session_state.chat_history:
        if isinstance(msg, HumanMessage):
            st.chat_message("user").write(msg.content)
        elif isinstance(msg, AIMessage):
            st.chat_message("assistant").write(msg.content)

    async def process():
        result = await run_agent(st.session_state.chat_history)
        st.session_state.chat_history.extend(result["messages"][len(st.session_state.chat_history):])

    asyncio.run(process())
