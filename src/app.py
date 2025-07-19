# app.py

import asyncio
import streamlit as st
from langchain_core.messages import HumanMessage, AIMessage
from graph import run_agent

st.set_page_config(page_title="ChemChatBot", page_icon="🧪")
st.title("🧪 ChemChatBot: Chemistry Assistant")
st.write("Ask about molecular weight, logP, or formula. Example: *What is the molecular weight of aspirin?*")

# Session state initialization
if "messages" not in st.session_state:
    st.session_state.messages = []

# Button to clear the conversation history
if st.button("Clear Chat"):
    st.session_state.messages = []


# Display chat history
for msg in st.session_state.messages:
    if isinstance(msg, HumanMessage):
        with st.chat_message("user"):
            st.markdown(msg.content)
    elif isinstance(msg, AIMessage):
        with st.chat_message("assistant"):
            st.markdown(msg.content)

# Chat input and response
user_input = st.chat_input("Ask me something about a molecule...")
if user_input:
    user_msg = HumanMessage(content=user_input)
    st.session_state.messages.append(user_msg)

    with st.chat_message("user"):
        st.markdown(user_input)

    with st.chat_message("assistant"):
        with st.spinner("Analyzing..."):
            result = asyncio.run(run_agent(st.session_state.messages))
            new_messages = result["messages"][len(st.session_state.messages):]
            st.session_state.messages.extend(new_messages)
            for msg in new_messages:
                st.markdown(msg.content)
