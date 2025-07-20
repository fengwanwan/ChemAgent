import asyncio
import streamlit as st
from langchain_core.messages import HumanMessage, AIMessage
from graph import run_agent

# Set Streamlit page configuration
st.set_page_config(page_title="ChemChatBot", page_icon="🧪")
st.title("🧪 ChemChatBot: Chemistry Assistant")
st.write("Ask about molecular weight, logP, or formula. Example: *What is the molecular weight of aspirin?*")

# --- Initialize session state ---
if "messages" not in st.session_state:
    st.session_state.messages = []
if "preserve_history" not in st.session_state:
    st.session_state.preserve_history = True  # Default to preserving history

# --- Sidebar toggle for preserving chat history ---
st.sidebar.title("Settings ⚙️")
st.session_state.preserve_history = st.sidebar.checkbox("Preserve chat history", value=True)

# --- Clear chat button ---
if st.sidebar.button("Clear Chat"):
    st.session_state.messages = []

# --- Display history if the user opted to preserve it ---
if st.session_state.preserve_history:
    for msg in st.session_state.messages:
        if isinstance(msg, HumanMessage):
            with st.chat_message("user"):
                st.markdown(msg.content)
        elif isinstance(msg, AIMessage):
            with st.chat_message("assistant"):
                st.markdown(msg.content)

# --- User input handling ---
user_input = st.chat_input("Ask me something about a molecule...")
if user_input:
    user_msg = HumanMessage(content=user_input)

    # Display user message (even if not preserving history)
    if st.session_state.preserve_history:
        st.session_state.messages.append(user_msg)
        with st.chat_message("user"):
            st.markdown(user_input)
    else:
        with st.chat_message("user"):
            st.markdown(user_input)

    with st.chat_message("assistant"):
        with st.spinner("Analyzing..."):
            # Build the input message list based on whether history is preserved
            input_messages = (
                st.session_state.messages + [user_msg]
                if st.session_state.preserve_history
                else [user_msg]
            )

            # Run the agent asynchronously
            result = asyncio.run(run_agent(input_messages))

            # Extract only the new assistant messages
            new_messages = result["messages"][len(input_messages):]

            # Filter for the last AIMessage (in case multiple assistant messages are returned)
            ai_messages = [msg for msg in new_messages if isinstance(msg, AIMessage)]
            if ai_messages:
                last_msg = ai_messages[-1]
                st.markdown(last_msg.content)

                # Save assistant response only if preserving history
                if st.session_state.preserve_history:
                    st.session_state.messages.extend(new_messages)