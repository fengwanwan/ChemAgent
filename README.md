
# 🧪 Chemistry Chatbot Agent

A LangGraph + Streamlit powered intelligent chemistry assistant that helps researchers with:

- 🧬 Molecular lookup (PubChem)
- 🔬 Molecular property calculations (RDKit)
- 🧫 Protein metadata (PDB)
- 🧪 Chemical similarity search (SmallWorld)
- 🌐 Chemistry-related web search (SerpAPI)
- 📄 Patent lookup for drug targets (Google Patents via SerpAPI)

## 🚀 Features

- 🔧 Async tools with robust error handling
- 🧠 Multi-turn memory (last 4 rounds)
- 📡 Tool routing via LangGraph agent
- 💬 Streamlit Web UI

## 📁 Project Structure

```
.
├── requirements.txt          # Python dependencies
├── src/
│   ├── app.py                # Streamlit UI
│   ├── graph.py              # LangGraph agent & memory
│   ├── tools.py              # All chemistry tools (async)
│   ├── prompts.py            # System prompt template
│   └── .env                  # API key template
            
```

## 🛠️ Setup

### 1. Clone the repo

```bash
git clone https://github.com/fengwanwan/ChemAgent.git
cd ChemAgent
```

### 2. Install dependencies

We recommend using a virtual environment:

```bash
conda create -n agent python=3.11
conda activate agent 
pip install -r requirements.txt
```

### 3. Setup environment variables

Copy `.env.example` to `.env` and add your API keys:

```bash
OPENAI_API_KEY=your-openai-key
SERPAPI_API_KEY=your-serpapi-key
```

### 4. Run the app

```bash
streamlit run src/app.py
```

## 🧪 Example Queries

- What is the molecular formula of aspirin?
- What's its logP and molecular weight?
- Show similar molecules to ibuprofen
- What is aspirin used for?
- Search patents related to GLP-1
- What protein data is available for PDB ID 1D66?

## 🧠 Learnings

This project showcases:
- Tool-augmented LLM reasoning
- Async multi-source API integration
- Memory + conversational context with LangGraph

---

Built for the **AI Scientist Take-Home Project**.
