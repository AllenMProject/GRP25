# scripts/getenv.py
from dotenv import load_dotenv
import os
from pathlib import Path

# look for .env in the project root (one level up from scripts/)
dotenv_path = Path(__file__).parent.parent / '.env'
load_dotenv(dotenv_path=dotenv_path)

MAPI_KEY = os.getenv("MAPI_KEY")
