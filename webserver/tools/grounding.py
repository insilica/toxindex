from langchain_google_community import GoogleSearchAPIWrapper
import os
# Check for API key in environment
cse_api_key = os.environ.get("CSE_API_KEY")
google_cse_id = os.environ.get("GOOGLE_CSE_ID")
search_tool = GoogleSearchAPIWrapper(google_api_key=cse_api_key, google_cse_id=google_cse_id)

serp_text = search_tool.run('is gentamicin nephrotoxic?')
print(serp_text)

from google import genai
from google.genai import types

# Configure the client
client = genai.Client()

# Define the grounding tool
grounding_tool = types.Tool(
    google_search=types.GoogleSearch()
)

# Configure generation settings
config = types.GenerateContentConfig(
    tools=[grounding_tool]
)

# Make the request
response = client.models.generate_content(
    model="gemini-2.5-flash",
    contents="Who won the euro 2024?",
    config=config,
)

# Print the grounded response
print(response.text)