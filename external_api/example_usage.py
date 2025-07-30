#!/usr/bin/env python3
"""
Example usage of the ToxIndex ToxRAP API client
"""

import json
from datetime import datetime
from test_toxrap_api import query_toxrap

# Simple usage with extended timeout
prompt = "Is PFOA carcinogenic?"
result = query_toxrap(f"{prompt} Output as Json", timeout_minutes=30)

if result['success']:
    print(f"✅ Success! ({len(result['result'])} chars)")
    print(result['result'][:200] + "...")
    
    # Save result as JSON
    timestamp = datetime.now().isoformat().replace(':', '-')
    filename = f"toxrap_result_{timestamp}.json"
    
    # Try to parse as JSON if it's a string, else save as text
    try:
        if isinstance(result['result'], str):
            parsed_result = json.loads(result['result'])
        else:
            parsed_result = result['result']
        
        with open(filename, "w", encoding="utf-8") as f:
            json.dump(parsed_result, f, indent=2, ensure_ascii=False)
        print(f"📄 Saved JSON result to: {filename}")
        
    except json.JSONDecodeError:
        # If not valid JSON, save as text file
        filename = f"toxrap_result_{timestamp}.txt"
        with open(filename, "w", encoding="utf-8") as f:
            f.write(result['result'])
        print(f"📄 Saved text result to: {filename}")
        
else:
    print(f"❌ Error: {result['error']}") 