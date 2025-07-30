#!/usr/bin/env python3
"""
ToxIndex ToxRAP API Client
Simple function to get responses from the toxrap task API

Usage:
    from test_toxrap_api import query_toxrap
    result = query_toxrap("Is Gentamicin nephrotoxic?")
"""

import requests
import time
from urllib.parse import urlencode

# Configuration
BASE_URL = "https://toxindex.com"
TEST_EMAIL = "test@test.com"  # Update with real credentials
TEST_PASSWORD = "test"         # Update with real password

def query_toxrap(prompt, timeout_minutes=20):
    """
    Query the toxrap API with a prompt and return the response.
    
    Args:
        prompt (str): The query/prompt to send to toxrap
        timeout_minutes (int): Maximum time to wait for response (default: 20)
    
    Returns:
        dict: {
            'success': bool,
            'result': str (if successful),
            'error': str (if failed)
        }
    """
    print(f"🚀 Querying ToxRAP: {prompt}")
    
    # Create session and login
    session = requests.Session()
    form_data = urlencode({'email': TEST_EMAIL, 'password': TEST_PASSWORD})
    
    try:
        login_response = session.post(f"{BASE_URL}/api/auth/login", 
            headers={'Content-Type': 'application/x-www-form-urlencoded'},
            data=form_data
        )
        
        if login_response.status_code != 200 or not login_response.json().get('success'):
            return {'success': False, 'error': 'Login failed'}
        
        # Create task
        task_response = session.post(f"{BASE_URL}/api/tasks", json={
            "message": prompt,
            "workflow": 1,
            "environment_id": None
        })
        
        if task_response.status_code != 200:
            return {'success': False, 'error': 'Task creation failed'}
        
        task_id = task_response.json()["task_id"]
        print(f"✅ Task created: {task_id}")
        
        # Poll for completion
        start_time = time.time()
        timeout_seconds = timeout_minutes * 60
        
        while time.time() - start_time < timeout_seconds:
            status_response = session.get(f"{BASE_URL}/api/tasks/{task_id}")
            if status_response.status_code != 200:
                return {'success': False, 'error': 'Failed to get task status'}
            
            status = status_response.json().get("status", "unknown")
            elapsed = time.time() - start_time
            
            print(f"   Status: {status} ({elapsed:.1f}s)")
            
            if status == "done":
                break
            elif status == "error":
                error_msg = status_response.json().get('error', 'Unknown error')
                return {'success': False, 'error': f'Task failed: {error_msg}'}
            elif status in ["processing", "starting", "running agent", "checking cache", 
                          "agent complete", "sending message", "uploading file to GCS", "file uploaded"]:
                time.sleep(30)  # Wait 30 seconds before next poll
            else:
                time.sleep(30)
        else:
            return {'success': False, 'error': f'Timeout after {timeout_minutes} minutes'}
        
        # Get results
        messages_response = session.get(f"{BASE_URL}/api/tasks/{task_id}/messages")
        if messages_response.status_code != 200:
            return {'success': False, 'error': 'Failed to get messages'}
        
        messages = messages_response.json().get("messages", [])
        assistant_messages = [m for m in messages if m.get("role") == "assistant"]
        
        if not assistant_messages:
            return {'success': False, 'error': 'No assistant messages found'}
        
        result = assistant_messages[-1].get("content", "")
        print(f"✅ Response received ({len(result)} characters)")
        
        return {'success': True, 'result': result}
        
    except Exception as e:
        return {'success': False, 'error': f'Exception: {str(e)}'}

def quick_test():
    """Quick test function"""
    result = query_toxrap("Is Gentamicin nephrotoxic?")
    
    if result['success']:
        print("\n📄 Response:")
        print("-" * 40)
        print(result['result'][:500] + "..." if len(result['result']) > 500 else result['result'])
        print("-" * 40)
    else:
        print(f"\n❌ Error: {result['error']}")

if __name__ == "__main__":
    # Example usage
    print("🧪 ToxIndex ToxRAP API Test")
    print("=" * 50)
    
    # You can either run a quick test
    quick_test()
    
    # Or use it interactively
    print("\n💡 To use this in your code:")
    print("from test_toxrap_api import query_toxrap")
    print("result = query_toxrap('Your prompt here')")
    print("if result['success']:")
    print("    print(result['result'])")
    print("else:")
    print("    print(f'Error: {result[\"error\"]}')") 