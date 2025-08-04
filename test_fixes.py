#!/usr/bin/env python3
"""
Test script to verify the fixes for the Flask application issues.
"""

import requests
import json
import time

def test_csrf_token():
    """Test CSRF token endpoint."""
    try:
        response = requests.get("http://localhost:6513/api/csrf-token")
        if response.status_code == 200:
            data = response.json()
            print(f"✅ CSRF token endpoint working: {data.get('csrf_token', 'No token')[:20]}...")
            return True
        else:
            print(f"❌ CSRF token endpoint failed: {response.status_code}")
            return False
    except Exception as e:
        print(f"❌ CSRF token test failed: {e}")
        return False

def test_health_endpoint():
    """Test health endpoint."""
    try:
        response = requests.get("http://localhost:6513/api/healthz")
        if response.status_code == 200:
            data = response.json()
            print(f"✅ Health endpoint working: {data}")
            return True
        else:
            print(f"❌ Health endpoint failed: {response.status_code}")
            return False
    except Exception as e:
        print(f"❌ Health test failed: {e}")
        return False

def test_socketio_connection():
    """Test Socket.IO connection."""
    try:
        import socketio
        sio = socketio.Client()
        
        @sio.event
        def connect():
            print("✅ Socket.IO connection successful")
            sio.disconnect()
        
        @sio.event
        def disconnect():
            print("✅ Socket.IO disconnect successful")
        
        sio.connect("http://localhost:6513")
        time.sleep(2)  # Give time for connection
        return True
    except Exception as e:
        print(f"❌ Socket.IO test failed: {e}")
        return False

def main():
    """Run all tests."""
    print("🧪 Testing Flask application fixes...")
    print("=" * 50)
    
    tests = [
        ("Health Endpoint", test_health_endpoint),
        ("CSRF Token", test_csrf_token),
        ("Socket.IO Connection", test_socketio_connection),
    ]
    
    results = []
    for test_name, test_func in tests:
        print(f"\n🔍 Testing {test_name}...")
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"❌ {test_name} test crashed: {e}")
            results.append((test_name, False))
    
    print("\n" + "=" * 50)
    print("📊 Test Results:")
    for test_name, result in results:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"  {test_name}: {status}")
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    print(f"\n🎯 Overall: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests passed! The fixes are working correctly.")
    else:
        print("⚠️  Some tests failed. Check the application logs for more details.")

if __name__ == "__main__":
    main() 