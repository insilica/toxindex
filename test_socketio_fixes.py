#!/usr/bin/env python3
"""
Test script to verify the Socket.IO fixes.
"""

import requests
import json
import time
import socketio
import threading

def test_health_endpoint():
    """Test health endpoint."""
    try:
        response = requests.get("http://localhost:6513/api/healthz")
        if response.status_code == 200:
            print("✅ Health endpoint working")
            return True
        else:
            print(f"❌ Health endpoint failed: {response.status_code}")
            return False
    except Exception as e:
        print(f"❌ Health endpoint test failed: {e}")
        return False

def test_socketio_connection():
    """Test Socket.IO connection."""
    try:
        # Create Socket.IO client
        sio = socketio.Client()
        
        connected = False
        error_occurred = False
        
        @sio.event
        def connect():
            nonlocal connected
            connected = True
            print("✅ Socket.IO connected successfully")
        
        @sio.event
        def disconnect():
            print("Socket.IO disconnected")
        
        @sio.event
        def error(data):
            nonlocal error_occurred
            error_occurred = True
            print(f"❌ Socket.IO error: {data}")
        
        # Connect to server
        sio.connect('http://localhost:6513')
        
        # Wait a bit for connection
        time.sleep(2)
        
        if connected and not error_occurred:
            print("✅ Socket.IO connection test passed")
            sio.disconnect()
            return True
        else:
            print("❌ Socket.IO connection test failed")
            sio.disconnect()
            return False
            
    except Exception as e:
        print(f"❌ Socket.IO test failed: {e}")
        return False

def main():
    """Run all tests."""
    print("🧪 Testing Socket.IO fixes...")
    print("=" * 50)
    
    tests = [
        ("Health Endpoint", test_health_endpoint),
        ("Socket.IO Connection", test_socketio_connection),
    ]
    
    results = []
    for test_name, test_func in tests:
        print(f"\n🔍 Testing: {test_name}")
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"❌ {test_name} test crashed: {e}")
            results.append((test_name, False))
    
    print("\n" + "=" * 50)
    print("📊 Test Results:")
    print("=" * 50)
    
    passed = 0
    total = len(results)
    
    for test_name, result in results:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{status}: {test_name}")
        if result:
            passed += 1
    
    print(f"\n🎯 Summary: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests passed! The Socket.IO fixes are working correctly.")
    else:
        print("⚠️  Some tests failed. Please check the issues above.")

if __name__ == "__main__":
    main() 